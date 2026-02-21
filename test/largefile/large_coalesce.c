/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Given two large nonblocking requests and they are actually contiguous in
 * fileview and buffer type, this program tests whether filetype and buftype
 * coalescing is skipped when calling MPI_Type_create_hindexed internally.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

#define FOUR_G 4294967296LL
#define TWO_G  2147483648LL
#define ONE_G  1073741824LL

static
int test_io_nc5(const char *out_path,
                MPI_Info    global_info)
{
    unsigned char *buf;
    int rank, nprocs, color, err, nerrs=0;
    int ncid, varid, dimid[2], req[3], st[3];
    MPI_Offset start[2], count[2];
    MPI_Info info;
    size_t i;
    MPI_Comm comm=MPI_COMM_NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    color = 1;

    if (nprocs > 2) {
        /* run on 2 ranks only, as this test allocates memory > 4GB per rank */
        /* split MPI_COMM_WORLD based on 'color' and use the same rank order */
        color = (rank < 2) ? 1 : 0;
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    }
    else
        comm = MPI_COMM_WORLD;

    if (!color) goto err_out;

    buf = (unsigned char*) calloc(TWO_G+1024,1);
    if (buf == NULL) {
        printf("malloc failed for size "OFFFMT"\n", TWO_G+1024);
        MPI_Finalize();
        return 1;
    }

    MPI_Info_dup(global_info, &info);
    MPI_Info_set(info, "romio_cb_write", "enable");
    MPI_Info_set(info, "romio_ds_read", "disable"); /* run slow without it */

    /* silence internal debug messages */
    setenv("PNETCDF_SAFE_MODE", "0", 1);

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "NPROCS", nprocs, &dimid[0]);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "X", TWO_G+1024, &dimid[1]);
    CHECK_ERR

    /* define a big 1D variable of ubyte type */
    err = ncmpi_def_var(ncid, "big_var", NC_UBYTE, 2, dimid, &varid);
    CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    /* now we are in data mode */
    for (i=0; i<20; i++) buf[ONE_G-10+i] = 'a'+i;
    for (i=0; i<20; i++) buf[TWO_G-10+i] = 'A'+i;

    start[0] = rank;
    count[0] = 1;

    start[1] = 0;
    count[1] = 10;
    err = ncmpi_iput_vara_uchar(ncid, varid, start, count, buf, &req[0]);
    CHECK_ERR

    /* 2nd request is not contiguous from the first */
    start[1] = 1024;
    count[1] = ONE_G-1024;
    err = ncmpi_iput_vara_uchar(ncid, varid, start, count, buf+1024, &req[1]);
    CHECK_ERR

    /* make file access and write buffer of 3rd request contiguous from the 2nd
     * request to check whether the internal fileview and buftype coalescing
     * are skipped */
    start[1] = ONE_G;
    count[1] = ONE_G+1024;
    err = ncmpi_iput_vara_uchar(ncid, varid, start, count, buf+ONE_G, &req[2]);
    CHECK_ERR

    err = ncmpi_wait_all(ncid, 3, req, st);
    CHECK_ERR

    /* read back to check contents */
    start[1] = ONE_G-10;
    count[1] = 20;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf);
    CHECK_ERR
    for (i=0; i<20; i++) {
        if (buf[i] != 'a'+i) {
            printf("%d (at line %d): expect buf["OFFFMT"]=%zd but got %d\n",
                   rank, __LINE__, ONE_G-10+i, i+'a', buf[i]);
            nerrs++;
        }
    }

    /* read back to check contents */
    start[1] = TWO_G-10;
    count[1] = 20;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf);
    CHECK_ERR
    for (i=0; i<20; i++) {
        if (buf[i] != 'A'+i) {
            printf("%d (at line %d): expect buf["OFFFMT"]=%zd but got %d\n",
                   rank, __LINE__, TWO_G-10+i, i+'A', buf[i]);
            nerrs++;
        }
    }

    /* test the same pattern but for iget */
    for (i=0; i<TWO_G+1024; i++) buf[i] = 0;

    start[1] = 0;
    count[1] = 10;
    err = ncmpi_iget_vara_uchar(ncid, varid, start, count, buf, &req[0]);
    CHECK_ERR

    start[1] = 1024;
    count[1] = ONE_G-1024;
    err = ncmpi_iget_vara_uchar(ncid, varid, start, count, buf+1024, &req[1]);
    CHECK_ERR

    start[1] = ONE_G;
    count[1] = ONE_G+1024;
    err = ncmpi_iget_vara_uchar(ncid, varid, start, count, buf+ONE_G, &req[2]);
    CHECK_ERR

    err = ncmpi_wait_all(ncid, 3, req, st);
    CHECK_ERR

    for (i=0; i<20; i++) {
        if (buf[ONE_G-10+i] != 'a'+i) {
            printf("%d (at line %d): expect buf["OFFFMT"]=%zd but got %d\n",
                   rank, __LINE__, ONE_G-10+i, i+'a', buf[i]);
            nerrs++;
        }
    }

    for (i=0; i<20; i++) {
        if (buf[TWO_G-10+i] != 'A'+i) {
            printf("%d (at line %d): expect buf["OFFFMT"]=%zd but got %d\n",
                   rank, __LINE__, TWO_G-10+i, i+'A', buf[i]);
            nerrs++;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if open for reading header */
    err = ncmpi_open(comm, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    free(buf);

    MPI_Info_free(&info);

err_out:
    if (comm != MPI_COMM_WORLD && comm != MPI_COMM_NULL)
        MPI_Comm_free(&comm);

    return nerrs;
}

static
int test_io_nc4(const char *out_path,
                MPI_Info    global_info)
{
    unsigned char *buf;
    int rank, nprocs, color, err, nerrs=0;
    int ncid, varid, dimid[2];
    MPI_Offset start[2], count[2];
    MPI_Info info;
    size_t i;
    MPI_Comm comm=MPI_COMM_NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    color = 1;

    if (nprocs > 2) {
        /* run on 2 ranks only, as this test allocates memory > 4GB per rank */
        /* split MPI_COMM_WORLD based on 'color' and use the same rank order */
        color = (rank < 2) ? 1 : 0;
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    }
    else
        comm = MPI_COMM_WORLD;

    if (!color) goto err_out;

    buf = (unsigned char*) calloc(TWO_G+1024,1);
    if (buf == NULL) {
        printf("malloc failed for size "OFFFMT"\n", TWO_G+1024);
        MPI_Finalize();
        return 1;
    }

    MPI_Info_dup(global_info, &info);
    MPI_Info_set(info, "romio_cb_write", "enable");
    MPI_Info_set(info, "romio_ds_read", "disable"); /* run slow without it */

    /* silence internal debug messages */
    setenv("PNETCDF_SAFE_MODE", "0", 1);

    /* Test for NetCDF 4 first as ncvalidator checks only read classic files */

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "NPROCS", nprocs, &dimid[0]);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "X", TWO_G+1024, &dimid[1]);
    CHECK_ERR

    /* define a big 1D variable of ubyte type */
    err = ncmpi_def_var(ncid, "big_var", NC_UBYTE, 2, dimid, &varid);
    CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    /* now we are in data mode */
    for (i=0; i<20; i++) buf[ONE_G-10+i] = 'a'+i;
    for (i=0; i<20; i++) buf[TWO_G-10+i] = 'A'+i;

    start[0] = rank;
    count[0] = 1;

    start[1] = 0;
    count[1] = 10;
    err = ncmpi_put_vara_uchar_all(ncid, varid, start, count, buf);
    CHECK_ERR

    /* 2nd request is not contiguous from the first */
    start[1] = 1024;
    count[1] = ONE_G-1024;
    err = ncmpi_put_vara_uchar_all(ncid, varid, start, count, buf+1024);
    CHECK_ERR

    /* make file access and write buffer of 3rd request contiguous from the 2nd
     * request to check whether the internal fileview and buftype coalescing
     * are skipped */
    start[1] = ONE_G;
    count[1] = ONE_G+1024;
    err = ncmpi_put_vara_uchar_all(ncid, varid, start, count, buf+ONE_G);
    CHECK_ERR

    start[1] = 0;
    count[1] = 10;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf);
    CHECK_ERR

    start[1] = 1024;
    count[1] = ONE_G-1024;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf+1024);
    CHECK_ERR

    start[1] = ONE_G;
    count[1] = ONE_G+1024;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf+ONE_G);
    CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if open to read header fine */
    err = ncmpi_open(comm, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    free(buf);

    MPI_Info_free(&info);

err_out:
    if (comm != MPI_COMM_WORLD && comm != MPI_COMM_NULL)
        MPI_Comm_free(&comm);

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int err;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    if (format == NC_FORMAT_NETCDF4)
        return test_io_nc4(out_path, info);
    else
        return test_io_nc5(out_path, info);
}

int main(int argc, char **argv) {

    int err;
#ifdef ENABLE_NETCDF4
    int formats[] = {NC_FORMAT_NETCDF4, NC_FORMAT_64BIT_DATA};
#else
    int formats[] = {NC_FORMAT_64BIT_DATA};
#endif
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 0; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "skip filetype buftype coalesce", opt, test_io);

    MPI_Finalize();

    return err;
}
