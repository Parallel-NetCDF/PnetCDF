/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests the use of nonblocking API.
 * The write buffer is a 2D array of size NY x NX
 * It writes the 2nd row of the memory buffer to the 1st row of the variable
 * array in file. Then it writes the 1st row of the memory buffer to the
 * 2nd row of the variable array in file.
 *
 * The expected reults from the output file contents are:
 * (when running on 1 MPI process)
 *
 *  % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 *         Y = UNLIMITED ; // (2 currently)
 *         X = 5 ;
 *    variables:
 *         int VAR(Y, X) ;
 *    data:
 *
 *    var =
 *      1, 1, 1, 1, 1,
 *      0, 0, 0, 0, 0 ;
 *    }
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 4
#define NX 5

int tst_iput(const char *out_path,
             const char *in_path, /* ignored */
             int         format,
             int         coll_io,
             MPI_Info    info)
{
    int i, j, err, ncid, varid, dimids[2], req[2], st[2], nerrs=0;
    int rank, nprocs, buf[NY+1][NX];
    MPI_Offset start[2], count[2];

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_FATAL_ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX,    &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_sync(ncid); CHECK_ERR
    }

    /* initialize the contents of the array */
    for (j=0; j<NY+1; j++) for (i=0; i<NX; i++) buf[j][i] = j;

    start[0] = 2*rank; start[1] = 0;
    count[0] = 1;      count[1] = NX;

    /* call nonblocking API */
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf[1], &req[0]);
    CHECK_ERR

    start[0] += 1;
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf[0], &req[1]);
    CHECK_ERR

    st[0] = st[1] = NC_NOERR;

    if (coll_io) {
        err = ncmpi_wait_all(ncid, 2, req, st);
        CHECK_ERR
    }
    else {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
        err = ncmpi_wait(ncid, 2, req, st);
        CHECK_ERR
    }
    err = st[0]; CHECK_ERR
    err = st[1]; CHECK_ERR

    /* check if the contents of buf are altered */
    for (j=0; j<NY; j++)
        for (i=0; i<NX; i++)
            if (buf[j][i] != j) {
                printf("Error at line %d in %s: buf[%d][%d]=%d != %d\n",
                __LINE__,__FILE__,j,i,buf[j][i],j);
                nerrs++;
                goto fn_exit;
            }

    /* check if root process can write to file header in data mode */
    err = ncmpi_rename_var(ncid, varid, "VAR"); CHECK_ERR

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_FATAL_ERR

    err = ncmpi_inq_varid(ncid, "VAR", &varid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* initialize the contents of the array to a different value */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = -1;

    /* read back variable */
    start[0] = 2*rank; start[1] = 0;
    count[0] = 2;      count[1] = NX;
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf[0]);
    else
        err = ncmpi_get_vara_int(ncid, varid, start, count, buf[0]);
    CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if the contents of buf are expected */
    for (j=0; j<2; j++) {
        int val = (j == 0) ? 1 : 0;
        for (i=0; i<NX; i++)
            if (buf[j][i] != val) {
                printf("Error: unexpected read buf[%d][%d]=%d, should be %d\n",
                       j,i,buf[j][i],val);
                nerrs++;
                goto fn_exit;
            }
    }

fn_exit:
    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int nerrs=0;

    nerrs = tst_iput(out_path, in_path, format, coll_io, info);
    if (nerrs > 0) goto err_out;

    /* disable PnetCDF internal buffering */
    MPI_Info_set(info, "nc_ibuf_size", "0");

    nerrs = tst_iput(out_path, in_path, format, coll_io, info);
    if (nerrs > 0) goto err_out;

err_out:
    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};

    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "ncmpi_iput_vara_int()", opt, test_io);

    MPI_Finalize();

    return err;
}
