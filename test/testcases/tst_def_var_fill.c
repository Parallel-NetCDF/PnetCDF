/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests fill mode for individual variables.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_def_var_fill tst_def_var_fill.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 tst_def_var_fill testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 8
#define NX 5

static int
tst_fmt(const char *out_path, int format, int coll_io, MPI_Info info)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, fmt, varid[2], dimid[2], expect, *buf;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* define dimension */
    err = ncmpi_def_dim(ncid, "Y", NY,        &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs, &dimid[1]); CHECK_ERR

    /* define 2 variables of size NY x NX */
    err = ncmpi_def_var(ncid, "var_nofill", NC_INT, 2, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_fill",   NC_INT, 2, dimid, &varid[1]); CHECK_ERR
    /* set fill modes for the variables */
    err = ncmpi_def_var_fill(ncid, NC_GLOBAL,1, NULL); EXP_ERR(NC_EGLOBAL)
    err = ncmpi_def_var_fill(ncid, varid[0], 1, NULL); CHECK_ERR
    err = ncmpi_def_var_fill(ncid, varid[1], 0, NULL); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

#ifdef STRONGER_CONSISTENCY
    err = ncmpi_sync(ncid); CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);
    err = ncmpi_sync(ncid); CHECK_ERR
#endif

    /* initialize I/O buffer */
    buf = (int*) malloc(sizeof(int) * NY*NX);

    for (i=0; i<NY*NX; i++) buf[i] = rank+5;

    /* write a subarray to each variable */
    start[0] = 0;
    start[1] = NX*rank+2;
    count[0] = NY;
    count[1] = 2;
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf);
    else
        err = ncmpi_put_vara_int(ncid, varid[0], start, count, buf);
    CHECK_ERR
    /* check if user put buffer contents altered */
    for (i=0; i<NY*NX; i++) {
        if (buf[i] != rank+5) {
            printf("Error in %s line %d: put buf[%d] altered from %d to %d\n",
                   __FILE__,__LINE__, i, rank+5, buf[i]);
            nerrs++;
        }
    }
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf);
    else
        err = ncmpi_put_vara_int(ncid, varid[1], start, count, buf);
    CHECK_ERR
    /* check if user put buffer contents altered */
    for (i=0; i<NY*NX; i++) {
        if (buf[i] != rank+5) {
            printf("Error in %s line %d: put buf[%d] altered from %d to %d\n",
                   __FILE__,__LINE__, i, rank+5, buf[i]);
            nerrs++;
        }
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* reopen the file and read data back */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_format(ncid, &fmt); CHECK_ERR
    if (fmt != format) {
        printf("Error at line %d of %s: expect %s but got %s\n",
               __LINE__,__FILE__,pnc_fmt_string(format),pnc_fmt_string(fmt));
        nerrs++;
    }

    /* inquire variabe IDs */
    err = ncmpi_inq_varid(ncid, "var_nofill", &varid[0]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "var_fill",   &varid[1]); CHECK_ERR

    /* read the subarray written by process (rank+1)%nproc */
    start[0] = 0;
    start[1] = NX*((rank+1)%nprocs);
    count[0] = NY;
    count[1] = NX;
    for (i=0; i<NY*NX; i++) buf[i] = -1;
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid[0], start, count, buf);
    else
        err = ncmpi_get_vara_int(ncid, varid[0], start, count, buf);
    CHECK_ERR

    /* check contents of variable var_nofill */
    expect = (rank+1)%nprocs + 5;
    for (i=0; i<NY; i++) for (j=0; j<NX; j++) {
        if (2 <= j && j < 4) {
            if (buf[i*NX+j] != expect) {
                printf("Error in %s line %d: expect get buf[%d]=%d but got %d\n",
                       __FILE__,__LINE__,i*NX+j, expect, buf[i*NX+j]);
                nerrs++;
            }
        }
        else if (buf[i*NX+j] == NC_FILL_INT) /* not an error */
            /* content of buf[i*NX+j] can be any value */
            printf("Warning in %s line %d: get buf[%d] same as NC_FILL_INT\n",
                   __FILE__,__LINE__,i*NX+j);
    }

    /* read the subarray written by process (rank+1)%nproc */
    for (i=0; i<NY*NX; i++) buf[i] = -1;
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid[1], start, count, buf);
    else
        err = ncmpi_get_vara_int(ncid, varid[1], start, count, buf);
    CHECK_ERR

    /* check contents of variable var_fill */
    for (i=0; i<NY; i++) for (j=0; j<NX; j++) {
        expect = NC_FILL_INT;
        if (2 <= j && j< 4) expect = (rank+1)%nprocs + 5;
        if (buf[i*NX+j] != expect) {
            printf("Error in %s line %d: expect get buf[%d]=%d but got %d\n",
                   __FILE__,__LINE__,i*NX+j, expect, buf[i*NX+j]);
            nerrs++;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char val[MPI_MAX_INFO_VAL];
    int nerrs=0, flag;

    /* check whether burst buffering is enabled */
    MPI_Info_get(info, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, val, &flag);
    if (flag && strcasecmp(val, "enable") == 0 &&
        (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC))
        /* does not work for NetCDF4 files when burst-buffering is enabled */
        return 0;

    nerrs = tst_fmt(out_path, format, coll_io, info);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "def_var_fill", opt, test_io);

    MPI_Finalize();

    return err;
}
