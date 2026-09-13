/*
 *  Copyright (C) 2026, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 */

/*
 * This program tests collective write to a record variable while one of the
 * process makes a zero-length request.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

// #define STAND_ALONE

#ifndef STAND_ALONE
#include <testutils.h>
#else
#include <assert.h>
#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        fprintf(stderr,"Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,ncmpi_strerrno(err)); \
        assert(0); \
    } \
}
#endif

#define DIM_X 24

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, err, nerrs=0, rank, nprocs, verbose=0;
    int ncid, dimids[2], varid, buf[DIM_X], unlimdimid, req=0;
    MPI_Offset start[2], count[2], num_rec=0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err  = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, dimids); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", DIM_X, &dimids[1]); CHECK_ERR

    /* create a record variable of type NC_INT */
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* Write some records of var data. */
    count[0] = 1;
    count[1] = DIM_X / nprocs;
    start[0] = 0;
    start[1] = count[1] * rank;

    /* The last process makes a zero-length request */
    if (rank == nprocs - 1)
        count[0] = count[1] = 0;

    if (verbose)
        printf("%d: start %lld %lld count %lld %lld\n", rank,
               start[0], start[1], count[0], count[1]);

    for (i=0; i<DIM_X; i++) buf[i] = rank;

    /* test nonblocking API */
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf, &req); CHECK_ERR

    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, NULL);
    else
        err = ncmpi_wait(ncid, 1, &req, NULL);
    CHECK_ERR

    /* test blocking API */
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf);
    else
        err = ncmpi_put_vara_int(ncid, varid, start, count, buf);
    CHECK_ERR

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    err  = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }
    err = ncmpi_inq_varid(ncid, "var", &varid); CHECK_ERR

    /* Check the number of records in the file */
    err = ncmpi_inq_unlimdim(ncid, &unlimdimid); CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, unlimdimid, &num_rec); CHECK_ERR

    if (num_rec != 1) {
        fprintf(stderr, "Error: expect 1 record, but got %lld\n", num_rec);
        nerrs++;
    }

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

#ifdef STAND_ALONE
int main(int argc, char **argv)
{
    int err;

    MPI_Init(&argc, &argv);

    err = test_io("testfile.nc", NULL, NC_FORMAT_CLASSIC, 1, MPI_INFO_NULL);

    MPI_Finalize();

    return err;
}
#else
int main(int argc, char **argv) {

    int err;
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 2;    /* enable and disable intra-node aggregation */
    opt.drv      = 2;    /* test GIO and MPI-IO driver */
    opt.ibuf     = 2;    /* enable and disable hint nc_ibuf_size */
    opt.bb       = 2;    /* enable and disable burst-buffering feature */
    opt.mod      = 2;    /* collective and independent data mode */
    opt.hdr_diff = true; /* run ncmpidiff for file header */
    opt.var_diff = true; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "only one record variable", opt, test_io);

    MPI_Finalize();

    return err;
}
#endif
