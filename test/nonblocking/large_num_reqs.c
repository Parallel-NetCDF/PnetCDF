/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * This program tests a large number of nonblocking requests (larger than
 * NC_REQUEST_CHUNK, the constant used to grow nonblocking put and get
 * queues.
 *
 *********************************************************************/
/*  $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <assert.h>
#include <pnetcdf.h>

#include <testutils.h>

#define FILE_NAME "testfile.nc"
#define NUM_REQS 1100   /* a number greater than NC_REQUEST_CHUNK */

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, ncid, dimid[2], varid, err, nerrs=0, rank, nprocs;
    int *buf, *req, *status;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
#define STRESS_ROMIO
#ifdef STRESS_ROMIO
    /* all processes write to the same file region. This can vigorously test
     * ROMIO for receiving write requests into the same local buffer at each
     * I/O aggregator.
     */
    err = ncmpi_def_dim(ncid, "X", 2,            &dimid[1]); CHECK_ERR
#else
    /* Writes from all processes are not overlapped */
    err = ncmpi_def_dim(ncid, "X", 2 * nprocs,   &dimid[1]); CHECK_ERR
#endif
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    req = (int*) malloc(sizeof(int) * NUM_REQS * 2);
    status = req + NUM_REQS;

    buf = (int*) calloc(NUM_REQS * 6,  sizeof(int));

    count[0] = 3; count[1] = 2;
    start[0] = 0;
#ifdef STRESS_ROMIO
    start[1] = 0;
#else
    start[1] = rank * 2;
#endif

    for (i=0; i<NUM_REQS; i++) {
        err = ncmpi_iput_vara_int(ncid, varid, start, count,
                                  &buf[i*6], &req[i]); CHECK_ERR
        start[0] += 3;
    }

    if (coll_io)
        err = ncmpi_wait_all(ncid, NUM_REQS, req, status);
    else
        err = ncmpi_wait(ncid, NUM_REQS, req, status);
    CHECK_ERR

    /* check each iput status */
    for (i=0; i<NUM_REQS; i++) {
        err = status[i];
        CHECK_ERR
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid); CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    start[0] = 0;
    for (i=0; i<NUM_REQS; i++) {
        err = ncmpi_iget_vara_int(ncid, varid, start, count,
                                  &buf[i*6], &req[i]); CHECK_ERR
        start[0] += 3;
    }

    if (coll_io)
        err = ncmpi_wait_all(ncid, NUM_REQS, req, status);
    else
        err = ncmpi_wait(ncid, NUM_REQS, req, status);
    CHECK_ERR

    /* check each iget status */
    for (i=0; i<NUM_REQS; i++) {
        err = status[i];
        CHECK_ERR
    }

    free(buf);
    free(req);
    err = ncmpi_close(ncid); CHECK_ERR

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
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "large number of iput/iget", opt, test_io);

    MPI_Finalize();

    return err;
}
