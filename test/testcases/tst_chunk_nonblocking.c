/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests chunking feature when using nonblocking APIs and one of
 * the processes makes no call to the API.
 *
 * Contributed by Danqing Wu.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */

#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>


#define DIM_LEN 8

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int err, nerrs=0, ncid, dimid, varid, rank, req, verbose=0;
    int vals[DIM_LEN] = {-1, -2, -3, -4, -5, -6, -7, -8};
    MPI_Offset start, count;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* enable chunking */
    MPI_Info_set(info, "nc_chunking", "enable");

    /* chunking is supported only when MPI-IO driver is used */
    MPI_Info_set(info, "nc_pncio", "disable");

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "x", DIM_LEN, &dimid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 1, &dimid, &varid);
    CHECK_ERR

    err = ncmpi_enddef(ncid);
    CHECK_ERR

    if (rank == 0)
    {
        start = 0;
        count = DIM_LEN;
        err = ncmpi_iput_vara_int(ncid, varid, &start, &count, vals, &req);
        CHECK_ERR
    }
    else
        req = NC_REQ_NULL;

    if (verbose) printf("rank = %d, before ncmpi_wait_all\n", rank);
    err = ncmpi_wait_all(ncid, 1, &req, NULL);
    CHECK_ERR
    if (verbose) printf("rank = %d, after ncmpi_wait_all\n", rank);

    err = ncmpi_close(ncid);

    return nerrs;
}

int main(int argc, char **argv)
{
    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};

    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 0; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 0; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "chunking using nonblocking APIs", opt, test_io);

    MPI_Finalize();

    return err;
}
