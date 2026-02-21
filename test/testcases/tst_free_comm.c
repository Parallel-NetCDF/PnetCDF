/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests whether PnetCDF duplicates MPI communicator and MPI info
 * object correctly, so the user supplied communicator and info objects can be
 * freed by users after ncmpi_create and ncmpi_open.
 *
 *  % mpiexec -n 4 tst_free_comm
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h> /* strcpy(), strncpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int nerrs=0, err, ncid;
    MPI_Comm comm=MPI_COMM_NULL;

    /* duplicate MPI_COMM_WORLD */
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a file */
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    MPI_Comm_free(&comm); comm = MPI_COMM_NULL;

    err = ncmpi_close(ncid); CHECK_ERR

    /* duplicate MPI_COMM_WORLD */
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    /* open the file */
    err = ncmpi_open(comm, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR

    MPI_Comm_free(&comm); comm = MPI_COMM_NULL;

    err = ncmpi_close(ncid); CHECK_ERR

    if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);

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
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 0; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "freeing MPI communicator", opt, test_io);

    MPI_Finalize();

    return err;
}
