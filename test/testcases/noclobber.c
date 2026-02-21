/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests if PnetCDF can return the right error code NC_EEXIST
 * when create mode NC_NOCLOBBER is used and the file exists.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strdup() */
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int err, nerrs=0, ncid;

    MPI_Barrier(MPI_COMM_WORLD);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a file if it does not exist */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* now the file exists, test if PnetCDF can return correct error code */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_NOCLOBBER, info, &ncid);
    EXP_ERR(NC_EEXIST) /* err == NC_EOFILE */

    MPI_Barrier(MPI_COMM_WORLD);

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
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "NC_NOCLOBBER and NC_EEXIST", opt, test_io);

    MPI_Finalize();

    return err;
}
