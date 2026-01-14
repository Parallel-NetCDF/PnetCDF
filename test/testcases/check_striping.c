/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests whether the file striping size and count retrieved from
 * MPI-IO hints are consistent among all MPI processes.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o get_striping get_striping.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 get_striping testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <assert.h>
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
    int err, nerrs=0, ncid, fmt;
    int striping_size=0, striping_count=0, root_striping_size, root_striping_count;

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_format(ncid, &fmt); CHECK_ERR
    err = ncmpi_inq_striping(ncid, &striping_size, &striping_count);
    if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC)
        EXP_ERR(NC_ENOTSUPPORT)
    else
        CHECK_ERR

    root_striping_size  = striping_size;
    root_striping_count = striping_count;
    err = MPI_Bcast(&root_striping_size,  1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_ERR(err)
    err = MPI_Bcast(&root_striping_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_ERR(err)
    if (root_striping_size != striping_size) {
        printf("Error at line %d in %s: inconsistent striping_size (root=%d local=%d)\n",
               __LINE__,__FILE__, root_striping_size, striping_size);
        nerrs++;
    }
    if (root_striping_count != striping_count) {
        printf("Error at line %d in %s: inconsistent striping_count (root=%d local=%d)\n",
               __LINE__,__FILE__, root_striping_count, striping_count);
        nerrs++;
    }

    err = ncmpi_close(ncid); CHECK_ERR

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
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 0; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "inquire striping info", opt, test_io);

    MPI_Finalize();

    return err;
}
