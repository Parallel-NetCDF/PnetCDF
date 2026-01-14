/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program check if error code NC_EMAXDIMS can be returned correctly, when
 * defining a variable with more than NC_MAX_VAR_DIMS dimensions.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_max_var_dims tst_max_var_dims.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_max_var_dims testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <limits.h> /* INT_MAX */
#include <pnetcdf.h>

#include <testutils.h>

#if NC_MAX_VAR_DIMS < INT_MAX

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int err, nerrs=0, ncid, varid, *dimid;
    size_t i;

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* define dimensions */
    dimid = (int*) malloc(sizeof(int) * (NC_MAX_VAR_DIMS+2));
    err = ncmpi_def_dim(ncid, "dim0", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", 1, &dimid[1]); CHECK_ERR

    for (i=2; i<(size_t)NC_MAX_VAR_DIMS+2; i++) dimid[i] = dimid[1];

    /* define variables */
    err = ncmpi_def_var(ncid, "v0", NC_INT, NC_MAX_VAR_DIMS+1, &dimid[0], &varid);
    EXP_ERR(NC_EMAXDIMS)

    err = ncmpi_def_var(ncid, "v1", NC_INT, NC_MAX_VAR_DIMS+1, &dimid[1], &varid);
    EXP_ERR(NC_EMAXDIMS)

    err = ncmpi_set_fill(ncid, NC_NOFILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    free(dimid);

    return nerrs;
}

int main(int argc, char **argv) {

    int err=0;
    loop_opts opt;

    MPI_Init(&argc, &argv);


    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "checking NC_MAX_VAR_DIMS", opt, test_io);

    MPI_Finalize();

    return err;
}
#else
int main(int argc, char **argv) {
    return 0;
}
#endif
