/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program adds two new variables to an existing netCDF file.
 * It is used to test if PnetCDF can correctly calculate the file offsets
 * for the two new variables, in particular for files that align the
 * fixed-size variables to a boundary larger than 4 bytes, for instance
 * a file created by PnetCDF with defaut alignment of 512 bytes.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o add_var add_var.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 add_var testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char var_name[256];
    int i, nvars, err, nerrs=0;
    int ncid, varid, dimid[2];
    MPI_Offset prev_off;

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "dim_1", 5, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim_2", 4, &dimid[1]); CHECK_ERR

    /* define a bunch of variables */
    for (i=0; i<10; i++) {
        sprintf(var_name, "var_%d", i);
        err = ncmpi_def_var(ncid, var_name, NC_INT, 2, dimid, &varid); CHECK_ERR
        err = ncmpi_def_var_fill(ncid, varid, 0, NULL); CHECK_ERR
    }
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* re-enter define mode */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* add 2 new dimensions */
    err = ncmpi_def_dim(ncid, "new_dim_1", 5, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "new_dim_2", 4, &dimid[1]); CHECK_ERR

    /* add 2 new variables */
    err = ncmpi_def_var(ncid, "new_var1", NC_INT,   2, dimid, &varid); CHECK_ERR
    err = ncmpi_def_var_fill(ncid, varid, 0, NULL); CHECK_ERR
    err = ncmpi_def_var(ncid, "new_var2", NC_FLOAT, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_def_var_fill(ncid, varid, 0, NULL); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_nvars(ncid, &nvars); CHECK_ERR
    if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC) {
        err = ncmpi_inq_varoffset(ncid, 0, &prev_off); EXP_ERR(NC_ENOTSUPPORT)
    }
    else {
        err = ncmpi_inq_varoffset(ncid, 0, &prev_off); CHECK_ERR
        for (i=1; i<nvars; i++) {
            MPI_Offset off;
            err = ncmpi_inq_varoffset(ncid, i, &off); CHECK_ERR
            if (off < prev_off + 5*4*4) { /* each variable is of size 5*4*4 bytes */
                err = ncmpi_inq_varname(ncid, i, var_name); CHECK_ERR
                printf("Error at line %d in %s: variable %s offset is set incorrectly\n",
                       __LINE__,__FILE__,var_name);
                nerrs++;
            }
            prev_off = off;
        }
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
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 1; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "checking offsets of new variables", opt, test_io);

    MPI_Finalize();

    return err;
}
