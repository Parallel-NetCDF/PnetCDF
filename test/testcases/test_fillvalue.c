/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests
 * 1. if PnetCDF allows put attribute _FillValue to global variable.
 * 2. if PnetCDF can return the right error code NC_EBADTYPE when put attribute
 *    _FillValue to a non-global variable with a different data type to the
 *    variable's.
 *
 * Expected results from running command ncdump on the output file:
 *
 * % ncdump testfile.nc
 * netcdf testfile {
 * variables:
 * 	int var ;
 * 		var:_FillValue = 5678 ;
 *
 * // global attributes:
 *		:_FillValue = 1.234f ;
 * data:
 *
 *  var = _ ;
 * }
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
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
    char val[MPI_MAX_INFO_VAL];
    int err, nerrs=0, flag, ncid, varid, int_buf;
    float flt_buf;

    /* check whether burst buffering is enabled */
    MPI_Info_get(info, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, val, &flag);
    if (flag && strcasecmp(val, "enable") == 0 &&
        (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC))
        /* does not work for NetCDF4 files when burst-buffering is enabled */
        return 0;

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a file */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    flt_buf = 1.234;
    err = ncmpi_put_att(ncid, NC_GLOBAL, "_FillValue", NC_FLOAT, 1, &flt_buf);
    CHECK_ERR

    err = ncmpi_def_var(ncid, "var", NC_INT, 0, NULL, &varid);
    CHECK_ERR

    err = ncmpi_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, &flt_buf);
    EXP_ERR(NC_EBADTYPE)

    int_buf = 5678;
    err = ncmpi_put_att(ncid, varid, "_FillValue", NC_INT, 1, &int_buf);
    CHECK_ERR

    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

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
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "_FillValue for NC_GLOBAL", opt, test_io);

    MPI_Finalize();

    return err;
}
