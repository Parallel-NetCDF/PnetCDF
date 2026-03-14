/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if one can get the number of record variables and fixed-
 * sized variables correctly. It first defines some number of fixed-size and
 * record variables and then calls the APIs
 *     ncmpi_inq_num_rec_vars() and ncmpi_inq_num_fix_vars()
 * to verify if the numbers are correct.
 *
 * The compile and run commands are given below. This program is to be run on
 * one MPI process.
 *
 *    % mpicc -g -o inq_num_vars inq_num_vars.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 inq_num_vars testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

static
int check_num_vars(int  ncid,
                   int  expected_nvars,
                   int  expected_num_rec_vars,
                   int  expected_num_fix_vars)
{
    int err, nerrs=0, nvars, num_rec_vars, num_fix_vars;

    /* NULL argument test */
    err = ncmpi_inq_nvars(ncid, NULL); CHECK_ERR
    err = ncmpi_inq_num_rec_vars(ncid, NULL); CHECK_ERR
    err = ncmpi_inq_num_fix_vars(ncid, NULL); CHECK_ERR

    err = ncmpi_inq_nvars(ncid, &nvars); CHECK_ERR
    err = ncmpi_inq_num_rec_vars(ncid, &num_rec_vars); CHECK_ERR
    err = ncmpi_inq_num_fix_vars(ncid, &num_fix_vars); CHECK_ERR

    if (nvars != expected_nvars) {
        printf("Error at line %d in %s: expecting %d number of variables defined, but got %d\n",
        __LINE__,__FILE__,expected_nvars, nvars);
        nerrs++;
    }
    if (num_rec_vars != expected_num_rec_vars) {
        printf("Error at line %d in %s: expecting %d number of record variables defined, but got %d\n",
        __LINE__,__FILE__,expected_num_rec_vars, num_rec_vars);
        nerrs++;
    }
    if (num_fix_vars != expected_num_fix_vars) {
        printf("Error at line %d in %s: expecting %d number of fixed-size variables defined, but got %d\n",
        __LINE__,__FILE__,expected_num_fix_vars, num_fix_vars);
        nerrs++;
    }
    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    char val[MPI_MAX_INFO_VAL];
    int nerrs=0, err, flag, ncid, varid[7], dimid[3];

    /* check whether burst buffering is enabled */
    MPI_Info_get(info, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, val, &flag);
    if (flag && strcasecmp(val, "enable") == 0 &&
        (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC))
        /* does not work for NetCDF4 files when burst-buffering is enabled */
        return 0;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* define dimension and variable */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y",       2,            &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X",       10,           &dimid[2]); CHECK_ERR

    nerrs = 0;

    err = ncmpi_def_var(ncid, "REC_VAR_1", NC_INT, 1, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "REC_VAR_2", NC_INT, 3, dimid, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "REC_VAR_3", NC_INT, 2, dimid, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "REC_VAR_4", NC_INT, 1, dimid, &varid[3]); CHECK_ERR

    nerrs += check_num_vars(ncid, 4, 4, 0);

    err = ncmpi_def_var(ncid, "FIX_VAR_1", NC_INT, 2, dimid+1, &varid[4]); CHECK_ERR
    err = ncmpi_def_var(ncid, "FIX_VAR_2", NC_INT, 1, dimid+1, &varid[5]); CHECK_ERR
    err = ncmpi_def_var(ncid, "FIX_VAR_3", NC_INT, 1, dimid+2, &varid[6]); CHECK_ERR

    nerrs += check_num_vars(ncid, 7, 4, 3);

    /* set fill mode, so ncmpidiff can compare 2 output files without error */
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* NULL argument test */
    err = ncmpi_inq_ndims(ncid, NULL); CHECK_ERR
    err = ncmpi_inq_nvars(ncid, NULL); CHECK_ERR
    err = ncmpi_inq_natts(ncid, NULL); CHECK_ERR
    err = ncmpi_inq_unlimdim(ncid, NULL); CHECK_ERR

    nerrs += check_num_vars(ncid, 7, 4, 3);

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the file for reading --------------------------------------------*/
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    nerrs += check_num_vars(ncid, 7, 4, 3);

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
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "number of variables", opt, test_io);

    MPI_Finalize();

    return err;
}
