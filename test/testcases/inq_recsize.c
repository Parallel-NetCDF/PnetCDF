/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if one can get the size of record block correctly. The
 * record block size is the sum of individual record of all record variables.
 * It first defines some number of record and fixed-size variables and then
 * calls the API ncmpi_inq_recsize() and verify if the numbers are correct.
 *
 * The compile and run commands are given below. This program is to be run on
 * one MPI process.
 *
 *    % mpicc -g -o inq_recsize inq_recsize.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 inq_recsize testfile.nc
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
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char val[MPI_MAX_INFO_VAL];
    int nerrs=0, err, flag, ncid, varid[7], dimid[3];
    MPI_Offset expected_recsize, recsize;

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
    expected_recsize = 0;

    /* define some record variables */
    err = ncmpi_def_var(ncid, "REC_VAR_1", NC_INT, 1, dimid, &varid[0]); CHECK_ERR
    expected_recsize += sizeof(int);
    err = ncmpi_def_var(ncid, "REC_VAR_2", NC_INT, 3, dimid, &varid[1]); CHECK_ERR
    expected_recsize += 2 * 10 * sizeof(int);
    err = ncmpi_def_var(ncid, "REC_VAR_3", NC_INT, 2, dimid, &varid[2]); CHECK_ERR
    expected_recsize += 2 * sizeof(int);
    err = ncmpi_def_var(ncid, "REC_VAR_4", NC_INT, 1, dimid, &varid[3]); CHECK_ERR
    expected_recsize += sizeof(int);

    /* define some fixed-size variables */
    err = ncmpi_def_var(ncid, "FIX_VAR_1", NC_INT, 2, dimid+1, &varid[4]); CHECK_ERR
    err = ncmpi_def_var(ncid, "FIX_VAR_2", NC_INT, 1, dimid+1, &varid[5]); CHECK_ERR
    err = ncmpi_def_var(ncid, "FIX_VAR_3", NC_INT, 1, dimid+2, &varid[6]); CHECK_ERR

    /* set fill mode, so ncmpidiff can compare 2 output files without error */
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* NULL argument test */
    err = ncmpi_inq_recsize(ncid, NULL); CHECK_ERR

    err = ncmpi_inq_recsize(ncid, &recsize); CHECK_ERR
    if (expected_recsize != recsize) {
        printf("Error at line %d in %s: expecting record size "OFFFMT" but got "OFFFMT"\n",
        __LINE__,__FILE__,expected_recsize, recsize);
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
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "inquiring record size", opt, test_io);

    MPI_Finalize();

    return err;
}
