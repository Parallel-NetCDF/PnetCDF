/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if the number of records is updated correctly. It first
 * writes to the 2nd record of a variable, followed by writing to the 1st
 * record. After the 2nd write, a call to ncmpi_inq_dimlen() or ncmpi_inq_dim()
 * should report seeing 2 records. Then, it does a similar test under
 * independent data mode. The same test will repeat for 3 cases:
 * 1) there is only one 1D record variable
 * 2) there is only one 3D record variable
 * 3) there are one 1D record variable and one 3D record variable
 *
 * The compile and run commands are given below. This program is to be run on
 * one MPI process.
 *
 *    % mpicc -g -o record record.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 record record.nc
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
int test_only_record_var_1D(const char *out_path, int format, MPI_Info info)
{
    int ncid, varid, dimid, buf[20], err, nerrs=0;
    MPI_Offset start[1], count[1], length;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_SELF, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* define dimension and variable */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimid); CHECK_ERR
    err = ncmpi_def_var(ncid, "REC_VAR_1D", NC_INT, 1, &dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* write the 2nd record first */
    buf[0] = 91;
    start[0] = 1; count[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* write the 1st record now */
    buf[0] = 90;
    start[0] = 0; count[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, dimid, &length); CHECK_ERR
    if (length != 2) {
        printf("Error at line %d in %s: expecting 2 records, but got "OFFFMT" record(s)\n",
        __LINE__,__FILE__,length);
        nerrs++;
    }

    /* test independent data mode */
    if (nerrs == 0 && format != NC_FORMAT_NETCDF4 &&
                      format != NC_FORMAT_NETCDF4_CLASSIC) {

        err = ncmpi_begin_indep_data(ncid); CHECK_ERR
        /* write the 4th record */
        buf[0] = 93;
        start[0] = 3; count[0] = 1;
        err = ncmpi_put_vara_int(ncid, varid, start, count, buf); CHECK_ERR

        /* write the 3rd record */
        buf[0] = 92; buf[1] = 93;
        start[0] = 2; count[0] = 2;
        err = ncmpi_put_vara_int(ncid, varid, start, count, buf); CHECK_ERR

        err = ncmpi_inq_dimlen(ncid, dimid, &length); CHECK_ERR
        if (length != 4) {
            printf("Error at line %d in %s: expecting 4 records, but got "OFFFMT" record(s)\n",
                   __LINE__,__FILE__,length);
            nerrs++;
        }
        err = ncmpi_end_indep_data(ncid); CHECK_ERR
    }
    err = ncmpi_close(ncid); CHECK_ERR
    return nerrs;
}

static
int test_only_record_var_3D(const char *out_path, int format, MPI_Info info)
{
    int i, ncid, varid, dimid[3], buf[20], err, nerrs=0;
    MPI_Offset start[3], count[3], length;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_SELF, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* define dimension and variable */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "FIX_DIM_Y", 2,          &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "FIX_DIM_X", 10,         &dimid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "REC_VAR_3D", NC_INT, 3, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    start[1] = 0; start[2] = 0; count[0] = 1; count[1] = 2; count[2] = 10;

    /* write the 2nd record first */
    for (i=0; i<20; i++) buf[i] = 91;
    start[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* write the 1st record now */
    for (i=0; i<20; i++) buf[i] = 90;
    start[0] = 0;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, dimid[0], &length); CHECK_ERR
    if (length != 2) {
        printf("Error at line %d in %s: expecting 2 records, but got "OFFFMT" record(s)\n",
        __LINE__,__FILE__,length);
        nerrs++;
    }

    /* test independent data mode */
    if (nerrs == 0 && format != NC_FORMAT_NETCDF4 &&
                      format != NC_FORMAT_NETCDF4_CLASSIC) {

        err = ncmpi_begin_indep_data(ncid); CHECK_ERR
        /* write the 4th record */
        for (i=0; i<20; i++) buf[i] = 93;
        start[0] = 3;
        err = ncmpi_put_vara_int(ncid, varid, start, count, buf); CHECK_ERR

        /* write the 3rd record */
        for (i=0; i<20; i++) buf[i] = 92;
        start[0] = 2;
        err = ncmpi_put_vara_int(ncid, varid, start, count, buf); CHECK_ERR

        err = ncmpi_inq_dimlen(ncid, dimid[0], &length); CHECK_ERR
        if (length != 4) {
            printf("Error at line %d in %s: expecting 4 records, but got "OFFFMT" record(s)\n",
                   __LINE__,__FILE__,length);
            nerrs++;
        }
        err = ncmpi_end_indep_data(ncid); CHECK_ERR
    }
    err = ncmpi_close(ncid); CHECK_ERR
    return nerrs;
}

static
int test_two_record_var(const char *out_path, int format, MPI_Info info)
{
    int i, ncid, varid[2], dimid[3], buf[20], err, nerrs=0;
    MPI_Offset start[3], count[3], length;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_SELF, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* define dimension and variable */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "FIX_DIM_Y", 2,          &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "FIX_DIM_X", 10,         &dimid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "REC_VAR_1D", NC_INT, 1, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "REC_VAR_3D", NC_INT, 3, dimid, &varid[1]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* REC_VAR_1D: write the 2nd record first */
    buf[0] = 91;
    start[0] = 1; count[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR

    /* write the 1st record now */
    buf[0] = 90;
    start[0] = 0; count[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, dimid[0], &length); CHECK_ERR
    if (length != 2) {
        printf("Error at line %d in %s: expecting 2 records, but got "OFFFMT" record(s)\n",
        __LINE__,__FILE__,length);
        nerrs++;
    }

    if (nerrs == 0) { /* test independent data mode */
        /* writing new records, HDF5 requires collective I/O */
        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            err = ncmpi_begin_indep_data(ncid); CHECK_ERR
        }

        /* write the 4th record */
        buf[0] = 93;
        start[0] = 3; count[0] = 1;
        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            err = ncmpi_put_vara_int(ncid, varid[0], start, count, buf); CHECK_ERR
        } else {
            err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR
        }

        /* write the 3rd and 4th records */
        buf[0] = 92; buf[1] = 93;
        start[0] = 2; count[0] = 2;
        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            err = ncmpi_put_vara_int(ncid, varid[0], start, count, buf); CHECK_ERR
        } else {
            err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR
        }

        err = ncmpi_inq_dimlen(ncid, dimid[0], &length); CHECK_ERR
        if (length != 4) {
            printf("Error at line %d in %s: expecting 4 records, but got "OFFFMT" record(s)\n",
                   __LINE__,__FILE__,length);
            nerrs++;
        }
        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            err = ncmpi_end_indep_data(ncid); CHECK_ERR
        }
    }

    /* REC_VAR_3D: write the 2nd record first */
    start[1] = 0; start[2] = 0; count[0] = 1; count[1] = 2; count[2] = 10;

    for (i=0; i<20; i++) buf[i] = 91;
    start[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR

    /* write the 1st record now */
    for (i=0; i<20; i++) buf[i] = 90;
    start[0] = 0;
    err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, dimid[0], &length); CHECK_ERR
    if (length != 4) {
        printf("Error at line %d in %s: expecting 4 records, but got "OFFFMT" record(s)\n",
        __LINE__,__FILE__,length);
        nerrs++;
    }

    if (nerrs == 0) { /* test independent data mode */
        /* writing new records, HDF5 requires collective I/O */
        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            err = ncmpi_begin_indep_data(ncid); CHECK_ERR
        }

        /* write the 4th record */
        for (i=0; i<20; i++) buf[i] = 93;
        start[0] = 3;
        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            err = ncmpi_put_vara_int(ncid, varid[1], start, count, buf); CHECK_ERR
        } else {
            err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR
        }

        /* write the 3rd record */
        for (i=0; i<20; i++) buf[i] = 92;
        start[0] = 2;
        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            err = ncmpi_put_vara_int(ncid, varid[1], start, count, buf); CHECK_ERR
        } else {
            err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR
        }

        err = ncmpi_inq_dimlen(ncid, dimid[0], &length); CHECK_ERR
        if (length != 4) {
            printf("Error at line %d in %s: expecting 4 records, but got "OFFFMT" record(s)\n",
                   __LINE__,__FILE__,length);
            nerrs++;
        }
        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            err = ncmpi_end_indep_data(ncid); CHECK_ERR
        }
    }
    err = ncmpi_close(ncid); CHECK_ERR
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
    int err=0, rank, flag;

    /* check whether burst buffering is enabled */
    MPI_Info_get(info, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, val, &flag);
    if (flag && strcasecmp(val, "enable") == 0 &&
        (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC))
        /* does not work for NetCDF4 files when burst-buffering is enabled */
        return 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank >= 1) return 0; /* this test is for running 1 process */

    err = test_only_record_var_1D(out_path, format, info);
    if (err > 0) return err;
    err = test_only_record_var_3D(out_path, format, info);
    if (err > 0) return err;
    err = test_two_record_var(out_path, format, info);
    if (err > 0) return err;

    return err;
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

    err = tst_main(argc, argv, "write records in reversed order", opt, test_io);

    MPI_Finalize();

    return err;
}
