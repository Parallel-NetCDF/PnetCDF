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
int test_only_record_var_1D(char *filename, int cmode)
{
    int ncid, varid, dimid, buf[20], err, nerrs=0;
    MPI_Offset start[1], count[1], length;
    MPI_Info info=MPI_INFO_NULL;

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_SELF, filename, cmode, info, &ncid); CHECK_ERR

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
        printf("Error at line %d in %s: expecting 2 records, but got %lld record(s)\n",
        __LINE__,__FILE__,length);
        nerrs++;
    }

    if (nerrs == 0 && !(cmode & NC_NETCDF4)) { /* test independent data mode */
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
            printf("Error at line %d in %s: expecting 4 records, but got %lld record(s)\n",
                   __LINE__,__FILE__,length);
            nerrs++;
        }
        err = ncmpi_end_indep_data(ncid); CHECK_ERR
    }
    err = ncmpi_close(ncid); CHECK_ERR
    return nerrs;
}

static
int test_only_record_var_3D(char *filename, int cmode)
{
    int i, ncid, varid, dimid[3], buf[20], err, nerrs=0;
    MPI_Offset start[3], count[3], length;
    MPI_Info info=MPI_INFO_NULL;

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_SELF, filename, cmode, info, &ncid); CHECK_ERR

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
        printf("Error at line %d in %s: expecting 2 records, but got %lld record(s)\n",
        __LINE__,__FILE__,length);
        nerrs++;
    }

    if (nerrs == 0 && !(cmode & NC_NETCDF4)) { /* test independent data mode */
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
            printf("Error at line %d in %s: expecting 4 records, but got %lld record(s)\n",
                   __LINE__,__FILE__,length);
            nerrs++;
        }
        err = ncmpi_end_indep_data(ncid); CHECK_ERR
    }
    err = ncmpi_close(ncid); CHECK_ERR
    return nerrs;
}

static
int test_two_record_var(char *filename, int cmode)
{
    int i, ncid, varid[2], dimid[3], buf[20], err, nerrs=0;
    MPI_Offset start[3], count[3], length;
    MPI_Info info=MPI_INFO_NULL;

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_SELF, filename, cmode, info, &ncid); CHECK_ERR

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
        printf("Error at line %d in %s: expecting 2 records, but got %lld record(s)\n",
        __LINE__,__FILE__,length);
        nerrs++;
    }

    if (nerrs == 0) { /* test independent data mode */
        /* writing new records, HDF5 requires collective I/O */
        if (!(cmode & NC_NETCDF4)) {
            err = ncmpi_begin_indep_data(ncid); CHECK_ERR
        }

        /* write the 4th record */
        buf[0] = 93;
        start[0] = 3; count[0] = 1;
        if (!(cmode & NC_NETCDF4)) {
            err = ncmpi_put_vara_int(ncid, varid[0], start, count, buf); CHECK_ERR
        } else {
            err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR
        }

        /* write the 3rd and 4th records */
        buf[0] = 92; buf[1] = 93;
        start[0] = 2; count[0] = 2;
        if (!(cmode & NC_NETCDF4)) {
            err = ncmpi_put_vara_int(ncid, varid[0], start, count, buf); CHECK_ERR
        } else {
            err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR
        }

        err = ncmpi_inq_dimlen(ncid, dimid[0], &length); CHECK_ERR
        if (length != 4) {
            printf("Error at line %d in %s: expecting 4 records, but got %lld record(s)\n",
                   __LINE__,__FILE__,length);
            nerrs++;
        }
        if (!(cmode & NC_NETCDF4)) {
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
        printf("Error at line %d in %s: expecting 4 records, but got %lld record(s)\n",
        __LINE__,__FILE__,length);
        nerrs++;
    }

    if (nerrs == 0) { /* test independent data mode */
        /* writing new records, HDF5 requires collective I/O */
        if (!(cmode & NC_NETCDF4)) {
            err = ncmpi_begin_indep_data(ncid); CHECK_ERR
        }

        /* write the 4th record */
        for (i=0; i<20; i++) buf[i] = 93;
        start[0] = 3;
        if (!(cmode & NC_NETCDF4)) {
            err = ncmpi_put_vara_int(ncid, varid[1], start, count, buf); CHECK_ERR
        } else {
            err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR
        }

        /* write the 3rd record */
        for (i=0; i<20; i++) buf[i] = 92;
        start[0] = 2;
        if (!(cmode & NC_NETCDF4)) {
            err = ncmpi_put_vara_int(ncid, varid[1], start, count, buf); CHECK_ERR
        } else {
            err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR
        }

        err = ncmpi_inq_dimlen(ncid, dimid[0], &length); CHECK_ERR
        if (length != 4) {
            printf("Error at line %d in %s: expecting 4 records, but got %lld record(s)\n",
                   __LINE__,__FILE__,length);
            nerrs++;
        }
        if (!(cmode & NC_NETCDF4)) {
            err = ncmpi_end_indep_data(ncid); CHECK_ERR
        }
    }
    err = ncmpi_close(ncid); CHECK_ERR
    return nerrs;
}

int main(int argc, char** argv) {
    char filename[256], *hint_value;
    int nerrs=0, rank, nprocs, err, bb_enabled=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        goto fn_exit;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for write records in reversed order", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    if (rank >= 1) goto fn_exit; /* this test is for running 1 process */

    /* check whether burst buffering is enabled */
    if (inq_env_hint("nc_burst_buf", &hint_value)) {
        if (strcasecmp(hint_value, "enable") == 0) bb_enabled = 1;
        free(hint_value);
    }

    /* CDF-1: test only one 1D record variable */
    nerrs += test_only_record_var_1D(filename, 0);
    /* CDF-1: test only one 3D record variable */
    nerrs += test_only_record_var_3D(filename, 0);
    /* CDF-1: test two record variables */
    nerrs += test_two_record_var(filename, 0);

    /* CDF-2: test only one 1D record variable */
    nerrs += test_only_record_var_1D(filename, NC_64BIT_OFFSET);
    /* CDF-2: test only one 3D record variable */
    nerrs += test_only_record_var_3D(filename, NC_64BIT_OFFSET);
    /* CDF-2: test two record variables */
    nerrs += test_two_record_var(filename, NC_64BIT_OFFSET);

    if (!bb_enabled) {
#ifdef USE_NETCDF4
        /* NETCDF4: test only one 1D record variable */
        nerrs += test_only_record_var_1D(filename, NC_NETCDF4);
        /* NETCDF4: test only one 3D record variable */
        nerrs += test_only_record_var_3D(filename, NC_NETCDF4);
        /* NETCDF4: test two record variables */
        nerrs += test_two_record_var(filename, NC_NETCDF4);

        /* NETCDF4_CLASSIC: test only one 1D record variable */
        nerrs += test_only_record_var_1D(filename, NC_NETCDF4|NC_CLASSIC_MODEL);
        /* NETCDF4_CLASSIC: test only one 3D record variable */
        nerrs += test_only_record_var_3D(filename, NC_NETCDF4|NC_CLASSIC_MODEL);
        /* NETCDF4_CLASSIC: test two record variables */
        nerrs += test_two_record_var(filename, NC_NETCDF4|NC_CLASSIC_MODEL);
#endif
    }

    /* CDF-5: test only one 1D record variable */
    nerrs += test_only_record_var_1D(filename, NC_64BIT_DATA);
    /* CDF-5: test only one 3D record variable */
    nerrs += test_only_record_var_3D(filename, NC_64BIT_DATA);
    /* CDF-5: test two record variables */
    nerrs += test_two_record_var(filename, NC_64BIT_DATA);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR && malloc_size > 0) { /* this test is for running 1 process */
        printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
               malloc_size);
        ncmpi_inq_malloc_list();
    }

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

