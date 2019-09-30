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

static int
tst_fmt(char *filename, int cmode)
{
    int nerrs=0, err, ncid, varid[7], dimid[3];
    MPI_Info info=MPI_INFO_NULL;

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid); CHECK_ERR

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

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char** argv) {
    char filename[256], *hint_value;
    int nerrs=0, rank, err, bb_enabled=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        goto fn_exit;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for no. record/fixed variables", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* printf("PnetCDF version string: \"%s\"\n", ncmpi_inq_libvers()); */

    /* check whether burst buffering is enabled */
    if (inq_env_hint("nc_burst_buf", &hint_value)) {
        if (strcasecmp(hint_value, "enable") == 0) bb_enabled = 1;
        free(hint_value);
    }

    nerrs += tst_fmt(filename, 0);
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += tst_fmt(filename, NC_NETCDF4);
        nerrs += tst_fmt(filename, NC_NETCDF4 | NC_CLASSIC_MODEL);
#endif
    }
    nerrs += tst_fmt(filename, NC_64BIT_DATA);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

fn_exit:
    MPI_Finalize();
    return (nerrs > 0);
}

