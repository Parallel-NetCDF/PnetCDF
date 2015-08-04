/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if one can get the number of record variables and fixed-
 * sized variables correctly. It first defines some number of fixed-size and
 * record variables and then calls the APIs
 *     ncmpi_inq_num_rec_vars() and ncmpi_inq_num_fix_vars()
 * to varify if the numbers are correct.
 *
 * The compile and run commands are given below. This program is to be run on
 * one MPI process.
 *
 *    % mpicc -g -o inq_num_vars inq_num_vars.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 inq_num_vars testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define FAIL_COLOR "\x1b[31mfail\x1b[0m\n"
#define PASS_COLOR "\x1b[32mpass\x1b[0m\n"

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}

void check_num_vars(int  ncid,
                    int  expected_nvars,
                    int  expected_num_rec_vars,
                    int  expected_num_fix_vars,
                    int *nfailed)
{
    int err, nvars, num_rec_vars, num_fix_vars;

    err = ncmpi_inq_nvars(ncid, &nvars); ERR
    err = ncmpi_inq_num_rec_vars(ncid, &num_rec_vars); ERR
    err = ncmpi_inq_num_fix_vars(ncid, &num_fix_vars); ERR

    if (nvars != expected_nvars) {
        printf("Error: expecting %d number of variables defined, but got %d\n", expected_nvars, nvars);
        (*nfailed)++;
    }
    if (num_rec_vars != expected_num_rec_vars) {
        printf("Error: expecting %d number of record variables defined, but got %d\n", expected_num_rec_vars, num_rec_vars);
        (*nfailed)++;
    }
    if (num_fix_vars != expected_num_fix_vars) {
        printf("Error: expecting %d number of fixed-size variables defined, but got %d\n", expected_num_fix_vars, num_fix_vars);
        (*nfailed)++;
    }
}

int main(int argc, char** argv) {
    char *filename="testfile.nc";
    int nfailed, nfailed_all, rank, nprocs, err;
    int ncid, cmode, varid[7], dimid[3];
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        goto fn_exit;
    }
    if (argc == 2) filename = argv[1];

    /* printf("PnetCDF version string: \"%s\"\n", ncmpi_inq_libvers()); */

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid); ERR

    /* define dimension and variable */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "Y",       2,            &dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "X",       10,           &dimid[2]); ERR

    nfailed = 0;

    err = ncmpi_def_var(ncid, "REC_VAR_1", NC_INT, 1, dimid, &varid[0]); ERR
    err = ncmpi_def_var(ncid, "REC_VAR_2", NC_INT, 3, dimid, &varid[1]); ERR
    err = ncmpi_def_var(ncid, "REC_VAR_3", NC_INT, 2, dimid, &varid[2]); ERR
    err = ncmpi_def_var(ncid, "REC_VAR_4", NC_INT, 1, dimid, &varid[3]); ERR

    check_num_vars(ncid, 4, 4, 0, &nfailed);

    err = ncmpi_def_var(ncid, "FIX_VAR_1", NC_INT, 2, dimid+1, &varid[4]); ERR
    err = ncmpi_def_var(ncid, "FIX_VAR_2", NC_INT, 1, dimid+1, &varid[5]); ERR
    err = ncmpi_def_var(ncid, "FIX_VAR_3", NC_INT, 1, dimid+2, &varid[6]); ERR

    check_num_vars(ncid, 7, 4, 3, &nfailed);

    err = ncmpi_enddef(ncid); ERR

    check_num_vars(ncid, 7, 4, 3, &nfailed);

    err = ncmpi_close(ncid); ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Reduce(&nfailed, &nfailed_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        char cmd_str[256];
        sprintf(cmd_str, "*** TESTING C   %s for no. record/fixed variables", argv[0]);
        if (nfailed_all > 0)
            printf("%s ------ "FAIL_COLOR" with %d mismatches\n",cmd_str,nfailed_all);
        else
            printf("%-66s ------ " PASS_COLOR, cmd_str);
    }

fn_exit:
    MPI_Finalize();
    return 0;
}

