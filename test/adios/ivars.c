/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program verify variable read capability of adios driver
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h> /* setenv() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include "pnetcdf.h"
#include <math.h>

#include <testutils.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "arrays.bp"
#define V1_NAME "var_double_2Darray"
#define V2_NAME "var_int_1Darray"
#define D1_NAME "NX"
#define D2_NAME "NY"

/* We are reading 2D data, a 6 x 12 grid. */
#define NX 10ULL
#define NY 100ULL
#define NDIMS 2

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main(int argc, char** argv) {
    char filename[256];
    int i, j, nerrs=0, rank, nprocs, err;
    int ncid;
    MPI_Offset start[2], count[2], stride[2];
    double data[2][2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        nerrs++;
        goto fn_exit;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, FILE_NAME);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str,
        "*** TESTING C   %s for using iget_vars to read a bp file",
        basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    start[0] = 0;
    start[1] = 0;
    count[0] = 2;
    count[1] = 2;
    stride[0] = 5;
    stride[1] = 50;
    err = ncmpi_iget_vars_double(ncid, 0, start, count, stride, (double*)data,
                                    NULL); CHECK_ERR

    err = ncmpi_wait_all(ncid, NC_GET_REQ_ALL, NULL, NULL); CHECK_ERR

    for(i = 0; i < 2; i++){
        for(j = 0; j < 2; j++){
            if (fabs(data[i][j] - (((double)(i * stride[0])) + ((double)(j *
                stride[1])) / 100)) > 0.0001){
                printf("Rank %d: Expect Var 0 [%d][%d] = %lf, but got %lf\n",
                        rank, i * (int)stride[0], j * (int)stride[1],
                        (((double)(i * stride[0])) + ((double)(j * stride[1]))
                         / 100), data[i][j]);
                nerrs++;
            }
        }
    }

    ncmpi_close(ncid);

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

