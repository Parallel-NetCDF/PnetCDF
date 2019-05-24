/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This example demonstrates how to read variable data in BP file using PnetCDF
 * 
 * We read the file located at test/adios/array.bp
 * 
 * The content of the file is:
 * netcdf arrays {
 * // file format: ADIOS BP Ver. 3
 *  dimensions:
 *          NX = 10 ;
 *          NY = 100 ;
 *  variables:
 *          double var_double_2Darray(NX, NY) ;
 *          int var_int_1Darray(NX) ;
 *  }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h> /* setenv() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include "pnetcdf.h"
#include <math.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "../../test/adios/arrays.bp"

/* We are reading 2D data, a 10 x 100 grid. */
#define NX 10ULL
#define NY 100ULL
#define NDIMS 2

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main(int argc, char** argv) {
    int i, rank, nprocs;
    int ncid;
    int reqids[2];
    MPI_Offset start[2], count[2];
    double data_double[NY];
    int data_int[NY];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    ncmpi_open(MPI_COMM_WORLD, FILE_NAME, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    
    /* The content of variable var_double_2Darray is 
     * var_double_2Darray[x][y] =  x + y / 100 
     */

    /* Collective read */
    start[0] = rank % NX;
    start[1] = 0;
    count[0] = 1;
    count[1] = NY;
    ncmpi_iget_vara_double(ncid, 0, start, count, data_double, reqids); 

    /* Read with different datatype */
    ncmpi_iget_vara_int(ncid, 0, start, count, data_int, reqids + 1); 

    /* Wait for the requests */
    ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);


    /* Check reqults */
    for(i = 0; i < NY; i++){
        if (fabs(data_double[i] - (((double)start[0]) + ((double)i) / 100)) 
        > 0.0001){
            printf("Rank %d: Expect Var 0 [%llu][%d] = %lf, but got %lf\n", 
                    rank, start[0], i, ((double)start[0]) + ((double)i) / 100, 
                    data_double[i]);
        }
    }

    for(i = 0; i < NY; i++){
        if (data_int[i] != (int)start[0]) {
            printf("Rank %d: Expect Var 0 [%llu][%d] = %d, but got %d\n", rank, 
                    start[0], i, (int)start[0], data_int[i]);
        }
    }

    ncmpi_close(ncid);

    MPI_Finalize();
    return 0;
}

