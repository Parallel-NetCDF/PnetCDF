/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
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
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <unistd.h> /* getopt() */

#include <mpi.h>
#include "pnetcdf.h"
#include <math.h>

static int verbose;

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

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] filename\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       filename - input BP file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv) {
    extern int optind;
    char filename[256];
    int i, rank, nprocs;
    int ncid;
    MPI_Offset start[2], count[2];
    double data_double[NY];
    int data_int[NY];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] != NULL)
        snprintf(filename, 256, "%s", argv[optind]);
    else {
        if (rank==0) {
            printf("Error: input file is required\n");
            usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);

    /* The content of variable var_double_2Darray is
     * var_double_2Darray[x][y] =  x + y / 100
     */

    /* Collective read */
    start[0] = rank % NX;
    start[1] = 0;
    count[0] = 1;
    count[1] = NY;
    ncmpi_get_vara_double_all(ncid, 0, start, count, data_double);
    for(i = 0; i < NY; i++){
        if (fabs(data_double[i] - (((double)start[0]) + ((double)i) / 100))
             > 0.0001){
            printf("Rank %d: Expect Var 0 [%llu][%d] = %lf, but got %lf\n",
                    rank, start[0], i, ((double)start[0]) + ((double)i) / 100,
                    data_double[i]);
        }
    }

    /* Read with different datatype */
    ncmpi_get_vara_int_all(ncid, 0, start, count, data_int);
    for(i = 0; i < NY; i++){
        if (data_int[i] != (int)start[0]) {
            printf("Rank %d: Expect Var 0 [%llu][%d] = %d, but got %d\n", rank,
                    start[0], i, (int)start[0], data_int[i]);
        }
    }

    /* Independent read */
    ncmpi_begin_indep_data(ncid);
    ncmpi_get_vara_double(ncid, 0, start, count, data_double);
    for(i = 0; i < NY; i++){
        if (fabs(data_double[i] - (((double)start[0]) + ((double)i) / 100))
            > 0.0001){
            printf("Rank %d: Expect Var 0 [%llu][%d] = %lf, but got %lf\n",
                    rank, start[0], i, ((double)start[0]) + ((double)i) / 100,
                    data_double[i]);
        }
    }

    ncmpi_close(ncid);

    MPI_Finalize();
    return 0;
}

