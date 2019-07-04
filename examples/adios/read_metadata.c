/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This examples demonstrates how to enumerate all variable, dimension,
 * and attributes in a BP file using PnetCDF.
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
#define FILE_NAME "../../test/adios/attributes.bp"

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

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
    int i, j, rank, nprocs, err, nerrs = 0;
    int ncid, natt, ndim, nudim, nvar;
    int dimids[1024];
    char name[NC_MAX_NAME + 1];
    nc_type vtype;
    MPI_Offset len;

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

    /* Open ADIOS BP file as if opening a netcdf file
     * PnetCDF can only read BP files for now, NC_NOWRITE must be set
     */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                        &ncid);
    ERR

    err = ncmpi_inq(ncid, &ndim, &nvar, &natt, &nudim);
    ERR
    if (verbose)
        printf("ndim: %d, nvar: %d, natt: %d, nudim: %d\n", ndim, nvar, natt,
               nudim);

    for(i = 0; i < ndim; i++){
        err = ncmpi_inq_dim(ncid, i, name, &len);
        if (verbose)
            printf("Dim %d: name = %s, length = %llu\n", i, name, len);
    }

    for(i = 0; i < natt; i++){
        err = ncmpi_inq_attname (ncid, NC_GLOBAL, i, name);
        ERR
        err = ncmpi_inq_att(ncid, NC_GLOBAL, name, NULL, &len);
        ERR
        if (verbose)
            printf("Att %d: name = %s, length = %llu\n", i, name, len);
    }

    for(i = 0; i < nvar; i++){
        err = ncmpi_inq_var(ncid, i, name, &vtype, &ndim, NULL, &natt);
        ERR
        if (verbose)
            printf("Var %d: name = %s, ndim = %d, natt = %d, dimmids = ", i,
                   name,  ndim, natt);
        if (ndim < 1024){
            err = ncmpi_inq_var(ncid, i, name, &vtype, &ndim, dimids, &natt);
            ERR
            if (verbose) {
                for (j=0; j<ndim; j++)
                    printf("%d, ", dimids[j]);
            }
        }
        if (verbose) printf("\n");
        for(j = 0; j < natt; j++){
            err = ncmpi_inq_attname (ncid, i, j, name);
            ERR
            err = ncmpi_inq_att(ncid, i, name, NULL, &len);
            ERR
            if (verbose)
                printf("\tAtt %d: name = %s, length = %llu\n", j, name, len);
        }
    }

    ncmpi_close(ncid);

    MPI_Finalize();

    return 0;
}

