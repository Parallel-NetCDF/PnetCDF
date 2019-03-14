/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This examples demonstrates how to enumerate all variable, dimension, 
 * and attributes in a BP file using PnetCDF.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h> /* setenv() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include "pnetcdf.h"
#include <math.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "../../test/adios/attributes.bp"

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

int main(int argc, char** argv) {
    int i, j, rank, nprocs, err, nerrs = 0;
    int ncid, natt, ndim, nudim, nvar;
    int dimids[1024];
    char name[NC_MAX_NAME + 1];
    nc_type vtype;
    MPI_Offset len;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Open ADIOS BP file as if opening a netcdf file
    // PnetCDF can only read BP files for now, NC_NOWRITE must be set
    err = ncmpi_open(MPI_COMM_WORLD, FILE_NAME, NC_NOWRITE, MPI_INFO_NULL, 
                        &ncid);
    ERR

    err = ncmpi_inq(ncid, &ndim, &nvar, &natt, &nudim);
    ERR
    printf("ndim: %d, nvar: %d, natt: %d, nudim: %d\n", ndim, nvar, natt, 
            nudim);

    for(i = 0; i < ndim; i++){
        err = ncmpi_inq_dim(ncid, i, name, &len);
        printf("Dim %d: name = %s, length = %llu\n", i, name, len);
    }

    for(i = 0; i < natt; i++){
        err = ncmpi_inq_attname (ncid, NC_GLOBAL, i, name);
        ERR
        err = ncmpi_inq_att(ncid, NC_GLOBAL, name, NULL, &len);
        ERR
        printf("Att %d: name = %s, length = %llu\n", i, name, len);
    }

    for(i = 0; i < nvar; i++){
        err = ncmpi_inq_var(ncid, i, name, &vtype, &ndim, NULL, &natt);
        ERR
        printf("Var %d: name = %s, ndim = %d, natt = %d, dimmids = ", i, name, 
                ndim, natt);
        if (ndim < 1024){
            err = ncmpi_inq_var(ncid, i, name, &vtype, &ndim, dimids, &natt);
            ERR
            for(j = 0; j < ndim; j++){
                printf("%d, ", dimids[j]);
            }
        }
        printf("\n");
        for(j = 0; j < natt; j++){
            err = ncmpi_inq_attname (ncid, i, j, name);
            ERR
            err = ncmpi_inq_att(ncid, i, name, NULL, &len);
            ERR
            printf("\tAtt %d: name = %s, length = %llu\n", j, name, len);
        }
    }

    ncmpi_close(ncid);

    MPI_Finalize();

    return 0;
}

