/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if NC_ENOTNC can be returned when opening a CDF-5 file
 * with corrupted header. In the corrupted file, bad_begin.nc5, the file
 * starting offset of second variable "var_small" is incorrectly set to
 * a value smaller than the end offset of its previous variable "var_big".
 * Command "ncoffsets bad_begin.nc5" shows the followings.
 *     ncoffsets bad_begin.nc5 
 *     netcdf bad_begin.nc5 {
 *     // file format: CDF-5
 *
 *     file header:
 *             size   = 220 bytes
 *             extent = 220 bytes
 *
 *     dimensions:
 *             dim0 = 4294967295
 *             dim1 = 10
 *
 *     fixed-size variables:
 *             short  var_big(dim0):
 *                    start file offset =         220
 *                    end   file offset =  8589934810
 *             short  var_small(dim1):
 *                    start file offset =  4294967515
 *                    end   file offset =  4294967535
 *     }
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_open_cdf5 tst_open_cdf5.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_open_cdf5
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h> /* setenv() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define FILE_NAME "simple_xy.nc4"

int main(int argc, char** argv) {
    char filename[256];
    int i, j, k, nerrs=0, rank, nprocs, err, ncid;
    int ndim, len, nvar;
    MPI_Offset dlen[1024];
    char tmp[1024];
    int dims[1024];
    int data[1024];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        goto fn_exit;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, FILE_NAME);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str,
        "*** TESTING C   %s for checking netcdf4 open",
        basename(argv[0]));
        printf("%-66s --- ", cmd_str); fflush(stdout);
        free(cmd_str);
    }


    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_ndims(ncid, &ndim); CHECK_ERR
    printf("number of dims: %d\n", ndim);

    for(i = 0; i < ndim; i++){
        err = ncmpi_inq_dim(ncid, i, tmp, dlen + i); CHECK_ERR
        printf("dim %s: %llu\n", tmp, dlen[i]);
    }

    err = ncmpi_inq_nvars(ncid, &nvar); CHECK_ERR
    for(i = 0; i < nvar; i++){
        err = ncmpi_inq_var(ncid, i, tmp, NULL, &ndim, NULL, NULL); CHECK_ERR
        err = ncmpi_inq_vardimid(ncid, i, dims); CHECK_ERR
        printf("var %s: (", tmp);
        for(j = 0; j < ndim; j++){
            err = ncmpi_inq_dimname (ncid, dims[j], tmp); CHECK_ERR
            printf("%s, ", tmp);
        }
        printf(")\n data: \n");
        err = ncmpi_get_var_int_all(ncid, i, data);
        if (ndim == 2){
            for(j = 0; j < dlen[dims[0]]; j++){
                for(k = 0; k < dlen[dims[1]]; k++){
                    printf("%d, ", data[j * dlen[dims[1]] + k]);
                }
                printf("\n"); 
            }
        }
    }
    
    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR && malloc_size > 0) /* this test is for running 1 process */
        printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
               malloc_size);

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

