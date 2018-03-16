/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program check if error code NC_EMAXDIMS can be returned correctly, when
 * defining a variable with more than NC_MAX_VAR_DIMS dimensions.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_max_var_dims tst_max_var_dims.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_max_var_dims testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <limits.h> /* INT_MAX */
#include <pnetcdf.h>

#include <testutils.h>

int main(int argc, char** argv) {
    char filename[256];
    int rank, nprocs, nerrs=0;
    int err, ncid;
#if NC_MAX_VAR_DIMS < INT_MAX
    int i, varid, *dimid;
#endif

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for checking NC_MAX_VAR_DIMS ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

#if NC_MAX_VAR_DIMS < INT_MAX
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* define dimensions */
    dimid = (int*) malloc((NC_MAX_VAR_DIMS+2) * sizeof(int));
    err = ncmpi_def_dim(ncid, "dim0", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", 1, &dimid[1]); CHECK_ERR

    for (i=2; i<NC_MAX_VAR_DIMS+2; i++) dimid[i] = dimid[1];

    /* define variables */
    err = ncmpi_def_var(ncid, "v0", NC_INT, NC_MAX_VAR_DIMS+1, &dimid[0], &varid);
    EXP_ERR(NC_EMAXDIMS)

    err = ncmpi_def_var(ncid, "v1", NC_INT, NC_MAX_VAR_DIMS+1, &dimid[1], &varid);
    EXP_ERR(NC_EMAXDIMS)

    err = ncmpi_set_fill(ncid, NC_NOFILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    free(dimid);

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
#else
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    if (rank == 0) printf(SKIP_STR);
#endif

    MPI_Finalize();
    return (nerrs > 0);
}

