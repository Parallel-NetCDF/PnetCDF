/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program adds two new variables to an existing netCDF file.
 * It is used to test if PnetCDF can correctly calculate the file offsets
 * for the two new variables, in particular for files that align the
 * fix-size variables to a boundary larger than 4 bytes, for instance
 * a file created by PnetCDF with defaut alignment of 512 bytes.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o add_var add_var.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 add_var testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

int main(int argc, char** argv) {
    char filename[256], var_name[NC_MAX_NAME];
    int i, nvars, rank, nprocs, err, nerrs=0;
    int ncid, varid, dimid[2];
    MPI_Offset prev_off, off;

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
        sprintf(cmd_str, "*** TESTING C   %s for checking offsets of new variables ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "dim_1", 5, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim_2", 4, &dimid[1]); CHECK_ERR

    /* define a bunch of variables */
    for (i=0; i<10; i++) {
        sprintf(var_name, "var_%d", i);
        err = ncmpi_def_var(ncid, var_name, NC_INT, 2, dimid, &varid); CHECK_ERR
    }
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* re-enter define mode */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* add 2 new dimensions */
    err = ncmpi_def_dim(ncid, "new_dim_1", 5, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "new_dim_2", 4, &dimid[1]); CHECK_ERR

    /* add 2 new variables */
    err = ncmpi_def_var(ncid, "new_var1", NC_INT,   2, dimid, &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "new_var2", NC_FLOAT, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_nvars(ncid, &nvars); CHECK_ERR
    err = ncmpi_inq_varoffset(ncid, 0, &prev_off); CHECK_ERR
    for (i=1; i<nvars; i++) {
        err = ncmpi_inq_varoffset(ncid, i, &off); CHECK_ERR
        if (off < prev_off + 5*4*4) { /* each variable is of size 5*4*4 bytes */
            err = ncmpi_inq_varname(ncid, i, var_name); CHECK_ERR
            printf("Error at line %d in %s: variable %s offset is set incorrectly\n",
                   __LINE__,__FILE__,var_name);
            nerrs++;
        }
        prev_off = off;
    }

    err = ncmpi_close(ncid); CHECK_ERR

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

    MPI_Finalize();
    return (nerrs > 0);
}

