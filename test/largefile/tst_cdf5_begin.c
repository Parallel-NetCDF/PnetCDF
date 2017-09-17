/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program creates one big variable followed by a small variable. It
 * writes some data to both variables and reads back to check the values.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_cdf5_begin tst_cdf5_begin.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_cdf5_begin
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

/* When using NetCDF 4.4.1 ad prior to create a CDF-5 file and defining a small
 * variable after a big variable (> 2^32-3 bytes), the file starting offset of
 * the small variable (and all variables defined after the big variable) is
 * calculated incorrectly. This test program detects this bug by checking the
 * contents of the possible overlaps between the two variables.
 */

int main(int argc, char** argv) {
    char filename[256];
    int i, err, rank, nprocs, nerrs=0, ncid, dimid[2], varid[2];
    short buf[10];
    MPI_Offset start[1], count[1];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        goto fn_exit;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for checking CDF-5 writes", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_DATA,
                       MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim0", NC_MAX_UINT, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", 10,          &dimid[1]); CHECK_ERR

    /* define one small variable after one big variable */
    err = ncmpi_def_var(ncid, "var_big",   NC_SHORT, 1, &dimid[0], &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_small", NC_SHORT, 1, &dimid[1], &varid[1]); CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_NOFILL, NULL); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* write to var_big in location overlapping with var_small when using
     * netCDF 4.4.x or prior */
    start[0] = NC_MAX_UINT/sizeof(short);
    count[0] = 10;
    for (i=0; i<10; i++) buf[i] = i;
    err = ncmpi_put_vara_short_all(ncid, varid[0], start, count, buf); CHECK_ERR

    /* write var_small */
    for (i=0; i<10; i++) buf[i] = -1;
    err = ncmpi_put_var_short_all(ncid, varid[1], buf); CHECK_ERR

    /* read back var_big and check contents */
    for (i=0; i<10; i++) buf[i] = -1;
    err = ncmpi_get_vara_short_all(ncid, varid[0], start, count,buf); CHECK_ERR
    for (i=0; i<10; i++) {
        if (buf[i] != i) {
            printf("Error at buf[%d] expect %d but got %hd\n",i,i,buf[i]);
            nerrs++;
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

