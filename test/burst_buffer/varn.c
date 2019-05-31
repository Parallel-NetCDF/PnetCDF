/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests nonblocking functionality of bb driver
 * Flushed requests can not be canceled, canceled request can't be waited
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <limits.h>
#include <testutils.h>
#include <libgen.h>

int main(int argc, char *argv[]) {
    int i;
    int err, nerrs = 0;
    int rank, np;
    int ncid, varid;
    int dimid[2];
    int buffer[10];
    char *filename;
    MPI_Offset starts[10][2], counts[10][2];
    MPI_Offset *Starts[10], *Counts[10];
    MPI_Info info;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    /* Determine test file name */
    if (argc > 1)
        filename = argv[1];
    else
        filename = "testfile.nc";

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for when requests are > buffer size", basename(argv[0]));
                printf("%-66s ------ ", cmd_str); fflush(stdout);
                free(cmd_str);
    }

     /* Initialize file info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_burst_buf", "enable");

    /* Create new netcdf file */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);    CHECK_ERR

    /* Define dimensions */
    err = ncmpi_def_dim(ncid, "X", np, dimid);    CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", 10, dimid + 1);    CHECK_ERR

    /* Define variable */
    err = ncmpi_def_var(ncid, "M", NC_INT, 2, dimid, &varid);    CHECK_ERR

    /* Switch to data mode */
    err = ncmpi_enddef(ncid);    CHECK_ERR

    for(i = 0; i < 10; i++){
        starts[i][0] = rank;
        starts[i][1] = i;
        counts[i][0] = 1;
        counts[i][1] = 1;
        Starts[i] = (MPI_Offset*)starts[i];
        Counts[i] = (MPI_Offset*)counts[i];
        buffer[i] = rank + i;
    }

    /* Standard varn */
    err = ncmpi_put_varn_int_all(ncid, varid, 10, Starts, Counts, buffer);    CHECK_ERR
    err = ncmpi_get_varn_int_all(ncid, varid, 10, Starts, Counts, buffer);    CHECK_ERR
    for(i = 0; i < 10; i++){
        if (buffer[i] != rank + i){
            nerrs++;
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n", __LINE__, __FILE__, i, rank + 1, buffer[i]);
        }
    }

    /* NULL counts */
    err = ncmpi_put_varn_int_all(ncid, varid, 10, Starts, NULL, buffer);    CHECK_ERR
    err = ncmpi_get_varn_int_all(ncid, varid, 10, Starts, NULL, buffer);    CHECK_ERR
    for(i = 0; i < 10; i++){
        if (buffer[i] != rank + i){
            nerrs++;
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n", __LINE__, __FILE__, i, rank + 1, buffer[i]);
        }
    }

    /* Partial NULL counts */
    for(i = 0; i < 10; i += 2){
        Counts[i] = (MPI_Offset*)counts[i];
    }
    err = ncmpi_put_varn_int_all(ncid, varid, 10, Starts, Counts, buffer);    CHECK_ERR
    err = ncmpi_get_varn_int_all(ncid, varid, 10, Starts, Counts, buffer);    CHECK_ERR
    for(i = 0; i < 10; i++){
        if (buffer[i] != rank + i){
            nerrs++;
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n", __LINE__, __FILE__, i, rank + 1, buffer[i]);
        }
    }

    /* Close the file */
    err = ncmpi_close(ncid);    CHECK_ERR

    MPI_Info_free(&info);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n", sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR, nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();

    return nerrs > 0;
}
