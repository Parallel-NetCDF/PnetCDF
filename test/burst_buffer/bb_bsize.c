/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests writing a variable exceeding the size of data buffer.
 * Each process writes a submatrix of size 1024 * 1024, a total of 1M cells
 * The submatrix form each processes is stacked along the first dimension so
 * that the variable dimensions are (1024 * np) * 1024 Each processes writes
 * it's rank to every cell
 * Each process starts by writing the first 1/8 of the rows at once, which
 * should caused an increase of the data buffer size to accommodate it Then each
 * process writes the remaining row one at a time
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <limits.h>
#include <testutils.h>
#include <libgen.h>

#define SIZE 1024

int buffer[SIZE * SIZE];
char bsize[32];

int main(int argc, char *argv[]) {
    int i, ret = NC_NOERR, nerr = 0;
    int rank, np;
    int ncid, varid;
    int dimid[2];
    char *filename;
    MPI_Offset start[2], count[2];
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
        sprintf(cmd_str, "*** TESTING C   %s for checking request > buffer size", basename(argv[0]));
                printf("%-66s ------ ", cmd_str); fflush(stdout);
                free(cmd_str);
    }

    /* Initialize file info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_burst_buf", "enable");
    /* Set default buffer size to 1/16 of the rows */
    sprintf(bsize, "%u", (unsigned int)(SIZE * SIZE / 16 * sizeof(int)));
    MPI_Info_set(info, "nc_burst_buf_flush_buffer_size", bsize);

    /* Create new netcdf file */
    ret = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_create: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Define dimensions */
    ret = ncmpi_def_dim(ncid, "X", SIZE * np, dimid);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_def_dim: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }
    ret = ncmpi_def_dim(ncid, "Y", SIZE, dimid + 1);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_def_dim: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Define variable */
    ret = ncmpi_def_var(ncid, "M", NC_INT, 2, dimid, &varid);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_def_var: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Switch to data mode */
    ret = ncmpi_enddef(ncid);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_enddef: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Initialize buffer */
    for (i = 0; i < SIZE * SIZE; i++) {
        buffer[i] = rank + 1;
    }

    /* Write first 1/8 of the rows */
    start[0] = SIZE * rank;
    start[1] = 0;
    count[0] = SIZE / 8;
    count[1] = SIZE;
    ret = ncmpi_put_vara_int_all(ncid, varid, start, count, buffer);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_put_vara_int: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Write remaining rows */
    start[0] = SIZE * rank + SIZE / 8;
    start[1] = 0;
    count[0] = 1;
    count[1] = SIZE;
    for (; start[0] < SIZE * (rank + 1); start[0]++) {
        ret = ncmpi_put_vara_int_all(ncid, varid, start, count, buffer);
        if (ret != NC_NOERR) {
            printf("Error at line %d in %s: ncmpi_put_vara_int: %d\n", __LINE__, __FILE__, ret);
            nerr++;
            goto ERROR;
        }
    }

    /*
     * Read it back
     * Flush on read is on by default so no additional action required
     */
    memset(buffer, 0, sizeof(buffer));
    start[0] = SIZE * rank;
    start[1] = 0;
    count[0] = SIZE;
    count[1] = SIZE;
    ret = ncmpi_get_vara_int_all(ncid, varid, start, count, buffer);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_get_vara_int: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Verify the result */
    for (i = 0; i < SIZE * SIZE; i++) {
        if (buffer[i] != rank + 1) {
            nerr++;
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n", __LINE__, __FILE__, i, rank + 1, buffer[i]);
        }
    }

    /* Close the file */
    ret = ncmpi_close(ncid);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_close: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    MPI_Info_free(&info);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    ret = ncmpi_inq_malloc_size(&malloc_size);
    if (ret == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n", sum_size);
    }

ERROR:
    MPI_Allreduce(MPI_IN_PLACE, &nerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerr) printf(FAIL_STR, nerr);
        else       printf(PASS_STR);
    }

    MPI_Finalize();

    return nerr > 0;
}
