/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example tests writing a high dimensional int variable when log io is
 * enabled.  Each process own a submatrix of size 2 * 2 * 2 * ... 1 * 1 * 1
 * that have as much size 2 dimensions as the buffer can accommodate.
 * All cell in the submatrix is it's rank + 1
 * The submatrix is combined by interleaving along the first dimension.
 * The variable will have size (2 * np) * 2 * 2 * ... 1 * 1 * 1
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file when dimension is set to 2.
 *
 *    % mpicc -O2 -o log_higndim.c -lpnetcdf
 *    % mpiexec -n 4 ./log_higndim [testfile.nc]
 *    % ncmpidump [testfile.nc]
 *    netcdf test {
 *    // file format: CDF-1
 *    dimensions:
 *            D0 = 8 ;
 *            D1 = 2 ;
 *    variables:
 *            int M(D0, D1) ;
 *    data:
 *
 *      M =
 *
 *       1, 1,
 *       2, 2,
 *       3, 3,
 *       4, 4,
 *       1, 1,
 *       2, 2,
 *       3, 3,
 *       4, 4 ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <limits.h>
#include <testutils.h>
#include <libgen.h>

#if NC_MAX_DIMS < 1024
#define DIM NC_MAX_DIMS
#else
#define DIM 1024
#endif

#define BSIZE 1024 * 1024

int main(int argc, char *argv[])
{
    char *filename=NULL, dimname[64];
    int i, ret=NC_NOERR, nerr=0;
    int rank, np, ndims;
    int ncid, varid;
    int *dimid=NULL, *buffer=NULL;
    long long j;
    MPI_Offset *start=NULL, *count=NULL, *stride=NULL;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (argc > 3) {
        if (!rank) printf("Usage: %s [filename]\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    /* Determine ndims and test file name */
    if (argc > 1)
        filename = strdup(argv[1]);
    else
        filename = strdup("testfile.nc");

    ndims = DIM;
    if (argc > 2)
        ndims = atoi(argv[2]);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for high dimensional variables", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        fflush(stdout);
        free(cmd_str);
    }

    /* Allocate buffers */
    dimid = (int*)malloc(sizeof(int) * ndims);
    start = (MPI_Offset*)malloc(sizeof(MPI_Offset) * ndims);
    count = (MPI_Offset*)malloc(sizeof(MPI_Offset) * ndims);
    stride = (MPI_Offset*)malloc(sizeof(MPI_Offset) * ndims);
    buffer = (int*)malloc(sizeof(int) * BSIZE);
    if (dimid == NULL || start == NULL || count == NULL || stride == NULL || buffer == NULL) {
        printf("Error at line %d in %s: malloc error\n", __LINE__, __FILE__);
        nerr++;
        goto ERROR;
    }

    /* Initialize buffers and calculate share among processes
     * Each process writes a hyper rectangle to the variable
     * The variable is formed by stacking the rectangle of every processes along first dimension
     */
    for (i = 0; i < BSIZE; i++) {
        buffer[i] = rank + 1;
    }
    for (i = 0, j = 2; i < ndims; i++) {
        start[i] = 0;
        stride[i] = 1;
        /* Most dimensions must be 1 for high dimensional variable
         * Set dimension size to 2 until we run out of buffer
         */
        if (j < BSIZE) {
            count[i] = 2;
            j <<= 1;
        }
        else{
            count[i] = 1;
        }
    }
    start[0] = rank;
    stride[0] = np;

    /* Create new netcdf file */
    ret = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_create: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Define dimensions */
    for (i = 0; i < ndims; i++) {
        sprintf(dimname, "D%d", i);
        /* Submatrix of each process stack along the first dimension */
        if (i == 0) {
            ret = ncmpi_def_dim(ncid, dimname, count[i] * np, dimid + i);
        }
        else{
            ret = ncmpi_def_dim(ncid, dimname, count[i], dimid + i);
        }
        if (ret != NC_NOERR) {
            printf("Error at line %d in %s: ncmpi_enddef: %d\n", __LINE__, __FILE__, ret);
            nerr++;
            goto ERROR;
        }
    }

    /* Define variable */
    ret = ncmpi_def_var(ncid, "M", NC_INT, ndims, dimid, &varid);
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

    /* Write variable */
    ret = ncmpi_put_vars_int_all(ncid, varid, start, count, stride, buffer);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_put_vars_int: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Read it back */
    ret = ncmpi_get_vars_int_all(ncid, varid, start, count, stride, buffer);
    if (ret != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_get_vars_int: %d\n", __LINE__, __FILE__, ret);
        nerr++;
        goto ERROR;
    }

    /* Verify the result */
    for (i = 0; i < BSIZE; i++) {
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

    if (start != NULL) free(start);
    if (count != NULL) free(count);
    if (stride != NULL) free(stride);
    if (dimid != NULL) free(dimid);
    if (buffer != NULL) free(buffer);
    if (filename != NULL) free(filename);

    MPI_Finalize();

    return nerr > 0;
}
