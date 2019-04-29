/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests vars APIs for variables with high dimensions.
 * In particular, this is to test fix in r3164.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

#define NVARS 2
#define NDIMS 16
#define DIMLEN 3
#define NRECS 4

int main(int argc, char** argv) {
    char filename[256], name[32];
    size_t nelms;
    short *buffer;
    int i, j, cmode, rank, nprocs, err, nerrs=0;
    int ncid, fvarid[NVARS], rvarid[NVARS], dimids[NDIMS], rdimids[NDIMS];
    MPI_Offset start[NDIMS], count[NDIMS], stride[NDIMS];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    memset(filename, 0, 256);
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for vars APIs on high-dim variables ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    cmode = NC_CLOBBER;
    cmode |= NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL,
                       &ncid); CHECK_ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "rdim", NC_UNLIMITED, &rdimids[0]); CHECK_ERR
    for (i=0; i<NDIMS; i++) {
        sprintf(name, "dim%d", i);
        err = ncmpi_def_dim(ncid, name, DIMLEN, &dimids[i]); CHECK_ERR
        if (i > 0) rdimids[i] = dimids[i];
    }
    for (i=0; i<NVARS; i++) {
        sprintf(name, "fix_var%d", i);
        err = ncmpi_def_var(ncid, name, NC_SHORT, NDIMS, dimids, &fvarid[i]); CHECK_ERR
        sprintf(name, "rec_var%d", i);
        err = ncmpi_def_var(ncid, name, NC_SHORT, NDIMS,rdimids, &rvarid[i]); CHECK_ERR
    }

    /* initialize the contents of fixed-size variables */
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    if (err != NC_NOERR) goto fn_exit;

    nelms = (NRECS > DIMLEN) ? NRECS : DIMLEN;
    for (i=1; i<NDIMS; i++) nelms *= DIMLEN;
    buffer = (short*) malloc(nelms * sizeof(short));
    if (buffer == NULL) {
        printf("Error %s at line %d: fail to allocate buffer of size %zu\n",
               argv[0], __LINE__, nelms * sizeof(int));
        goto fn_exit;
    }
    for (i=0; i<nelms; i++) buffer[i] = -1;


    /* initialize the contents of record variables */
    for (i=0; i<NVARS; i++) {
        for (j=0; j<NRECS; j++) {
            err = ncmpi_fill_var_rec(ncid, rvarid[i], j); CHECK_ERR
        }
    }

    for (i=0; i<nelms; i++) buffer[i] = i % 32768;

    for (i=0; i<NDIMS; i++) {
        start[i]  = 0;
        count[i]  = 2;
        stride[i] = 2;
    }
    /* only process 0 writes */
    if (rank > 0) for (i=0; i<NDIMS; i++) count[i] = 0;

    for (i=0; i<NVARS; i++) {
        start[0] = 0;
        err = ncmpi_put_vars_short_all(ncid, fvarid[i], start, count, stride,
                                       buffer); CHECK_ERR
        start[0] = 1;
        err = ncmpi_put_vars_short_all(ncid, rvarid[i], start, count, stride,
                                       buffer); CHECK_ERR
    }

    /* all processes read and verify */
    if (rank > 0) for (i=0; i<NDIMS; i++) count[i]  = 2;
    for (nelms=1,i=0; i<NDIMS; i++) nelms *= count[i];

    for (i=0; i<NVARS; i++) {
        for (j=0; j<nelms; j++) buffer[j] = -2;
        start[0] = 0;
        err = ncmpi_get_vars_short_all(ncid, fvarid[i], start, count, stride,
                                       buffer); CHECK_ERR
        for (j=0; j<nelms; j++) {
            if (buffer[j] != j%32768) {
                printf("Error at line %d: expect buffer[%d][%d]=%d but got %hd\n",
                       __LINE__, i, j, j%32768, buffer[j]);
                nerrs++;
                break;
            }
        }
        for (j=0; j<nelms; j++) buffer[j] = -2;
        start[0] = 1;
        err = ncmpi_get_vars_short_all(ncid, rvarid[i], start, count, stride,
                                       buffer); CHECK_ERR
        for (j=0; j<nelms; j++) {
            if (buffer[j] != j%32768) {
                printf("Error at line %d: expect buffer[%d][%d]=%d but got %hd\n",
                       __LINE__, i, j, j%32768, buffer[j]);
                nerrs++;
                break;
            }
        }
    }
    free(buffer);

fn_exit:
    err = ncmpi_close(ncid); CHECK_ERR

    /* check if open to read header fine */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}
