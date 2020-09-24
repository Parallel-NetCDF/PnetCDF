/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests nonblocking functionality of bb driver by making multiple
 * nonblocking request to test the driver's ability to handle large amount of
 * nonblocking requests.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <limits.h>
#include <testutils.h>
#include <libgen.h>

#define NREQ 2048
#define NROUND 4

int main(int argc, char *argv[]) {
    int i, j, err, nerrs = 0;
    int rank, np;
    int ncid, varid;
    int *buffer, *reqs, *stat;
    int dimid[2];
    char *filename;
    MPI_Offset start[2];
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
        sprintf(cmd_str, "*** TESTING C   %s for burst buffer big requests", basename(argv[0]));
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
    err = ncmpi_def_dim(ncid, "Y", NREQ, dimid + 1);    CHECK_ERR

    /* Define variable */
    err = ncmpi_def_var(ncid, "M", NC_INT, 2, dimid, &varid);    CHECK_ERR

    /* Switch to data mode */
    err = ncmpi_enddef(ncid);    CHECK_ERR

    /* Prepare Buffer */
    reqs = (int*)malloc(sizeof(int) * NREQ);
    stat = (int*)malloc(sizeof(int) * NREQ);
    buffer = (int*)malloc(sizeof(int) * NREQ);
    for (i = 0; i < NREQ; i++) {
        buffer[i] = rank + 1;
    }

    for (j = 0; j < NROUND; j++) {
        /* Test nonblocking put */
        for (i = 0; i < NREQ; i++) {
            start[0] = rank;
            start[1] = i;
            err = ncmpi_iput_var1_int(ncid, varid, start, buffer + i, reqs + i);    CHECK_ERR
        }
        err = ncmpi_wait_all(ncid, NREQ, reqs, stat);    CHECK_ERR
        for (i = 0; i < NREQ; i++) {
            err = stat[i];    CHECK_ERR
        }
        for (i = 0; i < NREQ; i++) {
            if (reqs[i] != NC_REQ_NULL) {
                printf("Error at line %d in %s: expecting reqs[%d] = NC_REQ_NULL but got %d\n", __LINE__, __FILE__, i, reqs[i]);
            }
        }

        /* Test nonblocking get */
        memset(buffer, 0, sizeof(int) * NREQ);
        for (i = 0; i < NREQ; i++ ) {
            start[0] = rank;
            start[1] = i;
            err = ncmpi_iget_var1_int(ncid, varid, start, buffer + i, reqs + i);    CHECK_ERR
        }
        err = ncmpi_wait_all(ncid, NREQ, reqs, stat);    CHECK_ERR
        for (i = 0; i < NREQ; i++) {
            err = stat[i];    CHECK_ERR
        }
        for (i = 0; i < NREQ; i++) {
            if (buffer[i] != rank + 1) {
                printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n", __LINE__, __FILE__, i, rank + 1, buffer[i]);
            }
        }
        for (i = 0; i < NREQ; i++) {
            if (reqs[i] != NC_REQ_NULL) {
                printf("Error at line %d in %s: expecting reqs[%d] = NC_REQ_NULL but got %d\n", __LINE__, __FILE__, i, reqs[i]);
            }
        }
    }

    /* Close the file */
    err = ncmpi_close(ncid);    CHECK_ERR

    MPI_Info_free(&info);
    free(buffer);
    free(reqs);
    free(stat);

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
