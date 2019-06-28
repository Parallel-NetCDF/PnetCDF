/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * This program tests a large number of nonblocking requests (larger than
 * NC_REQUEST_CHUNK, the constant used to grow nonblocking put and get
 * queues.
 *
 *********************************************************************/
/*  $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <assert.h>
#include <pnetcdf.h>

#include <testutils.h>

#define FILE_NAME "testfile.nc"
#define NUM_REQS 1100   /* a number greater than NC_REQUEST_CHUNK */

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    int i, ncid, dimid[2], varid, err, nerrs=0, rank, nprocs;
    int *buf, *req, *status;
    char filename[256];
    MPI_Offset start[2], count[2];
    MPI_Info info=MPI_INFO_NULL;

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

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for large number of iput/iget ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
#define STRESS_ROMIO
#ifdef STRESS_ROMIO
    /* all processes write to the same file region. This can vigorously test
     * ROMIO for receiving write requests into the same local buffer at each
     * I/O aggregator.
     */
    err = ncmpi_def_dim(ncid, "X", 2,            &dimid[1]); CHECK_ERR
#else
    /* Writes from all processes are not overlapped */
    err = ncmpi_def_dim(ncid, "X", 2 * nprocs,   &dimid[1]); CHECK_ERR
#endif
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    req = (int*) malloc(NUM_REQS * 2 * sizeof(int));
    status = req + NUM_REQS;

    buf = (int*) calloc(NUM_REQS * 6,  sizeof(int));

    count[0] = 3; count[1] = 2;
    start[0] = 0;
#ifdef STRESS_ROMIO
    start[1] = 0;
#else
    start[1] = rank * 2;
#endif

    for (i=0; i<NUM_REQS; i++) {
        err = ncmpi_iput_vara_int(ncid, varid, start, count,
                                  &buf[i*6], &req[i]); CHECK_ERR
        start[0] += 3;
    }

    err = ncmpi_wait_all(ncid, NUM_REQS, req, status); CHECK_ERR
    /* check each iput status */
    for (i=0; i<NUM_REQS; i++) {
        err = status[i];
        CHECK_ERR
    }

    start[0] = 0;
    for (i=0; i<NUM_REQS; i++) {
        err = ncmpi_iget_vara_int(ncid, varid, start, count,
                                  &buf[i*6], &req[i]); CHECK_ERR
        start[0] += 3;
    }

    err = ncmpi_wait_all(ncid, NUM_REQS, req, status); CHECK_ERR
    /* check each iget status */
    for (i=0; i<NUM_REQS; i++) {
        err = status[i];
        CHECK_ERR
    }

    free(buf);
    free(req);
    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0) {
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
            ncmpi_inq_malloc_list();
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

