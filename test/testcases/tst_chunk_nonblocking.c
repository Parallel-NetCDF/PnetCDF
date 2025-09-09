/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests chunking feature when using nonblocking APIs and one of
 * the processes makes no call to the API.
 *
 * Contributed by Danqing Wu.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */

#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>


#define DIM_LEN 8

int main(int argc, char **argv)
{
    char filename[256];
    int err, nerrs=0, ncid, dimid, varid, rank, req, verbose=0;
    int vals[DIM_LEN] = {-1, -2, -3, -4, -5, -6, -7, -8};
    MPI_Offset start, count;
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
        sprintf(cmd_str, "*** TESTING C   %s for chunking iput ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_chunking", "enable");

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);
    CHECK_ERR
    MPI_Info_free(&info);

    err = ncmpi_def_dim(ncid, "x", DIM_LEN, &dimid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 1, &dimid, &varid);
    CHECK_ERR

    err = ncmpi_enddef(ncid);
    CHECK_ERR

    if (rank == 0)
    {
        start = 0;
        count = DIM_LEN;
        err = ncmpi_iput_vara_int(ncid, varid, &start, &count, vals, &req);
        CHECK_ERR
    }
    else
        req = NC_REQ_NULL;

    if (verbose) printf("rank = %d, before ncmpi_wait_all\n", rank); fflush(stdout);
    err = ncmpi_wait_all(ncid, 1, &req, NULL);
    CHECK_ERR
    if (verbose) printf("rank = %d, after ncmpi_wait_all\n", rank); fflush(stdout);

    err = ncmpi_close(ncid);

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
