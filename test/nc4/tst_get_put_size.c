/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests whether get size and put size is counted accurately when 
 * accessing a NetCDF4 file.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define N 10

int main(int argc, char** argv) {
    int i;
    char filename[256];
    int rank, nprocs, err, nerrs=0;
    int ncid, cmode, varid[1], dimid[2], buf[1];
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Offset start[2];
    MPI_Offset size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

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
        sprintf(cmd_str, "*** TESTING C   %s for get size and put size ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* create a new file for writing ------------------------------------*/
    cmode = NC_CLOBBER | NC_NETCDF4;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR

    /* define dimension */
    err = ncmpi_def_dim(ncid, "X", nprocs, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", N, &dimid[0]); CHECK_ERR

    /* define 2 variables of size NY x NX */
    err = ncmpi_def_var(ncid, "V", NC_INT, 1, dimid, varid); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* write a subarray to each variable */
    start[0] = rank;
    for(i = 0; i < N; i++){
        start[1] = i;
        buf[0] = rank;
        err = ncmpi_put_var1_int_all(ncid, varid[0], start, buf); CHECK_ERR
        err = ncmpi_inq_put_size(ncid, &size); CHECK_ERR
        if (size != sizeof(int) * (i + 1)){
            printf("Error at line %d of %s: expect put_size = %ld but got %lld\n",
                __LINE__,__FILE__,sizeof(int) * (i + 1),size);
            nerrs++;
        }
        buf[0] = -1;
        err = ncmpi_get_var1_int_all(ncid, varid[0], start, buf); CHECK_ERR
        err = ncmpi_inq_get_size(ncid, &size); CHECK_ERR
        if (size != sizeof(int) * (i + 1)){
            printf("Error at line %d of %s: expect put_size = %ld but got %lld\n",
                __LINE__,__FILE__,sizeof(int) * (i + 1),size);
            nerrs++;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, comm);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

