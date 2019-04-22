/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
   This is an example program which writes a 1-D compressed array

   $Id$
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>
#include <mpi.h>
#include <testutils.h>

/* This is the name of the data file we will create. */
#define FILE_NAME "debug.nc"

#define N 5

int main(int argc, char **argv)
{
    /* IDs for the netCDF file, dimensions, and variables. */
    int i, j;
    int np, rank, nerrs = 0;
    int ncid, dimids[2], varid;
    int buf[N * N];
    int reqids[N];
    MPI_Offset start[2], count[2];
    MPI_Info info;

    /* Error handling. */
    int err;

    char *filename = FILE_NAME;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 3) {
        if (!rank)
            printf("Usage: %s [filename]\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    if (argc > 1) filename = argv[1];

    if (rank == 0) {
        char *cmd_str = (char *)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for writing compressed file", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    for(j = 0; j < 2; j++){
        /* Initialize file info */
        MPI_Info_create(&info);
        MPI_Info_set(info, "nc_compression", "enable");
        MPI_Info_set(info, "nc_zip_delay_init", "1");
        switch(j){
            case 0:
                MPI_Info_set(info, "nc_zip_comm_unit", "chunk");
                break;
            case 1:
                MPI_Info_set(info, "nc_zip_comm_unit", "proc");
                break;
        }

        /* Create the file. */
        err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);
        CHECK_ERR

        /* Define the dimension. */
        err = ncmpi_def_dim(ncid, "X", np * N, dimids);
        CHECK_ERR
        err = ncmpi_def_dim(ncid, "Y", np * N, dimids + 1);
        CHECK_ERR
        
        /* Define the variable. */
        err = ncmpi_def_var(ncid, "M", NC_INT, 2, dimids, &varid);
        CHECK_ERR

        /* End define mode. */
        err = ncmpi_enddef(ncid);
        CHECK_ERR

        // Init buffer
        for(i = 0; i < N * N; i++){
            buf[i] = rank + i + 1;
        }

        // Write variable
        count[0] = N;
        count[1] = N;
        for(i = 0; i < np; i++){
            start[0] = i * N ;
            start[1] = ((i + rank) % np) * N ;
            err = ncmpi_iput_vara_int(ncid, varid, start, count, buf, reqids + i); CHECK_ERR
        }

        /* Wait for all request */
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); CHECK_ERR

        /* Close the file. */
        err = ncmpi_close(ncid);
        CHECK_ERR

        /* Open the file. */
        err = ncmpi_open(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);
        CHECK_ERR

        // Free info
        MPI_Info_free(&info);

        // Read variable
        count[0] = N;
        count[1] = N;
        for(i = 0; i < np; i++){
            memset(buf, 0, sizeof(buf));
            start[0] = i * N ;
            start[1] = ((i + rank) % np) * N ;
            err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

            // Check results
            for(j = 0; j < N * N; j++){
                if (buf[j] != rank + j + 1){
                    printf("Error at %s:%d: expect round %d buf[%d]=%d but got %d\n",
                        __FILE__,__LINE__, i , rank + j + 1, buf[j]);
                }
            }
        }

        /* Close the file. */
        err = ncmpi_close(ncid);
        CHECK_ERR
    }

    /* check if there is any malloc residue */
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
        if (nerrs)
            printf(FAIL_STR, nerrs);
        else
            printf(PASS_STR);
    }

    MPI_Finalize();

    return (nerrs > 0);
}
