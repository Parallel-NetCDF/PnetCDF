/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
   This is an example program which writes a 2-D compressed array

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

#define NZIPDRIVER 4

#define N 10

int main(int argc, char **argv)
{
    /* IDs for the netCDF file, dimensions, and variables. */
    int i;
    int zipdriver, communit;
    int np, rank, nerrs = 0;
    int ncid, dimids[2], varid;
    int buf[N];
    MPI_Offset start[N][2], count[N][2];
    MPI_Offset *starts[N], *counts[N];
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

    for(zipdriver = 0; zipdriver < NZIPDRIVER; zipdriver++){
        for(communit = 0; communit < 2; communit++){
            /* Initialize file info */
            MPI_Info_create(&info);
            MPI_Info_set(info, "nc_compression", "enable");
            switch(communit){
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
            err = ncmpi_def_dim(ncid, "X", np, dimids);
            CHECK_ERR
            err = ncmpi_def_dim(ncid, "Y", N, dimids + 1);
            CHECK_ERR
            
            /* Define the variable. */
            err = ncmpi_def_var(ncid, "M", NC_INT, 2, dimids, &varid);
            CHECK_ERR

            /* Select compression driver. */
            buf[0] = zipdriver;
            err = ncmpi_put_att_int(ncid, varid, "_zipdriver", NC_INT, 1, buf);
            CHECK_ERR

            /* Define chunk size. */
            buf[0] = np;
            buf[1] = 5;
            err = ncmpi_put_att_int(ncid, varid, "_chunkdim", NC_INT, 2, buf);
            CHECK_ERR
            
            /* End define mode. */
            err = ncmpi_enddef(ncid);
            CHECK_ERR

            // Write variable
            for(i = 0; i < N; i++){
                start[i][0] = rank;
                start[i][1] = i;
                count[i][0] = 1;
                count[i][1] = 1;
                starts[i] = start[i];
                counts[i] = count[i];
                buf[i] = rank * N + i + 1;
            }
            err = ncmpi_iput_varn_int(ncid, varid, N, starts, counts, buf, NULL);

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
            memset(buf, 0, sizeof(buf));
            for(i = 0; i < N; i++){
                start[i][0] = rank;
                start[i][1] = i;
                count[i][0] = 1;
                count[i][1] = 1;
                starts[i] = start[i];
                counts[i] = count[i];
            }
            err = ncmpi_iget_varn_int(ncid, varid, N, starts, counts, buf, NULL);

            /* Wait for all request */
            err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); CHECK_ERR
            
            // Check results
            for(i = 0; i < N; i++){
                if (buf[i] != rank * N + i + 1){
                    printf("Error at %s:%d: expect buf[%d]=%d but got %d\n",
                        __FILE__,__LINE__, i , rank * N + i + 1, buf[i]);
                }
            }

            /* Close the file. */
            err = ncmpi_close(ncid);
            CHECK_ERR
        }
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
