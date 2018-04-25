/*********************************************************************
 *
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to use NC_REQ_ALL in nonblocking I/O operations.
 * The program writes 2 arrays by calling the nonblocking APIs with NULLs for
 * argument request ID. When calling ncmpi_wait_all(), NC_REQ_ALL is used to
 * commit all the pending requests without checking the individual statuses of
 * the requests.
 *
 *    To compile:
 *        mpicc -O2 req_all.c -o req_all -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * NC file produced by this example program:
 *
 *    % mpiexec -n 4 ./req_all /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    dimensions:
 *            Y = 8 ;
 *            X = 8 ;
 *    variables:
 *            int var_int(Y, X) ;
 *            float var_flt(Y, X) ;
 *    data:
 *
 *     var_int =
 *      10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13 ;
 *
 *     var_flt =
 *      10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5,
 *      10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5,
 *      10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5,
 *      10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5,
 *      10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5,
 *      10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5,
 *      10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5,
 *      10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 13.5, 13.5 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 8
#define NX 2

int main(int argc, char** argv)
{
    char filename[256];
    int i, j, rank, nprocs, nerrs=0, err;
    int ncid, cmode, varid[2], dimid[2], buf_int[NY][NX];
    float buf_flt[NY][NX];
    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[2], count[2];

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
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for NC_REQ_ALL ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* the global array is NY * (NX * nprocs) */
    global_ny = NY;
    global_nx = NX * nprocs;

    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++) {
             buf_int[i][j] = rank+10;
             buf_flt[i][j] = 10.5+rank;
        }

    /* define dimensions x and y */
    err = ncmpi_def_dim(ncid, "Y", global_ny, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", global_nx, &dimid[1]); CHECK_ERR

    /* define a 2D variable of integer type */
    err = ncmpi_def_var(ncid, "var_int", NC_INT, 2, dimid, &varid[0]); CHECK_ERR

    /* define a 2D variable of float type */
    err = ncmpi_def_var(ncid, "var_flt", NC_FLOAT, 2, dimid, &varid[1]); CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* now we are in data mode */
    start[0] = 0;
    start[1] = NX * rank;
    count[0] = NY;
    count[1] = NX;

    err = ncmpi_iput_vara_int(ncid, varid[0], start, count, &buf_int[0][0], NULL); CHECK_ERR
    err = ncmpi_iput_vara_float(ncid, varid[1], start, count, &buf_flt[0][0], NULL); CHECK_ERR

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); CHECK_ERR

    /* check if write buffer contents have been altered */
    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++) {
             if (buf_int[i][j] != rank+10) {
                 printf("Error at line %d in %s: expecting buffer[%d][%d]=%d but got %d\n",
                       __LINE__,__FILE__,i,j,rank+10,buf_int[i][j]);
                 nerrs++;
             }
             if (buf_flt[i][j] != 10.5 + rank) {
                 printf("Error at line %d in %s: expecting buffer[%d][%d]=%f but got %f\n",
                       __LINE__,__FILE__,i,j,10.5+rank,buf_flt[i][j]);
                 nerrs++;
             }
        }

    err = ncmpi_close(ncid);
    CHECK_ERR

    /* check if there is any PnetCDF internal malloc residue */
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

