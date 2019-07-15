/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to nonblocking IO with burst buffer driver.
 * It is same as using without the burst buffer driver with the only exception that
 * we may not be able to cancel nonblocking requests
 * In this example, every process write its rank to n-th cell in 3 1 X N variables.
 * N is the number of processes.
 * We will try to cancel the write operations on variable A, and C
 * While we can cancel the request for variable A, we can not do so for variable C
 * because it is already flushed to PFS when we wait on request of variable B
 *
 *    To compile:
 *        mpicc -O2 nonblocking.c -o nonblocking -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * netCDF file produced by this example program:
 *
 *    % mpiexec -n 4 ./nonblocking -b /scratch testfile.nc
 *    nonblocking.c: example of nonblocking IO
 *    Canceling write on variable A, get 0 (No error)
 *    Waiting on variable B, get 0 (No error)
 *    Canceling write on variable C, get -306 (Nonblocking requests already flushed.)
 *
 *    % ncmpidump testfile.nc
 *    netcdf test {
 *    // file format: CDF-1
 *    dimensions:
 *            X = 1 ;
 *            Y = 4 ;
 *    variables:
 *            int A(X, Y) ;
 *            int B(X, Y) ;
 *            int C(X, Y) ;
 *    data:
 *
 *     A =
 *      _, _, _, _ ;
 *
 *     B =
 *      0, 1, 2, 3 ;
 *
 *     C =
 *      0, 1, 2, 3 ;
 *    }
 *
 * Example batch script for running on Cori at NERSC using SLURM scheduler
 *
 * #!/bin/bash
 * #SBATCH -p debug
 * #SBATCH -N 1
 * #SBATCH -C haswell
 * #SBATCH -t 00:01:00
 * #SBATCH -o nonblocking_example.txt
 * #SBATCH -L scratch
 * #DW jobdw capacity=1289GiB access_mode=private type=scratch pool=sm_pool
 * srun -n 4 ./nonblocking -b $DW_JOB_PRIVATE ${SCRATCH}/testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy(), strdup() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [-b bb_dir] file_name\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-b bb_dir] Path to burst buffer\n"
    "       filename: output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256], *bb_dir;
    int i, rank, np, verbose=1, err, nerrs=0;
    int dimid[2];
    int varid[3];
    MPI_Offset start[2];
    int buffer[3];
    int req[3];
    int stat;
    int ncid, cmode;
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    /* get command-line arguments */
    bb_dir = strdup(".");
    while ((i = getopt(argc, argv, "hqb:")) != EOF)
        switch(i) {
            case 'b': free(bb_dir);
                      bb_dir = strdup(optarg);
                      break;
            case 'q': verbose = 0;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    if (verbose && rank == 0) printf("%s: example of nonblocking IO\n",__FILE__);

    /* Set up the hints for burst buffer driver in ncmpi_create
     * Note that the remaining part of the code remains unchanged
     * PnetCDF will warn if nc_burst_buf_dirname is not set.
     */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_burst_buf", "enable");
    MPI_Info_set(info, "nc_burst_buf_dirname", bb_dir);
    free(bb_dir);
    /* create a new file using clobber mode ----------------------------------*/
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    ERR

    /* Info can be freed after file creation */
    MPI_Info_free(&info);

    /* Fill up the variables with default value */
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); ERR

    /* Define dimensions */
    err = ncmpi_def_dim(ncid, "X", 1, dimid);    ERR
    err = ncmpi_def_dim(ncid, "Y", np, dimid + 1);    ERR

    /* Define variables */
    err = ncmpi_def_var(ncid, "A", NC_INT, 2, dimid, varid);    ERR
    err = ncmpi_def_var(ncid, "B", NC_INT, 2, dimid, varid + 1);    ERR
    err = ncmpi_def_var(ncid, "C", NC_INT, 2, dimid, varid + 2);    ERR

    /* burst buffer initialize log files on the first time we enters data mode */
    err = ncmpi_enddef(ncid);
    ERR

    /* Preparing the buffer */
    buffer[0] = buffer[1] = buffer[2] = rank;

    /* Do nonblocking write on variables */
    start[0] = 0;
    start[1] = rank;
    err = ncmpi_iput_var1_int(ncid, varid[0], start, buffer, req);    ERR
    err = ncmpi_iput_var1_int(ncid, varid[1], start, buffer + 1, req + 1);    ERR
    err = ncmpi_iput_var1_int(ncid, varid[2], start, buffer + 2, req + 2);    ERR

    /* Cancel first request */
    err = ncmpi_cancel(ncid, 1, req, &stat); ERR
    if (verbose && rank == 0)
        printf("Canceling write on variable A, get %d (%s)\n", stat, ncmpi_strerror(stat));

    /* Wait second request */
    err = ncmpi_wait_all(ncid, 1, req + 1, &stat); ERR
    if (verbose && rank == 0)
        printf("Waiting on variable B, get %d (%s)\n", stat, ncmpi_strerror(stat));

    /* Cancel third request */
    err = ncmpi_cancel(ncid, 1, req + 2, &stat); ERR
    if (verbose && rank == 0)
        printf("Canceling write on variable C, get %d (%s)\n", stat, ncmpi_strerror(stat));

    /* close file */
    err = ncmpi_close(ncid);
    ERR

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

