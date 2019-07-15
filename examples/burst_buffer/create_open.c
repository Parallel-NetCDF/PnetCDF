/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to create/open the file using the burst buffer driver.
 * It is a modified version of create_open.c under examples/C using the burst buffer driver.
 *
 *    To compile:
 *        mpicc -O2 create_open.c -o create_open -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * netCDF file produced by this example program:
 *
 *    % mpiexec -n 4 ./create_open -b /scratch testfile.nc
 *    create_open.c: example of file create and open
 *
 *    % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    }
 *
 * Example batch script for running on Cori at NERSC using SLURM scheduler
 *
 * #!/bin/bash
 * #SBATCH -p debug
 * #SBATCH -N 1
 * #SBATCH -C haswell
 * #SBATCH -t 00:01:00
 * #SBATCH -o create_open_example.txt
 * #SBATCH -L scratch
 * #DW jobdw capacity=1289GiB access_mode=private type=scratch pool=sm_pool
 * srun -n 4 ./create_open -b $DW_JOB_PRIVATE $SCRATCH/testfile.nc
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
    int i, rank, verbose=1, err, nerrs=0, ncid;
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* get command-line arguments */
    bb_dir = strdup(".");
    while ((i = getopt(argc, argv, "hqb:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'h':
            case 'b': free(bb_dir);
                      bb_dir = strdup(optarg);
                      break;
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    if (verbose && rank == 0) printf("%s: example of file create and open\n",__FILE__);

    /* Set up the hints for burst buffer driver in ncmpi_create
     * Note that the remaining part of the code remains unchanged
     * The burst buffer driver will not proceed if the log files already exists
     * to prevent overwriting existing files by accident
     * To open the file again, we need to delete the log file after file closing
     * The default value of nc_burst_buf_del_on_close is enable, we set it for the
     * purpose of demonstration.
     * PnetCDF will warn if nc_burst_buf_dirname is not set.
     */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_burst_buf", "enable");
    MPI_Info_set(info, "nc_burst_buf_del_on_close", "enable");
    MPI_Info_set(info, "nc_burst_buf_dirname", bb_dir);

    /* create a new file using clobber mode ----------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid); ERR

    /* Info can be freed after file creation */
    MPI_Info_free(&info);

    if (err != NC_NOERR) goto fn_exit;

    /* burst buffer initialize log files on the first time we enters data mode */
    err = ncmpi_enddef(ncid);
    ERR

    /* close file */
    err = ncmpi_close(ncid);
    ERR

    /* Set up the hints for burst buffer driver in ncmpi_create
     * Note that the remaining part of the code remains unchanged
     * PnetCDF will warn if nc_burst_buf_dirname is not set.
     */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_burst_buf", "enable");
    MPI_Info_set(info, "nc_burst_buf_dirname", bb_dir);
    free(bb_dir);

    /* open the newly created file for read only -----------------------------*/
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_WRITE, info, &ncid);
    ERR

    /* Info can be freed after file opening */
    MPI_Info_free(&info);

    if (err != NC_NOERR) goto fn_exit;

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

fn_exit:
    MPI_Finalize();
    return (nerrs > 0);
}

