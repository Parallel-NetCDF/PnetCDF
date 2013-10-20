/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* This program tests if PnetCDF can report correct file formats */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>

#define ERR {if(err!=NC_NOERR) {printf("Error(%d) at line %d: %s\n",err,__LINE__,ncmpi_strerror(err)); nerr++; }}

int main(int argc, char **argv) {
    int err, rank, nerr=0, format, ncid;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank > 0) {
        MPI_Finalize();
        return 0;
    }

    err = ncmpi_open(MPI_COMM_WORLD, "../data/test_int.nc", 0, MPI_INFO_NULL,
                     &ncid);
    ERR
 
    err = ncmpi_inq_format(ncid, &format); ERR
    if (format != 1) {
        printf("Error (line=%d): expecting CDF-1 format for file ../data/test.nc but got %d\n",__LINE__,format);
        nerr++;
    }
    err = ncmpi_close(ncid); ERR
  
    err = ncmpi_inq_file_format("../data/test_int_cdf5.nc", &format); ERR
    if (format != 5) {
        printf("Error (line=%d): expecting CDF-5 format for file ../data/test_int_cdf5.nc but got %d\n",__LINE__,format);
        nerr++;
    }

    char cmd_str[80];
    sprintf(cmd_str, "*** TESTING C   %s for inquiring CDF file format ", argv[0]);
    if (nerr == 0)
        printf("%-66s ------ pass\n", cmd_str);
    else
        printf("%-66s ------ failed\n", cmd_str);

    MPI_Finalize();
    return 0;
}
