/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* This program tests if PnetCDF can report correct file formats */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

int main(int argc, char **argv) {
    char dir_name[256], filename[512];
    int err, rank, nerrs=0, format, ncid;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s dir_name\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(dir_name, 256, "%s", argv[1]);
    else           strcpy(dir_name, ".");
    MPI_Bcast(dir_name, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for inquiring file formats ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    /* test CDF-1 -----------------------------------------------------------*/
    sprintf(filename,"%s/test_cdf1.nc",dir_name);
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_ERR

    /* test NULL argument */
    err = ncmpi_inq_format(ncid, NULL); CHECK_ERR

    err = ncmpi_inq_format(ncid, &format); CHECK_ERR
    if (format != NC_FORMAT_CLASSIC) {
        printf("Error at line %d in %s: expecting CDF-1 format for file %s but got %d\n",
               __LINE__,__FILE__,filename,format);
        nerrs++;
    }
    err = ncmpi_close(ncid); CHECK_ERR

    /* test NULL argument */
    err = ncmpi_inq_file_format(filename, NULL); CHECK_ERR

    err = ncmpi_inq_file_format(filename, &format); CHECK_ERR
    if (format != NC_FORMAT_CLASSIC) {
        printf("Error at line %d in %s: expecting CDF-1 format for file %s but got %d\n",
               __LINE__,__FILE__,filename,format);
        nerrs++;
    }

    /* test CDF-2 -----------------------------------------------------------*/
    sprintf(filename,"%s/test_cdf2.nc",dir_name);
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_ERR

    /* test NULL argument */
    err = ncmpi_inq_format(ncid, NULL); CHECK_ERR

    err = ncmpi_inq_format(ncid, &format); CHECK_ERR
    if (format != NC_FORMAT_CDF2) {
        printf("Error at line %d in %s: expecting CDF-2 format for file %s but got %d\n",
               __LINE__,__FILE__,filename,format);
        nerrs++;
    }
    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_inq_file_format(filename, &format); CHECK_ERR
    if (format != NC_FORMAT_CDF2) {
        printf("Error at line %d in %s: expecting CDF-2 format for file %s but got %d\n",
               __LINE__,__FILE__,filename,format);
        nerrs++;
    }

    /* test CDF-5 -----------------------------------------------------------*/
    sprintf(filename,"%s/test_cdf5.nc",dir_name);
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_ERR

    /* test NULL argument */
    err = ncmpi_inq_format(ncid, NULL); CHECK_ERR

    err = ncmpi_inq_format(ncid, &format); CHECK_ERR
    if (format != NC_FORMAT_CDF5) {
        printf("Error at line %d in %s: expecting CDF-5 format for file %s but got %d\n",
               __LINE__,__FILE__,filename,format);
        nerrs++;
    }
    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_inq_file_format(filename, &format); CHECK_ERR
    if (format != NC_FORMAT_CDF5) {
        printf("Error at line %d in %s: expecting CDF-5 format for file %s but got %d\n",
               __LINE__,__FILE__,filename,format);
        nerrs++;
    }

    /* test NetCDF4 --------------------------------------------------------*/
    sprintf(filename,"%s/test_netcdf4.nc",dir_name);
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    if (PNETCDF_DRIVER_NETCDF4 == 0)
        EXP_ERR(NC_ENOTBUILT)
    else { /* NetCDF-4 is enabled */
        CHECK_ERR

        /* test NULL argument */
        err = ncmpi_inq_format(ncid, NULL); CHECK_ERR

        err = ncmpi_inq_format(ncid, &format); CHECK_ERR
        if (format != NC_FORMAT_NETCDF4) {
            printf("Error at line %d in %s: expecting NetCDF-4 format for file %s but got %d\n",
                   __LINE__,__FILE__,filename,format);
            nerrs++;
        }
        err = ncmpi_close(ncid); CHECK_ERR

        /* test NULL argument */
        err = ncmpi_inq_file_format(filename, NULL); CHECK_ERR

        err = ncmpi_inq_file_format(filename, &format); CHECK_ERR
        if (format != NC_FORMAT_NETCDF4) {
            printf("Error at line %d in %s: expecting NETCDF4 format for file %s but got %d\n",
                   __LINE__,__FILE__,filename,format);
            nerrs++;
        }
    }

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
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}
