/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests if the correct error codes are returns given various
 * create/open modes.
 *
 * For example, NC_EINVAL_CMODE should be returned when creating a file using
 * comde with both NC_64BIT_OFFSET & NC_64BIT_DATA flags set.
 * When in safe mode, no file will be created.
 * When not in safe mode, a CDF-2 file will be created, as in this case the
 * flag NC_64BIT_OFFSET triumph NC_64BIT_DATA..
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <unistd.h> /* unlink(), access() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}

#define EXPECT_ERR(err_no) \
    if (err != err_no) { \
        nerrs++; \
        printf("Error at line %d: expect error code %s but got %s\n", \
               __LINE__,nc_err_code_name(err_no),nc_err_code_name(err)); \
    }

int main(int argc, char** argv)
{
    char filename[256];
    int rank, nprocs, err, nerrs=0, file_exist;
    int ncid, cmode, format;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    strcpy(filename, "testfile.nc");
    if (argc == 2) strcpy(filename, argv[1]);
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char cmd_str[256];
        sprintf(cmd_str, "*** TESTING C   %s for file create/open modes ", argv[0]);
        printf("%-66s ------ ", cmd_str); fflush(stdout);
    }

    /* create a new file and test various cmodes ----------------------------*/
    cmode = NC_CLOBBER;

    /* NC_64BIT_OFFSET and NC_64BIT_DATA should not appear together */
    cmode |= NC_64BIT_OFFSET | NC_64BIT_DATA;

    /* delete the file and ignore error */
    unlink(filename);

    /* test under safe mode enabled */
    setenv("PNETCDF_SAFE_MODE", "1", 1);
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    EXPECT_ERR(NC_EINVAL_CMODE)
    if (err == NC_NOERR) ncmpi_close(ncid);

    /* When opening a non-existing file for read, no file should be created */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    EXPECT_ERR(NC_ENOENT)

    file_exist = 0;
    if (rank == 0 && access(filename, F_OK) == 0) file_exist = 1;
    MPI_Bcast(&file_exist, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (file_exist) {
        printf("Error at line %d: opening a non-existing file (%s) creates the file by mistake\n", __LINE__, filename);
        nerrs++;
    }

    /* When opening a non-existing file for write, no file should be created */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_WRITE, MPI_INFO_NULL, &ncid);
    EXPECT_ERR(NC_ENOENT)

    file_exist = 0;
    if (rank == 0 && access(filename, F_OK) == 0) file_exist = 1;
    MPI_Bcast(&file_exist, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (file_exist) {
        printf("Error at line %d: opening a non-existing file (%s) creates the file by mistake\n", __LINE__, filename);
        nerrs++;
    }

    /* test under safe mode disabled */
    setenv("PNETCDF_SAFE_MODE", "0", 1);
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    EXPECT_ERR(NC_EINVAL_CMODE)

    /* file should be created and is in CDF-5 format */
    err = ncmpi_inq_format(ncid, &format);
    ERR
    if (format != 5) {
        printf("Error at line=%d: expecting CDF-5 format for file %s but got CDF-%d\n",__LINE__,filename,format);
        nerrs++;
    }
    err = ncmpi_close(ncid);
    ERR

    /* check if PnetCDF freed all internal malloc */
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
    return 0;
}

