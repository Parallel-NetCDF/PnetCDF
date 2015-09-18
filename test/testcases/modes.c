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
#include <unistd.h> /* unlink() */
#include <mpi.h>
#include <pnetcdf.h>

#define FAIL_COLOR "\x1b[31mfail\x1b[0m\n"
#define PASS_COLOR "\x1b[32mpass\x1b[0m\n"

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}

#define EXPECT_ERR(err_no) \
    if (err != err_no) { \
        nfails++; \
        printf("Error at line %d: expect error code %d but got %d\n", \
               __LINE__,err_no,err); \
    }

int main(int argc, char** argv)
{
    char filename[256];
    int rank, nprocs, err, nfails=0;
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

    /* no file should be created when in safe mode */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    EXPECT_ERR(NC_ENOENT)

    /* test under safe mode disabled */
    setenv("PNETCDF_SAFE_MODE", "0", 1);
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    EXPECT_ERR(NC_EINVAL_CMODE)

    /* file should be created and is in CDF-5 format */
    err = ncmpi_inq_format(ncid, &format);
    ERR
    if (format != 5) {
        printf("Error at line=%d: expecting CDF-5 format for file %s but got CDF-%d\n",__LINE__,filename,format);
        nfails++;
    }
    err = ncmpi_close(ncid);
    ERR

    MPI_Allreduce(MPI_IN_PLACE, &nfails, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    char cmd_str[256];
    sprintf(cmd_str, "*** TESTING C   %s for file create/open modes ", argv[0]);
    if (rank == 0) {
        if (nfails) printf("%-66s ------ " FAIL_COLOR, cmd_str);
        else        printf("%-66s ------ " PASS_COLOR, cmd_str);
    }


    MPI_Finalize();
    return 0;
}

