/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program (borrowed from netCDF library) tests defining and inquiring the
 * maximum allowable dimension size for CDF-1, 2, and 5 formats.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_dimsizes tst_dimsizes.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_dimsizes testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

#define DIMMAXCLASSIC (NC_MAX_INT - 3)
#define DIMMAX64OFFSET (NC_MAX_UINT - 3)
#define DIMMAX64DATA NC_MAX_INT64

/*
 * NC_CLASSIC => NC_INT_MAX - 3
 * NC_64BIT_OFFSET => NC_UINT_MAX - 3
 * NC_64BIT_DATA => NC_INT64_MAX
 * Note that for NC_64BIT_DATA, the max dimension size is different from netCDF
 * library. This is because PnetCDF uses MPI_Offset for dimension size and
 * MPI_Offset is a signed long long.
*/

int
main(int argc, char **argv)
{
    char filename[256];
    int rank, nprocs, err, nerrs=0;
    int ncid, dimid;
    MPI_Offset dimsize;

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
        sprintf(cmd_str, "*** TESTING C   %s for defining max dimension sizes ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* Writing Max Dimension Size For NC_CLASSIC */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR
    dimsize = DIMMAXCLASSIC;
    err = ncmpi_def_dim(ncid, "testdim", dimsize, &dimid); CHECK_ERR
    dimsize = -1;
    err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)
    dimsize = DIMMAXCLASSIC+1;
    err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)
    err = ncmpi_close(ncid); CHECK_ERR

    /* Reading Max Dimension Size For NC_CLASSIC */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOCLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_inq_dimid(ncid, "testdim", &dimid); CHECK_ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &dimsize); CHECK_ERR
    if (dimsize != DIMMAXCLASSIC) {
        printf("Error at line %d in %s: expecting dimsize %d but got %lld\n", __LINE__,__FILE__,DIMMAXCLASSIC,dimsize);
        nerrs++;
    }
    err = ncmpi_close(ncid); CHECK_ERR

    /* Writing Max Dimension Size For NC_64BIT_OFFSET */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER | NC_64BIT_OFFSET, MPI_INFO_NULL, &ncid); CHECK_ERR
    dimsize = DIMMAX64OFFSET;
    err = ncmpi_def_dim(ncid, "testdim", dimsize, &dimid); CHECK_ERR
    dimsize = -1;
    err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)
    dimsize = DIMMAX64OFFSET+1;
    err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)
    err = ncmpi_close(ncid); CHECK_ERR

    /* Reading Max Dimension Size For NC_64BIT_OFFSET */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOCLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_inq_dimid(ncid, "testdim", &dimid); CHECK_ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &dimsize); CHECK_ERR
    if (dimsize != DIMMAX64OFFSET) {
        printf("Error at line %d in %s: expecting dimsize %d but got %lld\n", __LINE__,__FILE__,DIMMAX64OFFSET,dimsize);
        nerrs++;
    }
    err = ncmpi_close(ncid); CHECK_ERR

    /* Writing Max Dimension Size For NC_64BIT_DATA */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER | NC_64BIT_DATA, MPI_INFO_NULL, &ncid); CHECK_ERR
    dimsize = DIMMAX64DATA;
    err = ncmpi_def_dim(ncid, "testdim", dimsize, &dimid); CHECK_ERR
    dimsize = -1;
    err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)
    err = ncmpi_close(ncid); CHECK_ERR

    /* Reading Max Dimension Size For NC_64BIT_DATA */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOCLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_inq_dimid(ncid, "testdim", &dimid); CHECK_ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &dimsize); CHECK_ERR
    if (dimsize != DIMMAX64DATA) {
        printf("Error at line %d in %s: expecting dimsize %lld but got %lld\n", __LINE__,__FILE__,DIMMAX64DATA,dimsize);
        nerrs++;
    }
    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}
