/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program try reading the output from the NetCDF example program
 * simple_xy_nc4_wr.c
 *
 * Code in simple_xy_nc4_wr.c is used to generate the file, we then read it
 * back with PnetCDF
 *
 * Original description of the example:
 *
 * Write the simple_xy file, with some of the features of netCDF-4.
 *
 * This is a very simple example which is based on the simple_xy example,
 * but whch uses netCDF-4 features, such as compression. Please see the
 * simple_xy example to learn more about the netCDF-3 API.
 *
 * Like simple_xy_wr.c, this program writes a 2D netCDF variable (calle
 * "data") and fills it with sample data.  It has two dimensions, "x" and
 * "y".
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h> /* setenv() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <pnetcdf.h>

#include <testutils.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "simple_xy.nc"

/* We are reading 2D data, a 6 x 12 grid. */
#define NX 6
#define NY 12
#define NDIMS 2

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main(int argc, char** argv) {
    char filename[256];
    int i, j, nerrs=0, rank, nprocs, err;
    int ncid, x_dimid, y_dimid, varid, ndim;
    int dimids[2];
    int data_out[NX][NY];
    int data_in[NX][NY];
    MPI_Offset dlen;
    char tmp[1024];
    int x, y;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        nerrs++;
        goto fn_exit;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, FILE_NAME);
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str,
        "*** TESTING C   %s for opening and reading a netcdf4 file",
        basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* Create some pretend data. If this wasn't an example program, we
     * would have some real data to write, for example, model output. */
    for (x = 0; x < NX; x++)
        for (y = 0; y < NY; y++)
            data_out[x][y] = x * NY + y;

    /* Create the file. The NC_NETCDF4 parameter tells netCDF to create
     * a file in netCDF-4/HDF5 standard. */

    /* Note NC_MPIIO is used in NetCDF 4.6.1 and earlier, but ignored in 4.6.2
     * and after. */
    if ((err = nc_create_par(filename, NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid)))
        CHECK_ERR

    /* Define the dimensions. */
    if ((err = nc_def_dim(ncid, "x", NX, &x_dimid)))
        CHECK_ERR
    if ((err = nc_def_dim(ncid, "y", NY, &y_dimid)))
        CHECK_ERR

    /* Set up variabe data. */
    dimids[0] = x_dimid;
    dimids[1] = y_dimid;

    /* Define the variable. */
    if ((err = nc_def_var(ncid, "data", NC_INT, NDIMS,
                                dimids, &varid)))
        CHECK_ERR

    /* No need to explicitly end define mode for netCDF-4 files. Write
        * the pretend data to the file. */
    if ((err = nc_put_var_int(ncid, varid, &data_out[0][0])))
        CHECK_ERR

    /* Close the file. */
    if ((err = nc_close(ncid)))
        CHECK_ERR

    /* Open with PnetCDF */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* Check if number of dimensions matches */
    err = ncmpi_inq_ndims(ncid, &ndim); CHECK_ERR
    if (ndim != 2) {
        printf("Error at line %d in %s: expect ndim = %d, but got %d\n",
        __LINE__, __FILE__, 2, ndim);
    }

    /* Check if dimension size matches */
    err = ncmpi_inq_dim(ncid, x_dimid, tmp, &dlen); CHECK_ERR
    if (dlen != NX) {
        printf("Error at line %d in %s: expect X dim size = %d, but got %d\n",
        __LINE__, __FILE__, NX, (int)dlen);
    }
    err = ncmpi_inq_dim(ncid, y_dimid, tmp, &dlen); CHECK_ERR
    if (dlen != NY) {
        printf("Error at line %d in %s: expect X dim size = %d, but got %d\n",
        __LINE__, __FILE__, NY, (int)dlen);
    }

    /* Check if data matches */
    err = ncmpi_get_var_int_all(ncid, varid, &data_in[0][0]);
    for(i = 0; i < NX; i++){
        for(j = 0; j < NY; j++){
            if (data_out[i][j] != data_in[i][j]){
                printf("Error at line %d in %s: expect data_in[%d][%d] = %d, but got %d\n",
                __LINE__, __FILE__, i, j, data_out[i][j], data_in[i][j]);
            }
        }
    }

    /* Close the file */
    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR && malloc_size > 0) /* this test is for running 1 process */
        printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
               malloc_size);

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

