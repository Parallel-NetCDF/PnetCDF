/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <time.h>   /* time() localtime(), asctime() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#define NY 10
#define NX 4

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to use ncmpi_put_vara_int_all() to write a 2D
 * 4-byte integer array in parallel. It first defines a netCDF variable of
 * size global_nx * global_ny where
 *    global_ny == NY and
 *    global_nx == (NX * number of MPI processes).
 * The data partitioning pattern is a column-wise partitioning across all
 * proceses. Each process writes a subarray of size ny * nx.
 *
 *    To compile:
 *        mpicc -O2 put_vara.c -o put_vara -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * NC file produced by this example program:
 *
 *    % mpiexec -n 4 ./put_vara /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 *            y = 10 ;
 *            x = 16 ;
 *    variables:
 *            int var(y, x) ;
 *                var:str_att_name = "example attribute of type text." ;
 *                var:float_att_name = 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f ;
 *    // global attributes:
 *                :history = "Wed Apr 30 11:18:58 2014\n",
 *       "" ;
 *    data:
 *
 *     var =
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 ;
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}

int main(int argc, char** argv) {
    char *filename="testfile.nc";
    int i, j, rank, nprocs, err;
    int ncid, cmode, varid, dimid[2], buf[NY][NX];
    char str_att[128];
    float float_att[100];
    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[2], count[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    if (argc == 2) filename = argv[1];

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* the global array is NY * (NX * nprocs) */
    global_ny = NY;
    global_nx = NX * nprocs;

    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++)
             buf[i][j] = rank;

    /* add a global attribute: a time stamp at rank 0 */
    time_t ltime = time(NULL); /* get the current calendar time */
    asctime_r(localtime(&ltime), str_att);

    /* make sure the time string are consistent among all processes */
    MPI_Bcast(str_att, strlen(str_att), MPI_CHAR, 0, MPI_COMM_WORLD);

    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "history", strlen(str_att),
                             &str_att[0]);
    ERR

    /* define dimensions x and y */
    err = ncmpi_def_dim(ncid, "Y", global_ny, &dimid[0]);
    ERR
    err = ncmpi_def_dim(ncid, "X", global_nx, &dimid[1]);
    ERR

    /* define a 2D variable of integer type */
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid);
    ERR

    /* add attributes to the variable */
    strcpy(str_att, "example attribute of type text.");
    err = ncmpi_put_att_text(ncid, varid, "str_att_name", strlen(str_att),
                             &str_att[0]);
    ERR

    for (i=0; i<8; i++) float_att[i] = i;
    err = ncmpi_put_att_float(ncid, varid, "float_att_name", NC_FLOAT, 8,
                              &float_att[0]);
    ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid);
    ERR

    /* now we are in data mode */
    start[0] = 0;
    start[1] = NX * rank;
    count[0] = NY;
    count[1] = NX;

    err = ncmpi_put_vara_int_all(ncid, varid, start, count, &buf[0][0]);
    ERR

    err = ncmpi_close(ncid);
    ERR

    MPI_Finalize();
    return 0;
}

