/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests if PnetCDF can return the right error code NC_EEXIST
 * when create mode NC_NOCLOBBER is used and the file exists.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <pnetcdf.h>

#define ERR if (err!=NC_NOERR) {printf("Error at line %d: %s\n", __LINE__,ncmpi_strerror(err)); exit(-1);}

int main(int argc, char **argv) {
    char filename[128];
    int  err, pass, ncid, cmode, rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    strcpy(filename, "testfile.nc");
    if (argc == 2) strcpy(filename, argv[1]);
    MPI_Bcast(filename, 128, MPI_CHAR, 0, MPI_COMM_WORLD);

    /* create a file if it does not exist */
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    ERR
    err = ncmpi_close(ncid); ERR

    /* now the file exists, test if PnetCDF can return correct error code */
    cmode = NC_NOCLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    if (err == NC_EEXIST) /* err == NC_EOFILE */
        pass = 1;
    else
        pass = 0;

    MPI_Allreduce(MPI_IN_PLACE, &pass, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    char cmd_str[80];
    sprintf(cmd_str, "*** TESTING C   %s for NC_NOCLOBBER and NC_EEXIST ", argv[0]);
    if (rank == 0) {
        if (pass) printf("%-66s ------ pass\n", cmd_str);
        else      printf("%-66s ------ failed\n", cmd_str);
    }

    MPI_Finalize();
    return 0;
}
