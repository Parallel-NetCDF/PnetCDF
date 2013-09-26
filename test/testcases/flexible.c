/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests the use of flexible API.
 * The write buffer is a 2D array of size NY x NX
 * The MPI data type for the buffer is defined by swapping the 1st and 2nd
 * rows of the array. It uses MPI_Type_create_hindex()
 *
 * The expected reults from the output file contents are:
 * (when running on 1 MPI process)
 *
 *  % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 * 	   Y = UNLIMITED ; // (2 currently)
 * 	   X = 5 ;
 *    variables:
 * 	   int var(Y, X) ;
 *    data:
 * 
 *    var =
 *      1, 1, 1, 1, 1,
 *      0, 0, 0, 0, 0 ;
 *    }
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <pnetcdf.h>

#define NY 2
#define NX 5
#define ERR if (err!=NC_NOERR) {printf("Error at line %d: %s\n", __LINE__,ncmpi_strerror(err)); exit(-1);}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {

    char         filename[128];
    int          i, j, err, ncid, varid, dimids[2], pass;
    int          rank, nprocs, blocklengths[2], buf[NY][NX], *bufptr;
    MPI_Offset   start[2], count[2];
    MPI_Aint     a0, a1, disps[2];
    MPI_Datatype buftype;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    strcpy(filename, "testfile.nc");
    if (rank == 0 && argc > 1) strcpy(filename, argv[1]);
    MPI_Bcast(filename, 128, MPI_CHAR, 0, MPI_COMM_WORLD);

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL,
                       &ncid); ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimids[0]); ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs,    &dimids[1]); ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimids, &varid); ERR
    err = ncmpi_enddef(ncid); ERR

    /* initialize the contents of the array */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = j;

    /* construct an MPI derived data type for swapping 1st row with 2nd row */
    blocklengths[0] = blocklengths[1] = NX;
    MPI_Get_address(buf[1], &a0);
    MPI_Get_address(buf[0], &a1);
    disps[0] = 0;
    disps[1] = a1 - a0;
    bufptr = buf[1];
    err = MPI_Type_create_hindexed(2, blocklengths, disps, MPI_INT, &buftype);
    if (err != MPI_SUCCESS) printf("MPI error MPI_Type_create_hindexed\n");
    MPI_Type_commit(&buftype);

    start[0] = 0; start[1] = NX*rank;
    count[0] = 2; count[1] = NX;

    /* call flexible API */
    err = ncmpi_put_vara_all(ncid, varid, start, count, bufptr, 1, buftype); ERR
    MPI_Type_free(&buftype);

    /* check if the contents of buf are altered */
    for (j=0; j<NY; j++)
        for (i=0; i<NX; i++)
            if (buf[j][i] != j)
                printf("buf[%d][%d] != %d\n",j,i,j,buf[j][i]);
 
    err = ncmpi_close(ncid); ERR

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); ERR

    err = ncmpi_inq_varid(ncid, "var", &varid); ERR

    /* initialize the contents of the array to a different value */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = -1;

    /* read back variable */
    start[0] = 0; start[1] = NX*rank;
    count[0] = 2; count[1] = NX;
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf[0]); ERR

    err = ncmpi_close(ncid); ERR

    /* check if the contents of buf are expected */
    pass = 1;
    for (j=0; j<2; j++) {
        int val = (j == 0) ? 1 : 0;
        for (i=0; i<NX; i++)
            if (buf[j][i] != val) {
                printf("Unexpected buf[%d][%d] != %d\n",j,i,val,buf[j][i]);
                pass = 0;
            }
    }

    MPI_Allreduce(MPI_IN_PLACE, &pass, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    char cmd_str[80];
    sprintf(cmd_str, "*** TESTING %s for using ncmpi_put_vara_all() ", argv[0]);
    if (rank == 0) {
        if (pass) printf("%-66s ------ pass\n", cmd_str);
        else      printf("%-66s ------ failed\n", cmd_str);
    }

    MPI_Finalize();
    return 0;
}
