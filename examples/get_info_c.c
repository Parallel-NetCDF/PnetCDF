/*********************************************************************
 *
 *  Copyright (C) 2010, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/*
 *    prints all MPI-IO hints used
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>


#define HANDLE_ERROR {                                \
    if (err != NC_NOERR)                              \
        printf("Error at line %d (%s)\n", __LINE__,   \
               ncmpi_strerror(err));                  \
}

/*----< print_info() >------------------------------------------------------*/
void print_info(MPI_Info *info_used)
{
    int  i, err, flag, nkeys;

    MPI_Info_get_nkeys(*info_used, &nkeys);
    printf("MPI File Info: nkeys = %d\n",nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(*info_used, i, key);
        MPI_Info_get_valuelen(*info_used, key, &valuelen, &flag);
        MPI_Info_get(*info_used, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %24s, value = %s\n",i,key,value);
    }
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int rank, ncid, err;
    MPI_Info info_used;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 2) {
        if (!rank) printf("Usage: %s filename\n",argv[0]);
        MPI_Finalize();
        return 1;
    }

    /* create the file */
    err = ncmpi_create(MPI_COMM_WORLD, argv[1], NC_CLOBBER|NC_64BIT_DATA,
                       MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        printf("Error: ncmpi_create() file %s (%s)\n",argv[1],ncmpi_strerror(err));
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(1);
    }

    /* exit the define mode */
    err = ncmpi_enddef(ncid);
    HANDLE_ERROR

    /* get all the hints used */
    err = ncmpi_get_file_info(ncid, &info_used);
    HANDLE_ERROR

    /* close the file */
    err = ncmpi_close(ncid);
    HANDLE_ERROR

    if (rank == 0) print_info(&info_used);
    MPI_Info_free(&info_used);

    MPI_Finalize();
    return 0;
}

