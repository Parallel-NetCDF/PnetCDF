/*********************************************************************
 *
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * To compile:
 *     % mpicc mpi_subarray.c -o mpi_subarray
 *
 * To run:
 *     % mpiexec -n 4 ./mpi_subarray output_file
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define LEN 16

#define ERROR(fname) \
    if (err != MPI_SUCCESS) { \
        int errorStringLen; \
        char errorString[MPI_MAX_ERROR_STRING]; \
        MPI_Error_string(err, errorString, &errorStringLen); \
        printf("Error at line %d when calling %s: %s\n",__LINE__,fname,errorString); \
    }

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char *filename;
    int err, rank, nprocs, verbose=0, psizes[2]={0,0}, buf[LEN*LEN];
    int gsizes[2], sizes[2], starts[2];
    MPI_File     fh;
    MPI_Datatype ftype;
    MPI_Status   status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    filename = "testfile";
    if (argc > 1) filename = argv[1];

    /* Creates a division of processors in a cartesian grid */
    err = MPI_Dims_create(nprocs, 2, psizes);
    ERROR("MPI_Dims_create");

    /* set global array sizes, local array sizes, starting offsets */
    starts[0] = LEN * (rank / psizes[1]);
    starts[1] = LEN * (rank % psizes[1]);
    gsizes[0] = LEN * psizes[0];
    gsizes[1] = LEN * psizes[1];
    sizes[0]  = sizes[1] = LEN;

    if (verbose) 
        printf("rank %d: gsizes=%d %d sizes=%d %d starts=%d %d\n", rank,
               gsizes[0],gsizes[1],sizes[0],sizes[1],starts[0],starts[1]);

    /* create file type: 2D subarray */
    err = MPI_Type_create_subarray(2, gsizes, sizes, starts, MPI_ORDER_C, MPI_INT, &ftype);
    ERROR("MPI_Type_create_subarray");
    err = MPI_Type_commit(&ftype);
    ERROR("MPI_Type_commit");

    /* open the file */
    err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                        MPI_INFO_NULL, &fh);
    ERROR("MPI_File_open");

    /* set the file view */
    err = MPI_File_set_view(fh, 0, MPI_INT, ftype, "native", MPI_INFO_NULL);
    ERROR("MPI_File_set_view");
    err = MPI_Type_free(&ftype);
    ERROR("MPI_Type_free");

    /* MPI collective write */
    err = MPI_File_write_all(fh, buf, LEN*LEN, MPI_INT, &status);
    ERROR("MPI_File_write_all");

    MPI_File_close(&fh);

    MPI_Finalize();
    return 0;
}
