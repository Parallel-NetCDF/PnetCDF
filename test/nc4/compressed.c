/*
   This program test reading a compressed netcdf 4 file
   The data is written using NetCDF and read by PnetCDF
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define FNAME "gzip_example.nc"

int main(int argc, char **argv) {
    int err, nerrs=0, rank, np;
    int ncid, ndims, nvars;
    char *dir_name=".", filename[512];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [dir_name]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) dir_name = argv[1];

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for reading compressed file", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    sprintf(filename, "%s/%s", dir_name, FNAME);

    /* Read with PnetCDF */
    /* Open the file */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    err = ncmpi_inq(ncid, &ndims, &nvars, NULL, NULL);
    CHECK_ERR

    /* Close file */
    ncmpi_close(ncid);

    /* check if there is any malloc residue */
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
