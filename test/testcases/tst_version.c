/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  Check whether PnetCDF version string returned from ncmpi_inq_libvers() 
 *  matches the constant PNETCDF_VERSION defined in header file pnetcdf.h.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy(), strtok(), strcmp() */
#include <libgen.h>  /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char *str, *pnetcdf_version_str;
    int err, nerrs=0, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for PnetCDF library version ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    str = (char*) malloc(strlen(ncmpi_inq_libvers())+1);
    strcpy(str, ncmpi_inq_libvers());
    pnetcdf_version_str = strtok(str, " ");

    if (pnetcdf_version_str == NULL) {
        printf("\nError: ncmpi_inq_libvers() returns ill form string %s\n",
               ncmpi_inq_libvers());
        nerrs++;
    }
    else if (strcmp(pnetcdf_version_str, PNETCDF_VERSION)) {
        printf("\nError: ncmpi_inq_libvers() returns %s does not match with PNETCDF_VERSION %s\n",
               pnetcdf_version_str, PNETCDF_VERSION);
        nerrs++;
    }
    free(str);

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

