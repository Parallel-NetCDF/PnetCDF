/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
   This program tests ncmpi_create(), ncmpi_open(), and ncmpi_close().
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */

#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

int main(int argc, char **argv)
{
    char filename[512];
    int i, err, nprocs, rank, nerrs=0, ncid;
    int format[3] = {0, NC_64BIT_OFFSET, NC_64BIT_DATA};
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 512, "%s", argv[1]);
    else           sprintf(filename, "%s.nc", argv[0]);

    if (rank == 0) {
        char *cmd_str = (char *)malloc(strlen(argv[0]) + 256);
            sprintf(cmd_str, "*** TESTING C   %s for file create", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_header_align_size", "100");
    MPI_Info_set(info, "nc_var_align_size", "200");
    MPI_Info_set(info, "nc_record_align_size", "300");

    for (i=0; i<3; i++) {
        /* Create a new file */
        int cmode = NC_CLOBBER | format[i];
        err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
        CHECK_ERR

        /* Close the file. */
        err = ncmpi_close(ncid);
        CHECK_ERR

        /* Open the file */
        err = ncmpi_open(MPI_COMM_WORLD, filename, NC_WRITE, info, &ncid);
        CHECK_ERR

        /* Close the file. */
        err = ncmpi_close(ncid);
        CHECK_ERR
    }

    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    /* check if there is any malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has "OFFFMT" bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs)
            printf(FAIL_STR, nerrs);
        else
            printf(PASS_STR);
    }

    MPI_Finalize();

    return (nerrs > 0);
}
