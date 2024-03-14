/*********************************************************************
 *
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program is to test
 *
 * large header size, i.e. > INT_MAX, i.e. 2 GiB
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <limits.h> /* INT_MAX */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

int main(int argc, char** argv)
{
    char filename[256];
    int rank, nprocs, err, nerrs=0;
    int ncid, cmode, dimid, varid, buf;
    MPI_Offset extent, start;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
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
        sprintf(cmd_str, "*** TESTING C   %s for large header ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* define a dimension of size = nprocs */
    err = ncmpi_def_dim(ncid, "dim", nprocs, &dimid);
    CHECK_ERR

    /* define a variable */
    err = ncmpi_def_var(ncid, "var0", NC_INT, 1, &dimid, &varid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, &dimid, &varid);
    CHECK_ERR

    /* make file header extent > 4 GiB */
    extent = (MPI_Offset)INT_MAX + 1024;
    err = ncmpi__enddef(ncid, 0, extent, 0, 0);
    CHECK_ERR

    /* write to the variable */
    start = rank;
    buf = rank;
    err = ncmpi_put_var1_int_all(ncid, varid, &start, &buf);
    CHECK_ERR

    err = ncmpi_close(ncid);
    CHECK_ERR

    if (err != NC_NOERR) goto err_out;

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_ERR

    /* inquire ID of the variable */
    err = ncmpi_inq_varid(ncid, "var1", &varid);
    CHECK_ERR

    /* read from the variable */
    buf = -1;
    err = ncmpi_get_var1_int_all(ncid, varid, &start, &buf);
    CHECK_ERR

    if (buf != rank) {
        nerrs++;
        printf("Error at line %d in %s: expecting read buf %d but got %d\n",
               __LINE__,__FILE__,rank,buf);
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0) {
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
            ncmpi_inq_malloc_list();
        }
    }

err_out:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

