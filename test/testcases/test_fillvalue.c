/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests
 * 1. if PnetCDF allows put attribute _FillValue to global variable.
 * 2. if PnetCDF can return the right error code NC_EBADTYPE when put attribute
 *    _FillValue to a non-global variable with a different data type to the
 *    variable's.
 *
 * Expected results from running command ncdump on the output file:
 *
 * % ncdump testfile.nc
 * netcdf testfile {
 * variables:
 * 	int var ;
 * 		var:_FillValue = 5678 ;
 *
 * // global attributes:
 *		:_FillValue = 1.234f ;
 * data:
 *	
 *  var = _ ;
 * }
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

int main(int argc, char **argv) {
    char filename[256];
    int  err, nerrs=0, ncid, cmode, varid, int_buf, rank, nprocs;
    float flt_buf;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
        sprintf(cmd_str, "*** TESTING C   %s for _FillValue for NC_GLOBAL ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* create a file */
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    flt_buf = 1.234;
    err = ncmpi_put_att(ncid, NC_GLOBAL, "_FillValue", NC_FLOAT, 1, &flt_buf);
    CHECK_ERR

    err = ncmpi_def_var(ncid, "var", NC_INT, 0, NULL, &varid);
    CHECK_ERR

    err = ncmpi_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, &flt_buf);
    EXP_ERR(NC_EBADTYPE)

    int_buf = 5678;
    err = ncmpi_put_att(ncid, varid, "_FillValue", NC_INT, 1, &int_buf);
    CHECK_ERR

    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

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
