/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests the special case of ONLY one record variable is defined
 * and the record size is not aligned with the 4-byte boundary. As defined in
 * CDF-1 and CDF-2 format specifications:
 *    "A special case: Where there is exactly one record variable, we drop the
 *    requirement that each record be four-byte aligned, so in this case there
 *    is no record padding."
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

#define STR_LEN 19
#define NUM_VALS 2

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char filename[256];
    int i, err, nerrs=0, rank, nprocs, cmode;
    int ncid, dimids[2], varid;
    char data[NUM_VALS][STR_LEN + 1], data_in[NUM_VALS*STR_LEN];
    MPI_Offset start[2];
    MPI_Offset count[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for only one record variable ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    strcpy(data[0], "2005-04-11_12:00:00"); /* 19 bytes not a multiply of 4 */
    strcpy(data[1], "2005-04-11_13:00:00");

    cmode = NC_CLOBBER;
    err  = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, dimids); CHECK_ERR
    err = ncmpi_def_dim(ncid, "text_dim", STR_LEN, &dimids[1]); CHECK_ERR

    /* create ONLY one record variable of type NC_CHAR and make sure each
     * record is of size not aligned with 4-byte boundary.
     */
    err = ncmpi_def_var(ncid, "text_var", NC_CHAR, 2, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* Write some records of var data. */
    count[0] = 1;
    count[1] = STR_LEN;
    start[0] = 0;
    start[1] = 0;
    for (i=0; i<NUM_VALS; i++) {
        err = ncmpi_put_vara_text_all(ncid, varid, start, count, data[start[0]]);
        CHECK_ERR
        start[0]++;
    }

    err = ncmpi_close(ncid); CHECK_ERR

    err  = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "text_var", &varid); CHECK_ERR

    /* read the entire data back */
    err = ncmpi_get_var_text_all(ncid, varid, data_in); CHECK_ERR

    /* check the contents */
    for (i=0; i<NUM_VALS; i++)
      if (strncmp(data[i], data_in+i*STR_LEN, STR_LEN)) {
          printf("Error at line %d in %s: expecting %s but got %s\n",
          __LINE__,__FILE__,data[i],data_in+i*STR_LEN);
          nerrs++;
      }

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

