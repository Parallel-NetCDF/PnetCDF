/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if the number of records is updated correctly. It first
 * writes to the 2nd record of a variable, followed by writing to the 1st
 * record. After the 2nd write, a call to ncmpi_inq_dimlen() or ncmpi_inq_dim()
 * should report seeing 2 records. Then, it does a similar test under
 * independent data mode.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o record record.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 record record.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}

int main(int argc, char** argv) {
    char *filename="testfile.nc";
    int pass=1, rank, nprocs, err;
    int ncid, cmode, varid, dimid, buf;
    MPI_Offset start, count, length;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        goto fn_exit;
    }
    if (argc == 2) filename = argv[1];

    if (rank >= 1) goto fn_exit; /* this test is for running 1 process */

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_SELF, filename, cmode, info, &ncid); ERR

    /* define dimension and variable */
    err = ncmpi_def_dim(ncid, "DIM", NC_UNLIMITED, &dimid); ERR
    err = ncmpi_def_var(ncid, "REC_VAR", NC_INT, 1, &dimid, &varid); ERR
    err = ncmpi_enddef(ncid); ERR

    /* write the 2nd record first */
    start = 1; count = 1;
    err = ncmpi_put_vara_int_all(ncid, varid, &start, &count, &buf); ERR

    /* write the 1st record now */
    start = 0; count = 1;
    err = ncmpi_put_vara_int_all(ncid, varid, &start, &count, &buf); ERR

    err = ncmpi_inq_dimlen(ncid, dimid, &length); ERR
    if (length != 2) {
        printf("Error: expecting 2 records, but got %lld record(s)\n",length);
        pass = 0;
    }

    if (pass) { /* test independent data mode */
        err = ncmpi_begin_indep_data(ncid); ERR
        /* write the 4th record */
        start = 3; count = 1;
        err = ncmpi_put_vara_int(ncid, varid, &start, &count, &buf); ERR

        /* write the 3rd record */
        start = 2; count = 1;
        err = ncmpi_put_vara_int(ncid, varid, &start, &count, &buf); ERR

        err = ncmpi_inq_dimlen(ncid, dimid, &length); ERR
        if (length != 4) {
            printf("Error: expecting 4 records, but got %lld record(s)\n",
                   length);
            pass = 0;
        }
        err = ncmpi_end_indep_data(ncid); ERR
    }

    err = ncmpi_close(ncid); ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR && malloc_size > 0) /* this test is for running 1 process */
        printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
               malloc_size);

    char cmd_str[80];
    sprintf(cmd_str, "*** TESTING C   %s for write records in reversed order", argv[0]);
    if (rank == 0) {
        if (pass) printf("%-66s ------ pass\n", cmd_str);
        else      printf("%-66s ------ failed\n", cmd_str);
    }

fn_exit:
    MPI_Finalize();
    return 0;
}

