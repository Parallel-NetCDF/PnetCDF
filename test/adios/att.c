/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program verify attribute read capability of adios driver
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h> /* setenv() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include "pnetcdf.h"
#include <math.h>

#include <testutils.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "attributes.bp"

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main(int argc, char** argv) {
    int nerrs=0, rank, nprocs, err;
    int ncid, vid, natt, int_attr;
    char filename[256], data[1024];
    double dbl_attr;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        nerrs++;
        goto fn_exit;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, FILE_NAME);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str,
        "*** TESTING C   %s for adios attribute read",
        basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    err = ncmpi_inq_natts(ncid, &natt); CHECK_ERR
    if (natt != 4){
        printf("Rank %d: Expect global atts = %d, but got %d\n", rank, 4, natt);
        nerrs++;
    }
    err = ncmpi_get_att_int(ncid, NC_GLOBAL, "temperature/number of levels",
                            &int_attr); CHECK_ERR
    if (int_attr != 1){
        printf("Rank %d: Expect global att 0 = %d, but got %d\n", rank, 1,
               int_attr);
        nerrs++;
    }
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, "temperature/description",  data);
    CHECK_ERR
    if (strcmp(data, "Global array written from 'size' processes") != 0){
        printf("Rank %d: Expect global att 1 = %s, but got %s\n", rank,
                "Global array written from 'size' processes", data);
        nerrs++;
    }
    err = ncmpi_get_att_double(ncid, NC_GLOBAL, "temperature/mean value",
                               &dbl_attr); CHECK_ERR
    if (dbl_attr != 4.5){
        printf("Rank %d: Expect global att 2 = %lf, but got %lf\n", rank, 4.5,
               dbl_attr);
        nerrs++;
    }
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, "temperature/date of coding",
                                data); CHECK_ERR
    if (strcmp(data, "Nov, 2009") != 0){
        printf("Rank %d: Expect global att 3 = %s, but got %s\n", rank,
                "Nov, 2009", data);
        nerrs++;
    }

    err = ncmpi_inq_varid(ncid, "temperature", &vid); CHECK_ERR

    err = ncmpi_inq_varnatts(ncid, vid, &natt); CHECK_ERR
    if (natt != 4){
        printf("Rank %d: Expect var %d atts = %d, but got %d\n", rank, vid, 4,
                natt);
        nerrs++;
    }
    err = ncmpi_get_att_int(ncid, vid, "number of levels", &int_attr);
    CHECK_ERR
    if (int_attr != 1){
        printf("Rank %d: Expect var %d att 0 = %d, but got %d\n", rank, vid, 1,
               int_attr);
        nerrs++;
    }
    err = ncmpi_get_att_text(ncid, vid, "description",  data); CHECK_ERR
    if (strcmp(data, "Global array written from 'size' processes") != 0){
        printf("Rank %d: Expect var %d att 1 = %s, but got %s\n", rank, vid,
                "Global array written from 'size' processes", data);
        nerrs++;
    }
    err = ncmpi_get_att_double(ncid, vid, "mean value", &dbl_attr);
    CHECK_ERR
    if (dbl_attr != 4.5){
        printf("Rank %d: Expect var %d att 2 = %lf, but got %lf\n", rank, vid,
                4.5, dbl_attr);
        nerrs++;
    }
    err = ncmpi_get_att_text(ncid, vid, "date of coding",  data); CHECK_ERR
    if (strcmp(data, "Nov, 2009") != 0){
        printf("Rank %d: Expect var %d att 3 = %s, but got %s\n", rank, vid,
                "Nov, 2009", data);
        nerrs++;
    }

    ncmpi_close(ncid);

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

