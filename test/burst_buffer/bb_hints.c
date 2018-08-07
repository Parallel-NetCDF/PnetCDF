/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program adds two new variables to an existing netCDF file.
 * It is used to test if PnetCDF can correctly calculate the file offsets
 * for the two new variables, in particular for files that align the
 * fixed-size variables to a boundary larger than 4 bytes, for instance
 * a file created by PnetCDF with default alignment of 512 bytes.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o add_var add_var.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 add_var testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

int main(int argc, char** argv) {
    char filename[256];
    int rank, nprocs, err, flag, nerrs=0;
    int log_enabled;
    int ncid;
    MPI_Info info, infoused;
    char hint[MPI_MAX_INFO_VAL];

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
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for checking offsets of new variables ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_burst_buf_del_on_close", "disable");
    MPI_Info_set(info, "nc_burst_buf_flush_buffer_size", "256");
    /* MPI_Info_set(info, "nc_burst_buf_dirname", "()@^$@!(_&$)@(#%%&)(*#$"); */

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid); CHECK_ERR
    err = ncmpi_inq_file_info(ncid, &infoused); CHECK_ERR

    MPI_Info_get(infoused, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, hint, &flag);
    if (flag && strcasecmp(hint, "enable") == 0)
        log_enabled = 1;
    else
        log_enabled = 0;

    if (log_enabled) {
        MPI_Info_get(infoused, "nc_burst_buf_del_on_close", MPI_MAX_INFO_VAL - 1, hint, &flag);
        if (flag) {
            if (strcmp(hint, "disable") != 0) {
                printf("Error at line %d: unexpected nc_burst_buf_del_on_close = %s, but got %s\n", __LINE__, "disable", hint);
                nerrs++;
            }
        }
        else{
            printf("Error at line %d: nc_burst_buf_del_on_close is not set\n", __LINE__);
            nerrs++;
        }

        MPI_Info_get(infoused, "nc_burst_buf_flush_buffer_size", MPI_MAX_INFO_VAL - 1, hint, &flag);
        if (flag) {
            if (strcmp(hint, "256") != 0) {
                printf("Error at line %d: unexpected nc_burst_buf_flush_buffer_size = %s, but got %s\n", __LINE__, "256", hint);
                nerrs++;
            }
        }
        else{
            printf("Error at line %d: nc_burst_buf_flush_buffer_size is not set\n", __LINE__);
            nerrs++;
        }
    }

    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    MPI_Info_free(&info);
    MPI_Info_free(&infoused);

    /* check if PnetCDF freed all internal malloc */
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

