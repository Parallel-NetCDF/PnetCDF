/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests whether the file striping size and count retrieved from
 * MPI-IO hints are consistent among all MPI processes.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o get_striping get_striping.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 get_striping testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

static int
tst_fmt(char *filename, int cmode)
{
    int err, nerrs=0, ncid, fmt;
    int striping_size=0, striping_count=0, root_striping_size, root_striping_count;

    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_format(ncid, &fmt); CHECK_ERR
    err = ncmpi_inq_striping(ncid, &striping_size, &striping_count);
    if (fmt == NC_FORMAT_NETCDF4 || fmt == NC_FORMAT_NETCDF4_CLASSIC)
        EXP_ERR(NC_ENOTSUPPORT)
    else
        CHECK_ERR

    root_striping_size  = striping_size;
    root_striping_count = striping_count;
    err = MPI_Bcast(&root_striping_size,  1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_ERR(err)
    err = MPI_Bcast(&root_striping_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_ERR(err)
    if (root_striping_size != striping_size) {
        printf("Error at line %d in %s: inconsistent striping_size (root=%d local=%d)\n",
               __LINE__,__FILE__, root_striping_size, striping_size);
        nerrs++;
    }
    if (root_striping_count != striping_count) {
        printf("Error at line %d in %s: inconsistent striping_count (root=%d local=%d)\n",
               __LINE__,__FILE__, root_striping_count, striping_count);
        nerrs++;
    }

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char** argv) {
    char filename[256], *hint_value;
    int rank, err, nerrs=0, bb_enabled=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for striping info ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* check whether burst buffering is enabled */
    if (inq_env_hint("nc_burst_buf", &hint_value)) {
        if (strcasecmp(hint_value, "enable") == 0) bb_enabled = 1;
        free(hint_value);
    }

    nerrs += tst_fmt(filename, 0);
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += tst_fmt(filename, NC_NETCDF4);
        nerrs += tst_fmt(filename, NC_NETCDF4 | NC_CLASSIC_MODEL);
#endif
    }
    nerrs += tst_fmt(filename, NC_64BIT_DATA);

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

