/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program opens corrupted NetCDF files and check whether the correct
 * error codes are returned.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h> /* setenv() */
#include <string.h> /* strlen() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: input file %s (%s)\n", \
        __LINE__,__FILE__,filename,ncmpi_strerrno(err)); \
    } \
}

#define EXP_ERR(exp) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expected_errno %s but got %s\n", \
        __LINE__,__FILE__,ncmpi_strerrno(exp), ncmpi_strerrno(err)); \
    } \
}

static int
str2err(char *err_str, int *expected_errno)
{
         if (!strcmp(err_str, "NC_ENOTNC"))    *expected_errno = NC_ENOTNC;
    else if (!strcmp(err_str, "NC_EMAXVARS"))  *expected_errno = NC_EMAXVARS;
    else if (!strcmp(err_str, "NC_EUNLIMIT"))  *expected_errno = NC_EUNLIMIT;
    else if (!strcmp(err_str, "NC_ENULLPAD"))  *expected_errno = NC_ENULLPAD;
    else if (!strcmp(err_str, "NC_EVARSIZE"))  *expected_errno = NC_EVARSIZE;
    else if (!strcmp(err_str, "NC_ENOTBUILT")) *expected_errno = NC_ENOTBUILT;
    else if (!strcmp(err_str, "NC_NOERR"))     *expected_errno = NC_NOERR;
    else return NC_EINVAL;
    return NC_NOERR;
}

int main(int argc, char** argv) {
    char filename[256];
    int nerrs=0, rank, err, ncid, expected_errno=NC_NOERR;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 3) {
        if (!rank) printf("Usage: %s [filename] expected_errno\n",argv[0]);
        goto fn_exit;
    }
    snprintf(filename, 256, "%s", argv[1]);
    err = str2err(argv[2], &expected_errno); CHECK_ERR

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** %s detecting corrupted file %s",
                basename(argv[0]), basename(argv[1]));
        printf("%-66s --- ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
#ifndef ENABLE_NULL_BYTE_HEADER_PADDING
    if (expected_errno == NC_ENULLPAD) CHECK_ERR
    else
#endif
    EXP_ERR(expected_errno)

    if (err == NC_NOERR || err == NC_ENULLPAD) {
        err = ncmpi_close(ncid); CHECK_ERR
    }

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR && malloc_size > 0) /* this test is for running 1 process */
        printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
               malloc_size);

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf("fail with %d mismatches\n",nerrs);
        else       printf("pass\n");
    }

    MPI_Finalize();
    return (nerrs > 0);
}

