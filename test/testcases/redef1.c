/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests redefine mode.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o redef1 redef1.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 redef1 testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

static int
tst_fmt(char *filename, int cmode)
{
    int i, j, k, rank, ncid, err, nerrs=0;
    int dim0id, dim1id, dim5id, dim9id, dim2id, dimsid[2], dims2id[2];
    int varid, var3id, var4id, var2id;
    int *data;
    double *dbl_data;
    MPI_Offset len0=10, len1=3, len5=5, len9=9, len2=10;
    MPI_Offset start[2], count[2];
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);

    /* Test NetCDF 4 first as ncvalidator checks only classic formats */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "dim0", len0, &dim0id); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", len1, &dim1id); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim5", len5, &dim5id); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim9", len9, &dim9id); CHECK_ERR

    dimsid[0] = dim0id;
    dimsid[1] = dim1id;
    err = ncmpi_def_var(ncid, "xyz", NC_INT, 2, dimsid, &varid); CHECK_ERR

    dimsid[0] = dim0id;
    dimsid[1] = dim5id;
    err = ncmpi_def_var(ncid, "connect", NC_INT, 2, dimsid, &var3id); CHECK_ERR

    dimsid[0] = dim0id;
    dimsid[1] = dim9id;
    err = ncmpi_def_var(ncid, "connect_exterior", NC_INT, 2, dimsid, &var4id); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* put data */
    start[0] = 0;
    start[1] = 0;
    count[0] = len0;
    count[1] = len1;

    data = (int*) malloc(len0*len1 * sizeof(int));
    k=0;
    for (i=0; i<len0; i++)
        for (j=0; j<len1; j++)
            data[i*len1+j] = k++;
    if (rank > 0) count[0] = count[1] = 0;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, &data[0]);
    CHECK_ERR
    free(data);

    count[0] = len0;
    count[1] = len5;
    data = (int*) malloc(len0*len5 * sizeof(int));
    k=0;
    for (i=0; i<len0; i++)
        for (j=0; j<len5; j++)
            data[i*len5+j] = k++;
    if (rank > 0) count[0] = count[1] = 0;
    err = ncmpi_put_vara_int_all(ncid, var3id, start, count, &data[0]);
    CHECK_ERR
    free(data);

    count[0] = len0;
    count[1] = len9;
    data = (int*) malloc(len0*len9 * sizeof(int));
    k=0;
    for (i=0; i<len0; i++)
        for (j=0; j<len9; j++)
            data[i*len9+j] = k++;
    if (rank > 0) count[0] = count[1] = 0;
    err = ncmpi_put_vara_int_all(ncid, var4id, start, count, &data[0]);
    CHECK_ERR
    free(data);

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(comm, filename, NC_WRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    err = ncmpi_redef(ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "dim2", len2, &dim2id); CHECK_ERR

    dims2id[0] = dim0id;
    dims2id[1] = dim2id;
    err = ncmpi_def_var(ncid, "xyz_r", NC_DOUBLE, 2, dims2id, &var2id);
    CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    start[0] = 0;
    start[1] = 0;
    count[0] = len0;
    count[1] = len2;
    k=0;
    dbl_data = (double*) malloc(len0*len2 * sizeof(double));
    for (i=0; i<len0; i++)
        for (j=0; j<len2; j++) {
            dbl_data[i*len2+j] = (k*k);
            k++;
        }
    if (rank > 0) count[0] = count[1] = 0;
    err = ncmpi_put_vara_double_all(ncid, var2id, start, count, &dbl_data[0]);
    CHECK_ERR
    free(dbl_data);

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char** argv)
{
    char filename[256], *hint_value;
    int commsize, rank, err, nerrs=0, bb_enabled=0;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &commsize);
    MPI_Comm_rank(comm, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "redef2.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for entering re-define mode ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

#ifdef DEBUG
    if (commsize > 1 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n",argv[0]);
#endif

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
