/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset() */
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

/* check if user put buffer contents altered */
#define CHECK_PUT_BUF \
    for (i=0; i<NY*NX; i++) { \
        if (buf[i] != rank+10) { \
            printf("Error at line %d in %s: user put buffer[%d] altered from %d to %d\n", \
                   __LINE__,__FILE__, i, rank+10, buf[i]); \
            nerrs++; \
        } \
    }

#define NDIMS 2
#define NY 4
#define NX 2

static int
tst_fmt(char *filename, int cmode)
{
    int err, nerrs=0, ncid, dimid[NDIMS], varid[5], ndims=NDIMS;
    int i, j, k, nprocs, rank, req, *buf;
    MPI_Offset start[NDIMS] = {0};
    MPI_Offset count[NDIMS] = {0};
    MPI_Offset stride[NDIMS] = {0};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", nprocs*NX, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v0", NC_INT, ndims, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v1", NC_INT, ndims, dimid, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v2", NC_INT, ndims, dimid, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v3", NC_INT, ndims, dimid, &varid[3]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v4", NC_INT, ndims, dimid, &varid[4]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    start[0] = 0;
    start[1] = rank*NX;
    count[0] = NY;
    count[1] = NX;
    buf = (int*) malloc((size_t)NY * NX * sizeof(int));
    for (i=0; i<NY*NX; i++) buf[i] = rank+10;

    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR
    CHECK_PUT_BUF

    err = ncmpi_put_vars_int_all(ncid, varid[1], start, count, NULL, buf); CHECK_ERR
    CHECK_PUT_BUF

    start[0] = 0;
    start[1] = rank;
    count[0] = NY;
    count[1] = NX;
    stride[0] = 1;
    stride[1] = nprocs;
    err = ncmpi_put_vars_int_all(ncid, varid[2], start, count, stride, buf); CHECK_ERR
    CHECK_PUT_BUF

    if (!(cmode & NC_NETCDF4)) {
        /* test bput_vars */
        err = ncmpi_buffer_attach(ncid, NY*NX*sizeof(int)); CHECK_ERR
    }

    start[0] = 0;
    start[1] = rank*NX;
    count[0] = NY;
    count[1] = NX;
    if (cmode & NC_NETCDF4) {
        err = ncmpi_put_vars_int_all(ncid, varid[3], start, count, NULL, buf); CHECK_ERR
    }
    else {
        err = ncmpi_bput_vars_int(ncid, varid[3], start, count, NULL, buf, &req); CHECK_ERR
        err = ncmpi_wait_all(ncid, 1, &req, NULL); CHECK_ERR
    }
    CHECK_PUT_BUF

    start[0] = 0;
    start[1] = rank;
    count[0] = NY;
    count[1] = NX;
    stride[0] = 1;
    stride[1] = nprocs;
    if (cmode & NC_NETCDF4) {
        err = ncmpi_put_vars_int_all(ncid, varid[4], start, count, stride, buf); CHECK_ERR
    }
    else {
        err = ncmpi_bput_vars_int(ncid, varid[4], start, count, stride, buf, &req); CHECK_ERR
        err = ncmpi_wait_all(ncid, 1, &req, NULL); CHECK_ERR
    }
    CHECK_PUT_BUF
    free(buf);

    if (!(cmode & NC_NETCDF4)) {
        err = ncmpi_buffer_detach(ncid); CHECK_ERR
    }

    buf = (int*) malloc((size_t)NY * NX * nprocs * sizeof(int));
    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    err = ncmpi_get_var_int_all(ncid, varid[0], buf); CHECK_ERR

    /* check read buffer contents */
    /*  v0 =
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13 ;
     */
    for (i=0; i<NY; i++) {
        for (j=0; j<nprocs; j++) {
            for (k=0; k<NX; k++) {
                if (buf[i*nprocs*NX+j*NX+k] != j+10) {
                    printf("Error at line %d in %s: expected buffer[%d]=%d but got %d\n",
                           __LINE__,__FILE__,i*nprocs*NX+j*NX+k, j+10, buf[i*nprocs*NX+j*NX+k]);
                    nerrs++;
                }
            }
        }
    }

    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    err = ncmpi_get_var_int_all(ncid, varid[1], buf); CHECK_ERR

    /* check read buffer contents */
    /*  v1 =
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13 ;
     */
    for (i=0; i<NY; i++) {
        for (j=0; j<nprocs; j++) {
            for (k=0; k<NX; k++) {
                if (buf[i*nprocs*NX+j*NX+k] != j+10) {
                    printf("Error at line %d in %s: expected buffer[%d]=%d but got %d\n",
                           __LINE__,__FILE__,i*nprocs*NX+j*NX+k, j+10, buf[i*nprocs*NX+j*NX+k]);
                    nerrs++;
                }
            }
        }
    }

    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    err = ncmpi_get_var_int_all(ncid, varid[2], buf); CHECK_ERR

    /* check read buffer contents */
    /*  v2 =
     *    10, 11, 12, 13, 10, 11, 12, 13,
     *    10, 11, 12, 13, 10, 11, 12, 13,
     *    10, 11, 12, 13, 10, 11, 12, 13,
     *    10, 11, 12, 13, 10, 11, 12, 13 ;
     */
    for (i=0; i<NY; i++) {
        for (k=0; k<NX; k++) {
            for (j=0; j<nprocs; j++) {
                if (buf[i*nprocs*NX+k*nprocs+j] != j+10) {
                    printf("Error at line %d in %s: expected buffer[%d]=%d but got %d\n",
                           __LINE__,__FILE__,i*nprocs*NX+k*nprocs+j, j+10, buf[i*nprocs*NX+k*nprocs+j]);
                    nerrs++;
                }
            }
        }
    }

    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    err = ncmpi_get_var_int_all(ncid, varid[3], buf); CHECK_ERR

    /* check read buffer contents */
    /*  v3 =
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13,
     *    10, 10, 11, 11, 12, 12, 13, 13 ;
     */
    for (i=0; i<NY; i++) {
        for (j=0; j<nprocs; j++) {
            for (k=0; k<NX; k++) {
                if (buf[i*nprocs*NX+j*NX+k] != j+10) {
                    printf("Error at line %d in %s: expected buffer[%d]=%d but got %d\n",
                           __LINE__,__FILE__,i*nprocs*NX+j*NX+k, j+10, buf[i*nprocs*NX+j*NX+k]);
                    nerrs++;
                }
            }
        }
    }

    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    err = ncmpi_get_var_int_all(ncid, varid[4], buf); CHECK_ERR

    /* check read buffer contents */
    /*  v4 =
     *    10, 11, 12, 13, 10, 11, 12, 13,
     *    10, 11, 12, 13, 10, 11, 12, 13,
     *    10, 11, 12, 13, 10, 11, 12, 13,
     *    10, 11, 12, 13, 10, 11, 12, 13 ;
     */
    for (i=0; i<NY; i++) {
        for (k=0; k<NX; k++) {
            for (j=0; j<nprocs; j++) {
                if (buf[i*nprocs*NX+k*nprocs+j] != j+10) {
                    printf("Error at line %d in %s: expected buffer[%d]=%d but got %d\n",
                           __LINE__,__FILE__,i*nprocs*NX+k*nprocs+j, j+10, buf[i*nprocs*NX+k*nprocs+j]);
                    nerrs++;
                }
            }
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);

    return nerrs;
}

int main(int argc, char **argv)
{
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
        sprintf(cmd_str, "*** TESTING C   %s for NULL stride ", basename(argv[0]));
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
