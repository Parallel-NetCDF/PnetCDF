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
        int exp = rank+10; \
        if (buf[i] != exp) { \
            printf("Error at line %d : put buf[%d]=%d alter from %d\n", \
                   __LINE__, i, buf[i], exp); \
            nerrs++; \
            goto err_out; \
        } \
    }

#define NDIMS 2
#define NY 4
#define NX 2

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int err, nerrs=0, ncid, dimid[NDIMS], varid[5], ndims=NDIMS;
    int i, j, k, nprocs, rank, req, *buf=NULL;
    MPI_Offset start[NDIMS] = {0};
    MPI_Offset count[NDIMS] = {0};
    MPI_Offset stride[NDIMS] = {0};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL); CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", nprocs*NX, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v0", NC_INT, ndims, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v1", NC_INT, ndims, dimid, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v2", NC_INT, ndims, dimid, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v3", NC_INT, ndims, dimid, &varid[3]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v4", NC_INT, ndims, dimid, &varid[4]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid); CHECK_ERR
    }

    buf = (int*) malloc(sizeof(int) * NY * NX);
    for (i=0; i<NY*NX; i++) buf[i] = rank+10;

    start[0] = 0;
    start[1] = rank*NX;
    count[0] = NY;
    count[1] = NX;

    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf);
    else
        err = ncmpi_put_vara_int(ncid, varid[0], start, count, buf);
    CHECK_ERR
    CHECK_PUT_BUF

    if (coll_io)
        err = ncmpi_put_vars_int_all(ncid, varid[1], start, count, NULL, buf);
    else
        err = ncmpi_put_vars_int(ncid, varid[1], start, count, NULL, buf);
    CHECK_ERR
    CHECK_PUT_BUF

    start[0] = 0;
    start[1] = rank;
    count[0] = NY;
    count[1] = NX;
    stride[0] = 1;
    stride[1] = nprocs;

    if (coll_io)
        err = ncmpi_put_vars_int_all(ncid, varid[2], start, count, stride, buf);
    else
        err = ncmpi_put_vars_int(ncid, varid[2], start, count, stride, buf);
    CHECK_ERR
    CHECK_PUT_BUF

    if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
        /* test bput_vars */
        err = ncmpi_buffer_attach(ncid, NY*NX*sizeof(int)); CHECK_ERR
    }

    start[0] = 0;
    start[1] = rank*NX;
    count[0] = NY;
    count[1] = NX;
    if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC) {
        if (coll_io)
            err = ncmpi_put_vars_int_all(ncid, varid[3], start, count, NULL, buf);
        else
            err = ncmpi_put_vars_int(ncid, varid[3], start, count, NULL, buf);
        CHECK_ERR
    }
    else {
        err = ncmpi_bput_vars_int(ncid, varid[3], start, count, NULL, buf, &req);
        CHECK_ERR
        if (coll_io)
            err = ncmpi_wait_all(ncid, 1, &req, NULL);
        else
            err = ncmpi_wait(ncid, 1, &req, NULL);
        CHECK_ERR
    }
    CHECK_PUT_BUF

    start[0] = 0;
    start[1] = rank;
    count[0] = NY;
    count[1] = NX;
    stride[0] = 1;
    stride[1] = nprocs;
    if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC) {
        if (coll_io)
            err = ncmpi_put_vars_int_all(ncid, varid[4], start, count, stride, buf);
        else
            err = ncmpi_put_vars_int(ncid, varid[4], start, count, stride, buf);
        CHECK_ERR
    }
    else {
        err = ncmpi_bput_vars_int(ncid, varid[4], start, count, stride, buf, &req);
        CHECK_ERR
        if (coll_io)
            err = ncmpi_wait_all(ncid, 1, &req, NULL);
        else
            err = ncmpi_wait(ncid, 1, &req, NULL);
        CHECK_ERR
    }
    CHECK_PUT_BUF
    free(buf);
    buf = NULL;

    if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
        err = ncmpi_buffer_detach(ncid); CHECK_ERR
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    buf = (int*) malloc(sizeof(int) * NY * NX * nprocs);
    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, varid[0], buf);
    else
        err = ncmpi_get_var_int(ncid, varid[0], buf);
    CHECK_ERR

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
                int idx = i*nprocs*NX+j*NX+k;
                if (buf[idx] != j+10) {
                    printf("Error at line %d: expect v0[%d]=%d but got %d\n",
                           __LINE__, idx, j+10, buf[idx]);
                    nerrs++;
                    goto err_out;
                }
            }
        }
    }

    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, varid[1], buf);
    else
        err = ncmpi_get_var_int(ncid, varid[1], buf);
    CHECK_ERR

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
                int idx = i*nprocs*NX+j*NX+k;
                if (buf[idx] != j+10) {
                    printf("Error at line %d: expect v1[%d]=%d but got %d\n",
                           __LINE__,idx, j+10, buf[idx]);
                    nerrs++;
                    goto err_out;
                }
            }
        }
    }

    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, varid[2], buf);
    else
        err = ncmpi_get_var_int(ncid, varid[2], buf);
    CHECK_ERR

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
                int idx = i*nprocs*NX+k*nprocs+j;
                if (buf[idx] != j+10) {
                    printf("Error at line %d: expect v2[%d]=%d but got %d\n",
                           __LINE__,idx, j+10, buf[idx]);
                    nerrs++;
                    goto err_out;
                }
            }
        }
    }

    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, varid[3], buf);
    else
        err = ncmpi_get_var_int(ncid, varid[3], buf);
    CHECK_ERR

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
                int idx = i*nprocs*NX+j*NX+k;
                if (buf[idx] != j+10) {
                    printf("Error at line %d: expect v3[%d]=%d but got %d\n",
                           __LINE__,idx, j+10, buf[idx]);
                    nerrs++;
                    goto err_out;
                }
            }
        }
    }

    memset(buf, 0, (size_t)NY * NX * nprocs * sizeof(int));
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, varid[4], buf);
    else
        err = ncmpi_get_var_int(ncid, varid[4], buf);
    CHECK_ERR

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
                int idx = i*nprocs*NX+k*nprocs+j;
                if (buf[idx] != j+10) {
                    printf("Error at line %d: expect v4[%d]=%d but got %d\n",
                           __LINE__,idx, j+10, buf[idx]);
                    nerrs++;
                    goto err_out;
                }
            }
        }
    }

err_out:
    err = ncmpi_close(ncid); CHECK_ERR
    if (buf != NULL) free(buf);

    return (nerrs > 0);
}

int main(int argc, char **argv) {

    int err;
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 2;    /* enable and disable intra-node aggregation */
    opt.drv      = 2;    /* test GIO and MPI-IO driver */
    opt.bb       = 2;    /* enable and disable burst-buffering feature */
    opt.mod      = 2;    /* collective and independent data mode */
    opt.hdr_diff = true; /* run ncmpidiff for file header */
    opt.var_diff = true; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "NULL stride", opt, test_io);

    MPI_Finalize();

    return err;
}

