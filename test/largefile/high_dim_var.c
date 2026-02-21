/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests vars APIs for variables with high dimensions.
 * In particular, this is to test fix in r3164.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

#define NVARS 2
#define NDIMS 16
#define DIMLEN 3
#define NRECS 4

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    char name[32];
    size_t nelms;
    short *buffer;
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, fvarid[NVARS], rvarid[NVARS], dimids[NDIMS], rdimids[NDIMS];
    MPI_Offset start[NDIMS], count[NDIMS], stride[NDIMS];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "rdim", NC_UNLIMITED, &rdimids[0]); CHECK_ERR
    for (i=0; i<NDIMS; i++) {
        sprintf(name, "dim%d", i);
        err = ncmpi_def_dim(ncid, name, DIMLEN, &dimids[i]); CHECK_ERR
        if (i > 0) rdimids[i] = dimids[i];
    }
    for (i=0; i<NVARS; i++) {
        sprintf(name, "fix_var%d", i);
        err = ncmpi_def_var(ncid, name, NC_SHORT, NDIMS, dimids, &fvarid[i]); CHECK_ERR
        sprintf(name, "rec_var%d", i);
        err = ncmpi_def_var(ncid, name, NC_SHORT, NDIMS,rdimids, &rvarid[i]); CHECK_ERR
    }

    /* initialize the contents of fixed-size variables */
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

#ifdef STRONGER_CONSISTENCY
    ncmpi_sync(ncid);
    MPI_Barrier(MPI_COMM_WORLD);
    ncmpi_sync(ncid);
#endif

    nelms = (NRECS > DIMLEN) ? NRECS : DIMLEN;
    for (i=1; i<NDIMS; i++) nelms *= DIMLEN;
    buffer = (short*) malloc(sizeof(short) * nelms);
    if (buffer == NULL) {
        printf("Error %s at line %d: fail to allocate buffer of size %zu\n",
               basename(__FILE__), __LINE__, nelms * sizeof(short));
        nerrs++;
        goto err_out;
    }
    for (i=0; i<nelms; i++) buffer[i] = -1;

    /* initialize the contents of record variables */
    for (i=0; i<NVARS; i++) {
        for (j=0; j<NRECS; j++) {
            err = ncmpi_fill_var_rec(ncid, rvarid[i], j);
            CHECK_ERR
        }
    }

#ifdef STRONGER_CONSISTENCY
    ncmpi_sync(ncid);
    MPI_Barrier(MPI_COMM_WORLD);
    ncmpi_sync(ncid);
#endif

    for (i=0; i<nelms; i++) buffer[i] = i % 32768;

    for (i=0; i<NDIMS; i++) {
        start[i]  = 0;
        count[i]  = 2;
        stride[i] = 2;
    }
    /* only process 0 writes */
    if (rank > 0) for (i=0; i<NDIMS; i++) count[i] = 0;

    for (i=0; i<NVARS; i++) {
        start[0] = 0;
        err = ncmpi_put_vars_short_all(ncid, fvarid[i], start, count, stride,
                                       buffer); CHECK_ERR
        start[0] = 1;
        err = ncmpi_put_vars_short_all(ncid, rvarid[i], start, count, stride,
                                       buffer); CHECK_ERR
    }

#ifdef STRONGER_CONSISTENCY
    ncmpi_sync(ncid);
    MPI_Barrier(MPI_COMM_WORLD);
    ncmpi_sync(ncid);
#endif

    /* all processes read and verify */
    if (rank > 0) for (i=0; i<NDIMS; i++) count[i]  = 2;
    for (nelms=1,i=0; i<NDIMS; i++) nelms *= count[i];

    for (i=0; i<NVARS; i++) {
        for (j=0; j<nelms; j++) buffer[j] = -2;
        start[0] = 0;
        err = ncmpi_get_vars_short_all(ncid, fvarid[i], start, count, stride,
                                       buffer); CHECK_ERR
        for (j=0; j<nelms; j++) {
            if (buffer[j] != j%32768) {
                printf("Error at line %d: expect buffer[%d][%d]=%d but got %hd\n",
                       __LINE__, i, j, j%32768, buffer[j]);
                nerrs++;
                break;
            }
        }
        for (j=0; j<nelms; j++) buffer[j] = -2;
        start[0] = 1;
        err = ncmpi_get_vars_short_all(ncid, rvarid[i], start, count, stride,
                                       buffer); CHECK_ERR
        for (j=0; j<nelms; j++) {
            if (buffer[j] != j%32768) {
                printf("Error at line %d: expect buffer[%d][%d]=%d but got %hd\n",
                       __LINE__, i, j, j%32768, buffer[j]);
                nerrs++;
                break;
            }
        }
    }
    free(buffer);

    err = ncmpi_close(ncid); CHECK_ERR

    /* check file header by just opening it */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

err_out:
    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 0; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "vars APIs on high-dim variables", opt, test_io);

    MPI_Finalize();

    return err;
}
