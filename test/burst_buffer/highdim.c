/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example tests writing a high dimensional int variable when log io is
 * enabled.  Each process own a submatrix of size 2 * 2 * 2 * ... 1 * 1 * 1
 * that have as much size 2 dimensions as the buffer can accommodate.
 * All cell in the submatrix is it's rank + 1
 * The submatrix is combined by interleaving along the first dimension.
 * The variable will have size (2 * np) * 2 * 2 * ... 1 * 1 * 1
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file when dimension is set to 2.
 *
 *    % mpicc -O2 -o log_higndim.c -lpnetcdf
 *    % mpiexec -n 4 ./log_higndim [testfile.nc]
 *    % ncmpidump [testfile.nc]
 *    netcdf test {
 *    // file format: CDF-1
 *    dimensions:
 *            D0 = 8 ;
 *            D1 = 2 ;
 *    variables:
 *            int M(D0, D1) ;
 *    data:
 *
 *      M =
 *
 *       1, 1,
 *       2, 2,
 *       3, 3,
 *       4, 4,
 *       1, 1,
 *       2, 2,
 *       3, 3,
 *       4, 4 ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <limits.h>
#include <testutils.h>
#include <libgen.h>

#if NC_MAX_DIMS < 1024
#define DIM NC_MAX_DIMS
#else
#define DIM 1024
#endif

#define BSIZE 1024 * 1024

static
int test_bb(const char *out_path,
            int         coll_io,
            MPI_Info    info)
{
    char *folder, *dup_out_path, dimname[64];
    int i, err=NC_NOERR, nerrs=0, rank, np, ndims, ncid, varid;
    int *dimid=NULL, *buffer=NULL;
    long long j;
    MPI_Offset *start=NULL, *count=NULL, *stride=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    ndims = DIM;

    /* Allocate buffers */
    dimid = (int*)malloc(sizeof(int) * ndims);
    start = (MPI_Offset*)malloc(sizeof(MPI_Offset) * ndims);
    count = (MPI_Offset*)malloc(sizeof(MPI_Offset) * ndims);
    stride = (MPI_Offset*)malloc(sizeof(MPI_Offset) * ndims);
    buffer = (int*)malloc(sizeof(int) * BSIZE);
    if (dimid == NULL || start == NULL || count == NULL || stride == NULL || buffer == NULL) {
        printf("Error at line %d in %s: malloc error\n", __LINE__, __FILE__);
        nerrs++;
        goto err_out;
    }

    /* Initialize buffers and calculate share among processes
     * Each process writes a hyper rectangle to the variable
     * The variable is formed by stacking the rectangle of every processes along first dimension
     */
    for (i = 0; i < BSIZE; i++) {
        buffer[i] = rank + 1;
    }
    for (i = 0, j = 2; i < ndims; i++) {
        start[i] = 0;
        stride[i] = 1;
        /* Most dimensions must be 1 for high dimensional variable
         * Set dimension size to 2 until we run out of buffer
         */
        if (j < BSIZE) {
            count[i] = 2;
            j <<= 1;
        }
        else{
            count[i] = 1;
        }
    }
    start[0] = rank;
    stride[0] = np;

    MPI_Info_set(info, "nc_burst_buf", "enable");
    MPI_Info_set(info, "nc_burst_buf_overwrite", "enable");

    dup_out_path = strdup(out_path);
    folder = dirname(dup_out_path);
    if (folder == NULL)
        MPI_Info_set(info, "nc_burst_buf_dirname", ".");
    else
        MPI_Info_set(info, "nc_burst_buf_dirname", folder);
    free(dup_out_path);

    /* Create new netcdf file */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* Define dimensions */
    for (i = 0; i < ndims; i++) {
        sprintf(dimname, "D%d", i);
        /* Submatrix of each process stack along the first dimension */
        if (i == 0)
            err = ncmpi_def_dim(ncid, dimname, count[i] * np, dimid + i);
        else
            err = ncmpi_def_dim(ncid, dimname, count[i], dimid + i);
        CHECK_ERR
    }

    /* Define variable */
    err = ncmpi_def_var(ncid, "M", NC_INT, ndims, dimid, &varid);
    CHECK_ERR

    /* Switch to data mode */
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* Write variable */
    if (coll_io)
        err = ncmpi_put_vars_int_all(ncid, varid, start, count, stride, buffer);
    else
        err = ncmpi_put_vars_int(ncid, varid, start, count, stride, buffer);
    CHECK_ERR

    /* Read it back */
    if (coll_io)
        err = ncmpi_get_vars_int_all(ncid, varid, start, count, stride, buffer);
    else
        err = ncmpi_get_vars_int(ncid, varid, start, count, stride, buffer);
    CHECK_ERR

    /* Verify the result */
    for (i = 0; i < BSIZE; i++) {
        if (buffer[i] != rank + 1) {
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n",
                    __LINE__, __FILE__, i, rank + 1, buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

    /* Close the file */
    err = ncmpi_close(ncid);
    CHECK_ERR

    if (start != NULL) free(start);
    if (count != NULL) free(count);
    if (stride != NULL) free(stride);
    if (dimid != NULL) free(dimid);
    if (buffer != NULL) free(buffer);

err_out:
    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int err=NC_NOERR;
    MPI_Info local_info;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    MPI_Info_dup(info, &local_info);
    err = test_bb(out_path, coll_io, local_info);
    MPI_Info_free(&local_info);

    MPI_Info_dup(info, &local_info);
    MPI_Info_set(local_info, "nc_burst_buf_shared_logs", "enable");
    err = test_bb(out_path, coll_io, local_info);
    MPI_Info_free(&local_info);

    return err;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "high dimensional variables", opt, test_io);

    MPI_Finalize();

    return err;
}
