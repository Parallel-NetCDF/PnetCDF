/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests PnetCDF flexible varm APIs, i.e. ncmpi_put_varm_all(),
 * ncmpi_get_varm_all(), and their nonblocking versions, to write a 2D array
 * double variable of size NY x NX*nproc in parallel.
 *
 * The data partitioning patterns on the variable is column-wise.
 * The local buffer has ghost cells surrounded along both dimensions.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -O2 -o flexible_varm flexible_varm.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./flexible_varm /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            Y = 6 ;
 *            X = 16 ;
 *    variables:
 *            double var(Y, X) ;
 *    data:
 *
 *    var =
 *      0,  6, 12, 18, 0,  6, 12, 18, 0,  6, 12, 18, 0,  6, 12, 18,
 *      1,  7, 13, 19, 1,  7, 13, 19, 1,  7, 13, 19, 1,  7, 13, 19,
 *      2,  8, 14, 20, 2,  8, 14, 20, 2,  8, 14, 20, 2,  8, 14, 20,
 *      3,  9, 15, 21, 3,  9, 15, 21, 3,  9, 15, 21, 3,  9, 15, 21,
 *      4, 10, 16, 22, 4, 10, 16, 22, 4, 10, 16, 22, 4, 10, 16, 22,
 *      5, 11, 17, 23, 5, 11, 17, 23, 5, 11, 17, 23, 5, 11, 17, 23 ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 32
#define NX 128
#define GHOST 2

#define INIT_PUT_BUF \
    for (i=0; i<array_of_sizes[0]; i++) { \
        for (j=0; j<array_of_sizes[1]; j++) { \
            if (i < GHOST || GHOST+array_of_subsizes[0] <= i || \
                j < GHOST || GHOST+array_of_subsizes[1] <= j) \
                buf[i][j] = -1; \
            else \
                buf[i][j] = (i-GHOST)*array_of_subsizes[1]+(j-GHOST); \
        } \
    }

#define CHECK_PUT_BUF \
    for (i=0; i<array_of_sizes[0]; i++) { \
        for (j=0; j<array_of_sizes[1]; j++) { \
            if (i < GHOST || GHOST+array_of_subsizes[0] <= i || \
                j < GHOST || GHOST+array_of_subsizes[1] <= j) { \
                if (buf[i][j] != -1) { \
                    printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%d\n", \
                           __LINE__,__FILE__,i,j,buf[i][j]); \
                    nerrs++; \
                    goto err_out; \
                } \
            } \
            else { \
                if (buf[i][j] != (i-GHOST)*array_of_subsizes[1]+(j-GHOST)) { \
                    printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%d\n", \
                           __LINE__,__FILE__,i,j,buf[i][j]); \
                    nerrs++; \
                    goto err_out; \
                } \
            } \
        } \
    }

#define INIT_GET_BUF \
    for (i=0; i<array_of_sizes[0]; i++) \
        for (j=0; j<array_of_sizes[1]; j++) \
            buf[i][j] = -2;

#define CHECK_GET_BUF \
    for (i=0; i<array_of_sizes[0]; i++) { \
        for (j=0; j<array_of_sizes[1]; j++) { \
            if (i < GHOST || GHOST+array_of_subsizes[0] <= i || \
                j < GHOST || GHOST+array_of_subsizes[1] <= j) { \
                if (buf[i][j] != -2) { \
                    printf("Unexpected get buffer[%d][%d]=%d\n", \
                           i,j,buf[i][j]); \
                    nerrs++; \
                    goto err_out; \
                } \
            } \
            else { \
                if (buf[i][j] != (i-GHOST)*array_of_subsizes[1]+(j-GHOST)) { \
                    printf("Unexpected get buffer[%d][%d]=%d\n", \
                           i,j,buf[i][j]); \
                    nerrs++; \
                    goto err_out; \
                } \
            } \
        } \
    }

#define INDEP_MODE 0
#define COLL_MODE 1

static
int tst_varm(const char *out_path,
             const char *in_path, /* ignored */
             int         format,
             int         coll_io,
             MPI_Info    info)
{
    int i, j, rank, nprocs, err, nerrs=0, req, status;
    int ncid, varid, dimid[2];
    int array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    int buf[NX+2*GHOST][NY+2*GHOST];
    MPI_Offset start[2], count[2], stride[2], imap[2];
    MPI_Datatype  subarray;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define 2 dimensions */
    err = ncmpi_def_dim(ncid, "Y", NY,        &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs, &dimid[1]); CHECK_ERR

    /* define a variable of size NY * (NX * nprocs) */
    err = ncmpi_def_var(ncid, "var", NC_DOUBLE, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

     start[0] = 0;  start[1] = NX * rank;
     count[0] = NY; count[1] = NX;
    stride[0] = 1; stride[1] = 1;
      imap[0] = 1;   imap[1] = NY; /* would be {NX, 1} if not transposing */

    /* var is partitioned along X dimension in a matrix transported way */
    array_of_sizes[0]    = NX + 2*GHOST;
    array_of_sizes[1]    = NY + 2*GHOST;
    array_of_subsizes[0] = NX;
    array_of_subsizes[1] = NY;
    array_of_starts[0]   = GHOST;
    array_of_starts[1]   = GHOST;
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C, MPI_INT, &subarray);
    MPI_Type_commit(&subarray);

    /* calling a blocking put_varm flexible API -----------------------------*/
    /* initiate put buffer contents */
    INIT_PUT_BUF
    if (coll_io)
        err = ncmpi_put_varm_all(ncid, varid, start, count, stride, imap, buf,
                                 1, subarray);
    else
        err = ncmpi_put_varm(ncid, varid, start, count, stride, imap, buf,
                             1, subarray);
    CHECK_ERR

    /* check the contents of put buffer */
    CHECK_PUT_BUF

    /* calling a nonblocking put_varm flexible API --------------------------*/
    /* initiate put buffer contents */
    INIT_PUT_BUF
    err = ncmpi_iput_varm(ncid, varid, start, count, stride, imap, buf,
                          1, subarray, &req);
    CHECK_ERR

    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, &status);
    else
        err = ncmpi_wait(ncid, 1, &req, &status);
    CHECK_ERR
    err = status; CHECK_ERR

    /* check the contents of put buffer */
    CHECK_PUT_BUF

   /* When running in independent data mode, flushing writes is necessary
    * before reading the data back.
    */
    MPI_Barrier(MPI_COMM_WORLD);
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    /* read back using a blocking get_varm flexible API ---------------------*/
    /* initiate get buffer contents */
    INIT_GET_BUF

    /* calling a blocking flexible API */
    if (coll_io)
        err = ncmpi_get_varm_all(ncid, varid, start, count, stride, imap, buf,
                                 1, subarray);
    else
        err = ncmpi_get_varm(ncid, varid, start, count, stride, imap, buf,
                             1, subarray);
    CHECK_ERR

    /* check the contents of get buffer */
    CHECK_GET_BUF

    /* read back using a non-blocking flexible API --------------------------*/
    /* initiate get buffer contents */
    INIT_GET_BUF

    /* calling a blocking flexible API */
    err = ncmpi_iget_varm(ncid, varid, start, count, stride, imap, buf,
                          1, subarray, &req);
    CHECK_ERR
    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, &status);
    else
        err = ncmpi_wait(ncid, 1, &req, &status);
    CHECK_ERR
    err = status; CHECK_ERR

    /* check the contents of get buffer */
    CHECK_GET_BUF

    MPI_Type_free(&subarray);

    err = ncmpi_close(ncid); CHECK_ERR

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
    int nerrs=0;

    nerrs = tst_varm(out_path, in_path, format, coll_io, info);
    if (nerrs > 0) goto err_out;

    /* disable PnetCDF internal buffering */
    MPI_Info_set(info, "nc_ibuf_size", "0");

    nerrs = tst_varm(out_path, in_path, format, coll_io, info);
    if (nerrs > 0) goto err_out;

err_out:
    return nerrs;
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
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "flexible varm APIs", opt, test_io);

    MPI_Finalize();

    return err;
}
