/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Similar to flexible.c, this program tests APIs with a need of type conversion.
 *
 * This program tests PnetCDF flexible APIs, ncmpi_put_vara_all(),
 * ncmpi_iput_vara() to write two 2D array variables (one is of 4-byte
 * integer byte and the other float type) in parallel. It then uses flexible
 * get/iget APIs to read data back and check the contents. The program first
 * defines 2 netCDF variables of sizes
 *    var_zy: NZ*nprocs x NY
 *    var_yx: NY x NX*nprocs
 *
 * The data partitioning patterns on the 2 variables are row-wise and
 * column-wise, respectively. Each process writes a subarray of size
 * NZ x NY and NY x NX to var_zy and var_yx, respectively.
 * Both local buffers have a ghost cell of length 3 surrounded along each
 * dimension.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -O2 -o flexible2 flexible2.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./flexible2 /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            Z = 20 ;
 *            Y = 5 ;
 *            X = 20 ;
 *    variables:
 *            int var_zy(Z, Y) ;
 *            float var_yx(Y, X) ;
 *    data:
 *
 *     var_zy =
 *      10, 10, 10, 10, 10,
 *      10, 10, 10, 10, 10,
 *      10, 10, 10, 10, 10,
 *      10, 10, 10, 10, 10,
 *      10, 10, 10, 10, 10,
 *      11, 11, 11, 11, 11,
 *      11, 11, 11, 11, 11,
 *      11, 11, 11, 11, 11,
 *      11, 11, 11, 11, 11,
 *      11, 11, 11, 11, 11,
 *      12, 12, 12, 12, 12,
 *      12, 12, 12, 12, 12,
 *      12, 12, 12, 12, 12,
 *      12, 12, 12, 12, 12,
 *      12, 12, 12, 12, 12,
 *      13, 13, 13, 13, 13,
 *      13, 13, 13, 13, 13,
 *      13, 13, 13, 13, 13,
 *      13, 13, 13, 13, 13,
 *      13, 13, 13, 13, 13 ;
 *
 *     var_yx =
 *      10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13,
 *      10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13,
 *      10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13,
 *      10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13,
 *      10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13 ;
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

#define NZ 5
#define NY 5
#define NX 70

/*----< tst_io() >-----------------------------------------------------------*/
static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, j, rank, nprocs, err, nerrs=0, req, status, ghost_len=3;
    int ncid, varid0, varid1, dimid[3], *buf_zy, verbose=0;
    int array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    double *buf_yx;
    MPI_Offset start[2], count[2];
    MPI_Datatype subarray;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define 3 dimensions */
    err = ncmpi_def_dim(ncid, "Z", NZ*nprocs, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NY,        &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs, &dimid[2]); CHECK_ERR

    /* define a variable of size (NZ * nprocs) * NY */
    err = ncmpi_def_var(ncid, "var_zy", NC_INT,   2, &dimid[0], &varid0); CHECK_ERR
    /* define a variable of size NY * (NX * nprocs) */
    err = ncmpi_def_var(ncid, "var_yx", NC_FLOAT, 2, &dimid[1], &varid1); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* var_zy is partitioned along Z dimension */
    array_of_sizes[0]    = NZ + 2*ghost_len;
    array_of_sizes[1]    = NY + 2*ghost_len;
    array_of_subsizes[0] = NZ;
    array_of_subsizes[1] = NY;
    array_of_starts[0]   = ghost_len;
    array_of_starts[1]   = ghost_len;
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C, MPI_INT, &subarray);
    MPI_Type_commit(&subarray);

    if (verbose && rank == 0) {
        printf("ghost_len = %d\n", ghost_len);
        printf("local array size    = %d x %d\n", array_of_sizes[0], array_of_sizes[1]);
        printf("local array subsize = %d x %d\n", array_of_subsizes[0], array_of_subsizes[1]);
        printf("local array start   = %d x %d\n", array_of_starts[0], array_of_starts[1]);
        printf("local array end     = %d x %d\n", array_of_starts[0]+array_of_subsizes[0], array_of_starts[1]+array_of_subsizes[1]);
    }

    int buffer_len = (NZ+2*ghost_len) * (NY+2*ghost_len);
    buf_zy = (int*) malloc(sizeof(int) * buffer_len);
    for (i=0; i<buffer_len; i++) buf_zy[i] = rank+10;

    start[0] = NZ * rank; start[1] = 0;
    count[0] = NZ;        count[1] = NY;
    /* calling a blocking flexible API */
    if (coll_io)
        err = ncmpi_put_vara_all(ncid, varid0, start, count, buf_zy, 1, subarray);
    else
        err = ncmpi_put_vara(ncid, varid0, start, count, buf_zy, 1, subarray);
    CHECK_ERR

    /* check the contents of put buffer */
    for (i=0; i<buffer_len; i++) {
        if (buf_zy[i] != rank+10) {
            printf("Error at line %d in %s: put buffer[%d] is altered\n",__LINE__,__FILE__,i);
            nerrs++;
            goto err_out;
        }
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    for (i=0; i<buffer_len; i++) buf_zy[i] = -1;
    /* calling a blocking flexible API */
    if (coll_io)
        err = ncmpi_get_vara_all(ncid, varid0, start, count, buf_zy, 1, subarray);
    else
        err = ncmpi_get_vara(ncid, varid0, start, count, buf_zy, 1, subarray);
    CHECK_ERR

    /* check the contents of get buffer */
    for (i=0; i<array_of_sizes[0]; i++) {
        for (j=0; j<array_of_sizes[1]; j++) {
            int index = i*array_of_sizes[1] + j;
            if (i < ghost_len || ghost_len+array_of_subsizes[0] <= i ||
                j < ghost_len || ghost_len+array_of_subsizes[1] <= j) {
                if (buf_zy[index] != -1) {
                    printf("Unexpected get buffer[%d][%d]=%d\n",
                           i,j,buf_zy[index]);
                    nerrs++;
                    goto err_out;
                }
            }
            else {
                if (buf_zy[index] != rank+10) {
                    printf("Unexpected get buffer[%d][%d]=%d\n",
                           i,j,buf_zy[index]);
                    nerrs++;
                    goto err_out;
                }
            }
        }
    }
    free(buf_zy);
    MPI_Type_free(&subarray);

    /* var_yx is partitioned along X dimension */
    array_of_sizes[0]    = NY + 2*ghost_len;
    array_of_sizes[1]    = NX + 2*ghost_len;
    array_of_subsizes[0] = NY;
    array_of_subsizes[1] = NX;
    array_of_starts[0]   = ghost_len;
    array_of_starts[1]   = ghost_len;
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C, MPI_DOUBLE,
                             &subarray);
    MPI_Type_commit(&subarray);

    buffer_len = (NY+2*ghost_len) * (NX+2*ghost_len);
    buf_yx = (double*) malloc(sizeof(double) * buffer_len);
    for (i=0; i<buffer_len; i++) buf_yx[i] = rank+10;

    start[0] = 0;  start[1] = NX * rank;
    count[0] = NY; count[1] = NX;

    /* calling a non-blocking flexible API */
    err = ncmpi_iput_vara(ncid, varid1, start, count, buf_yx, 1, subarray,&req);
    CHECK_ERR
    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, &status);
    else
        err = ncmpi_wait(ncid, 1, &req, &status);
    CHECK_ERR
    err = status; CHECK_ERR

    /* check the contents of put buffer */
    for (i=0; i<buffer_len; i++) {
        if (buf_yx[i] != rank+10) {
            printf("Error at line %d in %s: iput buffer[%d]=%f is altered\n",__LINE__,__FILE__,i,buf_yx[i]);
            nerrs++;
            goto err_out;
        }
    }

    for (i=0; i<buffer_len; i++) buf_yx[i] = -1;

    /* calling a non-blocking flexible API */
    err = ncmpi_iget_vara(ncid, varid1, start, count, buf_yx, 1, subarray,&req);
    CHECK_ERR
    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, &status);
    else
        err = ncmpi_wait(ncid, 1, &req, &status);
    CHECK_ERR
    err = status; CHECK_ERR

    /* check the contents of iget buffer */
    for (i=0; i<array_of_sizes[0]; i++) {
        for (j=0; j<array_of_sizes[1]; j++) {
            double exp;
            int index = i*array_of_sizes[1] + j;
            if (i < ghost_len || ghost_len+array_of_subsizes[0] <= i ||
                j < ghost_len || ghost_len+array_of_subsizes[1] <= j)
                exp = -1;
            else
                exp = rank+10;
            if (buf_yx[index] != exp) {
                printf("Error at %d: expect buffer[%d][%d]=%.1f but got %.1f\n",
                       __LINE__,i,j,exp,buf_yx[index]);
                nerrs++;
                goto err_out;
            }
        }
    }
    free(buf_yx);
    MPI_Type_free(&subarray);

    err = ncmpi_close(ncid); CHECK_ERR

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
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "flexible API + type conversion", opt, test_io);

    MPI_Finalize();

    return err;
}
