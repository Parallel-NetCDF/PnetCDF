/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example tests using a single call of ncmpi_put_varn_int_all() to
 * write a sequence of requests with arbitrary array indices and lengths.
 * Specifically, the fileview in each process is made a contiguous chunk.
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicc -O2 -o varn_contig varn_contig.c -lpnetcdf
 *    % mpiexec -n 4 ./varn_contig /pvfs2/wkliao/testfile.nc
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *             Y = 4 ;
 *             X = 10 ;
 *    variables:
 *             int var(Y, X) ;
 *    data:
 *
 *     var =
 *       200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
 *       300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
 *       400, 400, 400, 400, 400, 400, 400, 400, 400, 400,
 *       100, 100, 100, 100, 100, 100, 100, 100, 100, 100 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), memset() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 4
#define NX 10
#define NDIMS 2

static
int check_contents_for_fail(int *buffer)
{
    int i, nprocs;
    int expected[NY*NX] = {200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
                           300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
                           400, 400, 400, 400, 400, 400, 400, 400, 400, 400,
                           100, 100, 100, 100, 100, 100, 100, 100, 100, 100};

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* check if the contents of buf are expected */
    for (i=0; i<NY*NX; i++) {
        if (expected[i] >= nprocs) continue;
        if (buffer[i] != expected[i]) {
            printf("Expected read buf[%d]=%d, but got %d\n",
                   i,expected[i],buffer[i]);
            return 1;
        }
    }
    return 0;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, rank, nprocs, err, nerrs=0;
    int ncid, varid[3], dimid[2], num_reqs, *buffer, *r_buffer;
    MPI_Offset w_len, **starts=NULL, **counts=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#ifdef DEBUG
    if (nprocs != 4 && rank == 0)
        printf("Warning: %s is intended to run on 4 processes\n",argv[0]);
#endif

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* create a global array of size NY * NX */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, NDIMS, dimid, &varid[0]);
    CHECK_ERR
    if (nprocs < 4) { /* need 4 processes to fill the variables */
        err = ncmpi_set_fill(ncid, NC_FILL, NULL);
        CHECK_ERR
    }
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    /* make sure write fill requests sync-ed to the file before testing */
    err = ncmpi_sync(ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* pick arbitrary numbers of requests for 4 processes */
    num_reqs = 0;
    if (rank == 0)      num_reqs = 4;
    else if (rank == 1) num_reqs = 5;
    else if (rank == 2) num_reqs = 5;
    else if (rank == 3) num_reqs = 5;

    if (num_reqs > 0) {
        starts    = (MPI_Offset**) malloc(sizeof(MPI_Offset*) * num_reqs);
        counts    = (MPI_Offset**) malloc(sizeof(MPI_Offset*) * num_reqs);
        starts[0] = (MPI_Offset*)  calloc(num_reqs * NDIMS, sizeof(MPI_Offset));
        counts[0] = (MPI_Offset*)  calloc(num_reqs * NDIMS, sizeof(MPI_Offset));
        for (i=1; i<num_reqs; i++) {
            starts[i] = starts[i-1] + NDIMS;
            counts[i] = counts[i-1] + NDIMS;
        }
    }

    /* assign arbitrary starts and counts */
    const int y=0, x=1;
    if (rank == 0) {
        starts[0][y] = 3; starts[0][x] = 0; counts[0][y] = 1; counts[0][x] = 3;
        starts[1][y] = 3; starts[1][x] = 4; counts[1][y] = 1; counts[1][x] = 3;
        starts[2][y] = 3; starts[2][x] = 3; counts[2][y] = 1; counts[2][x] = 1;
        starts[3][y] = 3; starts[3][x] = 7; counts[3][y] = 1; counts[3][x] = 3;
        /*                  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                           100, 100, 100, 100, 100, 100, 100, 100, 100, 100
             req id:        0    0    0    2    1    1    1    3    3    3
         */
    } else if (rank ==1) {
        starts[0][y] = 0; starts[0][x] = 3; counts[0][y] = 1; counts[0][x] = 2;
        starts[1][y] = 0; starts[1][x] = 8; counts[1][y] = 1; counts[1][x] = 2;
        starts[2][y] = 0; starts[2][x] = 5; counts[2][y] = 1; counts[2][x] = 2;
        starts[3][y] = 0; starts[3][x] = 7; counts[3][y] = 1; counts[3][x] = 1;
        starts[4][y] = 0; starts[4][x] = 0; counts[4][y] = 1; counts[4][x] = 3;
        /*                 200, 200, 200, 200, 200, 200, 200, 200, 200, 200
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
             req id:        4    4    4    0    0    2    2    3    1    1
         */
    } else if (rank ==2) {
        starts[0][y] = 1; starts[0][x] = 1; counts[0][y] = 1; counts[0][x] = 3;
        starts[1][y] = 1; starts[1][x] = 0; counts[1][y] = 1; counts[1][x] = 1;
        starts[2][y] = 1; starts[2][x] = 5; counts[2][y] = 1; counts[2][x] = 2;
        starts[3][y] = 1; starts[3][x] = 7; counts[3][y] = 1; counts[3][x] = 3;
        starts[4][y] = 1; starts[4][x] = 4; counts[4][y] = 1; counts[4][x] = 1;
        /*                  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                           300, 300, 300, 300, 300, 300, 300, 300, 300, 300
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
             req id:        1    0    0    0    4    2    3    3    3    3
         */
    } else if (rank ==3) {
        starts[0][y] = 2; starts[0][x] = 3; counts[0][y] = 1; counts[0][x] = 3;
        starts[1][y] = 2; starts[1][x] = 6; counts[1][y] = 1; counts[1][x] = 2;
        starts[2][y] = 2; starts[2][x] = 0; counts[2][y] = 1; counts[2][x] = 2;
        starts[3][y] = 2; starts[3][x] = 8; counts[3][y] = 1; counts[3][x] = 2;
        starts[4][y] = 2; starts[4][x] = 2; counts[4][y] = 1; counts[4][x] = 1;
        /*                  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
                           400, 400, 400, 400, 400, 400, 400, 400, 400, 400
                            - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,  - ,
             req id:        2    2    4    0    0    0    1    1    3    3
         */
    }

    w_len = NX; /* total write length for this process */

    /* allocate I/O buffer and initialize its contents */
    r_buffer = (int*) malloc(sizeof(int) * NY*NX);
    buffer   = (int*) malloc(sizeof(int) * w_len);
    for (i=0; i<w_len; i++) buffer[i] = rank*100 + 100;

    /* check error code: NC_ENULLSTART */
    if (coll_io)
        err = ncmpi_put_varn_int_all(ncid, varid[0], 1, NULL, NULL, NULL);
    else
        err = ncmpi_put_varn_int(ncid, varid[0], 1, NULL, NULL, NULL);
    if (err != NC_ENULLSTART) {
        printf("expecting error code NC_ENULLSTART but got %s\n",ncmpi_strerrno(err));
        nerrs++;
    }

    /* write using varn API */
    if (coll_io)
        err = ncmpi_put_varn_int_all(ncid, varid[0], num_reqs, starts, counts, buffer);
    else
        err = ncmpi_put_varn_int(ncid, varid[0], num_reqs, starts, counts, buffer);
    CHECK_ERR

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR

    if (nprocs > 4) MPI_Barrier(MPI_COMM_WORLD);

    /* read back and check contents */
    memset(r_buffer, 0, NY*NX*sizeof(int));
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, varid[0], r_buffer);
    else
        err = ncmpi_get_var_int(ncid, varid[0], r_buffer);
    CHECK_ERR
    nerrs += check_contents_for_fail(r_buffer);

    err = ncmpi_close(ncid);
    CHECK_ERR

    free(buffer);
    free(r_buffer);
    if (num_reqs > 0) {
        free(starts[0]);
        free(counts[0]);
        free(starts);
        free(counts);
    }

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

    err = tst_main(argc, argv, "put_varn with contig fileview", opt, test_io);

    MPI_Finalize();

    return err;
}
