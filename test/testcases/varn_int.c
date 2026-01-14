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
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicc -O2 -o varn_int varn_int.c -lpnetcdf
 *    % mpiexec -n 4 ./varn_int /pvfs2/wkliao/testfile.nc
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *             Y = 4 ;
 *             X = 10 ;
 *             REC_DIM = UNLIMITED ; // (4 currently)
 *    variables:
 *             int fix_var(Y, X) ;
 *             int rec_var(REC_DIM, X) ;
 *    data:
 *
 *     fix_var =
 *       13, 13, 13, 11, 11, 10, 10, 12, 11, 11,
 *       10, 12, 12, 12, 13, 11, 11, 12, 12, 12,
 *       11, 11, 12, 13, 13, 13, 10, 10, 11, 11,
 *       10, 10, 10, 12, 11, 11, 11, 13, 13, 13 ;
 *
 *     rec_var =
 *       13, 13, 13, 11, 11, 10, 10, 12, 11, 11,
 *       10, 12, 12, 12, 13, 11, 11, 12, 12, 12,
 *       11, 11, 12, 13, 13, 13, 10, 10, 11, 11,
 *       10, 10, 10, 12, 11, 11, 11, 13, 13, 13 ;
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
int check_contents_for_fail(char *var_name, int *buffer)
{
    int i, err=0, nprocs;
    int expected[NY*NX] = {13, 13, 13, 11, 11, 10, 10, 12, 11, 11,
                           10, 12, 12, 12, 13, 11, 11, 12, 12, 12,
                           11, 11, 12, 13, 13, 13, 10, 10, 11, 11,
                           10, 10, 10, 12, 11, 11, 11, 13, 13, 13};

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* check if the contents of buf are expected */
    for (i=0; i<NY*NX; i++) {
        if (expected[i] >= nprocs+10) continue;
        if (buffer[i] != expected[i]) {
            printf("Error: var %s expect read buf[%d]=%d, but got %d\n",
                   var_name, i,expected[i],buffer[i]);
            err = 1;
            break;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    return (err > 0);
}

static
void permute(MPI_Offset a[NDIMS], MPI_Offset b[NDIMS])
{
    int i;
    MPI_Offset tmp;
    for (i=0; i<NDIMS; i++) {
        tmp = a[i]; a[i] = b[i]; b[i] = tmp;
    }
}

#define INDEP_MODE 0
#define COLL_MODE 1

static
int tst_io(const char *filename,
           int         coll_io,
           MPI_Info    info)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, varid[2], dimid[2], num_reqs=0, *buffer=NULL, *r_buffer=NULL;
    MPI_Offset w_len, **starts=NULL, **counts=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* create a global array of size NY * NX */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "fix_var", NC_INT, NDIMS, dimid, &varid[0]);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimid[0]);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "rec_var", NC_INT, NDIMS, dimid, &varid[1]);
    CHECK_ERR
    if (nprocs < 4) { /* need 4 processes to fill the variables */
        err = ncmpi_set_fill(ncid, NC_FILL, NULL);
        CHECK_ERR
    }
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    if (nprocs < 4) { /* need 4 processes to fill the variables */
        for (i=0; i<4; i++) { /* total 4 records */
            err = ncmpi_fill_var_rec(ncid, varid[1], i); CHECK_ERR
        }
    }

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* pick arbitrary numbers of requests for 4 processes */
    num_reqs = 0;
    if (rank == 0)      num_reqs = 4;
    else if (rank == 1) num_reqs = 6;
    else if (rank == 2) num_reqs = 5;
    else if (rank == 3) num_reqs = 4;

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
        starts[0][y] = 0; starts[0][x] = 5; counts[0][y] = 1; counts[0][x] = 2;
        starts[1][y] = 1; starts[1][x] = 0; counts[1][y] = 1; counts[1][x] = 1;
        starts[2][y] = 2; starts[2][x] = 6; counts[2][y] = 1; counts[2][x] = 2;
        starts[3][y] = 3; starts[3][x] = 0; counts[3][y] = 1; counts[3][x] = 3;
        /* rank 0 is writing the followings: ("-" means skip)
                  -  -  -  -  -  0  0  -  -  -
                  0  -  -  -  -  -  -  -  -  -
                  -  -  -  -  -  -  0  0  -  -
                  0  0  0  -  -  -  -  -  -  -
         */
    } else if (rank ==1) {
        starts[0][y] = 0; starts[0][x] = 3; counts[0][y] = 1; counts[0][x] = 2;
        starts[1][y] = 0; starts[1][x] = 8; counts[1][y] = 1; counts[1][x] = 2;
        starts[2][y] = 1; starts[2][x] = 5; counts[2][y] = 1; counts[2][x] = 2;
        starts[3][y] = 2; starts[3][x] = 0; counts[3][y] = 1; counts[3][x] = 2;
        starts[4][y] = 2; starts[4][x] = 8; counts[4][y] = 1; counts[4][x] = 2;
        starts[5][y] = 3; starts[5][x] = 4; counts[5][y] = 1; counts[5][x] = 3;
        /* rank 1 is writing the followings: ("-" means skip)
                  -  -  -  1  1  -  -  -  1  1
                  -  -  -  -  -  1  1  -  -  -
                  1  1  -  -  -  -  -  -  1  1
                  -  -  -  -  1  1  1  -  -  -
         */
    } else if (rank ==2) {
        starts[0][y] = 0; starts[0][x] = 7; counts[0][y] = 1; counts[0][x] = 1;
        starts[1][y] = 1; starts[1][x] = 1; counts[1][y] = 1; counts[1][x] = 3;
        starts[2][y] = 1; starts[2][x] = 7; counts[2][y] = 1; counts[2][x] = 3;
        starts[3][y] = 2; starts[3][x] = 2; counts[3][y] = 1; counts[3][x] = 1;
        starts[4][y] = 3; starts[4][x] = 3; counts[4][y] = 1; counts[4][x] = 1;
        /* rank 2 is writing the followings: ("-" means skip)
                  -  -  -  -  -  -  -  2  -  -
                  -  2  2  2  -  -  -  2  2  2
                  -  -  2  -  -  -  -  -  -  -
                  -  -  -  2  -  -  -  -  -  -
         */
    } else if (rank ==3) {
        starts[0][y] = 0; starts[0][x] = 0; counts[0][y] = 1; counts[0][x] = 3;
        starts[1][y] = 1; starts[1][x] = 4; counts[1][y] = 1; counts[1][x] = 1;
        starts[2][y] = 2; starts[2][x] = 3; counts[2][y] = 1; counts[2][x] = 3;
        starts[3][y] = 3; starts[3][x] = 7; counts[3][y] = 1; counts[3][x] = 3;
        /* rank 3 is writing the followings: ("-" means skip)
                  3  3  3  -  -  -  -  -  -  -
                  -  -  -  -  3  -  -  -  -  -
                  -  -  -  3  3  3  -  -  -  -
                  -  -  -  -  -  -  -  3  3  3
         */
    }

    w_len = 0; /* total write length for this process */
    for (i=0; i<num_reqs; i++) {
        MPI_Offset w_req_len=1;
        for (j=0; j<NDIMS; j++)
            w_req_len *= counts[i][j];
        w_len += w_req_len;
    }

    /* allocate I/O buffer and initialize its contents */
    r_buffer = (int*) malloc(sizeof(int) * NY*NX);
    buffer   = (int*) malloc(sizeof(int) * w_len);
    for (i=0; i<w_len; i++) buffer[i] = rank+10;

    /* check error code: NC_ENULLSTART */
    if (coll_io)
        err = ncmpi_put_varn_int_all(ncid, varid[0], 1, NULL, NULL, NULL);
    else
        err = ncmpi_put_varn_int(ncid, varid[0], 1, NULL, NULL, NULL);
    EXP_ERR(NC_ENULLSTART)

    /* write using varn API */
    if (coll_io)
        err = ncmpi_put_varn_int_all(ncid, varid[0], num_reqs, starts, counts, buffer);
    else
        err = ncmpi_put_varn_int(ncid, varid[0], num_reqs, starts, counts, buffer);
    CHECK_ERR

    /* check if user put buffer contents altered */
    for (i=0; i<w_len; i++) {
        if (buffer[i] != rank+10) {
            printf("Error at line %d in %s: user put buffer[%d] altered from %d to %d\n",
                   __LINE__,__FILE__,i, rank+10, buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

    if (nprocs > 4 || !coll_io) {
        /* This program was designed to run on 4 processes. If running on more
         * than 4, then we need to sync writes before reading, especially for
         * processes of rank >= 4. Similarly, when running in independent data
         * mode, flushing writes is necessary before reading the data back.
         */
        MPI_Barrier(MPI_COMM_WORLD);
        err = ncmpi_sync(ncid);
        CHECK_ERR
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* read back and check contents */
    memset(r_buffer, 0, NY*NX*sizeof(int));
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, varid[0], r_buffer);
    else
        err = ncmpi_get_var_int(ncid, varid[0], r_buffer);
    CHECK_ERR
    nerrs += check_contents_for_fail("fix_var", r_buffer);
    if (nerrs > 0) goto err_out;

    /* permute write order */
    if (num_reqs > 0) {
        permute(starts[1], starts[2]); permute(counts[1], counts[2]);
        permute(starts[2], starts[3]); permute(counts[2], counts[3]);
    }

    /* write using varn API */
    if (coll_io)
        err = ncmpi_put_varn_int_all(ncid, varid[1], num_reqs, starts, counts, buffer);
    else
        err = ncmpi_put_varn_int(ncid, varid[1], num_reqs, starts, counts, buffer);
    CHECK_ERR

    /* check if user put buffer contents altered */
    for (i=0; i<w_len; i++) {
        if (buffer[i] != rank+10) {
            printf("Error at line %d in %s: user put buffer[%d] altered from %d to %d\n",
                   __LINE__,__FILE__,i, rank+10, buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

    if (nprocs > 4 || !coll_io) {
        /* This program was designed to run on 4 processes. If running on more
         * than 4, then we need to sync writes before reading, especially for
         * processes of rank >= 4. Similarly, when running in independent data
         * mode, flushing writes is necessary before reading the data back.
         */
        MPI_Barrier(MPI_COMM_WORLD);
        err = ncmpi_sync(ncid);
        CHECK_ERR
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* read back using get_var API and check contents */
    memset(r_buffer, 0, NY*NX*sizeof(int));
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, varid[1], r_buffer);
    else
        err = ncmpi_get_var_int(ncid, varid[1], r_buffer);
    CHECK_ERR
    nerrs += check_contents_for_fail("rec_var", r_buffer);
    if (nerrs > 0) goto err_out;

    /* read back using get_varn API and check contents */
    for (i=0; i<w_len; i++) buffer[i] = -1;
    if (coll_io)
        err = ncmpi_get_varn_int_all(ncid, varid[1], num_reqs, starts, counts, buffer);
    else
        err = ncmpi_get_varn_int(ncid, varid[1], num_reqs, starts, counts, buffer);
    CHECK_ERR

    for (i=0; i<w_len; i++) {
        if (buffer[i] != rank+10) {
            printf("Error at line %d in %s: expecting rec_var buffer[%d]=%d but got %d\n",
                   __LINE__,__FILE__,i,rank+10,buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

    /* test flexible API, using a noncontiguous buftype */
    MPI_Datatype buftype;
    MPI_Type_vector((int)w_len, 1, 2, MPI_INT, &buftype);
    MPI_Type_commit(&buftype);
    free(buffer);
    buffer = (int*) malloc(sizeof(int) * w_len * 2);
    for (i=0; i<2*w_len; i++) buffer[i] = -1;
    if (coll_io)
        err = ncmpi_get_varn_all(ncid, varid[1], num_reqs, starts, counts, buffer, 1, buftype);
    else
        err = ncmpi_get_varn(ncid, varid[1], num_reqs, starts, counts, buffer, 1, buftype);
    CHECK_ERR
    MPI_Type_free(&buftype);

    for (i=0; i<w_len*2; i++) {
        if (i%2 && buffer[i] != -1) {
            printf("Error at line %d in %s: expecting rec_var buffer[%d]=-1 but got %d\n",
                   __LINE__,__FILE__,i,buffer[i]);
            nerrs++;
            goto err_out;
        }
        if (i%2 == 0 && buffer[i] != rank+10) {
            printf("Error at line %d in %s: expecting rec_var buffer[%d]=%d but got %d\n",
                   __LINE__,__FILE__,i,rank+10,buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

err_out:
    err = ncmpi_close(ncid);
    CHECK_ERR

    if (buffer != NULL) free(buffer);
    if (r_buffer != NULL) free(r_buffer);
    if (num_reqs > 0) {
        free(starts[0]);
        free(counts[0]);
        free(starts);
        free(counts);
    }

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int err, nerrs=0;
    MPI_Info info_dup;

    MPI_Info_dup(info, &info_dup);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

#ifdef DEBUG
    if (nprocs != 4 && rank == 0)
        printf("Warning: %s is intended to run on 4 processes\n",
                basename(__FILE__));
#endif

    nerrs = tst_io(out_path, coll_io, info_dup);
    if (nerrs > 0) goto err_out;

    /* disable PnetCDF internal buffering */
    MPI_Info_set(info, "nc_ibuf_size", "0");

    nerrs = tst_io(out_path, coll_io, info_dup);
    if (nerrs > 0) goto err_out;

err_out:
    MPI_Info_free(&info_dup);

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

    err = tst_main(argc, argv, "ncmpi_put_varn_int_all()", opt, test_io);

    MPI_Finalize();

    return err;
}
