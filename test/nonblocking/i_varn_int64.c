/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example tests nonblocking varn APIs, including
 * ncmpi_iput_varn_longlong(), ncmpi_iget_varn_longlong(), ncmpi_iput_varn(),
 * and ncmpi_iget_varn().
 * It first writes a sequence of requests with arbitrary array indices and
 * lengths to four variables of type NC_INT64, and reads back.
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicc -O2 -o i_varn_int64 i_varn_int64.c -lpnetcdf
 *    % mpiexec -n 4 ./i_varn_int64 /pvfs2/wkliao/testfile.nc
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *             Y = 4 ;
 *             X = 10 ;
 *    	       time = UNLIMITED ; // (4 currently)
 *    variables:
 *            int64 var0(Y, X) ;
 *            int64 var1(Y, X) ;
 *            int64 var2(Y, X) ;
 *            int64 var3(Y, X) ;
 *            int64 t_var0(time, X) ;
 *            int64 t_var1(time, X) ;
 *            int64 t_var2(time, X) ;
 *            int64 t_var3(time, X) ;
 *    data:
 *
 *     var0 =
 *      13, 13, 13, 11, 11, 10, 10, 12, 11, 11,
 *      10, 12, 12, 12, 13, 11, 11, 12, 12, 12,
 *      11, 11, 12, 13, 13, 13, 10, 10, 11, 11,
 *      10, 10, 10, 12, 11, 11, 11, 13, 13, 13 ;
 *
 *     var1 =
 *      12, 12, 12, 10, 10, 13, 13, 11, 10, 10,
 *      13, 11, 11, 11, 12, 10, 10, 11, 11, 11,
 *      10, 10, 11, 12, 12, 12, 13, 13, 10, 10,
 *      13, 13, 13, 11, 10, 10, 10, 12, 12, 12 ;
 *
 *     var2 =
 *      11, 11, 11, 13, 13, 12, 12, 10, 13, 13,
 *      12, 10, 10, 10, 11, 13, 13, 10, 10, 10,
 *      13, 13, 10, 11, 11, 11, 12, 12, 13, 13,
 *      12, 12, 12, 10, 13, 13, 13, 11, 11, 11 ;
 *
 *     var3 =
 *      10, 10, 10, 12, 12, 11, 11, 13, 12, 12,
 *      11, 13, 13, 13, 10, 12, 12, 13, 13, 13,
 *      12, 12, 13, 10, 10, 10, 11, 11, 12, 12,
 *      11, 11, 11, 13, 12, 12, 12, 10, 10, 10 ;
 *
 *     t_var0 =
 *      13, 13, 13, 11, 11, 10, 10, 12, 11, 11,
 *      10, 12, 12, 12, 13, 11, 11, 12, 12, 12,
 *      11, 11, 12, 13, 13, 13, 10, 10, 11, 11,
 *      10, 10, 10, 12, 11, 11, 11, 13, 13, 13 ;
 *
 *     t_var1 =
 *      12, 12, 12, 10, 10, 13, 13, 11, 10, 10,
 *      13, 11, 11, 11, 12, 10, 10, 11, 11, 11,
 *      10, 10, 11, 12, 12, 12, 13, 13, 10, 10,
 *      13, 13, 13, 11, 10, 10, 10, 12, 12, 12 ;
 *
 *     t_var2 =
 *      11, 11, 11, 13, 13, 12, 12, 10, 13, 13,
 *      12, 10, 10, 10, 11, 13, 13, 10, 10, 10,
 *      13, 13, 10, 11, 11, 11, 12, 12, 13, 13,
 *      12, 12, 12, 10, 13, 13, 13, 11, 11, 11 ;
 *
 *     t_var3 =
 *      10, 10, 10, 12, 12, 11, 11, 13, 12, 12,
 *      11, 13, 13, 13, 10, 12, 12, 13, 13, 13,
 *      12, 12, 13, 10, 10, 10, 11, 11, 12, 12,
 *      11, 11, 11, 13, 12, 12, 12, 10, 10, 10 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 4
#define NX 10
#define NDIMS 2

#define FATAL_ERR \
    if (err != NC_NOERR) { \
        printf("Error at line %d in %s: %s\n", __LINE__, __FILE__, ncmpi_strerrno(err)); \
        exit(1); \
    }

#define ERRS(n,a) { \
    int _i; \
    for (_i=0; _i<(n); _i++) { \
        if ((a)[_i] != NC_NOERR) { \
            printf("Error at line %d in %s: err[%d] %s\n", __LINE__, __FILE__, _i, \
                   ncmpi_strerrno((a)[_i])); \
            nerrs++; \
        } \
    } \
}

static
int clear_file_contents(int ncid, int *varid)
{
    int i, err, rank, nerrs=0;
    MPI_Offset start[2], count[2];
    long long *w_buffer = (long long*) malloc(NY*NX * sizeof(long long));

    for (i=0; i<NY*NX; i++) w_buffer[i] = -1;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    start[0] = start[1] = count[0] = count[1] = 0;
    if (rank == 0) { /* only rank 0 writes */
        count[0] = NY;
        count[1] = NX;
    }

    for (i=0; i<4; i++) {
        /* cannot use var APIs for record variables */
        err = ncmpi_iput_vara_longlong(ncid, varid[i], start, count, w_buffer, NULL);
        CHECK_ERR
    }
    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    free(w_buffer);

    /* When using burst buffering, flush the log to prevent new value being
     * skipped due to overlaping domain
     */
    err = ncmpi_flush(ncid); CHECK_ERR

    return nerrs;
}

static
int check_contents_for_fail(int ncid, int *varid)
{
    /* all processes read entire variables back and check contents */
    int i, j, nerrs=0, err, nprocs;
    long long expected[4][NY*NX] = {{13, 13, 13, 11, 11, 10, 10, 12, 11, 11,
                                     10, 12, 12, 12, 13, 11, 11, 12, 12, 12,
                                     11, 11, 12, 13, 13, 13, 10, 10, 11, 11,
                                     10, 10, 10, 12, 11, 11, 11, 13, 13, 13},
                                    {12, 12, 12, 10, 10, 13, 13, 11, 10, 10,
                                     13, 11, 11, 11, 12, 10, 10, 11, 11, 11,
                                     10, 10, 11, 12, 12, 12, 13, 13, 10, 10,
                                     13, 13, 13, 11, 10, 10, 10, 12, 12, 12},
                                    {11, 11, 11, 13, 13, 12, 12, 10, 13, 13,
                                     12, 10, 10, 10, 11, 13, 13, 10, 10, 10,
                                     13, 13, 10, 11, 11, 11, 12, 12, 13, 13,
                                     12, 12, 12, 10, 13, 13, 13, 11, 11, 11},
                                    {10, 10, 10, 12, 12, 11, 11, 13, 12, 12,
                                     11, 13, 13, 13, 10, 12, 12, 13, 13, 13,
                                     12, 12, 13, 10, 10, 10, 11, 11, 12, 12,
                                     11, 11, 11, 13, 12, 12, 12, 10, 10, 10}};

    long long *r_buffer = (long long*) malloc(NY*NX * sizeof(long long));

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs > 4) MPI_Barrier(MPI_COMM_WORLD);

    for (i=0; i<4; i++) {
        for (j=0; j<NY*NX; j++) r_buffer[j] = -1;
        err = ncmpi_get_var_longlong_all(ncid, varid[i], r_buffer);
        CHECK_ERR

        /* check if the contents of buf are expected */
        for (j=0; j<NY*NX; j++) {
            if (expected[i][j] >= nprocs) continue;
            if (r_buffer[j] != expected[i][j]) {
                printf("Error at line %d in %s: xxpect read buf[%d][%d]=%lld, but got %lld\n",
                       __LINE__,__FILE__,i,j,expected[i][j],r_buffer[j]);
                nerrs++;
            }
        }
    }
    free(r_buffer);
    return nerrs;
}

static int
check_num_pending_reqs(int ncid, int expected, int lineno)
/* check if PnetCDF can reports expected number of pending requests */
{
    int err, nerrs=0, n_pendings;
    err = ncmpi_inq_nreqs(ncid, &n_pendings);
    CHECK_ERR
    if (n_pendings != expected) {
        printf("Error at line %d in %s: expect %d pending requests but got %d\n",
               lineno, __FILE__, expected, n_pendings);
        nerrs++;
    }
    return nerrs;
}

/* swap two rows, a and b, of a 2D array */
static
void permute(MPI_Offset *a, MPI_Offset *b)
{
    int i;
    MPI_Offset tmp;
    for (i=0; i<NDIMS; i++) {
        tmp = a[i]; a[i] = b[i]; b[i] = tmp;
    }
}

static int
test_varn(int ncid, int rank, int *varid)
{
    int i, j, k, err, nerrs=0, bufsize=0;
    int nreqs, reqs[4], sts[4];
    long long *buffer[4], *cbuffer[4];
    int num_segs[4] = {4, 6, 5, 4};
    int req_lens[4], my_nsegs[4];
    MPI_Offset **starts[4], **counts[4];
    MPI_Offset n_starts[4][6][2] = {{{0,5}, {1,0}, {2,6}, {3,0}, {0,0}, {0,0}},
                                    {{0,3}, {0,8}, {1,5}, {2,0}, {2,8}, {3,4}},
                                    {{0,7}, {1,1}, {1,7}, {2,2}, {3,3}, {0,0}},
                                    {{0,0}, {1,4}, {2,3}, {3,7}, {0,0}, {0,0}}};
    MPI_Offset n_counts[4][6][2] = {{{1,2}, {1,1}, {1,2}, {1,3}, {0,0}, {0,0}},
                                    {{1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,3}},
                                    {{1,1}, {1,3}, {1,3}, {1,1}, {1,1}, {0,0}},
                                    {{1,3}, {1,1}, {1,3}, {1,3}, {0,0}, {0,0}}};

    /* n_starts[0][][] n_counts[0][][] indicate the following: ("-" means skip)
              -  -  -  -  -  X  X  -  -  -
              X  -  -  -  -  -  -  -  -  -
              -  -  -  -  -  -  X  X  -  -
              X  X  X  -  -  -  -  -  -  -
       n_starts[1][][] n_counts[1][][] indicate the following pattern.
              -  -  -  X  X  -  -  -  X  X
              -  -  -  -  -  X  X  -  -  -
              X  X  -  -  -  -  -  -  X  X
              -  -  -  -  X  X  X  -  -  -
       n_starts[2][][] n_counts[2][][] indicate the following pattern.
              -  -  -  -  -  -  -  X  -  -
              -  X  X  X  -  -  -  X  X  X
              -  -  X  -  -  -  -  -  -  -
              -  -  -  X  -  -  -  -  -  -
       n_starts[3][][] n_counts[3][][] indicate the following pattern.
              X  X  X  -  -  -  -  -  -  -
              -  -  -  -  X  -  -  -  -  -
              -  -  -  X  X  X  -  -  -  -
              -  -  -  -  -  -  -  X  X  X
     */

    /* allocate space for starts and counts */
    starts[0] = (MPI_Offset**) malloc(4 * 6 * sizeof(MPI_Offset*));
    counts[0] = (MPI_Offset**) malloc(4 * 6 * sizeof(MPI_Offset*));
    starts[0][0] = (MPI_Offset*) calloc(4 * 6 * NDIMS, sizeof(MPI_Offset));
    counts[0][0] = (MPI_Offset*) calloc(4 * 6 * NDIMS, sizeof(MPI_Offset));
    for (i=1; i<4; i++) {
        starts[i] = starts[i-1] + 6;
        counts[i] = counts[i-1] + 6;
        starts[i][0] = starts[i-1][0] + 6 * NDIMS;
        counts[i][0] = counts[i-1][0] + 6 * NDIMS;
    }
    for (i=0; i<4; i++) {
        for (j=1; j<6; j++) {
            starts[i][j] = starts[i][j-1] + NDIMS;
            counts[i][j] = counts[i][j-1] + NDIMS;
        }
    }

    /* set values for starts and counts */
    for (i=0; i<4; i++) {
        int n = (i + rank) % 4;
        my_nsegs[i] = num_segs[n]; /* number of segments for this request */
        for (j=0; j<6; j++) {
            for (k=0; k<NDIMS; k++) {
                starts[i][j][k] = n_starts[n][j][k];
                counts[i][j][k] = n_counts[n][j][k];
            }
        }
    }

    /* test error code: NC_ENULLSTART */
    err = ncmpi_iput_varn_longlong(ncid, varid[0], 1, NULL, NULL,
                                   NULL, &reqs[0]);
    if (err != NC_ENULLSTART) {
        printf("expecting error code NC_ENULLSTART but got %s\n",
               nc_err_code_name(err));
        nerrs++;
    }

    /* only rank 0, 1, 2, and 3 do I/O:
     * each of ranks 0 to 3 write 4 nonblocking requests */
    nreqs = 4;
    if (rank >= 4) nreqs = 0;

    /* calculate length of each varn request and allocate write buffer */
    for (i=0; i<nreqs; i++) {
        req_lens[i] = 0; /* total length this request */
        for (j=0; j<my_nsegs[i]; j++) {
            MPI_Offset req_len=1;
            for (k=0; k<NDIMS; k++)
                req_len *= counts[i][j][k];
            req_lens[i] += req_len;
        }

        /* allocate I/O buffer and initialize its contents */
        buffer[i] = (long long*) malloc(req_lens[i] * sizeof(long long));
        for (j=0; j<req_lens[i]; j++) buffer[i][j] = rank+10;
    }

    /* write using varn API */
    nerrs += clear_file_contents(ncid, varid);
    for (i=0; i<nreqs; i++) {
        err = ncmpi_iput_varn_longlong(ncid, varid[i], my_nsegs[i], starts[i],
                                       counts[i], buffer[i], &reqs[i]);
        CHECK_ERR
    }
    nerrs += check_num_pending_reqs(ncid, nreqs, __LINE__);
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    CHECK_ERR
    ERRS(nreqs, sts)

    /* check if write buffer contents have been altered */
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]; j++) {
            if (buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%lld\n",
                __LINE__,__FILE__,i,j,buffer[i][j]);
                nerrs++;
            }
        }
    }

    /* all processes read entire variables back and check contents */
    nerrs += check_contents_for_fail(ncid, varid);

    /* try with buffer being a single contiguous space */
    for (i=0; i<nreqs; i++) bufsize += req_lens[i];
    if (bufsize>0) {
        cbuffer[0] = (long long*) malloc(bufsize * sizeof(long long));
        for (i=1; i<nreqs; i++) cbuffer[i] = cbuffer[i-1] + req_lens[i-1];
        for (i=0; i<bufsize; i++) cbuffer[0][i] = rank+10;
    }
    else
        for (i=0; i<4; i++) cbuffer[i] = NULL;

    /* write usning varn API */
    nerrs += clear_file_contents(ncid, varid);
    for (i=0; i<nreqs; i++) {
        err = ncmpi_iput_varn_longlong(ncid, varid[i], my_nsegs[i], starts[i],
                                       counts[i], cbuffer[i], &reqs[i]);
        CHECK_ERR
    }
    nerrs += check_num_pending_reqs(ncid, nreqs, __LINE__);
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    CHECK_ERR
    ERRS(nreqs, sts)

    /* check if write buffer contents have been altered */
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]; j++) {
            if (cbuffer[i][j] != rank+10) {
                printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%lld\n",
                __LINE__,__FILE__,i,j,cbuffer[i][j]);
                nerrs++;
            }
        }
    }

    /* all processes read entire variables back and check contents */
    nerrs += check_contents_for_fail(ncid, varid);

    /* permute write order: so starts[*] are not in an increasing order:
     * swap segment 0 with segment 2 and swap segment 1 with segment 3
     */
    for (i=0; i<nreqs; i++) {
        permute(starts[i][0], starts[i][2]); permute(counts[i][0], counts[i][2]);
        permute(starts[i][1], starts[i][3]); permute(counts[i][1], counts[i][3]);
    }

    /* write usning varn API */
    nerrs += clear_file_contents(ncid, varid);
    for (i=0; i<nreqs; i++) {
        err = ncmpi_iput_varn_longlong(ncid, varid[i], my_nsegs[i], starts[i],
                                       counts[i], buffer[i], &reqs[i]);
        CHECK_ERR
    }
    nerrs += check_num_pending_reqs(ncid, nreqs, __LINE__);
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    CHECK_ERR
    ERRS(nreqs, sts)

    /* check if write buffer contents have been altered */
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]; j++) {
            if (buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%lld\n",
                __LINE__,__FILE__,i,j,buffer[i][j]);
                nerrs++;
            }
        }
    }

    /* all processes read entire variables back and check contents */
    nerrs += check_contents_for_fail(ncid, varid);

    /* read using get_varn API and check contents */
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]; j++) buffer[i][j] = -1;
        err = ncmpi_iget_varn_longlong(ncid, varid[i], my_nsegs[i], starts[i],
                                       counts[i], buffer[i], &reqs[i]);
        CHECK_ERR
    }

    nerrs += check_num_pending_reqs(ncid, nreqs, __LINE__);
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    CHECK_ERR
    ERRS(nreqs, sts)

    /* check if read buffer contents are expected */
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]; j++) {
            if (buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: expecting var %d buffer[%d][%d]=%d but got %lld\n",
                       __LINE__,__FILE__,varid[i],i,j,rank+10,buffer[i][j]);
                nerrs++;
            }
        }
    }

    for (i=0; i<nreqs; i++) free(buffer[i]);

    /* test flexible put API, using a noncontiguous buftype */
    nerrs += clear_file_contents(ncid, varid);
    for (i=0; i<nreqs; i++) {
        MPI_Datatype buftype;
        MPI_Type_vector(req_lens[i], 1, 2, MPI_LONG_LONG, &buftype);
        MPI_Type_commit(&buftype);
        buffer[i] = (long long*) malloc(req_lens[i] * 2 * sizeof(long long));
        for (j=0; j<req_lens[i]*2; j++) buffer[i][j] = rank+10;

        err = ncmpi_iput_varn(ncid, varid[i], my_nsegs[i], starts[i],
                              counts[i], buffer[i], 1, buftype, &reqs[i]);
        CHECK_ERR
        MPI_Type_free(&buftype);
    }
    nerrs += check_num_pending_reqs(ncid, nreqs, __LINE__);
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    CHECK_ERR
    ERRS(nreqs, sts)

    /* check if write buffer contents have been altered */
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]*2; j++) {
            if (buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%lld\n",
                __LINE__,__FILE__,i,j,buffer[i][j]);
                nerrs++;
            }
        }
    }

    /* all processes read entire variables back and check contents */
    nerrs += check_contents_for_fail(ncid, varid);

    /* test flexible get API, using a noncontiguous buftype */
    for (i=0; i<nreqs; i++) {
        MPI_Datatype buftype;
        MPI_Type_vector(req_lens[i], 1, 2, MPI_LONG_LONG, &buftype);
        MPI_Type_commit(&buftype);
        for (j=0; j<req_lens[i]*2; j++) buffer[i][j] = -1;
        err = ncmpi_iget_varn(ncid, varid[i], my_nsegs[i], starts[i],
                              counts[i], buffer[i], 1, buftype, &reqs[i]);
        CHECK_ERR
        MPI_Type_free(&buftype);
    }
    nerrs += check_num_pending_reqs(ncid, nreqs, __LINE__);
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    CHECK_ERR
    ERRS(nreqs, sts)

    /* check if read buffer contents are expected */
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]*2; j++) {
            if (j%2 && buffer[i][j] != -1) {
                printf("Error at line %d in %s: expecting buffer[%d][%d]=-1 but got %lld\n",
                       __LINE__,__FILE__,i,j,buffer[i][j]);
                nerrs++;
            }
            if (j%2 == 0 && buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: expecting buffer[%d][%d]=%d but got %lld\n",
                       __LINE__,__FILE__,i,j,rank+10,buffer[i][j]);
                nerrs++;
            }
        }
    }

    /* read back using a contiguous buffer. First swap back the starts[] and
     * counts[]. swap segment 0 with segment 2 and swap segment 1 with segment
     * 3
     */
    for (i=0; i<nreqs; i++) {
        permute(starts[i][0], starts[i][2]);
        permute(counts[i][0], counts[i][2]);
        permute(starts[i][1], starts[i][3]);
        permute(counts[i][1], counts[i][3]);
    }

    for (i=0; i<bufsize; i++) cbuffer[0][i] = -1;
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]; j++) buffer[i][j] = -1;
        err = ncmpi_iget_varn_longlong(ncid, varid[i], my_nsegs[i], starts[i],
                                       counts[i], cbuffer[i], &reqs[i]);
        CHECK_ERR
    }
    nerrs += check_num_pending_reqs(ncid, nreqs, __LINE__);
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    CHECK_ERR
    ERRS(nreqs, sts)

    /* check if read buffer contents are expected */
    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]; j++) {
            if (cbuffer[i][j] != rank+10) {
                printf("Error at line %d in %s: expecting buffer[%d][%d]=%d but got %lld\n",
                       __LINE__,__FILE__,i,j,rank+10,cbuffer[i][j]);
                nerrs++;
            }
        }
    }

    if (bufsize>0) free(cbuffer[0]);
    for (i=0; i<nreqs; i++) free(buffer[i]);
    free(starts[0][0]);
    free(counts[0][0]);
    free(starts[0]);
    free(counts[0]);

    return nerrs;
}

int main(int argc, char** argv)
{
    char filename[256];
    int rank, nprocs, err, nerrs=0;
    int ncid, cmode, varid[4], dimid[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for iput/iget varn ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

#ifdef DEBUG
    if (nprocs != 4 && rank == 0)
        printf("Warning: %s is intended to run on 4 processes\n",argv[0]);
#endif

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    FATAL_ERR

    /* create fixed-size variables of size NY * NX */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT64, NDIMS, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT64, NDIMS, dimid, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var2", NC_INT64, NDIMS, dimid, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var3", NC_INT64, NDIMS, dimid, &varid[3]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* test fixed-size variables */
    nerrs += test_varn(ncid, rank, varid);

    err = ncmpi_redef(ncid); CHECK_ERR

    /* add record variables */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "t_var0", NC_INT64, NDIMS, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "t_var1", NC_INT64, NDIMS, dimid, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "t_var2", NC_INT64, NDIMS, dimid, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "t_var3", NC_INT64, NDIMS, dimid, &varid[3]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* test record variables */
    nerrs += test_varn(ncid, rank, varid);

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0) {
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
            ncmpi_inq_malloc_list();
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

