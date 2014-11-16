/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example tests ncmpi_iput_varn_longlong(), ncmpi_iget_varn_longlong(),
 * and ncmpi_iget_varn().
 * It writes a sequence of requests with arbitrary array indices and lengths,
 * and reads back.
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
 *    variables:
 *            int64 var0(Y, X) ;
 *            int64 var1(Y, X) ;
 *            int64 var2(Y, X) ;
 *            int64 var3(Y, X) ;
 *    data:
 *
 *     var0 =
 *      3, 3, 3, 1, 1, 0, 0, 2, 1, 1,
 *      0, 2, 2, 2, 3, 1, 1, 2, 2, 2,
 *      1, 1, 2, 3, 3, 3, 0, 0, 1, 1,
 *      0, 0, 0, 2, 1, 1, 1, 3, 3, 3 ;
 *
 *     var1 =
 *      2, 2, 2, 0, 0, 3, 3, 1, 0, 0,
 *      3, 1, 1, 1, 2, 0, 0, 1, 1, 1,
 *      0, 0, 1, 2, 2, 2, 3, 3, 0, 0,
 *      3, 3, 3, 1, 0, 0, 0, 2, 2, 2 ;
 *
 *     var2 =
 *      1, 1, 1, 3, 3, 2, 2, 0, 3, 3,
 *      2, 0, 0, 0, 1, 3, 3, 0, 0, 0,
 *      3, 3, 0, 1, 1, 1, 2, 2, 3, 3,
 *      2, 2, 2, 0, 3, 3, 3, 1, 1, 1 ;
 *
 *     var3 =
 *      0, 0, 0, 2, 2, 1, 1, 3, 2, 2,
 *      1, 3, 3, 3, 0, 2, 2, 3, 3, 3,
 *      2, 2, 3, 0, 0, 0, 1, 1, 2, 2,
 *      1, 1, 1, 3, 2, 2, 2, 0, 0, 0 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <mpi.h>
#include <pnetcdf.h>

#define NY 4
#define NX 10
#define NDIMS 2

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}
#define ERRS(n,a) {int _i; for(_i=0;_i<(n);_i++)if(a[_i]!=NC_NOERR)printf("Error at line=%d: err[%d] %s\n", __LINE__, _i, ncmpi_strerror(a[_i]));}

int check_contents_for_fail(int n, long long *buffer)
{
    int i;
    long long expected[4][NY*NX] = {{3, 3, 3, 1, 1, 0, 0, 2, 1, 1,
                                     0, 2, 2, 2, 3, 1, 1, 2, 2, 2,
                                     1, 1, 2, 3, 3, 3, 0, 0, 1, 1,
                                     0, 0, 0, 2, 1, 1, 1, 3, 3, 3},
                                    {2, 2, 2, 0, 0, 3, 3, 1, 0, 0,
                                     3, 1, 1, 1, 2, 0, 0, 1, 1, 1,
                                     0, 0, 1, 2, 2, 2, 3, 3, 0, 0,
                                     3, 3, 3, 1, 0, 0, 0, 2, 2, 2},
                                    {1, 1, 1, 3, 3, 2, 2, 0, 3, 3,
                                     2, 0, 0, 0, 1, 3, 3, 0, 0, 0,
                                     3, 3, 0, 1, 1, 1, 2, 2, 3, 3,
                                     2, 2, 2, 0, 3, 3, 3, 1, 1, 1},
                                    {0, 0, 0, 2, 2, 1, 1, 3, 2, 2,
                                     1, 3, 3, 3, 0, 2, 2, 3, 3, 3,
                                     2, 2, 3, 0, 0, 0, 1, 1, 2, 2,
                                     1, 1, 1, 3, 2, 2, 2, 0, 0, 0}};

    /* check if the contents of buf are expected */
    for (i=0; i<NY*NX; i++) {
        if (buffer[i] != expected[n][i]) {
            printf("n=%d Expected read buf[%d]=%lld, but got %lld\n",
                   n, i,expected[n][i],buffer[i]);
            return 1;
        }
    }
    return 0;
}

void permute(MPI_Offset a[NDIMS], MPI_Offset b[NDIMS])
{
    int i;
    MPI_Offset tmp;
    for (i=0; i<NDIMS; i++) {
        tmp = a[i]; a[i] = b[i]; b[i] = tmp;
    }
}

int main(int argc, char** argv)
{
    char filename[128];
    int i, j, k, rank, nprocs, verbose=0, err, nfails=0;
    int ncid, cmode, varid[4], dimid[2], nreqs, reqs[4], sts[4];
    long long *buffer[4], *r_buffer;
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
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    strcpy(filename, "testfile.nc");
    if (argc == 2) strcpy(filename, argv[1]);
    MPI_Bcast(filename, 128, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (nprocs != 4 && rank == 0)
        printf("Warning: %s is intended to run on 4 processes\n",argv[0]);

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* create a global array of size NY * NX */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]); ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT64, NDIMS, dimid, &varid[0]); ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT64, NDIMS, dimid, &varid[1]); ERR
    err = ncmpi_def_var(ncid, "var2", NC_INT64, NDIMS, dimid, &varid[2]); ERR
    err = ncmpi_def_var(ncid, "var3", NC_INT64, NDIMS, dimid, &varid[3]); ERR
    err = ncmpi_enddef(ncid); ERR

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
        for (j=0; j<6; j++) {
            for (k=0; k<NDIMS; k++) {
                starts[i][j][k] = n_starts[i][j][k];
                counts[i][j][k] = n_counts[i][j][k];
            }
        }
    }

    /* check error code: NC_ENULLSTART */
    err = ncmpi_iput_varn_longlong(ncid, varid[0], 1, NULL, NULL,
                                   NULL, &reqs[0]);
    if (err != NC_ENULLSTART) {
        printf("expecting error code NC_ENULLSTART=%d but got %d\n",
               NC_ENULLSTART,err);
        nfails++;
    }

    nreqs = 4;
    if (rank >= 4) nreqs = 0;

    /* each of ranks 0 to 3 write 4 nonblocking requests */
    for (i=0; i<nreqs; i++) {
        int n = (i + rank) % 4;
        if (rank >= 4) { my_nsegs[i] = req_lens[i] = 0; continue; }
        req_lens[i] = 0; /* total length this request */
        my_nsegs[i] = num_segs[n]; /* number of segments for this request */
        for (j=0; j<my_nsegs[i]; j++) {
            MPI_Offset req_len=1;
            for (k=0; k<NDIMS; k++)
                req_len *= counts[n][j][k];
            req_lens[i] += req_len;
        }
        if (verbose) printf("req_lens[%d]=%d\n",i,req_lens[i]);

        /* allocate I/O buffer and initialize its contents */
        buffer[i] = (long long*) malloc(req_lens[i] * sizeof(long long));
        for (j=0; j<req_lens[i]; j++) buffer[i][j] = rank;

        /* write usning varn API */
        err = ncmpi_iput_varn_longlong(ncid, varid[i], my_nsegs[i], starts[n],
                                       counts[n], buffer[i], &reqs[i]);
        ERR
    }
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    ERRS(nreqs, sts)

    /* read entire variables back and check contents */
    r_buffer = (long long*) malloc(NY*NX * sizeof(long long));
    for (i=0; i<4; i++) {
        for (j=0; j<NY*NX; j++) r_buffer[j] = -1;
        err = ncmpi_get_var_longlong_all(ncid, varid[i], r_buffer);
        ERR
        if (nprocs >= 4) nfails += check_contents_for_fail(i, r_buffer);
    }

    /* permute write order: so starts[*] are not in an increasing order */
    for (i=0; i<4; i++) {
        permute(starts[i][1], starts[i][2]);
        permute(counts[i][1], counts[i][2]);
        permute(starts[i][2], starts[i][3]);
        permute(counts[i][2], counts[i][3]);
    }

    /* write usning varn API */
    for (i=0; i<nreqs; i++) {
        int n = (i + rank) % 4;
        err = ncmpi_iput_varn_longlong(ncid, varid[i], my_nsegs[i], starts[n],
                                       counts[n], buffer[i], &reqs[i]);
        ERR
    }
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    ERRS(nreqs, sts)

    /* read entire variables back and check contents */
    for (i=0; i<4; i++) {
        for (j=0; j<NY*NX; j++) r_buffer[j] = -1;
        err = ncmpi_get_var_longlong_all(ncid, varid[i], r_buffer);
        ERR
        if (nprocs >= 4) nfails += check_contents_for_fail(i, r_buffer);
    }

    /* read using get_varn API and check contents */
    for (i=0; i<nreqs; i++) {
        int n = (i + rank) % 4;
        for (j=0; j<req_lens[i]; j++) buffer[i][j] = -1;
        err = ncmpi_iget_varn_longlong(ncid, varid[i], my_nsegs[i], starts[n],
                                       counts[n], buffer[i], &reqs[i]);
        ERR
    }
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    ERRS(nreqs, sts)

    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]; j++) {
            if (buffer[i][j] != rank) {
                printf("Error at line %d: expecting buffer[%d][%d]=%d but got %lld\n",
                       __LINE__,i,j,rank,buffer[i][j]);
                nfails++;
            }
        }
    }

    /* test flexible API, using a noncontiguous buftype */
    for (i=0; i<nreqs; i++) {
        int n = (i + rank) % 4;
        MPI_Datatype buftype;
        MPI_Type_vector(req_lens[i], 1, 2, MPI_LONG_LONG, &buftype);
        MPI_Type_commit(&buftype);
        free(buffer[i]);
        buffer[i] = (long long*) malloc(req_lens[i] * 2 * sizeof(long long));
        for (j=0; j<req_lens[i]*2; j++) buffer[i][j] = -1;
        err = ncmpi_iget_varn(ncid, varid[i], my_nsegs[i], starts[n],
                              counts[n], buffer[i], 1, buftype, &reqs[i]);
        ERR
        MPI_Type_free(&buftype);
    }
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    ERRS(nreqs, sts)

    for (i=0; i<nreqs; i++) {
        for (j=0; j<req_lens[i]*2; j++) {
            if (j%2 && buffer[i][j] != -1) {
                printf("Error at line %d: expecting buffer[%d][%d]=-1 but got %lld\n",
                       __LINE__,i,j,buffer[i][j]);
                nfails++;
            }
            if (j%2 == 0 && buffer[i][j] != rank) {
                printf("Error at line %d: expecting buffer[%d][%d]=%d but got %lld\n",
                       __LINE__,i,j,rank,buffer[i][j]);
                nfails++;
            }
        }
    }

    err = ncmpi_close(ncid);
    ERR

    for (i=0; i<nreqs; i++) free(buffer[i]);
    free(r_buffer);
    free(starts[0][0]);
    free(counts[0][0]);
    free(starts[0]);
    free(counts[0]);

    MPI_Allreduce(MPI_IN_PLACE, &nfails, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    char cmd_str[80];
    sprintf(cmd_str, "*** TESTING C   %s for iput/iget varn ", argv[0]);
    if (rank == 0) {
        if (nfails) printf("%-66s ------ failed\n", cmd_str);
        else        printf("%-66s ------ pass\n", cmd_str);
    }


    MPI_Finalize();
    return 0;
}

