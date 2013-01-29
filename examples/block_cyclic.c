/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>

#define NY 10
#define NX 4

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example, generalized from column_wise.c, makes a number of nonblocking
 * API calls, each writes a block of columns into a 2D integer array. In other
 * words, the I/O pattern is a blocked cyclic along X dimension.
 *
 * Each process writes NX columns in total. The block length is controlled by
 * block_len. In this example, block_len is set to 2. Blocks are layout in a
 * cyclic fashion in the file. This example can test if PnetCDF can coalesce
 * file offsets and lengths when constructing a merged filetype.
 * 
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file. In this example, block_len = 2.
 *
 *    % mpicc -g -o column_wise column_wise.c -lpnetcdf
 *    % mpiexec -l -n 4 column_wise testfile.nc
 *    0:  0: myOff=  0 myNX=  4
 *    1:  1: myOff=  4 myNX=  4
 *    2:  2: myOff=  8 myNX=  4
 *    3:  3: myOff= 12 myNX=  4
 *
 *    0: [i=0] start=  0   0 count= 10   1
 *    0: [i=1] start=  0   1 count= 10   1
 *    0: [i=2] start=  0   8 count= 10   1
 *    0: [i=3] start=  0   9 count= 10   1
 *    1: [i=0] start=  0   2 count= 10   1
 *    1: [i=1] start=  0   3 count= 10   1
 *    1: [i=2] start=  0  10 count= 10   1
 *    1: [i=3] start=  0  11 count= 10   1
 *    2: [i=0] start=  0   4 count= 10   1
 *    2: [i=1] start=  0   5 count= 10   1
 *    2: [i=2] start=  0  12 count= 10   1
 *    2: [i=3] start=  0  13 count= 10   1
 *    3: [i=0] start=  0   6 count= 10   1
 *    3: [i=1] start=  0   7 count= 10   1
 *    3: [i=2] start=  0  14 count= 10   1
 *    3: [i=3] start=  0  15 count= 10   1
 *
 *    % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            Y = 10 ;
 *            X = 16 ;
 *    variables:
 *            int var(Y, X) ;
 *    data:
 *
 *     var =
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define MIN(a,b) (((a)<(b))?(a):(b))

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}

int main(int argc, char** argv) {
    int i, j, debug, rank, nprocs, err, myNX, G_NX, myOff;
    int ncid, cmode, varid, dimid[2], *reqs, *sts, **buf;
    int block_start, block_len;
    MPI_Offset start[2], count[2], stride[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    debug = 1;
    if (argc != 2) {
        if (!rank) printf("Usage: %s filename\n",argv[0]);
        MPI_Finalize();
        return 0;
    }

    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, argv[1], cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* the global array is NY * (NX * nprocs) */
    G_NX  = NX * nprocs;
    myOff = NX * rank;
    myNX   = NX;
    if (debug) printf("%2d: myOff=%3d myNX=%3d\n",rank,myOff,myNX);

    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]);
    ERR
    err = ncmpi_def_dim(ncid, "X", G_NX, &dimid[1]);
    ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid);
    ERR
    err = ncmpi_enddef(ncid);
    ERR

    /* First, fill the entire array with zeros, using a blocking I/O.
       Every process writes a subarray of size NY * myNX */
    buf    = (int**) malloc(myNX * sizeof(int*));
    buf[0] = (int*)  calloc(NY * myNX, sizeof(int));
    start[0] = 0;   start[1] = myOff;
    count[0] = NY;  count[1] = myNX;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf[0]);
    free(buf[0]);

    /* initialize the buffer with rank ID. Also make the case interesting,
       by allocatsing buffersd separately */
    for (i=0; i<myNX; i++) {
        buf[i] = (int*) malloc(NY * sizeof(int));
        for (j=0; j<NY; j++) buf[i][j] = rank+10;
    }

    reqs = (int*) malloc(myNX * sizeof(int));
    sts  = (int*) malloc(myNX * sizeof(int));

    /* each proc writes myNX columns of the 2D array, block_len controls the
       the number of contiguous columns at a time */
    block_start = 0;
    block_len   = 2;  /* can be 1, 2, 3, ..., myNX */
    if (block_len > myNX) block_len = myNX;

    start[0]  = 0;   start[1] = rank * block_len;
    count[0]  = NY;  count[1] = 1;
    for (i=0; i<myNX; i++) {
        err = ncmpi_iput_vara_int(ncid, varid, start, count, buf[i], &reqs[i]);
        ERR

        if (debug)
            printf("[i=%d] start=%3lld %3lld count=%3lld %3lld\n",
                   i, start[0],start[1], count[0],count[1]);

        if (i % block_len == block_len-1)  {
            int stride = MIN(myNX-1-i, block_len);
            block_start += block_len * nprocs;
            start[1] = block_start + stride * rank;
        }
        else
            start[1]++;
    }
    err = ncmpi_wait_all(ncid, myNX, reqs, sts);
    ERR
    err = ncmpi_close(ncid);
    ERR

    free(sts);
    free(reqs);
    for (i=0; i<myNX; i++) free(buf[i]);
    free(buf);

    MPI_Finalize();
    return 0;
}

