/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example makes a number of nonblocking API calls, each writes a single
 * column of a 2D integer array. Each process writes NX columns and any two
 * consecutive columns are of nprocs columns distance apart from each other.
 * In this case, the fileview of each process interleaves with all other
 * processes.
 * If simply concatenating fileviews of all the nonblocking calls will result
 * in a fileview that violates the MPI-IO requirement on the fileview of which
 * flattened file offsets must be monotonically non-decreasing. PnetCDF handles
 * this case by breaking down each nonblocking call into a lsit of offset-length
 * pairs, merging the pairs across multiple nonblocking calls, and sorting
 * them into an increasing order. The sorted pairs are used to construct a
 * fileview that meets the monotonically non-decreasing offset requirement,
 * and thus the nonblocking requests can be serviced by a single MPI-IO call. 
 * 
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicxx -O2 -o column_wise column_wise.cpp -lpnetcdf
 *    % mpiexec -l -n 4 ./column_wise /pvfs2/wkliao/testfile.nc
 *    0:  0: myOff=  0 myNX=  4
 *    1:  1: myOff=  4 myNX=  4
 *    2:  2: myOff=  8 myNX=  4
 *    3:  3: myOff= 12 myNX=  4
 *    0:  0: start=  0   0 count= 10   1
 *    1:  1: start=  0   1 count= 10   1
 *    2:  2: start=  0   2 count= 10   1
 *    3:  3: start=  0   3 count= 10   1
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
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
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
using namespace std;

#include <string.h>
#include <pnetcdf>

using namespace PnetCDF;
using namespace PnetCDF::exceptions;

#define NY 10
#define NX 4

int main(int argc, char** argv) {
    char filename[128];
    int i, j, debug, rank, nprocs, err, myNX, G_NX, myOff, num_reqs;
    int ncid, cmode, varid, dimid[2], *reqs, *sts, **buf;
    vector<MPI_Offset> start(2), count(2);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    debug = 1;
    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    if (argc == 2) strcpy(filename, argv[1]);
    else           strcpy(filename, "testfile.nc");

    try {
        NcmpiFile nc(MPI_COMM_WORLD, filename, NcmpiFile::replace,
                     NcmpiFile::data64bits);

        /* the global array is NY * (NX * nprocs) */
        G_NX  = NX * nprocs;
        myOff = NX * rank;
        myNX  = NX;
        if (debug) printf("%2d: myOff=%3d myNX=%3d\n",rank,myOff,myNX);

        /* define dimensions Y and X */
        vector<NcmpiDim> dimid(2);

        dimid[0] = nc.addDim("Y", NY);
        dimid[1] = nc.addDim("X", G_NX);

        /* define a 2D variable of integer type */
        NcmpiVar var = nc.addVar("var", ncmpiInt, dimid);

        /* First, fill the entire array with zeros, using a blocking I/O.
           Every process writes a subarray of size NY * myNX */
        buf    = (int**) malloc(myNX * sizeof(int*));
        buf[0] = (int*)  calloc(NY * myNX, sizeof(int));
        start[0] = 0;   start[1] = myOff;
        count[0] = NY;  count[1] = myNX;

        var.putVar_all(start, count, &buf[0][0]);

        free(buf[0]);

        /* initialize the buffer with rank ID. Also make the case interesting,
           by allocatsing buffersd separately */
        for (i=0; i<myNX; i++) {
            buf[i] = (int*) malloc(NY * sizeof(int));
            for (j=0; j<NY; j++) buf[i][j] = rank;
        }

        reqs = (int*) malloc(myNX * sizeof(int));
        sts  = (int*) malloc(myNX * sizeof(int));

        /* each proc writes myNX single columns of the 2D array */
        start[0]  = 0;   start[1] = rank;
        count[0]  = NY;  count[1] = 1;
        if (debug)
            printf("%2d: start=%3lld %3lld count=%3lld %3lld\n",
                   rank, start[0],start[1], count[0],count[1]);

        num_reqs = 0;
        for (i=0; i<myNX; i++) {
            var.iputVar(start, count, &buf[0][0], &reqs[num_reqs++]);
            start[1] += nprocs;
        }
        nc.Wait_all(num_reqs, reqs, sts);

        /* check status of all requests */
        for (i=0; i<num_reqs; i++)
            if (sts[i] != NC_NOERR)
                printf("Error: nonblocking write fails on request %d (%s)\n",
                       i, ncmpi_strerror(sts[i]));

        free(sts);
        free(reqs);
        for (i=0; i<myNX; i++) free(buf[i]);
        free(buf);

        /* file is close implicitly */
    }
    catch(NcmpiException& e) {
       cout << e.what() << " error code=" << e.errorCode() << " Error!\n";
       return 1;
    }

    MPI_Finalize();
    return 0;
}

