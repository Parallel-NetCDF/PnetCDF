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
 * this case by breaking down each nonblocking call into a list of offset-length
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

#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <pnetcdf>

using namespace PnetCDF;
using namespace PnetCDF::exceptions;

#define NY 10
#define NX 4

static void
usage(char *argv0)
{
    cerr <<
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [filename] output netCDF file name\n"
    << argv0;
}

int main(int argc, char** argv)
{
    extern int optind;
    char filename[256];
    int i, j, verbose=1, rank, nprocs, myNX, G_NX, myOff, num_reqs;
    int *reqs, *sts, **buf;
    vector<MPI_Offset> start(2), count(2);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    try {
        NcmpiFile nc(MPI_COMM_WORLD, filename, NcmpiFile::replace,
                     NcmpiFile::classic5);

        /* the global array is NY * (NX * nprocs) */
        G_NX  = NX * nprocs;
        myOff = NX * rank;
        myNX  = NX;
        if (verbose) printf("%2d: myOff=%3d myNX=%3d\n",rank,myOff,myNX);

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
           by allocating buffers separately */
        for (i=0; i<myNX; i++) {
            buf[i] = (int*) malloc(NY * sizeof(int));
            for (j=0; j<NY; j++) buf[i][j] = rank;
        }

        reqs = (int*) malloc(myNX * sizeof(int));
        sts  = (int*) malloc(myNX * sizeof(int));

        /* each proc writes myNX single columns of the 2D array */
        start[0]  = 0;   start[1] = rank;
        count[0]  = NY;  count[1] = 1;
        if (verbose)
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

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    int err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Finalize();
    return 0;
}

