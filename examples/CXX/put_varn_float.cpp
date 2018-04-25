/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to use a single call of NcmpiVar::putVarn_all()
 * to write a sequence of one-element requests with arbitrary array indices.
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicxx -O2 -o put_varn_float put_varn_float.cpp -lpnetcdf
 *    % mpiexec -n 4 ./put_varn_float /pvfs2/wkliao/testfile.nc
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
 *       3, 3, 3, 1, 1, 0, 0, 2, 1, 1,
 *       0, 2, 2, 2, 3, 1, 1, 2, 2, 2,
 *       1, 1, 2, 3, 3, 3, 0, 0, 1, 1,
 *       0, 0, 0, 2, 1, 1, 1, 3, 3, 3 ;
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

#define NY 4
#define NX 10
#define NDIMS 2

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
    int i, verbose=1, rank, nprocs, num_reqs;
    float *buffer;
    MPI_Offset **starts=NULL;

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

    if (nprocs != 4 && rank == 0 && verbose)
        printf("Warning: this program is intended to run on 4 processes\n");

    try {
        /* create a new file for writing -------------------------------------*/
        NcmpiFile nc(MPI_COMM_WORLD, filename, NcmpiFile::replace,
                     NcmpiFile::classic);

        /* define dimensions Y and X */
        vector<NcmpiDim> dimid(2);

        dimid[0] = nc.addDim("Y", NY);
        dimid[1] = nc.addDim("X", NX);

        /* define a 2D variable of integer type */
        NcmpiVar var = nc.addVar("var", ncmpiFloat, dimid);

        /* pick arbitrary numbers of requests for 4 processes */
        num_reqs = 0;
        if (rank == 0)      num_reqs = 8;
        else if (rank == 1) num_reqs = 13;
        else if (rank == 2) num_reqs = 9;
        else if (rank == 3) num_reqs = 10;

        if (num_reqs > 0) {
            starts    = (MPI_Offset**)malloc(num_reqs*       sizeof(MPI_Offset*));
            starts[0] = (MPI_Offset*) calloc(num_reqs*NDIMS, sizeof(MPI_Offset));
            for (i=1; i<num_reqs; i++)
                starts[i] = starts[i-1] + NDIMS;
        }

        /* assign arbitrary starts */
        const int y=0, x=1;
        if (rank == 0) {
            starts[0][y] = 0; starts[0][x] = 5;
            starts[1][y] = 1; starts[1][x] = 0;
            starts[2][y] = 2; starts[2][x] = 6;
            starts[3][y] = 3; starts[3][x] = 0;
            starts[4][y] = 0; starts[4][x] = 6;
            starts[5][y] = 2; starts[5][x] = 7;
            starts[6][y] = 3; starts[6][x] = 1;
            starts[7][y] = 3; starts[7][x] = 2;
            /* rank 0 is writing the following locations: ("-" means skip)
                      -  -  -  -  -  0  0  -  -  -
                      0  -  -  -  -  -  -  -  -  -
                      -  -  -  -  -  -  0  0  -  -
                      0  0  0  -  -  -  -  -  -  -
             */
        } else if (rank ==1) {
            starts[ 0][y] = 0; starts[ 0][x] = 3;
            starts[ 1][y] = 0; starts[ 1][x] = 8;
            starts[ 2][y] = 1; starts[ 2][x] = 5;
            starts[ 3][y] = 2; starts[ 3][x] = 0;
            starts[ 4][y] = 2; starts[ 4][x] = 8;
            starts[ 5][y] = 3; starts[ 5][x] = 4;
            starts[ 6][y] = 0; starts[ 6][x] = 4;
            starts[ 7][y] = 0; starts[ 7][x] = 9;
            starts[ 8][y] = 1; starts[ 8][x] = 6;
            starts[ 9][y] = 2; starts[ 9][x] = 1;
            starts[10][y] = 2; starts[10][x] = 9;
            starts[11][y] = 3; starts[11][x] = 5;
            starts[12][y] = 3; starts[12][x] = 6;
            /* rank 1 is writing the following locations: ("-" means skip)
                      -  -  -  1  1  -  -  -  1  1
                      -  -  -  -  -  1  1  -  -  -
                      1  1  -  -  -  -  -  -  1  1
                      -  -  -  -  1  1  1  -  -  -
             */
        } else if (rank ==2) {
            starts[0][y] = 0; starts[0][x] = 7;
            starts[1][y] = 1; starts[1][x] = 1;
            starts[2][y] = 1; starts[2][x] = 7;
            starts[3][y] = 2; starts[3][x] = 2;
            starts[4][y] = 3; starts[4][x] = 3;
            starts[5][y] = 1; starts[5][x] = 2;
            starts[6][y] = 1; starts[6][x] = 8;
            starts[7][y] = 1; starts[7][x] = 3;
            starts[8][y] = 1; starts[8][x] = 9;
            /* rank 2 is writing the following locations: ("-" means skip)
                      -  -  -  -  -  -  -  2  -  -
                      -  2  2  2  -  -  -  2  2  2
                      -  -  2  -  -  -  -  -  -  -
                      -  -  -  2  -  -  -  -  -  -
             */
        } else if (rank ==3) {
            starts[0][y] = 0; starts[0][x] = 0;
            starts[1][y] = 1; starts[1][x] = 4;
            starts[2][y] = 2; starts[2][x] = 3;
            starts[3][y] = 3; starts[3][x] = 7;
            starts[4][y] = 0; starts[4][x] = 1;
            starts[5][y] = 2; starts[5][x] = 4;
            starts[6][y] = 3; starts[6][x] = 8;
            starts[7][y] = 0; starts[7][x] = 2;
            starts[8][y] = 2; starts[8][x] = 5;
            starts[9][y] = 3; starts[9][x] = 9;
            /* rank 3 is writing the following locations: ("-" means skip)
                      3  3  3  -  -  -  -  -  -  -
                      -  -  -  -  3  -  -  -  -  -
                      -  -  -  3  3  3  -  -  -  -
                      -  -  -  -  -  -  -  3  3  3
             */
        }

        /* allocate I/O buffer and initialize its contents */
        buffer = (float*) malloc(num_reqs * sizeof(float));
        for (i=0; i<num_reqs; i++) buffer[i] =  (float)rank;

        /* set the buffer pointers to different offsets to the I/O buffer */
        var.putVarn_all(num_reqs, starts, NULL, buffer);

        free(buffer);
        if (num_reqs > 0) {
            free(starts[0]);
            free(starts);
        }

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

