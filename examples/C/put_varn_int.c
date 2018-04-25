/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to use a single call of ncmpi_put_varn_int_all() to
 * write a sequence of requests with arbitrary array indices and lengths.
 * Using ncmpi_put_varn_int_all() can achieve the same effect of HDF5 writing
 * a sequence of selected file locations through the following 2 APIs.
 *
 *   H5Sselect_elements(fid, H5S_SELECT_SET, NUMP, (const hssize_t **)coord);
 *   H5Dwrite(dataset, H5T_NATIVE_INT, mid, fid, H5P_DEFAULT, val);
 *
 * Note that in ncmpi_put_varn_int_all(), users can write more than one
 * element starting at each selected location.
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicc -O2 -o put_varn_int put_varn_int.c -lpnetcdf
 *    % mpiexec -n 4 ./put_varn_int /pvfs2/wkliao/testfile.nc
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
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

#define NY 4
#define NX 10
#define NDIMS 2

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [filename] output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv)
{
    extern int optind;
    char filename[256];
    int i, j, rank, nprocs, verbose=1, err, nerrs=0;
    int ncid, cmode, varid, dimid[2], num_reqs, *buffer;
    MPI_Offset w_len, **starts=NULL, **counts=NULL;

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

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* create a global array of size NY * NX */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]);
    ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]);
    ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, NDIMS, dimid, &varid);
    ERR
    err = ncmpi_enddef(ncid);
    ERR

    /* pick arbitrary numbers of requests for 4 processes */
    num_reqs = 0;
    if (rank == 0)      num_reqs = 4;
    else if (rank == 1) num_reqs = 6;
    else if (rank == 2) num_reqs = 5;
    else if (rank == 3) num_reqs = 4;

    if (num_reqs > 0) {
        starts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
        counts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
        starts[0] = (MPI_Offset*)  calloc(num_reqs*NDIMS, sizeof(MPI_Offset));
        counts[0] = (MPI_Offset*)  calloc(num_reqs*NDIMS, sizeof(MPI_Offset));
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
    buffer = (int*) malloc(w_len * sizeof(int));
    for (i=0; i<w_len; i++) buffer[i] = rank;

    /* set the buffer pointers to different offsets to the I/O buffer */
    err = ncmpi_put_varn_int_all(ncid, varid, num_reqs, starts,
                                 counts, buffer);
    ERR

    err = ncmpi_close(ncid);
    ERR

    free(buffer);
    if (num_reqs > 0) {
        free(starts[0]);
        free(counts[0]);
        free(starts);
        free(counts);
    }

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

