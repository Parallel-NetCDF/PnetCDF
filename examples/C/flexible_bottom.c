/*********************************************************************
 *
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This example shows how to use PnetCDF blocking and nonblocking flexible APIs
 * by supplying an MPI derived datatype and MPI_BOTTOM in the buffer argument.
 * The MPI derived datatype is constructed through MPI_Type_create_hindexed(),
 * which defines a noncontiguous memory space consisting of two separate blocks.
 *
 * A 2D global variable of size NY * NX * nprocs of type float is defined in
 * the file. Each process allocates two write buffers and constructs an MPI
 * derived datatype for the two buffers.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -O2 -o flexible_bottom flexible_bottom.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./flexible_bottom ./testfile.nc
 *
 *    % ncmpidump ./testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            Y = 4 ;
 *            X = 16 ;
 *    variables:
 *            float var(Y, X) ;
 *    data:
 *
 *     var =
 *      10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13,
 *      10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13,
 *      50, 50, 50, 50, 51, 51, 51, 51, 52, 52, 52, 52, 53, 53, 53, 53,
 *      50, 50, 50, 50, 51, 51, 51, 51, 52, 52, 52, 52, 53, 53, 53, 53 ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

#define NY 4
#define NX 4

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

#define CHECK_MPI_ERR \
    if (err != MPI_SUCCESS) { \
        int errorStringLen; \
        char errorString[MPI_MAX_ERROR_STRING]; \
        MPI_Error_string(err, errorString, &errorStringLen); \
        printf("Error at line %d: %s\n",__LINE__, errorString); \
    }

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

/*----< pnetcdf_check_mem_usage() >------------------------------------------*/
/* check PnetCDF library internal memory usage */
static int
pnetcdf_check_mem_usage(MPI_Comm comm)
{
    int err, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && verbose)
            printf("maximum heap memory allocated by PnetCDF internally is %lld bytes\n",
                   sum_size);

        /* check if there is any PnetCDF internal malloc residue */
        err = ncmpi_inq_malloc_size(&malloc_size);
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    char filename[256];
    int i, rank, nprocs, err, nerrs=0, req, status;
    int ncid, cmode, varid, dimid[2];
    float *buf[2];
    MPI_Offset start[2], count[2];
    MPI_Datatype btype;
    int array_of_blocklengths[2];
    MPI_Aint array_of_displacements[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    verbose = 1;

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

    if (verbose && rank == 0) printf("%s: example of using MPI_BOTTOM in flexible APIs\n",__FILE__);

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "Y", NY,        &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs, &dimid[1]); ERR

    /* define a variable of size NY * (NX * nprocs) */
    err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid, &varid); ERR

    /* exit define mode */
    err = ncmpi_enddef(ncid); ERR

    /* allocate two buffers */
    buf[0] = (float*) malloc(sizeof(float) * NY * NX / 2);
    buf[1] = (float*) malloc(sizeof(float) * NY * NX / 2);
    for (i=0; i<NY*NX/2; i++) {
        buf[0][i] = rank + 10;
        buf[1][i] = rank + 50;
    }

    /* create a noncontiguous buffer datatype consisting 2 blocks */
    MPI_Get_address(buf[0], &array_of_displacements[0]);
    MPI_Get_address(buf[1], &array_of_displacements[1]);
    array_of_blocklengths[0] = NY * NX / 2;
    array_of_blocklengths[1] = NY * NX / 2;

    err = MPI_Type_create_hindexed(2, array_of_blocklengths,
                                   array_of_displacements,
                                   MPI_FLOAT, &btype);
    CHECK_MPI_ERR

    err = MPI_Type_commit(&btype);
    CHECK_MPI_ERR

    start[0] = 0;
    start[1] = NX * rank;
    count[0] = NY;
    count[1] = NX;

    /* calling a blocking flexible API */
    err = ncmpi_put_vara_all(ncid, varid, start, count, MPI_BOTTOM, 1, btype);
    ERR

    /* calling a blocking flexible API */
    err = ncmpi_get_vara_all(ncid, varid, start, count, MPI_BOTTOM, 1, btype);
    ERR

    /* calling a nonblocking flexible API */
    err = ncmpi_iput_vara(ncid, varid, start, count, MPI_BOTTOM, 1, btype, &req);
    ERR

    err = ncmpi_wait_all(ncid, 1, &req, &status); ERR

    /* calling a nonblocking flexible API */
    err = ncmpi_iget_vara(ncid, varid, start, count, MPI_BOTTOM, 1, btype, &req);
    ERR

    err = ncmpi_wait_all(ncid, 1, &req, &status); ERR

    free(buf[0]);
    free(buf[1]);
    MPI_Type_free(&btype);

    err = ncmpi_close(ncid); ERR

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

