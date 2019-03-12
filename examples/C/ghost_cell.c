/*********************************************************************
 *
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/*
 * This example shows how to use varm API to write a 2D array buffer with ghost
 * cells to a variables. The size of ghost cells is nghosts and the ghost cells
 * cells appear on both ends of each dimension. The contents of ghost cells are
 * -8s and non-ghost cells are the process rank IDs.
 *
 * To compile:
 *        mpicc -O2 ghost_cell.c -o ghost_cell -lpnetcdf
 * To run:
 *        mpiexec -n num_processes ./ghost_cell [filename] [len]
 * where len decides the size of local array, which is len x len+1.
 * So, the variable is of size len*(len+1) * nprocs * sizeof(int)
 *
 * The screen outputs from running ncmpidump on the output file produced by
 * command:
 * % mpiexec -n 4 ./ghost_cell -q /pvfs2/wkliao/testfile.nc 4
 *
 * % ncmpidump /pvfs2/wkliao/testfile.nc
 *     netcdf testfile {
 *     // file format: CDF-5 (big variables)
 *     dimensions:
 *         Y = 8 ;
 *         X = 10 ;
 *     variables:
 *         int var(Y, X) ;
 *         data:
 *
 *     var =
 *         0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
 *         0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
 *         0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
 *         0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
 *         2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
 *         2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
 *         2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
 *         2, 2, 2, 2, 2, 3, 3, 3, 3, 3 ;
 *     }
 *
 * In this case, the contents of local buffers are shown below.
 *
 * rank 0:                                rank 1:
 *    -8, -8, -8, -8, -8, -8, -8, -8, -8     -8, -8, -8, -8, -8, -8, -8, -8, -8
 *    -8, -8, -8, -8, -8, -8, -8, -8, -8     -8, -8, -8, -8, -8, -8, -8, -8, -8
 *    -8, -8,  0,  0,  0,  0,  0, -8, -8     -8, -8,  1,  1,  1,  1,  1, -8, -8
 *    -8, -8,  0,  0,  0,  0,  0, -8, -8     -8, -8,  1,  1,  1,  1,  1, -8, -8
 *    -8, -8,  0,  0,  0,  0,  0, -8, -8     -8, -8,  1,  1,  1,  1,  1, -8, -8
 *    -8, -8,  0,  0,  0,  0,  0, -8, -8     -8, -8,  1,  1,  1,  1,  1, -8, -8
 *    -8, -8, -8, -8, -8, -8, -8, -8, -8     -8, -8, -8, -8, -8, -8, -8, -8, -8
 *    -8, -8, -8, -8, -8, -8, -8, -8, -8     -8, -8, -8, -8, -8, -8, -8, -8, -8
 *
 * rank 2:                                rank 3:
 *    -8, -8, -8, -8, -8, -8, -8, -8, -8     -8, -8, -8, -8, -8, -8, -8, -8, -8
 *    -8, -8, -8, -8, -8, -8, -8, -8, -8     -8, -8, -8, -8, -8, -8, -8, -8, -8
 *    -8, -8,  2,  2,  2,  2,  2, -8, -8     -8, -8,  3,  3,  3,  3,  3, -8, -8
 *    -8, -8,  2,  2,  2,  2,  2, -8, -8     -8, -8,  3,  3,  3,  3,  3, -8, -8
 *    -8, -8,  2,  2,  2,  2,  2, -8, -8     -8, -8,  3,  3,  3,  3,  3, -8, -8
 *    -8, -8,  2,  2,  2,  2,  2, -8, -8     -8, -8,  3,  3,  3,  3,  3, -8, -8
 *    -8, -8, -8, -8, -8, -8, -8, -8, -8     -8, -8, -8, -8, -8, -8, -8, -8, -8
 *    -8, -8, -8, -8, -8, -8, -8, -8, -8     -8, -8, -8, -8, -8, -8, -8, -8, -8
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [-k format] [-l len] [file_name]\n"
    "       [-h] Print this help\n"
    "       [-q] quiet mode\n"
    "       [-k format] file format: 1 for CDF-1, 2 for CDF-2, 3 for NetCDF4,\n"
    "                                4 for NetCDF4 classic model, 5 for CDF-5\n"
    "       [-l len] size of each dimension of the local array\n"
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

/*----< pnetcdf_io() >-------------------------------------------------------*/
static int
pnetcdf_io(MPI_Comm comm, char *filename, int cmode, int len)
{
    int i, j, rank, nprocs, ncid, bufsize, err, nerrs=0;
    int psizes[2], local_rank[2], dimids[2], varid, nghosts;
    int *buf, *buf_ptr;
    MPI_Offset gsizes[2], starts[2], counts[2], imap[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* calculate number of processes along each dimension */
    psizes[0] = psizes[1] = 0;
    MPI_Dims_create(nprocs, 2, psizes);
    if (verbose && rank == 0)
        printf("psizes=%d %d\n", psizes[0], psizes[1]);

    gsizes[0] = (MPI_Offset)len     * psizes[0]; /* global array size */
    gsizes[1] = (MPI_Offset)(len+1) * psizes[1];
    if (verbose && rank == 0)
        printf("global variable shape: %lld %lld\n", gsizes[0],gsizes[1]);

    /* find its local rank IDs along each dimension */
    local_rank[0] = rank / psizes[1];
    local_rank[1] = rank % psizes[1];
    if (verbose)
        printf("rank %d: dim rank=%d %d\n", rank,local_rank[0],local_rank[1]);

    counts[0] = len;
    counts[1] = len+1;
    starts[0] = local_rank[0] * counts[0];
    starts[1] = local_rank[1] * counts[1];
    if (verbose)
        printf("starts=%lld %lld counts=%lld %lld\n",starts[0],starts[1],counts[0],counts[1]);

    /* allocate and initialize buffer with ghost cells on both ends of each dim */
    nghosts = 2;
    bufsize = (counts[0] + 2 * nghosts) * (counts[1] + 2 * nghosts);
    buf = (int *) malloc(bufsize * sizeof(int));
    for (i=0; i<counts[0]+2*nghosts; i++)
    for (j=0; j<counts[1]+2*nghosts; j++) {
        if (nghosts <= i && i < counts[0]+nghosts &&
            nghosts <= j && j < counts[1]+nghosts)
            buf[i*(counts[1]+2*nghosts) + j] = rank;
        else
            buf[i*(counts[1]+2*nghosts) + j] = -8; /* all ghost cells have value -8 */
    }

    /* create the file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_create() file %s (%s)\n",
        __LINE__,__FILE__,filename,ncmpi_strerror(err));
        MPI_Abort(comm, -1);
        exit(1);
    }

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "Y", gsizes[0], &dimids[0]); ERR
    err = ncmpi_def_dim(ncid, "X", gsizes[1], &dimids[1]); ERR

    /* define variable */
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimids, &varid); ERR

    /* exit the define mode */
    err = ncmpi_enddef(ncid); ERR

    /* set up imap[] for excluding ghost cells */
    imap[1] = 1;
    imap[0] = counts[1] + 2 * nghosts;

    /* find the first non-ghost cell of the user buf */
    buf_ptr = buf + nghosts * (counts[1]+2*nghosts + 1);

    /* write the whole variable in file */
    err = ncmpi_put_varm_int_all(ncid, varid, starts, counts, NULL, imap, buf_ptr); ERR

    /* close the file */
    err = ncmpi_close(ncid); ERR

    free(buf);

    return nerrs;
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, rank, cmode, kind=0, len=0, nerrs=0;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqk:l:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'k': kind = atoi(optarg);
                      break;
            case 'l': len = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    len = (len <= 0) ? 4 : len;

    switch (kind) {
        case(2): cmode = NC_64BIT_OFFSET;             break;
        case(3): cmode = NC_NETCDF4;                  break;
        case(4): cmode = NC_NETCDF4|NC_CLASSIC_MODEL; break;
        case(5): cmode = NC_64BIT_DATA;               break;
        default: cmode = 0;
    }

#ifndef PNETCDF_DRIVER_NETCDF4
    /* netcdf4 driver is not enabled, skip */
    if (kind == 3 || kind == 4) {
        MPI_Finalize();
        return 0;
    }
#endif
    nerrs += pnetcdf_io(MPI_COMM_WORLD, filename, cmode, len);

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

