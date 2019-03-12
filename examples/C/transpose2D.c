/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/*
 *    This example shows how to use varm API to write a 2D integer array
 *    variable into a file. The variable in the file is a dimensional
 *    transposed array from the one stored in memory. In memory, a 2D array is
 *    partitioned among all processes in a block-block fashion in YX (i.e.
 *    row-major) order. The dimension structure of the transposed array is
 *
 *       int YX_var(Y, X) ;
 *       int XY_var(X, Y) ;
 *
 *    To compile:
 *        mpicc -O2 transpose2D.c -o transpose2D -lpnetcdf
 *    To run:
 *        mpiexec -n num_processes ./transpose [filename] [len]
 *    where len decides the size of local array, which is len x len+1.
 *    So, each variable is of size len*(len+1) * nprocs * sizeof(int)
 *
 *    % mpiexec -n 4 ./transpose2D
 *    % ncdump testfile.nc
 *    netcdf testfile {
 *    dimensions:
 *             Y = 4 ;
 *             X = 6 ;
 *    variables:
 *            int YX_var(Y, X) ;
 *            int XY_var(X, Y) ;
 *    data:
 *
 *     YX_var =
 *      0, 1, 2, 3, 4, 5,
 *      6, 7, 8, 9, 10, 11,
 *      12, 13, 14, 15, 16, 17,
 *      18, 19, 20, 21, 22, 23 ;
 *
 *     XY_var =
 *      0, 6, 12, 18,
 *      1, 7, 13, 19,
 *      2, 8, 14, 20,
 *      3, 9, 15, 21,
 *      4, 10, 16, 22,
 *      5, 11, 17, 23 ;
 *    }
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

#define NDIMS 2

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [-k format] [-l len] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
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
    char str[512];
    int i, j, rank, nprocs, ncid, bufsize, err, nerrs=0;
    int *buf, psizes[NDIMS], dimids[NDIMS], dimidsT[NDIMS];
    int XY_id, YX_id;
    MPI_Offset gsizes[NDIMS], start[NDIMS], count[NDIMS], imap[NDIMS];
    MPI_Offset startT[NDIMS], countT[NDIMS];
    MPI_Info info;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    for (i=0; i<NDIMS; i++)
        psizes[i] = 0;

    /* calculate number of processes along each dimension */
    MPI_Dims_create(nprocs, NDIMS, psizes);
    if (verbose && rank == 0) {
        sprintf(str, "psizes= ");
        for (i=0; i<NDIMS; i++) sprintf(str+strlen(str), "%d ",psizes[i]);
        printf("%s\n",str);
    }

    /* for each MPI rank, find its local rank IDs along each dimension in
     * start[] */
    int lower_dims=1;
    for (i=NDIMS-1; i>=0; i--) {
        start[i] = rank / lower_dims % psizes[i];
        lower_dims *= psizes[i];
    }
    if (verbose) {
        sprintf(str, "proc %d: dim rank= ", rank);
        for (i=0; i<NDIMS; i++) sprintf(str+strlen(str), "%lld ",start[i]);
        printf("%s\n",str);
    }

    bufsize = 1;
    for (i=0; i<NDIMS; i++) {
        gsizes[i]  = (MPI_Offset)(len + i) * psizes[i]; /* global array size */
        start[i] *= (MPI_Offset)(len + i);              /* start indices */
        count[i]  = (MPI_Offset)(len + i);              /* array elements */
        bufsize   *= (len + i);                         /* local buffer size */
    }

    /* allocate buffer and initialize buffer with contiguous numbers, in YX layout */
    buf = (int *) malloc(bufsize * sizeof(int));
    for (i=0; i<count[0]; i++)
    for (j=0; j<count[1]; j++)
        buf[i*count[1] + j] = (start[0]+i)*gsizes[1] + (start[1]+j);

    /* set an MPI-IO hint to disable file offset alignment for fixed-size
     * variables */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_var_align_size", "1");

    /* create the file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_create() file %s (%s)\n",
        __LINE__,__FILE__,filename,ncmpi_strerror(err));
        MPI_Abort(comm, -1);
        exit(1);
    }

    MPI_Info_free(&info);

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "Y", gsizes[0], &dimids[0]); ERR
    err = ncmpi_def_dim(ncid, "X", gsizes[1], &dimids[1]); ERR

    /* define variable with no transposed file layout: YX */
    err = ncmpi_def_var(ncid, "YX_var", NC_INT, NDIMS, dimids, &YX_id); ERR

    /* define variable with transposed file layout: XY */
    dimidsT[0] = dimids[1];
    dimidsT[1] = dimids[0];
    err = ncmpi_def_var(ncid, "XY_var", NC_INT, NDIMS, dimidsT, &XY_id); ERR

    /* exit the define mode */
    err = ncmpi_enddef(ncid); ERR

    /* write the whole variable in file: YX_var */
    err = ncmpi_put_vara_int_all(ncid, YX_id, start, count, buf); ERR

    /* transpose YX -> XY */
    imap[0] = 1;
    imap[1] = count[1];
    startT[0] = start[1];
    startT[1] = start[0];
    countT[0] = count[1];
    countT[1] = count[0];
    /* write the transposed variable XY_var */
    err = ncmpi_put_varm_int_all(ncid, XY_id, startT, countT, NULL, imap, buf);
    ERR

    /* erase buffer contents to prepare read */
    for (i=0; i<countT[0]; i++)
    for (j=0; j<countT[1]; j++)
        buf[i*countT[1] + j] = -1;

    /* read back the whole variable XY_var */
    err = ncmpi_get_vara_int_all(ncid, XY_id, startT, countT, buf); ERR

    /* check the contents */
    for (i=0; i<countT[0]; i++)
    for (j=0; j<countT[1]; j++) {
        int expect = (start[1]+i) + gsizes[1] * (start[0]+j);
        if (buf[i*countT[1] + j] !=  expect) {
            printf("Error at %s:%d: expect buf[%lld]=%d but got %d\n",
                   __FILE__,__LINE__,i*countT[1]+j,expect,buf[i*count[1]+j]);
            nerrs++;
            i = count[1]; /* break loop i */
            break;
        }
    }

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
    int i, rank, len=0, cmode=0, kind=0, nerrs=0;

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

    len = (len <= 0) ? 2 : len;

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

