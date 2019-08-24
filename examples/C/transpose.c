/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/*
 *    This example shows how to use varm API to write six 3D integer array
 *    variables into a file. Each variable in the file is a dimensional
 *    transposed array from the one stored in memory. In memory, a 3D array is
 *    partitioned among all processes in a block-block-block fashion and in
 *    ZYX (i.e. C) order. The dimension structures of the transposed six
 *    arrays are
 *       int ZYX_var(Z, Y, X) ;     ZYX -> ZYX
 *       int ZXY_var(Z, X, Y) ;     ZYX -> ZXY
 *       int YZX_var(Y, Z, X) ;     ZYX -> YZX
 *       int YXZ_var(Y, X, Z) ;     ZYX -> YXZ
 *       int XZY_var(X, Z, Y) ;     ZYX -> XZY
 *       int XYZ_var(X, Y, Z) ;     ZYX -> XYZ
 *
 *    To compile:
 *        mpicc -O2 transpose.c -o transpose -lpnetcdf
 *    To run:
 *        mpiexec -n num_processes ./transpose [filename] [len]
 *    where len decides the size of local array, which is len x len+1 x len+2.
 *    So, each variable is of size len*(len+1)*(len+2) * nprocs * sizeof(int)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

#define NDIMS 3

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
    int i, j, k, rank, nprocs, ncid, bufsize, err, nerrs=0;
    int *buf, psizes[NDIMS], dimids[NDIMS], dimidsT[NDIMS];
    int XYZ_id, XZY_id, YZX_id, YXZ_id, ZYX_id, ZXY_id;
    MPI_Offset gsizes[NDIMS], starts[NDIMS], counts[NDIMS], imap[NDIMS];
    MPI_Offset startsT[NDIMS], countsT[NDIMS];
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
     * starts[] */
    int lower_dims=1;
    for (i=NDIMS-1; i>=0; i--) {
        starts[i] = rank / lower_dims % psizes[i];
        lower_dims *= psizes[i];
    }
    if (verbose) {
        sprintf(str, "proc %d: dim rank= ", rank);
        for (i=0; i<NDIMS; i++) sprintf(str+strlen(str), "%lld ",starts[i]);
        printf("%s\n",str);
    }

    bufsize = 1;
    for (i=0; i<NDIMS; i++) {
        gsizes[i]  = (MPI_Offset)(len + i) * psizes[i]; /* global array size */
        starts[i] *= (MPI_Offset)(len + i);             /* start indices */
        counts[i]  = (MPI_Offset)(len + i);             /* array elements */
        bufsize   *= (len + i);
    }

    /* allocate buffer and initialize with contiguous numbers */
    buf = (int *) malloc(bufsize * sizeof(int));
    for (k=0; k<counts[0]; k++)
    for (j=0; j<counts[1]; j++)
    for (i=0; i<counts[2]; i++)
        buf[k*counts[1]*counts[2] +
                      j*counts[2] + i] = (starts[0]+k)*gsizes[1]*gsizes[2]
                                       + (starts[1]+j)*gsizes[2]
                                       + (starts[2]+i);

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
    for (i=0; i<NDIMS; i++) {
        sprintf(str, "%c", 'Z'-i);
        err = ncmpi_def_dim(ncid, str, gsizes[i], &dimids[i]);
        ERR
    }

    /* define variable with no transposed file layout: ZYX */
    err = ncmpi_def_var(ncid, "ZYX_var", NC_INT, NDIMS, dimids, &ZYX_id);
    ERR

    /* define variable with transposed file layout: ZYX -> ZXY */
    dimidsT[0] = dimids[0]; dimidsT[1] = dimids[2]; dimidsT[2] = dimids[1];
    err = ncmpi_def_var(ncid, "ZXY_var", NC_INT, NDIMS, dimidsT, &ZXY_id);
    ERR

    /* define variable with transposed file layout: ZYX -> YZX */
    dimidsT[0] = dimids[1]; dimidsT[1] = dimids[0]; dimidsT[2] = dimids[2];
    err = ncmpi_def_var(ncid, "YZX_var", NC_INT, NDIMS, dimidsT, &YZX_id);
    ERR

    /* define variable with transposed file layout: ZYX -> YXZ */
    dimidsT[0] = dimids[1]; dimidsT[1] = dimids[2]; dimidsT[2] = dimids[0];
    err = ncmpi_def_var(ncid, "YXZ_var", NC_INT, NDIMS, dimidsT, &YXZ_id);
    ERR

    /* define variable with transposed file layout: ZYX -> XZY */
    dimidsT[0] = dimids[2]; dimidsT[1] = dimids[0]; dimidsT[2] = dimids[1];
    err = ncmpi_def_var(ncid, "XZY_var", NC_INT, NDIMS, dimidsT, &XZY_id);
    ERR

    /* define variable with transposed file layout: ZYX -> XYZ */
    dimidsT[0] = dimids[2]; dimidsT[1] = dimids[1]; dimidsT[2] = dimids[0];
    err = ncmpi_def_var(ncid, "XYZ_var", NC_INT, NDIMS, dimidsT, &XYZ_id);
    ERR

    /* exit the define mode */
    err = ncmpi_enddef(ncid);
    ERR

    /* write the whole variable in file: ZYX */
    err = ncmpi_put_vara_int_all(ncid, ZYX_id, starts, counts, buf);
    ERR

    /* ZYX -> ZXY: */
    imap[1] = 1; imap[2] = counts[2]; imap[0] = counts[1]*counts[2];
    startsT[0] = starts[0]; startsT[1] = starts[2]; startsT[2] = starts[1];
    countsT[0] = counts[0]; countsT[1] = counts[2]; countsT[2] = counts[1];
    /* write the transposed variable */
    err = ncmpi_put_varm_int_all(ncid, ZXY_id, startsT, countsT, NULL, imap, buf);
    ERR

    /* ZYX -> YZX: */
    imap[2] = 1; imap[0] = counts[2]; imap[1] = counts[1]*counts[2];
    startsT[0] = starts[1]; startsT[1] = starts[0]; startsT[2] = starts[2];
    countsT[0] = counts[1]; countsT[1] = counts[0]; countsT[2] = counts[2];
    /* write the transposed variable */
    err = ncmpi_put_varm_int_all(ncid, YZX_id, startsT, countsT, NULL, imap, buf);
    ERR

    /* ZYX -> YXZ: */
    imap[1] = 1; imap[0] = counts[2]; imap[2] = counts[1]*counts[2];
    startsT[0] = starts[1]; startsT[1] = starts[2]; startsT[2] = starts[0];
    countsT[0] = counts[1]; countsT[1] = counts[2]; countsT[2] = counts[0];
    /* write the transposed variable */
    err = ncmpi_put_varm_int_all(ncid, YXZ_id, startsT, countsT, NULL, imap, buf);
    ERR

    /* ZYX -> XZY: */
    imap[0] = 1; imap[2] = counts[2]; imap[1] = counts[1]*counts[2];
    startsT[0] = starts[2]; startsT[1] = starts[0]; startsT[2] = starts[1];
    countsT[0] = counts[2]; countsT[1] = counts[0]; countsT[2] = counts[1];
    /* write the transposed variable */
    err = ncmpi_put_varm_int_all(ncid, XZY_id, startsT, countsT, NULL, imap, buf);
    ERR

    /* ZYX -> XYZ: */
    imap[0] = 1; imap[1] = counts[2]; imap[2] = counts[1]*counts[2];
    startsT[0] = starts[2]; startsT[1] = starts[1]; startsT[2] = starts[0];
    countsT[0] = counts[2]; countsT[1] = counts[1]; countsT[2] = counts[0];
    /* write the transposed variable */
    err = ncmpi_put_varm_int_all(ncid, XYZ_id, startsT, countsT, NULL, imap, buf);
    ERR

    /* close the file */
    err = ncmpi_close(ncid);
    ERR

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

    len = (len <= 0) ? 10 : len;

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

