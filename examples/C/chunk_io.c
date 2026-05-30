/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to use the chunking and compression features of
 * PnetCDF to write a 3D record variable of integer data type in parallel. It
 * first defines a netCDF variable of size
 *       (NTIMES * nprocs) x (NY * nprocs) x NX
 * where NTIMES, NY, and NX are predefined constant.
 * The data partitioning pattern for write is along Y dimension. Each process
 * writes a subarray of size (NY * NX) per record.
 *
 *    To compile:
 *        mpicc -O2 chunk_io.c -o chunk_io -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * output netCDF file produced by this example program:
 *
 *    % mpiexec -n 4 ./chunk_io testfile.nc
 *
 *    % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            time = UNLIMITED ; // (0 currently)  <-- Not used anymore
 *            Y = 8 ;
 *            X = 10 ;
 *           _datablock_dim_0 = 65721 ;
 *           _datablock_dim_1 = 185 ;
 *           _datablock_dim_2 = 185 ;
 *           _datablock_dim_3 = 185 ;
 *           _datablock_dim_4 = 185 ;
 *           _datablock_dim_5 = 185 ;
 *           _datablock_dim_6 = 185 ;
 *           _datablock_dim_7 = 185 ;
 *    variables:
 *            int var ;
 *                     var:_ndim = 3 ;
 *                     var:_dimids = 0, 1, 2 ;
 *                     var:_datatype = 4 ;
 *                     var:_varkind = 1 ;
 *                     var:_chunkdim = 1, 2, 10 ;
 *                     var:_filter = 2 ;
 *                     var:_metaoffset = 4LL ;
 *            byte _datablock_0(_datablock_dim_0) ;
 *                     _datablock_0:_varkind = 2 ;
 *            byte _datablock_1(_datablock_dim_1) ;
 *                     _datablock_1:_varkind = 2 ;
 *            byte _datablock_2(_datablock_dim_2) ;
 *                     _datablock_2:_varkind = 2 ;
 *            byte _datablock_3(_datablock_dim_3) ;
 *                     _datablock_3:_varkind = 2 ;
 *            byte _datablock_4(_datablock_dim_4) ;
 *                     _datablock_4:_varkind = 2 ;
 *            byte _datablock_5(_datablock_dim_5) ;
 *                     _datablock_5:_varkind = 2 ;
 *            byte _datablock_6(_datablock_dim_6) ;
 *                     _datablock_6:_varkind = 2 ;
 *            byte _datablock_7(_datablock_dim_7) ;
 *                     _datablock_7:_varkind = 2 ;
 *
 *    // global attributes:
 *                    :_comressed = 1 ;
 *                    :_nwrite = 8 ;
 *                    :_recsize = 8LL ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>

#define NTIMES 2
#define NY 2
#define NX 10

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [-k format] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-k format] file format: 1 for CDF-1, 2 for CDF-2, 5 for CDF-5\n"
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

/*----< compress() >--------------------------------------------------------*/
static int
compress(MPI_Comm comm, char *filename, int cmode)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, varid, dimid[3];
    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[3], count[3];
    MPI_Info info;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_chunking", "enable");
    MPI_Info_set(info, "nc_chunk_default_filter", "zlib");

    /* chunking is supported only when MPI-IO driver is used */
    MPI_Info_set(info, "nc_pncio", "disable");

    /* the global array is (NTIMES * nprocs) x (NY * nprocs) x NX */

    /* set chunking (1st dimension should always be 1 for record variable) */

    int *buf = (int*) malloc(sizeof(int) * NY * NX);
    for (i=0; i<NY*NX; i++) buf[i] = rank + i;

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR

    /* define dimensions x and y */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "Y",    NY * nprocs,  &dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "X",    NX,           &dimid[2]); ERR

    /* define a 3D variable of integer type */
    err = ncmpi_def_var(ncid, "var", NC_INT, 3, dimid, &varid); ERR

    /* set chunking */
    int chunk_dim[3] = {1, NY, NX};
    err = ncmpi_var_set_chunk(ncid, varid, chunk_dim);; ERR

    /* use default filter */
    err = ncmpi_var_set_filter(ncid, varid, NC_FILTER_DEFLATE); ERR

    /* exit define mode */
    err = ncmpi_enddef(ncid); ERR

    /* set subarray start and count. Each rank write a subarray of size NY *
     * NX. Patitioning pattern is along Y dimension.
     */
                  start[1] = NY*rank;  start[2] = 0;
    count[0] = 1; count[1] = NY;       count[2] = NX;

    /* write one time record at a time */
    for (i=0; i<NTIMES*nprocs; i++) {
        start[0] = i;
        err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf);
        ERR
    }

    /* check the current record dimension size */
    MPI_Offset dim_len;
    err = ncmpi_inq_dimlen(ncid, 0, &dim_len); ERR
    if (verbose && rank == 0)
        printf("Time dimension length (expect %lld) and got %lld\n",
                NTIMES * nprocs, dim_len);

    err = ncmpi_close(ncid); ERR

    MPI_Info_free(&info);

    free(buf);

    return nerrs;
}

/*----< decompress() >-------------------------------------------------------*/
static int
decompress(MPI_Comm comm, char *filename)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, varid, dimid, filter, chunk_dim[3];
    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[3], count[3];
    MPI_Info info;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_chunking", "enable");

    /* chunking is supported only when MPI-IO driver is used */
    MPI_Info_set(info, "nc_pncio", "disable");

    /* the global array is (NTIMES * nprocs) x (NY * nprocs) * NX */

    /* open the file for reading ----------------------------------------*/
    err = ncmpi_open(comm, filename, NC_NOWRITE, info, &ncid); ERR

    err = ncmpi_inq_varid(ncid, "var", &varid); ERR

    /* check the current record dimension size */
    MPI_Offset dim_len;
    err = ncmpi_inq_unlimdim(ncid, &dimid); ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &dim_len); ERR
    if (verbose && rank == 0)
        printf("Time dimension length = %lld\n", dim_len);

    /* get chunking */
    err = ncmpi_var_get_chunk(ncid, varid, chunk_dim);; ERR
    if (verbose && rank == 0)
        printf("chunk_dim[3]=%d %d %d\n",
               chunk_dim[0],chunk_dim[1],chunk_dim[2]);

    /* get filter */
    err = ncmpi_var_get_filter(ncid, varid, &filter); ERR
    if (verbose && rank == 0)
        printf("filter is %s\n", (filter == NC_FILTER_DEFLATE) ?
               "NC_FILTER_DEFLATE": (filter == NC_FILTER_SZ) ?
               "NC_FILTER_SZ" : "UNKNOWN");

    /* set subarray start and count. Each rank read a whole record at a time
     * for NTIMES times. Each process reads different records.
     */
    start[0] = rank; start[1] = 0;         start[2] = 0;
    count[0] = 1;    count[1] = NY*nprocs; count[2] = NX;

    if (verbose)
        printf("%d: start=%lld %lld %lld count=%lld %lld %lld\n", rank,
               start[0], start[1], start[2], count[0], count[1], count[2]);

    int *buf;
    buf = (int*) malloc(sizeof(int) * count[0] * count[1] * count[2]);
    for (j=0; j<count[0]*count[1]*count[2]; j++) buf[j] = -1;

    /* Each rank reads NTIMES records, one at a time. */
    for (i=0; i<NTIMES; i++) {
        err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); ERR
        start[0] += nprocs;
    }

    err = ncmpi_close(ncid); ERR

    MPI_Info_free(&info);

    free(buf);

    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, rank, kind=0, cmode=0, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 0;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqk:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'k': kind = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (verbose && rank == 0) printf("%s: example of using compression\n",__FILE__);

    switch (kind) {
        case(2): cmode = NC_64BIT_OFFSET; break;
        case(5): cmode = NC_64BIT_DATA;   break;
        default: cmode = 0;
    }

    nerrs += compress(MPI_COMM_WORLD, filename, cmode);

    nerrs += decompress(MPI_COMM_WORLD, filename);

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

