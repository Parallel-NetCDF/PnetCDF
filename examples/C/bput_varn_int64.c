/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example tests nonblocking buffered write varn APIs, including
 * ncmpi_bput_varn_longlong() and ncmpi_bput_varn(),
 * It first writes a sequence of requests with arbitrary array indices and
 * lengths to four variables of type NC_INT64, and reads back.
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicc -O2 -o i_varn_int64 i_varn_int64.c -lpnetcdf
 *    % mpiexec -n 4 ./i_varn_int64 /pvfs2/wkliao/testfile.nc
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *             Y = 10 ;
 *             X = 4 ;
 *    variables:
 *            int64 var0(Y, X) ;
 *            int64 var1(Y, X) ;
 *            int64 var2(Y, X) ;
 *            int64 var3(Y, X) ;
 *    data:
 *
 *     var0 =
 *      1, 0, 3, 0,
 *      1, 2, 3, 0,
 *      1, 2, 2, 0,
 *      3, 2, 1, 2,
 *      3, 1, 1, 3,
 *      0, 3, 1, 3,
 *      0, 3, 0, 3,
 *      2, 2, 0, 1,
 *      3, 2, 3, 1,
 *      3, 2, 3, 1 ;
 *
 *     var1 =
 *      2, 1, 0, 1,
 *      2, 3, 0, 1,
 *      2, 3, 3, 1,
 *      0, 3, 2, 3,
 *      0, 2, 2, 0,
 *      1, 0, 2, 0,
 *      1, 0, 1, 0,
 *      3, 3, 1, 2,
 *      0, 3, 0, 2,
 *      0, 3, 0, 2 ;
 *
 *     var2 =
 *      3, 2, 1, 2,
 *      3, 0, 1, 2,
 *      3, 0, 0, 2,
 *      1, 0, 3, 0,
 *      1, 3, 3, 1,
 *      2, 1, 3, 1,
 *      2, 1, 2, 1,
 *      0, 0, 2, 3,
 *      1, 0, 1, 3,
 *      1, 0, 1, 3 ;
 *
 *     var3 =
 *      0, 3, 2, 3,
 *      0, 1, 2, 3,
 *      0, 1, 1, 3,
 *      2, 1, 0, 1,
 *      2, 0, 0, 2,
 *      3, 2, 0, 2,
 *      3, 2, 3, 2,
 *      1, 1, 3, 0,
 *      2, 1, 2, 0,
 *      2, 1, 2, 0 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

#define NY 10
#define NX 4
#define NDIMS 2

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

#define ERRS(n,a) { \
    int _i; \
    for (_i=0; _i<(n); _i++) { \
        if ((a)[_i] != NC_NOERR) { \
            printf("Error at %s:%d : err[%d] %s\n", __FILE__,__LINE__, _i, \
                   ncmpi_strerror((a)[_i])); \
                   nerrs++; \
        } \
    } \
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

static MPI_Offset *** calloc_3D(size_t z, size_t y, size_t x)
{
    if (z*y*x == 0) return NULL;
    int _j, _k;
    MPI_Offset ***buf  = (MPI_Offset***) malloc(z     * sizeof(MPI_Offset**));
    MPI_Offset  **bufy = (MPI_Offset**)  malloc(z*y   * sizeof(MPI_Offset*));
    MPI_Offset   *bufx = (MPI_Offset*)   calloc(z*y*x,  sizeof(MPI_Offset));
    for (_k=0; _k<z; _k++, bufy+=y) {
        buf[_k] = bufy;
        for (_j=0; _j<y; _j++, bufx+=x)
            buf[_k][_j] = bufx;
    }
    return buf;
}

static void free_3D(MPI_Offset ***buf)
{
    free(buf[0][0]);
    free(buf[0]);
    free(buf);
}

static int check_contents(int ncid, int *varid)
{
    int i, j, err, nerrs=0;
    long long expected[4][NY*NX] = {{1, 0, 3, 0, 1, 2, 3, 0, 1, 2,
                                     2, 0, 3, 2, 1, 2, 3, 1, 1, 3,
                                     0, 3, 1, 3, 0, 3, 0, 3, 2, 2,
                                     0, 1, 3, 2, 3, 1, 3, 2, 3, 1},
                                    {2, 1, 0, 1, 2, 3, 0, 1, 2, 3,
                                     3, 1, 0, 3, 2, 3, 0, 2, 2, 0,
                                     1, 0, 2, 0, 1, 0, 1, 0, 3, 3,
                                     1, 2, 0, 3, 0, 2, 0, 3, 0, 2},
                                    {3, 2, 1, 2, 3, 0, 1, 2, 3, 0,
                                     0, 2, 1, 0, 3, 0, 1, 3, 3, 1,
                                     2, 1, 3, 1, 2, 1, 2, 1, 0, 0,
                                     2, 3, 1, 0, 1, 3, 1, 0, 1, 3},
                                    {0, 3, 2, 3, 0, 1, 2, 3, 0, 1,
                                     1, 3, 2, 1, 0, 1, 2, 0, 0, 2,
                                     3, 2, 0, 2, 3, 2, 3, 2, 1, 1,
                                     3, 0, 2, 1, 2, 0, 2, 1, 2, 0}};

    long long *r_buffer = (long long*) malloc(NY*NX*sizeof(long long));

    for (i=0; i<4; i++) {
        for (j=0; j<NY*NX; j++) r_buffer[j] = 99999;
        err = ncmpi_get_var_longlong_all(ncid, varid[i], r_buffer);
        ERR

        /* check if the contents of buf are expected */
        for (j=0; j<NY*NX; j++)
            if (r_buffer[j] != expected[i][j]) {
                printf("Expected file contents [%d][%d]=%lld, but got %lld\n",
                       i,j,expected[i][j],r_buffer[j]);
                nerrs++;
            }
    }
    free(r_buffer);
    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    char filename[256], *exec;
    int i, j, k, n, rank, nprocs, err, nerrs=0;
    int ncid, cmode, varid[4], dimid[2], nreqs, reqs[4], sts[4];
    long long *buffer[4];
    int num_segs[4], req_lens[4];
    MPI_Offset ***starts, ***counts;
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    exec = argv[0];

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

    if (nprocs != 4 && rank == 0 && verbose)
        printf("Warning: %s is intended to run on 4 processes\n",exec);

    /* set an MPI-IO hint to disable file offset alignment for fixed-size
     * variables */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_var_align_size", "1");

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    ERR

    MPI_Info_free(&info);

    /* create a global array of size NY * NX */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]); ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT64, NDIMS, dimid, &varid[0]); ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT64, NDIMS, dimid, &varid[1]); ERR
    err = ncmpi_def_var(ncid, "var2", NC_INT64, NDIMS, dimid, &varid[2]); ERR
    err = ncmpi_def_var(ncid, "var3", NC_INT64, NDIMS, dimid, &varid[3]); ERR
    if (nprocs < 4) { /* need 4 processes to fill the variables */
        err = ncmpi_set_fill(ncid, NC_FILL, NULL);
        ERR
    }
    err = ncmpi_enddef(ncid); ERR

    /* allocate space for starts and counts */
    starts = calloc_3D(4, 13, NDIMS);
    counts = calloc_3D(4, 13, NDIMS);

    n = rank % 4;
    num_segs[n] = 8; /* number of segments for this request */
    starts[n][0][0]=0; starts[n][0][1]=1; counts[n][0][0]=1; counts[n][0][1]=1;
    starts[n][1][0]=0; starts[n][1][1]=3; counts[n][1][0]=1; counts[n][1][1]=1;
    starts[n][2][0]=1; starts[n][2][1]=3; counts[n][2][0]=1; counts[n][2][1]=1;
    starts[n][3][0]=2; starts[n][3][1]=3; counts[n][3][0]=1; counts[n][3][1]=1;
    starts[n][4][0]=5; starts[n][4][1]=0; counts[n][4][0]=1; counts[n][4][1]=1;
    starts[n][5][0]=6; starts[n][5][1]=0; counts[n][5][0]=1; counts[n][5][1]=1;
    starts[n][6][0]=6; starts[n][6][1]=2; counts[n][6][0]=1; counts[n][6][1]=1;
    starts[n][7][0]=7; starts[n][7][1]=2; counts[n][7][0]=1; counts[n][7][1]=1;
    /* starts[n][][] n_counts[n][][] indicate the following: ("-" means skip)
              _  X  _  X
              _  _  _  X
              _  _  _  X
              _  _  _  _
              _  _  _  _
              X  _  _  _
              X  _  X  _
              _  _  X  _
              _  _  _  _
              _  _  _  _
     */
    n = (rank+1) % 4;
    num_segs[n] = 13; /* number of segments for this request */
    starts[n][ 0][0]=0;starts[n][ 0][1]=2;counts[n][ 0][0]=1;counts[n][ 0][1]=1;
    starts[n][ 1][0]=1;starts[n][ 1][1]=2;counts[n][ 1][0]=1;counts[n][ 1][1]=1;
    starts[n][ 2][0]=3;starts[n][ 2][1]=0;counts[n][ 2][0]=1;counts[n][ 2][1]=1;
    starts[n][ 3][0]=4;starts[n][ 3][1]=0;counts[n][ 3][0]=1;counts[n][ 3][1]=1;
    starts[n][ 4][0]=4;starts[n][ 4][1]=3;counts[n][ 4][0]=1;counts[n][ 4][1]=1;
    starts[n][ 5][0]=5;starts[n][ 5][1]=1;counts[n][ 5][0]=1;counts[n][ 5][1]=1;
    starts[n][ 6][0]=5;starts[n][ 6][1]=3;counts[n][ 6][0]=1;counts[n][ 6][1]=1;
    starts[n][ 7][0]=6;starts[n][ 7][1]=1;counts[n][ 7][0]=1;counts[n][ 7][1]=1;
    starts[n][ 8][0]=6;starts[n][ 8][1]=3;counts[n][ 8][0]=1;counts[n][ 8][1]=1;
    starts[n][ 9][0]=8;starts[n][ 9][1]=0;counts[n][ 9][0]=1;counts[n][ 9][1]=1;
    starts[n][10][0]=8;starts[n][10][1]=2;counts[n][10][0]=1;counts[n][10][1]=1;
    starts[n][11][0]=9;starts[n][11][1]=0;counts[n][11][0]=1;counts[n][11][1]=1;
    starts[n][12][0]=9;starts[n][12][1]=2;counts[n][12][0]=1;counts[n][12][1]=1;
    /* starts[n][][] counts[n][][] indicate the following pattern.
              _  _  X  _
              _  _  X  _
              _  _  _  _
              X  _  _  _
              X  _  _  X
              _  X  _  X
              _  X  _  X
              _  _  _  _
              X  _  X  _
              X  _  X  _
     */
    n = (rank+2) % 4;
    num_segs[n] = 7; /* number of segments for this request */
    starts[n][0][0]=1; starts[n][0][1]=1; counts[n][0][0]=1; counts[n][0][1]=1;
    starts[n][1][0]=2; starts[n][1][1]=1; counts[n][1][0]=1; counts[n][1][1]=2;
    starts[n][2][0]=3; starts[n][2][1]=1; counts[n][2][0]=1; counts[n][2][1]=1;
    starts[n][3][0]=3; starts[n][3][1]=3; counts[n][3][0]=1; counts[n][3][1]=1;
    starts[n][4][0]=7; starts[n][4][1]=0; counts[n][4][0]=1; counts[n][4][1]=2;
    starts[n][5][0]=8; starts[n][5][1]=1; counts[n][5][0]=1; counts[n][5][1]=1;
    starts[n][6][0]=9; starts[n][6][1]=1; counts[n][6][0]=1; counts[n][6][1]=1;
    /* starts[n][][] counts[n][][] indicate the following pattern.
              _  _  _  _
              _  X  _  _
              _  X  X  _
              _  X  _  X
              _  _  _  _
              _  _  _  _
              _  _  _  _
              X  X  _  _
              _  X  _  _
              _  X  _  _
     */
    n = (rank+3) % 4;
    num_segs[n] = 10; /* number of segments for this request */
    starts[n][0][0]=0; starts[n][0][1]=0; counts[n][0][0]=1; counts[n][0][1]=1;
    starts[n][1][0]=1; starts[n][1][1]=0; counts[n][1][0]=1; counts[n][1][1]=1;
    starts[n][2][0]=2; starts[n][2][1]=0; counts[n][2][0]=1; counts[n][2][1]=1;
    starts[n][3][0]=3; starts[n][3][1]=2; counts[n][3][0]=1; counts[n][3][1]=1;
    starts[n][4][0]=4; starts[n][4][1]=1; counts[n][4][0]=1; counts[n][4][1]=2;
    starts[n][5][0]=5; starts[n][5][1]=2; counts[n][5][0]=1; counts[n][5][1]=1;
    starts[n][6][0]=7; starts[n][6][1]=3; counts[n][6][0]=1; counts[n][6][1]=1;
    starts[n][7][0]=8; starts[n][7][1]=3; counts[n][7][0]=1; counts[n][7][1]=1;
    starts[n][8][0]=9; starts[n][8][1]=3; counts[n][8][0]=1; counts[n][8][1]=1;
     /*starts[n][][] counts[n][][] indicate the following pattern.
              X  _  _  _
              X  _  _  _
              X  _  _  _
              _  _  X  _
              _  X  X  _
              _  _  X  _
              _  _  _  _
              _  _  _  X
              _  _  _  X
              _  _  _  X
     */

    /* only rank 0, 1, 2, and 3 do I/O:
     * each of ranks 0 to 3 write 4 nonblocking requests */
    nreqs = 4;
    if (rank >= 4) nreqs = 0;

    /* bufsize must be max of data type converted before and after */
    MPI_Offset bufsize = 0;

    /* calculate length of each varn request, number of segments in each
     * varn request, and allocate write buffer */
    for (i=0; i<nreqs; i++) {
        req_lens[i] = 0; /* total length this request */
        for (j=0; j<num_segs[i]; j++) {
            MPI_Offset req_len=1;
            for (k=0; k<NDIMS; k++)
                req_len *= counts[i][j][k];
            req_lens[i] += req_len;
        }
        if (verbose) printf("req_lens[%d]=%d\n",i,req_lens[i]);

        /* allocate I/O buffer and initialize its contents */
        buffer[i] = (long long*) malloc(req_lens[i] * sizeof(long long));
        for (j=0; j<req_lens[i]; j++) buffer[i][j] = rank;

        bufsize += req_lens[i];
    }
    bufsize *= sizeof(long long);
    if (verbose) printf("%d: Attach buffer size %lld\n", rank, bufsize);

    /* give PnetCDF a space to buffer the nonblocking requests */
    if (bufsize > 0) {
        err = ncmpi_buffer_attach(ncid, bufsize); ERR
    }

    /* write using varn API */
    for (i=0; i<nreqs; i++) {
        err = ncmpi_bput_varn_longlong(ncid, varid[i], num_segs[i], starts[i],
                                       counts[i], buffer[i], &reqs[i]);
        ERR
    }
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    ERRS(nreqs, sts)

    /* check file contents */
    if (nprocs >= 4) nerrs += check_contents(ncid, varid);

    for (i=0; i<nreqs; i++) free(buffer[i]);

    /* test flexible put API, using a noncontiguous buftype */
    for (i=0; i<nreqs; i++) {
        MPI_Datatype buftype;
        MPI_Type_vector(req_lens[i], 1, 2, MPI_UNSIGNED, &buftype);
        MPI_Type_commit(&buftype);
        buffer[i] = (long long*) malloc(req_lens[i] * 2 * sizeof(long long));
        for (j=0; j<req_lens[i]*2; j++) buffer[i][j] = rank;

        err = ncmpi_bput_varn(ncid, varid[i], num_segs[i], starts[i],
                              counts[i], buffer[i], 1, buftype, &reqs[i]);
        ERR
        MPI_Type_free(&buftype);
    }
    err = ncmpi_wait_all(ncid, nreqs, reqs, sts);
    ERRS(nreqs, sts)

    /* check file contents */
    if (nprocs >= 4) nerrs += check_contents(ncid, varid);

    /* free the buffer space for bput */
    if (bufsize > 0) {
        err = ncmpi_buffer_detach(ncid); ERR
    }

    err = ncmpi_close(ncid);
    ERR

    for (i=0; i<nreqs; i++) free(buffer[i]);
    free_3D(starts);
    free_3D(counts);

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

