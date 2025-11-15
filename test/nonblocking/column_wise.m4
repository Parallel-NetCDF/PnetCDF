/*********************************************************************
 *
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example tests a number of nonblocking API calls, each writes a single
 * column of a 2D integer array. Each process writes NX columns and any two
 * consecutive columns are of nprocs columns distance apart from each other. In
 * this case, the fileview of each process interleaves with all other processes.
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
 *    % m4 column_wise.m4 > column_wise.c
 *    % mpicc -O2 -o column_wise column_wise.c -lpnetcdf
 *    % mpiexec -l -n 4 ./column_wise -l 4 /pvfs2/wkliao/testfile.nc
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
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13,
 *      10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13, 10, 11, 12, 13 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef HAVE_CONFIG_H
#include <config.h> /* output of 'configure' */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <unistd.h> /* getopt() */

#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 10
#define NX 70

typedef char text;

include(`foreach.m4')dnl
include(`utils.m4')dnl

define(`TEST_DATA_TYPE',`test_column_wise_$1(filename, cdf_formats[fmt], coll_io, len)')

define(`TEST_COLUMN_WISE',`dnl
static
int test_column_wise_$1(char *filename, int cdf, int coll_io, int len)
{
    int i, j, nerrs=0, rank, nprocs, err, myNX, G_NX, myOff, num_reqs;
    int ncid, cmode, varid, dimid[2], *reqs, *sts;
    $1 **buf;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    cmode = NC_CLOBBER;
    if (cdf == NC_FORMAT_CDF2)
        cmode |= NC_64BIT_OFFSET;
    else if (cdf == NC_FORMAT_CDF5)
        cmode |= NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* the global array is NY * (len * nprocs) */
    G_NX  = len * nprocs;
    myOff = len * rank;
    myNX  = len;

    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", G_NX, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_TYPE($1), 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* First, fill the entire array with zeros, using a blocking I/O.
       Every process writes a subarray of size NY * myNX */
    buf    = ($1**) malloc(sizeof($1*) * myNX);
    buf[0] = ($1*)  calloc(NY * myNX, sizeof($1));
    start[0] = 0;   start[1] = myOff;
    count[0] = NY;  count[1] = myNX;
    err = ncmpi_put_vara_`$1'_all(ncid, varid, start, count, buf[0]); CHECK_ERR
    free(buf[0]);

    /* When using burst buffering, flush the log to prevent new value being
     * skipped due to overlaping domain
     */
    err = ncmpi_flush(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* initialize the buffer with rank ID. Also make the case interesting,
       by allocatsing buffersd separately */
    for (i=0; i<myNX; i++) {
        buf[i] = ($1*) malloc(sizeof($1) * NY);
        for (j=0; j<NY; j++) buf[i][j] = ($1)rank+10;
    }

    reqs = (int*) malloc(sizeof(int) * myNX);
    sts  = (int*) malloc(sizeof(int) * myNX);

    /* each proc writes myNX single columns of the 2D array */
    start[0]  = 0;   start[1] = rank;
    count[0]  = NY;  count[1] = 1;

    num_reqs = 0;
    for (i=0; i<myNX; i++) {
        err = ncmpi_iput_vara_$1(ncid, varid, start, count, buf[i],
                                 &reqs[num_reqs++]); CHECK_ERR
        start[1] += nprocs;
    }

    /* try re-order the request list */
    for (i=0; i<myNX/2; i++) {
        int tmp = reqs[2*i];
        reqs[2*i] = reqs[2*i+1];
        reqs[2*i+1] = tmp;
    }

    /* test cancelling requests and see if the user write buffer is properly
     * byte-swapped back to it original form. To do this test, NY must be
     * changed to use a number > NC_BYTE_SWAP_BUFFER_SIZE/sizeof(int), say
     * 1025
     */
    err = ncmpi_cancel(ncid, num_reqs, reqs, sts); CHECK_ERR

    /* check if write buffer contents have been altered after cancelling */
    for (i=0; i<myNX; i++) {
        for (j=0; j<NY; j++) {
            if (buf[i][j] != ($1)rank+10) {
                printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=IFMT($1)\n",
                       __LINE__,__FILE__,i,j,buf[i][j]);
                nerrs++;
                i = myNX; break;
            }
        }
    }

    /* post iput requests again */
    start[1] = rank;
    num_reqs = 0;
    for (i=0; i<myNX; i++) {
        err = ncmpi_iput_vara_$1(ncid, varid, start, count, buf[i],
                                 &reqs[num_reqs++]); CHECK_ERR
        start[1] += nprocs;
    }

    /* try re-order the request list */
    for (i=0; i<myNX/2; i++) {
        int tmp = reqs[2*i];
        reqs[2*i] = reqs[2*i+1];
        reqs[2*i+1] = tmp;
    }

    if (coll_io)
        err = ncmpi_wait_all(ncid, num_reqs, reqs, sts);
    else
        err = ncmpi_wait(ncid, num_reqs, reqs, sts);
    CHECK_ERR

    /* check if write buffer contents have been altered after wait */
    for (i=0; i<myNX; i++) {
        for (j=0; j<NY; j++) {
            if (buf[i][j] != ($1)rank+10) {
                printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=IFMT($1)\n",
                       __LINE__,__FILE__,i,j,buf[i][j]);
                nerrs++;
                i = myNX; break;
            }
        }
    }

    /* check status of all requests */
    for (i=0; i<num_reqs; i++) {
        if (reqs[i] != NC_REQ_NULL) { /* add in PnetCDF v1.7.0 */
            printf("Error at line %d in %s: request ID %d fails to be set to NC_REQ_NULL\n",__LINE__,__FILE__,i);
            nerrs++;
            break;
        }
        if (sts[i] != NC_NOERR) {
            printf("Error at line %d in %s: nonblocking write fails on request %d (%s)\n",
                   __LINE__,__FILE__,i, ncmpi_strerror(sts[i]));
            nerrs++;
            break;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_varid(ncid, "var", &varid); CHECK_ERR

    /* read back using the same access pattern */
    for (i=0; i<myNX; i++)
        for (j=0; j<NY; j++) buf[i][j] = ($1)99;

    /* each proc reads myNX single columns of the 2D array */
    start[0]  = 0;   start[1] = rank;
    count[0]  = NY;  count[1] = 1;

    num_reqs = 0;
    for (i=0; i<myNX; i++) {
        err = ncmpi_iget_vara_$1(ncid, varid, start, count, buf[i],
                                 &reqs[num_reqs++]); CHECK_ERR
        start[1] += nprocs;
    }
    /* this test is to see if cancelling free up all the internal malloc */
    err = ncmpi_cancel(ncid, num_reqs, reqs, sts); CHECK_ERR

    /* post iget requests again */
    start[1] = rank;
    num_reqs = 0;
    for (i=0; i<myNX; i++) {
        err = ncmpi_iget_vara_$1(ncid, varid, start, count, buf[i],
                                 &reqs[num_reqs++]); CHECK_ERR
        start[1] += nprocs;
    }

    if (coll_io)
        err = ncmpi_wait_all(ncid, num_reqs, reqs, sts);
    else
        err = ncmpi_wait(ncid, num_reqs, reqs, sts);
    CHECK_ERR

    /* check status of all requests */
    for (i=0; i<num_reqs; i++)
        if (sts[i] != NC_NOERR) {
            printf("Error at line %d in %s: nonblocking write fails on request %d (%s)\n",
                   __LINE__,__FILE__,i, ncmpi_strerror(sts[i]));
            nerrs++;
            break;
        }

    for (i=0; i<myNX; i++) {
        for (j=0; j<NY; j++) {
            $1 expected = ($1)rank+10;
            if (buf[i][j] != expected) {
                printf("Error at line %d in %s: expect buf[%d][%d]=IFMT($1) but got IFMT($1)\n",
                       __LINE__,__FILE__,i,j,expected,buf[i][j]);
                nerrs++;
                i = myNX; break;
            }
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    free(sts);
    free(reqs);
    for (i=0; i<myNX; i++) free(buf[i]);
    free(buf);

    return nerrs;
}'
)dnl

TEST_COLUMN_WISE(text)
TEST_COLUMN_WISE(schar)
TEST_COLUMN_WISE(uchar)
TEST_COLUMN_WISE(short)
TEST_COLUMN_WISE(ushort)
TEST_COLUMN_WISE(int)
TEST_COLUMN_WISE(uint)
TEST_COLUMN_WISE(float)
TEST_COLUMN_WISE(double)
TEST_COLUMN_WISE(longlong)
TEST_COLUMN_WISE(ulonglong)

#define FILE_NAME "testfile.nc"

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTIONS]...[filename]\n"
    "       [-h] Print help\n"
    "       [-l num]: X dimension size of local array\n"
    "       [filename]: output netCDF file name (default: %s)\n";
    fprintf(stderr, help, argv0, FILE_NAME);
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, fmt, nerrs=0, rank, err, len, coll_io;
    int cdf_formats[3]={NC_FORMAT_CLASSIC, NC_FORMAT_CDF2, NC_FORMAT_CDF5};

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    len = NX;
    while ((i = getopt(argc, argv, "hl:")) != EOF)
        switch(i) {
            case 'l': len = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, FILE_NAME);
    else                      snprintf(filename, 256, "%s", argv[optind]);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for iput/iget interleaved access ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    for (coll_io=0; coll_io<2; coll_io++) {
        for (fmt=0; fmt<3; fmt++) {
            nerrs = TEST_DATA_TYPE(text);   if (nerrs > 0) goto err_out;
            nerrs = TEST_DATA_TYPE(schar);  if (nerrs > 0) goto err_out;
            nerrs = TEST_DATA_TYPE(short);  if (nerrs > 0) goto err_out;
            nerrs = TEST_DATA_TYPE(int);    if (nerrs > 0) goto err_out;
            nerrs = TEST_DATA_TYPE(float);  if (nerrs > 0) goto err_out;
            nerrs = TEST_DATA_TYPE(double); if (nerrs > 0) goto err_out;
            if (cdf_formats[fmt] == NC_FORMAT_CDF5) {
                nerrs = TEST_DATA_TYPE(uchar);     if (nerrs > 0) goto err_out;
                nerrs = TEST_DATA_TYPE(ushort);    if (nerrs > 0) goto err_out;
                nerrs = TEST_DATA_TYPE(uint);      if (nerrs > 0) goto err_out;
                nerrs = TEST_DATA_TYPE(longlong);  if (nerrs > 0) goto err_out;
                nerrs = TEST_DATA_TYPE(ulonglong); if (nerrs > 0) goto err_out;
            }
        }
    }

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has "OFFFMT" bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

err_out:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

