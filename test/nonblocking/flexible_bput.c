/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This example tests PnetCDF nonblocking buffered flexible varm API, i.e.
 * ncmpi_bput_varm() to write a 2D array double variable of size NY x NX*nproc
 * in parallel. In particular, we use a noncontiguous buffer type, a
 * noncontiguous imap[], and integer type in memory and double in file that
 * require a type conversion.
 *
 * The data partitioning patterns on the variable is column-wise.
 * The local buffer has ghost cells surrounded along both dimensions.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -O2 -o flexible_bput flexible_bput.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./flexible_bput -l 4 /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            Y = 6 ;
 *            X = 16 ;
 *    variables:
 *            double var(Y, X) ;
 *    data:
 *
 *    var =
 *      0,  6, 12, 18, 0,  6, 12, 18, 0,  6, 12, 18, 0,  6, 12, 18,
 *      1,  7, 13, 19, 1,  7, 13, 19, 1,  7, 13, 19, 1,  7, 13, 19,
 *      2,  8, 14, 20, 2,  8, 14, 20, 2,  8, 14, 20, 2,  8, 14, 20,
 *      3,  9, 15, 21, 3,  9, 15, 21, 3,  9, 15, 21, 3,  9, 15, 21,
 *      4, 10, 16, 22, 4, 10, 16, 22, 4, 10, 16, 22, 4, 10, 16, 22,
 *      5, 11, 17, 23, 5, 11, 17, 23, 5, 11, 17, 23, 5, 11, 17, 23 ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <unistd.h> /* getopt() */
#include <assert.h>

#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 6
#define NX 70
#define GHOST 2

#define INIT_PUT_BUF(buf) \
    for (i=0; i<array_of_sizes[0]; i++) { \
        for (j=0; j<array_of_sizes[1]; j++) { \
            if (i < GHOST || GHOST+array_of_subsizes[0] <= i || \
                j < GHOST || GHOST+array_of_subsizes[1] <= j) \
                buf[i][j] = -1; \
            else \
                buf[i][j] = (i-GHOST)*array_of_subsizes[1]+(j-GHOST); \
        } \
    }

#define CHECK_PUT_BUF(buf) \
    for (i=0; i<array_of_sizes[0]; i++) { \
        for (j=0; j<array_of_sizes[1]; j++) { \
            if (i < GHOST || GHOST+array_of_subsizes[0] <= i || \
                j < GHOST || GHOST+array_of_subsizes[1] <= j) { \
                if (buf[i][j] != -1) { \
                    printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%f\n", \
                           __LINE__,__FILE__,i,j,(double)buf[i][j]); \
                    nerrs++; \
                    goto err_out; \
                } \
            } \
            else { \
                if (buf[i][j] != (i-GHOST)*array_of_subsizes[1]+(j-GHOST)) { \
                    printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%f\n", \
                           __LINE__,__FILE__,i,j,(double)buf[i][j]); \
                    nerrs++; \
                    goto err_out; \
                } \
            } \
        } \
    }

#define INIT_GET_BUF(buf) \
    for (i=0; i<array_of_sizes[0]; i++) \
        for (j=0; j<array_of_sizes[1]; j++) \
            buf[i][j] = -2;

#define CHECK_GET_BUF(buf) \
    for (i=0; i<array_of_sizes[0]; i++) { \
        for (j=0; j<array_of_sizes[1]; j++) { \
            if (i < GHOST || GHOST+array_of_subsizes[0] <= i || \
                j < GHOST || GHOST+array_of_subsizes[1] <= j) { \
                if (buf[i][j] != -2) { \
                    printf("Error at line %d in %s: expect buffer[%d][%d] to be %d but got %f\n", \
                           __LINE__,__FILE__,i,j,-2,(double)buf[i][j]); \
                    nerrs++; \
                    goto err_out; \
                } \
            } \
            else { \
                int exp = (i-GHOST)*array_of_subsizes[1]+(j-GHOST); \
                if (buf[i][j] != exp) { \
                    printf("Error at line %d in %s: expect buffer[%d][%d] to be %d but got %f\n", \
                           __LINE__,__FILE__,i,j,exp,(double)buf[i][j]); \
                    nerrs++; \
                    goto err_out; \
                } \
            } \
        } \
    }

static
int tst_run(const char *filename, int coll_io, int len)
{
    int i, j, rank, nprocs, err, nerrs=0, req, status;
    int ncid, cmode, varid, dimid[2];
    int array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    int    **buf_int;
    double **buf_dbl;
    MPI_Offset start[2], count[2], stride[2], imap[2];
    MPI_Datatype  subarray;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    buf_int = (int**) malloc(sizeof(int*) * (len+2*GHOST));
    buf_int[0] = (int*) malloc(sizeof(int) * (len+2*GHOST) * (NY+2*GHOST));
    for (i=1; i<(len+2*GHOST); i++)
        buf_int[i] = buf_int[i-1] + (NY+2*GHOST);

    buf_dbl = (double**) malloc(sizeof(double*) * (len+2*GHOST));
    buf_dbl[0] = (double*) malloc(sizeof(double) * (len+2*GHOST) * (NY+2*GHOST));
    for (i=1; i<(len+2*GHOST); i++)
        buf_dbl[i] = buf_dbl[i-1] + (NY+2*GHOST);

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* define 2 dimensions */
    err = ncmpi_def_dim(ncid, "Y", NY,         &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", len*nprocs, &dimid[1]); CHECK_ERR

    /* define a variable of size NY * (len * nprocs) */
    err = ncmpi_def_var(ncid, "var", NC_DOUBLE, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

     start[0] = 0;  start[1] = len * rank;
     count[0] = NY; count[1] = len;
    stride[0] = 1; stride[1] = 1;
      imap[0] = 1;   imap[1] = NY; /* would be {len, 1} if not transposing */

    /* var is partitioned along X dimension in a matrix transported way */
    array_of_sizes[0]    = len + 2*GHOST;
    array_of_sizes[1]    = NY + 2*GHOST;
    array_of_subsizes[0] = len;
    array_of_subsizes[1] = NY;
    array_of_starts[0]   = GHOST;
    array_of_starts[1]   = GHOST;
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C, MPI_INT, &subarray);
    MPI_Type_commit(&subarray);

    /* calling a nonblocking bput_varm flexible API -------------------------*/
    /* initiate put buffer contents */
    INIT_PUT_BUF(buf_int)

    MPI_Offset bufsize = sizeof(double);
    for (i=0; i<2; i++) bufsize *= count[i];
    err = ncmpi_buffer_attach(ncid, bufsize); CHECK_ERR

    err = ncmpi_bput_varm(ncid, varid, start, count, stride, imap, buf_int[0],
                          1, subarray, &req);
    CHECK_ERR
    /* check if the contents of put buffer are altered */
    CHECK_PUT_BUF(buf_int)

    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, &status);
    else
        err = ncmpi_wait(ncid, 1, &req, &status);

    CHECK_ERR
    err = status; CHECK_ERR

    /* check the contents of put buffer are altered */
    CHECK_PUT_BUF(buf_int)

    err = ncmpi_buffer_detach(ncid); CHECK_ERR

    /* read back using a blocking get_varm flexible API ---------------------*/
    /* initiate get buffer contents */
    INIT_GET_BUF(buf_int)

    if (!coll_io) {
        /* calling a blocking flexible API */
        err = ncmpi_end_indep_data(ncid);
        CHECK_ERR
    }
    err = ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
                             buf_int[0], 1, subarray);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* check the contents of get buffer */
    CHECK_GET_BUF(buf_int)

    MPI_Type_free(&subarray);

    /* test case for no type conversion =====================================*/
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C, MPI_DOUBLE,
                             &subarray);
    MPI_Type_commit(&subarray);

    /* calling a nonblocking bput_varm flexible API -------------------------*/
    /* initiate put buffer contents */
    INIT_PUT_BUF(buf_dbl)

    err = ncmpi_buffer_attach(ncid, bufsize); CHECK_ERR

    err = ncmpi_bput_varm(ncid, varid, start, count, stride, imap, buf_dbl[0],
                          1, subarray, &req);
    CHECK_ERR

    /* check the contents of put buffer are altered */
    CHECK_PUT_BUF(buf_dbl)

    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, &status);
    else
        err = ncmpi_wait(ncid, 1, &req, &status);
    CHECK_ERR
    err = status; CHECK_ERR

    /* check the contents of put buffer are altered */
    CHECK_PUT_BUF(buf_dbl)

    err = ncmpi_buffer_detach(ncid); CHECK_ERR

    /* read back using a blocking get_varm flexible API ---------------------*/
    /* initiate get buffer contents */
    INIT_GET_BUF(buf_dbl)

    /* calling a blocking flexible API */
    if (coll_io)
        err = ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
                                 buf_dbl[0], 1, subarray);
    else
        err = ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                                 buf_dbl[0], 1, subarray);
    CHECK_ERR

    /* check the contents of get buffer */
    CHECK_GET_BUF(buf_dbl)

err_out:
    MPI_Type_free(&subarray);

    err = ncmpi_close(ncid); CHECK_ERR

    free(buf_int[0]);
    free(buf_int);
    free(buf_dbl[0]);
    free(buf_dbl);

    return nerrs;
}

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
    int i, len, rank, nprocs, nerrs=0, err=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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
        sprintf(cmd_str, "*** TESTING C   %s for flexible bput_varm ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    nerrs = tst_run(filename, 0, len);
    if (nerrs != 0) goto err_out;

    nerrs = tst_run(filename, 1, len);
    if (nerrs != 0) goto err_out;

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0) {
            printf("heap memory allocated by PnetCDF internally has "OFFFMT" bytes yet to be freed\n",
                   sum_size);
        }
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

