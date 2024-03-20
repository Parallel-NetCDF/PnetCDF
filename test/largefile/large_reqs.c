/*********************************************************************
 *
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program is to test writing and reading > 4GB in a single call to
 * MPI_File_write call. The user buffer is of size > 4GB per MPI rank.
 *
 * Two tests are includes:
 * 1. writing/reading one large variable (> 4GB)
 * 2. writing/reading multiple smaller variables of total size > 4GB.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

#define NY 1280
#define NX 1048576

static int verbose;

static
int tst_one_var(char *filename, MPI_Comm comm)
{
    size_t i, buf_len;
    int rank, nprocs, err, nerrs=0, ncid, cmode, varid, dimid[3], psize[2];
    int *buf;
    MPI_Offset start[3], count[3];

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    /* Creates a division of processors in a cartesian grid */
    psize[0] = psize[1] = 0;
    MPI_Dims_create(nprocs, 2, psize);

    /* Test classic CDF-5 format */
    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* define dimensions Z, Y, and X */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NY*psize[0], &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*psize[1], &dimid[2]); CHECK_ERR

    /* define a big 2D fixed-size variable of integer type */
    err = ncmpi_def_var(ncid, "var", NC_INT, 3, dimid, &varid); CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* now we are in data mode */
    start[0] = 0;
    start[1] = NY * (rank / psize[1]);
    start[2] = NX * (rank % psize[1]);
    count[0] = 1;

    if (verbose) {
        fflush(stdout);
        MPI_Barrier(comm);
        if (rank == 0) {
            float len = (float)NY*psize[0]*NX*psize[1]*sizeof(int);
            printf("\nglobal array is of size %d x %d = %.1f GiB\n",
                   NY*psize[0], NX*psize[1], len/1073741824);
        }
        printf("rank %d start=%lld %lld\n", rank, start[1],start[2]);
    }

    /* user buffer is contiguous */
    buf_len = (size_t)NY * NX;
    buf = (int*) malloc(buf_len * sizeof(int));
    for (i=0; i<buf_len; i++) buf[i] = (i + rank) % 128;

    /* write the entire variable */
    count[1] = NY;
    count[2] = NX;

    if (verbose)
        printf("write entire var - rank %d: write amount=%.1f GiB\n",
               rank, (float)count[1]*count[2]*sizeof(int)/1073741824);

    /* write */
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* read */
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* write the variable partially */
    count[1] = NY - 128;
    count[2] = NX - 128;

    if (verbose)
        printf("write partial var - rank %d: write amount=%.1f GiB\n",
               rank, (float)count[1]*count[2]*sizeof(int)/1073741824);

    /* write using contiguous user buffer */
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* read using contiguous user buffer */
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* Make user buffer non-contiguous */
    int gsize[2], lsize[2], lstart[2];
    MPI_Datatype buftype;

    gsize[0] = NY;
    gsize[1] = NX;
    lsize[0] = (int)count[1];
    lsize[1] = (int)count[2];
    lstart[0] = 0;
    lstart[1] = 0;
    MPI_Type_create_subarray(2, gsize, lsize, lstart, MPI_ORDER_C,
                             MPI_INT, &buftype);
    MPI_Type_commit(&buftype);

    if (verbose)
        printf("write noncontig buf - rank %d: write amount=%.1f GiB\n",
               rank, (float)count[1]*count[2]*sizeof(int)/1073741824);

    /* write */
    err = ncmpi_put_vara_all(ncid, varid, start, count, buf, 1, buftype);
    CHECK_ERR

    /* read */
    err = ncmpi_get_vara_all(ncid, varid, start, count, buf, 1, buftype);
    CHECK_ERR

    MPI_Type_free(&buftype);

    err = ncmpi_close(ncid); CHECK_ERR

    free(buf);

    return nerrs;
}

#define NVARS 1100
#define LEN 1024

static
int tst_vars(char *filename, MPI_Comm comm)
{
    size_t i, buf_len;
    int rank, nprocs, err, nerrs=0, *buf, *buf_ptr;
    int ncid, cmode, *varid, dimid[3], gap, psize[2]={0,0};
    MPI_Offset start[3], count[3];

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    /* Creates a division of processors in a cartesian grid */
    psize[0] = psize[1] = 0;
    MPI_Dims_create(nprocs, 2, psize);

    /* Test classic CDF-5 format */
    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* define dimensions Z, Y, and X */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", LEN*psize[0], &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", LEN*psize[1], &dimid[2]); CHECK_ERR

    varid = (int*) malloc(NVARS * sizeof(int));

    for (i=0; i<NVARS; i++) {
        /* define a 2D variables of integer type */
        char name[32];
        sprintf(name, "var.%zd", i);
        err = ncmpi_def_var(ncid, name, NC_INT, 3, dimid, &varid[i]);
        CHECK_ERR
    }

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (verbose && rank == 0) {
        float len = (float)LEN*psize[0]*LEN*psize[1]*sizeof(int);
        printf("Number of global arrays = %d\n", NVARS);
        printf("Each global array is of size %d x %d = %.1f MiB\n",
               LEN*psize[0], LEN*psize[1], len/1048576);
    }

    /* now we are in data mode */

    /* make user buffer noncontiguous */
    gap = 2;
    buf_len = (LEN * LEN + gap) * NVARS;
    buf = (int*) malloc(buf_len * sizeof(int));
    for (i=0; i<buf_len; i++) buf[i] = (i + rank) % 128;

    /* create a subarray datatype for user buffer */
    int gsize[2], lsize[2], lstart[2];
    MPI_Datatype buftype;

    gsize[0] = LEN;
    gsize[1] = LEN;
    lsize[0] = LEN - gap;
    lsize[1] = LEN - gap;
    lstart[0] = 0;
    lstart[1] = 0;
    MPI_Type_create_subarray(2, gsize, lsize, lstart, MPI_ORDER_C,
                             MPI_INT, &buftype);
    MPI_Type_commit(&buftype);

    /* set subarray offset and length */
    start[0] = 0;
    start[1] = LEN * (rank / psize[1]);
    start[2] = LEN * (rank % psize[1]);
    count[0] = 1;
    count[1] = lsize[0];
    count[2] = lsize[1];

    if (verbose)
        printf("rank %d start=%lld %lld count=%lld %lld\n",
               rank, start[1],start[2], count[1],count[2]);

    if (verbose)
        printf("%d: nonblocking write total amount = %.1f GiB\n",
               rank, (float)count[1]*count[2]*NVARS*sizeof(int)/1073741824);

    /* write */
    buf_ptr = buf;
    for (i=0; i<NVARS; i++) {
        /* write using non-contiguous user buffer */
        err = ncmpi_iput_vara(ncid, varid[i], start, count, buf_ptr, 1,
                              buftype, NULL);
        CHECK_ERR
        buf_ptr += LEN * LEN + gap;
    }

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    /* read */
    buf_ptr = buf;
    for (i=0; i<NVARS; i++) {
        /* read using non-contiguous user buffer */
        err = ncmpi_iget_vara(ncid, varid[i], start, count, buf_ptr, 1,
                              buftype, NULL);
        CHECK_ERR
        buf_ptr += LEN * LEN;
    }

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    MPI_Type_free(&buftype);

    free(varid);
    free(buf);

    return nerrs;
}

int main(int argc, char** argv)
{
    char filename[256];
    int rank, nprocs, err, nerrs=0, color;
    MPI_Comm comm;

    verbose = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for large requests ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    color = 1;

    if (nprocs > 2) {
        /* run on 2 ranks only, as this test allocates memory > 4GB per rank */
        /* split MPI_COMM_WORLD based on 'color' and use the same rank order */
        color = (rank < 2) ? 1 : 0;
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    }
    else
        comm = MPI_COMM_WORLD;

    if (color) {
        /* test one big variable */
        nerrs += tst_one_var(filename, comm);

        /* test a large number of smaller variables */
        nerrs += tst_vars(filename, comm);
    }

    if (comm != MPI_COMM_WORLD) MPI_Comm_free(&comm);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

