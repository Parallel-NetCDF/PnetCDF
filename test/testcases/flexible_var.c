/*********************************************************************
 *
 *  Copyright (C) 2022, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests PnetCDF flexible var APIs, i.e. ncmpi_put_var_all(),
 * ncmpi_get_var_all(), and their independent-mode, nonblocking versions, to
 * write a 2D array double variable of size NY x NX in parallel.
 *
 * The global 2D array is not partitioned.
 *
 * The local buffer on each process has ghost cells surrounded along both
 * dimensions.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -O2 -o flexible_var flexible_var.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./flexible_var /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            Y = 6 ;
 *            X = 4 ;
 *    variables:
 *            double var(Y, X) ;
 *    data:
 *
 *    var =
 *      0, 1, 2, 3,
 *      4, 5, 6, 7,
 *      8, 9, 10, 11,
 *      12, 13, 14, 15,
 *      16, 17, 18, 19,
 *      20, 21, 22, 23 ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

/* Smaller NY and NX are for demonstration, producing outputs as shown above.
#define NY 6
#define NX 4
*/
#define NY 32
#define NX 128
#define GHOST 2

#define PRINT_PUT_BUF(Y, X, buf) { \
    for (i=0; i<Y; i++) \
        for (j=0; j<X; j++) \
            printf("buf[%d][%d] = %d\n", i, j, buf[i][j]); \
}

#define INIT_PUT_BUF(Y, X, buf) { \
    for (i=0; i<Y; i++) \
        for (j=0; j<X; j++) \
            buf[i][j] = i*X + j; \
}

#define INIT_PUT_BUF_GHOST(Y, X, G, buf) \
    for (i=0; i<Y+2*G; i++) { \
        for (j=0; j<X+2*G; j++) { \
            if (i < G || G+Y <= i || j < G || G+X <= j) \
                buf[i][j] = -1; \
            else \
                buf[i][j] = (i-G)*X+(j-G); \
        } \
    }

#define CHECK_PUT_BUF(Y, X, buf) \
    for (i=0; i<Y; i++) { \
        for (j=0; j<X; j++) { \
            if (buf[i][j] != i*X+j) { \
                printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%d\n", \
                       __LINE__,__FILE__,i,j,buf[i][j]); \
                nerrs++; \
            } \
        } \
    }

#define CHECK_PUT_BUF_GHOST(Y, X, G, buf) \
    for (i=0; i<Y; i++) { \
        for (j=0; j<X; j++) { \
            if (i < G || G+Y <= i || j < G || G+X <= j) { \
                if (buf[i][j] != -1) { \
                    printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%d\n", \
                           __LINE__,__FILE__,i,j,buf[i][j]); \
                    nerrs++; \
                } \
            } \
            else if (buf[i][j] != (i-G)*X+(j-G)) { \
                printf("Error at line %d in %s: put buffer altered buffer[%d][%d]=%d\n", \
                       __LINE__,__FILE__,i,j,buf[i][j]); \
                nerrs++; \
            } \
        } \
    }

#define INIT_GET_BUF(Y, X, buf) {\
    for (i=0; i<Y; i++) \
        for (j=0; j<X; j++) \
            buf[i][j] = -2; \
}

#define CHECK_GET_BUF(Y, X, buf) \
    for (i=0; i<Y; i++) { \
        for (j=0; j<X; j++) { \
            if (buf[i][j] != i*X+j) { \
                printf("Error at line %d in %s: Unexpected get buffer[%d][%d]=%d\n", \
                       __LINE__,__FILE__,i,j,buf[i][j]); \
                nerrs++; \
            } \
        } \
    }
#define CHECK_GET_BUF_GHOST(Y, X, G, buf) \
    for (i=0; i<Y; i++) { \
        for (j=0; j<X; j++) { \
            if (i < G || G+Y <= i || \
                j < G || G+X <= j) { \
                if (buf[i][j] != -2) { \
                    printf("Error at line %d in %s: Unexpected get buffer[%d][%d]=%d\n", \
                           __LINE__,__FILE__,i,j,buf[i][j]); \
                    nerrs++; \
                } \
            } \
            else if (buf[i][j] != (i-G)*X+(j-G)) { \
                printf("Error at line %d in %s: Unexpected get buffer[%d][%d]=%d\n", \
                       __LINE__,__FILE__,i,j,buf[i][j]); \
                nerrs++; \
            } \
        } \
    }

static
int tst_io(const char *out_path,
           MPI_Info    info)
{
    int i, j, rank, nprocs, err, nerrs=0, req, status;
    int ncid, varid, dimid[2];
    int array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    int buf_ghost[NY+2*GHOST][NX+2*GHOST];
    int buf[NY][NX];
    MPI_Offset bufcount, start[2], count[2];
    MPI_Datatype subarray;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERROUT

    /* define 2 dimensions */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); CHECK_ERROUT
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]); CHECK_ERROUT

    /* define a variable of size NY * (NX * nprocs) */
    err = ncmpi_def_var(ncid, "var", NC_DOUBLE, 2, dimid, &varid); CHECK_ERROUT
    err = ncmpi_enddef(ncid); CHECK_ERROUT

    /* var is partitioned along X dimension in a matrix transported way */
    array_of_sizes[0]    = NY + 2*GHOST;
    array_of_sizes[1]    = NX + 2*GHOST;
    array_of_subsizes[0] = NY;
    array_of_subsizes[1] = NX;
    array_of_starts[0]   = GHOST;
    array_of_starts[1]   = GHOST;
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C, MPI_INT, &subarray);
    MPI_Type_commit(&subarray);

    start[0] = start[1] = 0;
    count[0] = count[1] = 0;

    /*----------------------------------------------------------------------*/
    /*---- test using bufcount == 1 with ghost cells -----------------------*/
    /*----------------------------------------------------------------------*/
    bufcount = 1;

    /* initiate put buffer contents (ghost cells are all -1) */
    INIT_PUT_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* calling a blocking put_var flexible API -----------------------------*/
    if (rank == 0) /* only rank 0 writes to the variable */
        err = ncmpi_put_var_all(ncid, varid, buf_ghost, bufcount, subarray);
    else /* other ranks write 0-sized data */
        err = ncmpi_put_vara_all(ncid, varid, start, count, buf_ghost, 0, MPI_INT);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    CHECK_ERROUT

    /* check the contents of put buffer. They should not be altered. */
    CHECK_PUT_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* read back and check the contents written in the file  ----------------*/
    INIT_GET_BUF(NY, NX, buf_ghost)
    err = ncmpi_get_var_all(ncid, varid, buf_ghost, bufcount, subarray);
    CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* read back using a non-blocking flexible API --------------------------*/
    INIT_GET_BUF(NY, NX, buf_ghost)
    err = ncmpi_iget_var(ncid, varid, buf_ghost, bufcount, subarray, &req);
    CHECK_ERROUT
    err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* initiate put buffer contents (ghost cells are all -1) */
    INIT_PUT_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* calling a nonblocking put_var flexible API --------------------------*/
    if (rank == 0) { /* only rank 0 writes to the variable */
        err = ncmpi_iput_var(ncid, varid, buf_ghost, bufcount, subarray, &req);
        CHECK_ERROUT
        /* check the contents of put buffer. They should not be altered. */
        CHECK_PUT_BUF_GHOST(NY, NX, GHOST, buf_ghost)
    }
    err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* read back and check the contents written in the file  ----------------*/
    INIT_GET_BUF(NY, NX, buf_ghost)
    err = ncmpi_get_var_all(ncid, varid, buf_ghost, bufcount, subarray);
    CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* read back using a non-blocking flexible API --------------------------*/
    INIT_GET_BUF(NY, NX, buf_ghost)
    err = ncmpi_iget_var(ncid, varid, buf_ghost, bufcount, subarray, &req);
    CHECK_ERROUT
    err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /*----------------------------------------------------------------------*/
    /*---- test using bufcount == NC_COUNT_IGNORE with no ghost cells ------*/
    /*----------------------------------------------------------------------*/
    bufcount = NC_COUNT_IGNORE;

    /* initiate put buffer contents, no ghost cells */
    INIT_PUT_BUF(NY, NX, buf)

    /* calling a blocking put_var flexible API -----------------------------*/
    if (rank == 0) /* only rank 0 writes to the variable */
        err = ncmpi_put_var_all(ncid, varid, buf, bufcount, MPI_INT);
    else /* other ranks write 0-sized data */
        err = ncmpi_put_vara_all(ncid, varid, start, count, buf, 0, MPI_INT);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    CHECK_ERROUT

    /* check the contents of put buffer. They should not be altered. */
    CHECK_PUT_BUF(NY, NX, buf)

    /* read back and check the contents written in the file  ----------------*/
    INIT_GET_BUF(NY, NX, buf)
    err = ncmpi_get_var_all(ncid, varid, buf, bufcount, MPI_INT);
    CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF(NY, NX, buf)

    /* read back using a non-blocking flexible API --------------------------*/
    INIT_GET_BUF(NY, NX, buf)
    err = ncmpi_iget_var(ncid, varid, buf, bufcount, MPI_INT, &req);
    CHECK_ERROUT
    err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF(NY, NX, buf)

    /* initiate put buffer contents, no ghost cells */
    INIT_PUT_BUF(NY, NX, buf)

    /* calling a nonblocking put_var flexible API --------------------------*/
    if (rank == 0) { /* only rank 0 writes to the variable */
        err = ncmpi_iput_var(ncid, varid, buf, bufcount, MPI_INT, &req);
        CHECK_ERROUT
        /* check the contents of put buffer. They should not be altered. */
        CHECK_PUT_BUF(NY, NX, buf)
    }
    err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* read back and check the contents written in the file  ----------------*/
    INIT_GET_BUF(NY, NX, buf)
    err = ncmpi_get_var_all(ncid, varid, buf, bufcount, MPI_INT);
    CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF(NY, NX, buf)

    /* read back using a non-blocking flexible API --------------------------*/
    INIT_GET_BUF(NY, NX, buf)
    err = ncmpi_iget_var(ncid, varid, buf, bufcount, MPI_INT, &req);
    CHECK_ERROUT
    err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF(NY, NX, buf)

    /*----------------------------------------------------------------------*/
    /*---- test independent I/O mode ---------------------------------------*/
    /*----------------------------------------------------------------------*/
    err = ncmpi_begin_indep_data(ncid);
    CHECK_ERROUT

    /*----------------------------------------------------------------------*/
    /*---- test using bufcount == 1 with ghost cells -----------------------*/
    /*----------------------------------------------------------------------*/
    bufcount = 1;

    /* initiate put buffer contents (ghost cells are all -1) */
    INIT_PUT_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* calling a blocking put_var flexible API -----------------------------*/
    if (rank == 0) { /* only rank 0 writes to the variable */
        err = ncmpi_put_var(ncid, varid, buf_ghost, bufcount, subarray);
        CHECK_ERROUT
        /* check the contents of put buffer. They should not be altered. */
        CHECK_PUT_BUF_GHOST(NY, NX, GHOST, buf_ghost)
    }

    /* file sync is required for non-zero ranks to see the data in file */
    err = ncmpi_sync(ncid);
    CHECK_ERROUT

    /* read back and check the contents written in the file  ----------------*/
    INIT_GET_BUF(NY, NX, buf_ghost)
    err = ncmpi_get_var(ncid, varid, buf_ghost, bufcount, subarray);
    CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* read back using a non-blocking flexible API --------------------------*/
    INIT_GET_BUF(NY, NX, buf_ghost)
    err = ncmpi_iget_var(ncid, varid, buf_ghost, bufcount, subarray, &req);
    CHECK_ERROUT
    err = ncmpi_wait(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* initiate put buffer contents (ghost cells are all -1) */
    INIT_PUT_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* calling a nonblocking put_var flexible API --------------------------*/
    if (rank == 0) { /* only rank 0 writes to the variable */
        err = ncmpi_iput_var(ncid, varid, buf_ghost, bufcount, subarray, &req);
        CHECK_ERROUT
        /* check the contents of put buffer. They should not be altered. */
        CHECK_PUT_BUF_GHOST(NY, NX, GHOST, buf_ghost)
    }
    err = ncmpi_wait(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* read back and check the contents written in the file  ----------------*/
    INIT_GET_BUF(NY, NX, buf_ghost)
    err = ncmpi_get_var(ncid, varid, buf_ghost, bufcount, subarray);
    CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /* read back using a non-blocking flexible API --------------------------*/
    INIT_GET_BUF(NY, NX, buf_ghost)
    err = ncmpi_iget_var(ncid, varid, buf_ghost, bufcount, subarray, &req);
    CHECK_ERROUT
    err = ncmpi_wait(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF_GHOST(NY, NX, GHOST, buf_ghost)

    /*----------------------------------------------------------------------*/
    /*---- test using bufcount == NC_COUNT_IGNORE with no ghost cells ------*/
    /*----------------------------------------------------------------------*/
    bufcount = NC_COUNT_IGNORE;

    /* initiate put buffer contents, no ghost cells */
    INIT_PUT_BUF(NY, NX, buf)

    /* calling a blocking put_var flexible API -----------------------------*/
    if (rank == 0) { /* only rank 0 writes to the variable */
        err = ncmpi_put_var(ncid, varid, buf, bufcount, MPI_INT);
        CHECK_ERROUT
        /* check the contents of put buffer. They should not be altered. */
        CHECK_PUT_BUF(NY, NX, buf)
    }

    /* file sync is required for non-zero ranks to see the data in file */
    err = ncmpi_sync(ncid);
    CHECK_ERROUT

    /* read back and check the contents written in the file  ----------------*/
    INIT_GET_BUF(NY, NX, buf)
    err = ncmpi_get_var(ncid, varid, buf, bufcount, MPI_INT);
    CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF(NY, NX, buf)

    /* read back using a non-blocking flexible API --------------------------*/
    INIT_GET_BUF(NY, NX, buf)
    err = ncmpi_iget_var(ncid, varid, buf, bufcount, MPI_INT, &req);
    CHECK_ERROUT
    err = ncmpi_wait(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF(NY, NX, buf)

    /* initiate put buffer contents, no ghost cells */
    INIT_PUT_BUF(NY, NX, buf)

    /* calling a nonblocking put_var flexible API --------------------------*/
    if (rank == 0) { /* only rank 0 writes to the variable */
        err = ncmpi_iput_var(ncid, varid, buf, bufcount, MPI_INT, &req);
        CHECK_ERROUT
        /* check the contents of put buffer. They should not be altered. */
        CHECK_PUT_BUF(NY, NX, buf)
    }
    err = ncmpi_wait(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* read back and check the contents written in the file  ----------------*/
    INIT_GET_BUF(NY, NX, buf)
    err = ncmpi_get_var(ncid, varid, buf, bufcount, MPI_INT);
    CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF(NY, NX, buf)

    /* read back using a non-blocking flexible API --------------------------*/
    INIT_GET_BUF(NY, NX, buf)
    err = ncmpi_iget_var(ncid, varid, buf, bufcount, MPI_INT, &req);
    CHECK_ERROUT
    err = ncmpi_wait(ncid, 1, &req, &status); CHECK_ERROUT
    err = status; CHECK_ERROUT

    /* check the contents of get buffer */
    CHECK_GET_BUF(NY, NX, buf)

    MPI_Type_free(&subarray);

    err = ncmpi_close(ncid); CHECK_ERROUT

err_out:
    return (nerrs > 0);
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int err, nerrs=0;
    MPI_Info info_dup;

    MPI_Info_dup(info, &info_dup);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    nerrs = tst_io(out_path, info_dup);
    if (nerrs > 0) return nerrs;

    /* disable PnetCDF internal buffering */
    MPI_Info_set(info_dup, "nc_ibuf_size", "0");

    nerrs = tst_io(out_path, info_dup);
    if (nerrs > 0) return nerrs;

    MPI_Info_free(&info_dup);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "flexible var APIs", opt, test_io);

    MPI_Finalize();

    return err;
}
