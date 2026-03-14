/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This example tests if PnetCDF duplicates the MPI derived data type supplied
 * by the user, when calling the flexible APIs. It tests a PnetCDF bug
 * prior to version 1.6.1.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 4
#define NX 100
#define NVARS 4
#define NGHOSTS 2

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, j, k, err, ncid, nerrs=0, rank, nprocs;
    int varid[NVARS], dimids[2], req[NVARS], st[NVARS], *buf[NVARS];
    int gsize[2], subsize[2], a_start[2], ghost;
    MPI_Offset start[2], count[2];
    MPI_Datatype buftype[NVARS];

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "Y", NY*nprocs, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX,        &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 2, dimids, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 2, dimids, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var2", NC_INT, 2, dimids, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var3", NC_INT, 2, dimids, &varid[3]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* initialize the contents of the array */
    start[0] = NY*rank; start[1] = 0;
    count[0] = NY;      count[1] = NX;

    for (i=0; i<NVARS; i++) {
        buf[i] = (int*) malloc(sizeof(int) * count[0] * count[1]);
        for (j=0; j<count[0]*count[1]; j++) buf[i][j] = j + rank*10;
    }

    if (coll_io) {
        err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf[0]);
        CHECK_ERR
        err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf[1]);
        CHECK_ERR
        err = ncmpi_put_vara_int_all(ncid, varid[2], start, count, buf[2]);
        CHECK_ERR
        err = ncmpi_put_vara_int_all(ncid, varid[3], start, count, buf[3]);
        CHECK_ERR
    }
    else {
        err = ncmpi_put_vara_int(ncid, varid[0], start, count, buf[0]);
        CHECK_ERR
        err = ncmpi_put_vara_int(ncid, varid[1], start, count, buf[1]);
        CHECK_ERR
        err = ncmpi_put_vara_int(ncid, varid[2], start, count, buf[2]);
        CHECK_ERR
        err = ncmpi_put_vara_int(ncid, varid[3], start, count, buf[3]);
        CHECK_ERR
    }

    /* check if user write buffer contents altered */
    for (i=0; i<NVARS; i++) {
        for (j=0; j<count[0]*count[1]; j++) {
            int exp = j + rank*10;
            if (buf[i][j] != exp) {
                printf("Error at line %d in %s: user put buffer[%d][%d] altered from %d to %d\n",
                       __LINE__,__FILE__,i,j, exp, buf[i][j]);
                nerrs++;
            }
        }
    }

    ghost    = NGHOSTS;
    gsize[1] = NX + 2 * ghost;
    gsize[0] = NY + 2 * ghost;

    /* reset buffer contents before reads */
    for (i=0; i<NVARS; i++) {
        free(buf[i]);
        buf[i] = (int*) malloc(sizeof(int) * gsize[0] * gsize[1]);
        for (j=0; j<gsize[0]*gsize[1]; j++)
            buf[i][j] = -1;
    }

    /* define an MPI datatype using MPI_Type_create_subarray() */
    subsize[1] = NX;
    subsize[0] = NY;
    a_start[1] = ghost;
    a_start[0] = ghost;

    for (i=0; i<NVARS; i++) {
        req[i] = NC_REQ_NULL;
        st[i]  = NC_NOERR;
        MPI_Type_create_subarray(2, gsize, subsize, a_start, MPI_ORDER_C,
                                 MPI_INT, &buftype[i]);
        MPI_Type_commit(&buftype[i]);

        err = ncmpi_iget_vara(ncid, varid[i], start, count, buf[i], 1,
                              buftype[i], &req[i]);
        CHECK_ERR
        MPI_Type_free(&buftype[i]);
    }

    if (coll_io)
        err = ncmpi_wait_all(ncid, NVARS, req, st);
    else
        err = ncmpi_wait(ncid, NVARS, req, st);
    CHECK_ERR

    for (i=0; i<NVARS; i++) {
        if (st[i] != NC_NOERR) {
            err = st[i];
            CHECK_ERR
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* check contents of read buffers */
    for (i=0; i<NVARS; i++) {
        for (j=0; j<gsize[0]; j++) {
            for (k=0; k<gsize[1]; k++) {
                int exp;
                if (j < ghost || j >= gsize[0]-ghost ||
                    k < ghost || k >= gsize[1]-ghost)
                    exp = -1;
                else
                    exp = (int)((j-ghost)*count[1]+(k-ghost) + rank*10);
                if (buf[i][j*gsize[1]+k] != exp) {
                    printf("Error at %d: var %d expect buf[%d][%d] = %d but got %d\n",
                           __LINE__, i, j, k, exp, buf[i][j*gsize[1]+k]);
                    nerrs++;
                    goto err_out;
                }
            }
        }
    }

    for (i=0; i<NVARS; i++)
        free(buf[i]);

err_out:
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
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "free buftype in flexible API", opt, test_io);

    MPI_Finalize();

    return err;
}
