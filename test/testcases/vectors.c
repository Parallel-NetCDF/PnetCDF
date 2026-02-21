/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 *
 *  Test flexible API ncmpi_put_vara()
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define VECCOUNT 4
#define BLOCKLEN 3
#define STRIDE   5

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int ncid, dimid, varid, rank, nprocs;
    MPI_Datatype vtype, rtype, usertype;
    MPI_Aint lb, extent;
    int *userbuf, *cmpbuf, i, err, errs=0, nerrs=0;
    size_t userbufsz;
    int count = 25;
    double pi = 3.14159;
    MPI_Offset start, acount;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#ifdef DEBUG
    if (nprocs > 2 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n",argv[0]);
#endif

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, "50k", 1024*50, &dimid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "vector", NC_DOUBLE, 1, &dimid, &varid);
    CHECK_ERR

    err = ncmpi_def_var_fill(ncid, varid, 0, NULL);
    CHECK_ERR

    err = ncmpi_enddef(ncid);
    CHECK_ERR

    MPI_Type_vector(VECCOUNT, BLOCKLEN, STRIDE, MPI_INT, &vtype);
    MPI_Type_create_resized(vtype, 0, STRIDE*VECCOUNT*sizeof(int), &rtype);
    MPI_Type_contiguous(count, rtype, &usertype);
    MPI_Type_commit(&usertype);

    MPI_Type_free(&vtype);
    MPI_Type_free(&rtype);

    MPI_Type_get_extent(usertype, &lb, &extent);
    userbufsz = extent;
    userbuf = (int*) malloc(userbufsz);
    cmpbuf = (int*) calloc(userbufsz, 1);
    for (i=0; i< userbufsz/sizeof(int); i++)
        userbuf[i] = pi*i;

    start = 10; acount = count*12;
    err = ncmpi_begin_indep_data(ncid);
    CHECK_ERR
    if (rank == 0) {
        err = ncmpi_put_vara(ncid, varid, &start, &acount, userbuf, 1, usertype);
        CHECK_ERR
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR

    err = ncmpi_close(ncid);
    CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR
    err = ncmpi_begin_indep_data(ncid);
    CHECK_ERR
    err = ncmpi_inq_varid(ncid, "vector", &varid);
    CHECK_ERR
    err = ncmpi_get_vara(ncid, varid, &start, &acount, cmpbuf, 1, usertype);
    CHECK_ERR

    MPI_Type_free(&usertype);

    err = ncmpi_close(ncid);
    CHECK_ERR

    for (i=0; errs < 10 &&  i < acount; i++) {
        /* vector of 4,3,5, so skip 4th and 5th items of every block */
        if (i%STRIDE >= BLOCKLEN) continue;
        if (userbuf[i] != cmpbuf[i]) {
            errs++;
            fprintf(stderr, "%d: expected 0x%x got 0x%x\n",
                    i, userbuf[i], cmpbuf[i]);
        }
    }
    nerrs += errs;
    free(userbuf);
    free(cmpbuf);

    return (nerrs > 0);
}

int main(int argc, char **argv) {

    int err;

    /* flexible APIs are not supported in NetCDF4 */
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
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "put_vara/get_vara", opt, test_io);

    MPI_Finalize();

    return err;
}
