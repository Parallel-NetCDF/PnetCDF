/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests redefine mode.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o redef1 redef1.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 redef1 testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, j, k, nprocs, rank, ncid, err, nerrs=0;
    int dim0id, dim1id, dim5id, dim9id, dim2id, dimsid[2], dims2id[2];
    int varid, var3id, var4id, var2id;
    int *data;
    double *dbl_data;
    MPI_Offset len0=10, len1=3, len5=5, len9=9, len2=10;
    MPI_Offset start[2], count[2];
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef DEBUG
    if (nprocs > 1 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n",
               basename(__FILE__));
#endif

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* Test NetCDF 4 first as ncvalidator checks only classic formats */
    err = ncmpi_create(comm, out_path, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "dim0", len0, &dim0id); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", len1, &dim1id); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim5", len5, &dim5id); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim9", len9, &dim9id); CHECK_ERR

    dimsid[0] = dim0id;
    dimsid[1] = dim1id;
    err = ncmpi_def_var(ncid, "xyz", NC_INT, 2, dimsid, &varid); CHECK_ERR

    dimsid[0] = dim0id;
    dimsid[1] = dim5id;
    err = ncmpi_def_var(ncid, "connect", NC_INT, 2, dimsid, &var3id); CHECK_ERR

    dimsid[0] = dim0id;
    dimsid[1] = dim9id;
    err = ncmpi_def_var(ncid, "connect_exterior", NC_INT, 2, dimsid, &var4id); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* put data */
    start[0] = 0;
    start[1] = 0;
    count[0] = len0;
    count[1] = len1;

    data = (int*) malloc(sizeof(int) * len0*len1);
    k=0;
    for (i=0; i<len0; i++)
        for (j=0; j<len1; j++)
            data[i*len1+j] = k++;
    if (rank > 0) count[0] = count[1] = 0;
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid, start, count, &data[0]);
    else
        err = ncmpi_put_vara_int(ncid, varid, start, count, &data[0]);
    CHECK_ERR
    free(data);

    count[0] = len0;
    count[1] = len5;
    data = (int*) malloc(sizeof(int) * len0*len5);
    k=0;
    for (i=0; i<len0; i++)
        for (j=0; j<len5; j++)
            data[i*len5+j] = k++;
    if (rank > 0) count[0] = count[1] = 0;
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, var3id, start, count, &data[0]);
    else
        err = ncmpi_put_vara_int(ncid, var3id, start, count, &data[0]);
    CHECK_ERR
    free(data);

    count[0] = len0;
    count[1] = len9;
    data = (int*) malloc(sizeof(int) * len0*len9);
    k=0;
    for (i=0; i<len0; i++)
        for (j=0; j<len9; j++)
            data[i*len9+j] = k++;
    if (rank > 0) count[0] = count[1] = 0;
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, var4id, start, count, &data[0]);
    else
        err = ncmpi_put_vara_int(ncid, var4id, start, count, &data[0]);
    CHECK_ERR
    free(data);

    /* file sync before file close and re-open it */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(comm, out_path, NC_WRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    err = ncmpi_redef(ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "dim2", len2, &dim2id); CHECK_ERR

    dims2id[0] = dim0id;
    dims2id[1] = dim2id;
    err = ncmpi_def_var(ncid, "xyz_r", NC_DOUBLE, 2, dims2id, &var2id);
    CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    start[0] = 0;
    start[1] = 0;
    count[0] = len0;
    count[1] = len2;
    k=0;
    dbl_data = (double*) malloc(sizeof(double) * len0*len2);
    for (i=0; i<len0; i++)
        for (j=0; j<len2; j++) {
            dbl_data[i*len2+j] = (k*k);
            k++;
        }
    if (rank > 0) count[0] = count[1] = 0;
    if (coll_io)
        err = ncmpi_put_vara_double_all(ncid, var2id, start, count, &dbl_data[0]);
    else
        err = ncmpi_put_vara_double(ncid, var2id, start, count, &dbl_data[0]);
    CHECK_ERR
    free(dbl_data);

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 1; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "entering re-define mode", opt, test_io);

    MPI_Finalize();

    return err;
}
