/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program is to test
 *
 * large number of dimensions per variable
 * large number of dimensions per file
 * large number of attributes per file
 * large number of variables per file
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

#define LARGE_NUM 102400

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    char str[32];
    int i, rank, nprocs, color, err, nerrs=0;
    int ncid, *varid, *dimids, intBuf[1];
    MPI_Comm comm=MPI_COMM_NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    color = 1;

    if (nprocs > 2) {
        /* run on 2 ranks only, as this test allocates memory > 4GB per rank */
        /* split MPI_COMM_WORLD based on 'color' and use the same rank order */
        color = (rank < 2) ? 1 : 0;
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    }
    else
        comm = MPI_COMM_WORLD;

    if (!color) goto err_out;

    dimids = (int*) malloc(sizeof(int) * LARGE_NUM);
    varid = (int*) malloc(sizeof(int) * LARGE_NUM);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    for (i=0; i<LARGE_NUM; i++) {
        sprintf(str, "dim%d", i);
        err = ncmpi_def_dim(ncid, str, 1, &dimids[i]);
        CHECK_ERR
    }

    err = ncmpi_def_var(ncid, "var", NC_INT, LARGE_NUM, dimids, &varid[0]);
    CHECK_ERR

    for (i=0; i<LARGE_NUM; i++) {
        sprintf(str, "attr%d", i);
        err = ncmpi_put_att(ncid, varid[0], str, NC_INT, 1, &i);
        CHECK_ERR
    }

    for (i=1; i<LARGE_NUM; i++) {
        signed char attrBuf[3]={1,2,3};
        sprintf(str, "var%d", i);
        err = ncmpi_def_var(ncid, str, NC_INT, 1, dimids, &varid[i]);
        CHECK_ERR
        err = ncmpi_put_att_text(ncid, varid[i], "attr text", 9, "some text");
        CHECK_ERR
        err = ncmpi_put_att_schar(ncid, varid[i], "attr short", NC_SHORT, 3, attrBuf);
        CHECK_ERR
    }

    err = ncmpi_enddef(ncid);
    CHECK_ERR

    intBuf[0] = rank;
    for (i=0; i<LARGE_NUM; i++) {
        err = ncmpi_put_var_int_all(ncid, varid[i], intBuf);
        CHECK_ERR
    }

    err = ncmpi_close(ncid);
    CHECK_ERR

    err = ncmpi_open(comm, out_path, NC_WRITE, info, &ncid);
    CHECK_ERR

    err = ncmpi_redef(ncid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    for (i=0; i<LARGE_NUM; i++) {
        MPI_Offset len;
        err = ncmpi_inq_dim(ncid, i, str, &len);
        CHECK_ERR
    }

    for (i=0; i<LARGE_NUM; i++) {
        int buf;
        sprintf(str, "attr%d", i);
        err = ncmpi_get_att(ncid, 0, str, &buf);
        CHECK_ERR
    }

    for (i=0; i<LARGE_NUM; i++) {
        nc_type xtype;
        int ndims, natts;
        err = ncmpi_inq_var(ncid, i, str, &xtype, &ndims, dimids, &natts);
        CHECK_ERR
    }

    err = ncmpi_close(ncid); CHECK_ERR

    free(varid);
    free(dimids);

err_out:
    if (comm != MPI_COMM_WORLD && comm != MPI_COMM_NULL)
        MPI_Comm_free(&comm);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 0; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "large DIMS, VARS, ATTRS", opt, test_io);

    MPI_Finalize();

    return err;
}
