/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This example tests if ncmpi_end_indep_data() works properly when nonblocking
 * APIs are called in the independent data mode, but the wait call is made later
 * in collective data mode.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 4
#define NX 10
#define NDIMS 2

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, varid, dimid[2], req, st;
    MPI_Offset start[2], count[2], stride[2];
    signed char buffer[NY][NX];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs,    &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_BYTE, NDIMS, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    for (i=0; i<NY; i++) for (j=0; j<NX; j++) buffer[i][j] = rank+10;

     start[0] = 0;   start[1] = rank;
     count[0] = NY;  count[1] = NX;
    stride[0] = 1;  stride[1] = nprocs;
    err = ncmpi_buffer_attach(ncid, NY*NX); CHECK_ERR

    err = ncmpi_begin_indep_data(ncid); CHECK_ERR
    err = ncmpi_bput_vars_schar(ncid, varid, start, count, stride,
                                &buffer[0][0], &req);
    CHECK_ERR

    /* check if write buffer contents have been altered */
    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++) {
            if (buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: put buffer[%d][%d]=%hhu altered, should be %d\n",
                       __LINE__,__FILE__,i,j,buffer[i][j],rank+10);
                nerrs++;
            }
        }

    err = ncmpi_end_indep_data(ncid); CHECK_ERR

    /* calling wait API after exiting independent data mode on purpose */
    err = ncmpi_wait_all(ncid, 1, &req, &st); CHECK_ERR
    err = st; CHECK_ERR

    /* check if write buffer contents have been altered */
    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++) {
            if (buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: put buffer[%d][%d]=%hhu altered, should be %d\n",
                       __LINE__,__FILE__,i,j,buffer[i][j],rank+10);
                nerrs++;
            }
        }

    err = ncmpi_buffer_detach(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

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
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "ncmpi_end_indep_data()", opt, test_io);

    MPI_Finalize();

    return err;
}
