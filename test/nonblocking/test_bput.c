/*********************************************************************
 *
 *  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/*  $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <assert.h>
#include <pnetcdf.h>

#include <testutils.h>

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char hint[MPI_MAX_INFO_VAL];
    int i, j, ncid, dimid[2], varid, err, nerrs=0, rank, bb_enabled;
    int flag, req[2], status[2];
    float  var[4][6];
    MPI_Offset bufsize,  start[2], count[2], stride[2], imap[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef DEBUG
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs > 1 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n", argv[0]);
#endif

    if (rank) goto err_out;

    MPI_Info_get(info, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, hint, &flag);
    if (flag && strcasecmp(hint, "enable") == 0)
        bb_enabled = 1;
    else
        bb_enabled = 0;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_SELF, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define a variable of a 6 x 4 integer array in the nc file */
    err = ncmpi_def_dim(ncid, "Y", 6, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 4, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* set the contents of the write buffer var, a 4 x 6 float array
          50, 51, 52, 53, 54, 55,
          56, 57, 58, 59, 60, 61,
          62, 63, 64, 65, 66, 67,
          68, 69, 70, 71, 72, 73
     */
    for (j=0; j<4; j++) for (i=0; i<6; i++) var[j][i] = 50.5 + j*6+i;

    /* bufsize must be max of data type converted before and after */
    bufsize = 4*6*sizeof(int);
    err = ncmpi_buffer_attach(ncid, bufsize); CHECK_ERR

    /* write var to the NC variable in the matrix transposed way */
    count[0]  = 6; count[1]  = 2;
    stride[0] = 1; stride[1] = 1;
    imap[0]   = 1; imap[1]   = 6;   /* would be {4, 1} if not transposing */

    if (rank > 0) /* non-root processes just participate the call */
        count[0] = count[1] = 0;

    /* write the first two columns of the NC variable in the matrix transposed way */
    start[0]  = 0; start[1]  = 0;
    err = ncmpi_bput_varm_float(ncid, varid, start, count, stride, imap, &var[0][0], &req[0]); CHECK_ERR

    /* check if write buffer contents have been altered */
    for (j=0; j<4; j++)
        for (i=0; i<6; i++) {
            if (var[j][i] != 50.5 + j*6+i) {
                printf("Error at line %d in %s: put buffer[%d][%d]=%f altered, should be %f\n",
                       __LINE__,__FILE__,j,i,var[j][i],50.5+j*6+i);
                assert(0);
            }
        }

    /* write the second two columns of the NC variable in the matrix transposed way */
    start[0]  = 0; start[1]  = 2;
    err = ncmpi_bput_varm_float(ncid, varid, start, count, stride, imap, &var[2][0], &req[1]); CHECK_ERR

    /* check if write buffer contents have been altered */
    for (j=0; j<4; j++)
        for (i=0; i<6; i++) {
            if (var[j][i] != 50.5 + j*6+i) {
                printf("Error at line %d in %s: put buffer[%d][%d]=%f altered, should be %f\n",
                       __LINE__,__FILE__,j,i,var[j][i],50.5+j*6+i);
                assert(0);
            }
        }

    if (coll_io)
        err = ncmpi_wait_all(ncid, 2, req, status);
    else
        err = ncmpi_wait(ncid, 2, req, status);
    CHECK_ERR

    /* check each bput status */
    for (i=0; i<2; i++) {
        err = status[i];
        CHECK_ERR
    }

    err = ncmpi_buffer_detach(ncid); CHECK_ERR

    /* the output from command "ncmpidump -v var test.nc" should be:
           var =
            50, 56, 62, 68,
            51, 57, 63, 69,
            52, 58, 64, 70,
            53, 59, 65, 71,
            54, 60, 66, 72,
            55, 61, 67, 73 ;
     */

    /* check if the contents of write buffer have been altered (should not be) */
    for (j=0; j<4; j++) {
        for (i=0; i<6; i++) {
            if (var[j][i] != 50.5+j*6+i) {
                /* this error is a pnetcdf internal error, if occurs */
                printf("Error at line %d in %s: put buffer[%d][%d]=%f altered, should be %f\n",
                       __LINE__,__FILE__,j,i,var[j][i],50.5+j*6+i);
                assert(0);
            }
        }
    }

    if (!bb_enabled) {
        /* Try calling a bput after buffer detached. Expecting error */
        start[0] = 0; start[1] = 0;
        count[0] = 1; count[1] = 1;
        err = ncmpi_bput_vara_float(ncid, varid, start, count, &var[0][0], &req[0]);
        EXP_ERR(NC_ENULLABUF)
    }

    err = ncmpi_close(ncid); CHECK_ERR

err_out:
    MPI_Barrier(MPI_COMM_WORLD);

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

    err = tst_main(argc, argv, "bput API", opt, test_io);

    MPI_Finalize();

    return err;
}
