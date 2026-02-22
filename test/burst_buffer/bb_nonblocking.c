/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests nonblocking functionality of bb driver
 * Flushed requests can not be canceled, canceled request can't be waited
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <limits.h>
#include <testutils.h>
#include <libgen.h>

static
int test_bb(const char *out_path,
            MPI_Info    info)
{
    char *folder, *dup_out_path;
    int err, tmp, nerrs = 0, rank, np, ncid, varid, dimid[2];
    int buffer, req1, req2, stat;
    MPI_Offset start[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    /* add file info */
    MPI_Info_set(info, "nc_burst_buf", "enable");
    MPI_Info_set(info, "nc_burst_buf_overwrite", "enable");

    dup_out_path = strdup(out_path);
    folder = dirname(dup_out_path);
    if (folder == NULL)
        MPI_Info_set(info, "nc_burst_buf_dirname", ".");
    else
        MPI_Info_set(info, "nc_burst_buf_dirname", folder);
    free(dup_out_path);

    /* Create new netcdf file */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* Define dimensions */
    err = ncmpi_def_dim(ncid, "X", np, dimid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", 4, dimid + 1); CHECK_ERR

    /* Define variable */
    err = ncmpi_def_var(ncid, "M", NC_INT, 2, dimid, &varid); CHECK_ERR

    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

    /* Switch to data mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    buffer = rank + 1;

    /* each process writes to a different row */
    start[0] = rank;
    start[1] = 0;
    err = ncmpi_iput_var1_int(ncid, varid, start, &buffer, &req1); CHECK_ERR
    start[1] = 1;
    err = ncmpi_iput_var1_int(ncid, varid, start, &buffer, &req2); CHECK_ERR
    start[1] = 0;
    err = ncmpi_get_var1_int_all(ncid, varid, start, &buffer); CHECK_ERR
    start[1] = 1;
    err = ncmpi_get_var1_int_all(ncid, varid, start, &buffer); CHECK_ERR
    err = ncmpi_cancel(ncid, 1, &req1, &stat); CHECK_ERR
    err = stat;    EXP_ERR(NC_EFLUSHED)
    err = ncmpi_wait_all(ncid, 1, &req2, &stat); CHECK_ERR
    err = stat; CHECK_ERR

    start[1] = 2;
    err = ncmpi_iput_var1_int(ncid, varid, start, &buffer, &req1); CHECK_ERR
    start[1] = 3;
    err = ncmpi_iput_var1_int(ncid, varid, start, &buffer, &req2); CHECK_ERR
    tmp = req1;
    err = ncmpi_cancel(ncid, 1, &req1, &stat); CHECK_ERR
    err = stat; CHECK_ERR
    start[1] = 2;
    err = ncmpi_get_var1_int_all(ncid, varid, start, &buffer); CHECK_ERR
    start[1] = 3;
    err = ncmpi_get_var1_int_all(ncid, varid, start, &buffer); CHECK_ERR
    req1 = tmp;
    err = ncmpi_wait_all(ncid, 1, &req1, &stat); CHECK_ERR
    err = stat;    EXP_ERR(NC_EINVAL_REQUEST)
    err = ncmpi_wait_all(ncid, 1, &req2, &stat); CHECK_ERR
    err = stat; CHECK_ERR

    /* Close the file */
    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int err=NC_NOERR;
    MPI_Info local_info;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    MPI_Info_dup(info, &local_info);
    err = test_bb(out_path, local_info);
    MPI_Info_free(&local_info);

    MPI_Info_dup(info, &local_info);
    MPI_Info_set(local_info, "nc_burst_buf_shared_logs", "enable");
    err = test_bb(out_path, local_info);
    MPI_Info_free(&local_info);

    MPI_Barrier(MPI_COMM_WORLD);
    return err;
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
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "nonblocking APIs", opt, test_io);

    MPI_Finalize();

    return err;
}
