/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program adds two new variables to an existing netCDF file.
 * It is used to test if PnetCDF can correctly calculate the file offsets
 * for the two new variables, in particular for files that align the
 * fixed-size variables to a boundary larger than 4 bytes, for instance
 * a file created by PnetCDF with default alignment of 512 bytes.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o add_var add_var.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 add_var testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* dirname() */
#include <pnetcdf.h>

#include <testutils.h>

static
int test_bb(const char *out_path,
            MPI_Info    info)
{
    char *folder, *dup_out_path;
    int err, flag, nerrs=0, ncid;
    MPI_Info infoused;
    char hint[MPI_MAX_INFO_VAL];

    MPI_Info_set(info, "nc_burst_buf", "enable");
    MPI_Info_set(info, "nc_burst_buf_del_on_close", "disable");
    MPI_Info_set(info, "nc_burst_buf_flush_buffer_size", "256");
    /* MPI_Info_set(info, "nc_burst_buf_dirname", "()@^$@!(_&$)@(#%%&)(*#$"); */

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

    err = ncmpi_inq_file_info(ncid, &infoused); CHECK_ERR

    MPI_Info_get(infoused, "nc_burst_buf_del_on_close", MPI_MAX_INFO_VAL - 1, hint, &flag);
    if (flag) {
        if (strcmp(hint, "disable") != 0) {
            printf("Error at line %d: unexpected nc_burst_buf_del_on_close = %s, but got %s\n",
                    __LINE__, "disable", hint);
            nerrs++;
            goto err_out;
        }
    }
    else{
        printf("Error at line %d: nc_burst_buf_del_on_close is not set\n", __LINE__);
        nerrs++;
        goto err_out;
    }

    MPI_Info_get(infoused, "nc_burst_buf_flush_buffer_size", MPI_MAX_INFO_VAL - 1, hint, &flag);
    if (flag) {
        if (strcmp(hint, "256") != 0) {
            printf("Error at line %d: unexpected nc_burst_buf_flush_buffer_size = %s, but got %s\n",
                    __LINE__, "256", hint);
            nerrs++;
            goto err_out;
        }
    }
    else{
        printf("Error at line %d: nc_burst_buf_flush_buffer_size is not set\n", __LINE__);
        nerrs++;
        goto err_out;
    }

    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    MPI_Info_free(&infoused);

err_out:
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
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "burt buffering hints", opt, test_io);

    MPI_Finalize();

    return err;
}
