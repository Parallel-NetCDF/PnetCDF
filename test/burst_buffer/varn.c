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
            int         coll_io,
            MPI_Info    info)
{
    char *folder, *dup_out_path;
    int i, err, nerrs = 0, rank, np, ncid, varid, dimid[2], buffer[10];
    MPI_Offset starts[10][2], counts[10][2];
    MPI_Offset *Starts[10], *Counts[10];

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
    err = ncmpi_def_dim(ncid, "Y", 10, dimid + 1); CHECK_ERR

    /* Define variable */
    err = ncmpi_def_var(ncid, "M", NC_INT, 2, dimid, &varid); CHECK_ERR

    /* Switch to data mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    for(i = 0; i < 10; i++){
        starts[i][0] = rank;
        starts[i][1] = i;
        counts[i][0] = 1;
        counts[i][1] = 1;
        Starts[i] = (MPI_Offset*)starts[i];
        Counts[i] = (MPI_Offset*)counts[i];
        buffer[i] = rank + i;
    }

    /* Standard varn */
    if (coll_io)
        err = ncmpi_put_varn_int_all(ncid, varid, 10, Starts, Counts, buffer);
    else
        err = ncmpi_put_varn_int(ncid, varid, 10, Starts, Counts, buffer);
    CHECK_ERR

    for (i=0; i<10; i++) buffer[0] = -1;
    if (coll_io)
        err = ncmpi_get_varn_int_all(ncid, varid, 10, Starts, Counts, buffer);
    else
        err = ncmpi_get_varn_int(ncid, varid, 10, Starts, Counts, buffer);
    CHECK_ERR
    for(i = 0; i < 10; i++){
        if (buffer[i] != rank + i){
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n",
                    __LINE__, __FILE__, i, rank + 1, buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

    /* NULL counts */
    if (coll_io)
        err = ncmpi_put_varn_int_all(ncid, varid, 10, Starts, NULL, buffer);
    else
        err = ncmpi_put_varn_int(ncid, varid, 10, Starts, NULL, buffer);
    CHECK_ERR

    for (i=0; i<10; i++) buffer[0] = -1;
    if (coll_io)
        err = ncmpi_get_varn_int_all(ncid, varid, 10, Starts, NULL, buffer);
    else
        err = ncmpi_get_varn_int(ncid, varid, 10, Starts, NULL, buffer);
    CHECK_ERR
    for(i = 0; i < 10; i++){
        if (buffer[i] != rank + i){
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n",
                    __LINE__, __FILE__, i, rank + 1, buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

    /* Partial NULL counts */
    for(i = 0; i < 10; i += 2){
        Counts[i] = (MPI_Offset*)counts[i];
    }
    if (coll_io)
        err = ncmpi_put_varn_int_all(ncid, varid, 10, Starts, Counts, buffer);
    else
        err = ncmpi_put_varn_int(ncid, varid, 10, Starts, Counts, buffer);
    CHECK_ERR

    for (i=0; i<10; i++) buffer[0] = -1;
    if (coll_io)
        err = ncmpi_get_varn_int_all(ncid, varid, 10, Starts, Counts, buffer);
    else
        err = ncmpi_get_varn_int(ncid, varid, 10, Starts, Counts, buffer);
    CHECK_ERR
    for(i = 0; i < 10; i++){
        if (buffer[i] != rank + i){
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n",
                    __LINE__, __FILE__, i, rank + 1, buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

    /* Close the file */
    err = ncmpi_close(ncid); CHECK_ERR

err_out:
    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int err=NC_NOERR;
    MPI_Info local_info;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    MPI_Info_dup(info, &local_info);
    err = test_bb(out_path, coll_io, local_info);
    MPI_Info_free(&local_info);

    MPI_Info_dup(info, &local_info);
    MPI_Info_set(local_info, "nc_burst_buf_shared_logs", "enable");
    err = test_bb(out_path, coll_io, local_info);
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

    err = tst_main(argc, argv, "varn API", opt, test_io);

    MPI_Finalize();

    return err;
}
