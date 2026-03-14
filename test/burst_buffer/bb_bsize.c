/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests writing a variable exceeding the size of data buffer.
 * Each process writes a submatrix of size 1024 * 1024, a total of 1M cells
 * The submatrix form each processes is stacked along the first dimension so
 * that the variable dimensions are (1024 * np) * 1024 Each processes writes
 * it's rank to every cell
 * Each process starts by writing the first 1/8 of the rows at once, which
 * should caused an increase of the data buffer size to accommodate it Then each
 * process writes the remaining row one at a time
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <limits.h>
#include <testutils.h>
#include <libgen.h> /* dirname() */

#define SIZE 1024

int buffer[SIZE * SIZE];
char bsize[32];

static
int test_bb(const char *out_path,
            int         coll_io,
            MPI_Info    info)
{
    char *folder, *dup_out_path;
    int i, err = NC_NOERR, nerrs = 0;
    int rank, np, ncid, varid, dimid[2];
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    /* add file info */
    MPI_Info_set(info, "nc_burst_buf", "enable");

    /* Set default buffer size to 1/16 of the rows */
    sprintf(bsize, "%u", (unsigned int)(SIZE * SIZE / 16 * sizeof(int)));
    MPI_Info_set(info, "nc_burst_buf_flush_buffer_size", bsize);

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
    err = ncmpi_def_dim(ncid, "X", SIZE * np, dimid);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", SIZE, dimid + 1);
    CHECK_ERR

    /* Define variable */
    err = ncmpi_def_var(ncid, "M", NC_INT, 2, dimid, &varid);
    CHECK_ERR

    /* Switch to data mode */
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* Initialize buffer */
    for (i = 0; i < SIZE * SIZE; i++) {
        buffer[i] = rank + 1;
    }

    /* Write first 1/8 of the rows */
    start[0] = SIZE * rank;
    start[1] = 0;
    count[0] = SIZE / 8;
    count[1] = SIZE;
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid, start, count, buffer);
    else
        err = ncmpi_put_vara_int(ncid, varid, start, count, buffer);
    CHECK_ERR

    /* Write remaining rows */
    start[0] = SIZE * rank + SIZE / 8;
    start[1] = 0;
    count[0] = 1;
    count[1] = SIZE;
    for (; start[0] < SIZE * (rank + 1); start[0]++) {
        if (coll_io)
            err = ncmpi_put_vara_int_all(ncid, varid, start, count, buffer);
        else
            err = ncmpi_put_vara_int(ncid, varid, start, count, buffer);
        CHECK_ERR
    }

    /*
     * Read it back
     * Flush on read is on by default so no additional action required
     */
    memset(buffer, 0, sizeof(buffer));
    start[0] = SIZE * rank;
    start[1] = 0;
    count[0] = SIZE;
    count[1] = SIZE;
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid, start, count, buffer);
    else
        err = ncmpi_get_vara_int(ncid, varid, start, count, buffer);
    CHECK_ERR

    /* Verify the result */
    for (i = 0; i < SIZE * SIZE; i++) {
        if (buffer[i] != rank + 1) {
            printf("Error at line %d in %s: expecting buffer[%d] = %d but got %d\n",
                    __LINE__, __FILE__, i, rank + 1, buffer[i]);
            nerrs++;
            goto err_out;
        }
    }

    /* Close the file */
    err = ncmpi_close(ncid);
    CHECK_ERR

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

    err = tst_main(argc, argv, "request size > buffer size", opt, test_io);

    MPI_Finalize();

    return err;
}
