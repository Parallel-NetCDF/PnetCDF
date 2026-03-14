/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 *
 *  Check if arguments start, count, stride, and imap are properly ignored
 *  when get/put a scalar variable.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>


static int
tst_fmt(const char *out_path, int format, int coll_io, MPI_Info info)
{
    int err, ncid, varid, buf;
    MPI_Offset start[1], count[1], stride[1], imap[1];

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define a scalar variable of integer type */
    err = ncmpi_def_var(ncid, "scalar_var", NC_INT, 0, NULL, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    buf = 1;
    start[0] = 1;
    count[0] = 2;
    stride[0] = 2;
    imap[0] = 2;

    /* put */
    if (coll_io) {
        err = ncmpi_put_var1_int_all(ncid, varid, NULL,  &buf); CHECK_ERR
        err = ncmpi_put_var1_int_all(ncid, varid, start, &buf); CHECK_ERR

        err = ncmpi_put_vara_int_all(ncid, varid, start, count, &buf); CHECK_ERR
        err = ncmpi_put_vara_int_all(ncid, varid, NULL, count, &buf); CHECK_ERR
        err = ncmpi_put_vara_int_all(ncid, varid, start, NULL, &buf); CHECK_ERR
        err = ncmpi_put_vara_int_all(ncid, varid, NULL, NULL, &buf); CHECK_ERR

        err = ncmpi_put_vars_int_all(ncid, varid, start, count, stride, &buf); CHECK_ERR
        err = ncmpi_put_vars_int_all(ncid, varid, NULL, count, stride, &buf); CHECK_ERR
        err = ncmpi_put_vars_int_all(ncid, varid, start, NULL, stride, &buf); CHECK_ERR
        err = ncmpi_put_vars_int_all(ncid, varid, start, count, NULL, &buf); CHECK_ERR
        err = ncmpi_put_vars_int_all(ncid, varid, NULL, NULL, NULL, &buf); CHECK_ERR

        err = ncmpi_put_varm_int_all(ncid, varid, start, count, stride, imap, &buf); CHECK_ERR
        err = ncmpi_put_varm_int_all(ncid, varid, NULL, NULL, NULL, NULL, &buf); CHECK_ERR
    }
    else {
        err = ncmpi_put_var1_int(ncid, varid, NULL,  &buf); CHECK_ERR
        err = ncmpi_put_var1_int(ncid, varid, start, &buf); CHECK_ERR

        err = ncmpi_put_vara_int(ncid, varid, start, count, &buf); CHECK_ERR
        err = ncmpi_put_vara_int(ncid, varid, NULL, count, &buf); CHECK_ERR
        err = ncmpi_put_vara_int(ncid, varid, start, NULL, &buf); CHECK_ERR
        err = ncmpi_put_vara_int(ncid, varid, NULL, NULL, &buf); CHECK_ERR

        err = ncmpi_put_vars_int(ncid, varid, start, count, stride, &buf); CHECK_ERR
        err = ncmpi_put_vars_int(ncid, varid, NULL, count, stride, &buf); CHECK_ERR
        err = ncmpi_put_vars_int(ncid, varid, start, NULL, stride, &buf); CHECK_ERR
        err = ncmpi_put_vars_int(ncid, varid, start, count, NULL, &buf); CHECK_ERR
        err = ncmpi_put_vars_int(ncid, varid, NULL, NULL, NULL, &buf); CHECK_ERR

        err = ncmpi_put_varm_int(ncid, varid, start, count, stride, imap, &buf); CHECK_ERR
        err = ncmpi_put_varm_int(ncid, varid, NULL, NULL, NULL, NULL, &buf); CHECK_ERR
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_varid(ncid, "scalar_var", &varid); CHECK_ERR

    /* get */
    if (coll_io) {
        err = ncmpi_get_var1_int_all(ncid, varid, NULL,  &buf); CHECK_ERR
        err = ncmpi_get_var1_int_all(ncid, varid, start, &buf); CHECK_ERR

        err = ncmpi_get_vara_int_all(ncid, varid, start, count, &buf); CHECK_ERR
        err = ncmpi_get_vara_int_all(ncid, varid, NULL, count, &buf); CHECK_ERR
        err = ncmpi_get_vara_int_all(ncid, varid, start, NULL, &buf); CHECK_ERR
        err = ncmpi_get_vara_int_all(ncid, varid, NULL, NULL, &buf); CHECK_ERR

        err = ncmpi_get_vars_int_all(ncid, varid, start, count, stride, &buf); CHECK_ERR
        err = ncmpi_get_vars_int_all(ncid, varid, NULL, count, stride, &buf); CHECK_ERR
        err = ncmpi_get_vars_int_all(ncid, varid, start, NULL, stride, &buf); CHECK_ERR
        err = ncmpi_get_vars_int_all(ncid, varid, start, count, NULL, &buf); CHECK_ERR
        err = ncmpi_get_vars_int_all(ncid, varid, NULL, NULL, NULL, &buf); CHECK_ERR

        err = ncmpi_get_varm_int_all(ncid, varid, start, count, stride, imap, &buf); CHECK_ERR
        err = ncmpi_get_varm_int_all(ncid, varid, NULL, NULL, NULL, NULL, &buf); CHECK_ERR
    }
    else {
        err = ncmpi_get_var1_int(ncid, varid, NULL,  &buf); CHECK_ERR
        err = ncmpi_get_var1_int(ncid, varid, start, &buf); CHECK_ERR

        err = ncmpi_get_vara_int(ncid, varid, start, count, &buf); CHECK_ERR
        err = ncmpi_get_vara_int(ncid, varid, NULL, count, &buf); CHECK_ERR
        err = ncmpi_get_vara_int(ncid, varid, start, NULL, &buf); CHECK_ERR
        err = ncmpi_get_vara_int(ncid, varid, NULL, NULL, &buf); CHECK_ERR

        err = ncmpi_get_vars_int(ncid, varid, start, count, stride, &buf); CHECK_ERR
        err = ncmpi_get_vars_int(ncid, varid, NULL, count, stride, &buf); CHECK_ERR
        err = ncmpi_get_vars_int(ncid, varid, start, NULL, stride, &buf); CHECK_ERR
        err = ncmpi_get_vars_int(ncid, varid, start, count, NULL, &buf); CHECK_ERR
        err = ncmpi_get_vars_int(ncid, varid, NULL, NULL, NULL, &buf); CHECK_ERR

        err = ncmpi_get_varm_int(ncid, varid, start, count, stride, imap, &buf); CHECK_ERR
        err = ncmpi_get_varm_int(ncid, varid, NULL, NULL, NULL, NULL, &buf); CHECK_ERR
    }

    err = ncmpi_close(ncid); CHECK_ERR

    return 0;
}

#define WAIT_CHECK { \
    CHECK_ERR \
    if (coll_io) \
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); \
    else \
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL); \
    CHECK_ERR \
}

static int
tst_fmt_nb(const char *out_path, int format, int coll_io, MPI_Info info)
{
    int err, ncid, varid, buf;
    MPI_Offset start[1], count[1], stride[1], imap[1];

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define a scalar variable of integer type */
    err = ncmpi_def_var(ncid, "scalar_var", NC_INT, 0, NULL, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    buf = 1;
    start[0] = 1;
    count[0] = 2;
    stride[0] = 2;
    imap[0] = 2;

    /* put */
    err = ncmpi_iput_var1_int(ncid, varid, NULL,  &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_var1_int(ncid, varid, start, &buf, NULL); WAIT_CHECK

    err = ncmpi_iput_vara_int(ncid, varid, start, count, &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_vara_int(ncid, varid, NULL, count, &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_vara_int(ncid, varid, start, NULL, &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_vara_int(ncid, varid, NULL, NULL, &buf, NULL); WAIT_CHECK

    err = ncmpi_iput_vars_int(ncid, varid, start, count, stride, &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_vars_int(ncid, varid, NULL, count, stride, &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_vars_int(ncid, varid, start, NULL, stride, &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_vars_int(ncid, varid, start, count, NULL, &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_vars_int(ncid, varid, NULL, NULL, NULL, &buf, NULL); WAIT_CHECK

    err = ncmpi_iput_varm_int(ncid, varid, start, count, stride, imap, &buf, NULL); WAIT_CHECK
    err = ncmpi_iput_varm_int(ncid, varid, NULL, NULL, NULL, NULL, &buf, NULL); WAIT_CHECK

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_varid(ncid, "scalar_var", &varid); CHECK_ERR

    /* get */
    err = ncmpi_iget_var1_int(ncid, varid, NULL,  &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_var1_int(ncid, varid, start, &buf, NULL); WAIT_CHECK

    err = ncmpi_iget_vara_int(ncid, varid, start, count, &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_vara_int(ncid, varid, NULL, count, &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_vara_int(ncid, varid, start, NULL, &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_vara_int(ncid, varid, NULL, NULL, &buf, NULL); WAIT_CHECK

    err = ncmpi_iget_vars_int(ncid, varid, start, count, stride, &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_vars_int(ncid, varid, NULL, count, stride, &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_vars_int(ncid, varid, start, NULL, stride, &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_vars_int(ncid, varid, start, count, NULL, &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_vars_int(ncid, varid, NULL, NULL, NULL, &buf, NULL); WAIT_CHECK

    err = ncmpi_iget_varm_int(ncid, varid, start, count, stride, imap, &buf, NULL); WAIT_CHECK
    err = ncmpi_iget_varm_int(ncid, varid, NULL, NULL, NULL, NULL, &buf, NULL); WAIT_CHECK

    err = ncmpi_close(ncid); CHECK_ERR

    return 0;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char val[MPI_MAX_INFO_VAL];
    int nerrs=0, flag;

#ifdef DEBUG
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs > 1 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n", argv[0]);
#endif

    /* check whether burst buffering is enabled */
    MPI_Info_get(info, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, val, &flag);
    if (flag && strcasecmp(val, "enable") == 0 &&
        (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC))
        /* does not work for NetCDF4 files when burst-buffering is enabled */
        return 0;

    /* test blocking APIs */
    nerrs += tst_fmt(out_path, format, coll_io, info);

    /* test nonblocking APIs */
    if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC)
        nerrs += tst_fmt_nb(out_path, format, coll_io, info);

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
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "get/put scalar variables", opt, test_io);

    MPI_Finalize();

    return err;
}
