/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This program tests a call to ncmpi_inq_header_size() when in the define
 * mode, which should calculate and return the latest file header size. This
 * can be useful for application users to decide how much free space to be
 * preserved in the file header section, i.e. by setting argument h_minfree
 * and/or v_align when calling ncmpi__enddef().
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

static int debug;

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    char *str;
    int  err, nerrs=0, rank, ncid, dimids[2], varid, int_buf;
    float flt_buf;
    double *dbl_buf;
    MPI_Offset old_h_size, old_h_extent, new_h_size, new_h_extent;

    debug = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a file */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    flt_buf = 1.234;
    err = ncmpi_put_att(ncid, NC_GLOBAL, "_FillValue", NC_FLOAT, 1, &flt_buf);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "X", 10, &dimids[0]); CHECK_ERR
    CHECK_ERR

    err = ncmpi_def_var(ncid, "int_var", NC_INT, 1, dimids, &varid);
    CHECK_ERR

    err = ncmpi_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, &flt_buf);
    EXP_ERR(NC_EBADTYPE)

    int_buf = 5678;
    err = ncmpi_put_att(ncid, varid, "_FillValue", NC_INT, 1, &int_buf);
    CHECK_ERR

    err = ncmpi_def_var(ncid, "dbl_var", NC_DOUBLE, 1, dimids, &varid);
    CHECK_ERR

    err = ncmpi_def_var(ncid, "short_var", NC_SHORT, 1, dimids, &varid);
    CHECK_ERR

    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &old_h_size); CHECK_ERR

    err = ncmpi_inq_header_extent(ncid, &old_h_extent); CHECK_ERR

    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld\n", __FILE__,__LINE__,
               old_h_size, old_h_extent);

    if (old_h_extent != 0) {
        printf("Error at %d: expect file extent size to be 0 but got %lld\n",
               __LINE__, old_h_extent);
        nerrs++;
        goto err_out;
    }

    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &new_h_size); CHECK_ERR

    err = ncmpi_inq_header_extent(ncid, &new_h_extent); CHECK_ERR

    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld\n", __FILE__,__LINE__,
               new_h_size, new_h_extent);

    if (new_h_size != old_h_size) {
        printf("Error at %d: expect file header size %lld but got %lld\n",
               __LINE__, old_h_size, new_h_size);
        nerrs++;
        goto err_out;
    }

    if (new_h_extent <= old_h_extent) {
        printf("Error at %d: expect file extent size > %lld but got %lld\n",
               __LINE__, old_h_extent, new_h_extent);
        nerrs++;
        goto err_out;
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    old_h_size = new_h_size;
    old_h_extent = new_h_extent;

    /* open the file */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_WRITE, info, &ncid);
    CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &new_h_size); CHECK_ERR

    err = ncmpi_inq_header_extent(ncid, &new_h_extent); CHECK_ERR

    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld\n", __FILE__,__LINE__,
               new_h_size, new_h_extent);

    if (new_h_size != old_h_size) {
        printf("Error at %d: expect file header size %lld but got %lld\n",
               __LINE__, old_h_size, new_h_size);
        nerrs++;
        goto err_out;
    }

    if (new_h_extent != old_h_extent) {
        printf("Error at %d: expect file extent size %lld but got %lld\n",
               __LINE__, old_h_extent, new_h_extent);
        nerrs++;
        goto err_out;
    }

    old_h_size = new_h_size;
    old_h_extent = new_h_extent;

    /* enter define mode and add new a dimension and a variable */
    err = ncmpi_redef(ncid); CHECK_ERR

    str = "new global attribute of text data type";
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "global_attr", strlen(str), str);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "Y", 10, &dimids[1]); CHECK_ERR
    CHECK_ERR

    err = ncmpi_def_var(ncid, "new_int_var", NC_INT, 2, dimids, &varid);
    CHECK_ERR

    dbl_buf = (double*) calloc(16, sizeof(double));
    err = ncmpi_put_att_double(ncid, varid, "attr", NC_DOUBLE, 16, dbl_buf);
    CHECK_ERR
    free(dbl_buf);

    err = ncmpi_inq_header_size(ncid, &new_h_size); CHECK_ERR

    err = ncmpi_inq_header_extent(ncid, &new_h_extent); CHECK_ERR

    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld\n", __FILE__,__LINE__,
               new_h_size, new_h_extent);

    if (new_h_size <= old_h_size) {
        printf("Error at %d: expect file header size > %lld but got %lld\n",
               __LINE__, old_h_size, new_h_size);
        nerrs++;
        goto err_out;
    }

    if (new_h_extent != old_h_extent) {
        printf("Error at %d: expect file extent size %lld but got %lld\n",
               __LINE__, old_h_extent, new_h_extent);
        nerrs++;
        goto err_out;
    }

    if (new_h_size > old_h_extent)
        /* header size grows beyond the current file extent size */
        err = ncmpi__enddef(ncid, 0, 512, 0, 0);
    else
        err = ncmpi_enddef(ncid);
    CHECK_ERR

    old_h_size = new_h_size;
    old_h_extent = new_h_extent;

    err = ncmpi_inq_header_size(ncid, &new_h_size); CHECK_ERR

    err = ncmpi_inq_header_extent(ncid, &new_h_extent); CHECK_ERR

    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld\n", __FILE__,__LINE__,
               new_h_size, new_h_extent);

    if (new_h_size != old_h_size) {
        printf("Error at %d: expect file header size %lld but got %lld\n",
               __LINE__, old_h_size, new_h_size);
        nerrs++;
        goto err_out;
    }

    if (new_h_extent < old_h_extent) {
        printf("Error at %d: expect file extent size >= %lld but got %lld\n",
               __LINE__, old_h_extent, new_h_extent);
        nerrs++;
        goto err_out;
    }

    err = ncmpi_close(ncid); CHECK_ERR

err_out:
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
    opt.chk      = 1; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "in define mode", opt, test_io);

    MPI_Finalize();

    return err;
}
