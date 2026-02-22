/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
   This program tests ncmpi_create(), ncmpi_open(), and ncmpi_close().
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

static
int file_op(const char *out_path,
            MPI_Info    info)
{
    int err, nerrs=0, ncid;

    /* Create a new file */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* Close the file. */
    err = ncmpi_close(ncid);
    CHECK_ERR

    /* Open the file */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_WRITE, info, &ncid);
    CHECK_ERR

    /* Close the file. */
    err = ncmpi_close(ncid);
    CHECK_ERR

    return (nerrs > 0);
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int err, nerrs=0;
    MPI_Info info_dup=MPI_INFO_NULL;

    MPI_Info_dup(info, &info_dup);

    MPI_Info_set(info_dup, "nc_header_align_size", "100");
    MPI_Info_set(info_dup, "nc_var_align_size", "200");
    MPI_Info_set(info_dup, "nc_record_align_size", "300");

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = file_op(out_path, info_dup);
    if (err != 0) nerrs++;

#ifdef MPICH_VERSION
    /* MPICH recognizes file system type acronym prefixed to the file name */
    if (strncmp("ufs:", out_path, 4)) {
        char *prefix_fname;
        prefix_fname = (char*) malloc(strlen(out_path) + 10);
        sprintf(prefix_fname, "ufs:%s", out_path);

        err = file_op(prefix_fname, info_dup);
        if (err != 0) nerrs++;

        free(prefix_fname);
    }
#endif

    if (info != MPI_INFO_NULL) MPI_Info_free(&info_dup);

    return (nerrs > 0);
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
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "file create/open", opt, test_io);

    MPI_Finalize();

    return err;
}

