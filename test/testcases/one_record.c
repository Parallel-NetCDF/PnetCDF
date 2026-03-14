/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests the special case of ONLY one record variable is defined
 * and the record size is not aligned with the 4-byte boundary. As defined in
 * CDF-1 and CDF-2 format specifications:
 *    "A special case: Where there is exactly one record variable, we drop the
 *    requirement that each record be four-byte aligned, so in this case there
 *    is no record padding."
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

#define STR_LEN 19
#define NUM_VALS 2

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, err, nerrs=0, rank, nprocs;
    int ncid, dimids[2], varid;
    char data[NUM_VALS][STR_LEN + 1], data_in[NUM_VALS*STR_LEN];
    MPI_Offset start[2];
    MPI_Offset count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    strcpy(data[0], "2005-04-11_12:00:00"); /* 19 bytes not a multiply of 4 */
    strcpy(data[1], "2005-04-11_13:00:00");

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err  = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, dimids); CHECK_ERR
    err = ncmpi_def_dim(ncid, "text_dim", STR_LEN, &dimids[1]); CHECK_ERR

    /* create ONLY one record variable of type NC_CHAR and make sure each
     * record is of size not aligned with 4-byte boundary.
     */
    err = ncmpi_def_var(ncid, "text_var", NC_CHAR, 2, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }
    /* Write some records of var data. */
    count[0] = 1;
    count[1] = STR_LEN;
    start[0] = 0;
    start[1] = 0;
    for (i=0; i<NUM_VALS; i++) {
        if (coll_io)
            err = ncmpi_put_vara_text_all(ncid, varid, start, count, data[start[0]]);
        else
            err = ncmpi_put_vara_text(ncid, varid, start, count, data[start[0]]);
        CHECK_ERR
        start[0]++;
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    err  = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }
    err = ncmpi_inq_varid(ncid, "text_var", &varid); CHECK_ERR

    /* read the entire data back */
    if (coll_io)
        err = ncmpi_get_var_text_all(ncid, varid, data_in);
    else
        err = ncmpi_get_var_text(ncid, varid, data_in);
    CHECK_ERR

    /* check the contents */
    for (i=0; i<NUM_VALS; i++)
      if (strncmp(data[i], data_in+i*STR_LEN, STR_LEN)) {
          printf("Error at line %d in %s: expecting %s but got %s\n",
          __LINE__,__FILE__,data[i],data_in+i*STR_LEN);
          nerrs++;
      }

    err = ncmpi_close(ncid); CHECK_ERR

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

    err = tst_main(argc, argv, "only one record variable", opt, test_io);

    MPI_Finalize();

    return err;
}
