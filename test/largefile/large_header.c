/*********************************************************************
 *
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program is to test
 *
 * large header size, i.e. > INT_MAX, i.e. 2 GiB
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <limits.h> /* INT_MAX */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int rank, nprocs, err, nerrs=0;
    int ncid, dimid, varid, buf;
    MPI_Offset extent, start;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define a dimension of size = nprocs */
    err = ncmpi_def_dim(ncid, "dim", nprocs, &dimid);
    CHECK_ERR

    /* define a variable */
    err = ncmpi_def_var(ncid, "var0", NC_INT, 1, &dimid, &varid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, &dimid, &varid);
    CHECK_ERR

    /* make file header extent > 4 GiB */
    extent = (MPI_Offset)INT_MAX + 1024;
    err = ncmpi__enddef(ncid, 0, extent, 0, 0);
    CHECK_ERR

    /* write to the variable */
    start = rank;
    buf = rank;
    err = ncmpi_put_var1_int_all(ncid, varid, &start, &buf);
    CHECK_ERR

    err = ncmpi_close(ncid);
    CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    /* inquire ID of the variable */
    err = ncmpi_inq_varid(ncid, "var1", &varid);
    CHECK_ERR

    /* read from the variable */
    buf = -1;
    err = ncmpi_get_var1_int_all(ncid, varid, &start, &buf);
    CHECK_ERR

    if (buf != rank) {
        printf("Error at line %d in %s: expecting read buf %d but got %d\n",
               __LINE__,__FILE__,rank,buf);
        nerrs++;
    }

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_64BIT_DATA};
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
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "large header", opt, test_io);

    MPI_Finalize();

    return err;
}
