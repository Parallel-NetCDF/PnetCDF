/******************************************************************************
 *
 *  Copyright (C) 2026, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *****************************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests a collective write to a record variable in which one of
 * the processes has nothing to contribute, i.e. it passes count == 0. Every
 * process still calls the collective API, so the call is legal.
 *
 * Such a request must not opt out of the collective synchronization of numrecs
 * that follows the I/O. If it does, the processes that do have data to write
 * are left calling an MPI_Allreduce() alone and the communicator goes out of
 * step, which shows up as a wrong error code, as corruption of a later
 * collective call, or as a hang.
 *
 * Fixed-size variables take no part in the numrecs synchronization and are
 * therefore not affected; the record variable is the case worth guarding.
 *
 *    To compile:
 *        mpicc -O2 tst_zero_len_recvar.c -o tst_zero_len_recvar -lpnetcdf
 *
 *    % mpiexec -n 4 ./tst_zero_len_recvar testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

#define NDBL 4  /* number of doubles written per process */

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, err, nerrs=0, ncid, dimid[2], varid, rank, nprocs, is_empty;
    MPI_Offset start[2], count[2], numrecs;
    double buf[NDBL], rbuf[NDBL];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* the last process contributes nothing; with a single process there is
     * nobody left to write, so keep its request non-empty
     */
    is_empty = (nprocs > 1 && rank == nprocs - 1);

    for (i=0; i<NDBL; i++) buf[i] = (double)(rank * NDBL + i);

    err = ncmpi_set_default_format(format, NULL); CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "time",  NC_UNLIMITED,   &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "nelem", nprocs * NDBL,  &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "v", NC_DOUBLE, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_NOFILL, NULL); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* all processes write into record 0 */
    start[0] = 0;
    count[0] = 1;
    start[1] = (MPI_Offset)rank * NDBL;
    count[1] = NDBL;
    if (is_empty) { start[1] = 0; count[1] = 0; }

    if (coll_io) {
        err = ncmpi_put_vara_double_all(ncid, varid, start, count, buf);
        CHECK_ERR
    }
    else {
        err = ncmpi_begin_indep_data(ncid); CHECK_ERR
        err = ncmpi_put_vara_double(ncid, varid, start, count, buf); CHECK_ERR
        err = ncmpi_end_indep_data(ncid); CHECK_ERR
    }

    /* one record must have been created, and all processes must agree */
    err = ncmpi_inq_dimlen(ncid, dimid[0], &numrecs); CHECK_ERR
    if (numrecs != 1) {
        printf("Error at %s:%d: expecting numrecs 1 but got %lld\n",
               __FILE__, __LINE__, numrecs);
        nerrs++;
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* read back and check the contents written by the non-empty processes */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "v", &varid); CHECK_ERR

    if (!is_empty) {
        for (i=0; i<NDBL; i++) rbuf[i] = -1.0;
        err = ncmpi_get_vara_double_all(ncid, varid, start, count, rbuf); CHECK_ERR
        for (i=0; i<NDBL; i++) {
            if (rbuf[i] != buf[i]) {
                printf("Error at %s:%d: expecting v[%lld] %g but got %g\n",
                       __FILE__, __LINE__, start[1]+i, buf[i], rbuf[i]);
                nerrs++;
                break;
            }
        }
    }
    else {
        /* still a collective call, with nothing to read */
        err = ncmpi_get_vara_double_all(ncid, varid, start, count, rbuf); CHECK_ERR
    }

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET,
                     NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 2;    /* enable and disable intra-node aggregation */
    opt.drv      = 2;    /* test GIO and MPI-IO driver */
    opt.ibuf     = 0;    /* test default hint nc_ibuf_size only */
    opt.bb       = 2;    /* enable and disable burst-buffering feature */
    opt.mod      = 2;    /* test both collective and independent data mode */
    opt.hdr_diff = true; /* run ncmpidiff for file header */
    opt.var_diff = false;/* skip ncmpidiff for variables */

    err = tst_main(argc, argv, "zero-length request to a record variable", opt,
                   test_io);

    MPI_Finalize();
    return err;
}
