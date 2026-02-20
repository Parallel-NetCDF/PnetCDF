/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* This test program checks whether a collective API can be nicely aborted
 * without causing the program to hang. It runs on 2 processes.
 * One process deliberately produces an error (using an illegal start argument),
 * while the other does not.
 *
 * Note when in safe mode, all processes obtain the same error code and hence
 * can terminate nicely. However, when not in safe mode, collective APIs can
 * hang for the following erros: NC_EBADID, NC_EPERM, NC_EINDEFINE, NC_EINDEP,
 * and NC_ENOTINDEP. These errors are considered fatal and program should stop.
 * For other kinds of errors, the processes causing the error will continue to
 * participate the collective I/O with zero-length requests, so that the MPI
 * collective I/O can complete with all processes participating the call.
 *
 * This program tests error code NC_EINVALCOORDS, not in the above list.
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <assert.h>
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
    char *val;
    int rank, nproc, ncid, err, nerrs=0, varid, dimids[1], req, status, exp;
    int safe_mode=0;
    double buf[2];
    MPI_Offset start[1], count[1];
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    val = getenv("PNETCDF_SAFE_MODE");
    if (val != NULL && atoi(val) == 1) safe_mode = 1;
#ifdef PNETCDF_DEBUG
    if (val == NULL && !safe_mode) safe_mode = 1;
#endif

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* Create a 2 element vector of doubles */
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR_ALL
    err = ncmpi_def_dim(ncid, "dim", 2, &dimids[0]); CHECK_ERR_ALL
    err = ncmpi_def_var(ncid, "var", NC_DOUBLE, 1, dimids, &varid); CHECK_ERR_ALL
    err = ncmpi_enddef(ncid); CHECK_ERR_ALL

    if (rank == 0) {
        start[0] = 0;
        count[0] = 2;
    } else if (rank == 1) {
        if (is_relax_coord_bound())
            start[0] = 3; /* illegal for a start > defined shape */
        else
            start[0] = 2; /* illegal for a start >= defined shape */
        count[0] = 0;
    }
    else {
        start[0] = 0;
        count[0] = 0;
    }

    buf[0] = 1.0;
    buf[1] = 2.0;

    err = ncmpi_put_vara_all(ncid, varid, start, count, buf, count[0],
                             MPI_DOUBLE);
    if ((safe_mode && nproc > 1) || rank == 1) exp = NC_EINVALCOORDS;
    else                                       exp = NC_NOERR;
    CHECK_EXP_ERR_ALL(exp)

    /* check if user put buffer contents altered */
    if (buf[0] != 1.0) {
        printf("Error at line %d in %s: user put buffer[%d] altered from %f to %f\n",
               __LINE__,__FILE__,0, 1.0, buf[0]);
        nerrs++;
    }
    if (buf[1] != 2.0) {
        printf("Error at line %d in %s: user put buffer[%d] altered from %f to %f\n",
               __LINE__,__FILE__,1, 2.0, buf[1]);
        nerrs++;
    }

    err = ncmpi_put_vara_double_all(ncid, varid, start, count, buf);
    if ((safe_mode && nproc > 1) || rank == 1) exp = NC_EINVALCOORDS;
    else                                       exp = NC_NOERR;
    CHECK_EXP_ERR_ALL(exp)

    /* check if user put buffer contents altered */
    if (buf[0] != 1.0) {
        printf("Error at line %d in %s: user put buffer[%d] altered from %f to %f\n",
               __LINE__,__FILE__,0, 1.0, buf[0]);
        nerrs++;
    }
    if (buf[1] != 2.0) {
        printf("Error at line %d in %s: user put buffer[%d] altered from %f to %f\n",
               __LINE__,__FILE__,1, 2.0, buf[1]);
        nerrs++;
    }

    if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
        err = ncmpi_iput_vara_double(ncid, varid, start, count, buf, &req);
        exp = (rank == 1) ? NC_EINVALCOORDS : NC_NOERR;
        EXP_ERR(exp)

        err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERR_ALL

        /* check if user put buffer contents altered */
        if (buf[0] != 1.0) {
            printf("Error at line %d in %s: user put buffer[%d] altered from %f to %f\n",
                   __LINE__,__FILE__,0, 1.0, buf[0]);
            nerrs++;
        }
        if (buf[1] != 2.0) {
            printf("Error at line %d in %s: user put buffer[%d] altered from %f to %f\n",
                   __LINE__,__FILE__,1, 2.0, buf[1]);
            nerrs++;
        }
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_get_vara_all(ncid, varid, start, count,
			     buf, count[0], MPI_DOUBLE);
    if ((safe_mode && nproc > 1) || rank == 1) exp = NC_EINVALCOORDS;
    else                                       exp = NC_NOERR;
    CHECK_EXP_ERR_ALL(exp)

    err = ncmpi_get_vara_double_all(ncid, varid, start, count, buf);
    if ((safe_mode && nproc > 1) || rank == 1) exp = NC_EINVALCOORDS;
    else                                       exp = NC_NOERR;
    CHECK_EXP_ERR_ALL(exp)

    if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
        err = ncmpi_iget_vara_double(ncid, varid, start, count, buf, &req);
        exp = (rank == 1) ? NC_EINVALCOORDS : NC_NOERR;
        EXP_ERR(exp)

        err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERR_ALL
    }

    err = ncmpi_close(ncid); CHECK_ERR_ALL

err_out:
    return nerrs;
}

#if 0
{
    int err, nerrs=0;

    /* test in non-safe mode */
    setenv("PNETCDF_SAFE_MODE", "0", 1);
    nerrs += test_collective_error(filename, 0, 0);
    if (nerrs) goto err_out;
    nerrs += test_collective_error(filename, 0, NC_64BIT_OFFSET);
    if (nerrs) goto err_out;
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += test_collective_error(filename, 0, NC_NETCDF4);
        if (nerrs) goto err_out;
        nerrs += test_collective_error(filename, 0, NC_NETCDF4 | NC_CLASSIC_MODEL);
        if (nerrs) goto err_out;
#endif
    }
    nerrs += test_collective_error(filename, 0, NC_64BIT_DATA);
    if (nerrs) goto err_out;

    /* test in safe mode */
    setenv("PNETCDF_SAFE_MODE", "1", 1);
    nerrs += test_collective_error(filename, 1, 0);
    if (nerrs) goto err_out;
    nerrs += test_collective_error(filename, 1, NC_64BIT_OFFSET);
    if (nerrs) goto err_out;
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += test_collective_error(filename, 1, NC_NETCDF4);
        if (nerrs) goto err_out;
        nerrs += test_collective_error(filename, 1, NC_NETCDF4 | NC_CLASSIC_MODEL);
        if (nerrs) goto err_out;
#endif
    }
    nerrs += test_collective_error(filename, 1, NC_64BIT_DATA);
    if (nerrs) goto err_out;

    return (nerrs > 0);
}
#endif

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
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "collective abort", opt, test_io);

    MPI_Finalize();

    return err;
}
