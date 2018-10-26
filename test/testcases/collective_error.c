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
int test_collective_error(char *filename, int safe_mode, int cmode)
{
    int rank, nproc, ncid, err, nerrs=0, varid, dimids[1], req, status;
    double buf[2];
    MPI_Offset start[1], count[1];
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    /* Create a 2 element vector of doubles */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim", 2, &dimids[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_DOUBLE, 1, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (rank == 0) {
        start[0] = 0;
        count[0] = 2;
    } else if (rank == 1) {
#if defined(PNETCDF_RELAX_COORD_BOUND) && PNETCDF_RELAX_COORD_BOUND==1
        start[0] = 3; /* illegal for a start > defined shape */
#else
        start[0] = 2; /* illegal for a start >= defined shape */
#endif
        count[0] = 0;
    }
    else {
        start[0] = 0;
        count[0] = 0;
    }

    buf[0] = 1.0;
    buf[1] = 2.0;

    err = ncmpi_put_vara_all(ncid, varid, start, count,
			     buf, count[0], MPI_DOUBLE);
    if ((safe_mode && nproc > 1) || rank == 1) EXP_ERR(NC_EINVALCOORDS)
    else                                       EXP_ERR(NC_NOERR)

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
    if ((safe_mode && nproc > 1) || rank == 1) EXP_ERR(NC_EINVALCOORDS)
    else                                       EXP_ERR(NC_NOERR)

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

    if (!(cmode & NC_NETCDF4)) {
        err = ncmpi_iput_vara_double(ncid, varid, start, count, buf, &req);
        if (rank == 1)
            EXP_ERR(NC_EINVALCOORDS)
        else
            EXP_ERR(NC_NOERR)

        err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERR

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

    err = ncmpi_get_vara_all(ncid, varid, start, count,
			     buf, count[0], MPI_DOUBLE);
    if ((safe_mode && nproc > 1) || rank == 1) EXP_ERR(NC_EINVALCOORDS)
    else                                       EXP_ERR(NC_NOERR)

    err = ncmpi_get_vara_double_all(ncid, varid, start, count, buf);
    if ((safe_mode && nproc > 1) || rank == 1) EXP_ERR(NC_EINVALCOORDS)
    else                                       EXP_ERR(NC_NOERR)

    if (!(cmode & NC_NETCDF4)) {
        err = ncmpi_iget_vara_double(ncid, varid, start, count, buf, &req);
        if (rank == 1)
            EXP_ERR(NC_EINVALCOORDS)
        else
            EXP_ERR(NC_NOERR)

        err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERR
    }

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char *argv[])
{
    char filename[256], *hint_value;
    int rank, err, nerrs=0, bb_enabled=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for collective abort ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* check whether burst buffering is enabled */
    if (inq_env_hint("nc_burst_buf", &hint_value)) {
        if (strcasecmp(hint_value, "enable") == 0) bb_enabled = 1;
        free(hint_value);
    }

    /* test in non-safe mode */
    setenv("PNETCDF_SAFE_MODE", "0", 1);
    nerrs += test_collective_error(filename, 0, 0);
    nerrs += test_collective_error(filename, 0, NC_64BIT_OFFSET);
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += test_collective_error(filename, 0, NC_NETCDF4);
        nerrs += test_collective_error(filename, 0, NC_NETCDF4 | NC_CLASSIC_MODEL);
#endif
    }
    nerrs += test_collective_error(filename, 0, NC_64BIT_DATA);

    /* test in safe mode */
    setenv("PNETCDF_SAFE_MODE", "1", 1);
    nerrs += test_collective_error(filename, 1, 0);
    nerrs += test_collective_error(filename, 1, NC_64BIT_OFFSET);
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += test_collective_error(filename, 1, NC_NETCDF4);
        nerrs += test_collective_error(filename, 1, NC_NETCDF4 | NC_CLASSIC_MODEL);
#endif
    }
    nerrs += test_collective_error(filename, 1, NC_64BIT_DATA);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}
