/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* This test program checks under safe mode whether a collective API can be
 * nicely aborted without causing the program to hang. It runs on 2 processes.
 * One process deliberately produces an error (using an illegal start argument),
 * while the other does not.
 * Note when not in safe mode, the program will hang.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define ERR { if (err!=NC_NOERR){printf("PE %d: error at line %d (%s)\n",rank,__LINE__,ncmpi_strerror(err)); nerrs++;}}

#define EXP_ERR(e) { \
    if (err!=e) { \
        printf("PE %d: error at line %d expecting %s but got %s\n",rank,__LINE__,nc_err_code_name(e),nc_err_code_name(err)); \
        nerrs++; \
    } \
}

static
int test_collective_error(char *filename, int safe_mode)
{
    int rank, nproc, ncid, err, nerrs=0, varid, dimids[1], req, status;
    double buf[2];
    MPI_Offset start[1], count[1];
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    /* Create a 2 element vector of doubles */
    err = ncmpi_create(comm, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim", 2, &dimids[0]); ERR
    err = ncmpi_def_var(ncid, "var", NC_DOUBLE, 1, dimids, &varid); ERR
    err = ncmpi_enddef(ncid); ERR

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
    if (nproc > 1) EXP_ERR(NC_EINVALCOORDS)
    else           EXP_ERR(NC_NOERR)

    /* check if user put buffer contents altered */
    if (buf[0] != 1.0) {
        printf("Error: user put buffer[%d] altered from %f to %f\n",
               0, 1.0, buf[0]);
        nerrs++;
    }
    if (buf[1] != 2.0) {
        printf("Error: user put buffer[%d] altered from %f to %f\n",
               1, 2.0, buf[1]);
        nerrs++;
    }

    err = ncmpi_put_vara_double_all(ncid, varid, start, count, buf);
    if (nproc > 1) EXP_ERR(NC_EINVALCOORDS)
    else           EXP_ERR(NC_NOERR)

    /* check if user put buffer contents altered */
    if (buf[0] != 1.0) {
        printf("Error: user put buffer[%d] altered from %f to %f\n",
               0, 1.0, buf[0]);
        nerrs++;
    }
    if (buf[1] != 2.0) {
        printf("Error: user put buffer[%d] altered from %f to %f\n",
               1, 2.0, buf[1]);
        nerrs++;
    }

    err = ncmpi_iput_vara_double(ncid, varid, start, count, buf, &req);
    if (rank == 1)
        EXP_ERR(NC_EINVALCOORDS)
    else
        EXP_ERR(NC_NOERR)

    err = ncmpi_wait_all(ncid, 1, &req, &status); ERR

    /* check if user put buffer contents altered */
    if (buf[0] != 1.0) {
        printf("Error: user put buffer[%d] altered from %f to %f\n",
               0, 1.0, buf[0]);
        nerrs++;
    }
    if (buf[1] != 2.0) {
        printf("Error: user put buffer[%d] altered from %f to %f\n",
               1, 2.0, buf[1]);
        nerrs++;
    }

    err = ncmpi_get_vara_all(ncid, varid, start, count,
			     buf, count[0], MPI_DOUBLE);
    if (nproc > 1) EXP_ERR(NC_EINVALCOORDS)
    else           EXP_ERR(NC_NOERR)

    err = ncmpi_get_vara_double_all(ncid, varid, start, count, buf);
    if (nproc > 1) EXP_ERR(NC_EINVALCOORDS)
    else           EXP_ERR(NC_NOERR)

    err = ncmpi_iget_vara_double(ncid, varid, start, count, buf, &req);
    if (rank == 1)
        EXP_ERR(NC_EINVALCOORDS)
    else
        EXP_ERR(NC_NOERR)

    err = ncmpi_wait_all(ncid, 1, &req, &status); ERR

    err = ncmpi_close(ncid); ERR

    return nerrs;
}

int main(int argc, char *argv[])
{
    char *filename="testfile.nc";
    int rank, err, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    if (argc == 2) filename = argv[1];
    assert(filename != NULL);

    if (rank == 0) {
        char cmd_str[256];
        sprintf(cmd_str, "*** TESTING C   %s for collective abort ", argv[0]);
        printf("%-66s ------ ", cmd_str); fflush(stdout);
    }

    /* Must test in safe mode, otherwise the program will hang if nproc > 1 */
    setenv("PNETCDF_SAFE_MODE", "1", 1);
    nerrs += test_collective_error(filename, 1);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs == 0) ? 0 : 1;
}
