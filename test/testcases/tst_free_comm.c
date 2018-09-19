/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests whether PnetCDF duplicates MPI communicator and MPI info
 * object correctly, so the user supplied communicator and info objects can be
 * freed by users after ncmpi_create and ncmpi_open.
 *
 *  % mpiexec -n 4 tst_free_comm
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h> /* strcpy(), strncpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

static int
tst_fmt(char *fname, int cmode)
{
    int nerrs=0, err, exp_err=NC_NOERR, ncid;
    MPI_Comm comm=MPI_COMM_NULL;
    MPI_Info info=MPI_INFO_NULL;

#ifndef ENABLE_NETCDF4
    if (cmode & NC_NETCDF4)
        exp_err = NC_ENOTBUILT;
#endif

    /* duplicate MPI_COMM_WORLD */
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    /* create MPI I/O hints */
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_no_indep_rw", "true");

    /* create a file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, fname, cmode, info, &ncid); EXP_ERR(exp_err)
    if (err == NC_ENOTBUILT) goto fn_exit;

    MPI_Comm_free(&comm);
    MPI_Info_free(&info);

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the file */

    /* duplicate MPI_COMM_WORLD */
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    /* create MPI I/O hints */
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_no_indep_rw", "true");

    err = ncmpi_open(comm, fname, NC_NOWRITE, info, &ncid); CHECK_ERR

    MPI_Comm_free(&comm); comm = MPI_COMM_NULL;
    MPI_Info_free(&info); info = MPI_INFO_NULL;

    err = ncmpi_close(ncid); CHECK_ERR

fn_exit:
    if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);
    return nerrs;
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {

    char filename[256];
    int  err, nerrs=0, rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for freeing MPI communicator ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    nerrs += tst_fmt(filename, 0);
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
    nerrs += tst_fmt(filename, NC_NETCDF4);
    nerrs += tst_fmt(filename, NC_NETCDF4|NC_CLASSIC_MODEL);
    nerrs += tst_fmt(filename, NC_64BIT_DATA);

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

