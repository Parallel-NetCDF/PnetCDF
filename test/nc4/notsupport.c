/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests fill mode for individual variables.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_def_var_fill tst_def_var_fill.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 tst_def_var_fill testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define EXPECT_ERR(err_no1) \
    if (err != err_no1) { \
        nerrs++; \
        printf("Error at line %d in %s: expect error code %s but got %s\n", \
               __LINE__,__FILE__,ncmpi_strerrno(err_no1),ncmpi_strerrno(err)); \
    }

#define NY 8
#define NX 5

int main(int argc, char** argv) {
    char filename[256];
    int rank, nprocs, err, nerrs=0;
    int ncid, cmode, varid, dimid, buf;
    int v1, v2;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Offset start[1] = {0};
    MPI_Datatype btype;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

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
        sprintf(cmd_str, "*** TESTING C   %s for error NC_ENOTSUPPORT ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* create a new file for writing ------------------------------------*/
    cmode = NC_CLOBBER | NC_NETCDF4;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR

    /* define dimension */
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid); CHECK_ERR

    /* define variable */
    err = ncmpi_def_var(ncid, "M",   NC_INT, 1, &dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_striping(ncid, &v1, &v2); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_inq_put_size(ncid, start); CHECK_ERR

    err = ncmpi_inq_get_size(ncid, start); CHECK_ERR

    err = ncmpi_inq_header_size(ncid, start); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_inq_header_extent(ncid, start); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_inq_nreqs(ncid, &v1); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_inq_buffer_usage(ncid, start); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_inq_buffer_size(ncid, start); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_sync_numrecs(ncid); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_flush(ncid); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_fill_var_rec(ncid, varid, 0); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_inq_varoffset(ncid, varid, start); EXPECT_ERR(NC_ENOTSUPPORT);

    err = MPI_Type_vector(1, 2, 1, MPI_INT, &btype);
    if (err != MPI_SUCCESS ){
        nerrs++;
        printf("Error at line %d in %s: %d\n", __LINE__, __FILE__, err);
    }
    err = MPI_Type_commit(&btype);
    if (err != MPI_SUCCESS ){
        nerrs++;
        printf("Error at line %d in %s: %d\n", __LINE__, __FILE__, err);
    }
    else {
        err = ncmpi_put_var1_all(ncid, varid, start, &buf, 1, btype); EXPECT_ERR(NC_ENOTSUPPORT);

        err = ncmpi_put_vard_all(ncid, varid, btype, &buf, 1, btype); EXPECT_ERR(NC_ENOTSUPPORT);
        err = MPI_Type_free(&btype);
        if (err != MPI_SUCCESS ){
            nerrs++;
            printf("Error at line %d in %s: %d\n", __LINE__, __FILE__, err);
        }
    }

    err = ncmpi_iput_var1_int(ncid, varid, start, &buf, NULL); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_cancel(ncid, NC_REQ_ALL, NULL, NULL); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_buffer_attach(ncid, 1); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_buffer_detach(ncid); EXPECT_ERR(NC_ENOTSUPPORT);

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, comm);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

