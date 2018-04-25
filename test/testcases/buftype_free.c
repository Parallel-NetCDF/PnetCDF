/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This example tests if PnetCDF duplicates the MPI derived data type supplied
 * by the user, when calling the flexible APIs. It tests a PnetCDF bug
 * prior to version 1.6.1.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 4
#define NX 4

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {

    char filename[256];
    int  i, j, err, ncid, varid[4], dimids[2], req[4], st[4], nerrs=0;
    int  rank, nprocs, buf[4][(NY+4)*(NX+4)];
    int  gsize[2], subsize[2], a_start[2], ghost;
    MPI_Offset start[2], count[2];
    MPI_Datatype buftype[4];

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
        sprintf(cmd_str, "*** TESTING C   %s for free buftype in flexible API ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "Y", NY*nprocs, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX,        &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 2, dimids, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 2, dimids, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var2", NC_INT, 2, dimids, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var3", NC_INT, 2, dimids, &varid[3]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* initialize the contents of the array */
    for (i=0; i<4; i++) for (j=0; j<(NY+4)*(NX+4); j++) buf[i][j] = rank+10;

    start[0] = NY*rank; start[1] = 0;
    count[0] = NY;      count[1] = NX;

    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf[0]); CHECK_ERR
    err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf[1]); CHECK_ERR
    err = ncmpi_put_vara_int_all(ncid, varid[2], start, count, buf[2]); CHECK_ERR
    err = ncmpi_put_vara_int_all(ncid, varid[3], start, count, buf[3]); CHECK_ERR

    /* check if user write buffer contents altered */
    for (i=0; i<4; i++) {
        for (j=0; j<(NY+4)*(NX+4); j++) {
            if (buf[i][j] != rank+10) {
                printf("Error at line %d in %s: user put buffer[%d][%d] altered from %d to %d\n",
                       __LINE__,__FILE__,i,j, rank+10, buf[i][j]);
                nerrs++;
            }
        }
    }

    /* define an MPI datatype using MPI_Type_create_subarray() */
    ghost      = 2;
    gsize[1]   = NX + 2 * ghost;
    gsize[0]   = NY + 2 * ghost;
    subsize[1] = NX;
    subsize[0] = NY;
    a_start[1] = ghost;
    a_start[0] = ghost;

    for (i=0; i<4; i++) {
        req[i] = NC_REQ_NULL;
        st[i]  = NC_NOERR;
        MPI_Type_create_subarray(2, gsize, subsize, a_start, MPI_ORDER_C, MPI_INT, &buftype[i]);
        MPI_Type_commit(&buftype[i]);

        err = ncmpi_iget_vara(ncid, varid[i], start, count, buf[i], 1, buftype[i], &req[i]); CHECK_ERR
        MPI_Type_free(&buftype[i]);
    }

    err = ncmpi_wait_all(ncid, 4, req, st); CHECK_ERR
    for (i=0; i<4; i++) {
        if (st[i] != NC_NOERR) {
            err = st[i];
            CHECK_ERR
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

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

