/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 *
 *  Test flexible API ncmpi_put_vara()
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define VECCOUNT 4
#define BLOCKLEN 3
#define STRIDE   5

int main(int argc, char ** argv)
{
    int ncid, dimid, varid, rank, nprocs;
    MPI_Datatype vtype, rtype, usertype;
    MPI_Aint lb, extent;
    int userbufsz, *userbuf, *cmpbuf, i, err, errs=0, nerrs=0;
    int count = 25;
    double pi = 3.14159;
    MPI_Offset start, acount;
    char filename[256];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for put_vara/get_vara ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

#ifdef DEBUG
    if (nprocs > 2 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n",argv[0]);
#endif

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, "50k", 1024*50, &dimid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "vector", NC_DOUBLE, 1, &dimid, &varid);
    CHECK_ERR

    err = ncmpi_enddef(ncid);
    CHECK_ERR

    MPI_Type_vector(VECCOUNT, BLOCKLEN, STRIDE, MPI_INT, &vtype);
    MPI_Type_create_resized(vtype, 0, STRIDE*VECCOUNT*sizeof(int), &rtype);
    MPI_Type_contiguous(count, rtype, &usertype);
    MPI_Type_commit(&usertype);

    MPI_Type_free(&vtype);
    MPI_Type_free(&rtype);

    MPI_Type_get_extent(usertype, &lb, &extent);
    userbufsz = extent;
    userbuf = (int*) malloc(userbufsz);
    cmpbuf = (int*) calloc(userbufsz, 1);
    for (i=0; i< userbufsz/sizeof(int); i++)
        userbuf[i] = pi*i;

    start = 10; acount = count*12;
    err = ncmpi_begin_indep_data(ncid);
    CHECK_ERR
    if (rank == 0) {
        err = ncmpi_put_vara(ncid, varid, &start, &acount, userbuf, 1, usertype);
        CHECK_ERR
    }

    err = ncmpi_close(ncid);
    CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR
    err = ncmpi_begin_indep_data(ncid);
    CHECK_ERR
    err = ncmpi_inq_varid(ncid, "vector", &varid);
    CHECK_ERR
    err = ncmpi_get_vara(ncid, varid, &start, &acount, cmpbuf, 1, usertype);
    CHECK_ERR
    err = ncmpi_close(ncid);
    CHECK_ERR

    for (i=0; errs < 10 &&  i < acount; i++) {
        /* vector of 4,3,5, so skip 4th and 5th items of every block */
        if (i%STRIDE >= BLOCKLEN) continue;
        if (userbuf[i] != cmpbuf[i]) {
            errs++;
            fprintf(stderr, "%d: expected 0x%x got 0x%x\n",
                    i, userbuf[i], cmpbuf[i]);
        }
    }
    nerrs += errs;
    free(userbuf);
    free(cmpbuf);
    MPI_Type_free(&usertype);

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
