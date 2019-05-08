/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program is to test
 *
 * large number of dimensions per variable
 * large number of dimensions per file
 * large number of attributes per file
 * large number of variables per file
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

#define LARGE_NUM 10240

int main(int argc, char** argv)
{
    char filename[256], str[32];
    int i, rank, nprocs, err, nerrs=0;
    int ncid, cmode, *varid, *dimids, intBuf[1];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
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
        sprintf(cmd_str, "*** TESTING C   %s for large DIMS, VARS, ATTRS ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    dimids = (int*) malloc(LARGE_NUM * sizeof(int));
    varid = (int*) malloc(LARGE_NUM * sizeof(int));

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    for (i=0; i<LARGE_NUM; i++) {
        sprintf(str, "dim%d", i);
        err = ncmpi_def_dim(ncid, str, 1, &dimids[i]);
        CHECK_ERR
    }

    err = ncmpi_def_var(ncid, "var", NC_INT, LARGE_NUM, dimids, &varid[0]);
    CHECK_ERR

    for (i=0; i<LARGE_NUM; i++) {
        sprintf(str, "attr%d", i);
        err = ncmpi_put_att(ncid, varid[0], str, NC_INT, 1, &i);
        CHECK_ERR
    }

    for (i=1; i<LARGE_NUM; i++) {
        signed char attrBuf[3]={1,2,3};
        sprintf(str, "var%d", i);
        err = ncmpi_def_var(ncid, str, NC_INT, 1, dimids, &varid[i]);
        CHECK_ERR
        err = ncmpi_put_att_text(ncid, varid[i], "attr text", 9, "some text");
        CHECK_ERR
        err = ncmpi_put_att_schar(ncid, varid[i], "attr short", NC_SHORT, 3, attrBuf);
        CHECK_ERR
    }

    err = ncmpi_enddef(ncid);
    CHECK_ERR

    intBuf[0] = rank;
    for (i=0; i<LARGE_NUM; i++) {
        err = ncmpi_put_var_int_all(ncid, varid[i], intBuf);
        CHECK_ERR
    }

    err = ncmpi_close(ncid);
    CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_ERR

    for (i=0; i<LARGE_NUM; i++) {
        MPI_Offset len;
        err = ncmpi_inq_dim(ncid, i, str, &len);
        CHECK_ERR
    }

    for (i=0; i<LARGE_NUM; i++) {
        int buf;
        sprintf(str, "attr%d", i);
        err = ncmpi_get_att(ncid, 0, str, &buf);
        CHECK_ERR
    }

    for (i=0; i<LARGE_NUM; i++) {
        nc_type xtype;
        int ndims, natts;
        err = ncmpi_inq_var(ncid, i, str, &xtype, &ndims, dimids, &natts);
        CHECK_ERR
    }

    err = ncmpi_close(ncid); CHECK_ERR

    free(varid);
    free(dimids);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0) {
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
            ncmpi_inq_malloc_list();
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

