/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
<<<<<<< HEAD
<<<<<<< HEAD
 * This program check whether BP file header information are synbced accross 
 * processes properly
=======
 * This program adds two new variables to an existing netCDF file.
 * It is used to test if PnetCDF can correctly calculate the file offsets
 * for the two new variables, in particular for files that align the
 * fixed-size variables to a boundary larger than 4 bytes, for instance
 * a file created by PnetCDF with defaut alignment of 512 bytes.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o add_var add_var.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 add_var testfile.nc
>>>>>>> 272f00a... header sync test
=======
 * This program check whether BP file header information are synbced accross 
 * processes properly
>>>>>>> 23ad299... parse var attributes
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

int cmp_int(int *in, int n, char *lbl){
    int i, j;
    int rank, np;
<<<<<<< HEAD
<<<<<<< HEAD
    int nerrs = 0;
=======
    int err, nerrs = 0;
>>>>>>> 272f00a... header sync test
=======
    int nerrs = 0;
>>>>>>> 3180e74... remove unused variables on test programs
    int *all;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    all = malloc(n * sizeof(int));

    MPI_Gather(&n, 1, MPI_INT, all, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for(i = 1; i < np; i++){
            if (all[i] != all[0]){
<<<<<<< HEAD
<<<<<<< HEAD
                printf("Number of element at rank %d = %d, expect %d\n", i, 
                        all[i], all[0]);
=======
                printf("Number of element at rank %d = %d, expect %d\n", i, all[i], all[0]);
>>>>>>> 272f00a... header sync test
=======
                printf("Number of element at rank %d = %d, expect %d\n", i, 
                        all[i], all[0]);
>>>>>>> ee05fef... wrap to 80 columns
                nerrs++;
            }
        }
    }

    for(j = 0; j < n; j++){
        MPI_Gather(in + j, 1, MPI_INT, all, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            for(i = 1; i < np; i++){
                if (all[i] != all[0]){
<<<<<<< HEAD
<<<<<<< HEAD
                    printf("%s[%d] at rank %d = %d, expect %d\n", lbl, j, i, 
                            all[i], all[0]);
=======
                    printf("%s[%d] at rank %d = %d, expect %d\n", lbl, j, i, all[i], all[0]);
>>>>>>> 272f00a... header sync test
=======
                    printf("%s[%d] at rank %d = %d, expect %d\n", lbl, j, i, 
                            all[i], all[0]);
>>>>>>> ee05fef... wrap to 80 columns
                    nerrs++;
                }
            }
        }
    }

    free(all);

    return nerrs;
}

int main(int argc, char** argv) {
<<<<<<< HEAD
<<<<<<< HEAD
    char filename[256], lbl[256];
    int rank, np, err, nerrs=0;
=======
    int i;
    char filename[256], *hint_value, lbl[256];
    int rank, np, err, nerrs=0, bb_enabled=0;
>>>>>>> 272f00a... header sync test
=======
    char filename[256], lbl[256];
    int rank, np, err, nerrs=0;
>>>>>>> 3180e74... remove unused variables on test programs
    int ncid, nvar, ndim;
    int varid, dimid, len;
    MPI_Offset dimlen;
    int dimids[1024];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "arrays.bp");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
<<<<<<< HEAD
<<<<<<< HEAD
        sprintf(cmd_str, 
                "*** TESTING C   %s for checking offsets of new variables ", 
                basename(argv[0]));
=======
        sprintf(cmd_str, "*** TESTING C   %s for checking offsets of new variables ", basename(argv[0]));
>>>>>>> 272f00a... header sync test
=======
        sprintf(cmd_str, 
                "*** TESTING C   %s for checking offsets of new variables ", 
                basename(argv[0]));
>>>>>>> ee05fef... wrap to 80 columns
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

<<<<<<< HEAD
<<<<<<< HEAD
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, 
                        &ncid); CHECK_ERR
=======
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR
>>>>>>> 272f00a... header sync test
=======
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, 
                        &ncid); CHECK_ERR
>>>>>>> ee05fef... wrap to 80 columns
    err = ncmpi_inq(ncid, &ndim, &nvar, NULL, NULL); CHECK_ERR

    nerrs += cmp_int(&ndim, 1, "ndim");
    nerrs += cmp_int(&nvar, 1, "nvar");

    for(dimid = 0; dimid < ndim; dimid++){
        err = ncmpi_inq_dimlen(ncid, dimid, &dimlen); CHECK_ERR
        len = (int)dimlen;
        sprintf(lbl, "size of dim %d", dimid);
        nerrs += cmp_int(&len, 1, lbl);
    }

    for(varid = 0; varid < nvar; varid++){
        err = ncmpi_inq_varndims(ncid, varid, &ndim); CHECK_ERR
        sprintf(lbl, "ndim of var %d", varid);
        nerrs += cmp_int(&ndim, 1, lbl);
        if (ndim < 1024){
            err = ncmpi_inq_vardimid(ncid, varid, dimids); CHECK_ERR
        }
        sprintf(lbl, "dims of var %d", varid);
        nerrs += cmp_int(dimids, ndim, lbl);
    }

    ncmpi_close(ncid);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
<<<<<<< HEAD
<<<<<<< HEAD
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, 
                    MPI_COMM_WORLD);
=======
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
>>>>>>> 272f00a... header sync test
=======
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, 
                    MPI_COMM_WORLD);
>>>>>>> ee05fef... wrap to 80 columns
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

