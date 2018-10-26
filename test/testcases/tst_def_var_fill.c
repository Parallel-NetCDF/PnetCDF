/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
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
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 8
#define NX 5

static int
tst_fmt(char *filename, int cmode)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, format, varid[2], dimid[2], expect, *buf;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* allocate I/O buffer */
    buf = (int*) malloc(NY*NX * sizeof(int));

    /* create a new file for writing ------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* define dimension */
    err = ncmpi_def_dim(ncid, "Y", NY,        &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs, &dimid[1]); CHECK_ERR

    /* define 2 variables of size NY x NX */
    err = ncmpi_def_var(ncid, "var_nofill", NC_INT, 2, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_fill",   NC_INT, 2, dimid, &varid[1]); CHECK_ERR
    /* set fill modes for the variables */
    err = ncmpi_def_var_fill(ncid, NC_GLOBAL,1, NULL); EXP_ERR(NC_EGLOBAL)
    err = ncmpi_def_var_fill(ncid, varid[0], 1, NULL); CHECK_ERR
    err = ncmpi_def_var_fill(ncid, varid[1], 0, NULL); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* initialize I/O buffer */
    for (i=0; i<NY*NX; i++) buf[i] = rank+5;

    /* write a subarray to each variable */
    start[0] = 0;
    start[1] = NX*rank+2;
    count[0] = NY;
    count[1] = 2;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR
    /* check if user put buffer contents altered */
    for (i=0; i<NY*NX; i++) {
        if (buf[i] != rank+5) {
            printf("Error in %s line %d: put buf[%d] altered from %d to %d\n",
                   __FILE__,__LINE__, i, rank+5, buf[i]);
            nerrs++;
        }
    }
    err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR
    /* check if user put buffer contents altered */
    for (i=0; i<NY*NX; i++) {
        if (buf[i] != rank+5) {
            printf("Error in %s line %d: put buf[%d] altered from %d to %d\n",
                   __FILE__,__LINE__, i, rank+5, buf[i]);
            nerrs++;
        }
    }
    err = ncmpi_close(ncid); CHECK_ERR

    /* reopen the file and read data back */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_format(ncid, &format); CHECK_ERR
    if (cmode == 0) {
        if (format != NC_FORMAT_CLASSIC) {
            printf("Error at line %d of %s: expect NC_FORMAT_CLASSIC but got %d\n",
                   __LINE__,__FILE__,format);
            nerrs++;
        }
    } else if (cmode == NC_64BIT_OFFSET) {
        if (format != NC_FORMAT_64BIT_OFFSET) {
            printf("Error at line %d of %s: expect NC_FORMAT_64BIT_OFFSET but got %d\n",
                   __LINE__,__FILE__,format);
            nerrs++;
        }
    } else if (cmode == NC_64BIT_DATA) {
        if (format != NC_FORMAT_64BIT_DATA) {
            printf("Error at line %d of %s: expect NC_FORMAT_64BIT_DATA but got %d\n",
                   __LINE__,__FILE__,format);
            nerrs++;
        }
    } else if (cmode == NC_NETCDF4) {
        if (format != NC_FORMAT_NETCDF4) {
            printf("Error at line %d of %s: expect NC_FORMAT_NETCDF4 but got %d\n",
                   __LINE__,__FILE__,format);
            nerrs++;
        }
    } else if (cmode == (NC_NETCDF4|NC_CLASSIC_MODEL)) {
        if (format != NC_FORMAT_NETCDF4_CLASSIC) {
            printf("Error at line %d of %s: expect NC_FORMAT_NETCDF4_CLASSIC but got %d\n",
                   __LINE__,__FILE__,format);
            nerrs++;
        }
    }

    /* inquire variabe IDs */
    err = ncmpi_inq_varid(ncid, "var_nofill", &varid[0]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "var_fill",   &varid[1]); CHECK_ERR

    /* read the subarray written by process (rank+1)%nproc */
    start[0] = 0;
    start[1] = NX*((rank+1)%nprocs);
    count[0] = NY;
    count[1] = NX;
    for (i=0; i<NY*NX; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR

    /* check contents of variable var_nofill */
    expect = (rank+1)%nprocs + 5;
    for (i=0; i<NY; i++) for (j=0; j<NX; j++) {
        if (2 <= j && j < 4) {
            if (buf[i*NX+j] != expect) {
                printf("Error in %s line %d: expect get buf[%d]=%d but got %d\n",
                       __FILE__,__LINE__,i*NX+j, expect, buf[i*NX+j]);
                nerrs++;
            }
        }
        else if (buf[i*NX+j] == NC_FILL_INT) /* not an error */
            /* content of buf[i*NX+j] can be any value */
            printf("Warning in %s line %d: get buf[%d] same as NC_FILL_INT\n",
                   __FILE__,__LINE__,i*NX+j);
    }

    /* read the subarray written by process (rank+1)%nproc */
    for (i=0; i<NY*NX; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR

    /* check contents of variable var_fill */
    for (i=0; i<NY; i++) for (j=0; j<NX; j++) {
        expect = NC_FILL_INT;
        if (2 <= j && j< 4) expect = (rank+1)%nprocs + 5;
        if (buf[i*NX+j] != expect) {
            printf("Error in %s line %d: expect get buf[%d]=%d but got %d\n",
                   __FILE__,__LINE__,i*NX+j, expect, buf[i*NX+j]);
            nerrs++;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);

    return nerrs;
}

int main(int argc, char** argv) {
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
        sprintf(cmd_str, "*** TESTING C   %s for def_var_fill ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* check whether burst buffering is enabled */
    if (inq_env_hint("nc_burst_buf", &hint_value)) {
        if (strcasecmp(hint_value, "enable") == 0) bb_enabled = 1;
        free(hint_value);
    }

    nerrs += tst_fmt(filename, 0);
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += tst_fmt(filename, NC_NETCDF4);
        nerrs += tst_fmt(filename, NC_NETCDF4 | NC_CLASSIC_MODEL);
#endif
    }
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

