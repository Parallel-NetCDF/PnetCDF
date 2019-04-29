/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program creates one big variable followed by a small variable. It
 * writes some data to both variables and reads back to check the values.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_cdf5_begin tst_cdf5_begin.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_cdf5_begin
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() strlen() */
#include <libgen.h> /* basename() */

/* This program can also be used to test NetCDF.
 * Add #define TEST_NETCDF and compile with command:
 * gcc -I/netcdf/path/include last_large_var.c -o last_large_var -L/netcdf/path/lib -lnetcdf
 */
#ifdef TEST_NETCDF
#include <netcdf.h>
#include <netcdf_meta.h>
#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,nc_strerror(err)); \
    } \
}
#define EXP_ERR(exp) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expecting %d but got %d\n", \
        __LINE__,__FILE__,exp, err); \
    } \
}
#define FileCreate(a,b,c,d,e)   nc_create(b,c,e)
#define FileOpen(a,b,c,d,e)     nc_open(b,c,e)
#define DefDim                  nc_def_dim
#define DefVar                  nc_def_var
#define SetFill                 nc_set_fill
#define EndDef                  nc_enddef
#define PutVarShort             nc_put_var_short
#define PutVaraShort            nc_put_vara_short
#define GetVaraShort            nc_get_vara_short
#define FileClose               nc_close
#define MPI_Init(a,b)
#define MPI_Comm_rank(a,b)
#define MPI_Comm_size(a,b)
#define MPI_Finalize()
typedef size_t len_t;
#else
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>
#define FileCreate              ncmpi_create
#define FileOpen                ncmpi_open
#define DefDim                  ncmpi_def_dim
#define DefVar                  ncmpi_def_var
#define SetFill                 ncmpi_set_fill
#define EndDef                  ncmpi_enddef
#define PutVarShort             ncmpi_put_var_short_all
#define PutVaraShort            ncmpi_put_vara_short_all
#define GetVaraShort            ncmpi_get_vara_short_all
#define FileClose               ncmpi_close
typedef MPI_Offset len_t;
#endif

/* When using NetCDF 4.4.1 ad prior to create a CDF-5 file and defining a small
 * variable after a big variable (> 2^32-3 bytes), the file starting offset of
 * the small variable (and all variables defined after the big variable) is
 * calculated incorrectly. This test program detects this bug by checking the
 * contents of the possible overlaps between the two variables.
 */

int main(int argc, char** argv) {
    char filename[256];
    int i, err, rank, nprocs, nerrs=0, ncid, dimid[2], varid[2];
    short buf[10];
    len_t start[1], count[1];

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
        sprintf(cmd_str, "*** TESTING C   %s for checking CDF-5 writes", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = FileCreate(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_DATA,
                     MPI_INFO_NULL, &ncid); CHECK_ERR
    err = DefDim(ncid, "dim0", NC_MAX_UINT, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "dim1", 10,          &dimid[1]); CHECK_ERR

    /* define one small variable after one big variable */
    err = DefVar(ncid, "var_big",   NC_SHORT, 1, &dimid[0], &varid[0]); CHECK_ERR
    err = DefVar(ncid, "var_small", NC_SHORT, 1, &dimid[1], &varid[1]); CHECK_ERR
    err = SetFill(ncid, NC_NOFILL, NULL); CHECK_ERR
    err = EndDef(ncid); CHECK_ERR

    /* write to var_big in location overlapping with var_small when using
     * netCDF 4.4.x or prior */
    start[0] = NC_MAX_UINT/sizeof(short);
    count[0] = 10;
    for (i=0; i<10; i++) buf[i] = i;
    err = PutVaraShort(ncid, varid[0], start, count, buf); CHECK_ERR

    /* write var_small */
    for (i=0; i<10; i++) buf[i] = -1;
    err = PutVarShort(ncid, varid[1], buf); CHECK_ERR

    /* read back var_big and check contents */
    for (i=0; i<10; i++) buf[i] = -1;
    err = GetVaraShort(ncid, varid[0], start, count,buf); CHECK_ERR
    for (i=0; i<10; i++) {
        if (buf[i] != i) {
            printf("Error at buf[%d] expect %d but got %hd\n",i,i,buf[i]);
            nerrs++;
        }
    }
    err = FileClose(ncid); CHECK_ERR

    /* check if open to read header fine */
    err = FileOpen(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = FileClose(ncid); CHECK_ERR

#ifdef TEST_NETCDF
    if (nerrs) printf("fail with %d mismatches\n",nerrs);
    else       printf("pass\n");
#else
    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR && malloc_size > 0) /* this test is for running 1 process */
        printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
               malloc_size);

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }
#endif

    MPI_Finalize();
    return (nerrs > 0);
}

