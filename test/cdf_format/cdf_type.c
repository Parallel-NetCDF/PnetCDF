/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* This program tests if PnetCDF can detect CDF-5 data types that are not
 * available in CDF-1 and CDF-2 formats
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#define DefDim                  nc_def_dim
#define DefVar                  nc_def_var
#define EndDef                  nc_enddef
#define PutAttInt               nc_put_att_int
#define FileClose               nc_close
#define NC_ERR_CODE_NAME        nc_strerrno
#define MPI_Init(a,b)
#define MPI_Comm_rank(a,b)
#define MPI_Comm_size(a,b)
#define MPI_Finalize()
#else
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>
#define FileCreate              ncmpi_create
#define DefDim                  ncmpi_def_dim
#define DefVar                  ncmpi_def_var
#define EndDef                  ncmpi_enddef
#define PutAttInt               ncmpi_put_att_int
#define FileClose               ncmpi_close
#define NC_ERR_CODE_NAME        ncmpi_strerrno
#endif

/*----< test_attr_types() >---------------------------------------------------*/
static
int test_attr_types(char *filename,
                    int   format)
{
    int i, err, rank, ncid, cmode, nerrs=0, attr=0;
    nc_type xtype[5]={NC_UBYTE, NC_USHORT, NC_UINT, NC_INT64, NC_UINT64};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cmode = NC_CLOBBER|format;

    /* create a file in CDF-1 or CDF-2 format */
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR
    for (i=0; i<5; i++) {
        char name[32];
        sprintf(name, "gattr_%d", i);
        err = PutAttInt(ncid, NC_GLOBAL, name, xtype[i], 1, &attr);
#ifdef TEST_NETCDF
        if (err != NC_EBADTYPE) {
            printf("Error at line %d in %s: expect NC_EBADTYPE but got %d\n",
                   __LINE__, __FILE__, err);
            nerrs++;
        }
#else
        if (err != NC_ESTRICTCDF2) {
            printf("Error at line %d in %s: expect NC_ESTRICTCDF2 but got %s\n",
                   __LINE__, __FILE__,NC_ERR_CODE_NAME(err));
            nerrs++;
        }
#endif
    }

    err = FileClose(ncid); CHECK_ERR

    return nerrs;
}

/*----< test_var_types() >----------------------------------------------------*/
static
int test_var_types(char *filename,
                   int   format)
{
    int i, err, rank, ncid, cmode, nerrs=0;
    int dimid, varid[5];
    nc_type xtype[5]={NC_UBYTE, NC_USHORT, NC_UINT, NC_INT64, NC_UINT64};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cmode = NC_CLOBBER|format;

    /* create a file in CDF-1 or CDF-2 format */
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = DefDim(ncid, "dim", NC_UNLIMITED, &dimid); CHECK_ERR
    for (i=0; i<5; i++) {
        char name[32];
        sprintf(name, "var_%d", i);
        err = DefVar(ncid, name, xtype[i], 1, &dimid, &varid[i]);
#ifdef TEST_NETCDF
        if (err != NC_EBADTYPE) {
            printf("Error at line %d in %s: expect NC_EBADTYPE but got %d\n",
                   __LINE__, __FILE__, err);
            nerrs++;
        }
#else
        if (err != NC_ESTRICTCDF2) {
            printf("Error at line %d in %s: expect NC_ESTRICTCDF2 but got %s\n",
                   __LINE__, __FILE__,NC_ERR_CODE_NAME(err));
            nerrs++;
        }
#endif
    }
    err = FileClose(ncid); CHECK_ERR

    return nerrs;
}

/*----< main() >--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char filename[256];
    int rank, nerrs=0;

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
        sprintf(cmd_str, "*** TESTING C   %s for CDF-5 type in CDF-1 and 2 ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    nerrs += test_attr_types(filename, 0);
    nerrs += test_attr_types(filename, NC_64BIT_OFFSET);

    nerrs += test_var_types(filename, 0);
    nerrs += test_var_types(filename, NC_64BIT_OFFSET);

#ifdef TEST_NETCDF
    if (nerrs) printf("fail with %d mismatches\n",nerrs);
    else       printf("pass\n");
#else
    MPI_Offset malloc_size, sum_size;
    int err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }
#endif

    MPI_Finalize();
    return (nerrs > 0);
}

