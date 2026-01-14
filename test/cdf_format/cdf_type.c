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
int test_attr_types(const char *filename,
                    int         format)
{
    int i, err, rank, ncid, cmode, nerrs=0, attr=0;
    nc_type xtype[5]={NC_UBYTE, NC_USHORT, NC_UINT, NC_INT64, NC_UINT64};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cmode = NC_CLOBBER|format;

    /* create a file in CDF-1 or CDF-2 format */
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (err != NC_NOERR) return 1;

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

    err = EndDef(ncid); CHECK_ERR
    err = FileClose(ncid); CHECK_ERR

    return nerrs;
}

/*----< test_var_types() >----------------------------------------------------*/
static
int test_var_types(const char *filename,
                   int         format)
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

/*----< test_io() >----------------------------------------------------------*/
static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,  /* ignored */
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int nerrs;

    nerrs = test_attr_types(out_path, 0); /* CDF-1 */
    if (nerrs > 0) return nerrs;

    nerrs = test_attr_types(out_path, NC_64BIT_OFFSET); /* CDF-2 */
    if (nerrs > 0) return nerrs;

    nerrs = test_var_types(out_path, 0); /* CDF-1 */
    if (nerrs > 0) return nerrs;

    nerrs = test_var_types(out_path, NC_64BIT_OFFSET); /* CDF-2 */
    if (nerrs > 0) return nerrs;

    return 0;
}

#ifndef TEST_NETCDF
/*----< main() >--------------------------------------------------------------*/
int main(int argc, char **argv) {

    int err, formats[] = {0};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "CDF-5 dtype in CDF-1 and 2", opt, test_io);

    MPI_Finalize();

    return err;
}
#else
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

    if (nerrs) printf("fail with %d mismatches\n",nerrs);
    else       printf("pass\n");

    MPI_Finalize();
    return (nerrs > 0);
}
#endif
