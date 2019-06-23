/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

#include <libgen.h> /* basename() */
#include <unistd.h>
#include "tests.h"

/*
 * Test driver for netCDF-3 interface.  This program performs tests against
 * the netCDF-3 specification for all user-level functions in an
 * implementation of the netCDF library.
 *
 * Unless invoked with "-r" (read_only) option, must be invoked from a
 * directory in which the invoker has write permission.
 *
 * Files:
 * The read-only tests read files:
 *     test.nc (see below)
 *     test_get.c (used merely as an example of a non-netCDF file)
 *
 * The write tests
 *     read test.nc (see below)
 *     write scratch.nc (deleted after each test)
 *
 * The file test.nc is created by running nc_test with the -c (create) option.
 * It is described by the following global variables.
 */


/*
 * global variables (defined by function init_gvars) describing file test.nc
 */
char dim_name[NDIMS][3];
MPI_Offset dim_len[NDIMS];
char var_name[NVARS][2+MAX_RANK];
nc_type var_type[NVARS];
int var_rank[NVARS];
int var_dimid[NVARS][MAX_RANK];
MPI_Offset var_shape[NVARS][MAX_RANK];
MPI_Offset var_nels[NVARS];
int  var_natts[NVARS];
char att_name[NVARS][MAX_NATTS][2];
char gatt_name[NGATTS][3];
nc_type att_type[NVARS][NGATTS];
nc_type gatt_type[NGATTS];
MPI_Offset att_len[NVARS][MAX_NATTS];
MPI_Offset gatt_len[NGATTS];

/*
 * command-line options
 */
int  read_only;                /* if 1, don't try to change files */
int  verbose;                /* if 1, print details of tests */
int  max_nmpt;                /* max. number of messages per test */

/*
 * Misc. global variables
 */
int  nfails;                /* number of failures in specific test */
char testfile[256];
char scratch[256];
MPI_Comm comm = MPI_COMM_WORLD; /* mpi communicator */
MPI_Info info;

static void
usage(char *progname)
{
#ifdef ENABLE_NETCDF4
    error("%s [-hrv245] [-n <MAX_NMPT>]\n", progname);
#else
    error("%s [-hrv25] [-n <MAX_NMPT>]\n", progname);
#endif
    error("   [-h] Print help\n" );
    error("   [-r] Just do read-only tests\n" );
    error("   [-v] Verbose mode\n" );
    error("   [-2] (with -c) create file with CDF-2 format\n" );
    error("   [-5] (with -c) create file with CDF-5 format\n" );
#ifdef ENABLE_NETCDF4
    error("   [-4] (with -c) create file with NetCDF-4 classic-model format\n" );
#endif
    error("   [-n <MAX_NMPT>] max. number of messages per test (Default: 8)\n");
    error("   [-d directory] directory for storing input/output files\n");
}

#define NC_CHECK_AND_PRINT                                               \
    nfailsTotal += nfails;                                               \
    if (verbose) print( "*** Testing %-30s ... ",func_name);             \
    if (nfails == 0 && verbose) {                                        \
        if (noks > 0) printf("%4d good comparisons. ok", noks);          \
        printf("\n");                                                    \
    }                                                                    \
    if (nfails > 0) {                                                    \
        print("\n\t### %d FAILURES TESTING %s! Stop ... ###\n",          \
              nfails,func_name);                                         \
        goto fn_exit;                                                    \
    }

#define NC_TEST(func) {                                                  \
    int noks;                                                            \
    char func_name[64];                                                  \
    nfails = 0;                                                          \
    sprintf(func_name, "test_%s",#func);                                 \
    noks = test_ ## func();                                              \
    NC_CHECK_AND_PRINT                                                   \
}

#define NC_TEST1(func, arg) {                                            \
    int noks;                                                            \
    char func_name[64];                                                  \
    sprintf(func_name, "test_%s",#func);                                 \
    noks = test_ ## func(arg);                                           \
    NC_CHECK_AND_PRINT                                                   \
}

#define NC_TEST2(func, arg1, arg2) {                                     \
    int noks;                                                            \
    char func_name[64];                                                  \
    sprintf(func_name, "test_%s",#func);                                 \
    noks = test_ ## func(arg1, arg2);                                    \
    NC_CHECK_AND_PRINT                                                   \
}

#if 1                /* both CRAY MPP and OSF/1 Alpha systems need this */
#include <signal.h>
#endif /* T90 */

int
main(int argc, char *argv[])
{
    extern char *optarg;
    char *cmd_str;
    int cdf_format, c;
    int numGatts, numTypes, numVars;
    int nfailsTotal = 0;        /* total number of failures */

#if 1                /* both CRAY MPP and OSF/1 Alpha systems need this */
        /*
         * Some of the extreme test assignments in this program trigger
         * floating point exceptions on CRAY T90
         */
        (void) signal(SIGFPE, SIG_IGN);
#endif

    MPI_Init(&argc, &argv);

    cdf_format = 1;         /* 1: CDF-1, 2: CDF-2 5: CDF-5 */
    read_only = 0;               /* assume may write in test dir as default */
    verbose = 0;
    max_nmpt = 8;
    strcpy(testfile, "test.nc");    /* read-only testfile */
    strcpy(scratch, "scratch.nc");  /* writable scratch file */

    while ((c = getopt(argc, argv, "245hrn:d:v")) != -1)
      switch(c) {
        case 'r':                /* just perform read-only tests */
          read_only = 1;
          break;
        case 'v':                /* verbose mode */
          verbose = 1;
          break;
        case 'n':                /* max. number of messages per test */
          max_nmpt = (int)strtol(optarg,NULL,10);
          break;
        case '2':
          cdf_format = 2;
          break;
        case '4':
          cdf_format = 4;
          break;
        case '5':
          cdf_format = 5;
          break;
        case 'd':
          sprintf(testfile, "%s/test.nc", optarg);
          sprintf(scratch, "%s/scratch.nc", optarg);
          break;
        case 'h':
        case '?':
          usage(argv[0]);
          MPI_Finalize();
          return 1;
      }

#ifndef ENABLE_NETCDF4
    if (cdf_format == 4) {
        printf("Error: NetCDF-4 support is not enabled at configure time\n");
        MPI_Finalize();
        return 1;
    }
#endif

    MPI_Info_create(&info);
    /* MPI_Info_set(info, "romio_pvfs2_posix_write", "enable"); */
    /* disable MPI-IO data sieving */
    MPI_Info_set(info, "romio_ds_write", "disable");
    MPI_Info_set(info, "romio_lustre_ds_in_coll", "disable");

    numGatts = 6;
    numVars  = 136;
    numTypes = 6;
    if (cdf_format == 5) {
        numGatts = NGATTS;
        numVars  = NVARS;
        numTypes = NTYPES;
    }

    if (cdf_format == 2)
        ncmpi_set_default_format(NC_FORMAT_CDF2, NULL);
    else if (cdf_format == 4)
        ncmpi_set_default_format(NC_FORMAT_NETCDF4_CLASSIC, NULL);
    else if (cdf_format == 5)
        ncmpi_set_default_format(NC_FORMAT_CDF5, NULL);
    else
        ncmpi_set_default_format(NC_FORMAT_CLASSIC, NULL);

    /* Initialize global variables defining test file */
    init_gvars(numGatts, numTypes, numVars);

    /* delete testfile file and ignore the error if not exist */
    unlink(testfile);

    /* create file test.nc for testing read operations */
    write_file(testfile, numGatts, numVars);
    if (nfailsTotal > 0) goto fn_exit;

    cmd_str = (char*)malloc(strlen(argv[0]) + 256);
    if (cdf_format == 4)
        sprintf(cmd_str, "*** TESTING C   %s for NetCDF4 classic-model format ", basename(argv[0]));
    else
        sprintf(cmd_str, "*** TESTING C   %s for format CDF-%d ", basename(argv[0]), cdf_format);
    printf("%-66s ------ ",cmd_str);
    free(cmd_str);

    /* Test read-only functions, using pregenerated test-file */
    NC_TEST(ncmpi_strerror);
    NC_TEST(ncmpi_open);
    NC_TEST(ncmpi_close);
    NC_TEST2(ncmpi_inq, numGatts, numVars);
    NC_TEST(ncmpi_inq_dimid);
    NC_TEST(ncmpi_inq_dim);
    NC_TEST(ncmpi_inq_dimlen);
    NC_TEST(ncmpi_inq_dimname);
    NC_TEST1(ncmpi_inq_varid, numVars);
    NC_TEST1(ncmpi_inq_var, numVars);
    NC_TEST1(ncmpi_inq_natts, numGatts);
    NC_TEST(ncmpi_inq_ndims);
    NC_TEST1(ncmpi_inq_nvars, numVars);
    NC_TEST(ncmpi_inq_unlimdim);
    NC_TEST1(ncmpi_inq_vardimid, numVars);
    NC_TEST1(ncmpi_inq_varname, numVars);
    NC_TEST2(ncmpi_inq_varnatts, numGatts, numVars);
    NC_TEST1(ncmpi_inq_varndims, numVars);
    NC_TEST1(ncmpi_inq_vartype, numVars);
    NC_TEST1(ncmpi_get_var_text, numVars);
    NC_TEST1(ncmpi_get_var_schar, numVars);
    NC_TEST1(ncmpi_get_var_short, numVars);
    NC_TEST1(ncmpi_get_var_int, numVars);
    NC_TEST1(ncmpi_get_var_long, numVars);
    NC_TEST1(ncmpi_get_var_float, numVars);
    NC_TEST1(ncmpi_get_var_double, numVars);
    NC_TEST1(ncmpi_get_var_uchar, numVars);
    NC_TEST1(ncmpi_get_var_ushort, numVars);
    NC_TEST1(ncmpi_get_var_uint, numVars);
    NC_TEST1(ncmpi_get_var_longlong, numVars);
    NC_TEST1(ncmpi_get_var_ulonglong, numVars);
    NC_TEST1(ncmpi_get_var1_text, numVars);
    NC_TEST1(ncmpi_get_var1_schar, numVars);
    NC_TEST1(ncmpi_get_var1_short, numVars);
    NC_TEST1(ncmpi_get_var1_int, numVars);
    NC_TEST1(ncmpi_get_var1_long, numVars);
    NC_TEST1(ncmpi_get_var1_float, numVars);
    NC_TEST1(ncmpi_get_var1_double, numVars);
    NC_TEST1(ncmpi_get_var1_uchar, numVars);
    NC_TEST1(ncmpi_get_var1_ushort, numVars);
    NC_TEST1(ncmpi_get_var1_uint, numVars);
    NC_TEST1(ncmpi_get_var1_longlong, numVars);
    NC_TEST1(ncmpi_get_var1_ulonglong, numVars);
    NC_TEST1(ncmpi_get_var1, numVars);
    NC_TEST1(ncmpi_get_vara_text, numVars);
    NC_TEST1(ncmpi_get_vara_schar, numVars);
    NC_TEST1(ncmpi_get_vara_short, numVars);
    NC_TEST1(ncmpi_get_vara_int, numVars);
    NC_TEST1(ncmpi_get_vara_long, numVars);
    NC_TEST1(ncmpi_get_vara_float, numVars);
    NC_TEST1(ncmpi_get_vara_double, numVars);
    NC_TEST1(ncmpi_get_vara_uchar, numVars);
    NC_TEST1(ncmpi_get_vara_ushort, numVars);
    NC_TEST1(ncmpi_get_vara_uint, numVars);
    NC_TEST1(ncmpi_get_vara_longlong, numVars);
    NC_TEST1(ncmpi_get_vara_ulonglong, numVars);
    NC_TEST1(ncmpi_get_vara, numVars);
    NC_TEST1(ncmpi_get_vars_text, numVars);
    NC_TEST1(ncmpi_get_vars_schar, numVars);
    NC_TEST1(ncmpi_get_vars_short, numVars);
    NC_TEST1(ncmpi_get_vars_int, numVars);
    NC_TEST1(ncmpi_get_vars_long, numVars);
    NC_TEST1(ncmpi_get_vars_float, numVars);
    NC_TEST1(ncmpi_get_vars_double, numVars);
    NC_TEST1(ncmpi_get_vars_uchar, numVars);
    NC_TEST1(ncmpi_get_vars_ushort, numVars);
    NC_TEST1(ncmpi_get_vars_uint, numVars);
    NC_TEST1(ncmpi_get_vars_longlong, numVars);
    NC_TEST1(ncmpi_get_vars_ulonglong, numVars);
    NC_TEST1(ncmpi_get_vars, numVars);
    NC_TEST1(ncmpi_get_varm_text, numVars);
    NC_TEST1(ncmpi_get_varm_schar, numVars);
    NC_TEST1(ncmpi_get_varm_short, numVars);
    NC_TEST1(ncmpi_get_varm_int, numVars);
    NC_TEST1(ncmpi_get_varm_long, numVars);
    NC_TEST1(ncmpi_get_varm_float, numVars);
    NC_TEST1(ncmpi_get_varm_double, numVars);
    NC_TEST1(ncmpi_get_varm_uchar, numVars);
    NC_TEST1(ncmpi_get_varm_ushort, numVars);
    NC_TEST1(ncmpi_get_varm_uint, numVars);
    NC_TEST1(ncmpi_get_varm_longlong, numVars);
    NC_TEST1(ncmpi_get_varm_ulonglong, numVars);
    NC_TEST1(ncmpi_get_varm, numVars);
    NC_TEST2(ncmpi_get_att_text, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_schar, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_short, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_int, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_long, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_float, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_double, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_uchar, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_ushort, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_uint, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_longlong, numGatts, numVars);
    NC_TEST2(ncmpi_get_att_ulonglong, numGatts, numVars);
    NC_TEST2(ncmpi_get_att,     numGatts, numVars);
    NC_TEST2(ncmpi_inq_att,     numGatts, numVars);
    NC_TEST2(ncmpi_inq_attname, numGatts, numVars);
    NC_TEST2(ncmpi_inq_attid,   numGatts, numVars);
    NC_TEST2(ncmpi_inq_attlen,  numGatts, numVars);
    NC_TEST2(ncmpi_inq_atttype, numGatts, numVars);

    /* nonblocking I/O (not supported in netCDF4) */
    if (cdf_format != 4) {
        NC_TEST1(ncmpi_iget_var_text, numVars);
        NC_TEST1(ncmpi_iget_var_schar, numVars);
        NC_TEST1(ncmpi_iget_var_short, numVars);
        NC_TEST1(ncmpi_iget_var_int, numVars);
        NC_TEST1(ncmpi_iget_var_long, numVars);
        NC_TEST1(ncmpi_iget_var_float, numVars);
        NC_TEST1(ncmpi_iget_var_double, numVars);
        NC_TEST1(ncmpi_iget_var_uchar, numVars);
        NC_TEST1(ncmpi_iget_var_ushort, numVars);
        NC_TEST1(ncmpi_iget_var_uint, numVars);
        NC_TEST1(ncmpi_iget_var_longlong, numVars);
        NC_TEST1(ncmpi_iget_var_ulonglong, numVars);
        NC_TEST1(ncmpi_iget_var, numVars);
        NC_TEST1(ncmpi_iget_var1_text, numVars);
        NC_TEST1(ncmpi_iget_var1_schar, numVars);
        NC_TEST1(ncmpi_iget_var1_short, numVars);
        NC_TEST1(ncmpi_iget_var1_int, numVars);
        NC_TEST1(ncmpi_iget_var1_long, numVars);
        NC_TEST1(ncmpi_iget_var1_float, numVars);
        NC_TEST1(ncmpi_iget_var1_double, numVars);
        NC_TEST1(ncmpi_iget_var1_uchar, numVars);
        NC_TEST1(ncmpi_iget_var1_ushort, numVars);
        NC_TEST1(ncmpi_iget_var1_uint, numVars);
        NC_TEST1(ncmpi_iget_var1_longlong, numVars);
        NC_TEST1(ncmpi_iget_var1_ulonglong, numVars);
        NC_TEST1(ncmpi_iget_var1, numVars);
        NC_TEST1(ncmpi_iget_vara_text, numVars);
        NC_TEST1(ncmpi_iget_vara_schar, numVars);
        NC_TEST1(ncmpi_iget_vara_short, numVars);
        NC_TEST1(ncmpi_iget_vara_int, numVars);
        NC_TEST1(ncmpi_iget_vara_long, numVars);
        NC_TEST1(ncmpi_iget_vara_float, numVars);
        NC_TEST1(ncmpi_iget_vara_double, numVars);
        NC_TEST1(ncmpi_iget_vara_uchar, numVars);
        NC_TEST1(ncmpi_iget_vara_ushort, numVars);
        NC_TEST1(ncmpi_iget_vara_uint, numVars);
        NC_TEST1(ncmpi_iget_vara_longlong, numVars);
        NC_TEST1(ncmpi_iget_vara_ulonglong, numVars);
        NC_TEST1(ncmpi_iget_vara, numVars);
        NC_TEST1(ncmpi_iget_vars_text, numVars);
        NC_TEST1(ncmpi_iget_vars_schar, numVars);
        NC_TEST1(ncmpi_iget_vars_short, numVars);
        NC_TEST1(ncmpi_iget_vars_int, numVars);
        NC_TEST1(ncmpi_iget_vars_long, numVars);
        NC_TEST1(ncmpi_iget_vars_float, numVars);
        NC_TEST1(ncmpi_iget_vars_double, numVars);
        NC_TEST1(ncmpi_iget_vars_uchar, numVars);
        NC_TEST1(ncmpi_iget_vars_ushort, numVars);
        NC_TEST1(ncmpi_iget_vars_uint, numVars);
        NC_TEST1(ncmpi_iget_vars_longlong, numVars);
        NC_TEST1(ncmpi_iget_vars_ulonglong, numVars);
        NC_TEST1(ncmpi_iget_vars, numVars);
        NC_TEST1(ncmpi_iget_varm_text, numVars);
        NC_TEST1(ncmpi_iget_varm_schar, numVars);
        NC_TEST1(ncmpi_iget_varm_short, numVars);
        NC_TEST1(ncmpi_iget_varm_int, numVars);
        NC_TEST1(ncmpi_iget_varm_long, numVars);
        NC_TEST1(ncmpi_iget_varm_float, numVars);
        NC_TEST1(ncmpi_iget_varm_double, numVars);
        NC_TEST1(ncmpi_iget_varm_uchar, numVars);
        NC_TEST1(ncmpi_iget_varm_ushort, numVars);
        NC_TEST1(ncmpi_iget_varm_uint, numVars);
        NC_TEST1(ncmpi_iget_varm_longlong, numVars);
        NC_TEST1(ncmpi_iget_varm_ulonglong, numVars);
        NC_TEST1(ncmpi_iget_varm, numVars);
    }
        /* Test write functions */
    if (! read_only) {
        /* delete scratch file and ignore the error if not exist */
        unlink(scratch);

        NC_TEST(ncmpi_create);
        NC_TEST2(ncmpi_redef, numGatts, numVars);
        /* NC_TEST(ncmpi_enddef);  redundant, as it calls test_ncmpi_redef() */
        NC_TEST2(ncmpi_sync, numGatts, numVars);
        if (cdf_format != 4)
            NC_TEST2(ncmpi_flush, numGatts, numVars);
        NC_TEST2(ncmpi_abort, numGatts, numVars);
        NC_TEST1(ncmpi_def_dim, numVars);
        NC_TEST(ncmpi_rename_dim);
        NC_TEST1(ncmpi_def_var, numVars);
        NC_TEST1(ncmpi_put_var_text, numVars);
        NC_TEST1(ncmpi_put_var_schar, numVars);
        NC_TEST1(ncmpi_put_var_short, numVars);
        NC_TEST1(ncmpi_put_var_int, numVars);
        NC_TEST1(ncmpi_put_var_long, numVars);
        NC_TEST1(ncmpi_put_var_float, numVars);
        NC_TEST1(ncmpi_put_var_double, numVars);
        NC_TEST1(ncmpi_put_var_uchar, numVars);
        NC_TEST1(ncmpi_put_var_ushort, numVars);
        NC_TEST1(ncmpi_put_var_uint, numVars);
        NC_TEST1(ncmpi_put_var_longlong, numVars);
        NC_TEST1(ncmpi_put_var_ulonglong, numVars);
        NC_TEST1(ncmpi_put_var1_text, numVars);
        NC_TEST1(ncmpi_put_var1_schar, numVars);
        NC_TEST1(ncmpi_put_var1_short, numVars);
        NC_TEST1(ncmpi_put_var1_int, numVars);
        NC_TEST1(ncmpi_put_var1_long, numVars);
        NC_TEST1(ncmpi_put_var1_float, numVars);
        NC_TEST1(ncmpi_put_var1_double, numVars);
        NC_TEST1(ncmpi_put_var1_uchar, numVars);
        NC_TEST1(ncmpi_put_var1_ushort, numVars);
        NC_TEST1(ncmpi_put_var1_uint, numVars);
        NC_TEST1(ncmpi_put_var1_longlong, numVars);
        NC_TEST1(ncmpi_put_var1_ulonglong, numVars);
        NC_TEST1(ncmpi_put_var1, numVars);
        NC_TEST1(ncmpi_put_vara_text, numVars);
        NC_TEST1(ncmpi_put_vara_schar, numVars);
        NC_TEST1(ncmpi_put_vara_short, numVars);
        NC_TEST1(ncmpi_put_vara_int, numVars);
        NC_TEST1(ncmpi_put_vara_long, numVars);
        NC_TEST1(ncmpi_put_vara_float, numVars);
        NC_TEST1(ncmpi_put_vara_double, numVars);
        NC_TEST1(ncmpi_put_vara_uchar, numVars);
        NC_TEST1(ncmpi_put_vara_ushort, numVars);
        NC_TEST1(ncmpi_put_vara_uint, numVars);
        NC_TEST1(ncmpi_put_vara_longlong, numVars);
        NC_TEST1(ncmpi_put_vara_ulonglong, numVars);
        NC_TEST1(ncmpi_put_vara, numVars);
        NC_TEST1(ncmpi_put_vars_text, numVars);
        NC_TEST1(ncmpi_put_vars_schar, numVars);
        NC_TEST1(ncmpi_put_vars_short, numVars);
        NC_TEST1(ncmpi_put_vars_int, numVars);
        NC_TEST1(ncmpi_put_vars_long, numVars);
        NC_TEST1(ncmpi_put_vars_float, numVars);
        NC_TEST1(ncmpi_put_vars_double, numVars);
        NC_TEST1(ncmpi_put_vars_uchar, numVars);
        NC_TEST1(ncmpi_put_vars_ushort, numVars);
        NC_TEST1(ncmpi_put_vars_uint, numVars);
        NC_TEST1(ncmpi_put_vars_longlong, numVars);
        NC_TEST1(ncmpi_put_vars_ulonglong, numVars);
        NC_TEST1(ncmpi_put_vars, numVars);
        NC_TEST1(ncmpi_put_varm_text, numVars);
        NC_TEST1(ncmpi_put_varm_schar, numVars);
        NC_TEST1(ncmpi_put_varm_short, numVars);
        NC_TEST1(ncmpi_put_varm_int, numVars);
        NC_TEST1(ncmpi_put_varm_long, numVars);
        NC_TEST1(ncmpi_put_varm_float, numVars);
        NC_TEST1(ncmpi_put_varm_double, numVars);
        NC_TEST1(ncmpi_put_varm_uchar, numVars);
        NC_TEST1(ncmpi_put_varm_ushort, numVars);
        NC_TEST1(ncmpi_put_varm_uint, numVars);
        NC_TEST1(ncmpi_put_varm_longlong, numVars);
        NC_TEST1(ncmpi_put_varm_ulonglong, numVars);
        NC_TEST1(ncmpi_put_varm, numVars);
        NC_TEST1(ncmpi_rename_var, numVars);
        NC_TEST2(ncmpi_put_att_text, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_schar, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_short, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_int, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_long, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_float, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_double, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_uchar, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_ushort, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_uint, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_longlong, numGatts, numVars);
        NC_TEST2(ncmpi_put_att_ulonglong, numGatts, numVars);
        NC_TEST2(ncmpi_put_att, numGatts, numVars);
        NC_TEST2(ncmpi_copy_att, numGatts, numVars);
        NC_TEST2(ncmpi_rename_att, numGatts, numVars);
        NC_TEST2(ncmpi_del_att, numGatts, numVars);
        NC_TEST1(ncmpi_set_fill, numVars);
        NC_TEST(ncmpi_delete);

        /* test nonblocking APIs */
        if (cdf_format != 4) {
            NC_TEST1(ncmpi_iput_var_text, numVars);
            NC_TEST1(ncmpi_iput_var_schar, numVars);
            NC_TEST1(ncmpi_iput_var_short, numVars);
            NC_TEST1(ncmpi_iput_var_int, numVars);
            NC_TEST1(ncmpi_iput_var_long, numVars);
            NC_TEST1(ncmpi_iput_var_float, numVars);
            NC_TEST1(ncmpi_iput_var_double, numVars);
            NC_TEST1(ncmpi_iput_var_uchar, numVars);
            NC_TEST1(ncmpi_iput_var_ushort, numVars);
            NC_TEST1(ncmpi_iput_var_uint, numVars);
            NC_TEST1(ncmpi_iput_var_longlong, numVars);
            NC_TEST1(ncmpi_iput_var_ulonglong, numVars);
            NC_TEST1(ncmpi_iput_var, numVars);
            NC_TEST1(ncmpi_iput_var1_text, numVars);
            NC_TEST1(ncmpi_iput_var1_schar, numVars);
            NC_TEST1(ncmpi_iput_var1_short, numVars);
            NC_TEST1(ncmpi_iput_var1_int, numVars);
            NC_TEST1(ncmpi_iput_var1_long, numVars);
            NC_TEST1(ncmpi_iput_var1_float, numVars);
            NC_TEST1(ncmpi_iput_var1_double, numVars);
            NC_TEST1(ncmpi_iput_var1_uchar, numVars);
            NC_TEST1(ncmpi_iput_var1_ushort, numVars);
            NC_TEST1(ncmpi_iput_var1_uint, numVars);
            NC_TEST1(ncmpi_iput_var1_longlong, numVars);
            NC_TEST1(ncmpi_iput_var1_ulonglong, numVars);
            NC_TEST1(ncmpi_iput_var1, numVars);
            NC_TEST1(ncmpi_iput_vara_text, numVars);
            NC_TEST1(ncmpi_iput_vara_schar, numVars);
            NC_TEST1(ncmpi_iput_vara_short, numVars);
            NC_TEST1(ncmpi_iput_vara_int, numVars);
            NC_TEST1(ncmpi_iput_vara_long, numVars);
            NC_TEST1(ncmpi_iput_vara_float, numVars);
            NC_TEST1(ncmpi_iput_vara_double, numVars);
            NC_TEST1(ncmpi_iput_vara_uchar, numVars);
            NC_TEST1(ncmpi_iput_vara_ushort, numVars);
            NC_TEST1(ncmpi_iput_vara_uint, numVars);
            NC_TEST1(ncmpi_iput_vara_longlong, numVars);
            NC_TEST1(ncmpi_iput_vara_ulonglong, numVars);
            NC_TEST1(ncmpi_iput_vara, numVars);
            NC_TEST1(ncmpi_iput_vars_text, numVars);
            NC_TEST1(ncmpi_iput_vars_schar, numVars);
            NC_TEST1(ncmpi_iput_vars_short, numVars);
            NC_TEST1(ncmpi_iput_vars_int, numVars);
            NC_TEST1(ncmpi_iput_vars_long, numVars);
            NC_TEST1(ncmpi_iput_vars_float, numVars);
            NC_TEST1(ncmpi_iput_vars_double, numVars);
            NC_TEST1(ncmpi_iput_vars_uchar, numVars);
            NC_TEST1(ncmpi_iput_vars_ushort, numVars);
            NC_TEST1(ncmpi_iput_vars_uint, numVars);
            NC_TEST1(ncmpi_iput_vars_longlong, numVars);
            NC_TEST1(ncmpi_iput_vars_ulonglong, numVars);
            NC_TEST1(ncmpi_iput_vars, numVars);
            NC_TEST1(ncmpi_iput_varm_text, numVars);
            NC_TEST1(ncmpi_iput_varm_schar, numVars);
            NC_TEST1(ncmpi_iput_varm_short, numVars);
            NC_TEST1(ncmpi_iput_varm_int, numVars);
            NC_TEST1(ncmpi_iput_varm_long, numVars);
            NC_TEST1(ncmpi_iput_varm_float, numVars);
            NC_TEST1(ncmpi_iput_varm_double, numVars);
            NC_TEST1(ncmpi_iput_varm_uchar, numVars);
            NC_TEST1(ncmpi_iput_varm_ushort, numVars);
            NC_TEST1(ncmpi_iput_varm_uint, numVars);
            NC_TEST1(ncmpi_iput_varm_longlong, numVars);
            NC_TEST1(ncmpi_iput_varm_ulonglong, numVars);
            NC_TEST1(ncmpi_iput_varm, numVars);
        }
        NC_TEST(ncmpi_set_default_format);
    }

fn_exit:
    MPI_Info_free(&info);

    if (nfailsTotal == 0)  {
        printf(PASS_STR);
    }
    else {
        print("\n%s: expects 0 failures ... ",argv[0]);
        printf(FAIL_STR, nfailsTotal);
    }
    MPI_Finalize();
    return nfailsTotal > 0;
}
