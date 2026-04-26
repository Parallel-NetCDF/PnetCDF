/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */


#ifndef H_UTILS
#define H_UTILS

#ifdef HAVE_CONFIG_H
#include <config.h> /* output of 'configure' */
#endif

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#if defined(HAVE_STDBOOL_H) && HAVE_STDBOOL_H == 1
#include <stdbool.h> /* type false and true */
typedef bool boolean;
#else
typedef int boolean;
#if !defined(__STDC_VERSION__) || __STDC_VERSION__ < 202311L
enum {false=0, true=1};
#endif
#endif

#define MODE_COLL  1
#define MODE_INDEP 0

#ifndef OFFFMT
#define OFFFMT "%lld"
#endif

#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,ncmpi_strerrno(err)); \
        assert(0); \
    } \
}

#define CHECK_ERR_ALL { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,ncmpi_strerrno(err)); \
    } \
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); \
    if (nerrs > 0) goto err_out; \
}

#define CHECK_ERROUT { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,ncmpi_strerrno(err)); \
        goto err_out; \
    } \
}

#define CHECK_FATAL_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,ncmpi_strerrno(err)); \
        assert(0); \
    } \
}

#define EXP_ERR(exp) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expecting %s but got %s\n", \
        __LINE__,__FILE__,ncmpi_strerrno(exp), ncmpi_strerrno(err)); \
    } \
}

#define CHECK_EXP_ERR_ALL(exp) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expecting %s but got %s\n", \
        __LINE__,__FILE__,ncmpi_strerrno(exp), ncmpi_strerrno(err)); \
    } \
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); \
    if (nerrs > 0) goto err_out; \
}

#define CHECK_NERRS_ALL if (nerrs != 0) assert(0);
/*
#define CHECK_NERRS_ALL { \
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); \
    if (nerrs > 0) assert(0) err_out; \
}
*/

int inq_env_hint(char *hint_key, char **hint_value);

#if PNETCDF_DEBUG_MODE == 1
#define PASS_STR "\x1b[32mpass\x1b[0m (%4.1fs)\n"
#define SKIP_STR "\x1b[32mskip\x1b[0m\n"
#define FAIL_STR "\x1b[31mfail\x1b[0m with %d mismatches\n"
#else
#define PASS_STR "pass (%4.1fs)\n"
#define SKIP_STR "skip\n"
#define FAIL_STR "fail with %d mismatches\n"
#endif

#define MPI_ERR(err) \
    if (err != MPI_SUCCESS) { \
        char err_string[MPI_MAX_ERROR_STRING+1]; \
        int  err_len; \
        MPI_Error_string(err, err_string, &err_len); \
        printf("MPI Error at file %s line %d (%s)\n",__FILE__,__LINE__,err_string); \
    }

#if defined(PNETCDF_DRIVER_NETCDF4) && PNETCDF_DRIVER_NETCDF4 == 1
extern int nc_formats[5];
#else
extern int nc_formats[3];
#endif

extern char* pnc_fmt_string(int format);

extern char* nc_err_code_name(int err);

#if MPI_VERSION < 3
/* MPI_OFFSET was first defined in MPI standard 2.2 */
#ifndef HAVE_DECL_MPI_OFFSET
#define MPI_OFFSET MPI_LONG_LONG_INT
#endif
#endif

#ifndef HAVE_STRDUP
extern char *strdup(const char *s);
#endif
#ifndef HAVE_STRCASECMP
extern int strcasecmp(const char *s1, const char *s2);
#endif

extern
char* remove_file_system_type_prefix(const char *filename);

extern
int is_relax_coord_bound(void);

extern
void tst_usage(char *argv0);

typedef struct {
    char *in_path;  /* input file path for read tests */
    int   num_fmts; /* number of file formats: CDF 1/2/3/4/5 */
    int  *formats;  /* [num_fmts] max are {NC_FORMAT_CLASSIC,
                                           NC_FORMAT_64BIT_OFFSET,
                                           NC_FORMAT_NETCDF4,
                                           NC_FORMAT_NETCDF4_CLASSIC,
                                           NC_FORMAT_64BIT_DATA}; */
    /* below 4 options:
     *      2 for testing both non-default and defaulte settings.
     *      1 for testing non-default setting only.
     *      0 for testing     defaulte setting only.
     */
    int  ina;  /* test of intra-node aggregation */
    int  drv;  /* test of GIO and MPI-IO drivers */
    int  ibuf; /* test of hint nc_ibuf_size */
    int  bb;   /* test of burst-buffering feature */
    int  mod;  /* test of independent data mode */

    boolean  hdr_diff; /* run ncmpidiff for file header */
    boolean  var_diff; /* run ncmpidiff for variables
                        * (disabled automatically when hdr_diff == false) */
} loop_opts;

extern
int tst_main(int argc, char **argv, char *msg, loop_opts opt,
             int (*tst_body)(const char*, const char*, int, int, MPI_Info));

#endif
