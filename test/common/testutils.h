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

#define MODE_COLL  1
#define MODE_INDEP 0

#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,ncmpi_strerrno(err)); \
    } \
}

#define CHECK_ERR_ALL { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,ncmpi_strerrno(err)); \
    } \
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); \
    if (nerrs > 0) goto fn_exit; \
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
        goto fn_exit; \
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
    if (nerrs > 0) goto fn_exit; \
}

#define CHECK_NERRS_ALL { \
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); \
    if (nerrs > 0) goto fn_exit; \
}

int inq_env_hint(char *hint_key, char **hint_value);

#ifdef PNETCDF_DEBUG
#define PASS_STR "\x1b[32mpass\x1b[0m\n"
#define SKIP_STR "\x1b[32mskip\x1b[0m\n"
#define FAIL_STR "\x1b[31mfail\x1b[0m with %d mismatches\n"
#else
#define PASS_STR "pass\n"
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
    char *in_path; /* input file path */
    int  num_fmts; /* number of file formats: CDF 1/2/4/5 */
    int *formats;  /* [num_fmts] max are {NC_FORMAT_CLASSIC,
                                          NC_FORMAT_64BIT_OFFSET,
                                          NC_FORMAT_NETCDF4_CLASSIC,
                                          NC_FORMAT_64BIT_DATA}; */
    int  ina;      /* add test of intra-node arrgatation */
    int  drv;      /* add test of MPI-IO driver */
    int  ind;      /* add test of hint romio_no_indep_rw */
    int  chk;      /* add test of hint nc_data_move_chunk_size */
    int  bb;       /* add test of burst-buffering feature */
    int  mod;      /* add test of independent data mode */
    int  hdr_diff; /* run ncmpidiff for file header only */
    int  var_diff; /* run ncmpidiff for variables */
} loop_opts;

extern
int tst_main(int argc, char **argv, char *msg, loop_opts opt,
             int (*tst_body)(const char*, const char*, int, int, MPI_Info));

#endif
