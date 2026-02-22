/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* This core subroutines of utility program ncmpidiff. It compares header and
 * variables of two files regardless the define order of the variables and
 * attributes.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h> /* stat() */
#include <sys/stat.h>  /* stat() */
#include <unistd.h>    /* stat() */
#include <errno.h>     /* errno */
#include <math.h>      /* INFINITY */
#include <assert.h>

#include <mpi.h>
#include <pnetcdf.h>

#include <ncmpidiff_core.h>

#ifndef ubyte
#define ubyte unsigned char
#endif
#ifndef ushort
#define ushort unsigned short
#endif
#ifndef uint
#define uint unsigned int
#endif
#ifndef int64
#define int64 long long
#endif
#ifndef uint64
#define uint64 unsigned long long
#endif

#define PRINT_CMD_OPTS \
    if (first_diff && cmd_opts != NULL) { \
        printf("%s\n", cmd_opts); \
        first_diff = 0; \
    }

#define OOM_ERROR { \
    fprintf(stderr, "Error: calloc() out of memory at line %d\n",__LINE__);  \
    exit(1); \
}

#define HANDLE_ERROR {                                                       \
    if (err != NC_NOERR) {                                                   \
        fprintf(stderr, "Error at line %d of file %s (%s)\n", __LINE__,      \
               __FILE__, ncmpi_strerror(err));                               \
        MPI_Abort(MPI_COMM_WORLD, -1);                                       \
        exit(-1);                                                            \
    }                                                                        \
}

#define HANDLE_FILE_ERR(filename) {                                          \
    if (err != NC_NOERR) {                                                   \
        fprintf(stderr, "Error at line %d: input file %s (%s)\n", __LINE__,  \
               filename, ncmpi_strerror(err));                               \
        MPI_Abort(MPI_COMM_WORLD, -1);                                       \
        exit(-1);                                                            \
    }                                                                        \
}

#define CHECK_GLOBAL_ATT_DIFF_CHAR {                                         \
    int pos;                                                                 \
    char *b1 = (char *)calloc((attlen[0] + 1) * 2, sizeof(char));            \
    char *b2 = b1 + attlen[0] + 1;                                           \
    if (!b1) OOM_ERROR                                                       \
    err = ncmpi_get_att_text(ncid[0], NC_GLOBAL, name[0], b1);               \
    HANDLE_ERROR                                                             \
    err = ncmpi_get_att_text(ncid[1], NC_GLOBAL, name[0], b2);               \
    HANDLE_ERROR                                                             \
    for (pos=0; pos<attlen[0]; pos++) {                                      \
        if (b1[pos] != b2[pos]) break;                                       \
    }                                                                        \
    if (pos != attlen[0]) {                                                  \
        char str[128], msg[1024];                                            \
        sprintf(msg, "DIFF: global ");                                       \
        sprintf(str, "attribute \"%s\" of type NC_CHAR at element %d of ",   \
                name[0], pos);                                               \
        strcat(msg, str);                                                    \
        sprintf(str, "value \"%s\" vs \"%s\"\n", b1, b2);                    \
        strcat(msg, str);                                                    \
        PRINT_CMD_OPTS \
        printf("%s", msg);                                                   \
        numHeadDIFF++;                                                       \
    }                                                                        \
    else if (verbose)                                                        \
        printf("\t\tSAME: attribute contents\n");                            \
    free(b1);                                                                \
    break;                                                                   \
}

#define CHECK_GLOBAL_ATT_DIFF(type, func) {                                  \
    int pos;                                                                 \
    type *b1 = (type *)calloc((attlen[0] + 1) * 2, sizeof(type));            \
    type *b2 = b1 + attlen[0] + 1;                                           \
    if (!b1) OOM_ERROR                                                       \
    err = func(ncid[0], NC_GLOBAL, name[0], b1);                             \
    HANDLE_ERROR                                                             \
    err = func(ncid[1], NC_GLOBAL, name[0], b2);                             \
    HANDLE_ERROR                                                             \
    for (pos=0; pos<attlen[0]; pos++) {                                      \
        if (b1[pos] != b2[pos]) break;                                       \
    }                                                                        \
    if (pos != attlen[0]) {                                                  \
        char str[128], msg[1024];                                            \
        sprintf(msg, "DIFF: global ");                                       \
        sprintf(str, "attribute \"%s\" of type \"%s\" at element %d of ",    \
                name[0], get_type(xtype[0]), pos);                           \
        strcat(msg, str);                                                    \
        sprintf(str, "value %g vs %g (difference = %e)\n",                   \
                (double)b1[pos],(double)b2[pos],(double)(b1[pos]-b2[pos]));  \
        strcat(msg, str);                                                    \
        PRINT_CMD_OPTS \
        printf("%s", msg);                                                   \
        numHeadDIFF++;                                                       \
    }                                                                        \
    else if (verbose)                                                        \
        printf("\t\tSAME: attribute contents\n");                            \
    free(b1);                                                                \
    break;                                                                   \
}

#define CHECK_VAR_ATT_DIFF_CHAR {                                            \
    int pos;                                                                 \
    char *b1 = (char *)calloc((attlen[0] + 1) * 2, sizeof(char));            \
    char *b2 = b1 + attlen[0] + 1;                                           \
    if (!b1) OOM_ERROR                                                       \
    err = ncmpi_get_att_text(ncid[0], varid[0], attrname, b1);               \
    HANDLE_ERROR                                                             \
    err = ncmpi_get_att_text(ncid[1], varid[1], attrname, b2);               \
    HANDLE_ERROR                                                             \
    for (pos=0; pos<attlen[0]; pos++) {                                      \
        if (b1[pos] != b2[pos]) break;                                       \
    }                                                                        \
    if (pos != attlen[0]) {                                                  \
        char str[1024], msg[1024];                                           \
        sprintf(msg, "DIFF: variable \"%s\" ", name[0]);                     \
        sprintf(str, "attribute \"%s\" of type NC_CHAR at element %d of ",   \
                attrname, pos);                                              \
        strcat(msg, str);                                                    \
        sprintf(str, "value \"%s\" vs \"%s\"\n", b1, b2);                    \
        strcat(msg, str);                                                    \
        PRINT_CMD_OPTS \
        printf("%s", msg);                                                   \
        numHeadDIFF++;                                                       \
    }                                                                        \
    else if (verbose)                                                        \
        printf("\t\tSAME: attribute contents\n");                            \
    free(b1);                                                                \
    break;                                                                   \
}

#define CHECK_VAR_ATT_DIFF(type, func) {                                     \
    int pos;                                                                 \
    type *b1 = (type *)calloc((attlen[0] + 1) * 2, sizeof(type));            \
    type *b2 = b1 + attlen[0] + 1;                                           \
    if (!b1) OOM_ERROR                                                       \
    err = func(ncid[0], varid[0], attrname, b1);                             \
    HANDLE_ERROR                                                             \
    err = func(ncid[1], varid[1], attrname, b2);                             \
    HANDLE_ERROR                                                             \
    for (pos=0; pos<attlen[0]; pos++) {                                      \
        if (b1[pos] != b2[pos]) break;                                       \
    }                                                                        \
    if (pos != attlen[0]) {                                                  \
        char str[1024], msg[1024];                                           \
        sprintf(msg, "DIFF: variable \"%s\" ", name[0]);                     \
        sprintf(str, "attribute \"%s\" of type \"%s\" at element %d of ",    \
                attrname, get_type(xtype[0]), pos);                          \
        strcat(msg, str);                                                    \
        sprintf(str, "value %g vs %g (difference = %e)\n",                   \
                (double)b1[pos],(double)b2[pos],(double)(b1[pos]-b2[pos]));  \
        strcat(msg, str);                                                    \
        PRINT_CMD_OPTS \
        printf("%s", msg);                                                   \
        numHeadDIFF++;                                                       \
    }                                                                        \
    else if (verbose)                                                        \
        printf("\t\tSAME: attribute contents\n");                            \
    free(b1);                                                                \
    break;                                                                   \
}

#define  ABS(x) ((x) >= 0) ? (x) : (-x)
#define UABS(x) (x)

#define CHECK_VAR_DIFF(type, func, xabs) {                                   \
    int pos, isDiff, worst = -1;                                             \
    type *b1, *b2;                                                           \
    b1 = (type *)calloc(varsize * 2, sizeof(type));                          \
    if (!b1) OOM_ERROR                                                       \
    b2 = b1 + varsize;                                                       \
    err = ncmpi_get_vara_##func(ncid[0], varid1, start, shape, b1);          \
    HANDLE_ERROR                                                             \
    err = ncmpi_get_vara_##func(ncid[1], varid2, start, shape, b2);          \
    HANDLE_ERROR                                                             \
    if (!check_tolerance) {                                                  \
        for (pos=0; pos<varsize; pos++) {                                    \
            if (b1[pos] != b2[pos])                                          \
                break;                                                       \
        }                                                                    \
    } else {                                                                 \
        for (pos=0; pos<varsize; pos++) {                                    \
            double abs_b1, abs_b2, abs_max, diff, ratio;                     \
            if ( b1[pos] == b2[pos] ) continue;                              \
            abs_b1 = xabs(b1[pos]);                                          \
            abs_b2 = xabs(b2[pos]);                                          \
            abs_max = (abs_b1 > abs_b2) ? abs_b1 : abs_b2;                   \
            diff = b1[pos] - b2[pos];                                        \
            diff = (diff >= 0) ? diff : -diff;                               \
            ratio = diff /  abs_max;                                         \
            if (diff <= tolerance_difference || ratio <= tolerance_ratio)    \
                continue;                                                    \
            /* fail to meet both tolerance errors */                         \
            worst = pos;                                                     \
            break;                                                           \
        }                                                                    \
    }                                                                        \
    if (pos != varsize || worst != -1) { /* diff is found */                 \
        double v1, v2;                                                       \
        if (ndims[0] == 0) { /* scalar variable */                           \
            PRINT_CMD_OPTS \
            if (worst == -1)                                                 \
                printf("DIFF: scalar variable \"%s\" of type \"%s\"\n",      \
                       name[0], get_type(xtype[0]));                         \
            else {                                                           \
                v1 = b1[worst];                                              \
                v2 = b2[worst];                                              \
                printf("DIFF (tolerance): scalar variable \"%s\" of type \"%s\" of value %g vs %g (difference = %e)\n", \
                       name[0], get_type(xtype[0]), v1, v2, v1-v2);          \
            }                                                                \
        } else {                                                             \
            int _i;                                                          \
            MPI_Offset *diffStart;                                           \
            diffStart = (MPI_Offset*) malloc(sizeof(MPI_Offset) * ndims[0]); \
            if (worst != -1) pos = worst;                                    \
            v1 = b1[pos];                                                    \
            v2 = b2[pos];                                                    \
            for (_i=ndims[0]-1; _i>=0; _i--) {                               \
                diffStart[_i] = pos % shape[_i] + start[_i];                 \
                pos /= shape[_i];                                            \
            }                                                                \
            PRINT_CMD_OPTS \
            if (worst == -1)                                                 \
                printf("DIFF: variable \"%s\" of type \"%s\" at element ["OFFFMT, \
                       name[0], get_type(xtype[0]), diffStart[0]);           \
            else                                                             \
                printf("DIFF (tolerance): variable \"%s\" of type \"%s\" at element ["OFFFMT, \
                       name[0], get_type(xtype[0]), diffStart[0]);           \
            for (_i=1; _i<ndims[0]; _i++)                                    \
                printf(", "OFFFMT, diffStart[_i]);                           \
            printf("] of value %g vs %g (difference = %e)\n", v1,v2,v1-v2);  \
            free(diffStart);                                                 \
        }                                                                    \
        numVarDIFF++;                                                        \
        pos = 1;                                                             \
    } else                                                                   \
        pos = 0;                                                             \
    MPI_Allreduce(&pos, &isDiff, 1, MPI_INT, MPI_MAX, comm);                 \
    if (isDiff == 0 && !rank && verbose)                                     \
        printf("\tSAME: variable \"%s\" contents\n",name[0]);                \
    free(b1);                                                                \
    break;                                                                   \
}

#ifndef EXIT_FAILURE
#ifndef vms
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#else
/* In OpenVMS, success is indicated by odd values and failure by even values. */
#define EXIT_SUCCESS 1
#define EXIT_FAILURE 0
#endif
#endif

/*----< get_type() >----------------------------------------------------------*/
static char*
get_type(int type)
{
    switch (type) {
        case NC_BYTE:   return "NC_BYTE";
        case NC_CHAR:   return "NC_CHAR";
        case NC_SHORT:  return "NC_SHORT";
        case NC_INT:    return "NC_INT";
        case NC_FLOAT:  return "NC_FLOAT";
        case NC_DOUBLE: return "NC_DOUBLE";
        case NC_UBYTE:  return "NC_UBYTE";
        case NC_USHORT: return "NC_USHORT";
        case NC_UINT:   return "NC_UINT";
        case NC_INT64:  return "NC_INT64";
        case NC_UINT64: return "NC_UINT64";
    }
    return "NC_NAT";
}

/*----< ncmpidiff_core() >---------------------------------------------------*/
MPI_Offset ncmpidiff_core(const char  *file1,
                          const char  *file2,
                          MPI_Comm     comm,
                          MPI_Info     info,
                          int          verbose,
                          int          quiet,
                          int          check_header,
                          int          check_variable_list,
                          int          check_entire_file,
                          int          num_vars,
                          char       **var_names,
                          int          check_tolerance,
                          int          first_diff,
                          char        *cmd_opts,
                          double       tolerance_difference,
                          double       tolerance_ratio)
{
    char *name[2];
    int i, j, err, rank, nprocs;
    int ncid[2], ndims[2], nvars[2], natts[2], recdim[2], *dimids[2], fmt[2];
    int cmp_nvars;
    long long numVarDIFF=0, numHeadDIFF=0, varDIFF, numDIFF;
    MPI_Offset *shape=NULL, varsize, *start=NULL;
    MPI_Offset attlen[2], dimlen[2];
    nc_type xtype[2];

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    ncid[0] = ncid[1] = -1;

    if (verbose && rank == 0) {
        printf("First  file: %s\n", file1);
        printf("Second file: %s\n", file2);
    }

    if (strcmp(file1, file2) == 0) {
        if (rank == 0)
            printf("Error: two input file names are identical (%s) ... exit\n", file1);
        return 1;
    }

    /* compare file format */
    err = ncmpi_inq_file_format(file1, &fmt[0]);
    HANDLE_FILE_ERR(file1)
    err = ncmpi_inq_file_format(file2, &fmt[1]);
    HANDLE_FILE_ERR(file2)

    if (fmt[0] != fmt[1]) {
        if (!quiet && rank == 0)
            printf("DIFF: file format (CDF-%d) != (CDF-%d)\n",fmt[0], fmt[1]);
        numHeadDIFF++;
        /* even formats are different, we continue to compare the contents
         * of the files (headers and variables).
         */
    }

    /* open files and retrieve headers into memory buffers */
    name[0] = (char*) calloc(NC_MAX_NAME, 1);
    if (!name[0]) OOM_ERROR
    name[1] = (char*) calloc(NC_MAX_NAME, 1);
    if (!name[1]) OOM_ERROR

    /* open files */
    err = ncmpi_open(comm, file1, NC_NOWRITE, info, &ncid[0]);
    HANDLE_ERROR
    err = ncmpi_open(comm, file2, NC_NOWRITE, info, &ncid[1]);
    HANDLE_ERROR

    /* retrieve metadata */
    err = ncmpi_inq(ncid[0], &ndims[0], &nvars[0], &natts[0], &recdim[0]);
    HANDLE_ERROR
    err = ncmpi_inq(ncid[1], &ndims[1], &nvars[1], &natts[1], &recdim[1]);
    HANDLE_ERROR

    /* compare file header */
    if (check_header && rank == 0) { /* only root checks header */
        int attnump;

        /* compare number of dimensions defined */
        if (ndims[0] != ndims[1]) {
            if (!quiet)
                printf("DIFF: number of dimensions (%d) != (%d)\n",ndims[0], ndims[1]);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of dimensions (%d)\n",ndims[0]);

        /* compare number of variables defined */
        if (nvars[0] != nvars[1]) {
            if (!quiet)
                printf("DIFF: number of variables (%d) != (%d)\n",nvars[0], nvars[1]);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of variables (%d)\n",nvars[0]);

        /* compare number of global attributes defined */
        if (natts[0] != natts[1]) {
            if (!quiet)
                printf("DIFF: number of global attributes (%d) != (%d)\n",natts[0], natts[1]);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of global attributes (%d)\n",natts[0]);

        /* compare attributes defined in 1st file and also in 2nd file */
        for (i=0; i<natts[0]; i++) {
            err = ncmpi_inq_attname(ncid[0], NC_GLOBAL, i, name[0]);
            HANDLE_ERROR
            /* find the attr with the same name from ncid[1] */
            err = ncmpi_inq_attid(ncid[1], NC_GLOBAL, name[0], &attnump);
            if (err == NC_ENOTATT) {
                if (!quiet)
                    printf("DIFF: global attribute \"%s\" defined in %s not found in %s\n",
                           name[0],file1,file2);
                numHeadDIFF++;
                continue; /* loop i */
            }

            err = ncmpi_inq_att(ncid[0], NC_GLOBAL, name[0], &xtype[0], &attlen[0]);
            HANDLE_ERROR
            err = ncmpi_inq_att(ncid[1], NC_GLOBAL, name[0], &xtype[1], &attlen[1]);
            HANDLE_ERROR

            /* compare attribute xtype */
            if (xtype[0] != xtype[1]) {
                if (!quiet)
                    printf("DIFF: global attribute \"%s\" data type (%s) != (%s)\n",
                           name[0],get_type(xtype[0]),get_type(xtype[1]));
                numHeadDIFF++;
                continue; /* loop i */
            }
            else if (verbose) {
                printf("Global attribute \"%s\":\n",name[0]);
                printf("\tSAME: data type (%s)\n",get_type(xtype[0]));
            }

            /* compare attribute length */
            if (attlen[0] != attlen[1]) {
                if (!quiet)
                    printf("DIFF: global attribute \"%s\" length ("OFFFMT") != ("OFFFMT")\n",
                           name[0],attlen[0],attlen[1]);
                numHeadDIFF++;
                continue; /* loop i */
            }
            else if (verbose)
                printf("\tSAME: length ("OFFFMT")\n",attlen[0]);

            /* compare attribute contents */
            switch (xtype[0]) {
                case NC_CHAR:   CHECK_GLOBAL_ATT_DIFF_CHAR
                case NC_SHORT:  CHECK_GLOBAL_ATT_DIFF(short,  ncmpi_get_att_short);
                case NC_INT:    CHECK_GLOBAL_ATT_DIFF(int,    ncmpi_get_att_int);
                case NC_FLOAT:  CHECK_GLOBAL_ATT_DIFF(float,  ncmpi_get_att_float);
                case NC_DOUBLE: CHECK_GLOBAL_ATT_DIFF(double, ncmpi_get_att_double);
                case NC_UBYTE:  CHECK_GLOBAL_ATT_DIFF(ubyte,  ncmpi_get_att_uchar);
                case NC_USHORT: CHECK_GLOBAL_ATT_DIFF(ushort, ncmpi_get_att_ushort);
                case NC_UINT:   CHECK_GLOBAL_ATT_DIFF(uint,   ncmpi_get_att_uint);
                case NC_INT64:  CHECK_GLOBAL_ATT_DIFF(int64,  ncmpi_get_att_longlong);
                case NC_UINT64: CHECK_GLOBAL_ATT_DIFF(uint64, ncmpi_get_att_ulonglong);
                default: ; /* TODO: handle unexpected types */
            }
        }

        /* check global attributes defined in 2nd file but not in 1st file */
        for (i=0; i<natts[1]; i++) {
            err = ncmpi_inq_attname(ncid[1], NC_GLOBAL, i, name[1]);
            HANDLE_ERROR
            /* find the attr with the same name from ncid[0] */
            if (ncmpi_inq_attid(ncid[0], NC_GLOBAL, name[1], &attnump) == NC_ENOTATT) {
                if (!quiet)
                    printf("DIFF: global attribute \"%s\" defined in %s not found in %s\n",
                           name[1],file2,file1);
                numHeadDIFF++;
            }
        }

        /* Compare dimensions */
        if (ndims[0] > 0 && ndims[1] > 0) {
            if (verbose)
                printf("Dimension:\n");
        } else
            goto cmp_vars;

        /* check dimensions in 1st file also appear in 2nd file */
        for (i=0; i<ndims[0]; i++) {
            int dimid;
            err = ncmpi_inq_dim(ncid[0], i, name[0], &dimlen[0]);
            HANDLE_ERROR
            /* find the dim with the same name from ncid[1] */
            err = ncmpi_inq_dimid(ncid[1], name[0], &dimid);
            if (err == NC_EBADDIM) {
                if (!quiet)
                    printf("DIFF: dimension \"%s\" defined in %s not found in %s\n",
                           name[0],file1,file2);
                numHeadDIFF++;
                continue; /* loop i */
            }

            /* compare dimension length */
            err = ncmpi_inq_dimlen(ncid[1], dimid, &dimlen[1]);
            HANDLE_ERROR
            if (dimlen[0] != dimlen[1]) {
                /* cast to quiet warning on 32 bit platforms */
                if (!quiet)
                    printf("DIFF: dimension \"%s\" length ("OFFFMT") != ("OFFFMT")\n",
                           name[0],dimlen[0],dimlen[1]);
                numHeadDIFF++;
            }
            else if (verbose)
                printf("\tSAME: dimension \"%s\" length ("OFFFMT")\n",
                       name[0],dimlen[0]);
        }

        /* check dimensions in 2nd file but not in 1st file */
        for (i=0; i<ndims[1]; i++) {
            int dimid;
            /* inquire dimension name from 2nd file */
            err = ncmpi_inq_dimname(ncid[1], i, name[1]);
            HANDLE_ERROR
            /* check if the dim name exists in the 1st file */
            if (ncmpi_inq_dimid(ncid[0], name[1], &dimid) == NC_EBADDIM) {
                if (!quiet)
                    printf("DIFF: dimension \"%s\" defined in %s not found in %s\n",
                           name[1],file2,file1);
                numHeadDIFF++;
            }
        }

        /* Compare variables' metadata */
cmp_vars:
        if (nvars[0] > 0 && nvars[1] > 0) {
            if (verbose)
                printf("Variables:\n");
        } else
            goto cmp_exit;

        /* check variables defined in 1st file and also in 2nd file */
        for (i=0; i<nvars[0]; i++) {
            int varid[2];

            varid[0] = i;
            err = ncmpi_inq_varndims(ncid[0], i, &ndims[0]); HANDLE_ERROR
            dimids[0] = (int*) calloc((size_t)ndims[0], sizeof(int));
            if (!dimids[0]) OOM_ERROR
            err = ncmpi_inq_var(ncid[0], i, name[0], &xtype[0], &ndims[0], dimids[0], &natts[0]);
            HANDLE_ERROR
            /* find the variable with the same name from ncid[1] */
            err = ncmpi_inq_varid(ncid[1], name[0], &varid[1]);
            if (err == NC_ENOTVAR) {
                if (!quiet)
                    printf("DIFF: variable \"%s\"defined in %s not found in %s\n",
                           name[0],file1,file2);
                numHeadDIFF++;
                numVarDIFF++;
                continue;
            }

            /* inquire variable metadata for varid[1] */
            err = ncmpi_inq_varndims(ncid[1], varid[1], &ndims[1]); HANDLE_ERROR
            dimids[1] = (int*) calloc((size_t)ndims[1], sizeof(int));
            if (!dimids[1]) OOM_ERROR
            err = ncmpi_inq_var(ncid[1], varid[1], name[1], &xtype[1], &ndims[1], dimids[1], &natts[1]);
            HANDLE_ERROR

            /* compare variable xtype */
            if (xtype[0] != xtype[1]) {
                if (!quiet)
                    printf("DIFF: variable \"%s\" data type (%s) != (%s)\n",
                           name[0],get_type(xtype[0]),get_type(xtype[1]));
                numHeadDIFF++;
            }
            else if (verbose) {
                printf("Variable \"%s\":\n",name[0]);
                printf("\tSAME: data type (%s)\n",get_type(xtype[0]));
            }

            /* compare variable ndims */
            if (ndims[0] != ndims[1]) {
                if (!quiet)
                    printf("DIFF: variable \"%s\" number of dimensions (%d) != (%d)\n",
                           name[0],ndims[0],ndims[1]);
                numHeadDIFF++;
            }
            else {
                if (verbose)
                    printf("\tSAME: number of dimensions (%d)\n",ndims[0]);

                /* compare variable's dimensionality */
                for (j=0; j<ndims[0]; j++) {
                    char dimname[2][NC_MAX_NAME];
                    /* get dim name for each dim ID */
                    err = ncmpi_inq_dim(ncid[0], dimids[0][j], dimname[0], &dimlen[0]);
                    HANDLE_ERROR
                    err = ncmpi_inq_dim(ncid[1], dimids[1][j], dimname[1], &dimlen[1]);
                    HANDLE_ERROR
                    if (verbose)
                        printf("\tdimension %d:\n",j);
                    if (strcmp(dimname[0], dimname[1]) != 0) {
                        if (!quiet)
                            printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's name (%s) != (%s)\n",
                                   name[0],get_type(xtype[0]),j,dimname[0],dimname[1]);
                        numHeadDIFF++;
                    }
                    else if (verbose)
                        printf("\t\tSAME: name (%s)\n",dimname[0]);

                    /* compare variable dimension j's length */
                    if (dimlen[0] != dimlen[1]) {
                        if (!quiet)
                            printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's length ("OFFFMT") != ("OFFFMT")\n",
                                   name[0],get_type(xtype[0]),j,dimlen[0],dimlen[1]);
                        numHeadDIFF++;
                    }
                    else if (verbose)
                        printf("\t\tSAME: length ("OFFFMT")\n",dimlen[0]);
                }
            }

            /* compare number of attributes of this variable */
            if (natts[0] != natts[1]) {
                if (!quiet)
                    printf("DIFF: variable \"%s\" number of attributes (%d) != (%d)\n",
                           name[0],natts[0],natts[1]);
                numHeadDIFF++;
            }
            else if (verbose)
                printf("\tSAME: number of attributes (%d)\n",natts[0]);

            /* var attributes in 1st file also appear in 2nd file */
            for (j=0; j<natts[0]; j++) {
                char attrname[NC_MAX_NAME];
                err = ncmpi_inq_attname(ncid[0], i, j, attrname);
                HANDLE_ERROR
                err = ncmpi_inq_att(ncid[0], i, attrname, &xtype[0], &attlen[0]);
                HANDLE_ERROR
                /* find the variable attr with the same name from ncid[1] */
                err = ncmpi_inq_att(ncid[1], varid[1], attrname, &xtype[1], &attlen[1]);
                if (err == NC_ENOTATT) {
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" defined in %s not found in %s\n",
                               name[0],attrname,file1,file2);
                    numHeadDIFF++;
                    continue;
                }
                if (verbose)
                    printf("\tattribute \"%s\":\n",attrname);

                /* compare attribute xtype */
                if (xtype[0] != xtype[1]) {
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" data type (%s) != (%s)\n",
                               name[0],attrname,get_type(xtype[0]),get_type(xtype[1]));
                    numHeadDIFF++;
                    continue; /* skip this attribute */
                }
                else if (verbose)
                    printf("\t\tSAME: data type (%s)\n",get_type(xtype[0]));

                /* compare attribute nelems */
                if (attlen[0] != attlen[1]) {
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" length ("OFFFMT") != ("OFFFMT")\n",
                               name[0],attrname,attlen[0],attlen[1]);
                    numHeadDIFF++;
                    continue; /* skip this attribute */
                }
                else if (verbose)
                    printf("\t\tSAME: length ("OFFFMT")\n",attlen[0]);

                /* compare attribute contents */
                switch (xtype[0]) {
                    case NC_CHAR:   CHECK_VAR_ATT_DIFF_CHAR
                    case NC_SHORT:  CHECK_VAR_ATT_DIFF(short,  ncmpi_get_att_short);
                    case NC_INT:    CHECK_VAR_ATT_DIFF(int,    ncmpi_get_att_int);
                    case NC_FLOAT:  CHECK_VAR_ATT_DIFF(float,  ncmpi_get_att_float);
                    case NC_DOUBLE: CHECK_VAR_ATT_DIFF(double, ncmpi_get_att_double);
                    case NC_UBYTE:  CHECK_VAR_ATT_DIFF(ubyte,  ncmpi_get_att_uchar);
                    case NC_USHORT: CHECK_VAR_ATT_DIFF(ushort, ncmpi_get_att_ushort);
                    case NC_UINT:   CHECK_VAR_ATT_DIFF(uint,   ncmpi_get_att_uint);
                    case NC_INT64:  CHECK_VAR_ATT_DIFF(int64,  ncmpi_get_att_longlong);
                    case NC_UINT64: CHECK_VAR_ATT_DIFF(uint64, ncmpi_get_att_ulonglong);
                    default: ; /* TODO: handle unexpected types */
                }
            }

            /* check attributes in 2nd file but not in 1st file */
            for (j=0; j<natts[1]; j++) {
                char attrname[NC_MAX_NAME];
                err = ncmpi_inq_attname(ncid[1], varid[1], j, attrname);
                HANDLE_ERROR
                /* find the variable attr with the same name from ncid[0] */
                err = ncmpi_inq_att(ncid[0], i, attrname, &xtype[0], &attlen[0]);
                if (err == NC_ENOTATT) {
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" defined in %s not found in %s\n",
                               name[0],attrname,file2,file1);
                    numHeadDIFF++;
                }
            }
            free(dimids[0]);
            free(dimids[1]);
        }

        /* check variables defined in 2nd file but not in 1st file */
        for (i=0; i<nvars[1]; i++) { /* check variables in file2 but not in file1 */
            int varid;
            err = ncmpi_inq_varname(ncid[1], i, name[1]);
            HANDLE_ERROR
            /* find the variable with the same name from ncid[0] */
            err = ncmpi_inq_varid(ncid[0], name[1], &varid);
            if (err == NC_ENOTVAR) {
                if (!quiet)
                    printf("DIFF: variable \"%s\" defined in %s not found in %s\n",
                           name[1],file2,file1);
                numHeadDIFF++;
                numVarDIFF++;
            }
        }
    }

    /* compare variable contents */
    cmp_nvars = 0;
    if (check_variable_list) /* variable list is given at command line */
        cmp_nvars = num_vars;

    if (check_entire_file) { /* In this case, header has been checked */
        /* var_names is initialized to NULL */
        ncmpi_inq_nvars(ncid[0], &cmp_nvars);
        num_vars = cmp_nvars;
        var_names = (char**) calloc((size_t)cmp_nvars, sizeof(char*));
        if (!var_names) OOM_ERROR
        /* collect all the variable names from 1st file */
        for (i=0; i<cmp_nvars; i++) {
            ncmpi_inq_varname(ncid[0], i, name[0]);
            var_names[i] = strdup(name[0]);
            if (!var_names[i]) OOM_ERROR
        }
    }
    if (!rank && verbose) printf("number of variables to be compared = %d\n",cmp_nvars);

    /* compare variables, one at a time */
    for (i=0; i<cmp_nvars; i++) {
        int varid1, varid2;

        /* find variable ID in 1st file corresponding to var_names[i] */
        err = ncmpi_inq_varid(ncid[0], var_names[i], &varid1);
        if (err == NC_ENOTVAR) {
            if (!check_header) {
                if (!rank && !quiet)
                    printf("WARN: variable \"%s\" defined in %s not found in %s\n",
                           var_names[i],file2,file1);
                numVarDIFF++;
            }
            continue;
        }

        /* find variable ID in 2nd file corresponding to var_names[i] */
        err = ncmpi_inq_varid(ncid[1], var_names[i], &varid2);
        if (err == NC_ENOTVAR) {
            if (!check_header) {
                if (!rank && !quiet)
                    printf("WARN: variable \"%s\" defined in %s not found in %s\n",
                           var_names[i],file1,file2);
                numVarDIFF++;
            }
            continue;
        }

        /* Header comparison may have been skipped. Even if file headers have
         * been compared, we still need to compare variable's xtype and
         * dimensions to skip variables when their structures are different.
         */

        err = ncmpi_inq_varndims(ncid[0], varid1, &ndims[0]); HANDLE_ERROR
        dimids[0] = (int*) calloc((size_t)ndims[0], sizeof(int));
        if (!dimids[0]) OOM_ERROR
        err = ncmpi_inq_var(ncid[0], varid1, name[0], &xtype[0], &ndims[0], dimids[0], &natts[0]);
        HANDLE_ERROR
        err = ncmpi_inq_varndims(ncid[1], varid2, &ndims[1]); HANDLE_ERROR
        dimids[1] = (int*) calloc((size_t)ndims[1], sizeof(int));
        if (!dimids[1]) OOM_ERROR
        err = ncmpi_inq_var(ncid[1], varid2, name[1], &xtype[1], &ndims[1], dimids[1], &natts[1]);
        HANDLE_ERROR

        /* compare variable's NC data type */
        if (xtype[0] != xtype[1]) {
            if (!check_header) { /* if header has not been checked */
                if (!rank && !quiet)
                    printf("DIFF: variable \"%s\" data type (%s) != (%s)\n",
                           name[0],get_type(xtype[0]),get_type(xtype[1]));
                numHeadDIFF++;
                numVarDIFF++;
            }
            continue; /* skip this variable */
        }
        else if (!check_header && !rank && verbose) {
            printf("Variable \"%s\":\n",name[0]);
            printf("\tSAME: data type (%s)\n",get_type(xtype[0]));
        }

        /* compare variable's number of dimensions */
        if (ndims[0] != ndims[1]) {
            if (!check_header) { /* if header has not been checked */
                if (!rank && !quiet)
                    printf("DIFF: variable \"%s\" number of dimensions (%d) != (%d)\n",
                           name[0],ndims[0],ndims[1]);
                numHeadDIFF++;
                numVarDIFF++;
            }
            continue; /* skip this variable */
        }
        else if (!check_header && !rank && verbose)
            printf("\tSAME: number of dimensions (%d)\n",ndims[0]);

        shape = (MPI_Offset*) calloc((size_t)ndims[0] * 2, sizeof(MPI_Offset));
        if (!shape) OOM_ERROR
        start = shape + ndims[0];

        /* compare variable's dimension sizes, not dimension's names */
        for (j=0; j<ndims[0]; j++) {
            err = ncmpi_inq_dimlen(ncid[0], dimids[0][j], &dimlen[0]);
            HANDLE_ERROR
            err = ncmpi_inq_dimlen(ncid[1], dimids[1][j], &dimlen[1]);
            HANDLE_ERROR
            if (!check_header && !rank && verbose)
                printf("\tDimension %d:\n",j);
            if (dimlen[0] != dimlen[1]) {
                if (!check_header) { /* if header has not been checked */
                    if (!rank && !quiet)
                        printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's length ("OFFFMT") != ("OFFFMT")\n",
                               name[0],get_type(xtype[0]),j,dimlen[0],dimlen[1]);
                    numHeadDIFF++;
                    numVarDIFF++;
                }
                break; /* skip this variable */
            }
            else if (!check_header && !rank && verbose)
                printf("\t\tSAME: length ("OFFFMT")\n",dimlen[0]);
            shape[j] = dimlen[0];
        }
        if (j != ndims[0]) {
            free(shape);
            free(dimids[0]);
            free(dimids[1]);
            continue; /* skip this variable */
        }

        if (ndims[0] > 0 && dimids[0][0] == recdim[0]) { /* record variable */
            err = ncmpi_inq_dimlen(ncid[0], recdim[0], &shape[0]);
            HANDLE_ERROR
            if (shape[0] == 0) {
                /* No record has been written to the file, skip comparison */
                free(shape);
                free(dimids[0]);
                free(dimids[1]);
                continue;
            }
        }

        /* calculate read amount of this process in start[] and shape[] */
        for (j=0; j<ndims[0]; j++) {
            /* partition along dimension j among processes */
            if (shape[j] >= nprocs) {
                MPI_Offset dimLen = shape[j];
                shape[j] = dimLen / nprocs;
                start[j] = shape[j] * rank;
                if (rank < dimLen % nprocs) {
                    start[j] += rank;
                    shape[j]++;
                }
                else
                    start[j] += dimLen % nprocs;
                break;
            }
        }
        /* if none of shape[*] >= nprocs, then let all processes compare the
         * whole variable
         */

        varsize = 1;
        /* block partition the variable along the 1st dimension */
        for (j=0; j<ndims[0]; j++) varsize *= shape[j];

        /* compare the variable contents */
        switch (xtype[0]) {
            case NC_CHAR:   CHECK_VAR_DIFF(char,   text_all,       ABS);
            case NC_SHORT:  CHECK_VAR_DIFF(short,  short_all,      ABS);
            case NC_INT:    CHECK_VAR_DIFF(int,    int_all,        ABS);
            case NC_FLOAT:  CHECK_VAR_DIFF(float,  float_all,      ABS);
            case NC_DOUBLE: CHECK_VAR_DIFF(double, double_all,     ABS);
            case NC_UBYTE:  CHECK_VAR_DIFF(ubyte,  uchar_all,     UABS);
            case NC_USHORT: CHECK_VAR_DIFF(ushort, ushort_all,    UABS);
            case NC_UINT:   CHECK_VAR_DIFF(uint,   uint_all,      UABS);
            case NC_INT64:  CHECK_VAR_DIFF(int64,  longlong_all,   ABS);
            case NC_UINT64: CHECK_VAR_DIFF(uint64, ulonglong_all, UABS);
            default: ; /* TODO: handle unexpected types */
        }
        free(shape);
        free(dimids[0]);
        free(dimids[1]);
    }

cmp_exit:
    /* free up the memory previously allocated */
    if (check_entire_file) {
        for (i=0; i<num_vars; i++)
            free(var_names[i]);
        free(var_names);
    }

    free(name[0]);
    free(name[1]);

    for (i=0; i<2; i++) {
        /* close files */
        if (ncid[i] >= 0) {
            err = ncmpi_close(ncid[i]);
            HANDLE_ERROR
        }
    }

    /* summary of the difference */
    MPI_Reduce(&numVarDIFF, &varDIFF, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, comm);
    if (rank == 0 && !quiet) {
        if (check_header) {
            if (numHeadDIFF == 0)
                printf("Headers of two files are the same\n");
            else
                printf("Number of differences in header %lld\n",numHeadDIFF);
        }
        if (check_variable_list) {
            if (varDIFF == 0)
                printf("Compared variable(s) are the same\n");
            else
                printf("Compared variables(s) has %lld differences\n",varDIFF);
        }
        if (check_entire_file) {
            if (varDIFF == 0)
                printf("All variables of two files are the same\n");
            else
                printf("Number of differences in variables %lld\n",varDIFF);
        }
    }

    if (rank == 0) numDIFF = varDIFF + numHeadDIFF;
    MPI_Bcast(&numDIFF, 1, MPI_LONG_LONG_INT, 0, comm);

    return numDIFF;
}
