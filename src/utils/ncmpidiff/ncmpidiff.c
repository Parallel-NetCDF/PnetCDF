/*
 *  Copyright (C) 2010, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* This utility program compares header and variables of two files regardless
 * the define order of the variables and attributes. It can also compare a
 * subset of the variables, for example
 *    mpiexec -n 8 ncmpidiff -v var1,var2 file1.nc file2.nc
 *
 * or compare the header only, for example,
 *    mpiexec -n 8 ncmpidiff -h file1.nc file2.nc
 *
 * or compare header + a subset of variables, for example,
 *    mpiexec -n 8 ncmpidiff -h -v var1,var2 file1.nc file2.nc
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>   /* INFINITY */

#include <mpi.h>
#include <pnetcdf.h>

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

#define CHECK_GLOBAL_ATT_DIFF(type, func) {                                  \
    int pos;                                                                 \
    type *b1, *b2;                                                           \
    b1 = (type *)calloc((attlen[0] + 1) * 2, sizeof(type));                  \
    if (!b1) OOM_ERROR                                                       \
    b2 = b1 + attlen[0] + 1;                                                 \
    err = func(ncid[0], NC_GLOBAL, name[0], b1);                             \
    HANDLE_ERROR                                                             \
    err = func(ncid[1], NC_GLOBAL, name[0], b2);                             \
    HANDLE_ERROR                                                             \
    for (pos=0; pos<attlen[0]; pos++) {                                      \
        if (b1[pos] != b2[pos]) {                                            \
            char str[128], msg[1024];                                        \
            sprintf(msg, "DIFF: global attribute \"%s\" of type \"%s\" at element %d of value ", \
                    name[0], get_type(xtype[0]), pos);                       \
            if (xtype[0] == NC_CHAR)                                         \
                sprintf(str, "\"%s\" vs \"%s\"\n", b1, b2);                  \
            else                                                             \
                sprintf(str, "%g vs %g (difference = %e)\n", b1,b2,b1-b2);   \
            strcat(msg, str);                                                \
            printf("%s", msg);                                               \
            numHeadDIFF++;                                                   \
            break;                                                           \
        }                                                                    \
    }                                                                        \
    if (pos == attlen[0] && verbose)                                         \
        printf("\tSAME: attribute contents\n");                              \
    free(b1);                                                                \
    break;                                                                   \
}

#define CHECK_VAR_ATT_DIFF(type, func) {                                     \
    int pos;                                                                 \
    type *b1, *b2;                                                           \
    b1 = (type *)calloc(attlen[0] * 2, sizeof(type));                        \
    if (!b1) OOM_ERROR                                                       \
    b2 = b1 + attlen[0];                                                     \
    err = func(ncid[0], varid[0], attrname, b1);                             \
    HANDLE_ERROR                                                             \
    err = func(ncid[1], varid[1], attrname, b2);                             \
    HANDLE_ERROR                                                             \
    for (pos=0; pos<attlen[0]; pos++) {                                      \
        if (b1[pos] != b2[pos]) {                                            \
            char str[128], msg[1024];                                        \
            sprintf(msg, "DIFF: variable \"%s\" attribute \"%s\" of type \"%s\" at element %d of value ", \
                    name[0], attrname, get_type(xtype[0]), pos);             \
            if (xtype[0] == NC_CHAR)                                         \
                sprintf(str, "\"%s\" vs \"%s\"\n", b1, b2);                  \
            else                                                             \
                sprintf(str, "%g vs %g (difference = %e)\n", b1,b2,b1-b2);   \
            strcat(msg, str);                                                \
            printf("%s", msg);                                               \
            numHeadDIFF++;                                                   \
            break;                                                           \
        }                                                                    \
    }                                                                        \
    if (pos == attlen[0] && verbose)                                         \
        printf("\t\tSAME: attribute contents\n");                            \
    free(b1);                                                                \
    break;                                                                   \
}

#define CHECK_VAR_DIFF(type, func) {                                         \
    int pos, isDiff, worst = -1;                                             \
    type *b1, *b2;                                                           \
    b1 = (type *)calloc(varsize * 2, sizeof(type));                          \
    if (!b1) OOM_ERROR                                                       \
    b2 = b1 + varsize;                                                       \
    err = func(ncid[0], varid1, start, shape, b1);                           \
    HANDLE_ERROR                                                             \
    err = func(ncid[1], varid2, start, shape, b2);                           \
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
            abs_b1 = (b1[pos] >= 0) ? b1[pos] : -b1[pos];                    \
            abs_b2 = (b2[pos] >= 0) ? b2[pos] : -b2[pos];                    \
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
            diffStart = (MPI_Offset*) malloc(ndims[0] * sizeof(MPI_Offset)); \
            if (worst != -1) pos = worst;                                    \
            v1 = b1[pos];                                                    \
            v2 = b2[pos];                                                    \
            for (_i=ndims[0]-1; _i>=0; _i--) {                               \
                diffStart[_i] = pos % shape[_i] + start[_i];                 \
                pos /= shape[_i];                                            \
            }                                                                \
            if (worst == -1)                                                 \
                printf("DIFF: variable \"%s\" of type \"%s\" at element [%lld", \
                       name[0], get_type(xtype[0]), diffStart[0]);           \
            else                                                             \
                printf("DIFF (tolerance): variable \"%s\" of type \"%s\" at element [%lld", \
                       name[0], get_type(xtype[0]), diffStart[0]);           \
            for (_i=1; _i<ndims[0]; _i++)                                    \
                printf(", %lld", diffStart[_i]);                             \
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

/*----< usage() >-------------------------------------------------------------*/
static void
usage(int rank, char *progname)
{
#define USAGE   "\
  Compare the contents of two netCDF files.\n\
  [-b]             Verbose output\n\
  [-q]             quiet mode (no output if two files are the same)\n\
  [-h]             Compare header information only, no variables\n\
  [-v var1[,...]]  Compare variable(s) <var1>,... only\n\
  [-t diff,ratio]  Tolerance: diff is absolute element-wise difference\n\
                   and ratio is relative element-wise difference defined\n\
                   as |x - y|/max(|x|, |y|)\n\
  file1 file2      File names of two input netCDF files to be compared\n"

    if (rank == 0) {
        printf("  %s [-b] [-q] [-h] [-v ...] [-t diff,ratio] file1 file2\n%s",
               progname, USAGE);
        printf("*PnetCDF library version %s\n", ncmpi_inq_libvers());
    }
    MPI_Finalize();
    exit(1);
}

struct vspec {
    int    nvars;
    char **names; /* [nvars] */
};

/*----< get_var_names() >-----------------------------------------------------*/
static void
get_var_names(char *optarg, struct vspec* vspecp)
{
    char *cp=optarg, **cpp;
    int nvars = 1;

    /* compute number of variable names in comma-delimited list */
    vspecp->nvars = 1;
    while (*cp++)
        if (*cp == ',')
            nvars++;

    vspecp->names = (char **) calloc((size_t)nvars, sizeof(char*));
    if (!vspecp->names) OOM_ERROR

    cpp = vspecp->names;
    /* copy variable names into list */
    for (cp = strtok(optarg, ",");
         cp != NULL;
         cp = strtok((char *) NULL, ",")) {

        *cpp = strdup(cp);
        if (!*cpp) OOM_ERROR
        cpp++;
    }
    vspecp->nvars = nvars;
}

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

/*----< main() >--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    extern char *optarg;
    extern int optind;
    char *name[2];
    int i, j, c, err, rank, nprocs, verbose, quiet, check_tolerance;
    int ncid[2], ndims[2], nvars[2], natts[2], recdim[2], *dimids[2], fmt[2];
    int cmp_nvars, check_header, check_variable_list, check_entire_file;
    long long numVarDIFF=0, numHeadDIFF=0, varDIFF, numDIFF;
    double tolerance_ratio, tolerance_difference;
    MPI_Offset *shape=NULL, varsize, *start=NULL;
    MPI_Offset attlen[2], dimlen[2];
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    nc_type xtype[2];
    struct vspec var_list;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    verbose             = 0;
    quiet               = 0;
    check_header        = 0;
    check_variable_list = 0;
    check_entire_file   = 0;
    var_list.names      = NULL;
    var_list.nvars      = 0;
    check_tolerance     = 0;

    while ((c = getopt(argc, argv, "bhqt:v:")) != -1) {
        char *str, *ptr;
        switch(c) {
            case 'h':               /* compare header only */
                check_header = 1;
                break;
            case 'v':               /* variable names */
                /* make list of names of variables specified */
                get_var_names(optarg, &var_list);
                check_variable_list = 1;
                break;
            case 'b':
                verbose = 1;
                break;
            case 'q':
                quiet = 1;
                break;
            case 't':
                str = strdup(optarg);
                ptr = strtok(str, ",");
                if (ptr == NULL) {
                    usage(rank, argv[0]);
                    break;
                } else
                    sscanf(ptr, "%lf", &tolerance_difference);
                ptr = strtok(NULL, ",");
                if (ptr == NULL) {
                    usage(rank, argv[0]);
                    break;
                } else
                    sscanf(ptr, "%lf", &tolerance_ratio);
                check_tolerance = 1;
                free(str);
                break;
            case '?':
                usage(rank, argv[0]);
                break;
        }
    }

    /* quiet mode overwrites verbose */
    if (quiet) verbose = 0;

    if (argc - optind != 2) usage(rank, argv[0]);

    if (verbose && check_tolerance && rank == 0) {
        printf("Tolerance absolute difference = %e\n", tolerance_difference);
        printf("Tolerance ratio    difference = %e\n", tolerance_ratio);
    }

    if (check_header == 0 && check_variable_list == 0) {
        /* variable list is not provided, check header and all variables */
        check_entire_file = 1;
        check_header      = 1;
    }

    /* Nov. 18, 2014 -- disable subfiling as it does not correctly handle the
     * cases when  nprocs < num_subfiles */
    MPI_Info_create (&info);
    MPI_Info_set (info, "pnetcdf_subfiling", "disable");

    ncid[0] = ncid[1] = -1;

    /* open files and retrieve headers into memory buffers */
    for (i=0; i<2; i++) { /* i=0 for 1st file, i=1 for 2nd file */
        if (verbose && rank == 0) {
            if (i == 0) printf("First  file: %s\n", argv[optind+i]);
            if (i == 1) printf("Second file: %s\n", argv[optind+i]);
        }

        /* file format version */
        err = ncmpi_inq_file_format(argv[optind+i], &fmt[i]);
        HANDLE_ERROR

        if (fmt[i] == NC_FORMAT_NETCDF4 || fmt[i] == NC_FORMAT_NETCDF4_CLASSIC) {
#ifndef ENABLE_NETCDF4
            /* HDF5 files are not supported when --enable-netcdf4 is not used */
            if (rank == 0)
                fprintf(stderr, "Error: HDF5 based NetCDF4 file %s is not supported\n",
                        argv[optind+i]);
            numHeadDIFF++;
            quiet = 1;
            goto cmp_exit;
#endif
        } else if (fmt[i] == NC_FORMAT_BP) {
            /* BP files are not supported */
            if (rank == 0)
                fprintf(stderr, "Error: BP file %s is not supported\n",
                        argv[optind+i]);
            numHeadDIFF++;
            quiet = 1;
            goto cmp_exit;
        } else if (fmt[i] != NC_FORMAT_CLASSIC &&
                   fmt[i] != NC_FORMAT_CDF2 &&
                   fmt[i] != NC_FORMAT_CDF5) {
            /* valid classic NetCDF files are CDF-1, CDF-2, and CDF-5 */
            if (rank == 0)
                fprintf(stderr, "Error: %s is not a classic NetCDF file\n",
                        argv[optind+i]);
            numHeadDIFF++;
            quiet = 1;
            goto cmp_exit;
        }

        name[i] = (char*) calloc(NC_MAX_NAME, 1);
        if (!name[i]) OOM_ERROR

        /* open files */
        err = ncmpi_open(comm, argv[optind+i], NC_NOWRITE, info, &ncid[i]);
        HANDLE_ERROR

        err = ncmpi_inq(ncid[i], &ndims[i], &nvars[i], &natts[i], &recdim[i]);
        HANDLE_ERROR
    }

    /* compare file format */
    if (fmt[0] != fmt[1]) {
        if (!quiet && rank == 0)
            printf("DIFF: file format (CDF-%d) != (CDF-%d)\n",fmt[0], fmt[1]);
        numHeadDIFF++;
        /* even formats are different, we continue to compare the contents
         * of the files (headers and variables).
         */
    }

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
                           name[0],argv[optind],argv[optind+1]);
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
                    printf("DIFF: global attribute \"%s\" length (%lld) != (%lld)\n",
                           name[0],attlen[0],attlen[1]);
                numHeadDIFF++;
                continue; /* loop i */
            }
            else if (verbose)
                printf("\tSAME: length (%lld)\n",attlen[0]);

            /* compare attribute contents */
            switch (xtype[0]) {
                case NC_CHAR:   CHECK_GLOBAL_ATT_DIFF(char,   ncmpi_get_att_text);
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
                           name[1],argv[optind+1],argv[optind]);
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
                           name[0],argv[optind],argv[optind+1]);
                numHeadDIFF++;
                continue; /* loop i */
            }

            /* compare dimension length */
            err = ncmpi_inq_dimlen(ncid[1], dimid, &dimlen[1]);
            HANDLE_ERROR
            if (dimlen[0] != dimlen[1]) {
                /* cast to quiet warning on 32 bit platforms */
                if (!quiet)
                    printf("DIFF: dimension \"%s\" length (%lld) != (%lld)\n",
                           name[0],(long long int)dimlen[0],(long long int)dimlen[1]);
                numHeadDIFF++;
            }
            else if (verbose)
                printf("\tSAME: dimension \"%s\" length (%lld)\n",
                       name[0],(long long int)dimlen[0]);
        }

        /* check dimensions in 2nd file but not in 1st file */
        for (i=0; i<ndims[1]; i++) {
            int dimid;
            err = ncmpi_inq_dim(ncid[1], i, name[1], &dimlen[1]);
            HANDLE_ERROR
            /* find the dim with the same name from ncid[0] */
            if (ncmpi_inq_dimid(ncid[1], name[0], &dimid) == NC_EBADDIM) {
                if (!quiet)
                    printf("DIFF: dimension \"%s\" defined in %s not found in %s\n",
                           name[0],argv[optind+1],argv[optind]);
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
            dimids[0] = (int*) calloc((size_t)ndims[0], SIZEOF_INT);
            if (!dimids[0]) OOM_ERROR
            err = ncmpi_inq_var(ncid[0], i, name[0], &xtype[0], &ndims[0], dimids[0], &natts[0]);
            HANDLE_ERROR
            /* find the variable with the same name from ncid[1] */
            err = ncmpi_inq_varid(ncid[1], name[0], &varid[1]);
            if (err == NC_ENOTVAR) {
                if (!quiet)
                    printf("DIFF: variable \"%s\"defined in %s not found in %s\n",
                           name[0],argv[optind],argv[optind+1]);
                numHeadDIFF++;
                numVarDIFF++;
                continue;
            }

            /* inquire variable metadata for varid[1] */
            err = ncmpi_inq_varndims(ncid[1], varid[1], &ndims[1]); HANDLE_ERROR
            dimids[1] = (int*) calloc((size_t)ndims[1], SIZEOF_INT);
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
                            printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's length (%lld) != (%lld)\n",
                                   name[0],get_type(xtype[0]),j,(long long int)dimlen[0],(long long int)dimlen[1]);
                        numHeadDIFF++;
                    }
                    else if (verbose)
                        printf("\t\tSAME: length (%lld)\n",(long long int)dimlen[0]);
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
                               name[0],attrname,argv[optind],argv[optind+1]);
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
                        printf("DIFF: variable \"%s\" attribute \"%s\" length (%lld) != (%lld)\n",
                               name[0],attrname,(long long int)attlen[0],(long long int)attlen[1]);
                    numHeadDIFF++;
                    continue; /* skip this attribute */
                }
                else if (verbose)
                    printf("\t\tSAME: length (%lld)\n",(long long int)attlen[0]);

                /* compare attribute contents */
                switch (xtype[0]) {
                    case NC_CHAR:   CHECK_VAR_ATT_DIFF(char,   ncmpi_get_att_text);
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
                               name[0],attrname,argv[optind+1],argv[optind]);
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
                           name[1],argv[optind+1],argv[optind]);
                numHeadDIFF++;
                numVarDIFF++;
            }
        }
    }

    /* compare variable contents */
    cmp_nvars = 0;
    if (check_variable_list) /* variable list is given at command line */
        cmp_nvars = var_list.nvars;

    if (check_entire_file) { /* In this case, header has been checked */
        /* var_list.names is initialized to NULL */
        ncmpi_inq_nvars(ncid[0], &cmp_nvars);
        var_list.nvars = cmp_nvars;
        var_list.names = (char**) calloc((size_t)cmp_nvars, sizeof(char*));
        if (!var_list.names) OOM_ERROR
        /* collect all the variable names from 1st file */
        for (i=0; i<cmp_nvars; i++) {
            ncmpi_inq_varname(ncid[0], i, name[0]);
            var_list.names[i] = strdup(name[0]);
            if (!var_list.names[i]) OOM_ERROR
        }
    }
    if (!rank && verbose) printf("number of variables to be compared = %d\n",cmp_nvars);

    /* compare variables, one at a time */
    for (i=0; i<cmp_nvars; i++) {
        int varid1, varid2;

        /* find variable ID in 1st file corresponding to var_list.names[i] */
        err = ncmpi_inq_varid(ncid[0], var_list.names[i], &varid1);
        if (err == NC_ENOTVAR) {
            if (!check_header) {
                if (!rank && !quiet)
                    printf("WARN: variable \"%s\" defined in %s not found in %s\n",
                           var_list.names[i],argv[optind+1],argv[optind]);
                numVarDIFF++;
            }
            continue;
        }

        /* find variable ID in 2nd file corresponding to var_list.names[i] */
        err = ncmpi_inq_varid(ncid[1], var_list.names[i], &varid2);
        if (err == NC_ENOTVAR) {
            if (!check_header) {
                if (!rank && !quiet)
                    printf("WARN: variable \"%s\" defined in %s not found in %s\n",
                           var_list.names[i],argv[optind],argv[optind+1]);
                numVarDIFF++;
            }
            continue;
        }

        /* Header comparison may have been skipped. Even if file headers have
         * been compared, we still need to compare variable's xtype and
         * dimensions to skip variables when their structures are different.
         */

        err = ncmpi_inq_varndims(ncid[0], varid1, &ndims[0]); HANDLE_ERROR
        dimids[0] = (int*) calloc((size_t)ndims[0], SIZEOF_INT);
        if (!dimids[0]) OOM_ERROR
        err = ncmpi_inq_var(ncid[0], varid1, name[0], &xtype[0], &ndims[0], dimids[0], &natts[0]);
        HANDLE_ERROR
        err = ncmpi_inq_varndims(ncid[1], varid2, &ndims[1]); HANDLE_ERROR
        dimids[1] = (int*) calloc((size_t)ndims[1], SIZEOF_INT);
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

        shape = (MPI_Offset*) calloc((size_t)ndims[0] * 2, SIZEOF_MPI_OFFSET);
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
                        printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's length (%lld) != (%lld)\n",
                               name[0],get_type(xtype[0]),j,(long long int)dimlen[0],(long long int)dimlen[1]);
                    numHeadDIFF++;
                    numVarDIFF++;
                }
                break; /* skip this variable */
            }
            else if (!check_header && !rank && verbose)
                printf("\t\tSAME: length (%lld)\n",(long long int)dimlen[0]);
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
            case NC_CHAR:   CHECK_VAR_DIFF(char,   ncmpi_get_vara_text_all);
            case NC_SHORT:  CHECK_VAR_DIFF(short,  ncmpi_get_vara_short_all);
            case NC_INT:    CHECK_VAR_DIFF(int,    ncmpi_get_vara_int_all);
            case NC_FLOAT:  CHECK_VAR_DIFF(float,  ncmpi_get_vara_float_all);
            case NC_DOUBLE: CHECK_VAR_DIFF(double, ncmpi_get_vara_double_all);
            case NC_UBYTE:  CHECK_VAR_DIFF(ubyte,  ncmpi_get_vara_uchar_all);
            case NC_USHORT: CHECK_VAR_DIFF(ushort, ncmpi_get_vara_ushort_all);
            case NC_UINT:   CHECK_VAR_DIFF(uint,   ncmpi_get_vara_uint_all);
            case NC_INT64:  CHECK_VAR_DIFF(int64,  ncmpi_get_vara_longlong_all);
            case NC_UINT64: CHECK_VAR_DIFF(uint64, ncmpi_get_vara_ulonglong_all);
            default: ; /* TODO: handle unexpected types */
        }
        free(shape);
        free(dimids[0]);
        free(dimids[1]);
    }
    free(name[0]);
    free(name[1]);

    /* free up the memory previously allocated */
    if (var_list.nvars) {
        for (i=0; i<var_list.nvars; i++)
            free(var_list.names[i]);
        free(var_list.names);
    }

cmp_exit:
    for (i=0; i<2; i++) {
        /* close files */
        if (ncid[i] >= 0) {
            err = ncmpi_close(ncid[i]);
            HANDLE_ERROR
        }
    }
    if (info != MPI_INFO_NULL)
        MPI_Info_free(&info);

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

    MPI_Finalize();
    exit((numDIFF == 0) ? EXIT_SUCCESS : EXIT_FAILURE);
}
