/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* This utility program compares header and variables of two classic cdf-based
 * NetCDF files regardless the define order of the variables and attributes. It
 * can also compare a subset of the variables, for example
 *     cdfdiff -v var1,var2 file1.nc file2.nc
 *
 * or compare the header only, for example,
 *     cdfdiff -h file1.nc file2.nc
 *
 * or * compare header + a subset of variables, for example,
 *     cdfdiff -h -v var1,var2 file1.nc file2.nc
 */

#include <stdio.h>
#include <stdlib.h> /* malloc(), calloc(), free() */
#include <string.h> /* strtok(), strdup(), strlen(), strerror(), strcmp(), memcmp() */
#include <sys/types.h> /* lseek() */
#include <unistd.h> /* getopt(), lseek() */
#include <math.h>   /* INFINITY */

/* include subroutines from ncvalidator.c for reading file header */
#define BUILD_CDFDIFF
#include "ncvalidator.c"

/* divide variable into chunks, so diff is done in chunks */
#define READ_CHUNK_SIZE 4194304

#define OOM_ERROR { \
    fprintf(stderr, "Error: calloc() out of memory at line %d\n",__LINE__); \
    exit(1); \
}

#define CHECK_VAR_DIFF(itype) {                                          \
    size_t indx, esize, nelems;                                          \
    itype *b1 = (itype*)buf[0];                                          \
    itype *b2 = (itype*)buf[1];                                          \
    esize = sizeof(itype);                                               \
    nelems = rdLen[0] / esize;                                           \
    if (esize > 1 && is_little_endian) {                                 \
        swapn(b1, nelems, esize);                                        \
        swapn(b2, nelems, esize);                                        \
    }                                                                    \
    for (indx=0; indx<nelems; indx++) {                                  \
        double abs_b1, abs_b2, abs_max, diff, ratio;                     \
        if ( b1[indx] == b2[indx] ) continue;                            \
        abs_b1 = (b1[indx] >= 0) ? b1[indx] : -b1[indx];                 \
        abs_b2 = (b2[indx] >= 0) ? b2[indx] : -b2[indx];                 \
        abs_max = (abs_b1 > abs_b2) ? abs_b1 : abs_b2;                   \
        diff = b1[indx] - b2[indx];                                      \
        diff = (diff >= 0) ? diff : -diff;                               \
        ratio = diff /  abs_max;                                         \
        if (diff <= tolerance_difference || ratio <= tolerance_ratio)    \
            continue;                                                    \
        /* fail to meet both tolerance errors */                         \
        worst = chunk_off + indx; /* mark a difference is found */       \
        worst1 = b1[indx];                                               \
        worst2 = b2[indx];                                               \
        break;                                                           \
    }                                                                    \
    chunk_off += nelems;                                                 \
}

#define CHECK_VAR_DIFF_UNSIGNED(itype) {                                 \
    size_t indx, esize, nelems;                                          \
    itype *b1 = (itype*)buf[0];                                          \
    itype *b2 = (itype*)buf[1];                                          \
    esize = sizeof(itype);                                               \
    nelems = rdLen[0] / esize;                                           \
    if (esize > 1 && is_little_endian) {                                 \
        swapn(b1, nelems, esize);                                        \
        swapn(b2, nelems, esize);                                        \
    }                                                                    \
    for (indx=0; indx<nelems; indx++) {                                  \
        double abs_max, diff, ratio;                                     \
        if ( b1[indx] == b2[indx] ) continue;                            \
        abs_max = (b1[indx] > b2[indx]) ? b1[indx] : b2[indx];           \
        diff = b1[indx] - b2[indx];                                      \
        diff = (diff >= 0) ? diff : -diff;                               \
        ratio = diff /  abs_max;                                         \
        if (diff <= tolerance_difference || ratio <= tolerance_ratio)    \
            continue;                                                    \
        /* fail to meet both tolerance errors */                         \
        worst = chunk_off + indx; /* mark a difference is found */       \
        worst1 = b1[indx];                                               \
        worst2 = b2[indx];                                               \
        break;                                                           \
    }                                                                    \
    chunk_off += nelems;                                                 \
}

#define PRINTF_VAR_DIFF(itype) {                                         \
    int esize = sizeof(itype);                                           \
    itype *b1 = (itype*)str[0];                                          \
    itype *b2 = (itype*)str[1];                                          \
    if (esize > 1 && is_little_endian) {                                 \
        swapn(b1, 1, esize);                                             \
        swapn(b2, 1, esize);                                             \
    }                                                                    \
    worst1 = (double)b1[0];                                              \
    worst2 = (double)b2[0];                                              \
    printf("%g vs %g (difference = %e)\n", worst1,worst2,worst1-worst2); \
}

#define SWAP2B(a) ( (((a) & 0xff) << 8) | \
                    (((a) >> 8) & 0xff) )

/*----< swapn() >--------------------------------------------------*/
/* in-place byte swap, used when Endianness mismatches */
static void
swapn(void   *buf,
      size_t  nelems,  /* number of elements in buf[] */
      size_t  esize)   /* byte size of each element */
{
    size_t i;

    if (esize <= 1 || nelems <= 0) return;  /* no need */

    if (esize == 4) { /* this is the most common case */
        uint32_t *dest = (uint32_t*) buf;
        for (i=0; i<nelems; i++)
            dest[i] =  ((dest[i]) << 24)
                    | (((dest[i]) & 0x0000ff00) << 8)
                    | (((dest[i]) & 0x00ff0000) >> 8)
                    | (((dest[i]) >> 24));
    }
    else if (esize == 8) {
        uint64_t *dest = (uint64_t*) buf;
        for (i=0; i<nelems; i++)
            dest[i] = ((dest[i] & 0x00000000000000FFULL) << 56) |
                      ((dest[i] & 0x000000000000FF00ULL) << 40) |
                      ((dest[i] & 0x0000000000FF0000ULL) << 24) |
                      ((dest[i] & 0x00000000FF000000ULL) <<  8) |
                      ((dest[i] & 0x000000FF00000000ULL) >>  8) |
                      ((dest[i] & 0x0000FF0000000000ULL) >> 24) |
                      ((dest[i] & 0x00FF000000000000ULL) >> 40) |
                      ((dest[i] & 0xFF00000000000000ULL) >> 56);
    }
    else if (esize == 2) {
        uint16_t *dest = (uint16_t*) buf;
        for (i=0; i<nelems; i++)
            dest[i] = (uint16_t)(((dest[i] & 0xff) << 8) |
                                 ((dest[i] >> 8) & 0xff));
    }
    else {
        unsigned char tmp, *op = (unsigned char*)buf;
        /* for esize is not 1, 2, or 4 */
        while (nelems-- > 0) {
            for (i=0; i<esize/2; i++) {
                tmp           = op[i];
                op[i]         = op[esize-1-i];
                op[esize-1-i] = tmp;
            }
            op += esize;
        }
    }
}

/*----< usage() >-------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  Compare the contents of two files in classic netCDF formats.\n\
  [-b]             Verbose output\n\
  [-q]             quiet mode (no output if two files are the same)\n\
  [-h]             Compare header information only, no variables\n\
  [-v var1[,...]]  Compare variable(s) <var1>,... only\n\
  [-t diff,ratio]  Tolerance: diff is absolute element-wise difference\n\
                   and ratio is relative element-wise difference defined\n\
                   as |x - y|/max(|x|, |y|)\n\
  file1 file2      File names of two input netCDF files to be compared\n\
  *PnetCDF library version PNETCDF_RELEASE_VERSION of PNETCDF_RELEASE_DATE\n"

    printf("  %s [-b] [-q] [-h] [-v ...] [-t diff,ratio] file1 file2\n%s",
           progname, USAGE);
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
    char *str, *ptr;
    size_t nbytes;
    int i, j, k, m, n, c, err, verbose, quiet, isDiff;
    int fd[2], nvars[2], ndims[2], nattrs[2], check_tolerance;
    int cmp_nvars, check_header, check_variable_list, check_entire_file;
    long long numVarDIFF=0, numHeadDIFF=0, numDIFF;
    long long r, worst, chunk_off, numrecs;
    double tolerance_difference, tolerance_ratio;
    double worst1, worst2;
    void *buf[2] = {NULL, NULL};
    struct vspec var_list;
    NC *ncp[2];

    /* find Endianness of the running machine */
    int is_little_endian = check_little_endian();

    verbose             = 0;
    quiet               = 0;
    check_header        = 0;
    check_variable_list = 0;
    check_entire_file   = 0;
    var_list.names      = NULL;
    var_list.nvars      = 0;
    check_tolerance     = 0;

    while ((c = getopt(argc, argv, "bhqv:t:")) != -1)
        switch(c) {
            case 'h':               /* compare header only */
                check_header = 1;
                break;
            case 'v':               /* variable names */
                /* make list of names of variables specified */
                get_var_names(optarg, &var_list);
                check_variable_list = 1;
                break;
            case 't':
                str = strdup(optarg);
                ptr = strtok(str, ",");
                if (ptr == NULL) {
                    usage(argv[0]);
                    break;
                } else
                    sscanf(ptr, "%lf", &tolerance_difference);
                ptr = strtok(NULL, ",");
                if (ptr == NULL) {
                    usage(argv[0]);
                    break;
                } else
                    sscanf(ptr, "%lf", &tolerance_ratio);
                check_tolerance = 1;
                free(str);
                break;
            case 'b':
                verbose = 1;
                break;
            case 'q':
                quiet = 1;
                break;
            case '?':
                usage(argv[0]);
                break;
        }

    /* quiet mode overwrites verbose */
    if (quiet) verbose = 0;

    if (argc - optind != 2) usage(argv[0]);

    if (verbose && check_tolerance) {
        printf("Tolerance absolute difference = %e\n", tolerance_difference);
        printf("Tolerance ratio    difference = %e\n", tolerance_ratio);
    }

    if (check_header == 0 && check_variable_list == 0) {
        /* variable list is not provided, check header and all variables */
        check_entire_file = 1;
        check_header      = 1;
    }

    /* open files and retrieve headers into memory buffers */
    for (i=0; i<2; i++) { /* i=0 for 1st file, i=1 for 2nd file */
        const char *cdf_signature="CDF";
        const char *hdf5_signature="\211HDF\r\n\032\n";
        char signature[8];
        ssize_t rlen;

        /* open file */
        fd[i] = open(argv[optind+i], O_RDONLY);
        if (fd[i] == -1) {
            fprintf(stderr, "Error at line %d: open file %s (%s)\n",
                    __LINE__, argv[optind+i], strerror(errno));
            exit(1);
        }
        if (verbose && i == 0) printf("First  file: %s\n", argv[optind+i]);
        if (verbose && i == 1) printf("Second file: %s\n", argv[optind+i]);

        /* read file format signature, first 8 bytes of file */
        rlen = read(fd[i], signature, 8);
        if (rlen != 8) {
            close(fd[i]); /* ignore error */
            fprintf(stderr, "Error at line %d: reading file %s (%s)\n",
                    __LINE__, argv[optind+i], strerror(errno));
            exit(1);
        }

        /* HDF5 files are not supported */
        if (memcmp(signature, hdf5_signature, 8) == 0) {
            close(fd[i]); /* ignore error */
            fprintf(stderr, "Error: HDF5 based NetCDF4 file %s is not supported\n",
                    argv[optind+i]);
            exit(1);
        }
        else if (memcmp(signature, cdf_signature, 3) == 0) {
            /* classic NetCDF files */
            if (signature[3] != 1 && signature[3] != 2 && signature[3] != 5) {
                close(fd[i]); /* ignore error */
                fprintf(stderr, "Error: %s is not a classic NetCDF file\n",
                        argv[optind+i]);
                exit(1);
            }
        } else {
            close(fd[i]); /* ignore error */
            fprintf(stderr, "Error: %s is not a NetCDF file\n",argv[optind+i]);
            exit(1);
        }

        /* Allocate NC object which stores the entire file header */
        ncp[i] = (NC*) calloc(1, sizeof(NC));
        if (ncp[i] == NULL) OOM_ERROR

        /* read and validate the header */
        err = val_get_NC(fd[i], ncp[i]);
        if (err != NC_NOERR && err != NC_ENULLPAD && err != -1) {
            printf("File \"%s\" fails to conform with classic CDF-%d format specifications\n",
                   argv[optind+i], ncp[i]->format);
            exit(1);
        }

        nvars[i] = ncp[i]->vars.ndefined;
        ndims[i] = ncp[i]->dims.ndefined;
    }

    /* compare file format */
    if (ncp[0]->format != ncp[1]->format) {
        if (!quiet)
            printf("DIFF: file format (CDF-%d) != (CDF-%d)\n",
                   ncp[0]->format, ncp[1]->format);
        numHeadDIFF++;
        /* even formats are different, we continue to compare the contents
         * of the files (headers and variables).
         */
    }

    /* compare file header */
    if (check_header) {
        NC_attr *attr[2];
        NC_dim  *dim[2];
        NC_var  *var[2];

        /* compare number of dimensions defined */
        if (ndims[0] != ndims[1]) {
            if (!quiet)
                printf("DIFF: number of dimensions (%d) != (%d)\n",
                       ndims[0], ndims[1]);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of dimensions (%d)\n",ndims[0]);

        /* compare number of variables defined */
        if (nvars[0] != nvars[1]) {
            if (!quiet)
                printf("DIFF: number of variables (%d) != (%d)\n",
                       nvars[0], nvars[1]);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of variables (%d)\n", nvars[0]);

        /* compare number of global attributes defined */
        nattrs[0] = ncp[0]->attrs.ndefined;
        nattrs[1] = ncp[1]->attrs.ndefined;
        if (nattrs[0] != nattrs[1]) {
            if (!quiet)
                printf("DIFF: number of global attributes (%d) != (%d)\n",
                       nattrs[0], nattrs[1]);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of global attributes (%d)\n", nattrs[0]);

        /* compare attributes defined in 1st file and also in 2nd file */
        for (i=0; i<nattrs[0]; i++) {
            /* compare attribute name */
            attr[0] = ncp[0]->attrs.value[i];
            k = i % nattrs[1];
            for (j=0; j<nattrs[1]; j++) {
                attr[1] = ncp[1]->attrs.value[k];
                if (strcmp(attr[0]->name, attr[1]->name))
                    k = (k + 1) % nattrs[1];
                else
                    break; /* loop j */
            }
            if (j == nattrs[1]) { /* not found in 2nd file */
                if (!quiet)
                    printf("DIFF: global attribute \"%s\" defined in %s not found in %s\n",
                           attr[0]->name, argv[optind], argv[optind+1]);
                numHeadDIFF++;
                continue; /* loop i */
            }

            /* compare attribute xtype */
            if (attr[0]->xtype != attr[1]->xtype) {
                if (!quiet)
                    printf("DIFF: global attribute \"%s\" data type (%s) != (%s)\n",
                           attr[0]->name,
                           get_type(attr[0]->xtype), get_type(attr[1]->xtype));
                numHeadDIFF++;
                continue; /* loop i */
            }
            else if (verbose) {
                printf("Global attribute \"%s\":\n",attr[0]->name);
                printf("\tSAME: data type (%s)\n",get_type(attr[0]->xtype));
            }

            /* compare attribute length */
            if (attr[0]->nelems != attr[1]->nelems) {
                if (!quiet)
                    printf("DIFF: global attribute \"%s\" length (%lld) != (%lld)\n",
                           attr[0]->name, attr[0]->nelems, attr[1]->nelems);
                numHeadDIFF++;
                continue; /* loop i */
            }
            else if (verbose)
                printf("\tSAME: length (%lld)\n", attr[0]->nelems);

            /* compare attribute contents */
            nbytes = attr[0]->nelems * len_nctype(attr[0]->xtype);
            if ((isDiff = memcmp(attr[0]->xvalue, attr[1]->xvalue, nbytes)) != 0) {
                char *str[2];
                str[0] = (char*) attr[0]->xvalue;
                str[1] = (char*) attr[1]->xvalue;
                /* find the array index of first element in difference */
                for (m=0; m<nbytes; m++)
                    if (str[0][m] != str[1][m])
                        break;
                isDiff = m / len_nctype(attr[0]->xtype);
                if (!quiet) {
                    printf("DIFF: global attribute \"%s\" of type \"%s\" at element %d of value ",
                            attr[0]->name, get_type(attr[0]->xtype), isDiff);
                    if (attr[0]->xtype == NC_CHAR) {
                        char *_val[2];
                        _val[0] = (char*) calloc(attr[0]->nelems + 1, 1);
                        _val[1] = (char*) calloc(attr[1]->nelems + 1, 1);
                        strncpy(_val[0], str[0], attr[0]->nelems);
                        strncpy(_val[1], str[1], attr[1]->nelems);
                        printf("\"%s\" vs \"%s\"\n", _val[0], _val[1]);
                        free(_val[0]);
                        free(_val[1]);
                    }
                    else if (attr[0]->xtype == NC_BYTE)   PRINTF_VAR_DIFF(signed char)
                    else if (attr[0]->xtype == NC_UBYTE)  PRINTF_VAR_DIFF(unsigned char)
                    else if (attr[0]->xtype == NC_SHORT)  PRINTF_VAR_DIFF(short)
                    else if (attr[0]->xtype == NC_USHORT) PRINTF_VAR_DIFF(unsigned short)
                    else if (attr[0]->xtype == NC_INT)    PRINTF_VAR_DIFF(int)
                    else if (attr[0]->xtype == NC_UINT)   PRINTF_VAR_DIFF(unsigned int)
                    else if (attr[0]->xtype == NC_FLOAT)  PRINTF_VAR_DIFF(float)
                    else if (attr[0]->xtype == NC_DOUBLE) PRINTF_VAR_DIFF(double)
                    else if (attr[0]->xtype == NC_INT64)  PRINTF_VAR_DIFF(long long)
                    else if (attr[0]->xtype == NC_UINT64) PRINTF_VAR_DIFF(unsigned long long)
                }
                numHeadDIFF++;
            }
            else if (verbose)
                printf("\tSAME: attribute contents\n");
        }

        /* check global attributes defined in 2nd file but not in 1st file */
        for (i=0; i<nattrs[1]; i++) {
            attr[1] = ncp[1]->attrs.value[i];
            /* compare attribute name */
            k = i % nattrs[0];
            for (j=0; j<nattrs[0]; j++) {
                attr[0] = ncp[0]->attrs.value[k];
                if (strcmp(attr[1]->name, attr[0]->name))
                    k = (k + 1) % nattrs[0];
                else
                    break; /* loop j */
            }
            if (j == nattrs[0]) { /* not found in 1st file */
                if (!quiet)
                    printf("DIFF: global attribute \"%s\" defined in %s not found in %s\n",
                           attr[1]->name, argv[optind+1], argv[optind]);
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
            dim[0] = ncp[0]->dims.value[i];

            /* compare dimension name */
            k = i % ndims[1];
            for (j=0; j<ndims[1]; j++) {
                dim[1] = ncp[1]->dims.value[k];
                if (strcmp(dim[0]->name, dim[1]->name))
                    k = (k + 1) % ndims[1];
                else
                    break; /* loop j */
            }
            if (j == ndims[1]) { /* not found in 2nd file */
                if (!quiet)
                    printf("DIFF: dimension \"%s\" defined in %s not found in %s\n",
                           dim[1]->name, argv[optind+1], argv[optind]);
                numHeadDIFF++;
                continue; /* loop i */
            }

            /* compare dimension length */
            if (dim[0]->size != dim[1]->size) {
                if (!quiet)
                    printf("DIFF: dimension \"%s\" length (%lld) != (%lld)\n",
                           dim[0]->name, dim[0]->size, dim[1]->size);
                numHeadDIFF++;
            }
            else if (verbose)
                printf("\tSAME: dimension \"%s\" length (%lld)\n",
                       dim[0]->name, dim[0]->size);
        }

        /* check dimensions in 2nd file but not in 1st file */
        for (i=0; i<ndims[1]; i++) {
            dim[1] = ncp[1]->dims.value[i];
            /* compare dimension name */
            k = i % ndims[0];
            for (j=0; j<ndims[0]; j++) {
                dim[0] = ncp[0]->dims.value[k];
                if (strcmp(dim[1]->name, dim[0]->name))
                    k = (k + 1) % ndims[0];
                else
                    break; /* loop j */
            }
            if (j == ndims[0]) { /* not found in 1st file */
                if (!quiet)
                    printf("DIFF: dimension \"%s\" defined in %s not found in %s\n",
                           dim[1]->name, argv[optind+1], argv[optind]);
                numHeadDIFF++;
            }
        }

        /* Compare variables' metadata */
cmp_vars:
        if (nvars[0] > 0 && nvars[1] > 0) {
            if (verbose)
                printf("Variables:\n");
        } else
            goto fn_exit;

        /* check variables defined in 1st file and also in 2nd file */
        for (i=0; i<nvars[0]; i++) {
            var[0] = ncp[0]->vars.value[i];

            /* compare variable name */
            k = i % nvars[1];
            for (j=0; j<nvars[1]; j++) {
                var[1] = ncp[1]->vars.value[k];
                if (strcmp(var[0]->name, var[1]->name))
                    k = (k + 1) % nvars[1];
                else
                    break; /* loop j */
            }
            if (j == nvars[1]) { /* not found in 2nd file */
                if (!quiet)
                    printf("DIFF: variable \"%s\" defined in %s not found in %s\n",
                           var[0]->name, argv[optind], argv[optind+1]);
                numHeadDIFF++;
                numVarDIFF++;
                continue; /* loop i */
            }

            /* compare variable xtype */
            if (var[0]->xtype != var[1]->xtype) {
                if (!quiet)
                    printf("DIFF: variable \"%s\" data type (%s) != (%s)\n",
                           var[0]->name,
                           get_type(var[0]->xtype), get_type(var[1]->xtype));
                numHeadDIFF++;
            }
            else if (verbose) {
                printf("Variable \"%s\":\n", var[0]->name);
                printf("\tSAME: data type (%s)\n",get_type(var[0]->xtype));
            }

            /* compare variable ndims */
            if (var[0]->ndims != var[1]->ndims) {
                if (!quiet)
                    printf("DIFF: variable \"%s\" number of dimensions (%d) != (%d)\n",
                           var[0]->name, var[0]->ndims, var[1]->ndims);
                numHeadDIFF++;
            }
            else {
                if (verbose)
                    printf("\tSAME: number of dimensions (%d)\n", var[0]->ndims);

                /* compare variable's dimensionality */
                for (j=0; j<var[0]->ndims; j++) {
                    dim[0] = ncp[0]->dims.value[var[0]->dimids[j]];
                    dim[1] = ncp[1]->dims.value[var[1]->dimids[j]];

                    if (verbose)
                        printf("\tdimension %d:\n",j);

                    /* compare variable dimension j's name */
                    if (strcmp(dim[0]->name, dim[1]->name) != 0) {
                        if (!quiet)
                            printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's name (%s) != (%s)\n",
                                   var[0]->name, get_type(var[0]->xtype), j,
                                   dim[0]->name, dim[1]->name);
                        numHeadDIFF++;
                    }
                    else if (verbose)
                        printf("\t\tSAME: name (%s)\n", dim[0]->name);

                    /* compare variable dimension j's length */
                    if (dim[0]->size != dim[1]->size) {
                        if (!quiet)
                            printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's length (%lld) != (%lld)\n",
                                   var[0]->name, get_type(var[0]->xtype), j,
                                   dim[0]->size, dim[1]->size);
                        numHeadDIFF++;
                    }
                    else if (verbose)
                        printf("\t\tSAME: length (%lld)\n", dim[0]->size);
                }
            }

            /* compare variable's attributes */
            nattrs[0] = var[0]->attrs.ndefined;
            nattrs[1] = var[1]->attrs.ndefined;

            /* compare number of attributes of this variable */
            if (nattrs[0] != nattrs[1]) {
                if (!quiet)
                    printf("DIFF: variable \"%s\" number of attributes (%d) != (%d)\n",
                           var[0]->name, nattrs[0], nattrs[1]);
                numHeadDIFF++;
            }
            else if (verbose)
                printf("\tSAME: number of attributes (%d)\n", nattrs[0]);

            /* attributes in 1st file also appear in 2nd file */
            for (j=0; j<nattrs[0]; j++) {
                attr[0] = var[0]->attrs.value[j];

                /* find the variable attr with the same name from 2nd file */
                n = j % nattrs[1];
                for (m=0; m<nattrs[1]; m++) {
                    attr[1] = var[1]->attrs.value[n];
                    if (strcmp(attr[0]->name, attr[1]->name) != 0)
                        n = (n + 1) % nattrs[1];
                    else
                        break;
                }
                if (m == nattrs[1]) { /* not found in 2nd file */
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" defined in %s not found in %s\n",
                               var[0]->name, attr[0]->name, argv[optind], argv[optind+1]);
                    numHeadDIFF++;
                    continue; /* skip this attribute */
                }
                if (verbose)
                    printf("\tattribute \"%s\":\n", attr[0]->name);

                /* compare attribute xtype */
                if (attr[0]->xtype != attr[1]->xtype) {
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" data type (%s) != (%s)\n",
                               var[0]->name, attr[0]->name,
                               get_type(attr[0]->xtype), get_type(attr[1]->xtype));
                    numHeadDIFF++;
                    continue; /* skip this attribute */
                }
                else if (verbose)
                    printf("\t\tSAME: data type (%s)\n", get_type(attr[0]->xtype));

                /* compare attribute nelems */
                if (attr[0]->nelems != attr[1]->nelems) {
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" length (%lld) != (%lld)\n",
                               var[0]->name, attr[0]->name,
                               attr[0]->nelems, attr[1]->nelems);
                    numHeadDIFF++;
                    continue; /* skip this attribute */
                }
                else if (verbose)
                    printf("\t\tSAME: length (%lld)\n", attr[0]->nelems);

                /* compare attribute contents */
                nbytes = attr[0]->nelems * len_nctype(attr[0]->xtype);
                if ((isDiff = memcmp(attr[0]->xvalue, attr[1]->xvalue, nbytes)) != 0) {
                    char *str[2];
                    str[0] = (char*) attr[0]->xvalue;
                    str[1] = (char*) attr[1]->xvalue;
                    /* find the array index of first element in difference */
                    for (m=0; m<nbytes; m++)
                        if (str[0][m] != str[1][m])
                            break;
                    isDiff = m / len_nctype(attr[0]->xtype);
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" of type \"%s\" at element %d\n",
                               var[0]->name, attr[0]->name, get_type(attr[0]->xtype), isDiff);
                    numHeadDIFF++;
                }
                else if (verbose)
                    printf("\t\tSAME: attribute contents\n");
            }

            /* check attributes in 2nd file but not in 1st file */
            for (j=0; j<nattrs[1]; j++) {
                attr[1] = var[1]->attrs.value[j];

                /* find the variable attr with the same name from 1st file */
                n = j % nattrs[0];
                for (m=0; m<nattrs[0]; m++) {
                    attr[0] = var[0]->attrs.value[n];
                    if (strcmp(attr[0]->name, attr[1]->name) != 0)
                        n = (n + 1) % nattrs[0];
                    else
                        break;
                }
                if (m == nattrs[0]) { /* not found in 1st file */
                    if (!quiet)
                        printf("DIFF: variable \"%s\" attribute \"%s\" defined in %s not found in %s\n",
                               var[1]->name, attr[1]->name, argv[optind+1],argv[optind]);
                    numHeadDIFF++;
                }
            }
        }

        /* check variables defined in 2nd file but not in 1st file */
        for (i=0; i<nvars[1]; i++) {
            var[1] = ncp[1]->vars.value[i];
            /* compare variable name */
            k = i % nvars[0];
            for (j=0; j<nvars[0]; j++) {
                var[0] = ncp[0]->vars.value[k];
                if (strcmp(var[1]->name, var[0]->name))
                    k = (k + 1) % nvars[0];
                else
                    break; /* loop j */
            }
            if (j == nvars[0]) { /* not found in 1st file */
                if (!quiet)
                    printf("DIFF: variable \"%s\" defined in %s not found in %s\n",
                           var[1]->name, argv[optind+1], argv[optind]);
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
        cmp_nvars = nvars[0];
        var_list.nvars = nvars[0];
        var_list.names = (char**) calloc((size_t)cmp_nvars, sizeof(char*));
        if (!var_list.names) OOM_ERROR
        /* collect names of all variables from 1st file */
        for (i=0; i<cmp_nvars; i++)
            var_list.names[i] = ncp[0]->vars.value[i]->name;
    }
    if (verbose) printf("number of variables to be compared = %d\n",cmp_nvars);

    /* allocate read buffers */
    buf[0] = (void*) malloc(READ_CHUNK_SIZE);
    buf[1] = (void*) malloc(READ_CHUNK_SIZE);

    /* compare variables, one at a time */
    for (i=0; i<cmp_nvars; i++) {
        char *var_name;
        int varid[2], *dimids[2], isRecVar[2];
        long long remainLen, nChunks, varsize[2], offset[2];
        off_t seek_ret;
        nc_type xtype[2];
        NC_var  *var[2];

        /* find variable ID in 1st file corresponding to var_list.names[i] */
        for (j=0; j<nvars[0]; j++)
            if (strcmp(ncp[0]->vars.value[j]->name, var_list.names[i]) == 0)
                break;
        if (j == nvars[0]) {
            if (!check_header) {
                if (!quiet)
                    printf("WARN: variable \"%s\" s not found in %s\n",
                           var_list.names[i],argv[optind]);
                numVarDIFF++;
            }
            continue;
        } else
            varid[0] = j;

        /* find variable ID in 2nd file corresponding to var_list.names[i] */
        for (j=0; j<nvars[1]; j++)
            if (strcmp(ncp[1]->vars.value[j]->name, var_list.names[i]) == 0)
                break;
        if (j == nvars[1]) {
            if (!check_header) {
                if (!quiet)
                    printf("WARN: variable \"%s\" s not found in %s\n",
                           var_list.names[i],argv[optind+1]);
                numVarDIFF++;
            }
            continue;
        } else
            varid[1] = j;

        /* these below pointers are just short cuts */
        for (j=0; j<2; j++) {
                 var[j] = ncp[j]->vars.value[varid[j]];
              dimids[j] = var[j]->dimids;
               xtype[j] = var[j]->xtype;
               ndims[j] = var[j]->ndims;
            isRecVar[j] = 0;
            if (ndims[j] > 0 && dimids[j][0] == ncp[j]->dims.unlimited_id)
                isRecVar[j] = 1;
        }
        var_name = var_list.names[i];

	/* Header comparison may have been skipped. Even if file headers have
	 * been compared, we still need to compare variable's xtype and
	 * dimensions to skip variables when their structures are different.
	 */

	/* compare variable's NC data type */
        if (xtype[0] != xtype[1]) {
            if (!check_header) { /* if header has not been checked */
                if (!quiet)
                    printf("DIFF: variable \"%s\" data type (%s) != (%s)\n",
                           var_name,get_type(xtype[0]),get_type(xtype[1]));
                numHeadDIFF++;
                numVarDIFF++;
            }
            continue; /* skip this variable */
        }
        else if (!check_header && verbose) {
            printf("Variable \"%s\":\n",var_name);
            printf("\tSAME: data type (%s)\n",get_type(xtype[0]));
        }

        /* compare variable's number of dimensions */
        if (ndims[0] != ndims[1]) {
            if (!check_header) { /* if header has not been checked */
                if (!quiet)
                    printf("DIFF: variable \"%s\" number of dimensions (%d) != (%d)\n",
                           var_name,ndims[0],ndims[1]);
                numHeadDIFF++;
                numVarDIFF++;
            }
            continue; /* skip this variable */
        }
        else if (!check_header && verbose)
            printf("\tSAME: number of dimensions (%d)\n",ndims[0]);

        /* compare variable's dimension sizes, not dimension's names */
        for (j=0; j<ndims[0]; j++) {
            long long dimlen[2];
            dimlen[0] = ncp[0]->dims.value[dimids[0][j]]->size;
            dimlen[1] = ncp[1]->dims.value[dimids[1][j]]->size;
            if (!check_header && verbose)
                printf("\tDimension %d:\n",j);
            if (dimlen[0] != dimlen[1]) {
                if (!check_header) { /* if header has not been checked */
                    if (!quiet)
                        printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's length (%lld) != (%lld)\n",
                               var_name,get_type(xtype[0]),j,dimlen[0],dimlen[1]);
                    numHeadDIFF++;
                    numVarDIFF++;
                }
                break; /* skip this variable */
            }
            else if (!check_header && verbose)
                printf("\t\tSAME: length (%lld)\n",dimlen[0]);
        }
        if (j != ndims[0])
            continue; /* diff found, skip this variable */

	/* variable total size should be the same */
	assert(var[0]->len == var[1]->len);

	/* calculate variable's size in bytes. Cannot use var[i]->len, as
         * it is upward aligned to a 4-byte boundary.
         */
        for (j=0; j<2; j++) {
	    varsize[j] = var[j]->xsz;
            for (k=0; k<var[j]->ndims; k++) {
                if (var[j]->shape[k] != NC_UNLIMITED)
                    varsize[j] *= var[j]->shape[k];
            }
        }
        assert(varsize[0] == varsize[1]);

        worst = -1;
        chunk_off = 0; /* array index offset for this chunk */
        numrecs = (isRecVar[0]) ? ncp[0]->numrecs : 1;
        /* numrecs = 0 for fixed-size variables */

        for (r=0; r<numrecs; r++) {
            for (j=0; j<2; j++) {
                long long recsize = (isRecVar[j]) ? ncp[j]->recsize : 0;
                /* starting file offsets can be different between 2 files */
                offset[j] = var[j]->begin + recsize * r;

                /* lseek to starting file offset of the next record */
                seek_ret = lseek(fd[j], (off_t)offset[j], SEEK_SET);
                if (seek_ret < 0) {
                    fprintf(stderr, "Error on lseek offset at %lld file %s (%s)\n",
                            offset[j], argv[optind+i], strerror(errno));
                    goto fn_exit; /* fatal file error */
                }
            }

	    /* for record variable, varsize is the size of one record of this
             * variable. For fixed-size variable, it is the whole size.
             */
	    remainLen = varsize[0];
            nChunks = remainLen / READ_CHUNK_SIZE;
            if (remainLen % READ_CHUNK_SIZE) nChunks++;

            /* compare the variable contents, one chunk at a time */
            for (k=0; k<nChunks; k++) {
                ssize_t rdLen[2];
                size_t  rSize = MIN(remainLen, READ_CHUNK_SIZE);
                /* read next chunk */
                rdLen[0] = read(fd[0], buf[0], rSize);
                rdLen[1] = read(fd[1], buf[1], rSize);
                if (rdLen[0] != rdLen[1]) {
                    if (!quiet)
                        printf("DIFF: variable \"%s\" of type \"%s\" record %lld read size (%zd) != (%zd)\n",
                                var_name, get_type(xtype[0]), r, rdLen[0], rdLen[1]);
                    numVarDIFF++;
                    break; /* loop k */
                }
                remainLen -= rdLen[0];

                /* compare contents of chunks */
                if (check_tolerance && worst == -1) {
                         if (xtype[0] == NC_CHAR)   CHECK_VAR_DIFF(char)
                    else if (xtype[0] == NC_BYTE)   CHECK_VAR_DIFF(signed char)
                    else if (xtype[0] == NC_UBYTE)  CHECK_VAR_DIFF_UNSIGNED(unsigned char)
                    else if (xtype[0] == NC_SHORT)  CHECK_VAR_DIFF(short)
                    else if (xtype[0] == NC_USHORT) CHECK_VAR_DIFF_UNSIGNED(unsigned short)
                    else if (xtype[0] == NC_INT)    CHECK_VAR_DIFF(int)
                    else if (xtype[0] == NC_UINT)   CHECK_VAR_DIFF_UNSIGNED(unsigned int)
                    else if (xtype[0] == NC_FLOAT)  CHECK_VAR_DIFF(float)
                    else if (xtype[0] == NC_DOUBLE) CHECK_VAR_DIFF(double)
                    else if (xtype[0] == NC_INT64)  CHECK_VAR_DIFF(long long)
                    else if (xtype[0] == NC_UINT64) CHECK_VAR_DIFF_UNSIGNED(unsigned long long)
                    if (worst != -1) {
                        r = numrecs; /* break both loops k and r */
                        break;
                    }
                }
                else { /* if (check_tolerance) */
                    if (memcmp(buf[0], buf[1], rdLen[0]) == 0) continue;

                    /* difference found */
                    if (!quiet) {
                        char *str[2];
                        str[0] = (char*) buf[0];
                        str[1] = (char*) buf[1];
                        if (ndims[0] == 0) { /* scalar variable */
                            printf("DIFF: scalar variable \"%s\" of type \"%s\" of value ",
                                   var_name, get_type(xtype[0]));
                        }
                        else {
                            int _i;
                            long long pos, *diffStart;

                            /* find the array index of first element in difference */
                            for (m=0; m<rdLen[0]; m++)
                                if (str[0][m] != str[1][m])
                                    break;
                            pos = ((long long)m + k * nChunks) / var[0]->xsz;

                            diffStart = (long long*) malloc(var[0]->ndims * sizeof(long long));
                            for (_i=var[0]->ndims-1; _i>=0; _i--) {
                                if (isRecVar[0] && _i == 0) {
                                    diffStart[_i] = r; /* r-th record */
                                } else {
                                    diffStart[_i] = pos % var[0]->shape[_i];
                                    pos /= var[0]->shape[_i];
                                }
                            }
                            printf("DIFF: variable \"%s\" of type \"%s\" at element [%lld",
                                   var_name, get_type(xtype[0]), diffStart[0]);
                            for (_i=1; _i<var[0]->ndims; _i++)
                                printf(", %lld", diffStart[_i]);
                            printf("] of value ");
                            free(diffStart);

                            pos = m / var[0]->xsz;
                            pos *= var[0]->xsz;
                            str[0] += pos;
                            str[1] += pos;
                        }
                        if (xtype[0] == NC_CHAR) printf("%s vs %s\n", str[0], str[1]);
                        else if (xtype[0] == NC_BYTE)   PRINTF_VAR_DIFF(signed char)
                        else if (xtype[0] == NC_UBYTE)  PRINTF_VAR_DIFF(unsigned char)
                        else if (xtype[0] == NC_SHORT)  PRINTF_VAR_DIFF(short)
                        else if (xtype[0] == NC_USHORT) PRINTF_VAR_DIFF(unsigned short)
                        else if (xtype[0] == NC_INT)    PRINTF_VAR_DIFF(int)
                        else if (xtype[0] == NC_UINT)   PRINTF_VAR_DIFF(unsigned int)
                        else if (xtype[0] == NC_FLOAT)  PRINTF_VAR_DIFF(float)
                        else if (xtype[0] == NC_DOUBLE) PRINTF_VAR_DIFF(double)
                        else if (xtype[0] == NC_INT64)  PRINTF_VAR_DIFF(long long)
                        else if (xtype[0] == NC_UINT64) PRINTF_VAR_DIFF(unsigned long long)
                    }
                    worst = 0;
                    numVarDIFF++;
                    r = numrecs; /* break both loops k and r */
                    break;
                }
            } /* loop k */
        } /* loop r */

        if (!check_tolerance && worst == 0)
            continue; /* go to next variable */

        if (worst == -1) {
            if (verbose)
                printf("\tSAME: variable \"%s\" contents\n",var_name);
            continue; /* go to next variable */
        }

        /* diff is found when check_tolerance = 1 */
        if (ndims[0] == 0) { /* scalar variable */
            printf("DIFF (tolerance): scalar variable \"%s\" of type \"%s\" of value %g vs %g (difference = %e)\n",
                   var_name, get_type(xtype[0]), worst1, worst2, worst1-worst2);
        } else {
            int _i;
            long long pos, *diffStart;
            diffStart = (long long*) malloc(ndims[0] * sizeof(long long));
            if (worst != -1) pos = worst;
            for (_i=ndims[0]-1; _i>=0; _i--) {
                long long dimLen = var[0]->shape[_i];
                if (_i == 0 && isRecVar[0]) dimLen = numrecs;
                diffStart[_i] = pos % dimLen;
                pos /= dimLen;
            }
            if (worst == -1)
                printf("DIFF: variable \"%s\" of type \"%s\" first at element [%lld",
                       var_name, get_type(xtype[0]), diffStart[0]);
            else
                printf("DIFF (tolerance): variable \"%s\" of type \"%s\" at element [%lld",
                       var_name, get_type(xtype[0]), diffStart[0]);
            for (_i=1; _i<ndims[0]; _i++)
                printf(", %lld", diffStart[_i]);
            if (worst == -1)
                printf("]\n");
            else
                printf("] of value %g vs %g (difference = %e)\n", worst1, worst2, worst1-worst2);
            free(diffStart);
        }
        numVarDIFF++;
    } /* loop i */

fn_exit:
    /* close files and free up the memory previously allocated */
    for (i=0; i<2; i++) {
        if (-1 == close(fd[i]))
            fprintf(stderr, "Error on close file %s (%s)\n",
                    argv[optind+i], strerror(errno));
        if (buf[i] != NULL) free(buf[i]);
        free_NC_dimarray(&ncp[i]->dims);
        free_NC_attrarray(&ncp[i]->attrs);
        free_NC_vararray(&ncp[i]->vars);
        free(ncp[i]);
    }
    if (var_list.nvars)
        free(var_list.names);

    /* summary of the difference */
    if (!quiet) {
        if (check_header) {
            if (numHeadDIFF == 0)
                printf("Headers of two files are the same\n");
            else
                printf("Number of differences in header: %lld\n",numHeadDIFF);
        }
        if (check_variable_list) {
            if (numVarDIFF == 0)
                printf("Compared variable(s) are the same\n");
            else
                printf("Compared variables(s) has %lld differences\n",numVarDIFF);
        }
        if (check_entire_file) {
            if (numVarDIFF == 0)
                printf("All variables of two files are the same\n");
            else
                printf("Number of differences in variables: %lld\n",numVarDIFF);
        }
    }

    numDIFF = numVarDIFF + numHeadDIFF;

    exit((numDIFF == 0) ? EXIT_SUCCESS : EXIT_FAILURE);
}
