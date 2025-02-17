/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* This file contains the implementation of CDL header file APIs, which are
 * public to PnetCDF users, i.e. visible in file pnetcdf.h.
 *    cdl_hdr_open()       opens and parses the CDL file's header
 *    cdl_hdr_inq_format() returns file format version
 *    cdl_hdr_inq_ndims()  returns number of dimensions defined in CDL file
 *    cdl_hdr_inq_dim()    returns metadata of a dimension
 *    cdl_hdr_inq_nvars()  returns number of variables
 *    cdl_hdr_inq_var()    returns metadata of a variable defined in CDL file
 *    cdl_hdr_inq_nattrs() returns number of attributes of a given variable
 *    cdl_hdr_inq_attr()   returns metadata of an attribute
 *    cdl_hdr_close()      closes the CDL file
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#ifndef TEST_RUN
#include <mpi.h>
#include <pnetcdf.h>

#else
typedef int nc_type;
typedef long long MPI_Offset;

/* below define C macros are duplicated from pnetcdf.h */
#define NC_GLOBAL -1
#define NC_UNLIMITED 0L

#define NC_NAT          0       /**< Not A Type */
#define NC_BYTE         1       /**< signed 1 byte integer */
#define NC_CHAR         2       /**< ISO/ASCII character */
#define NC_SHORT        3       /**< signed 2 byte integer */
#define NC_INT          4       /**< signed 4 byte integer */
#define NC_LONG         NC_INT  /**< \deprecated required for backward compatibility. */
#define NC_FLOAT        5       /**< single precision floating point number */
#define NC_DOUBLE       6       /**< double precision floating point number */
#define NC_UBYTE        7       /**< unsigned 1 byte int */
#define NC_USHORT       8       /**< unsigned 2-byte int */
#define NC_UINT         9       /**< unsigned 4-byte int */
#define NC_INT64        10      /**< signed 8-byte int */
#define NC_UINT64       11      /**< unsigned 8-byte int */
#define NC_STRING       12      /**< string */

#define NC_MAX_NFILES   1024

#define NC_NOERR    0      /**< No Error */
#define NC_EBADID   (-33)
#define NC_ENFILE   (-34)  /**< Too many netcdfs open */
#define NC_EBADDIM  (-46)  /**< Invalid dimension id or name */
#define NC_ENOTNC   (-51)  /**< Not a netcdf file */
#define NC_ENOMEM   (-61)  /**< Memory allocation (malloc) failure */
#define NC_EFILE    (-204) /**< Unknown error in file operation */
#define NC_ENOENT   (-220) /**< File does not exist */

#define ERR { \
    if (err != NC_NOERR) { \
        char *err_str; \
        switch(err) { \
            case NC_EBADID:  err_str = "Bad ID"; \
                             break; \
            case NC_ENFILE:  err_str = "Too many CDL file open"; \
                             break; \
            case NC_EBADDIM: err_str = "Invalid dimension id or name"; \
                             break; \
            case NC_ENOTNC:  err_str = "Not a CDL file"; \
                             break; \
            case NC_EFILE:   err_str = "Unknown error in file operation"; \
                             break; \
            case NC_ENOENT:  err_str = "File does not exist"; \
                             break; \
            default: err_str = "Other unknown error"; \
        } \
        printf("Error in %s and %d : %s\n", __FILE__,__LINE__, err_str); \
        err = 1; \
        goto err_out; \
    } \
}
#endif

#define LINE_SIZE 1024

#define ERR_FORMAT(msg) { \
    printf("Error in %s at %d: input file format at %s\n", \
           __func__,__LINE__,msg); \
    err = NC_ENOTNC; \
    goto err_out; \
}

#define DIM_ARRAY_GROWBY   32
#define VAR_ARRAY_GROWBY   128
#define ATTR_ARRAY_GROWBY  8
#define GATTR_ARRAY_GROWBY 64

typedef struct {
    char       *name;
    MPI_Offset  size;
} CDL_dim;

typedef struct {
    int      nelems;      /* number of defined dimensions */
    CDL_dim *value;
} CDL_dimarray;

typedef struct {
    char       *name;    /* name of the attributes */
    nc_type     xtype;   /* external NC data type of the attribute */
    MPI_Offset  nelems;  /* number of attribute elements */
    void       *value;   /* attribute contents */
} CDL_attr;

typedef struct {
    int       nelems;  /* number of defined attributes */
    CDL_attr *value;
} CDL_attrarray;

typedef struct {
    char          *name;    /* name of the variable */
    nc_type        xtype;   /* variable's external NC data type */
    int            ndims;   /* number of dimensions */
    int           *dimids;  /* [ndims] array of dimension IDs */
    CDL_attrarray  attrs;   /* attribute array */
} CDL_var;

typedef struct {
    int      nelems;    /* number of defined variables */
    CDL_var *value;
} CDL_vararray;

typedef struct {
    int           format;   /* 1, 2, or 5 corresponding to CDF-1, 2, or 5 */
    CDL_dimarray  dims;     /* dimensions defined */
    CDL_attrarray attrs;    /* global attributes defined */
    CDL_vararray  vars;     /* variables defined */
} CDL_header;

static CDL_header *cdl_filelist[NC_MAX_NFILES];
static int cdl_nfiles;

static
int c_to_xtype(char *ctype)
{
    if (!strcmp(ctype, "byte"))   return NC_BYTE;
    if (!strcmp(ctype, "char"))   return NC_CHAR;
    if (!strcmp(ctype, "short"))  return NC_SHORT;
    if (!strcmp(ctype, "int"))    return NC_INT;
    if (!strcmp(ctype, "float"))  return NC_FLOAT;
    if (!strcmp(ctype, "double")) return NC_DOUBLE;
    if (!strcmp(ctype, "ubyte"))  return NC_UBYTE;
    if (!strcmp(ctype, "ushort")) return NC_USHORT;
    if (!strcmp(ctype, "uint"))   return NC_UINT;
    if (!strcmp(ctype, "int64"))  return NC_INT64;
    if (!strcmp(ctype, "uint64")) return NC_UINT64;
    return NC_NAT;
}

static
int get_dimid(CDL_header *header,
              char       *name,
              int        *dimid)
{
    int i;
    CDL_dim *dimp = header->dims.value;

    for (i=0; i<header->dims.nelems; i++)
        if (!strcmp(dimp[i].name, name)) {
            *dimid = i;
            return NC_NOERR;
        }

    printf("Error in %s at %d: failed to find dim ID for %s\n",
           __func__,__LINE__,name);
    return NC_EBADDIM;
}

static
int parse_signature(char **bptr)
{
    char *line;
    int err=NC_NOERR;

    /* first line is signature "netcdf" and file name */
    line = strtok(*bptr, "\n");
    if (line[0] == '}') return 1; /* end of file */

    *bptr += strlen(line) + 1;
    while (**bptr == '\n') (*bptr)++;
    if (line != NULL) {
        if (strlen(line) == 0) ERR_FORMAT("signature line")
        char *signature = strtok(line, " ");
        if (signature == NULL || strcmp(signature, "netcdf"))
            ERR_FORMAT("signature")
#ifdef TEST_RUN
        printf("LINE %d signature: %s\n",__LINE__, signature);
#endif
        char *filename = strtok(NULL, " ");
        if (filename == NULL)
            ERR_FORMAT("input file name")
#ifdef TEST_RUN
        printf("LINE %d input file name: %s\n",__LINE__, filename);
#endif
    }
err_out:
    return err;
}

static
int parse_format(char       **bptr,
                 CDL_header  *header,
                 int         *dim_sec,
                 int         *var_sec,
                 int         *gattr_sec)
{
    char *line, *key, *val;
    int err=NC_NOERR;

    /* 2nd line is a comment containing file format version */
    header->format = 1;

    while ((line = strtok(*bptr, "\n")) != NULL) {
        if (line[0] == '}') return 1; /* end of file */

        *bptr += strlen(line) + 1;
        while (**bptr == '\n') (*bptr)++;
        key = strtok(line, " \t\n");
        if (key == NULL) continue; /* blank line */
        if (!strncmp(key, "data:", 5)) return 1; /* data section */
        if (!strcmp(key, "//")) { /* comment line */
            val = strtok(NULL, "\n");
            if (val == NULL) continue;
            char *fmt_str="file format: CDF-";
            if (!strncmp(val, fmt_str, strlen(fmt_str))) {
                val += strlen(fmt_str);
                val[1] = '\0';
                header->format = atoi(val);
                if (*val != '1' && *val != '2' && *val != '5')
                    return NC_ENOTNC;
#ifdef TEST_RUN
                printf("LINE %d CDF file format: CDF-%d\n",
                    __LINE__,header->format);
#endif
            }
            else if (!strcmp(key+3, "global attributes:"))
                *gattr_sec = 1;
        }
        else {
            if (!strcmp(key+3, "global attributes:"))
                *gattr_sec = 1;
            else if (!strcmp(key, "dimensions:"))
                *dim_sec = 1;
            else if (!strcmp(key, "variables:"))
                *var_sec = 1;
            break;
        }
    }
    return err;
}

static
int parse_dims(char         **bptr,
               CDL_dimarray  *dims,
               int           *var_sec,
               int           *gattr_sec)
{
    char *line, *key, *val;
    int err=0;

    while ((line = strtok(*bptr, "\n")) != NULL) {
        if (line[0] == '}') return 1; /* end of file */

        *bptr += strlen(line) + 1;
        while (**bptr == '\n') (*bptr)++;
        key = strtok(line, "\n");
        if (key == NULL || !strcmp(key, "//")) continue;
        if (!strncmp(key, "data:", 5)) return 1; /* data section */
        if (!strcmp(key, "variables:")) {
            *var_sec = 1;
            break;
        }
        if (!strcmp(key, "// global attributes:")) {
            *gattr_sec = 1;
            break;
        }
        key = strtok(key, " \t");
        val = strtok(NULL, " =");

        if (dims->nelems % DIM_ARRAY_GROWBY == 0) {
            size_t len = dims->nelems + DIM_ARRAY_GROWBY;
            dims->value = (CDL_dim*) realloc(dims->value,
                                             sizeof(CDL_dim) * len);
            if (dims->value == NULL) {
                printf("Error in %s at %d: fail to realloc of size %zd (%s)\n",
                       __func__,__LINE__,sizeof(CDL_dim)*len,strerror(errno));
                return NC_ENOMEM;
            }
        }
        CDL_dim *dimp = &dims->value[dims->nelems];
        dimp->name = strdup(key);
        if (!strcmp(val, "UNLIMITED"))
            dimp->size = NC_UNLIMITED;
        else
            dimp->size = atoi(val);
        dims->nelems++;

#ifdef TEST_RUN
        printf("LINE %d DIM name %s val %s\n",__LINE__, key, val);
#endif
        if (key[0] == '}') return 1;
        if (!strncmp(key, "data:", 5)) return 1; /* data section */
    }

    return err;
}

#define PARSE_ATTR(itype, tail, conv) { \
    itype *buf = (itype*) attrp->value; \
    attrp->nelems=0; \
    key = strtok(str, ", ;"); \
    while (key != NULL) { \
        if (tail) key[strlen(key) - tail] = '\0'; \
        buf[attrp->nelems++] = conv(key); \
        key = strtok(NULL, ", ;"); \
    } \
}

static
int parse_attr_value(CDL_attr *attrp,
                     char     *val)
{
    int err=0;

    if (val[1] == '\"') {
        attrp->xtype = NC_CHAR;
        char *tail = strrchr(val, '\"');
        *tail = '\0';
        attrp->value = strdup(val+2);
#ifdef TEST_RUN
        printf("LINE %d ATTR name %s type NC_CHAR %s\n",__LINE__,
            attrp->name, (char*)attrp->value);
#endif
        attrp->nelems = strlen(attrp->value);
    }
    else { /* not text */
        int xtype, nelems=0;
        char *str=strdup(val);
        char *key = strtok(str, ", ;");
        xtype = NC_INT;
        while (key != NULL) {
            nelems++;
            size_t key_len = strlen(key);
            if (key[key_len-1] == 'b') xtype = NC_BYTE;
            else if (key[key_len-1] == 's') xtype = NC_SHORT;
            else if (key[key_len-1] == 'f') xtype = NC_FLOAT;
            else if (key[key_len-1] == 'd') xtype = NC_DOUBLE;
            else if (key[key_len-1] == 'U') xtype = NC_UINT;
            else if (key[key_len-2] == 'U') {
                if (key[key_len-1] == 'B')      xtype = NC_UBYTE;
                else if (key[key_len-1] == 'S') xtype = NC_USHORT;
            }
            else if (key[key_len-1] == 'L') {
                if (key[key_len-3] == 'U') xtype = NC_UINT64;
                else                       xtype = NC_INT64;
            }
            key = strtok(NULL, ", ;");
        }
        attrp->xtype = xtype;
#ifdef TEST_RUN
        printf("LINE %d ATTR type %d name %s nelems=%d\n",
            __LINE__, xtype, attrp->name, nelems);
#endif

        strcpy(str, val);

        attrp->value = (void*) malloc(sizeof(double) * nelems);
        if (attrp->value == NULL) {
            printf("Error in %s at %d: fail to malloc of size %zd (%s)\n",
                    __func__,__LINE__,sizeof(double)*nelems,strerror(errno));
            return NC_ENOMEM;
        }

             if (xtype == NC_INT)    PARSE_ATTR(int, 0, atoi)
        else if (xtype == NC_BYTE)   PARSE_ATTR(signed char, 1, atoi)
        else if (xtype == NC_SHORT)  PARSE_ATTR(short, 1, atoi)
        else if (xtype == NC_FLOAT)  PARSE_ATTR(float, 1, atof)
        else if (xtype == NC_DOUBLE) PARSE_ATTR(double, 1, atof)
        else if (xtype == NC_UBYTE)  PARSE_ATTR(unsigned char, 2, atoi)
        else if (xtype == NC_USHORT) PARSE_ATTR(unsigned short, 2, atoi)
        else if (xtype == NC_UINT)   PARSE_ATTR(unsigned int, 1, atoi)
        else if (xtype == NC_INT64)  PARSE_ATTR(long long, 2, atoll)
        else if (xtype == NC_UINT64) PARSE_ATTR(unsigned long long, 2, atoll)

        free(str);
    }
    return err;
}

static
int parse_attr(char          **bptr,
               char           *var_name,
               CDL_attrarray  *attrs)
{
    char *line, line_cpy[LINE_SIZE], *key, *val;
    int err=0;
    size_t prefix_len;

    /* parse attributes */
    while ((line = strtok(*bptr, "\n")) != NULL) {
        if (line[0] == '}') return 1; /* end of file */

        *bptr += strlen(line) + 1;
        while (**bptr == '\n') (*bptr)++;

        strncpy(line_cpy, line, LINE_SIZE-1);

        prefix_len = strspn(line_cpy, " \t");
        key = line_cpy + prefix_len;
        if (*var_name != '\0' && !strncmp(key, "// global attributes:", 21)) {
            *bptr = line;
            break;
        }
        if (!strncmp(key, "//", 2)) {
            continue;
        }
        if (!strncmp(key, "data:", 5)) return 1; /* data section */

        char *name;
        if (key[0] == ':') { /* global attribute */
            name = "";
            val = key + 1;
        }
        else {
            char *colon = strchr(key, ':');
            if (colon == NULL) { /* not attribute */
                *bptr = line;
                break;
            }
            *colon = '\0';
            name = key;
            if (strcmp(name, var_name)) ERR_FORMAT("attribute name");
            val = colon + 1;
        }

        /* check if it is an attribute of this variable */
        if (strcmp(name, var_name)) {
            *bptr = line;
            break;
        }

        if (*var_name == '\0') { /* global attributes */
            if (attrs->nelems % GATTR_ARRAY_GROWBY == 0) {
                size_t len = attrs->nelems + GATTR_ARRAY_GROWBY;
                attrs->value = (CDL_attr*) realloc(attrs->value,
                                                   sizeof(CDL_attr) * len);
                if (attrs->value == NULL) {
                    printf("Error in %s at %d: fail to realloc of size %zd (%s)\n",
                           __func__,__LINE__,sizeof(CDL_attr)*len,strerror(errno));
                    return NC_ENOMEM;
                }
            }
        }
        else {
            if (attrs->nelems % ATTR_ARRAY_GROWBY == 0) {
                size_t len = attrs->nelems + ATTR_ARRAY_GROWBY;
                attrs->value = (CDL_attr*) realloc(attrs->value,
                                                   sizeof(CDL_attr) * len);
                if (attrs->value == NULL) {
                    printf("Error in %s at %d: fail to realloc of size %zd (%s)\n",
                        __func__,__LINE__,sizeof(CDL_attr)*len,strerror(errno));
                    return NC_ENOMEM;
                }
            }
        }
        CDL_attr *attrp = &attrs->value[attrs->nelems];

        char *next;

        next = strchr(val, ' ');
        if (next == NULL) ERR_FORMAT("attribute name");
        *next = '\0';

        char *attr_name = val;
        attrp->name = strdup(attr_name);

        next = strchr(next+1, '=');
        val = next + 1;

        attrs->nelems++;
        parse_attr_value(attrp, val);
    }
err_out:
    return err;
}

static
int parse_vars(char       **bptr,
               CDL_header  *header,
               int         *gattr_sec)
{
    char *type, *var_name;
    char *line, orig_line[LINE_SIZE], *key, *val;
    int err=0;
    size_t prefix_len;
    CDL_vararray *vars = &header->vars;

    while ((line = strtok(*bptr, "\n")) != NULL) {
        if (line[0] == '}') return 1; /* end of file */

        *bptr += strlen(line) + 1;
        while (**bptr == '\n') (*bptr)++;

        strncpy(orig_line, line, LINE_SIZE-1);

        prefix_len = strspn(line, " \t");
        line += prefix_len;
        if (!strncmp(line, "// global attributes:", 21)) {
            *bptr = line;
            *gattr_sec = 1;
            break;
        }
        if (!strncmp(line, "//", 2)) {
            continue;
        }
        if (!strncmp(line, "data:", 5)) return 1; /* data section */

        type = strtok(line, " \t\n");

        if (vars->nelems % VAR_ARRAY_GROWBY == 0) {
            size_t len = vars->nelems + VAR_ARRAY_GROWBY;
            vars->value = (CDL_var*) realloc(vars->value,
                                             sizeof(CDL_var) * len);
            if (vars->value == NULL) {
                printf("Error in %s at %d: fail to realloc of size %zd (%s)\n",
                        __func__,__LINE__,sizeof(CDL_var)*len,strerror(errno));
                return NC_ENOMEM;
            }
        }
        CDL_var *varp = &vars->value[vars->nelems];
        var_name = strtok(NULL, "( ");
        varp->name = strdup(var_name);
        varp->xtype = c_to_xtype(type);
        varp->ndims = 0;
        varp->attrs.nelems = 0;
        varp->attrs.value = NULL;

        val = strtok(NULL, ");");
        if (val == NULL) {
            varp->dimids = NULL;
#ifdef TEST_RUN
            printf("LINE %d VAR type %s name %s IS SCALAR\n",
                __LINE__, type, varp->name);
#endif
        }
        else {
            char *str=strdup(val);
            key = strtok(str, ", ");
            while (key != NULL) {
                varp->ndims++;
                key = strtok(NULL, ", ");
            }
#ifdef TEST_RUN
            printf("LINE %d VAR type %s name %s val %s ndims=%d\n",
                __LINE__, type, varp->name, val, varp->ndims);
#endif

            varp->dimids = (int*) malloc(sizeof(int) * varp->ndims);
            if (varp->dimids == NULL) {
                printf("Error in %s at %d: fail to malloc of size %zd (%s)\n",
                        __func__,__LINE__,sizeof(int)*varp->ndims,strerror(errno));
                return NC_ENOMEM;
            }
            strcpy(str, val);
            key = strtok(str, ", ");
            int i = 0;
            while (key != NULL) {
                int dimid;
                err = get_dimid(header, key, &dimid);
                if (err != NC_NOERR) {
                    free(str);
                    free(varp->name);
                    varp->name = NULL;
                    free(varp->dimids);
                    varp->dimids = NULL;
                    return err;
                }
                varp->dimids[i++] = dimid;
                key = strtok(NULL, ", ");
            }
            free(str);
        }
        vars->nelems++;

        err = parse_attr(bptr, varp->name, &varp->attrs);
        if (err)
            goto err_out;
    }
err_out:
    return err;
}

static
int hdr_free(CDL_header *header)
{
    int i, j, err=NC_NOERR;

    /* free dimensions */
    for (i=0; i<header->dims.nelems; i++)
        free(header->dims.value[i].name);

    if (header->dims.value != NULL)
        free(header->dims.value);

    /* free global attributes */
    for (i=0; i<header->attrs.nelems; i++) {
        free(header->attrs.value[i].name);
        free(header->attrs.value[i].value);
    }

    if (header->attrs.value != NULL)
        free(header->attrs.value);

    /* free variables */
    for (i=0; i<header->vars.nelems; i++) {
        CDL_var *varp = &header->vars.value[i];

        if (varp->name != NULL) free(varp->name);
        if (varp->dimids != NULL) free(varp->dimids);

        /* free attributes */
        for (j=0; j<varp->attrs.nelems; j++) {
            free(varp->attrs.value[j].name);
            free(varp->attrs.value[j].value);
        }
        if (varp->attrs.value != NULL)
            free(varp->attrs.value);
    }

    if (header->vars.value != NULL)
        free(header->vars.value);

    free(header);

    return err;
}

int cdl_hdr_close(int hid)
{
    int err=NC_NOERR;
    CDL_header *header;

    /* check if hid is valid */
    if (hid >= NC_MAX_NFILES || cdl_filelist[hid] == NULL)
        return NC_EBADID;

    header = cdl_filelist[hid];
    err = hdr_free(header);

    cdl_filelist[hid] = NULL;

    while (cdl_nfiles > 0 && cdl_filelist[cdl_nfiles-1] == NULL)
        cdl_nfiles--;

    return err;
}

/*----< cdl_hdr_open() >-----------------------------------------------------*/
/* Reads a CDL header file, parses it, and stores the metadata in an internal
 * buffer. Returns an ID for future inquery APIs to retrive the metadata.
 */
int cdl_hdr_open(const char *filename,
                 int        *hid)
{
    FILE *fptr;
    char *fbuf, *bptr;
    int err, dim_sec, var_sec, gattr_sec;
    size_t rlen;
    CDL_header *header=NULL;

    if (cdl_nfiles + 1 > NC_MAX_NFILES)
        return NC_ENFILE;

    fptr = fopen(filename, "r");
    if (fptr == NULL) {
        printf("Error in %s at %d: fail to open file %s (%s)\n",
               __func__,__LINE__,filename,strerror(errno));
        return NC_ENOENT;
    }
    err = fseek(fptr, 0, SEEK_END);
    if (err < 0) {
        printf("Error in %s at %d: fail to fseek SEEK_END file %s (%s)\n",
               __func__,__LINE__,filename,strerror(errno));
        return NC_EFILE;
    }

    /* read the entire file into a buffer */
    long file_size = ftell(fptr);
    fbuf = (char *) malloc(file_size);
    if (fbuf == NULL) {
        printf("Error in %s at %d: fail to malloc of size %zd (%s)\n",
               __func__,__LINE__,file_size,strerror(errno));
        return NC_ENOMEM;
    }
    err = fseek(fptr, 0, SEEK_SET);
    if (err < 0) {
        printf("Error in %s at %d: fail to fseek SEEK_SET file %s (%s)\n",
               __func__,__LINE__,filename,strerror(errno));
        return NC_EFILE;
    }
    rlen = fread(fbuf, 1, file_size, fptr);
    if (file_size > 0 && rlen == 0) {
        printf("Error in %s at %d: fail to fread file %s (%s)\n",
               __func__,__LINE__,filename,strerror(errno));
        return NC_EFILE;
    }
    bptr = fbuf;
    err = fclose(fptr);
    if (err == EOF) {
        printf("Error in %s at %d: fail to fclose file %s (%s)\n",
               __func__,__LINE__,filename,strerror(errno));
        return NC_EFILE;
    }

    /* check netCDF file signature */
    err = parse_signature(&bptr);
    if (err != NC_NOERR) goto err_out;

    header = (CDL_header*) calloc(1, sizeof(CDL_header));
    if (header == NULL) {
        printf("Error in %s at %d: fail to calloc of size %zd (%s)\n",
               __func__,__LINE__,sizeof(CDL_header),strerror(errno));
        return NC_ENOMEM;
    }

    dim_sec = var_sec = gattr_sec = 0;

    /* check netCDF file format */
    err = parse_format(&bptr, header, &dim_sec, &var_sec, &gattr_sec);
    if (err < NC_NOERR) goto err_out;

    /* Dimension section */
    if (dim_sec) {
        err = parse_dims(&bptr, &header->dims, &var_sec, &gattr_sec);
        if (err < NC_NOERR) goto err_out;
    }

    /* Variable section */
    if (var_sec) {
        err = parse_vars(&bptr, header, &gattr_sec);
        if (err < NC_NOERR) goto err_out;
    }

    /* Global attribute section */
    if (gattr_sec) {
        err = parse_attr(&bptr, "", &header->attrs);
        if (err < NC_NOERR) goto err_out;
    }

err_out:
    free(fbuf);

    if (err < NC_NOERR) {
        if (header != NULL)
            hdr_free(header);
        return err;
    }

    *hid = cdl_nfiles;
    cdl_filelist[cdl_nfiles++] = header;

    return NC_NOERR;
}

/*----< cdl_hdr_inq_format() >-----------------------------------------------*/
/* Returns the format ID of an opened CDL header file ID.
 * 1: classic CDF-1 format
 * 2: classic CDF-2 format
 * 5: classic CDF-5 format
 */
int cdl_hdr_inq_format(int  hid,
                       int *format)
{
    CDL_header *header;

    /* check if hid is valid */
    if (hid >= NC_MAX_NFILES || cdl_filelist[hid] == NULL)
        return NC_EBADID;

   if (format == NULL)
       return NC_NOERR;

    header = cdl_filelist[hid];

    *format = header->format;
    return NC_NOERR;
}

/*----< cdl_hdr_inq_ndims() >------------------------------------------------*/
/* Returns the number of dimensions defined in the CDL header file. */
int cdl_hdr_inq_ndims(int  hid,
                      int *ndims)
{
    CDL_header *header;

    /* check if hid is valid */
    if (hid >= NC_MAX_NFILES || cdl_filelist[hid] == NULL)
        return NC_EBADID;

    if (ndims == NULL)
        return NC_NOERR;

    header = cdl_filelist[hid];

    *ndims = header->dims.nelems;
    return NC_NOERR;
}

/*----< cdl_hdr_inq_dim() >--------------------------------------------------*/
/* Returns the metadata of a given dimension, specified by its ID, 'dimid',
 * including its name and size. Note users should not alter the contents of
 * character string pointed by 'name'.
 */
int cdl_hdr_inq_dim(int          hid,
                    int          dimid,
                    char       **name,
                    MPI_Offset  *size)
{
    CDL_header *header;

    /* check if hid is valid */
    if (hid >= NC_MAX_NFILES || cdl_filelist[hid] == NULL)
        return NC_EBADID;

    header = cdl_filelist[hid];

    /* check dimid if valid */
    if (dimid < 0 || dimid >= header->dims.nelems)
        return NC_EBADID;

    if (name != NULL) *name = header->dims.value[dimid].name;
    if (size != NULL) *size = header->dims.value[dimid].size;
    return NC_NOERR;
}

/*----< cdl_hdr_inq_nattrs() >-----------------------------------------------*/
/* Returns the number of attributes of a given variable, specified by its ID,
 * 'varid'.
 */
int cdl_hdr_inq_nattrs(int  hid,
                       int  varid,
                       int *nattrs)
{
    CDL_header *header;

    /* check if hid is valid */
    if (hid >= NC_MAX_NFILES || cdl_filelist[hid] == NULL)
        return NC_EBADID;

    header = cdl_filelist[hid];

    /* check varid if valid */
    if (varid == NC_GLOBAL)
        *nattrs = header->attrs.nelems;
    else if (varid < NC_GLOBAL || varid >= header->vars.nelems)
        return NC_EBADID;
    else
        *nattrs = header->vars.value[varid].attrs.nelems;

    return NC_NOERR;
}

/*----< cdl_hdr_inq_attr() >-------------------------------------------------*/
/* Returns the metadata of an attribute of a given attribute ID, specified by
 * it ID 'attrid' of a given variable, specified by its ID, 'varid'. The
 * returned metadata includes the attribute's name, external netCDF data type,
 * number of elements, and a pointer to its contents in the external format,
 * i.e. in Big-Endian. Note users should not alter the contents pointed by
 * arguments 'name' and 'value'.
 */
int cdl_hdr_inq_attr(int          hid,
                     int          varid,
                     int          attrid,
                     char       **name,
                     nc_type     *xtype,
                     MPI_Offset  *nelems,
                     void       **value)
{
    CDL_attr *attrp;
    CDL_attrarray *attr_array;
    CDL_header *header;

    /* check if hid is valid */
    if (hid >= NC_MAX_NFILES || cdl_filelist[hid] == NULL)
        return NC_EBADID;

    header = cdl_filelist[hid];

    /* check varid if valid */
    if (varid == NC_GLOBAL)
        attr_array = &header->attrs;
    else if (varid < NC_GLOBAL || varid >= header->vars.nelems)
        return NC_EBADID;
    else
        attr_array = &header->vars.value[varid].attrs;

    /* check attrid if valid */
    if (attrid < 0 || attrid >= attr_array->nelems)
        return NC_EBADID;

    attrp = &attr_array->value[attrid];

    if (name   != NULL) *name   = attrp->name;
    if (xtype  != NULL) *xtype  = attrp->xtype;
    if (nelems != NULL) *nelems = attrp->nelems;
    if (value  != NULL) *value  = attrp->value;

    return NC_NOERR;
}

/*----< cdl_hdr_inq_nvars() >------------------------------------------------*/
/* Returns the number of variables defined in the CDL file, specified by its
 * header ID, 'hid'.
 */
int cdl_hdr_inq_nvars(int  hid,
                      int *nvars)
{
    CDL_header *header;

    /* check if hid is valid */
    if (hid >= NC_MAX_NFILES || cdl_filelist[hid] == NULL)
        return NC_EBADID;

    header = cdl_filelist[hid];

    /* check header if valid */
    if (header == NULL) return 1;

    if(nvars != NULL) *nvars = header->vars.nelems;
    return NC_NOERR;
}

/*----< cdl_hdr_inq_var() >--------------------------------------------------*/
/* Returns the metadata of a variable, specified by its ID, 'varid'. The
 * returned metadata includes the variable's name, external netCDF data type,
 * number of dimensions, and a pointer to an arrays storing all the dimension
 * IDs.  Note users should not alter the contents pointed by arguments 'name'
 * and 'dimids'.
 */
int cdl_hdr_inq_var(int       hid,
                    int       varid,
                    char    **name,
                    nc_type  *xtype,
                    int      *ndims,
                    int     **dimids)
{
    CDL_header *header;

    /* check if hid is valid */
    if (hid >= NC_MAX_NFILES || cdl_filelist[hid] == NULL)
        return NC_EBADID;

    header = cdl_filelist[hid];

    /* check varid if valid */
    if (varid < 0 || varid >= header->vars.nelems)
        return NC_EBADID;

    if (name   != NULL) *name   = header->vars.value[varid].name;
    if (xtype  != NULL) *xtype  = header->vars.value[varid].xtype;
    if (ndims  != NULL) *ndims  = header->vars.value[varid].ndims;
    if (dimids != NULL) *dimids = header->vars.value[varid].dimids;

    return NC_NOERR;
}


#ifdef TEST_RUN
int main(int argc, char **argv)
{
    char *name;
    int i, j, err=0, hid, format, ndims, *dimids, nvars, nattrs;
    void *value;
    nc_type xtype;
    MPI_Offset size, nelems;

    if (argc != 2) {
        printf("Usage: %s <CDL file>\n",argv[0]);
        exit(1);
    }

    /* open the input file in CDL format */
    err = cdl_hdr_open(argv[1], &hid);
    if (err != NC_NOERR) exit(1);
    printf("Input CDF file : %s\n", argv[1]);

    /* retrieve file format information of the input file */
    err = cdl_hdr_inq_format(hid, &format); ERR
    printf("CDF file format: CDF-%d\n", format);

    /* retrieve the number of dimensions defined in the CDL file */
    err = cdl_hdr_inq_ndims(hid, &ndims); ERR
    printf("Number of dimensions: %d\n", ndims);

    /* retrieve metadata of each dimension defined in the CDL file */
    for (i=0; i<ndims; i++) {
        err = cdl_hdr_inq_dim(hid, i, &name, &size); ERR
        printf("\t name %s size " OFFFMT "\n",name, size);
    }

    /* retrieve number of variables defined in the CDL file */
    err = cdl_hdr_inq_nvars(hid, &nvars); ERR
    printf("Number of variables: %d\n", nvars);

    /* retrieve metadata of each variable defined in the CDL file */
    for (i=0; i<nvars; i++) {
        err = cdl_hdr_inq_var(hid, i, &name, &xtype, &ndims, &dimids); ERR

        /* retrieve the number of attributes associated with variable i */
        err = cdl_hdr_inq_nattrs(hid, i, &nattrs); ERR

        printf("\t name %s type %d ndims %d nattr %d\n",
               name, xtype, ndims, nattrs);
        for (j=0; j<ndims; j++)
             printf("\t\tdimid %d\n",dimids[j]);

        for (j=0; j<nattrs; j++) {
            /* retrieve metadata of attribute j associated with variable j */
            err = cdl_hdr_inq_attr(hid, i, j, &name, &xtype, &nelems, &value);
            ERR
            if (xtype == NC_CHAR)
                printf("\t\tattr %s type %d nelems "OFFFMT" (%s)\n",
                       name, xtype,nelems,(char*)value);
            else
                printf("\t\tattr %s type %d nelems "OFFFMT"\n",
                       name, xtype, nelems);
        }
    }

    /* retrieve the number of global attributes */
    err = cdl_hdr_inq_nattrs(hid, NC_GLOBAL, &nattrs); ERR
#ifdef TEST_RUN
    printf("global attrs: nattrs %d\n", nattrs);
#endif

    /* retrieve metadata of each global attribute */
    for (i=0; i<nattrs; i++) {
        err = cdl_hdr_inq_attr(hid, NC_GLOBAL, i, &name, &xtype, &nelems, &value);
        ERR
        if (xtype == NC_CHAR)
            printf("\t name %s type %d nelems "OFFFMT" (%s)\n",
                    name, xtype, nelems,(char*)value);
        else
            printf("\t name %s type %d nelems "OFFFMT"\n",
                    name, xtype, nelems);
    }

err_out:
    /* close CDL file */
    err = cdl_hdr_close(hid); ERR

    return err;
}
#endif

