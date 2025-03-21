#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include <mpi.h>
#include <pnetcdf.h>

static int debug;

#define LINE_SIZE 1024

#define ERR_FORMAT(msg) { \
    printf("Error line %d: input file format at %s\n", __LINE__,msg); \
    err = -1; \
    goto err_out; \
}

#define DIM_ARRAY_GROWBY 2
#define VAR_ARRAY_GROWBY 2
#define ATTR_ARRAY_GROWBY 2

typedef struct {
    char       *name;
    MPI_Offset  size;
} HR_dim;

typedef struct {
    int     nelems;      /* number of defined dimensions */
    HR_dim *value;
} HR_dimarray;

typedef struct {
    char       *name;     /* name of the attributes */
    nc_type     xtype;    /* external NC data type of the attribute */
    MPI_Offset  nelems;   /* number of attribute elements */
    void       *xvalue;   /* the actual data, in external representation */
} HR_attr;

typedef struct {
    int      nelems;  /* number of defined attributes */
    HR_attr *value;
} HR_attrarray;

typedef struct {
    char         *name;    /* name of the variable */
    nc_type       xtype;   /* variable's external NC data type */
    int           ndims;   /* number of dimensions */
    int          *dimids;  /* [ndims] array of dimension IDs */
    HR_attrarray  attrs;   /* attribute array */
} HR_var;

typedef struct {
    int     nelems;    /* number of defined variables */
    HR_var *value;
} HR_vararray;

typedef struct {
    int           format;   /* 1, 2, or 5 corresponding to CDF-1, 2, or 5 */
    HR_dimarray   dims;     /* dimensions defined */
    HR_attrarray  attrs;    /* global attributes defined */
    HR_vararray   vars;     /* variables defined */
} NC_header;

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
int get_dimid(NC_header   *header,
              char *name)
{
    int i;
    HR_dim *dimp = header->dims.value;

    for (i=0; i<header->dims.nelems; i++)
        if (!strcmp(dimp[i].name, name))
            return i;

    printf("Error: failed to find dim ID for %s\n", name);

    return -1;
}

static
int parse_signature(char **bptr)
{
    char *line;
    int err=0;

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
        if (debug) printf("LINE %d signature: %s\n",__LINE__, signature);
        char *filename = strtok(NULL, " ");
        if (filename == NULL)
            ERR_FORMAT("input file name")
        if (debug) printf("LINE %d input file name: %s\n",__LINE__, filename);
    }
err_out:
    return err;
}

static
int parse_format(char      **bptr,
                 NC_header  *header,
                 int        *dim_sec,
                 int        *var_sec,
                 int        *gattr_sec)
{
    char *line, *key, *val;
    int err=0;

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
                if (debug) printf("LINE %d CDF file format: CDF-%d\n", __LINE__,header->format);
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
int parse_dims(char        **bptr,
               HR_dimarray  *dims,
               int          *var_sec,
               int          *gattr_sec)
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
            dims->value = (HR_dim*) realloc(dims->value,
                                            sizeof(HR_dim) * len);
        }
        HR_dim *dimp = &dims->value[dims->nelems];
        dimp->name = strdup(key);
        if (!strcmp(val, "UNLIMITED"))
            dimp->size = NC_UNLIMITED;
        else
            dimp->size = atoi(val);
        dims->nelems++;

        if (debug) printf("LINE %d DIM name %s val %s\n",__LINE__, key, val);
        if (key[0] == '}') return 1;
        if (!strncmp(key, "data:", 5)) return 1; /* data section */
    }

    return err;
}

#define PARSE_ATTR(itype, tail, conv) { \
    itype *buf = (itype*) attrp->xvalue; \
    attrp->nelems=0; \
    key = strtok(str, ", ;"); \
    while (key != NULL) { \
        if (tail) key[strlen(key) - tail] = '\0'; \
        buf[attrp->nelems++] = conv(key); \
        key = strtok(NULL, ", ;"); \
    } \
}

static
int parse_attr_value(HR_attr *attrp,
                     char    *val)
{
    int err=0;

    if (val[1] == '\"') {
        attrp->xtype = NC_CHAR;
        char *tail = strrchr(val, '\"');
        *tail = '\0';
        attrp->xvalue = strdup(val+2);
        if (debug) printf("LINE %d ATTR name %s type NC_CHAR %s\n",__LINE__, attrp->name, (char*)attrp->xvalue);
        attrp->nelems = strlen(attrp->xvalue);
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
        if (debug) printf("LINE %d ATTR type %d name %s nelems=%d\n",__LINE__, xtype, attrp->name, nelems);

        strcpy(str, val);

        attrp->xvalue = (void*) malloc(sizeof(double) * nelems);

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
int parse_attr(char         **bptr,
               char          *var_name,
               HR_attrarray  *attrs)
{
    char *line, line_cpy[LINE_SIZE], *key, *val;
    int err=0;
    size_t prefix_len;

    /* parse attributes */
    while ((line = strtok(*bptr, "\n")) != NULL) {
        if (line[0] == '}') return 1; /* end of file */

        *bptr += strlen(line) + 1;
        while (**bptr == '\n') (*bptr)++;

        strcpy(line_cpy, line);

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

        /* check if an attribute of this variable */
        if (strcmp(name, var_name)) {
            *bptr = line;
            break;
        }
        if (attrs->nelems % ATTR_ARRAY_GROWBY == 0) {
            size_t len = attrs->nelems + ATTR_ARRAY_GROWBY;
            attrs->value = (HR_attr*) realloc(attrs->value,
                                              sizeof(HR_attr) * len);
        }
        HR_attr *attrp = &attrs->value[attrs->nelems];

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
int parse_vars(char    **bptr,
               NC_header *header,
               int *gattr_sec)
{
    char *type, *var_name;
    char *line, orig_line[LINE_SIZE], *key, *val;
    int err=0;
    size_t prefix_len;
    HR_vararray *vars = &header->vars;

    while ((line = strtok(*bptr, "\n")) != NULL) {
        if (line[0] == '}') return 1; /* end of file */

        *bptr += strlen(line) + 1;
        while (**bptr == '\n') (*bptr)++;

        strncpy(orig_line, line, LINE_SIZE);

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
            vars->value = (HR_var*) realloc(vars->value,
                                            sizeof(HR_var) * len);
        }
        HR_var *varp = &vars->value[vars->nelems];
        var_name = strtok(NULL, "( ");
        varp->name = strdup(var_name);
        varp->xtype = c_to_xtype(type);
        varp->ndims = 0;
        varp->attrs.nelems = 0;
        varp->attrs.value = NULL;

        val = strtok(NULL, ");");
        if (val == NULL) {
            varp->dimids = NULL;
            if (debug) printf("LINE %d VAR type %s name %s IS SCALAR\n",__LINE__, type, varp->name);
        }
        else {
            char *str=strdup(val);
            key = strtok(str, ", ");
            while (key != NULL) {
                varp->ndims++;
                key = strtok(NULL, ", ");
            }
            if (debug) printf("LINE %d VAR type %s name %s val %s ndims=%d\n",__LINE__, type, varp->name, val, varp->ndims);

            varp->dimids = (int*) malloc(sizeof(int) * varp->ndims);
            strcpy(str, val);
            key = strtok(str, ", ");
            int i = 0;
            while (key != NULL) {
                int dimid = get_dimid(header, key);
                if (dimid < 0) {
                    free(str);
                    free(varp->name);
                    varp->name = NULL;
                    free(varp->dimids);
                    varp->dimids = NULL;
                    return -1;
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

int nc_header_free(NC_header *header)
{
    int i, j, err=0;

    /* free dimensions */
    for (i=0; i<header->dims.nelems; i++)
        free(header->dims.value[i].name);

    if (header->dims.value != NULL)
        free(header->dims.value);

    /* free global attributes */
    for (i=0; i<header->attrs.nelems; i++) {
        free(header->attrs.value[i].name);
        free(header->attrs.value[i].xvalue);
    }

    if (header->attrs.value != NULL)
        free(header->attrs.value);

    /* free variables */
    for (i=0; i<header->vars.nelems; i++) {
        HR_var *varp = &header->vars.value[i];

        if (varp->name != NULL) free(varp->name);
        if (varp->dimids != NULL) free(varp->dimids);

        /* free attributes */
        for (j=0; j<varp->attrs.nelems; j++) {
            free(varp->attrs.value[j].name);
            free(varp->attrs.value[j].xvalue);
        }
        if (varp->attrs.value != NULL)
            free(varp->attrs.value);
    }

    if (header->vars.value != NULL)
        free(header->vars.value);

    free(header);

    return err;
}

NC_header *nc_header_parse(char *filename)
{
    FILE *fptr;
    char *fbuf, *bptr;
    int err=0, dim_sec, var_sec, gattr_sec;
    NC_header *header=NULL;

    fptr = fopen(filename, "r");
    if (fptr == NULL) {
        printf("Error: fail to open file %s (%s)\n",
               filename,strerror(errno));
        return NULL;
    }
    err = fseek(fptr, 0, SEEK_END);
    if (err < 0) {
        printf("Error: fail to fseek file %s (%s)\n",
               filename,strerror(errno));
        return NULL;
    }

    /* read the entire file into a buffer */
    long file_size = ftell(fptr);
    fbuf = (char *) malloc(file_size);
    err = fseek(fptr, 0, SEEK_SET);
    fread(fbuf, 1, file_size, fptr);
    bptr = fbuf;
    fclose(fptr);

    /* check netCDF file signature */
    err = parse_signature(&bptr);
    if (err) goto err_out;

    header = (NC_header*) calloc(1, sizeof(NC_header));

    dim_sec = var_sec = gattr_sec = 0;

    /* check netCDF file format */
    err = parse_format(&bptr, header, &dim_sec, &var_sec, &gattr_sec);
    if (err) goto err_out;

    /* Dimension section */
    if (dim_sec) {
        err = parse_dims(&bptr, &header->dims, &var_sec, &gattr_sec);
        if (err) goto err_out;
    }

    /* Variable section */
    if (var_sec) {
        err = parse_vars(&bptr, header, &gattr_sec);
        if (err) goto err_out;
    }

    /* Global attribute section */
    if (gattr_sec) {
        err = parse_attr(&bptr, "", &header->attrs);
        if (err) goto err_out;
    }

err_out:
    free(fbuf);

    if (err == -1) {
        nc_header_free(header);
        return NULL;
    }

    return header;
}

#define TEST_RUN

#ifdef TEST_RUN
#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error at %s:%d : %s\n", __FILE__,__LINE__, \
               ncmpi_strerrno(err)); \
        err = 1; \
        goto err_out; \
    } \
}


int main(int argc, char **argv)
{
    int i, j, err=0;
    NC_header *header;

    MPI_Init(&argc, &argv);

    if (argc != 2) {
        printf("Usage: %s header_file\n",argv[0]);
        exit(1);
    }

    debug=0;

    header = nc_header_parse(argv[1]);
    if (header == NULL) exit(1);

    if (debug) {
        printf("==================================================\n");
        printf("CDF file format: CDF-%d\n", header->format);
        printf("dim: nelems %d\n", header->dims.nelems);
        for (i=0; i<header->dims.nelems; i++) {
            HR_dim *dimp = header->dims.value+i;
            printf("\t name %s size %lld\n",dimp->name, dimp->size);
        }

        printf("var: nelems %d\n", header->vars.nelems);
        for (i=0; i<header->vars.nelems; i++) {
            HR_var *varp = header->vars.value+i;
            printf("\t name %s type %d ndims %d nattr %d\n",
                    varp->name, varp->xtype, varp->ndims,varp->attrs.nelems);
            for (j=0; j<varp->ndims; j++)
                printf("\t\tdimid %d\n",varp->dimids[j]);
            for (j=0; j<varp->attrs.nelems; j++) {
                HR_attr *attrp = varp->attrs.value + j;
                if (attrp->xtype == NC_CHAR)
                    printf("\t\tattr %s type %d nelems %lld (%s)\n",
                        attrp->name, attrp->xtype,attrp->nelems,(char*)attrp->xvalue);
                else
                    printf("\t\tattr %s type %d nelems %lld\n",
                        attrp->name, attrp->xtype,attrp->nelems);
            }
        }
        printf("global attrs: nelems %d\n", header->attrs.nelems);
        for (i=0; i<header->attrs.nelems; i++) {
            HR_attr *attrp = header->attrs.value+i;
            if (attrp->xtype == NC_CHAR)
                printf("\t name %s type %d nelems %lld (%s)\n",
                    attrp->name, attrp->xtype, attrp->nelems,(char*)attrp->xvalue);
            else
                printf("\t name %s type %d nelems %lld\n",
                    attrp->name, attrp->xtype, attrp->nelems);
        }
    }

    /* create a new netcdf file */
    int ncid, cmode;

    cmode = NC_CLOBBER;
    if (header->format == 2) cmode |= NC_64BIT_OFFSET;
    else if (header->format == 5) cmode |= NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, "testfile.nc", cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* define dimensions */
    for (i=0; i<header->dims.nelems; i++) {
        int dimid;
        HR_dim *dimp = header->dims.value+i;
        err = ncmpi_def_dim(ncid, dimp->name, dimp->size, &dimid); ERR
    }

    /* define variables */
    for (i=0; i<header->vars.nelems; i++) {
        int varid;
        HR_var *varp = header->vars.value+i;
        err = ncmpi_def_var(ncid, varp->name, varp->xtype, varp->ndims,
                            varp->dimids, &varid); ERR

        /* define local attributes */
        for (j=0; j<varp->attrs.nelems; j++) {
            HR_attr *attrp = varp->attrs.value+j;
            err = ncmpi_put_att(ncid, varid, attrp->name, attrp->xtype,
                                attrp->nelems, attrp->xvalue); ERR
        }
    }

    /* define global attributes */
    for (i=0; i<header->attrs.nelems; i++) {
        HR_attr *attrp = header->attrs.value+i;
        err = ncmpi_put_att(ncid, NC_GLOBAL, attrp->name, attrp->xtype,
                                    attrp->nelems, attrp->xvalue); ERR
    }
    err = ncmpi_close(ncid); ERR

err_out:
    err = nc_header_free(header);

    MPI_Finalize();
    return err;
}
#endif

