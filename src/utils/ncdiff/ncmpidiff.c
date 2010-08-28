/*
 *  Copyright (C) 2010, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* wkliao: This diff utility compares header and variables of two files
 *         regardless the define order of the variables and attributes.
 *
 *         It can also compare a subset of the variables, for example
 *           mpiexec -n 8 ncmpidiff -v var1,var2 file1.nc file2.nc
 *
 *         or compare the header only, for example,
 *           mpiexec -n 8 ncmpidiff -h file1.nc file2.nc
 *
 *         or compare header + a subset of variables, for example,
 *           mpiexec -n 8 ncmpidiff -h -v var1,var2 file1.nc file2.nc
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>
#include <pnetcdf.h>

#define HANDLE_ERROR {                                         \
    if (err != NC_NOERR) {                                     \
        fprintf(stderr, "Error at line %d (%s)\n", __LINE__,   \
               ncmpi_strerror(err));                           \
        MPI_Abort(MPI_COMM_WORLD, -1);                         \
        exit(-1);                                              \
    }                                                          \
}

#define CHECK_GLOBAL_ATT_DIFF(type, func, nctype) {                    \
    int   pos, len = attlen1 * sizeof(type);                           \
    type *b1 = (type *)malloc(len);                                    \
    type *b2 = (type *)malloc(len);                                    \
    err = func(ncid1, NC_GLOBAL, name1, b1);                           \
    HANDLE_ERROR                                                       \
    err = func(ncid2, NC_GLOBAL, name1, b2);                           \
    HANDLE_ERROR                                                       \
    if ((pos = memcmp(b1, b2, len)) != 0) {                            \
        printf("DIFF: global attribute \"%s\" of type \"%s\" (%d)\n",  \
                name1,#nctype,pos);                                    \
        numHeadDIFF++;                                                 \
    }                                                                  \
    else if (verbose)                                                  \
        printf("\tSAME: attribute contents\n");                        \
    free(b1);                                                          \
    free(b2);                                                          \
    break;                                                             \
}

#define CHECK_VAR_ATT_DIFF(type, func, nctype) {                              \
    int   pos, len = attlen1 * sizeof(type);                                  \
    type *b1 = (type *)malloc(len);                                           \
    type *b2 = (type *)malloc(len);                                           \
    err = func(ncid1, i, attrname, b1);                                       \
    HANDLE_ERROR                                                              \
    err = func(ncid2, varid, attrname, b2);                                   \
    HANDLE_ERROR                                                              \
    if ((pos = memcmp(b1, b2, len)) != 0) {                                   \
        printf("DIFF: variable \"%s\" attribute \"%s\" of type \"%s\" (%d)\n",\
               name1,attrname,#nctype,pos);                                   \
        numHeadDIFF++;                                                        \
    }                                                                         \
    else if (verbose)                                                         \
        printf("\t\tSAME: attribute contents\n");                             \
    free(b1);                                                                 \
    free(b2);                                                                 \
    break;                                                                    \
}

#define CHECK_VAR_DIFF(type, func, nctype) {                                 \
    int   pos, isDiff, len = varsize * sizeof(type);                         \
    type *b1 = (type *)malloc(len);                                          \
    type *b2 = (type *)malloc(len);                                          \
    err = func(ncid1, varid1, start, shape, b1);                             \
    HANDLE_ERROR                                                             \
    err = func(ncid2, varid2, start, shape, b2);                             \
    HANDLE_ERROR                                                             \
    if ((pos = memcmp(b1, b2, len)) != 0) {                                  \
        printf("DIFF: variable \"%s\" of type \"%s\" (%d)\n",                \
               name1,#nctype,pos);                                           \
        numVarDIFF++;                                                        \
    }                                                                        \
    MPI_Allreduce(&pos, &isDiff, 1, MPI_INT, MPI_MAX, comm);                 \
    if (isDiff == 0 && !rank && verbose)                                     \
        printf("\tSAME: variable \"%s\" contents\n",name1);                  \
    free(b1);                                                                \
    free(b2);                                                                \
    break;                                                                   \
}

/*----< usage() >-------------------------------------------------------------*/
static void
usage(int rank, char *progname)
{
#define USAGE   "\
  Compare the contents of two netCDF files.\n\n\
  [-b]             Verbose output\n\
  [-h]             Compare header information only, no data\n\
  [-v var1[,...]]  Compare variable(s) <var1>,... only\n\
  file1 file2      File names of two input netCDF files to be compared\n"

    if (rank == 0) {
        printf("  %s [-h] [-v ...] [-b] file1 file2\n%s", progname, USAGE);
        printf("  netcdf library version %s\n", ncmpi_inq_libvers());
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

    vspecp->names = (char **) malloc(nvars * sizeof(char*));
    if (!vspecp->names) error("out of memory");

    cpp = vspecp->names;
    /* copy variable names into list */
    for (cp = strtok(optarg, ",");
         cp != NULL;
         cp = strtok((char *) NULL, ",")) {

        *cpp = (char *) malloc(strlen(cp) + 1);
        if (!*cpp) error("out of memory");
        strcpy(*cpp, cp);
        cpp++;
    }
    vspecp->nvars = nvars;
}

/*----< get_type() >----------------------------------------------------------*/
static char*
get_type(int type)
{
    switch (type) {
        case NC_CHAR:   return "char";
        case NC_SHORT:  return "short";
        case NC_INT:    return "int";
        case NC_FLOAT:  return "float";
        case NC_DOUBLE: return "double";
    }
    return "";
}

/*----< main() >--------------------------------------------------------------*/
int main(int argc, char **argv) {
    int i, j, err, isRecvar, rank, nprocs;
    int ncid1, ndims1, nvars1, natts1, unlimdimid1, *dimids1;
    int ncid2, ndims2, nvars2, natts2, unlimdimid2, *dimids2;
    char c, str[512], *name1, *name2;
    MPI_Offset *shape, varsize, *start;
    MPI_Offset attlen1, dimlen1, attlen2, dimlen2;
    ncmpi_type type1, type2;
    MPI_Comm comm=MPI_COMM_WORLD;
    int nvars, verbose, check_header, check_variable_list, check_entire_file;
    int numVarDIFF=0, numHeadDIFF=0, varDIFF;
    struct vspec var_list;
    extern char *optarg;
    extern int optind;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    verbose             = 0;
    check_header        = 0;
    check_variable_list = 0;
    check_entire_file   = 0;
    var_list.names      = NULL;
    var_list.nvars      = 0;

    while ((c = getopt(argc, argv, "bhv:")) != EOF)
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
            case '?':
                usage(rank, argv[0]);
                break;
        }

    if (argc - optind != 2) usage(rank, argv[0]);

    if (check_header == 0 && check_variable_list == 0) {
        check_entire_file = 1;
        check_header      = 1;
    }

    str[0] = '\0';
    dimids1 = (int*) malloc(NC_MAX_DIMS * sizeof(int));
    dimids2 = (int*) malloc(NC_MAX_DIMS * sizeof(int));
    name1   = (char*) malloc(NC_MAX_NAME);
    name2   = (char*) malloc(NC_MAX_NAME);
    shape   = (MPI_Offset*) malloc(NC_MAX_VAR_DIMS * sizeof(MPI_Offset));
    start   = (MPI_Offset*) malloc(NC_MAX_VAR_DIMS * sizeof(MPI_Offset));

    /* open files */
    err = ncmpi_open(comm, argv[optind], NC_NOWRITE, MPI_INFO_NULL, &ncid1);
    HANDLE_ERROR
    err = ncmpi_open(comm, argv[optind+1], NC_NOWRITE, MPI_INFO_NULL, &ncid2);
    HANDLE_ERROR

    /* check header */
    if (check_header && rank == 0) { /* only root checks header */
        int attnump;

        err = ncmpi_inq(ncid1, &ndims1, &nvars1, &natts1, &unlimdimid1);
        HANDLE_ERROR
        err = ncmpi_inq(ncid2, &ndims2, &nvars2, &natts2, &unlimdimid2);
        HANDLE_ERROR
        if (ndims1 != ndims2) { /* check number of dimensions if equal */
            printf("DIFF: number of dimensions (%d) != (%d)\n",ndims1, ndims2);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of dimensions (%d)\n",ndims1);
        if (nvars1 != nvars2) { /* check number of variables if equal */
            printf("DIFF: number of variables (%d) != (%d)\n",nvars1, nvars2);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of variables (%d)\n",nvars1);
        if (natts1 != natts2) { /* check number of global attributes if equal */
            printf("DIFF: number of global attributes (%d) != (%d)\n",natts1, natts2);
            numHeadDIFF++;
        }
        else if (verbose)
            printf("SAME: number of global attributes (%d)\n",natts1);

        /* Compare global attributes, assume CHAR attributes. */
        for (i=0; i<natts1; i++) { /* check what's in file1 also in file2 */
            err = ncmpi_inq_attname(ncid1, NC_GLOBAL, i, name1);
            HANDLE_ERROR
            /* find the attr with the same name from ncid2 */
            err = ncmpi_inq_attid(ncid2, NC_GLOBAL, name1, &attnump);
            if (err == NC_ENOTATT) {
                printf("DIFF: global attribute \"%s\" not found in file %s\n",
                       name1,argv[optind+1]);
                numHeadDIFF++;
                continue;
            }

            err = ncmpi_inq_att(ncid1, NC_GLOBAL, name1, &type1, &attlen1);
            HANDLE_ERROR
            err = ncmpi_inq_att(ncid2, NC_GLOBAL, name1, &type2, &attlen2);
            HANDLE_ERROR
            if (type1 != type2) {
                printf("DIFF: global attribute \"%s\" data type (%s) != (%s)\n",
                       name1,get_type(type1),get_type(type2));
                numHeadDIFF++;
                continue;
            }
            else if (verbose) {
                printf("Global attribute \"%s\":\n",name1);
                printf("\tSAME: data type (%s)\n",get_type(type1));
            }

            if (attlen1 != attlen2) {
                printf("DIFF: global attribute \"%s\" length (%lld) != (%lld)\n",
                       name1,(long long int)attlen1,(long long int)attlen2);
                numHeadDIFF++;
                continue;
            }
            else if (verbose)
                printf("\tSAME: length (%lld)\n",(long long int)attlen1);

            switch (type1) {
                case NC_CHAR:   CHECK_GLOBAL_ATT_DIFF(char,   ncmpi_get_att_text,   NC_CHAR)
                case NC_SHORT:  CHECK_GLOBAL_ATT_DIFF(short,  ncmpi_get_att_short,  NC_SHORT)
                case NC_INT:    CHECK_GLOBAL_ATT_DIFF(int,    ncmpi_get_att_int,    NC_INT)
                case NC_FLOAT:  CHECK_GLOBAL_ATT_DIFF(float,  ncmpi_get_att_float,  NC_FLOAT)
                case NC_DOUBLE: CHECK_GLOBAL_ATT_DIFF(double, ncmpi_get_att_double, NC_DOUBLE)
                default: ; /* TODO: handle unexpected types */
            }
        }
        for (i=0; i<natts2; i++) { /* check attributes in file2 but not in file1 */
            err = ncmpi_inq_attname(ncid2, NC_GLOBAL, i, name2);
            HANDLE_ERROR
            /* find the attr with the same name from ncid1 */
            if (ncmpi_inq_attid(ncid1, NC_GLOBAL, name2, &attnump) == NC_ENOTATT) {
                numHeadDIFF++;
                printf("DIFF: global attribute \"%s\" not found in file %s\n",
                       name1,argv[optind]);
            }
        }

        /* Compare dimension */
        if (ndims1 && verbose)
            printf("Dimension:\n");

        for (i=0; i<ndims1; i++) { /* check dimensions in file1 also in file2 */
            int dimid;
            err = ncmpi_inq_dim(ncid1, i, name1, &dimlen1);
            HANDLE_ERROR
            /* find the dim with the same name from ncid2 */
            err = ncmpi_inq_dimid(ncid2, name1, &dimid);
            if (err == NC_EBADDIM) {
                printf("DIFF: dimension \"%s\" not found in file %s\n",
                       name1,argv[optind+1]);
                numHeadDIFF++;
                continue;
            }

            err = ncmpi_inq_dimlen(ncid2, dimid, &dimlen2);
            HANDLE_ERROR
            if (dimlen1 != dimlen2) {
                /* cast to quiet warning on 32 bit platforms */
                printf("DIFF: dimension \"%s\" length (%lld) != (%lld)\n",
                       name1,(long long int)dimlen1,(long long int)dimlen2);
                numHeadDIFF++;
            }
            else if (verbose)
                printf("\tSAME: dimension \"%s\" length (%lld)\n",
                       name1,(long long int)dimlen1);
        }
        for (i=0; i<ndims2; i++) { /* check dimensions in file2 but not in file1 */
            int dimid;
            err = ncmpi_inq_dim(ncid2, i, name2, &dimlen2);
            HANDLE_ERROR
            /* find the dim with the same name from ncid1 */
            if (ncmpi_inq_dimid(ncid2, name1, &dimid) == NC_EBADDIM) {
                printf("DIFF: dimension \"%s\" not found in file %s\n",
                       name1,argv[optind]);
                numHeadDIFF++;
            }
        }

        /* Compare variables' metadata */
        for (i=0; i<nvars1; i++) {
            int varid;
            err = ncmpi_inq_var(ncid1, i, name1, &type1, &ndims1, dimids1, &natts1);
            HANDLE_ERROR
            /* find the variable with the same name from ncid2 */
            err = ncmpi_inq_varid(ncid2, name1, &varid);
            if (err == NC_ENOTVAR) {
                printf("DIFF: variable \"%s\" not found in file %s\n",
                       name1,argv[optind+1]);
                numHeadDIFF++;
                numVarDIFF++;
                continue;
            }
            err = ncmpi_inq_var(ncid2, varid, name2, &type2, &ndims2, dimids2, &natts2);
            HANDLE_ERROR

            if (type1 != type2) {
                printf("DIFF: variable \"%s\" data type (%s) != (%s)\n",
                       name1,get_type(type1),get_type(type2));
                numHeadDIFF++;
            }
            else if (verbose) {
                printf("Variable \"%s\":\n",name1);
                printf("\tSAME: data type (%s)\n",get_type(type1));
            }

            if (ndims1 != ndims2) {
                printf("DIFF: variable \"%s\" number of dimensions (%d) != (%d)\n",
                       name1,ndims1,ndims2);
                numHeadDIFF++;
            }
            else {
                if (verbose)
                    printf("\tSAME: number of dimensions (%d)\n",ndims1);

                for (j=0; j<ndims1; j++) { /* check variable's dimensionality */
                    char dimname1[NC_MAX_NAME], dimname2[NC_MAX_NAME];
                    /* get dim name for each dim ID */
                    err = ncmpi_inq_dim(ncid1, dimids1[j], dimname1, &dimlen1);
                    HANDLE_ERROR
                    err = ncmpi_inq_dim(ncid1, dimids2[j], dimname2, &dimlen2);
                    HANDLE_ERROR
                    if (verbose)
                        printf("\tdimension %d:\n",j);
                    if (strcmp(dimname1, dimname2) != 0) {
                        printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's name (%s) != (%s)\n",
                               name1,get_type(type1),j,dimname1,dimname2);
                        numHeadDIFF++;
                    }
                    else if (verbose)
                        printf("\t\tSAME: name (%s)\n",dimname1);
                    if (dimlen1 != dimlen2) {
                        printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's length (%lld) != (%lld)\n",
                               name1,get_type(type1),j,(long long int)dimlen1,(long long int)dimlen2);
                        numHeadDIFF++;
                    }
                    else if (verbose)
                        printf("\t\tSAME: length (%lld)\n",(long long int)dimlen1);
                }
            }

            if (natts1 != natts2) {
                printf("DIFF: variable \"%s\" number of attributes (%d) != (%d)\n",
                       name1,natts1,natts2);
                numHeadDIFF++;
            }
            else if (verbose)
                printf("\tSAME: number of attributes (%d)\n",natts1);

            /* var attributes, assume CHAR attributes */
            for (j=0; j<natts1; j++) {
                char attrname[NC_MAX_NAME];
                err = ncmpi_inq_attname(ncid1, i, j, attrname);
                HANDLE_ERROR
                err = ncmpi_inq_att(ncid1, i, attrname, &type1, &attlen1);
                HANDLE_ERROR
                /* find the variable attr with the same name from ncid2 */
                err = ncmpi_inq_att(ncid2, varid, attrname, &type2, &attlen2);
                if (err == NC_ENOTATT) {
                    printf("DIFF: variable \"%s\" attribute \"%s\" not found in file %s\n",
                           name1,attrname,argv[optind+1]);
                    numHeadDIFF++;
                    continue;
                }
                if (verbose)
                    printf("\tattribute \"%s\":\n",attrname);

                if (type1 != type2) {
                    printf("DIFF: variable \"%s\" attribute \"%s\" data type (%s) != (%s)\n",
                           name1,attrname,get_type(type1),get_type(type2));
                    numHeadDIFF++;
                    continue; /* skip this attribute */
                }
                else if (verbose)
                    printf("\t\tSAME: data type (%s)\n",get_type(type1));
                if (attlen1 != attlen2) {
                    printf("DIFF: variable \"%s\" attribute \"%s\" length (%lld) != (%lld)\n",
                           name1,attrname,(long long int)attlen1,(long long int)attlen2);
                    numHeadDIFF++;
                    continue; /* skip this attribute */
                }
                else if (verbose)
                    printf("\t\tSAME: length (%lld)\n",(long long int)attlen1);

                switch (type1) {
                    case NC_CHAR:   CHECK_VAR_ATT_DIFF(char,   ncmpi_get_att_text,   NC_CHAR)
                    case NC_SHORT:  CHECK_VAR_ATT_DIFF(short,  ncmpi_get_att_short,  NC_SHORT)
                    case NC_INT:    CHECK_VAR_ATT_DIFF(int,    ncmpi_get_att_int,    NC_INT)
                    case NC_FLOAT:  CHECK_VAR_ATT_DIFF(float,  ncmpi_get_att_float,  NC_FLOAT)
                    case NC_DOUBLE: CHECK_VAR_ATT_DIFF(double, ncmpi_get_att_double, NC_DOUBLE)
                    default: ; /* TODO: handle unexpected types */
                }
            }
            for (j=0; j<natts2; j++) {
                char attrname[NC_MAX_NAME];
                err = ncmpi_inq_attname(ncid2, varid, j, attrname);
                HANDLE_ERROR
                /* find the variable attr with the same name from ncid1 */
                err = ncmpi_inq_att(ncid1, i, attrname, &type1, &attlen1);
                if (err == NC_ENOTATT) {
                    printf("DIFF: variable \"%s\" attribute \"%s\" not found in file %s\n",
                           name1,attrname,argv[optind]);
                    numHeadDIFF++;
                }
            }
        }
        for (i=0; i<nvars2; i++) { /* check variables in file2 but not in file1 */
            int varid;
            err = ncmpi_inq_varname(ncid2, i, name2);
            HANDLE_ERROR
            /* find the variable with the same name from ncid1 */
            err = ncmpi_inq_varid(ncid1, name2, &varid);
            if (err == NC_ENOTVAR) {
                printf("DIFF: variable \"%s\" not found in file %s\n",
                       name2,argv[optind]);
                numHeadDIFF++;
                numVarDIFF++;
            }
        }
    }

    /* compare variable contents */
    for (i=0; i<NC_MAX_VAR_DIMS; i++) start[i] = 0;

    nvars = 0;
    if (check_variable_list) nvars = var_list.nvars;

    if (check_entire_file) { /* header has been checked */
        ncmpi_inq_nvars(ncid1, &nvars);
        var_list.nvars = nvars;
        var_list.names = (char**) malloc(nvars * sizeof(char*));
        /* get all the variable names from file1 */
        for (i=0; i<nvars; i++) {
            ncmpi_inq_varname(ncid1, i, name1);
            var_list.names[i] = (char *) malloc(strlen(name1) + 1);
            strcpy(var_list.names[i], name1);
        }
    }
    if (!rank && verbose) printf("number of variables to be compared = %d\n",nvars);

    for (i=0; i<nvars; i++) { /* compare one variable at a time */
        int varid1, varid2;

        err = ncmpi_inq_varid(ncid1, var_list.names[i], &varid1);
        if (err == NC_ENOTVAR) {
            if (!check_header && !rank) {
                printf("WARN: variable \"%s\" not found in file %s\n",
                       var_list.names[i],argv[optind]);
                numVarDIFF++;
            }
            continue;
        }
        err = ncmpi_inq_varid(ncid2, var_list.names[i], &varid2);
        if (err == NC_ENOTVAR) {
            if (!check_header && !rank) {
                printf("WARN: variable \"%s\" not found in file %s\n",
                       var_list.names[i],argv[optind+1]);
                numVarDIFF++;
            }
            continue;
        }
        err = ncmpi_inq_var(ncid1, varid1, name1, &type1, &ndims1, dimids1, &natts1);
        HANDLE_ERROR
        err = ncmpi_inq_var(ncid2, varid2, name2, &type2, &ndims2, dimids2, &natts2);
        HANDLE_ERROR

        /* check data type */
        if (type1 != type2) {
            if (!rank) {
                if (!check_header) /* if header has not been checked */
                    printf("DIFF: variable \"%s\" data type (%s) != (%s)\n",
                       name1,get_type(type1),get_type(type2));
                numHeadDIFF++;
                numVarDIFF++;
            }
            continue; /* skip this variable */
        }
        else if (!check_header && !rank && verbose) {
            printf("Variable \"%s\":\n",name1);
            printf("\tSAME: data type (%s)\n",get_type(type1));
        }

        /* check number of dimensions */
        if (ndims1 != ndims2) {
            if (!rank) {
                if (!check_header) /* if header has not been checked */
                    printf("DIFF: variable \"%s\" number of dimensions (%d) != (%d)\n",
                       name1,ndims1,ndims2);
                numHeadDIFF++;
                numVarDIFF++;
            }
            continue; /* skip this variable */
        }
        else if (!check_header && !rank && verbose)
            printf("\tSAME: number of dimensions (%d)\n",ndims1);

        /* check dimension length only */
        for (j=0; j<ndims1; j++) { /* check variable's dimensionality */
            err = ncmpi_inq_dimlen(ncid1, dimids1[j], &dimlen1);
            HANDLE_ERROR
            err = ncmpi_inq_dimlen(ncid2, dimids2[j], &dimlen2);
            HANDLE_ERROR
            if (!check_header && !rank && verbose)
                printf("\tDimension %d:\n",j);
            if (dimlen1 != dimlen2) {
                if (!rank) {
                    if (!check_header) /* if header has not been checked */
                        printf("DIFF: variable \"%s\" of type \"%s\" dimension %d's length (%lld) != (%lld)\n",
                               name1,get_type(type1),j,(long long int)dimlen1,(long long int)dimlen2);
                    numHeadDIFF++;
                    numVarDIFF++;
                }
                break; /* skip this variable */
            }
            else if (!check_header && !rank && verbose)
                printf("\t\tSAME: length (%lld)\n",(long long int)dimlen1);
            shape[j] = dimlen1;
        }
        if (j != ndims1)
            continue; /* skip this variable */

        isRecvar = 0;
        varsize  = 1;
        /* block partition the variable along the 1st dimension */
        for (j=0; j<ndims1; j++) {
            if (j == 0) {
                shape[j] /= nprocs;
                start[j] = shape[j] * rank;
            }
            varsize *= shape[j];
            if (dimids1[j] == unlimdimid1)
                isRecvar = 1;
        }
 
        /* compare the variable contents */
        switch (type1) {
            case NC_CHAR:   CHECK_VAR_DIFF(char,   ncmpi_get_vara_text_all,   NC_CHAR)
            case NC_SHORT:  CHECK_VAR_DIFF(short,  ncmpi_get_vara_short_all,  NC_SHORT)
            case NC_INT:    CHECK_VAR_DIFF(int,    ncmpi_get_vara_int_all,    NC_INT)
            case NC_FLOAT:  CHECK_VAR_DIFF(float,  ncmpi_get_vara_float_all,  NC_FLOAT)
            case NC_DOUBLE: CHECK_VAR_DIFF(double, ncmpi_get_vara_double_all, NC_DOUBLE)
            default: ; /* TODO: handle unexpected types */
        }
    }

    /* close files */
    err = ncmpi_close(ncid1);
    HANDLE_ERROR
    err = ncmpi_close(ncid2);
    HANDLE_ERROR

    /* summary of the difference */
    MPI_Reduce(&numVarDIFF, &varDIFF, 1, MPI_INT, MPI_SUM, 0, comm);
    if (rank == 0) {
        if (check_header) {
            if (numHeadDIFF == 0)
                printf("Headers of two files are the same\n");
            else
                printf("Number of differences in header %d\n",numHeadDIFF);
        }
        if (check_variable_list) {
            if (varDIFF == 0)
                printf("Compared variable(s) are the same\n");
            else
                printf("Compared variables(s) has %d differences\n",varDIFF);
        }
        if (check_entire_file) {
            if (varDIFF == 0)
                printf("All variables of two files are the same\n");
            else
                printf("Number of differences in variables %d\n",varDIFF);
        }
    }

    /* free up the memory previously allocated */
    if (var_list.nvars) {
        for (i=0; i<var_list.nvars; i++)
            free(var_list.names[i]);
        free(var_list.names);
    }
    free(dimids1);
    free(dimids2);
    free(name1);
    free(name2);
    free(shape);
    free(start);

    MPI_Finalize();
    return 0;
}
