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
#include <sys/types.h> /* stat() */
#include <sys/stat.h>  /* stat() */
#include <unistd.h>    /* stat() */
#include <errno.h>     /* errno */
#include <math.h>      /* INFINITY */

#include <mpi.h>
#include <pnetcdf.h>
#include <ncmpidiff_core.h>


#define OOM_ERROR { \
    fprintf(stderr, "Error: calloc() out of memory at line %d\n",__LINE__);  \
    exit(1); \
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
get_var_names(char    *opt_arg,
              int     *nvars,
              char ***names)
{
    char *cp=opt_arg, **cpp;
    *nvars = 1;

    /* compute number of variable names in comma-delimited list */
    while (*cp++)
        if (*cp == ',')
            (*nvars)++;

    *names = (char **) calloc((size_t)*nvars, sizeof(char*));
    if (!*names) OOM_ERROR

    cpp = *names;
    /* copy variable names into list */
    for (cp = strtok(opt_arg, ",");
         cp != NULL;
         cp = strtok((char *) NULL, ",")) {

        *cpp = strdup(cp);
        if (!*cpp) OOM_ERROR
        cpp++;
    }
}

/*----< main() >--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    extern char *optarg;
    extern int optind;
    char cmd_opts[1024], **var_names;
    int i, c, rank, nprocs, verbose, quiet, check_tolerance;
    int first_diff, ncid[2], num_vars;
    int check_header, check_variable_list, check_entire_file;
    long long numDIFF;
    double tolerance_ratio, tolerance_difference;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    if (nprocs == 1)
        strcpy(cmd_opts, "ncmpidiff");
    else
        sprintf(cmd_opts, "Rank %d: ncmpidiff", rank);

    for (i=1; i<argc; i++) {
        strcat(cmd_opts, " ");
        strcat(cmd_opts, argv[i]);
    }

    verbose             = 0;
    quiet               = 0;
    check_header        = 0;
    check_variable_list = 0;
    check_entire_file   = 0;
    var_names           = NULL;
    num_vars            = 0;
    check_tolerance     = 0;
    first_diff          = 1;

    while ((c = getopt(argc, argv, "bhqt:v:")) != -1) {
        char *str, *ptr;
        switch(c) {
            case 'h':               /* compare header only */
                check_header = 1;
                break;
            case 'v':               /* variable names */
                /* make list of names of variables specified */
                get_var_names(optarg, &num_vars, &var_names);
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
            default:
                usage(rank, argv[0]);
                break;
        }
    }

    /* quiet mode overwrites verbose */
    if (quiet) verbose = 0;

    if (argc - optind != 2) usage(rank, argv[0]);

    /* check stat of input files, e.g. exist or not */
    if (rank == 0) {
        struct stat sb;
        ncid[0] = (stat(argv[optind],   &sb) != 0) ? errno : 0;
        ncid[1] = (stat(argv[optind+1], &sb) != 0) ? errno : 0;
    }
    MPI_Bcast(ncid, 2, MPI_INT, 0, comm);
    if (rank == 0) {
        if (ncid[0] != 0)
            fprintf(stderr,"Error: ncmpidiff input file \"%s\" (%s)\n",
                    argv[optind], strerror(ncid[0]));
        if (ncid[1] != 0)
            fprintf(stderr,"Error: ncmpidiff input file \"%s\" (%s)\n",
                    argv[optind+1], strerror(ncid[1]));
    }
    if (ncid[0] != 0 || ncid[1] != 0) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

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
    MPI_Info_create(&info);
    MPI_Info_set(info, "pnetcdf_subfiling", "disable");

    numDIFF = ncmpidiff_core(argv[optind], argv[optind+1],
                             comm, info, verbose, quiet, check_header,
                             check_variable_list, check_entire_file,
                             num_vars, var_names, check_tolerance,
                             first_diff, cmd_opts, tolerance_difference,
                             tolerance_ratio);

    MPI_Info_free(&info);

    if (num_vars > 0) {
        for (i=0; i<num_vars; i++)
            free(var_names[i]);
        free(var_names);
    }

    MPI_Finalize();

    exit((numDIFF == 0) ? EXIT_SUCCESS : EXIT_FAILURE);
}
