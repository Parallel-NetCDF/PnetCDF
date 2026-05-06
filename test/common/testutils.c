/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */


#include <stdio.h>
#include <libgen.h> /* basename() */
#include <limits.h>
#include <string.h> /* strchr(), strerror(), strdup(), strcpy(), strlen() */
#include <unistd.h> /* getopt(), stat() */
#include <sys/types.h> /* stat() */
#include <sys/stat.h> /* stat() */

#include <mpi.h>

#include <pnetcdf.h>
#include "testutils.h"
#include <ncmpidiff_core.h>

#if defined(PNETCDF_DRIVER_NETCDF4) && PNETCDF_DRIVER_NETCDF4 == 1
int nc_formats[5] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_NETCDF4,
                     NC_FORMAT_NETCDF4_CLASSIC, NC_FORMAT_64BIT_DATA};
#else
int nc_formats[3] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
#endif

char* nc_err_code_name(int err)
{
    static char unknown_str[32];

    if (err > 0) { /* system error */
        const char *cp = (const char *) strerror(err);
        if (cp == NULL)
            sprintf(unknown_str,"Unknown error code %d",err);
        else
            sprintf(unknown_str,"Error code %d (%s)",err,cp);
        return unknown_str;
    }

#define ERR_CODE_STR(err) case (err): return #err;

    switch (err) {
        ERR_CODE_STR(NC_NOERR)
        ERR_CODE_STR(NC_EBADID)
        ERR_CODE_STR(NC_ENFILE)
        ERR_CODE_STR(NC_EEXIST)
        ERR_CODE_STR(NC_EINVAL)
        ERR_CODE_STR(NC_EPERM)
        ERR_CODE_STR(NC_ENOTINDEFINE)
        ERR_CODE_STR(NC_EINDEFINE)
        ERR_CODE_STR(NC_EINVALCOORDS)
        ERR_CODE_STR(NC_EMAXDIMS)
        ERR_CODE_STR(NC_ENAMEINUSE)
        ERR_CODE_STR(NC_ENOTATT)
        ERR_CODE_STR(NC_EMAXATTS)
        ERR_CODE_STR(NC_EBADTYPE)
        ERR_CODE_STR(NC_EBADDIM)
        ERR_CODE_STR(NC_EUNLIMPOS)
        ERR_CODE_STR(NC_EMAXVARS)
        ERR_CODE_STR(NC_ENOTVAR)
        ERR_CODE_STR(NC_EGLOBAL)
        ERR_CODE_STR(NC_ENOTNC)
        ERR_CODE_STR(NC_ESTS)
        ERR_CODE_STR(NC_EMAXNAME)
        ERR_CODE_STR(NC_EUNLIMIT)
        ERR_CODE_STR(NC_ENORECVARS)
        ERR_CODE_STR(NC_ECHAR)
        ERR_CODE_STR(NC_EEDGE)
        ERR_CODE_STR(NC_ESTRIDE)
        ERR_CODE_STR(NC_EBADNAME)
        ERR_CODE_STR(NC_ERANGE)
        ERR_CODE_STR(NC_ENOMEM)
        ERR_CODE_STR(NC_EVARSIZE)
        ERR_CODE_STR(NC_EDIMSIZE)
        ERR_CODE_STR(NC_ETRUNC)
        ERR_CODE_STR(NC_EAXISTYPE)
        ERR_CODE_STR(NC_EDAP)
        ERR_CODE_STR(NC_ECURL)
        ERR_CODE_STR(NC_EIO)
        ERR_CODE_STR(NC_ENODATA)
        ERR_CODE_STR(NC_EDAPSVC)
        ERR_CODE_STR(NC_EDAS)
        ERR_CODE_STR(NC_EDDS)
        ERR_CODE_STR(NC_EDATADDS)
        ERR_CODE_STR(NC_EDAPURL)
        ERR_CODE_STR(NC_EDAPCONSTRAINT)
        ERR_CODE_STR(NC_ETRANSLATION)
        ERR_CODE_STR(NC_EACCESS)
        ERR_CODE_STR(NC_EAUTH)
        ERR_CODE_STR(NC_ENOTFOUND)
        ERR_CODE_STR(NC_ECANTREMOVE)
        ERR_CODE_STR(NC_EHDFERR)
        ERR_CODE_STR(NC_ECANTREAD)
        ERR_CODE_STR(NC_ECANTWRITE)
        ERR_CODE_STR(NC_ECANTCREATE)
        ERR_CODE_STR(NC_EFILEMETA)
        ERR_CODE_STR(NC_EDIMMETA)
        ERR_CODE_STR(NC_EATTMETA)
        ERR_CODE_STR(NC_EVARMETA)
        ERR_CODE_STR(NC_ENOCOMPOUND)
        ERR_CODE_STR(NC_EATTEXISTS)
        ERR_CODE_STR(NC_ENOTNC4)
        ERR_CODE_STR(NC_ESTRICTNC3)
        ERR_CODE_STR(NC_ENOTNC3)
        ERR_CODE_STR(NC_ENOPAR)
        ERR_CODE_STR(NC_EPARINIT)
        ERR_CODE_STR(NC_EBADGRPID)
        ERR_CODE_STR(NC_EBADTYPID)
        ERR_CODE_STR(NC_ETYPDEFINED)
        ERR_CODE_STR(NC_EBADFIELD)
        ERR_CODE_STR(NC_EBADCLASS)
        ERR_CODE_STR(NC_EMAPTYPE)
        ERR_CODE_STR(NC_ELATEFILL)
        ERR_CODE_STR(NC_ELATEDEF)
        ERR_CODE_STR(NC_EDIMSCALE)
        ERR_CODE_STR(NC_ENOGRP)
        ERR_CODE_STR(NC_ESTORAGE)
        ERR_CODE_STR(NC_EBADCHUNK)
        ERR_CODE_STR(NC_ENOTBUILT)
        ERR_CODE_STR(NC_EDISKLESS)
        ERR_CODE_STR(NC_ECANTEXTEND)
        ERR_CODE_STR(NC_EMPI)
        // ERR_CODE_STR(NC_EURL)
        // ERR_CODE_STR(NC_ECONSTRAINT)
        ERR_CODE_STR(NC_ESMALL)
        ERR_CODE_STR(NC_ENOTINDEP)
        ERR_CODE_STR(NC_EINDEP)
        ERR_CODE_STR(NC_EFILE)
        ERR_CODE_STR(NC_EREAD)
        ERR_CODE_STR(NC_EWRITE)
        ERR_CODE_STR(NC_EOFILE)
        ERR_CODE_STR(NC_EMULTITYPES)
        ERR_CODE_STR(NC_EIOMISMATCH)
        ERR_CODE_STR(NC_ENEGATIVECNT)
        ERR_CODE_STR(NC_EUNSPTETYPE)
        ERR_CODE_STR(NC_EINVAL_REQUEST)
        ERR_CODE_STR(NC_EAINT_TOO_SMALL)
        ERR_CODE_STR(NC_ENOTSUPPORT)
        ERR_CODE_STR(NC_ENULLBUF)
        ERR_CODE_STR(NC_EPREVATTACHBUF)
        ERR_CODE_STR(NC_ENULLABUF)
        ERR_CODE_STR(NC_EPENDINGBPUT)
        ERR_CODE_STR(NC_EINSUFFBUF)
        ERR_CODE_STR(NC_ENOENT)
        ERR_CODE_STR(NC_EINTOVERFLOW)
        ERR_CODE_STR(NC_ENOTENABLED)
        ERR_CODE_STR(NC_EBAD_FILE)
        ERR_CODE_STR(NC_ENO_SPACE)
        ERR_CODE_STR(NC_EQUOTA)
        ERR_CODE_STR(NC_ENULLSTART)
        ERR_CODE_STR(NC_ENULLCOUNT)
        ERR_CODE_STR(NC_EINVAL_CMODE)
        ERR_CODE_STR(NC_EINVAL_OMODE)
        ERR_CODE_STR(NC_ETYPESIZE)
        ERR_CODE_STR(NC_ETYPE_MISMATCH)
        ERR_CODE_STR(NC_ETYPESIZE_MISMATCH)
        ERR_CODE_STR(NC_ESTRICTCDF2)
        ERR_CODE_STR(NC_ENOTRECVAR)
        ERR_CODE_STR(NC_ENOTFILL)
        ERR_CODE_STR(NC_EMULTIDEFINE)
        ERR_CODE_STR(NC_EMULTIDEFINE_OMODE)
        ERR_CODE_STR(NC_EMULTIDEFINE_CMODE)
        ERR_CODE_STR(NC_EMULTIDEFINE_DIM_NUM)
        ERR_CODE_STR(NC_EMULTIDEFINE_DIM_SIZE)
        ERR_CODE_STR(NC_EMULTIDEFINE_DIM_NAME)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_NUM)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_NAME)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_NDIMS)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_DIMIDS)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_TYPE)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_LEN)
        ERR_CODE_STR(NC_EMULTIDEFINE_NUMRECS)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_BEGIN)
        ERR_CODE_STR(NC_EMULTIDEFINE_ATTR_NUM)
        ERR_CODE_STR(NC_EMULTIDEFINE_ATTR_SIZE)
        ERR_CODE_STR(NC_EMULTIDEFINE_ATTR_NAME)
        ERR_CODE_STR(NC_EMULTIDEFINE_ATTR_TYPE)
        ERR_CODE_STR(NC_EMULTIDEFINE_ATTR_LEN)
        ERR_CODE_STR(NC_EMULTIDEFINE_ATTR_VAL)
        ERR_CODE_STR(NC_EMULTIDEFINE_FNC_ARGS)
        ERR_CODE_STR(NC_EMULTIDEFINE_FILL_MODE)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_FILL_MODE)
        ERR_CODE_STR(NC_EMULTIDEFINE_VAR_FILL_VALUE)
        default:
              sprintf(unknown_str,"Unknown code %d",err);
    }
    return unknown_str;
}

/*----< inq_env_hint() >-----------------------------------------------------*/
int
inq_env_hint(char *hint_key, char **hint_value)
{
    char *warn_str="Warning: skip ill-formed hint set in PNETCDF_HINTS";
    char *env_str;

    /* read hints set in the environment variable PNETCDF_HINTS, a string of
     * hints separated by ";" and each hint is in the form of hint=value. E.g.
     * "cb_nodes=16;cb_config_list=*:6". If this environment variable is set,
     * this subroutine allocates char array for hint_value, copy the hint
     * value to it, and return 1. Otherwise it returns 0 with *value set to
     * NULL.
     */

    *hint_value = NULL;

    /* get environment variable PNETCDF_HINTS */
    if ((env_str = getenv("PNETCDF_HINTS")) != NULL) {
        char *env_str_cpy, *hint, *next_hint, *key, *val, *deli;
        char *hint_saved=NULL;

        env_str_cpy = strdup(env_str);
        next_hint = env_str_cpy;

        do {
            hint = next_hint;
            deli = strchr(hint, ';');
            if (deli != NULL) {
                *deli = '\0'; /* add terminate char */
                next_hint = deli + 1;
            }
            else next_hint = "\0";
            if (hint_saved != NULL) free(hint_saved);

            /* skip all-blank hint */
            hint_saved = strdup(hint);
            if (strtok(hint, " \t") == NULL) continue;

            free(hint_saved);
            hint_saved = strdup(hint); /* save hint for error message */

            deli = strchr(hint, '=');
            if (deli == NULL) { /* ill-formed hint */
                printf("%s: '%s'\n", warn_str, hint_saved);
                continue;
            }
            *deli = '\0';

            /* hint key */
            key = strtok(hint, "= \t");
            if (key == NULL || NULL != strtok(NULL, "= \t")) {
                /* expect one token before = */
                printf("%s: '%s'\n", warn_str, hint_saved);
                continue;
            }

            /* hint value */
            val = strtok(deli+1, "= \t");
            if (NULL != strtok(NULL, "= \t")) { /* expect one token before = */
                printf("%s: '%s'\n", warn_str, hint_saved);
                continue;
            }
            if (strcasecmp(key,hint_key) == 0) {
                /* inquired hint is found */
                if (val != NULL) {
                    *hint_value = (char*) malloc(strlen(val)+1);
                    strcpy(*hint_value, val);
                }
                if (hint_saved != NULL) free(hint_saved);
                free(env_str_cpy);
                return (val == NULL) ? 0 : 1;
            }
        } while (*next_hint != '\0');

        if (hint_saved != NULL) free(hint_saved);
        free(env_str_cpy);
    }
    return 0;
}

/* File system types recognized by ROMIO in MPICH 4.0.0 */
static const char* fstypes[] = {"ufs", "nfs", "xfs", "pvfs2", "gpfs", "panfs", "lustre", "daos", "testfs", "ime", "quobyte", NULL};

/* Return a pointer to filename by removing the file system type prefix name if
 * there is any.  For example, when filename = "lustre:/home/foo/testfile.nc",
 * remove "lustre:" to return a pointer to "/home/foo/testfile.nc", so the name
 * can be used in POSIX open() calls.
 */
char* remove_file_system_type_prefix(const char *filename)
{
    char *ret_filename = (char*)filename;

    if (filename == NULL) return NULL;

    if (strchr(filename, ':') != NULL) { /* there is a prefix end with ':' */
        /* check if prefix is one of recognized file system types */
        int i=0;
        while (fstypes[i] != NULL) {
            size_t prefix_len = strlen(fstypes[i]);
            if (!strncmp(filename, fstypes[i], prefix_len)) { /* found */
                ret_filename += prefix_len + 1;
                break;
            }
            i++;
        }
    }

    return ret_filename;
}

int is_relax_coord_bound(void)
{
    char *env_str;
    int relax_coord_bound;

#if defined(PNETCDF_RELAX_COORD_BOUND) && PNETCDF_RELAX_COORD_BOUND == 1
    relax_coord_bound = 1;
#else
    relax_coord_bound = 0;
#endif
    if ((env_str = getenv("PNETCDF_RELAX_COORD_BOUND")) != NULL) {
        /* the env variable is set */
        if (*env_str == '0') relax_coord_bound = 0;
        else                 relax_coord_bound = 1;
    }

    return relax_coord_bound;
}

void
static tst_main_usage(char *argv0)
{
    char *base_name = basename(argv0);
    char *help =
    "Usage: %s [OPTIONS]...[filename]\n"
    "       [-h] Print help\n"
    "       [-q] quiet mode\n"
    "       [-k] Keep output files (default: no)\n"
    "       [-i  in_path]: input file path (default: NULL)\n"
    "       [-o out_path]: output netCDF file name (default: %s.nc)\n";
    fprintf(stderr, help, base_name, base_name);
}

int tst_main(int        argc,
             char      **argv,
             char       *msg,  /* short description about the test */
             loop_opts   opt,  /* test options */
             int       (*tst_body)(const char*,const char*,int,int,MPI_Info))
{
    extern int optind;
    extern char *optarg;
    char *in_path=NULL, *out_path=NULL, *out_dir;

    /* IDs for the netCDF file, dimensions, and variables. */
    int nprocs, rank, err, nerrs=0, keep_files, quiet;
    int i, a, d, b, s, m, u;
    int s_ina, e_ina, s_drv, e_drv, s_bb, e_bb, s_mod, e_mod, s_ds, e_ds;
    int s_ibuf, e_ibuf;

    MPI_Info info=MPI_INFO_NULL;
    double timing = MPI_Wtime();

#if defined(PNETCDF_PROFILING) && PNETCDF_PROFILING == 1
    double itiming[256]; int k=0;
#endif

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    keep_files = 0;
    quiet = 0;

    while ((i = getopt(argc, argv, "hqki:o:")) != EOF)
        switch(i) {
            case 'q':
                quiet = 1;
                break;
            case 'k':
                keep_files = 1;
                break;
            case 'i':
                in_path = strdup(optarg);
                break;
            case 'o':
                out_path = strdup(optarg);
                break;
            case 'h':
            default:  if (rank==0) tst_main_usage(argv[0]);
                      MPI_Finalize();
                      exit(1);
        }

#ifndef TESTOUTDIR
    out_dir = getenv("TESTOUTDIR");
    if (out_dir == NULL) out_dir = ".";
#else
    out_dir = TESTOUTDIR;
#endif

    if (out_path == NULL)
        out_path = strdup("testfile.nc");
#if 0
    else {
        /* check if filename is a directory */
        struct stat sb;

        snprintf(filename, 256, "%s", argv[optind]);
        if (stat(filename, &sb) == 0 && S_ISDIR(sb.st_mode))
            append_suffix = 0;
    }
#endif

    if (rank == 0) {
        char *cmd_str = (char *)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s - %s", basename(argv[0]), msg);
        printf("%-63s -- ", cmd_str);
        free(cmd_str);
    }

    char cmd_opts[64];
    sprintf(cmd_opts, "Rank %d: ncmpidiff", rank);

    char *ptr = strrchr(out_path, '.');
    if (ptr != NULL) *ptr = '\0';

    MPI_Info_create(&info);

    /* Set common I/O hints. Using smaller values for hints below can make
     * the tests more rigorous.
     */
    MPI_Info_set(info, "ind_wr_buffer_size", "60");
    MPI_Info_set(info, "ind_rd_buffer_size", "70");
    MPI_Info_set(info, "cb_buffer_size", "500");

#define SET_OPT(key) {               \
    if (opt.key == 2) {              \
        s_ ## key = 1; e_##key = 0;  \
    }                                \
    else if (opt.key == 1) {         \
        s_##key = 1; e_##key = 1;    \
    }                                \
    else { /* #key == 0 */           \
        s_##key = 0; e_##key = 0;    \
    }                                \
}

    SET_OPT(ina)    /* test intra-node aggregation */
    SET_OPT(drv)    /* test MPI-IO, GIO drivers */
    SET_OPT(ibuf)   /* test hint nc_ibuf_size */
    SET_OPT(mod)    /* test collective/independent data mode */
    SET_OPT(bb)     /* test of burst-buffering feature */

#if !defined(ENABLE_GIO) || ENABLE_GIO == 0
    s_drv = e_drv = 1; /* skip testing GIO driver */
#endif

#if !defined(PNETCDF_BURST_BUFFERING) || PNETCDF_BURST_BUFFERING == 0
    s_bb = e_bb = 0; /* skip testing burst buffering */
#endif

    s_ds = 1; e_ds = 0; /* test both date sieving enabled and disabled */

    for (i=0; i<opt.num_fmts; i++) {
        char out_filename[512], ext[16], *base_file;

        base_file = NULL;
        if (opt.formats[i] > 0)
            sprintf(ext, "nc%d", opt.formats[i]);
        else /* for tests that are not testing different CDF versions */
            strcpy(ext, "nc");

        /* For indice a,d,b,s,m, 0 is PnetCDF's default setting */
        for (a=s_ina;  a>=e_ina;  a--) {
        for (d=s_drv;  d>=e_drv;  d--) {
        for (u=s_ibuf; u>=e_ibuf; u--) {
        for (b=s_bb;   b>=e_bb;   b--) {
        for (s=s_ds;   s>=e_ds;   s--) {
        for (m=s_mod;  m>=e_mod;  m--) {

#if defined(PNETCDF_PROFILING) && PNETCDF_PROFILING == 1
            MPI_Barrier(MPI_COMM_WORLD);
            itiming[k] = MPI_Wtime();
#endif

            sprintf(out_filename, "%s.%s", out_path, ext);

            /* Whether or not to enable the intra-node aggregation */
            if (a == 0) { /* PnetCDF's default */
                MPI_Info_set(info, "nc_num_aggrs_per_node", "0");
                strcat(out_filename, ".noina");
            } else {
                MPI_Info_set(info, "nc_num_aggrs_per_node", "2");
                strcat(out_filename, ".ina");
            }

            /* Use MPI-IO or GIO driver */
            if (d == 0) { /* PnetCDF's default */
                MPI_Info_set(info, "nc_driver", "gio");
                strcat(out_filename, ".gio");
            } else {
                MPI_Info_set(info, "nc_driver", "mpiio");
                strcat(out_filename, ".mpio");
            }

            /* Set hint nc_ibuf_size */
            if (u == 0) { /* test nc_ibuf_size = 1 MB */
                /* PnetCDF's default PNC_DEFAULT_IBUF_SIZE 16777216 */
                MPI_Info_set(info, "nc_ibuf_size", "16777216");
                strcat(out_filename, ".ibuf");
            } else { /* do not use internal buffering */
                MPI_Info_set(info, "nc_ibuf_size", "0");
                strcat(out_filename, ".noibuf");
            }

            /* Whether or not to enable the burst buffering */
            if (b == 0) { /* PnetCDF's default */
                MPI_Info_set(info, "nc_burst_buf", "disable");
                strcat(out_filename, ".nobb");
            }
            else {
                MPI_Info_set(info, "nc_burst_buf", "enable");
                MPI_Info_set(info, "nc_burst_buf_dirname", out_dir);
                MPI_Info_set(info, "nc_burst_buf_overwrite", "enable");
                strcat(out_filename, ".bb");
            }

            /* Whether or not to enable data sieving */
            if (s == 0) { /* PnetCDF's default */
                MPI_Info_set(info, "romio_ds_read",  "automatic");
                MPI_Info_set(info, "romio_ds_write", "automatic");
                strcat(out_filename, ".ds");
            }
            else {
                MPI_Info_set(info, "romio_ds_read",  "disable");
                MPI_Info_set(info, "romio_ds_write", "disable");
                strcat(out_filename, ".nods");
            }

            /* Test PnetCDF collective or independent APIs */
            if (m == 0) /* collective data mode */
                strcat(out_filename, ".coll_mod");
            else /* independent data mode */
                strcat(out_filename, ".indep_mod");

            /* NetCDF4 does not allow to extend number of record numbers in
             * independent data mode. NC_ECANTEXTEND will be returned.
             */
            if (m == 1 &&
                (opt.formats[i] == NC_FORMAT_NETCDF4_CLASSIC ||
                 opt.formats[i] == NC_FORMAT_NETCDF4))
                continue;

#define RUN_ERR(msg, fname) {\
    printf("\n%s %-44s (INA=%s drv=%s ibuf=%s BB=%s DS=%s mode=%s)\n", \
           msg, fname, (a)?"yes":"no", (d)?"mpio":"gio", (u)?"no":"yes", \
           (b)?"yes":"no", (s)?"no":"auto", (m)?"indp":"coll"); \
}

#define DIFF_ERR(msg, f1, f2) {\
    printf("\n%s %-44s %s (INA=%s drv=%s ibuf=%s BB=%s DS=%s mode=%s)\n", \
           msg, f1, f2, (a)?"yes":"no", (d)?"mpio":"gio", (u)?"no":"yes", \
           (b)?"yes":"no", (s)?"no":"auto", (m)?"indp":"coll"); \
}

            double time_body = MPI_Wtime();
            if (!quiet && rank == 0)
                RUN_ERR("Testing", out_filename)

            /* tst_body() is the core of test program */
            int coll_io = (m == 0);
            nerrs = tst_body(out_filename, in_path, opt.formats[i], coll_io,
                             info);
            MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_MAX,
                          MPI_COMM_WORLD);
            if (nerrs > 0) {
                fflush(stdout);
                if (rank == 0)
                    RUN_ERR("\nFAILED", out_filename)
                goto err_out;
            }

            if (!quiet) {
                time_body = MPI_Wtime() - time_body;
                MPI_Allreduce(MPI_IN_PLACE, &time_body, 1, MPI_DOUBLE, MPI_MAX,
                              MPI_COMM_WORLD);
                if (rank == 0)
                    printf(" (%.2fs)\n", time_body);
            }

#if defined(PNETCDF_PROFILING) && PNETCDF_PROFILING == 1
            itiming[k] = MPI_Wtime() - itiming[k]; k++;
#endif
            /* wait for all processes to complete */
            MPI_Barrier(MPI_COMM_WORLD);

            /* run ncmpidiff to compare output files */
            if (base_file == NULL) { /* skip first file */
                base_file = strdup(out_filename);
                goto skip_diff;
            }

            if (!opt.hdr_diff) goto skip_diff;

            if (strcmp(base_file, out_filename)) {
                int check_header=1, check_entire_file, first_diff=1;

                check_entire_file = (opt.var_diff) ? 1 : 0;

                /* ncmpidiff does not support netCDF4 files */
                if (opt.formats[i] == NC_FORMAT_NETCDF4_CLASSIC ||
                    opt.formats[i] == NC_FORMAT_NETCDF4) {
                    goto skip_diff;
                }

                if (!quiet && rank == 0)
                    DIFF_ERR("ncmpidiff", out_filename, base_file)

#ifdef MIMIC_LUSTRE
                /* use a larger stripe size when running ncmpidiff is faster */
                setenv("MIMIC_STRIPE_SIZE", "1048576", 1);
#else
                unsetenv("MIMIC_STRIPE_SIZE");
#endif
                /* running ncmpidiff also validates the file header */
                MPI_Offset numDIFF;
                numDIFF = ncmpidiff_core(out_filename, base_file,
                                         MPI_COMM_WORLD, info, 0,
                                         quiet, check_header, 0,
                                         check_entire_file, 0, NULL, 0,
                                         first_diff, cmd_opts, 0, 0);

                /* check error so it can error out collectively */
                MPI_Allreduce(MPI_IN_PLACE, &numDIFF, 1, MPI_OFFSET, MPI_MAX,
                              MPI_COMM_WORLD);
                if (numDIFF > 0) {
                    if (rank == 0)
                        DIFF_ERR("FAILED: ncmpidiff", out_filename, base_file)

                    nerrs = 1;
                    goto err_out;
                }

                /* wait for all ranks to complete diff before file delete */
                MPI_Barrier(MPI_COMM_WORLD);
            }
skip_diff:
            if (!keep_files && base_file != NULL && strcmp(base_file, out_filename)) {
                if (rank == 0)
                    ncmpi_delete(out_filename, MPI_INFO_NULL);

                /* wait for deletion to complete before next iteration */
                MPI_Barrier(MPI_COMM_WORLD);
            }
        } /* loop m */
        } /* loop s */
        } /* loop b */
        } /* loop u */
        } /* loop d */
        } /* loop a */

        if (base_file != NULL) {
            if (!keep_files) {
                if (rank == 0)
                    ncmpi_delete(base_file, MPI_INFO_NULL);

                /* all ranks wait for root to complete file deletion */
                MPI_Barrier(MPI_COMM_WORLD);
            }
            free(base_file);
        }
    }
    MPI_Info_free(&info);

    /* check if there is any malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has "OFFFMT" bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

err_out:
    if (in_path  != NULL) free(in_path);
    if (out_path != NULL) free(out_path);

    timing = MPI_Wtime() - timing;
    MPI_Allreduce(MPI_IN_PLACE, &timing, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#if defined(PNETCDF_PROFILING) && PNETCDF_PROFILING == 1
    MPI_Allreduce(MPI_IN_PLACE, itiming, 256, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0 && !quiet) {
        for (i=0; i<k; i++) printf("k=%d timing[%3d]=%.4f\n",k,i,itiming[i]);
        printf("\n");
    }
#endif

    if (rank == 0) {
        if (nerrs)
            printf(FAIL_STR, nerrs);
        else
            printf(PASS_STR, timing);
    }

    return (nerrs > 0);
}

/*----< pnc_fmt_string() >---------------------------------------------------*/
char* pnc_fmt_string(int format)
{
    switch(format) {
        case NC_FORMAT_CLASSIC:         return "NC_FORMAT_CLASSIC";
        case NC_FORMAT_64BIT_OFFSET:    return "NC_FORMAT_64BIT_OFFSET";
        case NC_FORMAT_NETCDF4:         return "NC_FORMAT_NETCDF4";
        case NC_FORMAT_NETCDF4_CLASSIC: return "NC_FORMAT_NETCDF4_CLASSIC";
        case NC_FORMAT_64BIT_DATA:      return "NC_FORMAT_64BIT_DATA";
        default:                        return "UNKNOWN";
    }
}
