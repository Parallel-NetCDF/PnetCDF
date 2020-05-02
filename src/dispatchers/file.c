/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>     /* getenv() */
#include <string.h>     /* strtok(), strtok_r(), strchr(), strcpy(), strdup() */
#include <strings.h>    /* strcasecmp() */
#include <fcntl.h>      /* open() */
#include <sys/types.h>  /* lseek() */
#include <unistd.h>     /* read(), close(), lseek() */
#include <assert.h>     /* assert() */
#include <errno.h>      /* errno */

#ifdef ENABLE_THREAD_SAFE
#include<pthread.h>
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef ENABLE_NETCDF4
#include <netcdf.h>
#endif

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

#ifdef ENABLE_ADIOS
#include "adios_read.h"
#include <arpa/inet.h>
#define BP_MINIFOOTER_SIZE 28
#define BUFREAD64(buf,var) memcpy(&var, buf, 8); if (diff_endian) swap_64(&var);
#endif

/* TODO: the following 3 global variables make PnetCDF not thread safe */

/* static variables are initialized to NULLs */
static PNC *pnc_filelist[NC_MAX_NFILES];
static int  pnc_numfiles;

/* This is the default create format for ncmpi_create and nc__create.
 * The use of this file scope variable is not thread-safe.
 */
static int ncmpi_default_create_format = NC_FORMAT_CLASSIC;

#define NCMPII_HANDLE_ERROR(func) \
    if (mpireturn != MPI_SUCCESS) { \
        int errorStringLen; \
        char errorString[MPI_MAX_ERROR_STRING]; \
        MPI_Error_string(mpireturn, errorString, &errorStringLen); \
        printf("%s error at line %d file %s (%s)\n", func, __LINE__, __FILE__, errorString); \
    }

/*----< new_id_PNCList() >---------------------------------------------------*/
/* Return a new ID (array index) from the PNC list, pnc_filelist[] that is
 * not used. Note the used elements in pnc_filelist[] may not be contiguous.
 * For example, some files created/opened later may be closed earlier than
 * others, leaving those array elements NULL in the middle.
 */
static int
new_id_PNCList(int *new_id, PNC *pncp)
{
    int i, err;

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif
    *new_id = -1;
    if (pnc_numfiles == NC_MAX_NFILES) { /* Too many files open */
        DEBUG_ASSIGN_ERROR(err, NC_ENFILE)
    }
    else {
        err = NC_NOERR;
        for (i=0; i<NC_MAX_NFILES; i++) { /* find the first unused element */
            if (pnc_filelist[i] == NULL) {
                *new_id = i;
                pnc_filelist[i] = pncp;
                pnc_numfiles++; /* increment number of files opened */
                break;
            }
        }
    }
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif

    return err;
}

/*----< del_from_PNCList() >-------------------------------------------------*/
static void
del_from_PNCList(int ncid)
{
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif

    /* validity of ncid should have been checked already */
    pnc_filelist[ncid] = NULL;
    pnc_numfiles--;

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif
}

/*----< PNC_check_id() >-----------------------------------------------------*/
int
PNC_check_id(int ncid, PNC **pncp)
{
    int err=NC_NOERR;

    assert(pncp != NULL);

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif

    if (pnc_numfiles == 0 || ncid < 0 || ncid >= NC_MAX_NFILES)
        DEBUG_ASSIGN_ERROR(err, NC_EBADID)
    else
        *pncp = pnc_filelist[ncid];

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif

    return err;
}

/*----< construct_info() >---------------------------------------------------*/
static void
combine_env_hints(MPI_Info  user_info,  /* IN */
                  MPI_Info *new_info)   /* OUT: may be MPI_INFO_NULL */
{
    char *warn_str="Warning: skip ill-formed hint set in PNETCDF_HINTS";
    char *env_str;

    /* take hints from the environment variable PNETCDF_HINTS, a string of
     * hints separated by ";" and each hint is in the form of hint=value. E.g.
     * "cb_nodes=16;cb_config_list=*:6". If this environment variable is set,
     * it overrides the same hints that were set by MPI_Info_set() called in
     * the application program.
     */
    if (user_info != MPI_INFO_NULL)
        MPI_Info_dup(user_info, new_info); /* ignore error */
    else
        *new_info = MPI_INFO_NULL;

    /* get environment variable PNETCDF_HINTS */
    if ((env_str = getenv("PNETCDF_HINTS")) != NULL) {
#ifdef USE_STRTOK_R
        char *env_str_cpy, *env_str_saved, *hint, *key;
        env_str_cpy = strdup(env_str);
        env_str_saved = env_str_cpy;
        hint = strtok_r(env_str_cpy, ";", &env_str_saved);
        while (hint != NULL) {
            char *hint_saved = strdup(hint);
            char *val = strchr(hint, '=');
            if (val == NULL) { /* ill-formed hint */
                if (NULL != strtok(hint, " \t"))
                    printf("%s: '%s'\n", warn_str, hint_saved);
                /* else case: ignore white-spaced hints */
                free(hint_saved);
                hint = strtok_r(NULL, ";", &env_str_saved); /* get next hint */
                continue;
            }
            key = strtok(hint, "= \t");
            val = strtok(NULL, "= \t");
            if (NULL != strtok(NULL, "= \t")) /* expect no more token */
                printf("%s: '%s'\n", warn_str, hint_saved);
            else {
                if (*new_info == MPI_INFO_NULL)
                    MPI_Info_create(new_info); /* ignore error */
                MPI_Info_set(*new_info, key, val); /* override or add */
            }
            /* printf("env hint: key=%s val=%s\n",key,val); */
            hint = strtok_r(NULL, ";", &env_str_saved);
            free(hint_saved);
        }
        free(env_str_cpy);
#else
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
            if (*new_info == MPI_INFO_NULL)
                MPI_Info_create(new_info); /* ignore error */
            MPI_Info_set(*new_info, key, val); /* override or add */

        } while (*next_hint != '\0');

        if (hint_saved != NULL) free(hint_saved);
        free(env_str_cpy);
#endif
    }
    /* return no error as all hints are advisory */
}

/*----< ncmpi_create() >-----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_create(MPI_Comm    comm,
             const char *path,
             int         cmode,
             MPI_Info    info,
             int        *ncidp)
{
    int rank, nprocs, status=NC_NOERR, err;
    int safe_mode=0, mpireturn, relax_coord_bound, format;
    char *env_str;
    MPI_Info combined_info;
    void *ncp;
    PNC *pncp;
    PNC_driver *driver;
#ifdef BUILD_DRIVER_FOO
    int enable_foo_driver=0;
#endif
#ifdef ENABLE_BURST_BUFFER
    int enable_bb_driver=0;
#endif

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

#ifdef PNETCDF_DEBUG
    safe_mode = 1;
    /* When debug mode is enabled at the configure time, safe_mode is by
     * default enabled. This can be overwritten by the run-time environment
     * variable PNETCDF_SAFE_MODE */
#endif
    /* get environment variable PNETCDF_SAFE_MODE
     * if it is set to 1, then we perform a strict parameter consistent test
     */
    if ((env_str = getenv("PNETCDF_SAFE_MODE")) != NULL) {
        if (*env_str == '0') safe_mode = 0;
        else                 safe_mode = 1;
        /* if PNETCDF_SAFE_MODE is set but without a value, *env_str can
         * be '\0' (null character). In this case, safe_mode is enabled */
    }

    /* get environment variable PNETCDF_RELAX_COORD_BOUND
     * if it is set to 0, then we perform a strict start bound check
     */
#ifndef RELAX_COORD_BOUND
    relax_coord_bound = 0;
#else
    relax_coord_bound = 1;
#endif
    if ((env_str = getenv("PNETCDF_RELAX_COORD_BOUND")) != NULL) {
        if (*env_str == '0') relax_coord_bound = 0;
        else                 relax_coord_bound = 1;
        /* if PNETCDF_RELAX_COORD_BOUND is set but without a value, *env_str
         * can be '\0' (null character). This is equivalent to setting
         * relax_coord_bound to 1 */
    }

    /* path's validity is checked in MPI-IO with error code MPI_ERR_BAD_FILE
     * path consistency is checked in MPI-IO with error code MPI_ERR_NOT_SAME
     */
    if (path == NULL || *path == '\0') DEBUG_RETURN_ERROR(NC_EBAD_FILE)

    if (nprocs > 1) { /* Check cmode consistency */
        int root_cmode = cmode; /* only root's matters */
        TRACE_COMM(MPI_Bcast)(&root_cmode, 1, MPI_INT, 0, comm);
        NCMPII_HANDLE_ERROR("MPI_Bcast")

        /* Overwrite cmode with root's cmode */
        if (root_cmode != cmode) {
            cmode = root_cmode;
            DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE_CMODE)
        }

        if (safe_mode) { /* sync status among all processes */
            err = status;
            TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, comm);
            NCMPII_HANDLE_ERROR("MPI_Allreduce")
        }
        /* continue to use root's cmode to create the file, but will report
         * cmode inconsistency error, if there is any */
    }

    /* combine user's MPI info and PNETCDF_HINTS env variable */
    combine_env_hints(info, &combined_info);

#ifdef BUILD_DRIVER_FOO
    if (combined_info == MPI_INFO_NULL)
        MPI_Info_create(&combined_info);
    {
        char value[MPI_MAX_INFO_VAL];
        int flag;

        /* check if nc_foo_driver is enabled */
        MPI_Info_get(combined_info, "nc_foo_driver", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag && strcasecmp(value, "enable") == 0)
            enable_foo_driver = 1;
    }
#endif
#ifdef ENABLE_BURST_BUFFER
    if (combined_info == MPI_INFO_NULL)
        MPI_Info_create(&combined_info);
    {
        char value[MPI_MAX_INFO_VAL];
        int flag;

        /* check if nc_burst_buf is enabled */
        MPI_Info_get(combined_info, "nc_burst_buf", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag && strcasecmp(value, "enable") == 0)
            enable_bb_driver = 1;
    }
#endif

    /* Use environment variable and cmode to tell the file format
     * which is later used to select the right driver.
     */

#ifdef ENABLE_NETCDF4
    /* It is illegal to have NC_64BIT_OFFSET & NC_64BIT_DATA & NC_NETCDF4 */
    if ((cmode & (NC_64BIT_OFFSET|NC_NETCDF4)) ==
                 (NC_64BIT_OFFSET|NC_NETCDF4) ||
        (cmode & (NC_64BIT_DATA|NC_NETCDF4)) ==
                 (NC_64BIT_DATA|NC_NETCDF4)) {
        if (combined_info != MPI_INFO_NULL)
            MPI_Info_free(&combined_info);
        DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)
    }
#else
    if (cmode & NC_NETCDF4) {
        if (combined_info != MPI_INFO_NULL)
            MPI_Info_free(&combined_info);
        DEBUG_RETURN_ERROR(NC_ENOTBUILT)
    }
#endif

    /* It is illegal to have both NC_64BIT_OFFSET & NC_64BIT_DATA */
    if ((cmode & (NC_64BIT_OFFSET|NC_64BIT_DATA)) ==
                 (NC_64BIT_OFFSET|NC_64BIT_DATA)) {
        if (combined_info != MPI_INFO_NULL)
            MPI_Info_free(&combined_info);
        DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)
    }

    /* Check if cmode contains format specific flag */
    if (fIsSet(cmode, NC_64BIT_DATA))
        format = NC_FORMAT_CDF5;
    else if (fIsSet(cmode, NC_64BIT_OFFSET))
        format = NC_FORMAT_CDF2;
    else if (fIsSet(cmode, NC_NETCDF4)) {
        if (fIsSet(cmode, NC_CLASSIC_MODEL))
            format = NC_FORMAT_NETCDF4_CLASSIC;
        else
            format = NC_FORMAT_NETCDF4;
    }
    else if (fIsSet(cmode, NC_CLASSIC_MODEL))
        format = NC_FORMAT_CLASSIC;
    else {
        /* if no file format flag is set in cmode, use default */
        ncmpi_inq_default_format(&format);
        if (format == NC_FORMAT_CDF5)
            cmode |= NC_64BIT_DATA;
        else if (format == NC_FORMAT_CDF2)
            cmode |= NC_64BIT_OFFSET;
        else if (format == NC_FORMAT_NETCDF4)
            cmode |= NC_NETCDF4;
        else if (format == NC_FORMAT_NETCDF4_CLASSIC)
            cmode |= NC_NETCDF4 | NC_CLASSIC_MODEL;
    }

#ifdef ENABLE_NETCDF4
    if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC) {
        driver = nc4io_inq_driver();
#ifdef ENABLE_BURST_BUFFER
        /* Burst buffering does not support NetCDF-4 files yet.
         * If hint nc_burst_buf is enabled in combined_info, disable it.
         */
        if (enable_bb_driver == 1)
            MPI_Info_set(combined_info, "nc_burst_buf", "disable");
        enable_bb_driver = 0;
#endif
    }
    else
#endif
#ifdef BUILD_DRIVER_FOO
    if (enable_foo_driver)
        driver = ncfoo_inq_driver();
    else
#endif
#ifdef ENABLE_BURST_BUFFER
    if (enable_bb_driver)
        driver = ncbbio_inq_driver();
    else
#endif
        /* default is the driver built on top of MPI-IO */
        driver = ncmpio_inq_driver();

    /* allocate a new PNC object */
    pncp = (PNC*) NCI_Malloc(sizeof(PNC));
    if (pncp == NULL) {
        *ncidp = -1;
        if (combined_info != MPI_INFO_NULL)
            MPI_Info_free(&combined_info);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    /* generate a new nc file ID from NCPList */
    err = new_id_PNCList(ncidp, pncp);
    if (err != NC_NOERR) {
        if (combined_info != MPI_INFO_NULL)
            MPI_Info_free(&combined_info);
        return err;
    }

    /* Duplicate comm, because users may free it (though unlikely). Note
     * MPI_Comm_dup() is collective. We pass pncp->comm to drivers, so there
     * is no need for a driver to duplicate it again.
     */
    if (comm != MPI_COMM_WORLD && comm != MPI_COMM_SELF)
        MPI_Comm_dup(comm, &pncp->comm);
    else
        pncp->comm = comm;

    /* calling the driver's create subroutine */
    err = driver->create(pncp->comm, path, cmode, *ncidp, combined_info, &ncp);
    if (status == NC_NOERR) status = err;
    if (combined_info != MPI_INFO_NULL) MPI_Info_free(&combined_info);
    if (status != NC_NOERR && status != NC_EMULTIDEFINE_CMODE) {
        del_from_PNCList(*ncidp);
        if (pncp->comm != MPI_COMM_WORLD && pncp->comm != MPI_COMM_SELF)
            MPI_Comm_free(&pncp->comm); /* a collective call */
        NCI_Free(pncp);
        *ncidp = -1;
        return status;
    }

    /* fill in pncp members */
    pncp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (pncp->path == NULL) {
        driver->close(ncp); /* close file and ignore error */
        del_from_PNCList(*ncidp);
        if (pncp->comm != MPI_COMM_WORLD && pncp->comm != MPI_COMM_SELF)
            MPI_Comm_free(&pncp->comm); /* a collective call */
        NCI_Free(pncp);
        *ncidp = -1;
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(pncp->path, path);
    pncp->mode       = cmode;
    pncp->driver     = driver;
    pncp->ndims      = 0;
    pncp->unlimdimid = -1;
    pncp->nvars      = 0;
    pncp->nrec_vars  = 0;
    pncp->vars       = NULL;
    pncp->flag       = NC_MODE_DEF | NC_MODE_CREATE;
    pncp->ncp        = ncp;
    pncp->format     = format;

    if (safe_mode)          pncp->flag |= NC_MODE_SAFE;
    if (!relax_coord_bound) pncp->flag |= NC_MODE_STRICT_COORD_BOUND;

    return status;
}

#define _NDIMS_ 16

/*----< ncmpi_open() >-------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_open(MPI_Comm    comm,
           const char *path,
           int         omode,
           MPI_Info    info,
           int        *ncidp)  /* OUT */
{
    int i, j, nalloc, rank, nprocs, format, status=NC_NOERR, err;
    int safe_mode=0, mpireturn, relax_coord_bound, DIMIDS[_NDIMS_], *dimids;
    char *env_str;
    MPI_Info combined_info;
    void *ncp;
    PNC *pncp;
    PNC_driver *driver;
#ifdef BUILD_DRIVER_FOO
    int enable_foo_driver=0;
#endif
#ifdef ENABLE_BURST_BUFFER
    int enable_bb_driver=0;
#endif

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

#ifdef PNETCDF_DEBUG
    safe_mode = 1;
    /* When debug mode is enabled at the configure time, safe_mode is by
     * default enabled. This can be overwritten by the run-time environment
     * variable PNETCDF_SAFE_MODE */
#endif
    /* get environment variable PNETCDF_SAFE_MODE
     * if it is set to 1, then we perform a strict parameter consistent test
     */
    if ((env_str = getenv("PNETCDF_SAFE_MODE")) != NULL) {
        if (*env_str == '0') safe_mode = 0;
        else                 safe_mode = 1;
        /* if PNETCDF_SAFE_MODE is set but without a value, *env_str can
         * be '\0' (null character). In this case, safe_mode is enabled */
    }

    /* get environment variable PNETCDF_RELAX_COORD_BOUND
     * if it is set to 0, then we perform a strict start bound check
     */
#ifndef RELAX_COORD_BOUND
    relax_coord_bound = 0;
#else
    relax_coord_bound = 1;
#endif
    if ((env_str = getenv("PNETCDF_RELAX_COORD_BOUND")) != NULL) {
        if (*env_str == '0') relax_coord_bound = 0;
        else                 relax_coord_bound = 1;
        /* if PNETCDF_RELAX_COORD_BOUND is set but without a value, *env_str
         * can be '\0' (null character). This is equivalent to setting
         * relax_coord_bound to 1 */
    }

    /* path's validity is checked in MPI-IO with error code MPI_ERR_BAD_FILE
     * path consistency is checked in MPI-IO with error code MPI_ERR_NOT_SAME
     */
    if (path == NULL || *path == '\0') DEBUG_RETURN_ERROR(NC_EBAD_FILE)

    /* Check the file signature to tell the file format which is later used to
     * select the right driver.
     */
    format = NC_FORMAT_UNKNOWN;
    if (rank == 0) {
        err = ncmpi_inq_file_format(path, &format);
        if (err != NC_NOERR) {
            if (nprocs == 1) return err;
            format = err;
        }
        else if (format == NC_FORMAT_UNKNOWN) {
            if (nprocs == 1) DEBUG_RETURN_ERROR(NC_ENOTNC)
            format = NC_ENOTNC;
        }
#ifndef ENABLE_NETCDF4
        else if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC) {
            if (nprocs == 1) DEBUG_RETURN_ERROR(NC_ENOTBUILT)
            format = NC_ENOTBUILT;
        }
#endif
    }

    if (nprocs > 1) { /* root broadcasts format and omode */
        int root_omode, msg[2];

        msg[0] = format; /* only root's matters (format or error code) */

        /* Check omode consistency:
         * Note if omode contains NC_NOWRITE, it is equivalent to NC_CLOBBER.
         * In pnetcdf.h, they both are defined the same value, 0.
         * Only root's omode matters.
         */
        msg[1] = omode; /* only root's matters */

        TRACE_COMM(MPI_Bcast)(&msg, 2, MPI_INT, 0, comm);
        NCMPII_HANDLE_ERROR("MPI_Bcast")

        /* check format error (a fatal error, must return now) */
        format = msg[0];
        if (format < 0) return format; /* all netCDF errors are negative */

        /* check omode consistency */
        root_omode = msg[1];
        if (root_omode != omode) {
            omode = root_omode;
            DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE_OMODE)
        }

        if (safe_mode) { /* sync status among all processes */
            err = status;
            TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, comm);
            NCMPII_HANDLE_ERROR("MPI_Allreduce")
        }
        /* continue to use root's omode to open the file, but will report omode
         * inconsistency error, if there is any */
    }

    /* combine user's MPI info and PNETCDF_HINTS env variable */
    combine_env_hints(info, &combined_info);

#ifdef BUILD_DRIVER_FOO
    if (combined_info == MPI_INFO_NULL)
        MPI_Info_create(&combined_info);
    {
        char value[MPI_MAX_INFO_VAL];
        int flag;

        /* check if nc_foo_driver is enabled */
        MPI_Info_get(combined_info, "nc_foo_driver", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag && strcasecmp(value, "enable") == 0)
            enable_foo_driver = 1;

    }
#endif
#ifdef ENABLE_BURST_BUFFER
    if (combined_info == MPI_INFO_NULL)
        MPI_Info_create(&combined_info);
    {
        char value[MPI_MAX_INFO_VAL];
        int flag;

        /* check if nc_burst_buf is enabled */
        MPI_Info_get(combined_info, "nc_burst_buf", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag && strcasecmp(value, "enable") == 0)
            enable_bb_driver = 1;
    }
#endif

#ifdef ENABLE_NETCDF4
    if (format == NC_FORMAT_NETCDF4_CLASSIC || format == NC_FORMAT_NETCDF4) {
        driver = nc4io_inq_driver();
#ifdef ENABLE_BURST_BUFFER
        /* Burst buffering does not support NetCDF-4 files yet.
         * If hint nc_burst_buf is enabled in combined_info, disable it.
         */
        if (enable_bb_driver == 1)
            MPI_Info_set(combined_info, "nc_burst_buf", "disable");
        enable_bb_driver = 0;
#endif
    }
    else
#else
    if (format == NC_FORMAT_NETCDF4_CLASSIC || format == NC_FORMAT_NETCDF4)
        DEBUG_RETURN_ERROR(NC_ENOTBUILT)
    else
#endif
#ifdef BUILD_DRIVER_FOO
    if (enable_foo_driver)
        driver = ncfoo_inq_driver();
    else
#endif
#ifdef ENABLE_BURST_BUFFER
    if (enable_bb_driver)
        driver = ncbbio_inq_driver();
    else
#endif
    {
        /* ncmpio driver */
        if (format == NC_FORMAT_CLASSIC ||
            format == NC_FORMAT_CDF2 ||
            format == NC_FORMAT_CDF5) {
            driver = ncmpio_inq_driver();
        }
#ifdef ENABLE_ADIOS
        else if (format == NC_FORMAT_BP) {
            driver = ncadios_inq_driver();
        }
#endif
        else /* unrecognized file format */
            DEBUG_RETURN_ERROR(NC_ENOTNC)
    }

    /* allocate a PNC object */
    pncp = (PNC*) NCI_Malloc(sizeof(PNC));
    if (pncp == NULL) {
        *ncidp = -1;
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    /* generate a new nc file ID from NCPList */
    err = new_id_PNCList(ncidp, pncp);
    if (err != NC_NOERR) return err;

    /* Duplicate comm, because users may free it (though unlikely). Note
     * MPI_Comm_dup() is collective. We pass pncp->comm to drivers, so there
     * is no need for a driver to duplicate it again.
     */
    if (comm != MPI_COMM_WORLD && comm != MPI_COMM_SELF)
        MPI_Comm_dup(comm, &pncp->comm);
    else
        pncp->comm = comm;

    /* calling the driver's open subroutine */
    err = driver->open(pncp->comm, path, omode, *ncidp, combined_info, &ncp);
    if (status == NC_NOERR) status = err;
    if (combined_info != MPI_INFO_NULL) MPI_Info_free(&combined_info);
    if (status != NC_NOERR && status != NC_EMULTIDEFINE_OMODE &&
        status != NC_ENULLPAD) {
        /* NC_EMULTIDEFINE_OMODE and NC_ENULLPAD are not fatal error. We
         * continue the rest open procedure */
        del_from_PNCList(*ncidp);
        if (pncp->comm != MPI_COMM_WORLD && pncp->comm != MPI_COMM_SELF)
            MPI_Comm_free(&pncp->comm); /* a collective call */
        NCI_Free(pncp);
        *ncidp = -1;
        return status;
    }

    /* fill in pncp members */
    pncp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (pncp->path == NULL) {
        driver->close(ncp); /* close file and ignore error */
        del_from_PNCList(*ncidp);
        if (pncp->comm != MPI_COMM_WORLD && pncp->comm != MPI_COMM_SELF)
            MPI_Comm_free(&pncp->comm); /* a collective call */
        NCI_Free(pncp);
        *ncidp = -1;
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(pncp->path, path);
    pncp->mode       = omode;
    pncp->driver     = driver;
    pncp->ndims      = 0;
    pncp->unlimdimid = -1;
    pncp->nvars      = 0;
    pncp->vars       = NULL;
    pncp->flag       = 0;
    pncp->ncp        = ncp;
    pncp->format     = format;

    if (!fIsSet(omode, NC_WRITE)) pncp->flag |= NC_MODE_RDONLY;
    if (safe_mode)                pncp->flag |= NC_MODE_SAFE;
    if (!relax_coord_bound)       pncp->flag |= NC_MODE_STRICT_COORD_BOUND;

    /* inquire number of dimensions, variables defined and rec dim ID */
    err = driver->inq(pncp->ncp, &pncp->ndims, &pncp->nvars, NULL,
                      &pncp->unlimdimid);
    if (err != NC_NOERR) goto fn_exit;

    if (pncp->nvars == 0) return status; /* no variable defined in the file */

    /* make a copy of variable metadata at the dispatcher layer, because
     * sanity check is done at the dispatcher layer
     */

    /* allocate chunk size for pncp->vars[] */
    nalloc = _RNDUP(pncp->nvars, PNC_VARS_CHUNK);
    pncp->vars = NCI_Malloc(nalloc * sizeof(PNC_var));
    if (pncp->vars == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto fn_exit;
    }

    dimids = DIMIDS;

    /* construct array of PNC_var for all variables */
    for (i=0; i<pncp->nvars; i++) {
        int ndims, max_ndims=_NDIMS_;
        pncp->vars[i].shape  = NULL;
        pncp->vars[i].recdim = -1;   /* if fixed-size variable */
        err = driver->inq_var(pncp->ncp, i, NULL, &pncp->vars[i].xtype, &ndims,
                              NULL, NULL, NULL, NULL, NULL);
        if (err != NC_NOERR) break; /* loop i */
        pncp->vars[i].ndims = ndims;

        if (ndims > 0) {
            pncp->vars[i].shape = (MPI_Offset*)
                                  NCI_Malloc(ndims * SIZEOF_MPI_OFFSET);
            if (ndims > max_ndims) { /* avoid repeated malloc */
                if (dimids == DIMIDS) dimids = NULL;
                dimids = (int*) NCI_Realloc(dimids, ndims * SIZEOF_INT);
                max_ndims = ndims;
            }
            err = driver->inq_var(pncp->ncp, i, NULL, NULL, NULL,
                                  dimids, NULL, NULL, NULL, NULL);
            if (err != NC_NOERR) break; /* loop i */
            if (dimids[0] == pncp->unlimdimid)
                pncp->vars[i].recdim = pncp->unlimdimid;
            for (j=0; j<ndims; j++) {
                /* obtain size of dimension j */
                err = driver->inq_dim(pncp->ncp, dimids[j], NULL,
                                      pncp->vars[i].shape+j);
                if (err != NC_NOERR) break; /* loop i */
            }
        }
    }
    if (err != NC_NOERR) { /* error happens in loop i */
        assert(i < pncp->nvars);
        for (j=0; j<=i; j++) {
            if (pncp->vars[j].shape != NULL)
                NCI_Free(pncp->vars[j].shape);
        }
        NCI_Free(pncp->vars);
    }
    if (dimids != DIMIDS) NCI_Free(dimids);

fn_exit:
    if (err != NC_NOERR) {
        driver->close(ncp); /* close file and ignore error */
        if (pncp->comm != MPI_COMM_WORLD && pncp->comm != MPI_COMM_SELF)
            MPI_Comm_free(&pncp->comm); /* a collective call */
        del_from_PNCList(*ncidp);
        NCI_Free(pncp->path);
        NCI_Free(pncp);
        *ncidp = -1;
        if (status == NC_NOERR) status = err;
    }

    return status;
}

/*----< ncmpi_close() >------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_close(int ncid)
{
    int i, err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_close() */
    err = pncp->driver->close(pncp->ncp);

    /* Remove from the PNCList, even if err != NC_NOERR */
    del_from_PNCList(ncid);

    /* free the PNC object */
    if (pncp->comm != MPI_COMM_WORLD && pncp->comm != MPI_COMM_SELF)
        MPI_Comm_free(&pncp->comm); /* a collective call */

    NCI_Free(pncp->path);
    for (i=0; i<pncp->nvars; i++)
        if (pncp->vars[i].shape != NULL)
            NCI_Free(pncp->vars[i].shape);
    if (pncp->vars != NULL)
        NCI_Free(pncp->vars);
    NCI_Free(pncp);

    return err;
}

/*----< ncmpi_enddef() >-----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_enddef(int ncid) {
    int err=NC_NOERR;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (!(pncp->flag & NC_MODE_DEF)) DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)

    if (pncp->flag & NC_MODE_SAFE) { /* safe mode */
        int minE, mpireturn;
        /* check the error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }
    else if (err != NC_NOERR) return err; /* fatal error */

    /* calling the subroutine that implements ncmpi_enddef() */
    err = pncp->driver->enddef(pncp->ncp);
    if (err != NC_NOERR) return err;

    fClr(pncp->flag, NC_MODE_INDEP); /* default enters collective data mode */
    fClr(pncp->flag, NC_MODE_DEF);
    return NC_NOERR;
}

/*----< ncmpi__enddef() >----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi__enddef(int        ncid,
              MPI_Offset h_minfree,
              MPI_Offset v_align,
              MPI_Offset v_minfree,
              MPI_Offset r_align)
{
    int err=NC_NOERR;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (!(pncp->flag & NC_MODE_DEF)) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

    if (h_minfree < 0 || v_align < 0 || v_minfree < 0 || r_align < 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto err_check;
    }

err_check:
    if (pncp->flag & NC_MODE_SAFE) { /* safe mode */
        int minE, mpireturn;
        MPI_Offset root_args[4];

        /* first check the error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;

        /* check if h_minfree, v_align, v_minfree, and r_align are consistent
         * among all processes */
        root_args[0] = h_minfree;
        root_args[1] = v_align;
        root_args[2] = v_minfree;
        root_args[3] = r_align;
        TRACE_COMM(MPI_Bcast)(&root_args, 4, MPI_OFFSET, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");

        if (root_args[0] != h_minfree ||
            root_args[1] != v_align   ||
            root_args[2] != v_minfree ||
            root_args[3] != r_align)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }
    else if (err != NC_NOERR) return err; /* fatal error */

    /* calling the subroutine that implements ncmpi__enddef() */
    err = pncp->driver->_enddef(pncp->ncp, h_minfree, v_align,
                                           v_minfree, r_align);
    if (err != NC_NOERR) return err;

    fClr(pncp->flag, NC_MODE_INDEP); /* default enters collective data mode */
    fClr(pncp->flag, NC_MODE_DEF);
    return NC_NOERR;
}

/*----< ncmpi_redef() >------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_redef(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (fIsSet(pncp->flag, NC_MODE_RDONLY)) /* read-only */
        DEBUG_RETURN_ERROR(NC_EPERM)
    /* if open mode is inconsistent, then this return might cause parallel
     * program to hang */

    /* cannot be in define mode, must enter from data mode */
    if (fIsSet(pncp->flag, NC_MODE_DEF)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* calling the subroutine that implements ncmpi_redef() */
    err = pncp->driver->redef(pncp->ncp);
    if (err != NC_NOERR) return err;

    fSet(pncp->flag, NC_MODE_DEF);
    return NC_NOERR;
}

/*----< ncmpi_sync() >-------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpi_sync(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_sync() */
    return pncp->driver->sync(pncp->ncp);
}

/*----< ncmpi_flush() >-------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpi_flush(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_flush() */
    return pncp->driver->flush(pncp->ncp);
}

/*----< ncmpi_abort() >------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_abort(int ncid)
{
    int i, err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_abort() */
    err = pncp->driver->abort(pncp->ncp);

    /* Remove from the PNCList, even if err != NC_NOERR */
    del_from_PNCList(ncid);

    /* free the PNC object */
    if (pncp->comm != MPI_COMM_WORLD && pncp->comm != MPI_COMM_SELF)
        MPI_Comm_free(&pncp->comm); /* a collective call */

    NCI_Free(pncp->path);
    for (i=0; i<pncp->nvars; i++)
        if (pncp->vars[i].shape != NULL)
            NCI_Free(pncp->vars[i].shape);
    if (pncp->vars != NULL)
        NCI_Free(pncp->vars);
    NCI_Free(pncp);

    return err;
}

/*----< ncmpi_set_fill() >---------------------------------------------------*/
/* This is a collective subroutine.
 * This subroutine serves both purposes of setting and inquiring the fill mode.
 */
int
ncmpi_set_fill(int  ncid,
               int  fill_mode,     /* mode to be changed by user */
               int *old_fill_mode) /* current fill mode */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (fIsSet(pncp->flag, NC_MODE_RDONLY)) /* read-only */
        DEBUG_RETURN_ERROR(NC_EPERM)

    /* not allowed to call in data mode for classic formats */
    if ((pncp->format != NC_FORMAT_NETCDF4) && !(pncp->flag & NC_MODE_DEF))
        DEBUG_RETURN_ERROR(NC_ENOTINDEFINE)

    /* calling the subroutine that implements ncmpi_set_fill() */
    err = pncp->driver->set_fill(pncp->ncp, fill_mode, old_fill_mode);
    if (err != NC_NOERR) return err;

    if (fill_mode == NC_FILL)
        fSet(pncp->flag, NC_MODE_FILL);
    else /* NC_NOFILL */
        fClr(pncp->flag, NC_MODE_FILL);

    return NC_NOERR;
}

/*----< ncmpi_inq_format() >-------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_format(int  ncid,
                 int *formatp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (formatp != NULL) *formatp = pncp->format;

    return NC_NOERR;
}

#ifdef ENABLE_ADIOS
static void swap_64(void *data)
{
    uint64_t *dest = (uint64_t*) data;
    uint64_t tmp;
    memcpy(&tmp, dest, 8);
    *dest = ((tmp & 0x00000000000000FFULL) << 56) |
            ((tmp & 0x000000000000FF00ULL) << 40) |
            ((tmp & 0x0000000000FF0000ULL) << 24) |
            ((tmp & 0x00000000FF000000ULL) <<  8) |
            ((tmp & 0x000000FF00000000ULL) >>  8) |
            ((tmp & 0x0000FF0000000000ULL) >> 24) |
            ((tmp & 0x00FF000000000000ULL) >> 40) |
            ((tmp & 0xFF00000000000000ULL) >> 56);
}

static int adios_parse_endian(char *footer, int *diff_endianness) {
    unsigned int version;
    unsigned int test = 1; /* If high bit set, big endian */

    version = ntohl (*(uint32_t *) (footer + BP_MINIFOOTER_SIZE - 4));
    char *v = (char *) (&version);
    if ((*v && !*(char *) &test) /* Both writer and reader are big endian */
        || (!*(v+3) && *(char *) &test)){ /* Both are little endian */
        *diff_endianness = 0; /* No need to change endianness */
    }
    else{
        *diff_endianness = 1;
    }

    return 0;
}
#endif

/*----< ncmpi_inq_file_format() >--------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_file_format(const char *filename,
                      int        *formatp) /* out */
{
    const char *cdf_signature="CDF";
    const char *hdf5_signature="\211HDF\r\n\032\n";
    const char *path;
    char signature[8];
    int fd;
    ssize_t rlen;

    if (formatp == NULL) return NC_NOERR;

    *formatp = NC_FORMAT_UNKNOWN;

    /* remove the file system type prefix name if there is any.
     * For example, when filename = "lustre:/home/foo/testfile.nc", remove
     * "lustre:" to make path = "/home/foo/testfile.nc" in open() below
     */
    path = strchr(filename, ':');
    if (path == NULL) path = filename; /* no prefix */
    else              path++;

    /* must include config.h on 32-bit machines, as AC_SYS_LARGEFILE is called
     * at the configure time and it defines _FILE_OFFSET_BITS to 64 if large
     * file feature is supported.
     */
    if ((fd = open(path, O_RDONLY, 00400)) == -1) { /* open for read */
             if (errno == ENOENT)       DEBUG_RETURN_ERROR(NC_ENOENT)
        else if (errno == EACCES)       DEBUG_RETURN_ERROR(NC_EACCESS)
        else if (errno == ENAMETOOLONG) DEBUG_RETURN_ERROR(NC_EBAD_FILE)
        else {
            fprintf(stderr,"Error on opening file %s (%s)\n",
                    filename,strerror(errno));
            DEBUG_RETURN_ERROR(NC_EFILE)
        }
    }
    /* get first 8 bytes of file */
    rlen = read(fd, signature, 8);
    if (rlen != 8) {
        close(fd); /* ignore error */
        DEBUG_RETURN_ERROR(NC_EFILE)
    }
    if (close(fd) == -1) {
        DEBUG_RETURN_ERROR(NC_EFILE)
    }

    if (memcmp(signature, cdf_signature, 3) == 0) {
             if (signature[3] == 5)  *formatp = NC_FORMAT_CDF5;
        else if (signature[3] == 2)  *formatp = NC_FORMAT_CDF2;
        else if (signature[3] == 1)  *formatp = NC_FORMAT_CLASSIC;
    }

    /* check if the file is an HDF5. */
    if (*formatp == NC_FORMAT_UNKNOWN) {
        /* The HDF5 superblock is located by searching for the HDF5 format
         * signature at byte offset 0, byte offset 512, and at successive
         * locations in the file, each a multiple of two of the previous
         * location; in other words, at these byte offsets: 0, 512, 1024, 2048,
         * and so on. The space before the HDF5 superblock is referred as to
         * "user block".
         */
        off_t offset=0;

        fd = open(path, O_RDONLY, 00400); /* error check already done */
        /* get first 8 bytes of file */
        rlen = read(fd, signature, 8); /* error check already done */

        while (rlen == 8 && memcmp(signature, hdf5_signature, 8)) {
            offset = (offset == 0) ? 512 : offset * 2;
            lseek(fd, offset, SEEK_SET);
            rlen = read(fd, signature, 8);
        }
        close(fd); /* ignore error */

        if (rlen == 8) { /* HDF5 signature found */
            /* TODO: whether the file is NC_FORMAT_NETCDF4_CLASSIC is
             * determined by HDF5 attribute "_nc3_strict" which requires a call
             * to H5Aget_name(). For now, we do not distinguish
             * NC_CLASSIC_MODEL, but simply return NETCDF4 format.
             */
#ifdef ENABLE_NETCDF4
            int err, ncid;
            err = nc_open(path, NC_NOWRITE, &ncid);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
            err = nc_inq_format(ncid, formatp);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
            err = nc_close(ncid);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
#else
            *formatp = NC_FORMAT_NETCDF4;
#endif
        }
    }

#ifdef ENABLE_ADIOS
    /* check if the file is a BP. */
    if (*formatp == NC_FORMAT_UNKNOWN) {
        off_t fsize;
        int diff_endian;
        char footer[BP_MINIFOOTER_SIZE];
        off_t h1, h2, h3;

        /* test if the file footer follows BP specification */
        if ((fd = open(path, O_RDONLY, 00400)) == -1) {
                 if (errno == ENOENT)       DEBUG_RETURN_ERROR(NC_ENOENT)
            else if (errno == EACCES)       DEBUG_RETURN_ERROR(NC_EACCESS)
            else if (errno == ENAMETOOLONG) DEBUG_RETURN_ERROR(NC_EBAD_FILE)
            else {
                fprintf(stderr,"Error on opening file %s (%s)\n",
                        filename,strerror(errno));
                DEBUG_RETURN_ERROR(NC_EFILE)
            }
        }

        /* Seek to end of file */
        fsize = lseek(fd, (off_t)(-(BP_MINIFOOTER_SIZE)), SEEK_END);

        /* read footer */
        rlen = read(fd, footer, BP_MINIFOOTER_SIZE);
        if (rlen != BP_MINIFOOTER_SIZE) {
            close(fd);
            DEBUG_RETURN_ERROR(NC_EFILE)
        }
        if (close(fd) == -1) {
            DEBUG_RETURN_ERROR(NC_EFILE)
        }

        /* check endianness of file and this running system */
        adios_parse_endian(footer, &diff_endian);

        BUFREAD64(footer,      h1) /* file offset of process group index table */
        BUFREAD64(footer + 8,  h2) /* file offset of variable index table */
        BUFREAD64(footer + 16, h3) /* file offset of attribute index table */

        /* All index tables must fall within the range of file size.
         * Process group index table must comes before variable index table.
         * Variable index table must comes before attribute index table.
         */
        if (0 < h1 && h1 < fsize &&
            0 < h2 && h2 < fsize &&
            0 < h3 && h3 < fsize &&
            h1 < h2 && h2 < h3){
            /* basic footer check is passed, now we try to open the file with
             * ADIOS library to make sure it is indeed a BP formated file
             */
            ADIOS_FILE *fp;
            fp = adios_read_open_file(path, ADIOS_READ_METHOD_BP,
                                        MPI_COMM_SELF);
            if (fp != NULL) {
                *formatp = NC_FORMAT_BP;
                adios_read_close(fp);
            }
        }
    }
#endif

    return NC_NOERR;
}

/*----< ncmpi_inq_version() >------------------------------------------------*/
int
ncmpi_inq_version(int ncid, int *nc_mode)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (nc_mode == NULL) return NC_NOERR;

    if (pncp->format == NC_FORMAT_CDF5)
        *nc_mode = NC_64BIT_DATA;
    else if (pncp->format == NC_FORMAT_CDF2)
        *nc_mode = NC_64BIT_OFFSET;
    else if (pncp->format == NC_FORMAT_CLASSIC)
        *nc_mode = NC_CLASSIC_MODEL;

#ifdef ENABLE_NETCDF4
    else if (pncp->format == NC_FORMAT_NETCDF4)
        *nc_mode = NC_NETCDF4;
    else if (pncp->format == NC_FORMAT_NETCDF4_CLASSIC)
        *nc_mode = NC_NETCDF4 | NC_CLASSIC_MODEL;
#endif

#ifdef ENABLE_ADIOS
    else if (pncp->format == NC_FORMAT_BP)
        *nc_mode = NC_BP;
#endif

    return NC_NOERR;
}

/*----< ncmpi_inq() >--------------------------------------------------------*/
int
ncmpi_inq(int  ncid,
          int *ndimsp,
          int *nvarsp,
          int *nattsp,
          int *xtendimp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq() */
    return pncp->driver->inq(pncp->ncp, ndimsp, nvarsp, nattsp, xtendimp);
}

/*----< ncmpi_inq_ndims() >--------------------------------------------------*/
int
ncmpi_inq_ndims(int  ncid,
                int *ndimsp)
{
    return ncmpi_inq(ncid, ndimsp, NULL, NULL, NULL);
}

/*----< ncmpi_inq_nvars() >--------------------------------------------------*/
int
ncmpi_inq_nvars(int  ncid,
                int *nvarsp)
{
    return ncmpi_inq(ncid, NULL, nvarsp, NULL, NULL);
}

/*----< ncmpi_inq_natts() >--------------------------------------------------*/
int
ncmpi_inq_natts(int  ncid,
                int *nattsp)
{
    return ncmpi_inq(ncid, NULL, NULL, nattsp, NULL);
}

/*----< ncmpi_inq_unlimdim() >-----------------------------------------------*/
int
ncmpi_inq_unlimdim(int  ncid,
                   int *unlimdimidp)
{
    return ncmpi_inq(ncid, NULL, NULL, NULL, unlimdimidp);
}

/*----< ncmpi_inq_path() >---------------------------------------------------*/
/* Get the file pathname which was used to open/create the ncid's file.
 * pathlen and path must already be allocated. Ignored if NULL.
 * This is an independent subroutine.
 */
int
ncmpi_inq_path(int   ncid,
               int  *pathlen,/* Ignored if NULL */
               char *path)   /* must have already been allocated. Ignored if NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

#if 0
    /* calling the subroutine that implements ncmpi_inq_path() */
    return pncp->driver->inq_misc(pncp->ncp, pathlen, path, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
#endif
    if (pathlen != NULL) {
        if (pncp->path == NULL) *pathlen = 0;
        else                    *pathlen = (int)strlen(pncp->path);
    }
    if (path != NULL) {
        if (pncp->path == NULL) *path = '\0';
        else                    strcpy(path, pncp->path);
    }
    return NC_NOERR;
}

/*----< ncmpi_inq_num_fix_vars() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_num_fix_vars(int ncid, int *num_fix_varsp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (num_fix_varsp == NULL) return NC_NOERR;

#ifdef ENABLE_NETCDF4
    if (pncp->format == NC_FORMAT_NETCDF4 ||
        pncp->format == NC_FORMAT_NETCDF4_CLASSIC) {
        /* calling the subroutine that implements ncmpi_inq_num_fix_vars() */
        return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, num_fix_varsp,
                                      NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                      NULL, NULL, NULL, NULL, NULL);
    }
#endif

    *num_fix_varsp = pncp->nvars - pncp->nrec_vars;

    /* number of fixed-size variables can also be calculated below.
    int i;
    *num_fix_varsp = 0;
    for (i=0; i<pncp->nvars; i++) {
        if (pncp->vars[i].recdim < 0)
            (*num_fix_varsp)++;
    }
    */

    return NC_NOERR;
}

/*----< ncmpi_inq_num_rec_vars() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_num_rec_vars(int ncid, int *num_rec_varsp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (num_rec_varsp == NULL) return NC_NOERR;

#ifdef ENABLE_NETCDF4
    if (pncp->format == NC_FORMAT_NETCDF4 ||
        pncp->format == NC_FORMAT_NETCDF4_CLASSIC) {
        /* calling the subroutine that implements ncmpi_inq_num_rec_vars() */
        return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL,
                                      num_rec_varsp, NULL, NULL, NULL, NULL,
                                      NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        }
#endif

    *num_rec_varsp = pncp->nrec_vars;

    /* number of record variables can also be calculated below.
    int i;
    *num_rec_varsp = 0;
    for (i=0; i<pncp->nvars; i++) {
        if (pncp->vars[i].recdim >= 0)
            (*num_rec_varsp)++;
    }
    */

    return NC_NOERR;
}

/*----< ncmpi_inq_striping() >-----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_striping(int ncid, int *striping_size, int *striping_count)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_striping() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  striping_size, striping_count, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_header_size() >--------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_header_size(int ncid, MPI_Offset *header_size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (header_size == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_header_size() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, header_size, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_header_extent() >------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_header_extent(int ncid, MPI_Offset *header_extent)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (header_extent == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_header_extent() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, header_extent, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_recsize() >------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_recsize(int ncid, MPI_Offset *recsize)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (recsize == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_recsize() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, recsize, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_put_size() >-----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_put_size(int ncid, MPI_Offset *put_size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (put_size == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_put_size() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, put_size,
                                  NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_get_size() >-----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_get_size(int ncid, MPI_Offset *get_size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (get_size == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_get_size() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  get_size, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_file_info() >----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_file_info(int ncid, MPI_Info *info)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (info == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_file_info() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, info, NULL, NULL, NULL);
}

/* ncmpi_get_file_info() is now deprecated, replaced by ncmpi_inq_file_info() */
int
ncmpi_get_file_info(int ncid, MPI_Info *info)
{
    return ncmpi_inq_file_info(ncid, info);
}

/*----< ncmpi_begin_indep_data() >-------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_begin_indep_data(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_begin_indep_data() */
    err = pncp->driver->begin_indep_data(pncp->ncp);
    if (err != NC_NOERR) return err;

    fSet(pncp->flag, NC_MODE_INDEP);
    return NC_NOERR;
}

/*----< ncmpi_end_indep_data() >---------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_end_indep_data(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_end_indep_data() */
    err = pncp->driver->end_indep_data(pncp->ncp);
    if (err != NC_NOERR) return err;

    fClr(pncp->flag, NC_MODE_INDEP);
    return NC_NOERR;
}

/*----< ncmpi_sync_numrecs() >-----------------------------------------------*/
/* this API is collective, but can be called in independent data mode.
 * Note numrecs (number of records) is always sync-ed in memory and file in
 * collective data mode.
 */
int
ncmpi_sync_numrecs(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_sync_numrecs() */
    return pncp->driver->sync_numrecs(pncp->ncp);
}

/*----< ncmpi_set_default_format() >-----------------------------------------*/
/* This function sets a default create file format.
 * Valid formats are NC_FORMAT_CLASSIC, NC_FORMAT_CDF2, and NC_FORMAT_CDF5
 * This API is NOT collective, as there is no way to check against an MPI
 * communicator. It should be called by all MPI processes that intend to
 * create a file later. Consistency check will have to be done in other APIs.
 */
int
ncmpi_set_default_format(int format, int *old_formatp)
{
    int err;

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif

    /* Return existing format if desired. */
    if (old_formatp != NULL)
        *old_formatp = ncmpi_default_create_format;

    /* Make sure only valid format is set. */
    if (format != NC_FORMAT_CLASSIC &&
        format != NC_FORMAT_CDF2 &&
        format != NC_FORMAT_NETCDF4 &&
        format != NC_FORMAT_NETCDF4_CLASSIC &&
        format != NC_FORMAT_CDF5) {
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
    }
    else {
        ncmpi_default_create_format = format;
        err = NC_NOERR;
    }

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif

    return err;
}

/*----< ncmpi_inq_default_format() >-----------------------------------------*/
/* returns a value suitable for a create flag.  Will return one or more of the
 * following values OR-ed together:
 * NC_64BIT_OFFSET, NC_CLOBBER, NC_LOCK, NC_SHARE */
int
ncmpi_inq_default_format(int *formatp)
{
    if (formatp == NULL) DEBUG_RETURN_ERROR(NC_EINVAL)

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif

    *formatp = ncmpi_default_create_format;

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif

    return NC_NOERR;
}

/*----< ncmpi_inq_files_opened() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_files_opened(int *num,    /* cannot be NULL */
                       int *ncids)  /* can be NULL */
{
    int i;

    if (num == NULL) DEBUG_RETURN_ERROR(NC_EINVAL)

#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif

    *num = pnc_numfiles;

    if (ncids != NULL) { /* ncids can be NULL */
        *num = 0;
        for (i=0; i<NC_MAX_NFILES; i++) {
            if (pnc_filelist[i] != NULL) {
                ncids[*num] = i;
                (*num)++;
            }
        }
    }
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif

    return NC_NOERR;
}

/*----< ncmpi_inq_nreqs() >--------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_nreqs(int  ncid,
                int *nreqs) /* number of pending nonblocking requests */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (nreqs == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_nreqs() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, NULL, nreqs, NULL, NULL);
}

/*----< ncmpi_inq_buffer_usage() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_buffer_usage(int         ncid,
                       MPI_Offset *usage) /* amount of space used so far */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (usage == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_buffer_usage() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, usage, NULL);
}

/*----< ncmpi_inq_buffer_size() >--------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_buffer_size(int         ncid,
                      MPI_Offset *buf_size) /* amount of space attached */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (buf_size == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_buffer_size() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, buf_size);
}

/*----< ncmpi_buffer_attach() >----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_buffer_attach(int        ncid,
                    MPI_Offset bufsize) /* amount of memory space allowed for
                                           PnetCDF library to buffer the
                                           nonblocking requests */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_buffer_attach() */
    return pncp->driver->buffer_attach(pncp->ncp, bufsize);
}

/*----< ncmpi_buffer_detach() >----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_buffer_detach(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_buffer_detach() */
    return pncp->driver->buffer_detach(pncp->ncp);
}

/*----< ncmpi_delete() >-----------------------------------------------------*/
/*
 * filename: the name of the file we will remove.
 * info: MPI info object, in case underlying file system needs hints.
 *
 * This API is implemented in src/driver/ncmpio/ncmpio_file.c
 *
 */

/*----< ncmpi_wait() >-------------------------------------------------------*/
/* This API is an independent subroutine. */
int
ncmpi_wait(int  ncid,
           int  num_reqs, /* number of requests */
           int *req_ids,  /* [num_reqs]: IN/OUT */
           int *statuses) /* [num_reqs], can be NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_wait() */
    return pncp->driver->wait(pncp->ncp, num_reqs, req_ids, statuses,
                              NC_REQ_INDEP);
}

/*----< ncmpi_wait_all() >---------------------------------------------------*/
/* This API is a collective subroutine. */
int
ncmpi_wait_all(int  ncid,
               int  num_reqs, /* number of requests */
               int *req_ids,  /* [num_reqs]: IN/OUT */
               int *statuses) /* [num_reqs], can be NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_wait_all() */
    return pncp->driver->wait(pncp->ncp, num_reqs, req_ids, statuses,
                              NC_REQ_COLL);
}

/*----< ncmpi_cancel() >-----------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_cancel(int  ncid,
             int  num_reqs, /* number of requests */
             int *req_ids,  /* [num_reqs]: IN/OUT */
             int *statuses) /* [num_reqs], can be NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_cancel() */
    return pncp->driver->cancel(pncp->ncp, num_reqs, req_ids, statuses);
}

