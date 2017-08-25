/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>
#include <string.h>  /* strtok(), strcpy(), strchr() */
#include <strings.h> /* strcasecmp() */

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncx.h"
#ifdef ENABLE_SUBFILING
#include "subfile.h"
#endif

/* Begin Of Dataset Functions */

/*----< ncmpii_create() >----------------------------------------------------*/
/**  \ingroup datasets
Create a new netCDF file.

This function creates a new netCDF dataset, returning a netCDF ID that can
subsequently be used to refer to the netCDF dataset in other PnetCDF function
calls. The new netCDF dataset opened for write access and placed in define
mode, ready for you to add dimensions, variables, and attributes.

\param comm The MPI communicator. This API is a collective routine: all
processes must provide the same value for cmode, and all processes must provide
filenames that reference the same file. (Values for info may vary.) comm must
be an MPI intracommunicator.

\param path The file name of the new netCDF dataset.

\param cmode The creation mode flag. The following flags are available:
NC_NOCLOBBER, NC_SHARE, NC_64BIT_OFFSET, and NC_64BIT_DATA.

\param info MPI info object. It is used to provide file access hints,including
existing MPI hints as well as PnetCDF hints.  For MPI hints, users are referred
to MPI user guide for further information. For PnetCDF hints see below.

\param ncidp Pointer to location where returned netCDF ID is to be stored.

<h2>The cmode Flag</h2>

The cmode flag is used to control the type of file created, and some aspects of
how it may be used.

Setting NC_NOCLOBBER means you do not want to clobber (overwrite) an existing
dataset; an error (NC_EEXIST) is returned if the specified dataset already
exists.

The NC_SHARE flag in PnetCDF does not mean sharing the file with the processes
in this MPI program. Instead, it means the file will be concurrently shared
by a different MPI program. Hence, PnetCDF calls MPI_File_sync() right after
every time an MPI_File_write() call is made. This includes writing to metadata
(file header) as well as array data.

Setting NC_64BIT_OFFSET causes PnetCDF to create a 64-bit offset format file
(CDF-2), instead of a netCDF classic format file.  The 64-bit offset format
imposes far fewer restrictions on large (i.e. over 2 GB) data files.  See Large
File Support (The PnetCDF Users Guide).

Setting NC_64BIT_DATA causes PnetCDF to create a 64-bit data format file
(CDF-5). The 64-bit data format allows define variables with more than 4
billion array elements. See Large File Support (The PnetCDF Users Guide).

A zero value (defined for convenience as NC_CLOBBER) specifies the default
behavior: overwrite any existing dataset with the same file name.

<h2>The info object Flag</h2>

Starting from version 1.3.1, the following PnetCDF hints are available:

- nc_header_align_size: This hint allows some extra space between the end of
the header describing the entire file and the first variable. If you have an
application that periodically wishes to add more variables to an already
existing file, expanding the file header size may result in an expensive move
of the entire data file to make room for the definition of the new variables.
Hence, setting this hint to a value that is big enough to accommodate any
additional variables means you may leave your application code as-is and yet
still see tremendous performance improvements.

- nc_var_align_size: If you are writing to a block-based parallel file system,
such as IBM's GPFS or Lustre, then an application write becomes a block write
at the file system layer. If a write straddles two blocks, then locks must be
acquired for both blocks. Aligning the start of a variable to a block boundary,
combined with collective I/O optimization in the MPI-IO library can often
eliminate all unaligned file system accesses.

- nc_record_align_size: This hint aligns the starting file offset of the
record variable section.

- nc_header_read_chunk_size: PnetCDF reads the file headers in chunks. This
hint indicates the chunk size (in bytes). The default is 256 KB.

\returns ::NC_NOERR No error.
\returns ::NC_ENOMEM System out of memory.
\returns ::NC_EEXIST Specified file name exists when using NC_NOCLOBBER.
Can be use to check if the file exists.
\returns ::NC_EMULTIDEFINE_OMODE Bad file create/open mode or modes are
inconsistent across processes
\returns ::NC_EFILE: Unknown error in file operation

<h1>Examples</h1>

In this example we create a netCDF dataset named foo.nc; we want the dataset to
be created in the current directory only if a dataset with that name does not
already exist:

@code
     #include <mpi.h>
     #include <pnetcdf.h>
        ...
     int status;
     int ncid;
     MPI_Info info;
        ...
     MPI_Info_create (&info);
     MPI_Info_set (info, "romio_no_indep_rw",    "true");
     MPI_Info_set (info, "nc_header_align_size", "4194304");
     MPI_Info_set (info, "nc_var_align_size",    "1048576");
     MPI_Info_set (info, "nc_record_align_size", "1048576");

     status = ncmpi_create(MPI_COMM_WORLD, "foo.nc", NC_NOCLOBBER, info, &ncid);
     if (status != NC_NOERR) handle_error(status);
     MPI_Info_free(&info);
@endcode

In this example we create a netCDF dataset named foo.nc. It will
be in the CDF-5 format.

@code
     #include <mpi.h>
     #include <pnetcdf.h>
        ...
     int status;
     int ncid;
     int cmode = NC_NOCLOBBER | NC_64BIT_DATA;
     MPI_Info info = MPI_INFO_NULL;
        ...
     status = ncmpi_create(MPI_COMM_WORLD, "foo.nc", cmode, info, &ncid);
     if (status != NC_NOERR) handle_error(status);
@endcode
*/
int
ncmpii_create(MPI_Comm     comm,
              const char  *path,
              int          cmode,
              int          ncid,
              MPI_Info     info,
              void       **ncpp)
{
    int i, err, safe_mode=0, mpireturn, default_format;
    char *env_str;
    MPI_Info   env_info=MPI_INFO_NULL;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;
    NC *ncp=NULL;

#ifdef PNETCDF_DEBUG
    safe_mode = 1;
    /* this configure time setting will be overwritten by the run-time
     * environment variable PNETCDF_SAFE_MODE */
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

    /* path's validity and cmode consistency have been checked in ncmpi_create()
     * in src/dispatchers/file.c */

    /* check default format */
    ncmpi_inq_default_format(&default_format);

#if SIZEOF_MPI_OFFSET <  8
    /* check cmode */
    if (fIsSet(cmode, NC_64BIT_DATA)     ||
        fIsSet(cmode, NC_64BIT_OFFSET)   ||
        default_format == NC_FORMAT_CDF5 || 
        default_format == NC_FORMAT_CDF2) {
        /* unlike serial netcdf, we will not bother to support
         * NC_64BIT_OFFSET on systems with off_t smaller than 8 bytes.
         * serial netcdf has proven it's possible if datasets are small, but
         * that's a hassle we don't want to worry about */
        DEBUG_RETURN_ERROR(NC_ESMALL)
    }
#endif

    /* NC_DISKLESS is not supported yet */
    if (cmode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* NC_MMAP is not supported yet */
    if (cmode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* It is illegal to have both NC_64BIT_OFFSET & NC_64BIT_DATA */
    if ((cmode & (NC_64BIT_OFFSET|NC_64BIT_DATA)) ==
                 (NC_64BIT_OFFSET|NC_64BIT_DATA)) {
        DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)
    }

    /* take hints from the environment variable PNETCDF_HINTS
     * a string of hints separated by ";" and each hint is in the
     * form of hint=value. E.g. cb_nodes=16;cb_config_list=*:6
     * If this environment variable is set, it overrides any values that
     * were set by using calls to MPI_Info_set in the application code.
     */
    if (info != MPI_INFO_NULL) {
#ifdef HAVE_MPI_INFO_DUP
        mpireturn = MPI_Info_dup(info, &env_info);
        if (mpireturn != MPI_SUCCESS)
            DEBUG_RETURN_ERROR(ncmpii_handle_error(mpireturn, "MPI_Info_dup"))
#else
        printf("Warning: MPI info is ignored as MPI_Info_dup() is missing\n");
#endif
    }
    if ((env_str = getenv("PNETCDF_HINTS")) != NULL) {
        if (env_info == MPI_INFO_NULL)
            MPI_Info_create(&env_info); /* ignore error */

        char *env_str_cpy, *key;
        env_str_cpy = (char*) NCI_Malloc(strlen(env_str)+1);
        strcpy(env_str_cpy, env_str);
        key = strtok(env_str_cpy, ";");
        while (key != NULL) {
            char *val;
            val = strchr(key, '=');
            if (val == NULL) continue; /* ill-formed hint */
            *val = '\0';
            val++;
            /* printf("env hint: key=%s val=%s\n",key,val); */
            MPI_Info_set(env_info, key, val); /* override or add */
            key = strtok(NULL, ";");
        }
        NCI_Free(env_str_cpy);
    }

    /* get header chunk size from user info */
    if (env_info != MPI_INFO_NULL) {
        int flag;
        char value[MPI_MAX_INFO_VAL];
        MPI_Info_get(env_info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;
            chunksize = strtoll(value,NULL,10);
            if (errno != 0) chunksize = NC_DEFAULT_CHUNKSIZE;
        }
    }

    /* allocate buffer for header object NC */
    ncp = ncmpii_new_NC(&chunksize);
    if (ncp == NULL) {
        if (env_info != MPI_INFO_NULL) MPI_Info_free(&env_info);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    ncp->safe_mode = safe_mode;
    ncp->abuf      = NULL;
    ncp->old       = NULL;
#ifdef ENABLE_SUBFILING
    ncp->subfile_mode = 1;
    if (env_info != MPI_INFO_NULL) {
        int flag;
        char value[MPI_MAX_INFO_VAL];
        MPI_Info_get(env_info, "pnetcdf_subfiling", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag && strcasecmp(value, "disable") == 0)
            ncp->subfile_mode = 0;
    }
    ncp->ncid_sf = -1; /* subfile ncid; init to -1 */
    ncp->nc_num_subfiles = 0; /* num_subfiles; init to 0 */
#endif

    /* set the file format version based on the create mode, cmode */
    if (fIsSet(cmode, NC_64BIT_DATA)) {
        fSet(ncp->flags, NC_64BIT_DATA);
        ncp->format = 5;
    } else if (fIsSet(cmode, NC_64BIT_OFFSET)) {
        fSet(ncp->flags, NC_64BIT_OFFSET);
        ncp->format = 2;
    } else {
        if (default_format == NC_FORMAT_CDF5) {
            fSet(ncp->flags, NC_64BIT_DATA);
            ncp->format = 5;
        }
        else if (default_format == NC_FORMAT_CDF2) {
            fSet(ncp->flags, NC_64BIT_OFFSET);
            ncp->format = 2;
        }
        else {
            fSet(ncp->flags, NC_32BIT);
            ncp->format = 1;
        }
    }

    /* find the true header size (not-yet aligned) */
    ncp->xsz = ncmpii_hdr_len_NC(ncp);

    /* PnetCDF default fill mode is no fill */
    fSet(ncp->flags, NC_NOFILL);

    /* open file collectively */
    err = ncmpiio_create(comm, path, cmode, env_info, ncp);
    if (env_info != MPI_INFO_NULL) MPI_Info_free(&env_info);
    if (err != NC_NOERR) { /* fatal error */
        if (err == NC_EMULTIDEFINE_OMODE) err = NC_EMULTIDEFINE_CMODE;
        ncmpii_free_NC(ncp);
        DEBUG_RETURN_ERROR(err)
    }

    fSet(ncp->flags, NC_CREAT);

    /* initialize arrays storing pending non-blocking requests */
    ncp->numGetReqs = 0;
    ncp->numPutReqs = 0;
    ncp->get_list   = NULL;
    ncp->put_list   = NULL;

    /* initialize unlimited_id as no unlimited dimension yet defined */
    ncp->dims.unlimited_id = -1;

#ifndef SEARCH_NAME_LINEARLY
    for (i=0; i<HASH_TABLE_SIZE; i++) {
        /* initialize dim name lookup table */
        ncp->dims.nameT[i].num = 0;
        ncp->dims.nameT[i].list = NULL;
        /* initialize var name lookup table */
        ncp->vars.nameT[i].num = 0;
        ncp->vars.nameT[i].list = NULL;
    }
#endif

    ncp->ncid = ncid;

    *ncpp = (void*)ncp;

    return NC_NOERR;
}

/*----< ncmpii_open() >------------------------------------------------------*/
int
ncmpii_open(MPI_Comm    comm,
            const char *path,
            int         omode,
            int         ncid,
            MPI_Info    info,
            void       **ncpp)
{
    int i, err, status=NC_NOERR, safe_mode=0, mpireturn;
    char *env_str;
    MPI_Info   env_info=MPI_INFO_NULL;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;
    NC *ncp=NULL;
#ifndef SEARCH_NAME_LINEARLY
    NC_nametable *nameT;
#endif

#ifdef PNETCDF_DEBUG
    safe_mode = 1;
    /* this configure time setting will be overwritten by the run-time
     * environment variable PNETCDF_SAFE_MODE */
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

    /* path's validity and omode consistency have been checked in ncmpi_open()
     * in src/dispatchers/file.c */

    /* NC_DISKLESS is not supported yet */
    if (omode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* NC_MMAP is not supported yet */
    if (omode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* take hints from the environment variable PNETCDF_HINTS
     * a string of hints separated by ";" and each hint is in the
     * form of hint=value. E.g. cb_nodes=16;cb_config_list=*:6
     * If this environment variable is set, it  overrides any values that
     * were set by using calls to MPI_Info_set in the application code.
     */
    if (info != MPI_INFO_NULL) {
#ifdef HAVE_MPI_INFO_DUP
        mpireturn = MPI_Info_dup(info, &env_info);
        if (mpireturn != MPI_SUCCESS)
            DEBUG_RETURN_ERROR(ncmpii_handle_error(mpireturn, "MPI_Info_dup"))
#else
        printf("Warning: MPI info is ignored as MPI_Info_dup() is missing\n");
#endif
    }
    if ((env_str = getenv("PNETCDF_HINTS")) != NULL) {
        if (env_info == MPI_INFO_NULL)
            MPI_Info_create(&env_info); /* ignore error */

        char *env_str_cpy, *key;
        env_str_cpy = (char*) NCI_Malloc(strlen(env_str)+1);
        strcpy(env_str_cpy, env_str);
        key = strtok(env_str_cpy, ";");
        while (key != NULL) {
            char *val;
            val = strchr(key, '=');
            if (val == NULL) continue; /* ill-formed hint */
            *val = '\0';
            val++;
            /* printf("env hint: key=%s val=%s\n",key,val); */
            MPI_Info_set(env_info, key, val); /* override or add */
            key = strtok(NULL, ";");
        }
        NCI_Free(env_str_cpy);
    }

    /* get header chunk size from user info, if provided */
    if (env_info != MPI_INFO_NULL) {
        int flag;
        char value[MPI_MAX_INFO_VAL];
        MPI_Info_get(env_info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;
            chunksize = strtoll(value,NULL,10);
            if (errno != 0) chunksize = NC_DEFAULT_CHUNKSIZE;
        }
    }

    /* allocate NC file object */
    ncp = ncmpii_new_NC(&chunksize);
    if (ncp == NULL) {
        if (env_info != MPI_INFO_NULL) MPI_Info_free(&env_info);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    ncp->safe_mode = safe_mode;
    ncp->old       = NULL;
#ifdef ENABLE_SUBFILING
    ncp->subfile_mode = 1;
    if (env_info != MPI_INFO_NULL) {
        int flag;
        char value[MPI_MAX_INFO_VAL];
        MPI_Info_get(env_info, "pnetcdf_subfiling", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag && strcasecmp(value, "disable") == 0)
            ncp->subfile_mode = 0;
    }
    ncp->ncid_sf   = -1;
    ncp->nc_num_subfiles = 0;
#endif

    /* open the file in parallel */
    err = ncmpiio_open(comm, path, omode, env_info, ncp);
    if (err != NC_NOERR) { /* fatal error */
        if (env_info != MPI_INFO_NULL) MPI_Info_free(&env_info);
        ncmpii_free_NC(ncp);
        if (status == NC_NOERR) status = err;
        DEBUG_RETURN_ERROR(status) /* return the 1st error encountered */
    }

    /* PnetCDF's default mode is no fill */
    fSet(ncp->flags, NC_NOFILL);

    /* read header from file into an NC object pointed by ncp */
    err = ncmpii_hdr_get_NC(ncp);
    if (err != NC_NOERR) { /* fatal error */
        ncmpiio_close(ncp->nciop, 0);
        if (env_info != MPI_INFO_NULL) MPI_Info_free(&env_info);
        ncmpii_free_NC(ncp);
        if (status == NC_NOERR) status = err;
        DEBUG_RETURN_ERROR(status) /* return the 1st error encountered */
    }

    /* initialize arrays storing pending non-blocking requests */
    ncp->numGetReqs = 0;
    ncp->numPutReqs = 0;
    ncp->get_list   = NULL;
    ncp->put_list   = NULL;

#ifdef ENABLE_SUBFILING
    if (ncp->subfile_mode) {
        /* check attr for subfiles */
        err = ncmpii_get_att(ncp, NC_GLOBAL, "num_subfiles",
                             &ncp->nc_num_subfiles, NC_INT);
        if (err == NC_NOERR && ncp->nc_num_subfiles > 1) {
            /* ignore error NC_ENOTATT if this attribute is not defined */
            int nvars;

            err = ncmpii_inq(ncp, NULL, &nvars, NULL, NULL);

            if (status == NC_NOERR) status = err;

            for (i=0; i<nvars; i++) {
                err = ncmpii_get_att(ncp, i, "num_subfiles",
                                     &ncp->vars.value[i]->num_subfiles, NC_INT);
                if (err == NC_ENOTATT) continue;
                if (err != NC_NOERR && status == NC_NOERR) { /* other error */
                    status = err;
                    continue;
                }

                if (ncp->vars.value[i]->num_subfiles > 1) {
                    err = ncmpii_get_att(ncp, i, "ndims_org",
                                         &ncp->vars.value[i]->ndims_org, NC_INT);
                    if (status == NC_NOERR) status = err;
                }
            }

            if (ncp->nc_num_subfiles > 1) {
                err = ncmpii_subfile_open(ncp, &ncp->ncid_sf);
                if (status == NC_NOERR) status = err;
            }
        }
    }
    else
        ncp->nc_num_subfiles = 0;
#endif

    if (env_info != info) MPI_Info_free(&env_info);

    /* update the total number of record variables */
    ncp->vars.num_rec_vars = 0;
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.num_rec_vars += IS_RECVAR(ncp->vars.value[i]);

#ifndef SEARCH_NAME_LINEARLY
    /* initialize dim name lookup table */
    nameT = ncp->dims.nameT;
    memset(nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

    /* populate dim name lookup table */
    for (i=0; i<ncp->dims.ndefined; i++) {
        /* hash the dim name into a key for name lookup */
        int key = HASH_FUNC(ncp->dims.value[i]->name->cp);
        nameT = &ncp->dims.nameT[key];
        if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
            nameT->list = (int*) NCI_Realloc(nameT->list,
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) * SIZEOF_INT);
        nameT->list[nameT->num] = i;
        nameT->num++;
    }

    /* initialize var name lookup table */
    nameT = ncp->vars.nameT;
    memset(nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

    /* populate var name lookup table */
    for (i=0; i<ncp->vars.ndefined; i++) {
        /* hash the var name into a key for name lookup */
        int key = HASH_FUNC(ncp->vars.value[i]->name->cp);
        nameT = &ncp->vars.nameT[key];
        if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
            nameT->list = (int*) NCI_Realloc(nameT->list,
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) * SIZEOF_INT);
        nameT->list[nameT->num] = i;
        nameT->num++;
    }
#endif

    ncp->ncid = ncid;

    *ncpp = (void*)ncp;

    return status;
}

/*----< ncmpii_redef() >-----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii_redef(void *ncdp)
{
    int err=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (NC_readonly(ncp)) DEBUG_RETURN_ERROR(NC_EPERM) /* read-only */
    /* if open mode is inconsistent, then this return might cause parallel
     * program to hang */

    /* cannot be in define mode, must enter from data mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* sync all metadata, including numrecs, if changed in independent mode.
     * also ensure exiting define mode always entering collective data mode
     */
    if (NC_indep(ncp))
        ncmpiio_end_indep_data(ncp);

    if (NC_doFsync(ncp)) { /* re-read the header from file */
        err = ncmpii_read_NC(ncp);
        if (err != NC_NOERR) return err;
    }

    /* duplicate a header to be used in enddef() for checking if header grows */
    ncp->old = ncmpii_dup_NC(ncp);
    if (ncp->old == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* we are now entering define mode */
    fSet(ncp->flags, NC_INDEF);

    return NC_NOERR;
}

/*----< ncmpii_begin_indep_data() >------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii_begin_indep_data(void *ncdp)
{
    int err=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (NC_indef(ncp))  /* must not be in define mode */
        DEBUG_RETURN_ERROR(NC_EINDEFINE)

    if (NC_indep(ncp))  /* already in indep data mode */
        return NC_NOERR;

    /* we need no MPI_File_sync() here. If users want a stronger data
     * consistency, they can call ncmpi_sync()
     */
#if 0 && !defined(DISABLE_FILE_SYNC)
    if (!NC_readonly(ncp) && NC_collectiveFhOpened(ncp->nciop)) {
        /* calling file sync for those already open the file */
        int err, mpireturn;
        /* MPI_File_sync() is collective */
        TRACE_IO(MPI_File_sync)(ncp->nciop->collective_fh);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_File_sync");
            if (err == NC_NOERR) return err;
        }
        TRACE_COMM(MPI_Barrier)(ncp->nciop->comm);
    }
#endif

    fSet(ncp->flags, NC_INDEP);

    err = ncmpii_check_mpifh(ncp, 0);

    return err;
}

/*----< ncmpii_end_indep_data() >--------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii_end_indep_data(void *ncdp)
{
    NC *ncp = (NC*)ncdp;

    if (!NC_indep(ncp)) /* must be in independent data mode */
        DEBUG_RETURN_ERROR(NC_ENOTINDEP)

    return ncmpiio_end_indep_data(ncp);
}

/*----< ncmpiio_end_indep_data() >-------------------------------------------*/
/* this function is called when:
 * 1. ncmpi_end_indep_data()
 * 2. ncmpi_redef() from independent data mode entering to define more
 * 3. ncmpii_close() when closing the file
 * This function is collective.
 */
int
ncmpiio_end_indep_data(NC *ncp)
{
    int status=NC_NOERR;

    if (!NC_readonly(ncp)) {
        if (ncp->vars.num_rec_vars > 0) {
            /* numrecs dirty bit may not be the same across all processes.
             * force sync in memory no matter if dirty or not.
             */
            set_NC_ndirty(ncp);
            status = ncmpiio_sync_numrecs(ncp, ncp->numrecs);
            /* the only possible dirty part of the header is numrecs */
        }

#ifndef DISABLE_FILE_SYNC
        /* calling file sync for those already open the file */
        if (NC_doFsync(ncp) && NC_independentFhOpened(ncp->nciop)) {
            int mpireturn;
            /* MPI_File_sync() is collective */
            TRACE_IO(MPI_File_sync)(ncp->nciop->independent_fh);
            if (mpireturn != MPI_SUCCESS) {
                int err = ncmpii_handle_error(mpireturn, "MPI_File_sync");
                if (status == NC_NOERR) status = err;
            }
            TRACE_COMM(MPI_Barrier)(ncp->nciop->comm);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_handle_error(mpireturn, "MPI_Barrier");
        }
#endif
    }

    fClr(ncp->flags, NC_INDEP);

    return status;
}

/*----< ncmpii_enddef() >----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii_enddef(void *ncdp)
{
    NC *ncp = (NC*)ncdp;

    if (!NC_indef(ncp)) /* must currently in define mode */
        DEBUG_RETURN_ERROR(NC_ENOTINDEFINE)

    return ncmpiio_enddef(ncp);
}

/*----< ncmpii__enddef() >---------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpii__enddef(void       *ncdp,
               MPI_Offset  h_minfree,
               MPI_Offset  v_align,
               MPI_Offset  v_minfree,
               MPI_Offset  r_align)
{
    int err=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (!NC_indef(ncp)) /* must currently in define mode */
        DEBUG_RETURN_ERROR(NC_ENOTINDEFINE)

    if (ncp->safe_mode) {
        int status, mpireturn;
        MPI_Offset root_args[4];

        /* check if h_minfree, v_align, v_minfree, and r_align are consistent
         * among all processes */
        root_args[0] = h_minfree;
        root_args[1] = v_align;
        root_args[2] = v_minfree;
        root_args[3] = r_align;
        TRACE_COMM(MPI_Bcast)(&root_args, 4, MPI_OFFSET, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Bcast");

        if (root_args[0] != h_minfree ||
            root_args[1] != v_align   ||
            root_args[2] != v_minfree ||
            root_args[3] != r_align)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Allreduce");

        if (status != NC_NOERR) return status;
    }

    return ncmpiio__enddef(ncp, h_minfree, v_align, v_minfree, r_align);
}

/*----< ncmpii_sync_numrecs() >-----------------------------------------------*/
/* this API is collective, but can be called in independent data mode.
 * Note numrecs is always sync-ed in memory and update in file in collective
 * data mode.
 */
int
ncmpii_sync_numrecs(void *ncdp)
{
    int err=NC_NOERR;
    NC *ncp=(NC*)ncdp;

    /* cannot be in define mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* check if we have defined record variables */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    if (!NC_indep(ncp)) /* in collective data mode, numrecs is always sync-ed */
        return NC_NOERR;
    else /* if called in independent mode, we force sync in memory */
        set_NC_ndirty(ncp);

    /* sync numrecs in memory and file */
    err = ncmpiio_sync_numrecs(ncp, ncp->numrecs);

#ifndef DISABLE_FILE_SYNC
    if (NC_doFsync(ncp)) { /* NC_SHARE is set */
        int mpierr, mpireturn;
        if (NC_indep(ncp)) {
            TRACE_IO(MPI_File_sync)(ncp->nciop->independent_fh);
        }
        else {
            TRACE_IO(MPI_File_sync)(ncp->nciop->collective_fh);
        }
        if (mpireturn != MPI_SUCCESS) {
            mpierr = ncmpii_handle_error(mpireturn, "MPI_File_sync");
            if (err == NC_NOERR) err = mpierr;
        }
        TRACE_COMM(MPI_Barrier)(ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Barrier");
    }
#endif
    return err;
}

/*----< ncmpii_sync() >------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpii_sync(void *ncdp)
{
    int err;
    NC *ncp = (NC*)ncdp;

    /* cannot be in define mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    if (NC_readonly(ncp))
        /* calling sync for file opened for read only means re-read header */
        return ncmpii_read_NC(ncp);

    /* the only part of header that can be dirty is numrecs (caused only by
     * independent APIs) */
    if (ncp->vars.num_rec_vars > 0 && NC_indep(ncp)) {
        /* sync numrecs in memory among processes and in file */
        set_NC_ndirty(ncp);
        err = ncmpiio_sync_numrecs(ncp, ncp->numrecs);
        if (err != NC_NOERR) return err;
    }

    /* calling MPI_File_sync() on both collective and independent handlers */
    return ncmpiio_sync(ncp->nciop);
}

/*----< ncmpii_abort() >-----------------------------------------------------*/
/* This API is a collective subroutine */
int
ncmpii_abort(void *ncdp)
{
   /*
    * In data mode, same as ncmpiio_close.
    * In define mode, descard new definition.
    * If file is just created, remove the file.
    */
    int status=NC_NOERR, err, doUnlink = 0;
    NC *ncp = (NC*)ncdp;

    /* delete the file if it is newly created by ncmpi_create() */
    doUnlink = NC_IsNew(ncp);

    if (ncp->old != NULL) {
        /* a plain redef, not a create */
        assert(!NC_IsNew(ncp));
        assert(fIsSet(ncp->flags, NC_INDEF));
        ncmpii_free_NC(ncp->old);
        ncp->old = NULL;
        fClr(ncp->flags, NC_INDEF);
    }

    if (!doUnlink) {
        if (!NC_readonly(ncp) &&  /* file is open for write */
             NC_indep(ncp)) {     /* in independent data mode */
            status = ncmpiio_end_indep_data(ncp); /* sync header */
        }

        if (NC_doFsync(ncp)) {
            err = ncmpiio_sync(ncp->nciop); /* calling MPI_File_sync() */
            if (status == NC_NOERR ) status = err;
        }
    }

    /* close the file */
    err = ncmpiio_close(ncp->nciop, doUnlink);
    if (status == NC_NOERR ) status = err;

    ncp->nciop = NULL;

    /* free up space occupied by the header metadata */
    ncmpii_free_NC(ncp);

    return status;
}

/*----< ncmpi_delete() >-----------------------------------------------------*/
/* doesn't do anything to release resources, so call ncmpi_close before calling
 * this function.
 *
 * filename: the name of the file we will remove.
 * info: MPI info object, in case underlying file system needs hints.
 */
int
ncmpi_delete(const char *filename,
             MPI_Info    info)
{
    int err=NC_NOERR, mpireturn;

    TRACE_IO(MPI_File_delete)((char*)filename, info);
    if (mpireturn != MPI_SUCCESS)
        err = ncmpii_handle_error(mpireturn, "MPI_File_delete");
    return err;
}

/* End Of Dataset Functions */

/*----< ncmpii_check_mpifh() >-----------------------------------------------*/
int
ncmpii_check_mpifh(NC  *ncp,
                   int  collective)
{
    int mpireturn;

    if (collective && NC_indep(ncp)) /* collective handle but in indep mode */
        DEBUG_RETURN_ERROR(NC_EINDEP)

    if (!collective && !NC_indep(ncp)) /* indep handle but in collective mode */
        DEBUG_RETURN_ERROR(NC_ENOTINDEP)

    if (collective && !NC_collectiveFhOpened(ncp->nciop)) {
        TRACE_IO(MPI_File_open)(ncp->nciop->comm, (char*)ncp->nciop->path,
                                ncp->nciop->mpiomode, ncp->nciop->mpiinfo,
                                &ncp->nciop->collective_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_open");

        set_NC_collectiveFh(ncp->nciop);
    }
    else if (!collective && !NC_independentFhOpened(ncp->nciop)) {
        TRACE_IO(MPI_File_open)(MPI_COMM_SELF, (char*)ncp->nciop->path,
                                ncp->nciop->mpiomode, ncp->nciop->mpiinfo,
                                &ncp->nciop->independent_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_open");

        set_NC_independentFh(ncp->nciop);
    }

    return NC_NOERR;
}

/*----< ncmpii_inq_misc() >--------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpii_inq_misc(void       *ncdp,
                int        *pathlen,
                char       *path,
                int        *num_fix_varsp,
                int        *num_rec_varsp,
                int        *striping_size,
                int        *striping_count,
                MPI_Offset *header_size,
                MPI_Offset *header_extent,
                MPI_Offset *recsize,
                MPI_Offset *put_size,
                MPI_Offset *get_size,
                MPI_Info   *info_used,
                int        *nreqs,
                MPI_Offset *usage,
                MPI_Offset *buf_size)
{
    int i, flag, mpireturn;
    char value[MPI_MAX_INFO_VAL];
    NC *ncp=(NC*)ncdp;

    /* Get the file pathname which was used to open/create the ncid's file.
     * path must already be allocated. Ignored if NULL */
    if (ncp->nciop->path == NULL) {
        if (pathlen != NULL) *pathlen = 0;
        if (path    != NULL) *path = '\0';
    } else {
        if (pathlen != NULL) *pathlen = (int)strlen(ncp->nciop->path);
        if (path    != NULL) strcpy(path, ncp->nciop->path);
    }

    /* obtain the number of fixed-size variables */
    if (num_fix_varsp != NULL) {
        if (NC_indef(ncp)) {
            /* if in define mode, recalculate the number of record variables */
            *num_fix_varsp = 0;
            for (i=0; i<ncp->vars.ndefined; i++)
                *num_fix_varsp += IS_RECVAR(ncp->vars.value[i]);
        }
        else
            *num_fix_varsp = ncp->vars.num_rec_vars;

        /* no. fixed-size == ndefined - no. record variables */
        *num_fix_varsp = ncp->vars.ndefined - *num_fix_varsp;
    }

    /* obtain the number of record variables */
    if (num_rec_varsp != NULL) {
        if (NC_indef(ncp)) {
            /* if in define mode, recalculate the number of record variables */
            *num_rec_varsp = 0;
            for (i=0; i<ncp->vars.ndefined; i++)
                *num_rec_varsp += IS_RECVAR(ncp->vars.value[i]);
        }
        else
            *num_rec_varsp = ncp->vars.num_rec_vars;
    }

    /* obtain file (system) striping settings, striping size and count, if they
     * are available from MPI-IO hint. Otherwise, 0s are returned.
     */
    if (striping_size != NULL) {
        MPI_Info_get(ncp->nciop->mpiinfo, "striping_unit", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        *striping_size = 0;
        if (flag) {
            errno = 0;
            *striping_size = (int)strtol(value,NULL,10);
            if (errno != 0) *striping_size = 0;
        }
    }

    if (striping_count != NULL) {
        MPI_Info_get(ncp->nciop->mpiinfo, "striping_factor", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        *striping_count = 0;
        if (flag) {
            errno = 0;
            *striping_count = (int)strtol(value,NULL,10);
            if (errno != 0) *striping_count = 0;
        }
    }

    /* the amount of writes, in bytes, committed to file system so far */
    if (put_size != NULL) *put_size = ncp->nciop->put_size;

    /* the amount of reads, in bytes, obtained from file system so far */
    if (get_size != NULL) *get_size = ncp->nciop->get_size;

    if (recsize != NULL) *recsize = ncp->recsize;

    if (header_size != NULL) *header_size = ncp->xsz;

    if (header_extent != NULL) *header_extent = ncp->begin_var;

    if (info_used != NULL) {
#ifdef HAVE_MPI_INFO_DUP
        mpireturn = MPI_Info_dup(ncp->nciop->mpiinfo, info_used);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Info_dup");
#else
        mpireturn = MPI_File_get_info(ncp->nciop->collective_fh, info_used);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_File_get_info");
#endif

        sprintf(value, "%lld", ncp->nciop->hints.h_align);
        MPI_Info_set(*info_used, "nc_header_align_size", value);

        sprintf(value, "%lld", ncp->nciop->hints.v_align);
        MPI_Info_set(*info_used, "nc_var_align_size", value);

        sprintf(value, "%lld", ncp->nciop->hints.r_align);
        MPI_Info_set(*info_used, "nc_record_align_size", value);

        sprintf(value, "%lld", ncp->nciop->hints.header_read_chunk_size);
        MPI_Info_set(*info_used, "nc_header_read_chunk_size", value);

#ifdef ENABLE_SUBFILING
        sprintf(value, "%d", ncp->nciop->hints.subfile_mode);
        MPI_Info_set(*info_used, "pnetcdf_subfiling", value);
        sprintf(value, "%d", ncp->nciop->hints.num_subfiles);
        MPI_Info_set(*info_used, "nc_num_subfiles", value);
#endif
    }

    if (nreqs != NULL) {
        /* cannot just use *nreqs = ncp->numGetReqs + ncp->numPutReqs;
         * because some request IDs are repeated, such as record variables and
         * varn requests
         */
        *nreqs = 0;
        for (i=0; i<ncp->numGetReqs; i++) {
            if (i > 0 && ncp->get_list[i].id == ncp->get_list[i-1].id)
                continue;
            (*nreqs)++;
        }
        for (i=0; i<ncp->numPutReqs; i++) {
            if (i > 0 && ncp->put_list[i].id == ncp->put_list[i-1].id)
                continue;
            (*nreqs)++;
        }
    }

    if (usage != NULL) {
        /* check if the buffer has been previously attached */
        if (ncp->abuf == NULL) DEBUG_RETURN_ERROR(NC_ENULLABUF)
        /* return the current usage in bytes */
        *usage = ncp->abuf->size_used;
    }

    if (buf_size != NULL) {
        /* check if the buffer has been previously attached */
        if (ncp->abuf == NULL) DEBUG_RETURN_ERROR(NC_ENULLABUF)
        /* return the current usage in bytes */
        *buf_size = ncp->abuf->size_allocated;
    }

    return NC_NOERR;
}

