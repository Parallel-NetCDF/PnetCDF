/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_create() : dispatcher->create()
 * ncmpi_open()   : dispatcher->open()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy(), strchr() */
#ifdef HAVE_ACCESS
#include <unistd.h>  /* access() */
#endif
#include <errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< ncmpio_create() >----------------------------------------------------*/
int
ncmpio_create(MPI_Comm     comm,
              const char  *path,
              int          cmode,
              int          ncid,
              MPI_Info     info, /* user's info and env info combined */
              void       **ncpp)
{
    char *env_str;
    int i, rank, mpiomode, err, mpireturn, default_format;
    MPI_File fh;
    MPI_Comm dup_comm;
    MPI_Info info_used;
    NC *ncp=NULL;

    *ncpp = NULL;

    /* Note path's validity and cmode consistency have been checked in
     * ncmpi_create() in src/dispatchers/file.c and
     * path consistency will be done in MPI_File_open */

    /* First, check whether cmode is valid or supported ---------------------*/
    /* NC_DISKLESS is not supported yet */
    if (cmode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* NC_MMAP is not supported yet */
    if (cmode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* Get default format, in case cmode does not include either
     * NC_64BIT_OFFSET or NC_64BIT_DATA */
    ncmpi_inq_default_format(&default_format);

#if SIZEOF_MPI_OFFSET <  8
    /* check cmode */
    if (fIsSet(cmode, NC_64BIT_DATA)     || fIsSet(cmode, NC_64BIT_OFFSET)   ||
        default_format == NC_FORMAT_CDF5 || default_format == NC_FORMAT_CDF2) {
        /* unlike serial netcdf, we will not bother to support
         * NC_64BIT_OFFSET on systems with off_t smaller than 8 bytes.
         * serial netcdf has proven it's possible if datasets are small, but
         * that's a hassle we don't want to worry about */
        DEBUG_RETURN_ERROR(NC_ESMALL)
    }
#endif

    /* It is illegal to have both NC_64BIT_OFFSET & NC_64BIT_DATA */
    if ((cmode & (NC_64BIT_OFFSET|NC_64BIT_DATA)) ==
                 (NC_64BIT_OFFSET|NC_64BIT_DATA))
        DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* Handle file clobber --------------------------------------------------*/
    MPI_Comm_rank(comm, &rank);

    mpiomode = MPI_MODE_RDWR | MPI_MODE_CREATE;

    if (fIsSet(cmode, NC_NOCLOBBER)) {
        /* check if file exists: NC_EEXIST is returned if the file * already
         * exists and NC_NOCLOBBER mode is used in ncmpi_create */
#ifdef HAVE_ACCESS
        int file_exist;
        /* if access() is available, use it to check whether file already exists
         * rank 0 calls access() and broadcasts file_exist */
        if (rank == 0) {
            /* remove the file system type prefix name if there is any.
             * For example, path=="lustre:/home/foo/testfile.nc",
             * use "/home/foo/testfile.nc" when calling access()
             */
            char *filename = strchr(path, ':');
            if (filename == NULL) filename = (char*)path; /* no prefix */
            else                  filename++;

            if (access(filename, F_OK) == 0) file_exist = 1;
            else                             file_exist = 0;
        }
        TRACE_COMM(MPI_Bcast)(&file_exist, 1, MPI_INT, 0, comm);
        if (file_exist) DEBUG_RETURN_ERROR(NC_EEXIST)
#else
        /* add MPI_MODE_EXCL mode for MPI_File_open to check file existence */
        fSet(mpiomode, MPI_MODE_EXCL);
#endif
    }
    else { /* NC_CLOBBER is the default mode in create */
        /* rank 0 deletes the file and ignores error code.
         * Note calling MPI_File_set_size is expensive as it calls truncate()
         */
        if (rank == 0) {
#ifdef HAVE_UNLINK
            char *filename = strchr(path, ':');
            if (filename == NULL) filename = (char*)path; /* no prefix */
            else                  filename++;

            err = unlink(filename);
            if (err < 0 && errno != ENOENT) /* ignore ENOENT: file not exist */
                DEBUG_ASSIGN_ERROR(err, NC_EFILE) /* other error */
            else
                err = NC_NOERR;
#else
            err = NC_NOERR;
            TRACE_IO(MPI_File_delete)((char*)path, MPI_INFO_NULL);
            if (mpireturn != MPI_SUCCESS) {
                int errorclass;
                MPI_Error_class(mpireturn, &errorclass);
                if (errorclass != MPI_ERR_NO_SUCH_FILE) /* ignore this error */
                    err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_delete");
            }
#endif
        }
        /* all processes must wait here until file deletion is completed */
        TRACE_COMM(MPI_Bcast)(&err, 1, MPI_INT, 0, comm);
        if (err != NC_NOERR) return err;
    }

    /* create file collectively -------------------------------------------- */
    TRACE_IO(MPI_File_open)(comm, (char *)path, mpiomode, info, &fh);
    if (mpireturn != MPI_SUCCESS) {
#ifndef HAVE_ACCESS
        if (fIsSet(cmode, NC_NOCLOBBER)) {
            /* This is the case when NC_NOCLOBBER is used in file creation and
             * function access() is not available. MPI_MODE_EXCL is set in open
             * mode. When MPI_MODE_EXCL is used and the file already exists,
             * MPI-IO should return error class MPI_ERR_FILE_EXISTS. But, some
             * MPI-IO implementations (older ROMIO) do not correctly return
             * this error class. In this case, we can do the followings: check
             * errno to see if it set to EEXIST. Note usually rank 0 makes the
             * file open call and can be the only one having errno set.
             */
            TRACE_COMM(MPI_Bcast)(&errno, 1, MPI_INT, 0, comm);
            if (errno == EEXIST) DEBUG_RETURN_ERROR(NC_EEXIST)
        }
#endif
        return ncmpii_error_mpi2nc(mpireturn, "MPI_File_open");
        /* for NC_NOCLOBBER, MPI_MODE_EXCL was added to mpiomode. If the file
         * already exists, MPI-IO should return error class MPI_ERR_FILE_EXISTS
         * which PnetCDF will return error code NC_EEXIST. This is checked
         * inside of ncmpii_error_mpi2nc()
         */
    }

    /* duplicate communicator as user may free it later */
    mpireturn = MPI_Comm_dup(comm, &dup_comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_dup");

    /* get the file info actually used by MPI-IO (may alter user's info) */
    mpireturn = MPI_File_get_info(fh, &info_used);
    if (mpireturn != MPI_SUCCESS) {
        MPI_Comm_free(&dup_comm);
        return ncmpii_error_mpi2nc(mpireturn, "MPI_File_get_info");
    }

    /* Now the file has been successfully created, allocate/set NC object */

    /* allocate buffer for header object NC and initialize its contents */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* set the file format version based on the create mode, cmode */
         if (fIsSet(cmode, NC_64BIT_DATA))   ncp->format = 5;
    else if (fIsSet(cmode, NC_64BIT_OFFSET)) ncp->format = 2;
    else {
             if (default_format == NC_FORMAT_CDF5) ncp->format = 5;
        else if (default_format == NC_FORMAT_CDF2) ncp->format = 2;
        else                                       ncp->format = 1;
    }

    fClr(ncp->flags, NC_MODE_RDONLY); /* create automatically enter write mode */
    fSet(ncp->flags, NC_MODE_CREATE);
    fSet(ncp->flags, NC_MODE_DEF);    /* create automatically enter define mode */
    fClr(ncp->flags, NC_MODE_FILL);   /* PnetCDF default mode is no fill */

    ncp->ncid         = ncid;
    ncp->safe_mode    = 0;
    /* initialize arrays storing pending non-blocking requests */
    ncp->numGetReqs   = 0;
    ncp->numPutReqs   = 0;
#ifdef ENABLE_SUBFILING
    ncp->subfile_mode = 0;
    ncp->num_subfiles = 0;
    ncp->ncp_sf       = NULL; /* pointer to subfile NC object */
#endif

    ncp->chunk        = NC_DEFAULT_CHUNKSIZE;
    ncp->h_align      = 0; /* value 0 indicates the hint is not set */
    ncp->v_align      = 0;
    ncp->r_align      = 0;
    ncp->h_minfree    = 0;
    ncp->v_minfree    = 0;

    /* extract I/O hints from user info */
    ncmpio_set_pnetcdf_hints(ncp, info);

    /* find the true header size (not-yet aligned) */
    ncp->xsz          = ncmpio_hdr_len_NC(ncp);

    ncp->get_list     = NULL;
    ncp->put_list     = NULL;
    ncp->abuf         = NULL;
    ncp->old          = NULL;

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

    /* For file create, ignore if NC_NOWRITE set in cmode by user */
    ncp->iomode         = cmode | NC_WRITE;
    ncp->comm           = dup_comm;
    ncp->mpiinfo        = info_used;
    ncp->mpiomode       = mpiomode;
    ncp->put_size       = 0;
    ncp->get_size       = 0;
    ncp->collective_fh  = fh;
    ncp->independent_fh = MPI_FILE_NULL;
    ncp->path = (char*) NCI_Malloc(strlen(path) + 1);
    strcpy(ncp->path, path);

#ifdef PNETCDF_DEBUG
    /* PNETCDF_DEBUG is set at configure time, which will be overwritten by
     * the run-time environment variable PNETCDF_SAFE_MODE */
    ncp->safe_mode = 1;
#endif
    /* If environment variable PNETCDF_SAFE_MODE is set to 1, then we perform
     * a strict consistent test, i.e. arguments used in def_dim/def_var APIs
     */
    if ((env_str = getenv("PNETCDF_SAFE_MODE")) != NULL) {
        if (*env_str == '0') ncp->safe_mode = 0;
        else                 ncp->safe_mode = 1;
        /* if PNETCDF_SAFE_MODE is set but without a value, *env_str can
         * be '\0' (null character). In this case, safe_mode is enabled */
    }

    *ncpp = (void*)ncp;

    return NC_NOERR;
}

