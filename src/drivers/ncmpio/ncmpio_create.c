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
              MPI_Info     user_info, /* user's and env info combined */
              void       **ncpp)
{
    char *env_str, *filename;
    int rank, mpiomode, err, mpireturn, default_format, file_exist = 1;
    MPI_File fh;
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

    /* Check cmode for other illegal flags already done in dispatcher layer */

    /* Get default format, in case cmode does not include either
     * NC_64BIT_OFFSET or NC_64BIT_DATA */
    ncmpi_inq_default_format(&default_format);

    /* Handle file clobber --------------------------------------------------*/
    MPI_Comm_rank(comm, &rank);

    mpiomode = MPI_MODE_RDWR | MPI_MODE_CREATE;

    /* remove the file system type prefix name if there is any.  For example,
     * path=="lustre:/home/foo/testfile.nc", use "/home/foo/testfile.nc" when
     * calling access()
     */
    filename = strchr(path, ':');
    if (filename == NULL) filename = (char*)path; /* no prefix */
    else                  filename++;

#ifdef HAVE_ACCESS
    /* if access() is available, use it to check whether file already exists
     * rank 0 calls access() and broadcasts file_exist */
    if (rank == 0) {
        if (access(filename, F_OK) == -1) file_exist = 0;
        errno = 0; /* reset errno */
    }
#endif

    if (fIsSet(cmode, NC_NOCLOBBER)) {
        /* check if file exists: NC_EEXIST is returned if the file already
         * exists and NC_NOCLOBBER mode is used in ncmpi_create */
#ifdef HAVE_ACCESS
        TRACE_COMM(MPI_Bcast)(&file_exist, 1, MPI_INT, 0, comm);
        if (file_exist) DEBUG_RETURN_ERROR(NC_EEXIST)
#else
        /* add MPI_MODE_EXCL mode for MPI_File_open to check file existence */
        fSet(mpiomode, MPI_MODE_EXCL);
        errno = 0; /* reset errno, as MPI_File_open may change it */
#endif
    }
    else { /* NC_CLOBBER is the default mode in create */
        /* rank 0 truncates or deletes the file and ignores error code.
         * Note calling MPI_File_set_size is expensive as it calls truncate()
         */
        err = NC_NOERR;
        if (rank == 0 && file_exist) {
#ifdef HAVE_TRUNCATE
            err = truncate(filename, 0);
            if (err < 0 && errno != ENOENT) /* ignore ENOENT: file not exist */
                DEBUG_ASSIGN_ERROR(err, NC_EFILE) /* other error */
            else
                err = NC_NOERR;
#else
            /* call MPI_File_set_size() to truncate the file. Note this can
             * be expensive.
             */
            err = NC_NOERR;
            TRACE_IO(MPI_File_open)(MPI_COMM_SELF, (char *)path, MPI_MODE_RDWR,
                                    MPI_INFO_NULL, &fh);
            if (mpireturn != MPI_SUCCESS) {
                int errorclass;
                MPI_Error_class(mpireturn, &errorclass);
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_open");
            }
            else {
                TRACE_IO(MPI_File_set_size)(fh, 0);
                if (mpireturn != MPI_SUCCESS) {
                    int errorclass;
                    MPI_Error_class(mpireturn, &errorclass);
                    err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_set_size");
                }
                else {
                    TRACE_IO(MPI_File_close)(&fh);
                    if (mpireturn != MPI_SUCCESS) {
                        int errorclass;
                        MPI_Error_class(mpireturn, &errorclass);
                        err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_close");
                    }
                }
            }
#endif
            if (errno == ENOENT) errno = 0; /* reset errno */
        }
        /* all processes must wait here until file deletion is completed */
        TRACE_COMM(MPI_Bcast)(&err, 1, MPI_INT, 0, comm);
        if (err != NC_NOERR) return err;
    }

    /* create file collectively -------------------------------------------- */
    TRACE_IO(MPI_File_open)(comm, (char *)path, mpiomode, user_info, &fh);
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
    else
        /* reset errno, as MPI_File_open may change it, even for MPI_SUCCESS */
        errno = 0;

    /* get the I/O hints used/modified by MPI-IO */
    mpireturn = MPI_File_get_info(fh, &info_used);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_File_get_info");

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

    fSet(ncp->flags, NC_MODE_CREATE);
    /* create automatically enter write mode */
    fClr(ncp->flags, NC_MODE_RDONLY);
    /* create automatically enter define mode */
    fSet(ncp->flags, NC_MODE_DEF);
    /* PnetCDF default mode is no fill */
    fClr(ncp->flags, NC_MODE_FILL);

    ncp->ncid = ncid;

    /* chunk size for reading header, set to default before check hints */
    ncp->chunk = NC_DEFAULT_CHUNKSIZE;

    /* calculate the true header size (not-yet aligned)
     * No need to do this now.
     * ncp->xsz = ncmpio_hdr_len_NC(ncp);
     */

    /* initialize unlimited_id as no unlimited dimension yet defined */
    ncp->dims.unlimited_id = -1;

    /* buffer to pack noncontiguous user buffers when calling wait() */
    ncp->ibuf_size = NC_DEFAULT_IBUF_SIZE;

    /* Extract PnetCDF specific I/O hints from user_info and set default hint
     * values into info_used. Note some MPI libraries, such as MPICH 3.3.1 and
     * priors fail to preserve user hints that are not recogniozed by the MPI
     * libraries.
     */
    ncmpio_set_pnetcdf_hints(ncp, user_info, info_used);

    /* For file create, ignore if NC_NOWRITE set in cmode by user */
    ncp->iomode         = cmode | NC_WRITE;
    ncp->comm           = comm;  /* reuse comm duplicated in dispatch layer */
    ncp->mpiinfo        = info_used; /* is not MPI_INFO_NULL */
    ncp->mpiomode       = mpiomode;
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

