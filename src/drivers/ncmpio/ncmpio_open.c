/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_open() : dispatcher->open()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy() */
#ifdef HAVE_ACCESS
#include <unistd.h>  /* access() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"
#ifdef ENABLE_SUBFILING
#include "ncmpio_subfile.h"
#endif

/*----< ncmpio_open() >------------------------------------------------------*/
int
ncmpio_open(MPI_Comm     comm,
            const char  *path,
            int          omode,
            int          ncid,
            MPI_Info     user_info, /* user's and env info combined */
            void       **ncpp)
{
    char *env_str;
    int mpiomode, err, status=NC_NOERR, mpireturn;
    MPI_File fh;
    MPI_Info info_used;
    NC *ncp=NULL;

    *ncpp = NULL;

    /* Note path's validity and omode consistency have been checked in
     * ncmpi_open() in src/dispatchers/file.c and
     * path consistency will be done in MPI_File_open */

    /* First, check whether omode is valid or supported ---------------------*/
    /* NC_DISKLESS is not supported yet */
    if (omode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* NC_MMAP is not supported yet */
    if (omode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

#if 0 && defined(HAVE_ACCESS)
    if (mpiomode == MPI_MODE_RDONLY) { /* file should already exit */
        int rank, file_exist;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            if (access(path, F_OK) == 0) file_exist = 1;
            else                         file_exist = 0;
        }
        TRACE_COMM(MPI_Bcast)(&file_exist, 1, MPI_INT, 0, comm);
        if (!file_exist) DEBUG_RETURN_ERROR(NC_ENOENT)
    }
#endif

    /* open file collectively ---------------------------------------------- */
    mpiomode = fIsSet(omode, NC_WRITE) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;

    TRACE_IO(MPI_File_open)(comm, (char *)path, mpiomode, user_info, &fh);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_File_open");

    /* get the file info used/modified by MPI-IO */
    mpireturn = MPI_File_get_info(fh, &info_used);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_File_get_info");

    /* Now the file has been successfully opened, allocate/set NC object */

    /* path's validity and omode consistency have been checked in ncmpi_open()
     * in src/dispatchers/file.c */

    /* allocate buffer for header object NC */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* PnetCDF default fill mode is no fill */
    fClr(ncp->flags, NC_MODE_FILL);
    if (!fIsSet(omode, NC_WRITE)) fSet(ncp->flags, NC_MODE_RDONLY);

    ncp->ncid = ncid;

    /* chunk size for reading header (set default before check hints) */
    ncp->chunk = NC_DEFAULT_CHUNKSIZE;

    /* buffer to pack noncontiguous user buffers when calling wait() */
    ncp->ibuf_size = NC_DEFAULT_IBUF_SIZE;

    /* Extract PnetCDF specific I/O hints from user_info and set default hint
     * values into info_used. Note some MPI libraries, such as MPICH 3.3.1 and
     * priors fail to preserve user hints that are not recogniozed by the MPI
     * libraries.
     */
    ncmpio_set_pnetcdf_hints(ncp, user_info, info_used);

    ncp->iomode         = omode;
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

    /* read header from file into NC object pointed by ncp -------------------*/
    err = ncmpio_hdr_get_NC(ncp);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) { /* fatal error */
        ncmpio_close_files(ncp, 0);
        ncmpio_free_NC(ncp);
        return err;
    }

#ifdef ENABLE_SUBFILING
    if (ncp->subfile_mode) {
        /* check subfiling attribute */
        err = ncmpio_get_att(ncp, NC_GLOBAL, "_PnetCDF_SubFiling.num_subfiles",
                             &ncp->num_subfiles, MPI_INT);
        if (err == NC_NOERR && ncp->num_subfiles > 1) {
            int i;
            /* ignore error NC_ENOTATT if this attribute is not defined */
            for (i=0; i<ncp->vars.ndefined; i++) {
                /* variables may have different numbers of subfiles */
                err = ncmpio_get_att(ncp, i, "_PnetCDF_SubFiling.num_subfiles",
                             &ncp->vars.value[i]->num_subfiles,MPI_INT);
                if (err == NC_ENOTATT) continue;
                if (err != NC_NOERR) return err;
                if (ncp->vars.value[i]->num_subfiles > 1) {
                    /* find the orginal ndims of variable i */
                    err = ncmpio_get_att(ncp,i,"_PnetCDF_SubFiling.ndims_org",
                                 &ncp->vars.value[i]->ndims_org,MPI_INT);
                    if (err != NC_NOERR) return err;
                    ncp->vars.value[i]->dimids_org = (int*) NCI_Malloc(
                              ncp->vars.value[i]->ndims_org * SIZEOF_INT);
                    err = ncmpio_get_att(ncp,i,"_PnetCDF_SubFiling.dimids_org",
                              ncp->vars.value[i]->dimids_org, MPI_INT);
                    if (err != NC_NOERR) return err;
                }
            }
            /* open subfile */
            err = ncmpio_subfile_open(ncp);
            if (err != NC_NOERR) return err;
        }
        else ncp->num_subfiles = 0;
    }
    else
        ncp->num_subfiles = 0;
#endif

#ifndef SEARCH_NAME_LINEARLY
    /* initialize and populate name lookup tables ---------------------------*/
    ncmpio_hash_table_populate_NC_dim(&ncp->dims);
    ncmpio_hash_table_populate_NC_var(&ncp->vars);
    ncmpio_hash_table_populate_NC_attr(ncp);
#endif

    *ncpp = (void*)ncp;

    return status;
}

