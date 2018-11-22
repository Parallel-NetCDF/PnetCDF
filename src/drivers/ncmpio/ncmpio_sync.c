/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_sync()         : dispatcher->sync()
 * ncmpi_sync_numrecs() : dispatcher->sync_numrecs()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset() */
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncx.h>
#include "ncmpio_NC.h"

/*----< ncmpio_file_sync() >-------------------------------------------------*/
/* This function must be called collectively, no matter if it is in collective
 * or independent data mode.
 */
int
ncmpio_file_sync(NC *ncp) {
#ifndef DISABLE_FILE_SYNC
    int mpireturn;

    if (ncp->independent_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync)(ncp->independent_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_File_sync");
    }

    /* ncp->collective_fh is never MPI_FILE_NULL as collective mode is
     * default in PnetCDF */
    TRACE_IO(MPI_File_sync)(ncp->collective_fh);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_File_sync");

    TRACE_COMM(MPI_Barrier)(ncp->comm);
#endif
    return NC_NOERR;
}

#define NC_NUMRECS_OFFSET 4

/*----< ncmpio_write_numrecs() >---------------------------------------------*/
/* root process writes the new record number into file.
 * This function is called by:
 * 1. ncmpio_sync_numrecs
 * 2. collective nonblocking wait API, if the new number of records is bigger
 */
int
ncmpio_write_numrecs(NC         *ncp,
                     MPI_Offset  new_numrecs)
{
    int rank, mpireturn, err;
    MPI_File fh;

    /* root process writes numrecs in file */
    MPI_Comm_rank(ncp->comm, &rank);
    if (rank > 0) return NC_NOERR;

    /* return now if there is no record variabled defined */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    fh = ncp->collective_fh;
    if (NC_indep(ncp))
        fh = ncp->independent_fh;

    if (new_numrecs > ncp->numrecs || NC_ndirty(ncp)) {
        int len;
        char pos[8], *buf=pos;
        MPI_Status mpistatus;

        /* update ncp->numrecs */
        if (new_numrecs > ncp->numrecs) ncp->numrecs = new_numrecs;

        if (ncp->format < 5) {
            if (ncp->numrecs != (int)ncp->numrecs)
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            len = X_SIZEOF_SIZE_T;
            err = ncmpix_put_uint32((void**)&buf, (uint)ncp->numrecs);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
        }
        else {
            len = X_SIZEOF_INT64;
            err = ncmpix_put_uint64((void**)&buf, (uint64)ncp->numrecs);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
        }
        /* ncmpix_put_xxx advances the 1st argument with size len */

#ifdef _USE_MPI_GET_COUNT
        /* explicitly initialize mpistatus object to 0. For zero-length read,
         * MPI_Get_count may report incorrect result for some MPICH version,
         * due to the uninitialized MPI_Status object passed to MPI-IO calls.
         * Thus we initialize it above to work around.
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));
#endif
        /* root's file view always includes the entire file header */
        TRACE_IO(MPI_File_write_at)(fh, NC_NUMRECS_OFFSET, (void*)pos, len,
                                    MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_write_at");
            if (err == NC_EFILE) DEBUG_RETURN_ERROR(NC_EWRITE)
        }
        else {
#ifdef _USE_MPI_GET_COUNT
            int put_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
            ncp->put_size += put_size;
#else
            ncp->put_size += len;
#endif
        }
    }
    return NC_NOERR;
}

/*----< ncmpio_sync_numrecs() >-----------------------------------------------*/
/* Synchronize the number of records in memory among all processes and write
 * numrecs to file.
 * This function is called by:
 * 1. ncmpi_sync_numrecs(): by the user
 * 2. ncmpi_sync(): by the user
 * 3. ncmpi_end_indep_data(): exit from independent data mode
 * 4. all blocking collective put APIs when writing record variables
 * 5. ncmpi_close(): file close and currently in independent data mode
 *
 * This API is collective, but can be called in independent data mode.
 * Note numrecs is always sync-ed in memory and update in file in collective
 * data mode.
 */
int
ncmpio_sync_numrecs(void *ncdp)
{
    int status=NC_NOERR, mpireturn;
    MPI_Offset max_numrecs;
    NC *ncp=(NC*)ncdp;

    /* cannot be in define mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* check if we have defined record variables */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* check write permission */
    if (NC_readonly(ncp)) DEBUG_RETURN_ERROR(NC_EPERM)

    if (!NC_indep(ncp)) /* in collective data mode, numrecs is always sync-ed */
        return NC_NOERR;
    else /* if called in independent mode, we force sync in memory */
        set_NC_ndirty(ncp);

    /* return now if there is no record variabled defined */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* find the max numrecs among all processes
     * Note max numrecs may be smaller than some process's ncp->numrecs
     */
    TRACE_COMM(MPI_Allreduce)(&ncp->numrecs, &max_numrecs, 1, MPI_OFFSET,
                              MPI_MAX, ncp->comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");

    /* root process writes max_numrecs to file */
    status = ncmpio_write_numrecs(ncp, max_numrecs);

    if (ncp->safe_mode == 1) {
        /* broadcast root's status, because only root writes to the file */
        int root_status = status;
        TRACE_COMM(MPI_Bcast)(&root_status, 1, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        /* root's write has failed, which is serious */
        if (root_status == NC_EWRITE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
    }

    /* update numrecs in all processes's memory */
    ncp->numrecs = max_numrecs;

    /* clear numrecs dirty bit */
    fClr(ncp->flags, NC_NDIRTY);

#ifndef DISABLE_FILE_SYNC
    if (NC_doFsync(ncp)) { /* NC_SHARE is set */
        int mpierr, mpireturn;
        if (NC_indep(ncp)) {
            TRACE_IO(MPI_File_sync)(ncp->independent_fh);
        }
        else {
            TRACE_IO(MPI_File_sync)(ncp->collective_fh);
        }
        if (mpireturn != MPI_SUCCESS) {
            mpierr = ncmpii_error_mpi2nc(mpireturn, "MPI_File_sync");
            if (status == NC_NOERR) status = mpierr;
        }
        TRACE_COMM(MPI_Barrier)(ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Barrier");
    }
#endif
    return status;
}

/*----< ncmpio_sync() >------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpio_sync(void *ncdp)
{
    int err;
    NC *ncp = (NC*)ncdp;

    /* cannot be in define mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

#if 0
    /* In PnetCDF, header metadata is always sync-ed among all processes.
     * There is no need to re-read the header from file.
     */
    if (NC_readonly(ncp))
        /* calling sync for file opened for read only means re-read header */
        return ncmpio_read_NC(ncp);
#else
    if (NC_readonly(ncp)) return NC_NOERR;
#endif

    /* the only part of header that can be dirty is numrecs (caused only by
     * independent APIs) */
    if (ncp->vars.num_rec_vars > 0 && NC_indep(ncp)) {
        /* sync numrecs in memory among processes and in file */
        set_NC_ndirty(ncp);
        err = ncmpio_sync_numrecs(ncp);
        if (err != NC_NOERR) return err;
    }

    /* calling MPI_File_sync() on both collective and independent handlers */
    return ncmpio_file_sync(ncp);
}

