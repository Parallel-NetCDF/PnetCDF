/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_close() : dispatcher->close()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"
#ifdef ENABLE_SUBFILING
#include "ncmpio_subfile.h"
#endif

/*----< ncmpio_free_NC() >----------------------------------------------------*/
void
ncmpio_free_NC(NC *ncp)
{
    if (ncp == NULL) return;

    ncmpio_free_NC_dimarray(&ncp->dims);
    ncmpio_free_NC_attrarray(&ncp->attrs);
    ncmpio_free_NC_vararray(&ncp->vars);

    /* The only case that ncp->mpiinfo is MPI_INFO_NULL is when exiting endef
     * from a redef. All other cases reaching here are from ncmpi_close, in
     * which case ncp->mpiinfo is never MPI_INFO_NULL.
     */
    if (ncp->mpiinfo != MPI_INFO_NULL) MPI_Info_free(&ncp->mpiinfo);

    if (ncp->get_list != NULL) NCI_Free(ncp->get_list);
    if (ncp->put_list != NULL) NCI_Free(ncp->put_list);
    if (ncp->abuf     != NULL) NCI_Free(ncp->abuf);
    if (ncp->path     != NULL) NCI_Free(ncp->path);

    NCI_Free(ncp);
}

/*----< ncmpio_close_files() >-----------------------------------------------*/
int
ncmpio_close_files(NC *ncp, int doUnlink) {
    int mpireturn;

    assert(ncp != NULL); /* this should never occur */

    if (ncp->independent_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_close)(&ncp->independent_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_File_close");
    }

    if (ncp->collective_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_close)(&ncp->collective_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_File_close");
    }

    if (doUnlink) {
        /* called from ncmpi_abort, if the file is being created and is still
         * in define mode, the file is deleted */
        TRACE_IO(MPI_File_delete)((char *)ncp->path, ncp->mpiinfo);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_File_delete");
    }
    return NC_NOERR;
}

/*----< ncmpio_close() >------------------------------------------------------*/
/* This function is collective */
int
ncmpio_close(void *ncdp)
{
    int err=NC_NOERR, status=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (NC_indef(ncp)) { /* currently in define mode */
        status = ncmpio__enddef(ncp, 0, 0, 0, 0); /* TODO: defaults */

        if (status != NC_NOERR) {
            /* To do: Abort new definition, if any */
            if (ncp->old != NULL) {
                ncmpio_free_NC(ncp->old);
                ncp->old = NULL;
                fClr(ncp->flags, NC_MODE_DEF);
            }
        }
    }

    if (!NC_readonly(ncp) &&  /* file is open for write */
         NC_indep(ncp)) {     /* exit independent data mode will sync header */
        err = ncmpio_end_indep_data(ncp);
        if (status == NC_NOERR) status = err;
    }

    /* if entering this function in  collective data mode, we do not have to
     * update header in file, as file header is always up-to-date */

#ifdef ENABLE_SUBFILING
    /* ncmpio__enddef() will update ncp->num_subfiles */
    /* TODO: should check ncid_sf? */
    /* if the file has subfiles, close them first */
    if (ncp->num_subfiles > 1) {
        err = ncmpio_subfile_close(ncp);
        if (status == NC_NOERR) status = err;
    }
#endif

    /* We can cancel or complete all outstanding nonblocking I/O.
     * For now, cancelling makes more sense. */
#ifdef COMPLETE_NONBLOCKING_IO
    if (ncp->numLeadGetReqs > 0) {
        err = ncmpio_wait(ncp, NC_GET_REQ_ALL, NULL, NULL, NC_REQ_INDEP);
        if (status == NC_NOERR) status = err;
        if (status == NC_NOERR) status = NC_EPENDING;
    }
    if (ncp->numLeadPutReqs > 0) {
        err = ncmpio_wait(ncp, NC_PUT_REQ_ALL, NULL, NULL, NC_REQ_INDEP);
        if (status == NC_NOERR) status = err;
        if (status == NC_NOERR) status = NC_EPENDING;
    }
#else
    if (ncp->numLeadGetReqs > 0) {
        int rank;
        MPI_Comm_rank(ncp->comm, &rank);
        printf("PnetCDF warning: %d nonblocking get requests still pending on process %d. Cancelling ...\n",ncp->numLeadGetReqs,rank);
        err = ncmpio_cancel(ncp, NC_GET_REQ_ALL, NULL, NULL);
        if (status == NC_NOERR) status = err;
        if (status == NC_NOERR) status = NC_EPENDING;
    }
    if (ncp->numLeadPutReqs > 0) {
        int rank;
        MPI_Comm_rank(ncp->comm, &rank);
        printf("PnetCDF warning: %d nonblocking put requests still pending on process %d. Cancelling ...\n",ncp->numLeadPutReqs,rank);
        err = ncmpio_cancel(ncp, NC_PUT_REQ_ALL, NULL, NULL);
        if (status == NC_NOERR) status = err;
        if (status == NC_NOERR) status = NC_EPENDING;
    }
#endif

    /* If the user wants a stronger data consistency by setting NC_SHARE */
    if (NC_doFsync(ncp)) {
        err = ncmpio_file_sync(ncp); /* calling MPI_File_sync() */
        if (status == NC_NOERR) status = err;
    }

    /* calling MPI_File_close() */
    err = ncmpio_close_files(ncp, 0);
    if (status == NC_NOERR) status = err;

    /* free up space occupied by the header metadata */
    ncmpio_free_NC(ncp);

    return status;
}

