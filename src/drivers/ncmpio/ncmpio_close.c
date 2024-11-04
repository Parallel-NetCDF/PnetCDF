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
#include <sys/types.h> /* open(), lseek() */
#include <sys/stat.h>  /* open() */
#include <fcntl.h>     /* open() */
#include <unistd.h>    /* truncate(), lseek() */
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

    if (ncp->get_list      != NULL) NCI_Free(ncp->get_list);
    if (ncp->put_list      != NULL) NCI_Free(ncp->put_list);
    if (ncp->abuf          != NULL) NCI_Free(ncp->abuf);
    if (ncp->path          != NULL) NCI_Free(ncp->path);
    if (ncp->nonaggr_ranks != NULL) NCI_Free(ncp->nonaggr_ranks);

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

    if (ncp->nprocs > 1 && ncp->collective_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_close)(&ncp->collective_fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_File_close");
    }

    if (doUnlink) {
        /* called from ncmpi_abort, if the file is being created and is still
         * in define mode, the file is deleted */
        if (ncp->rank == 0) {
            TRACE_IO(MPI_File_delete)((char *)ncp->path, ncp->mpiinfo);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_error_mpi2nc(mpireturn, "MPI_File_delete");
        }
        if (ncp->nprocs > 1)
            MPI_Barrier(ncp->comm);
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
        printf("PnetCDF warning: %d nonblocking get requests still pending on process %d. Cancelling ...\n",ncp->numLeadGetReqs,ncp->rank);
        err = ncmpio_cancel(ncp, NC_GET_REQ_ALL, NULL, NULL);
        if (status == NC_NOERR) status = err;
        if (status == NC_NOERR) status = NC_EPENDING;
    }
    if (ncp->numLeadPutReqs > 0) {
        printf("PnetCDF warning: %d nonblocking put requests still pending on process %d. Cancelling ...\n",ncp->numLeadPutReqs,ncp->rank);
        err = ncmpio_cancel(ncp, NC_PUT_REQ_ALL, NULL, NULL);
        if (status == NC_NOERR) status = err;
        if (status == NC_NOERR) status = NC_EPENDING;
    }
#endif

    /* calling MPI_File_close() */
    err = ncmpio_close_files(ncp, 0);
    if (status == NC_NOERR) status = err;

    /* file is open for write and no variable has been defined */
    if (!NC_readonly(ncp) && ncp->vars.ndefined == 0) {
        /* wait until all processes close the file */
        if (ncp->nprocs > 1) MPI_Barrier(ncp->comm);

        if (ncp->rank == 0) {
            /* ignore all errors, as unexpected file size if not a fatal error */
#ifdef HAVE_TRUNCATE
            /* when calling POSIX I/O, remove file type prefix from file name */
            char *path = ncmpii_remove_file_system_type_prefix(ncp->path);
            int fd = open(path, O_RDWR, 0666);
            if (fd != -1) {
                /* obtain file size */
                off_t file_size = lseek(fd, 0, SEEK_END);
                /* truncate file size to header size, if larger than header */
                if (file_size > ncp->xsz && ftruncate(fd, ncp->xsz) < 0) {
                    err = ncmpii_error_posix2nc("ftruncate");
                    if (status == NC_NOERR) status = err;
                }
                close(fd);
            }
#else
            MPI_File fh;
            int mpireturn;
            mpireturn = MPI_File_open(MPI_COMM_SELF, ncp->path, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
            if (mpireturn == MPI_SUCCESS) {
                /* obtain file size */
                MPI_Offset *file_size;
                mpireturn = MPI_File_seek(fh, 0, MPI_SEEK_END);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_File_seek");
                    if (status == NC_NOERR) status = err;
                }
                mpireturn = MPI_File_get_position(fh, &file_size);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_File_get_position");
                    if (status == NC_NOERR) status = err;
                }
                /* truncate file size to header size, if larger than header */
                if (file_size > ncp->xsz) {
                    mpireturn = MPI_File_set_size(fh, ncp->xsz);
                    if (mpireturn != MPI_SUCCESS) {
                        err = ncmpii_error_mpi2nc(mpireturn,"MPI_File_set_size");
                        if (status == NC_NOERR) status = err;
                    }
                }
                mpireturn = MPI_File_close(&fh);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_File_close");
                    if (status == NC_NOERR) status = err;
                }
            }
            else {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_File_open");
                if (status == NC_NOERR) status = err;
            }
#endif
        }
        if (ncp->nprocs > 1) MPI_Barrier(ncp->comm);
    }

    /* free up space occupied by the header metadata */
    ncmpio_free_NC(ncp);

    return status;
}

