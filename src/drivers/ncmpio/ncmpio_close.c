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
    if (ncp->nonaggr_ranks != NULL) NCI_Free(ncp->nonaggr_ranks);

    NCI_Free(ncp);
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

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    int i, j, ntimers;
    double tt[16], max_t[16], put_time=0, get_time=0;
    MPI_Offset sizes[16], max_sizes[16], max_npairs_put=0, max_npairs_get=0;

    /* print intra-node aggregation timing breakdown */
    if (ncp->num_aggrs_per_node > 0) {
        j = 0;
        for (i=0; i<6; i++) sizes[j++] = ncp->maxmem_put[i];
        for (i=0; i<6; i++) sizes[j++] = ncp->maxmem_get[i];
        sizes[12] = ncp->ina_npairs_put;
        sizes[13] = ncp->ina_npairs_get;

        MPI_Allreduce(sizes, max_sizes, 14, MPI_OFFSET, MPI_MAX, ncp->comm);
        max_npairs_put = max_sizes[12];
        max_npairs_get = max_sizes[13];

        for (i=0; i<12; i++) tt[i] = (float)(max_sizes[i]) / 1048576.0; /* in MiB */
        if (ncp->rank == 0 && max_npairs_put > 0)
            printf("%s: INA put npairs=%lld mem=%.1f %.1f %.1f %.1f %.1f %.1f (MiB)\n",
                   __func__, max_sizes[12], tt[0],tt[1],tt[2],tt[3],tt[4],tt[5]);
        if (ncp->rank == 0 && max_npairs_get > 0)
            printf("%s: INA get npairs=%lld mem=%.1f %.1f %.1f %.1f %.1f %.1f (MiB)\n",
                   __func__, max_sizes[13], tt[6],tt[7],tt[8],tt[9],tt[10],tt[11]);

        if (max_npairs_put > 0) { /* put npairs > 0 */
            put_time = ncp->ina_time_init + ncp->ina_time_flatten;
            ntimers = 4;
            for (i=0; i<ntimers; i++) {
                tt[i]     = ncp->ina_time_put[i];
                put_time += tt[i];
            }
            tt[ntimers]   = ncp->ina_time_init;
            tt[ntimers+1] = ncp->ina_time_flatten;
            tt[ntimers+2] = put_time;

            MPI_Reduce(tt, max_t, ntimers+3, MPI_DOUBLE, MPI_MAX, 0, ncp->comm);
            put_time = max_t[ntimers+2];
            if (ncp->rank == 0)
                printf("%s: INA put timing %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f = %5.2f\n",
                __func__, max_t[ntimers],max_t[ntimers+1],max_t[0],max_t[1],max_t[2],max_t[3],put_time);
        }
        if (max_npairs_get > 0) { /* get npairs > 0 */
            get_time = ncp->ina_time_init + ncp->ina_time_flatten;
            ntimers = 4;
            for (i=0; i<ntimers; i++) {
                tt[i]     = ncp->ina_time_get[i];
                get_time += tt[i];
            }
            tt[ntimers]   = ncp->ina_time_init;
            tt[ntimers+1] = ncp->ina_time_flatten;
            tt[ntimers+2] = get_time;

            MPI_Reduce(tt, max_t, ntimers+3, MPI_DOUBLE, MPI_MAX, 0, ncp->comm);
            if (ncp->rank == 0)
                printf("%s: INA get timing %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f = %5.2f\n",
                __func__, max_t[ntimers],max_t[ntimers+1],max_t[0],max_t[1],max_t[2],max_t[3],max_t[ntimers+2]);
        }
    }
#endif

    /* close the file */
    err = ncmpio_file_close(ncp);
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
            TRACE_IO(MPI_File_open, (MPI_COMM_SELF, ncp->path, MPI_MODE_RDWR, MPI_INFO_NULL, &fh));
            if (mpireturn == MPI_SUCCESS) {
                /* obtain file size */
                MPI_Offset *file_size;
                TRACE_IO(MPI_File_get_size, (fh, &file_size));
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    if (status == NC_NOERR) status = err;
                }
                /* truncate file size to header size, if larger than header */
                if (file_size > ncp->xsz) {
                    TRACE_IO(MPI_File_set_size, (fh, ncp->xsz));
                    if (mpireturn != MPI_SUCCESS) {
                        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                        if (status == NC_NOERR) status = err;
                    }
                }
                TRACE_IO(MPI_File_close, (&fh));
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    if (status == NC_NOERR) status = err;
                }
            }
            else {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                if (status == NC_NOERR) status = err;
            }
#endif
        }
        if (ncp->nprocs > 1) MPI_Barrier(ncp->comm);
    }

    /* free the intra-node aggregation communicator */
    if (ncp->ina_comm != MPI_COMM_NULL)
        MPI_Comm_free(&ncp->ina_comm);

    /* collectively return the same error code */
    if (ncp->nprocs > 1)
        MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_MIN, ncp->comm);

    /* free up space occupied by the header metadata */
    ncmpio_free_NC(ncp);

    return status;
}

