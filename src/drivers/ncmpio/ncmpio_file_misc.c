/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_close()            : dispatcher->close()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_inq_xxx()          : dispatcher->inq_misc()
 * ncmpi_flush()            : dispatcher->flush()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy() */
#include <assert.h>
#include <errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< dup_NC() >-----------------------------------------------------------*/
static NC *
dup_NC(const NC *ref)
{
    NC *ncp;

    ncp = (NC *) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) return NULL;

    /* copy most of the NC members over */
    *ncp = *ref;

#ifndef SEARCH_NAME_LINEARLY
    /* set hash tables to NULL, indicating space not-yet allocated */
    ncp->dims.nameT = NULL;
    ncp->vars.nameT = NULL;
    ncp->attrs.nameT = NULL;
#endif

    if (ncmpio_dup_NC_dimarray(&ncp->dims,   &ref->dims)  != NC_NOERR ||
        ncmpio_dup_NC_attrarray(&ncp->attrs, &ref->attrs) != NC_NOERR ||
        ncmpio_dup_NC_vararray(&ncp->vars,   &ref->vars, ref->hash_size_attr) != NC_NOERR) {
        ncmpio_free_NC(ncp);
        return NULL;
    }

    if (ref->nonaggr_ranks != NULL) {
        size_t len = sizeof(int) * ncp->num_nonaggrs;
        ncp->nonaggr_ranks = (int*) NCI_Malloc(len);
        memcpy(ncp->nonaggr_ranks, ref->nonaggr_ranks, len);
    }

    /* fields below should not copied from ref */
    ncp->comm       = MPI_COMM_NULL;
    ncp->mpiinfo    = MPI_INFO_NULL;
    ncp->get_list   = NULL;
    ncp->put_list   = NULL;
    ncp->abuf       = NULL;
    ncp->path       = NULL;

    return ncp;
}

/*----< ncmpio_redef() >-----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpio_redef(void *ncdp)
{
    int err, status=NC_NOERR;
    NC *ncp = (NC*)ncdp;

#if 0
    if (NC_readonly(ncp)) DEBUG_RETURN_ERROR(NC_EPERM) /* read-only */
    /* if open mode is inconsistent, then this return might cause parallel
     * program to hang */

    /* cannot be in define mode, must enter from data mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* sync all metadata, including numrecs, if changed in independent mode.
     * also ensure exiting define mode always entering collective data mode
     */
#endif
    if (NC_indep(ncp)) /* exit independent mode, if in independent mode */
        ncmpio_end_indep_data(ncp);

    /* duplicate header to be used in enddef() for checking if header grows */
    ncp->old = dup_NC(ncp);
    if (ncp->old == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* we are now entering define mode */
    fSet(ncp->flags, NC_MODE_DEF);

    /* must reset fileview as header extent may later change in enddef() */
    err = ncmpio_file_set_view(ncp, 0, MPI_BYTE, 0, NULL, NULL);
    DEBUG_ASSIGN_ERROR(status, err)

    return status;
}

/*----< ncmpio_begin_indep_data() >------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpio_begin_indep_data(void *ncdp)
{
    NC *ncp = (NC*)ncdp;

    if (NC_indef(ncp))  /* must not be in define mode */
        DEBUG_RETURN_ERROR(NC_EINDEFINE)

    if (NC_indep(ncp))  /* already in indep data mode */
        return NC_NOERR;
        /* starting from 1.2.0, calling begin_indep_data() in independent data
         * mode is no longer considered illegal
        DEBUG_RETURN_ERROR(NC_EINDEP)
        */

    /* If users want a stronger data consistency, ncmpi_sync() should be called
     * following this subroutine. */

    /* raise independent flag */
    fSet(ncp->flags, NC_MODE_INDEP);

    /* Barrier is necessary to prevent non-aggregators from calling open()
     * before the file is being collectively created by the aggregators.
     */
    MPI_Barrier(ncp->comm);

    if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
        /* When using PnetCDF's PNCIO driver, there are 2 scenarios:
         * 1. When intra-node aggregation (INA) is enabled, at the end of
         *    ncmpi_create/ncmpi_open, non-aggregators' pncio_fh are NULL. Thus
         *    switching to independent data mode, we can re-use pncio_fh to
         *    store file handler of file opened with MPI_COMM_SELF. Note
         *    whether pncio_fh is NULL or not does not tell whether INA is
         *    enabled or not.
         * 2. When INA is disabled, all ranks calls PNCIO_File_open() and thus
         *    pncio_fh should not be NULL. In other word, this scenario should
         *    not reach here at all. Because PnetCDF's PNCIO driver relaxes
         *    File_setview subroutine to be able to called independently, the
         *    same pncio_fh can be used for both collective and independent I/O
         *    APIs. Note we cannot re-used pncio_fh for the above scenario 1,
         *    because in the collective data mode, all ranks must participate
         *    each collective I/O call,
         */
        int err;
        char *filename;

        if (ncp->pncio_fh != NULL)
            /* Only INA non-aggregators' pncio_fh can be NULL, because
             * aggregators open the file collectively and their pncio_fh can
             * never be NULL.
             */
            return NC_NOERR;

        filename = ncmpii_remove_file_system_type_prefix(ncp->path);

        ncp->pncio_fh = (PNCIO_File*) NCI_Calloc(1,sizeof(PNCIO_File));
        ncp->pncio_fh->file_system = ncp->fstype;
        ncp->pncio_fh->node_ids.num_nodes = 1;
        ncp->pncio_fh->node_ids.ids = (int*) NCI_Malloc(sizeof(int));
        ncp->pncio_fh->node_ids.ids[0] = 0;

        int omode = fClr(ncp->mpiomode, MPI_MODE_CREATE);

        err = PNCIO_File_open(MPI_COMM_SELF, filename, omode, ncp->mpiinfo,
                              ncp->pncio_fh);
        if (err != NC_NOERR)
            return err;

        /* Get the I/O hints used/modified by MPI-IO. Note ncp->mpiinfo may
         * have been populated. It can be discarded and replaced by the one
         * used by MPI-IO.
         */
        if (ncp->mpiinfo != MPI_INFO_NULL)
            MPI_Info_free(&ncp->mpiinfo);

        err = PNCIO_File_get_info(ncp->pncio_fh, &ncp->mpiinfo);
        if (err != NC_NOERR) return err;

        /* Add PnetCDF hints into ncp->mpiinfo */
        ncmpio_hint_set(ncp, ncp->mpiinfo);

        NCI_Free(ncp->pncio_fh->node_ids.ids);
        ncp->pncio_fh->node_ids.num_nodes = 0;
        ncp->pncio_fh->node_ids.ids = NULL;

        return NC_NOERR;
    }

    /* PnetCDF's default mode is collective. MPI file handle, collective_fh,
     * will never be MPI_FILE_NULL. We must use a separate MPI file handle
     * opened with MPI_COMM_SELF, because MPI_File_set_view is a collective
     * call and accessing a subarray requires a call to MPI_File_set_view.
     * In independent data mode, no collective MPI operation can be implicitly
     * called.
     */
    if (ncp->independent_fh == MPI_FILE_NULL) {
        char *mpi_name;
        int mpireturn;
        TRACE_IO(MPI_File_open, (MPI_COMM_SELF, ncp->path, ncp->mpiomode,
                                 ncp->mpiinfo, &ncp->independent_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);

        /* for those ranks whose mpiinfo is NULL, retrieve info */
        if (ncp->mpiinfo == MPI_INFO_NULL) {
            /* get the I/O hints used/modified by MPI-IO */
            mpireturn = MPI_File_get_info(ncp->independent_fh, &ncp->mpiinfo);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_error_mpi2nc(mpireturn, mpi_name);

            /* Copy MPI-IO hints into ncp->mpiinfo */
            ncmpio_hint_set(ncp, ncp->mpiinfo);
        }
    }
    return NC_NOERR;
}

/*----< ncmpio_end_indep_data() >--------------------------------------------*/
/* This is a collective subroutine.
 * It can be called from:
 *    1. ncmpi_end_indep_data()
 *    2. ncmpi_redef() from independent data mode entering to define more
 *    3. ncmpi_close() when closing the file
 */
int
ncmpio_end_indep_data(void *ncdp)
{
    int status=NC_NOERR;
    NC *ncp = (NC*)ncdp;

    if (NC_indef(ncp))  /* must not be in define mode */
        DEBUG_RETURN_ERROR(NC_EINDEFINE)

    if (!NC_indep(ncp)) /* already in collective ata mode */
        return NC_NOERR;
        /* starting from 1.9.0, calling end_indep_data() in collective data
         * mode is no longer considered illegal
        DEBUG_RETURN_ERROR(NC_ENOTINDEP)
        */

    if (!NC_readonly(ncp)) {
        if (ncp->vars.num_rec_vars > 0) {
            /* numrecs dirty bit may not be the same across all processes.
             * force sync in memory no matter if dirty or not.
             */
            set_NC_ndirty(ncp);
            status = ncmpio_sync_numrecs(ncp);
            /* the only possible dirty part of the header is numrecs */
        }
    }

    fClr(ncp->flags, NC_MODE_INDEP);

    return status;
}

/*----< ncmpio_abort() >-----------------------------------------------------*/
/* This API is a collective subroutine */
int
ncmpio_abort(void *ncdp)
{
   /*
    * In data mode, same as ncmpi_close().
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
        assert(fIsSet(ncp->flags, NC_MODE_DEF));
        ncmpio_free_NC(ncp->old);
        ncp->old = NULL;
        fClr(ncp->flags, NC_MODE_DEF);
    }

    if (!doUnlink) {
        if (!NC_readonly(ncp) &&  /* file is open for write */
             NC_indep(ncp)) {     /* in independent data mode */
            /* exit independent mode, if in independent mode */
            status = ncmpio_end_indep_data(ncp); /* will sync header */
        }
    }

    /* close the file */
    err = ncmpio_file_close(ncp);
    if (status == NC_NOERR ) status = err;

    if (doUnlink) {
        err = ncmpio_file_delete(ncp);
        status = (status == NC_NOERR) ? err : status;
    }

    /* free up space occupied by the header metadata */
    ncmpio_free_NC(ncp);

    return status;
}

/*----< ncmpio_inq() >-------------------------------------------------------*/
int
ncmpio_inq(void *ncdp,
           int  *ndimsp,
           int  *nvarsp,
           int  *nattsp,
           int  *xtendimp)
{
    NC *ncp = (NC*)ncdp;

    if (ndimsp   != NULL) *ndimsp   = ncp->dims.ndefined;
    if (nvarsp   != NULL) *nvarsp   = ncp->vars.ndefined;
    if (nattsp   != NULL) *nattsp   = ncp->attrs.ndefined;
    if (xtendimp != NULL) *xtendimp = ncp->dims.unlimited_id;

    return NC_NOERR;
}

/*----< ncmpio_inq_misc() >--------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpio_inq_misc(void       *ncdp,
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
    if (pathlen != NULL) {
        if (ncp->path == NULL) *pathlen = 0;
        else                   *pathlen = (int)strlen(ncp->path);
    }
    if (path != NULL) {
        if (ncp->path == NULL) *path = '\0';
        else                   strcpy(path, ncp->path);
    }

    /* obtain the number of fixed-size variables */
    if (num_fix_varsp != NULL) {
        int num_rec_vars;
        if (NC_indef(ncp)) {
            /* if in define mode, recalculate the number of record variables */
            num_rec_vars = 0;
            for (i=0; i<ncp->vars.ndefined; i++)
                num_rec_vars += IS_RECVAR(ncp->vars.value[i]);
        }
        else
            num_rec_vars = ncp->vars.num_rec_vars;

        /* no. fixed-size == ndefined - no. record variables */
        *num_fix_varsp = ncp->vars.ndefined - num_rec_vars;
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
        MPI_Info_get(ncp->mpiinfo, "striping_unit", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        *striping_size = 0;
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            *striping_size = (int)strtol(value,NULL,10);
            if (errno != 0) *striping_size = 0;
        }
    }

    if (striping_count != NULL) {
        MPI_Info_get(ncp->mpiinfo, "striping_factor", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        *striping_count = 0;
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            *striping_count = (int)strtol(value,NULL,10);
            if (errno != 0) *striping_count = 0;
        }
    }

    /* the amount of writes, in bytes, committed to file system so far */
    if (put_size != NULL) *put_size = ncp->put_size;

    /* the amount of reads, in bytes, obtained from file system so far */
    if (get_size != NULL) *get_size = ncp->get_size;

    if (recsize != NULL) *recsize = ncp->recsize;

    if (header_size != NULL) {
        if (NC_indef(ncp))
            /* When called in define mode, calculate and return the current
             * header size. Cannot do the same for header extent, as the empty
             * space depends on arguments h_minfree and v_align of
             * ncmpi__enddef().
             */
            *header_size = ncmpio_hdr_len_NC(ncp);
        else
            *header_size = ncp->xsz;
    }

    if (header_extent != NULL) *header_extent = ncp->begin_var;

    if (info_used != NULL) {
        mpireturn = MPI_Info_dup(ncp->mpiinfo, info_used);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Info_dup");

        /* PnetCDF hints have been added to ncp->mpiinfo at ncmpi_enddef.
         *
         * Note MPI implementations may choose to ignore unrecognized hints and
         * MPI_File_get_info() may returns no PnetCDF hints. We need to add the
         * PnbetCDF hints explicitly to the info object brefore returning it to
         * user.
         */

        sprintf(value, OFFFMT, ncp->v_align);
        MPI_Info_set(*info_used, "nc_var_align_size", value);

        sprintf(value, OFFFMT, ncp->r_align);
        MPI_Info_set(*info_used, "nc_record_align_size", value);

        sprintf(value, "%d", ncp->hdr_chunk);
        MPI_Info_set(*info_used, "nc_header_read_chunk_size", value);

        sprintf(value, "%d", ncp->data_chunk);
        MPI_Info_set(*info_used, "pnc_data_move_chunk_size", value);

        if (fIsSet(ncp->flags, NC_MODE_SWAP_ON))
            MPI_Info_set(*info_used, "nc_in_place_swap", "enable");
        else if (fIsSet(ncp->flags, NC_MODE_SWAP_OFF))
            MPI_Info_set(*info_used, "nc_in_place_swap", "disable");
        else
            MPI_Info_set(*info_used, "nc_in_place_swap", "auto");

        sprintf(value, OFFFMT, ncp->ibuf_size);
        MPI_Info_set(*info_used, "nc_ibuf_size", value);

#ifdef ENABLE_SUBFILING
        if (ncp->subfile_mode)
            MPI_Info_set(*info_used, "pnetcdf_subfiling", "enable");
        else
            MPI_Info_set(*info_used, "pnetcdf_subfiling", "disable");
        sprintf(value, "%d", ncp->num_subfiles);
        MPI_Info_set(*info_used, "nc_num_subfiles", value);
#else
        MPI_Info_set(*info_used, "pnetcdf_subfiling", "disable");
#endif
    }

    if (nreqs != NULL)
        *nreqs = ncp->numLeadGetReqs + ncp->numLeadPutReqs;

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

/*----< ncmpi_delete() >-----------------------------------------------------*/
/* doesn't do anything to release resources. Users are advised to call
 * ncmpi_close() before calling this function.
 *
 * filename: the name of the file we will remove.
 * info: MPI info object, in case underlying file system needs hints.
 */
int
ncmpi_delete(const char *filename,
             MPI_Info    info)
{
    int err = NC_NOERR;
#ifdef MIMIC_LUSTRE
    char *path = ncmpii_remove_file_system_type_prefix(filename);
    err = unlink(path);
    if (err != 0)
        err = ncmpii_error_posix2nc("unlink");
#else
    err = PNCIO_File_delete(filename);
#if 0
    char *mpi_name;
    int mpireturn;

    TRACE_IO(MPI_File_delete, (filename, info));
    if (mpireturn != MPI_SUCCESS)
        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
#endif
#endif
    return err;
}


/*----< ncmpio_flush() >------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode
 */
int
ncmpio_flush(void *ncdp)
{
    /* Flush has no effect in ncmpio */
    return NC_NOERR;
}

