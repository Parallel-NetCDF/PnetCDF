/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in
 * src/dispatchers/var_getput.m4
 *
 * ncmpi_get_vard() : dispatcher->get_vard()
 * ncmpi_put_vard() : dispatcher->put_vard()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h> /* memcpy() */
#include <limits.h> /* INT_MAX */
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< getput_vard() >------------------------------------------------------*/
/* Starting from 1.9.1, it is possible to read/write multiple variables in a
 * single call to vard. Because it is difficult to tell the variable boundaries
 * from a single filetype, unless it is flattened, we enforce the following
 * requirements for vard APIs.
 * 1. The element data type of filetype must be the conform with the NC
 *    external data type of the variable, i.e. MPI_FLOAT <--> NC_FLOAT.
 *    If the requirement is violated, NC_ETYPE_MISMATCH will be returned.
 * 2. All variables accessed by filetype must be of the same data type.
 *    NC_EMULTITYPES will be returned if violated.
 * 3. buftype must contain only one element data type. Otherwise NC_EMULTITYPES
 *    will be returned.
 */
static int
getput_vard(NC               *ncp,
            NC_var           *varp,
            MPI_Datatype      filetype, /* access layout in the file */
            void             *buf,
            MPI_Offset        bufcount,
            MPI_Datatype      buftype,  /* data type of the buffer */
            int               reqMode)
{
    void *xbuf=NULL;
    int mpireturn, status=NC_NOERR, err=NC_NOERR, xtype_is_contig=1;
    int el_size, buftype_is_contig=0, need_swap_back_buf=0;
    int nelems=0, need_convert=0, need_swap=0, coll_indep, rw_flag;
    MPI_File fh;
    MPI_Offset fnelems=0, bnelems=0, offset=0;
    MPI_Datatype etype=MPI_DATATYPE_NULL, xtype=MPI_BYTE;
#if MPI_VERSION >= 3
    MPI_Count filetype_size=0;
    MPI_Count true_lb=0, true_ub=0, true_extent=0;
#else
    int filetype_size=0;
    MPI_Aint true_lb=0, true_ub=0, true_extent=0;
#endif

#ifdef ENABLE_SUBFILING
    /* call a separate routine if variable is stored in subfiles */
    if (varp->num_subfiles > 1) {
        printf("This feature for subfiling is yet to implement\n");
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
    }
#endif

    /* First check filetype whether it meets the requirement */
    if (filetype == MPI_DATATYPE_NULL) {
        /* This is actually an invalid filetype that can cause
         * MPI_File_set_view to fail. In PnetCDF, we simply consider this as
         * a zero-length request.
         */
        filetype_size = 0;
        goto err_check;
    }

    /* obtain size of filetype and its true upper bound.
     * filetype's lb may not always be 0 (e.g. created by constructor
     * MPI_Type_create_hindexed), we need to find the true last byte accessed
     * by this request, true_ub, in order to calculate new_numrecs.
     */
#if MPI_VERSION >= 3
    /* MPI_Type_size_x is introduced in MPI 3.0 */
    mpireturn = MPI_Type_size_x(filetype, &filetype_size);
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size_x");
        xtype = MPI_BYTE;
        goto err_check;
    }
    /* MPI_Type_get_true_extent_x is introduced in MPI 3.0 */
    MPI_Type_get_true_extent_x(filetype, &true_lb, &true_extent);
    true_ub = true_lb + true_extent;
#else
    /* PROBLEM: In MPI_Type_size(), argument filetype_size is a 4-byte integer,
     * cannot be used for large filetypes. Prior to MPI 3.0 standard, argument
     * "size" of MPI_Type_size is of type int. When int overflows, the returned
     * value in argument "size" may be a negative. */
    mpireturn = MPI_Type_size(filetype, &filetype_size);
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size");
        xtype = MPI_BYTE;
        goto err_check;
    }
    if (filetype_size < 0) { /* int overflow */
        DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
        goto err_check;
    }
    MPI_Type_get_true_extent(filetype, &true_lb, &true_extent);
    true_ub = true_lb + true_extent;
#endif

    if (filetype_size == 0) { /* zero-length request */
        goto err_check;
    }

#ifndef ENABLE_LARGE_SINGLE_REQ
    /* Not all MPI-IO libraries support single requests larger than 2 GiB */
    if (filetype_size > INT_MAX) {
        DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
        goto err_check;
    }
#endif

    /* get the corresponding MPI datatype of variable external type */
    xtype = ncmpii_nc2mpitype(varp->xtype);

    /* find the element type of filetype. ncmpii_dtype_decode() checks
     * NC_EMULTITYPES */
    err = ncmpii_dtype_decode(filetype, &etype, NULL, &fnelems, NULL, NULL);
    if (err != NC_NOERR) goto err_check;

    /* element type of filetype must be the same as variable's NC type */
    if (etype != xtype) {
        DEBUG_ASSIGN_ERROR(err, NC_ETYPE_MISMATCH)
        goto err_check;
    }

    /* done with checking filetype, now check buftype */

    if (bufcount == 0 && buftype != MPI_DATATYPE_NULL) {
        /* if this process has nothing to read/write */
        goto err_check;
    }

    if (buftype == MPI_DATATYPE_NULL) {
        /* In this case, the request size is the same as filetype */
        buftype = etype = xtype;
        MPI_Type_size(buftype, &el_size);
        bufcount = filetype_size / el_size;
        buftype_is_contig = 1;
        bnelems = bufcount;
    }
    else {
        /* find the element type of filetype. ncmpii_dtype_decode() checks
         * NC_EMULTITYPES */
        err = ncmpii_dtype_decode(buftype, &etype, &el_size, &bnelems, NULL,
                                  &buftype_is_contig);
        if (err != NC_NOERR) goto err_check;

        /* type conversion between non-char and char is not allowed */
        err = NCMPII_ECHAR(varp->xtype, etype);
        if (err != NC_NOERR) goto err_check;

        bnelems *= bufcount;
#ifndef ENABLE_LARGE_SINGLE_REQ
        if (bnelems != (int)bnelems) {
            DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
            goto err_check;
        }
#endif

        /* filetype's number of elements must be equal to request's */
        if (fnelems != bnelems) {
            DEBUG_ASSIGN_ERROR(err, NC_EIOMISMATCH)
            goto err_check;
        }
    }
    xtype_is_contig = buftype_is_contig;

    /* check if type conversion and Endianness byte swap is needed */
    need_convert = ncmpii_need_convert(ncp->format, varp->xtype, etype);
    need_swap    = NEED_BYTE_SWAP(varp->xtype, etype);

    if (fIsSet(reqMode, NC_REQ_WR)) {
        int in_place_swap = 0;
        if (need_swap) {
            if (fIsSet(ncp->flags, NC_MODE_SWAP_ON))
                in_place_swap = 1;
            else if (! fIsSet(ncp->flags, NC_MODE_SWAP_OFF)) { /* auto mode */
                if (filetype_size > NC_BYTE_SWAP_BUFFER_SIZE)
                    in_place_swap = 1;
            }
        }
        if (!need_convert &&
            (!need_swap || (in_place_swap && buftype_is_contig))) {
            /* reuse buftype, bufcount, buf in later MPI file write */
            xbuf = buf;
            if (need_swap) {
                ncmpii_in_swapn(xbuf, bnelems, varp->xsz);
                need_swap_back_buf = 1;
            }
        }
        else { /* must allocate a buffer to convert/swap/pack */
            xbuf = NCI_Malloc((size_t)filetype_size);
            if (xbuf == NULL) {
                DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
                goto err_check;
            }
            need_swap_back_buf = 0;
            xtype_is_contig = 1;

            /* pack buf to xbuf, byte-swap and type-convert on xbuf, which
             * will later be used in MPI file write */
            err = ncmpio_pack_xbuf(ncp->format, varp, bufcount, buftype,
                                   buftype_is_contig, bnelems, etype, el_size,
                                   MPI_DATATYPE_NULL, need_convert, need_swap,
                                   filetype_size, buf, xbuf);
            if (err != NC_NOERR && err != NC_ERANGE) {
                if (xbuf != buf) NCI_Free(xbuf);
                xbuf = NULL;
                goto err_check;
            }
        }
    }
    else { /* read request */
        if (!need_convert && (!need_swap || buftype_is_contig)) {
            /* reuse buftype, bufcount, buf in later MPI file read */
            xbuf = buf;
        }
        else { /* must allocate a buffer to convert/swap/unpack */
            xbuf = NCI_Malloc((size_t)filetype_size);
            if (xbuf == NULL) {
                DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
                goto err_check;
            }
            xtype_is_contig = 1;
        }
    }

    /* Set nelems and xtype which will be used in MPI read/write */
    if (buf != xbuf) {
        /* xbuf is a malloc-ed contiguous buffer */
        nelems = (int)bnelems;
    }
    else {
        /* we can safely use bufcount and buftype in MPI File read/write.
         * Note buftype may be noncontiguous. */
        nelems = bufcount;
        xtype = buftype;
    }

    /* set fileview's displacement to the variable's starting file offset */
    offset = varp->begin;

err_check:
    /* check error before going into a collective call.
     * If an error has been detected on one or more processes, we'll still
     * conduct a zero-byte operation (everyone has to participate in the
     * collective I/O call) but return the error at the end. NC_ERANGE is not a
     * fatal error, we proceed with write request.
     */
    if ((err != NC_NOERR && err != NC_ERANGE) || bufcount == 0 ||
        filetype_size == 0) {
        if (fIsSet(reqMode, NC_REQ_INDEP)) {
            if (need_swap_back_buf)
                /* byte-swap back to buf's original contents */
                ncmpii_in_swapn(buf, bnelems, varp->xsz);
            if (xbuf != NULL && xbuf != buf) NCI_Free(xbuf);
            return err;
        }
        /* for NC_REQ_COLL, this process must participate successive collective
         * MPI-IO calls as a zero-length request.
         */
        offset   = 0;
        bufcount = 0;
        nelems   = 0;
        filetype_size = 0;
        filetype = MPI_BYTE;
        buftype  = MPI_BYTE;
        xtype    = MPI_BYTE;
        etype    = MPI_BYTE;
    }
    status = err;

    if (fIsSet(reqMode, NC_REQ_COLL)) {
        fh = ncp->collective_fh;
        coll_indep = NC_REQ_COLL;
    } else {
        fh = ncp->independent_fh;
        coll_indep = NC_REQ_INDEP;
    }

    /* set the MPI-IO fileview, this is a collective call */
    err = ncmpio_file_set_view(ncp, fh, &offset, filetype);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        nelems = 0; /* skip this request */
    }

    rw_flag = (fIsSet(reqMode, NC_REQ_RD)) ? NC_REQ_RD : NC_REQ_WR;

    err = ncmpio_read_write(ncp, rw_flag, coll_indep, offset, nelems,
                            xtype, xbuf, xtype_is_contig);
    if (status == NC_NOERR) status = err;

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
     */

    if (fIsSet(reqMode, NC_REQ_RD)) {
        if (filetype_size == 0) return status;

        /* swap/convert/unpack xbuf into user buffer, buf, when necessary */
        err = ncmpio_unpack_xbuf(ncp->format, varp, bufcount, buftype,
                                 buftype_is_contig, bnelems, etype,
                                 MPI_DATATYPE_NULL, need_convert, need_swap,
                                 buf, xbuf);
        if (status == NC_NOERR) status = err;
    }
    else { /* write request */
        if (need_swap_back_buf) /* byte-swap back to buf's original contents */
            ncmpii_in_swapn(buf, bnelems, varp->xsz);

        if (IS_RECVAR(varp)) {
            /* update header's number of records */
            MPI_Offset new_numrecs = true_ub / ncp->recsize;
            if (true_ub % ncp->recsize) new_numrecs++;

            if (fIsSet(reqMode, NC_REQ_INDEP)) {
                /* For independent put, we delay the sync for numrecs until
                 * the next collective call, such as end_indep(), sync(),
                 * enddef(), or close(). This is because if we update numrecs
                 * to file now, race condition can happen. Note numrecs in
                 * memory may be inconsistent and obsolete till then.
                 */
                if (ncp->numrecs < new_numrecs) {
                    ncp->numrecs = new_numrecs;
                    set_NC_ndirty(ncp);
                }
            }
            else { /* NC_REQ_COLL: sync numrecs in memory and file */
                /* new_numrecs may be different among processes.
                 * First, find the max numrecs among all processes.
                 */
                MPI_Offset max_numrecs;
                TRACE_COMM(MPI_Allreduce)(&new_numrecs, &max_numrecs, 1,
                                          MPI_OFFSET, MPI_MAX, ncp->comm);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
                    if (status == NC_NOERR) status = err;
                }
                /* In collective mode, ncp->numrecs is always sync-ed among
                   processes */
                if (ncp->numrecs < max_numrecs) {
                    err = ncmpio_write_numrecs(ncp, max_numrecs);
                    if (status == NC_NOERR) status = err;
                    ncp->numrecs = max_numrecs;
                }
            }
        }

        if (NC_doFsync(ncp)) { /* NC_SHARE is set */
            TRACE_IO(MPI_File_sync)(fh);
            if (fIsSet(reqMode, NC_REQ_COLL))
                TRACE_COMM(MPI_Barrier)(ncp->comm);
        }
    }

    if (xbuf != NULL && xbuf != buf) NCI_Free(xbuf);

    return status;
}

/*----< ncmpio_get_vard() >--------------------------------------------------*/
int
ncmpio_get_vard(void         *ncdp,
                int           varid,
                MPI_Datatype  filetype,  /* access layout in file */
                void         *buf,
                MPI_Offset    bufcount,
                MPI_Datatype  buftype,   /* data type of the buffer */
                int           reqMode)
{
    NC *ncp=(NC*)ncdp;

    if (fIsSet(reqMode, NC_REQ_ZERO) && fIsSet(reqMode, NC_REQ_COLL))
        /* this collective API has a zero-length request */
        return ncmpio_getput_zero_req(ncp, reqMode);

    /* Note sanity check for ncdp and varid has been done in dispatchers */

    return getput_vard(ncp, ncp->vars.value[varid], filetype, buf, bufcount,
                       buftype, reqMode);
}

/*----< ncmpio_put_vard() >--------------------------------------------------*/
int
ncmpio_put_vard(void         *ncdp,
                int           varid,
                MPI_Datatype  filetype, /* access layout in the file */
                const void   *buf,
                MPI_Offset    bufcount,
                MPI_Datatype  buftype,   /* data type of the buffer */
                int           reqMode)
{
    NC *ncp=(NC*)ncdp;

    if (fIsSet(reqMode, NC_REQ_ZERO) && fIsSet(reqMode, NC_REQ_COLL))
        /* this collective API has a zero-length request */
        return ncmpio_getput_zero_req(ncp, reqMode);

    /* Note sanity check for ncdp and varid has been done in dispatchers */

    return getput_vard(ncp, ncp->vars.value[varid], filetype, (void*)buf,
                       bufcount, buftype, reqMode);
}
