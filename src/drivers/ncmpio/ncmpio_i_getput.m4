dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in
 * src/dispatchers/var_getput.m4
 *
 * ncmpi_iget_var<kind>()        : dispatcher->iget_var()
 * ncmpi_iput_var<kind>()        : dispatcher->iput_var()
 * ncmpi_iget_var<kind>_<type>() : dispatcher->iget_var()
 * ncmpi_iput_var<kind>_<type>() : dispatcher->iput_var()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <limits.h> /* INT_MAX */
#include <assert.h>

#include <string.h> /* memcpy() */
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< abuf_malloc() >------------------------------------------------------*/
/* allocate memory space from the attached buffer pool */
static int
abuf_malloc(NC *ncp, MPI_Offset nbytes, void **buf, int *abuf_index)
{
    /* extend the table size if more entries are needed */
    if (ncp->abuf->tail + 1 == ncp->abuf->table_size) {
        ncp->abuf->table_size += NC_ABUF_DEFAULT_TABLE_SIZE;
        ncp->abuf->occupy_table = (NC_buf_status*)
                   NCI_Realloc(ncp->abuf->occupy_table,
                   (size_t)ncp->abuf->table_size * sizeof(NC_buf_status));
    }
    /* mark the new entry is used and store the requested buffer size */
    ncp->abuf->occupy_table[ncp->abuf->tail].is_used  = 1;
    ncp->abuf->occupy_table[ncp->abuf->tail].req_size = nbytes;
    *abuf_index = ncp->abuf->tail;

    *buf = (char*)ncp->abuf->buf + ncp->abuf->size_used;
    ncp->abuf->size_used += nbytes;
    ncp->abuf->tail++;

    return NC_NOERR;
}

/*----< abuf_dealloc() >-----------------------------------------------------*/
/* deallocate (actually un-register) memory space from the attached buffer
 * pool
 */
static int
abuf_dealloc(NC *ncp, int abuf_index)
{
    assert(abuf_index == ncp->abuf->tail - 1);

    /* mark the tail entry un-used */
    ncp->abuf->size_used -= ncp->abuf->occupy_table[abuf_index].req_size;
    ncp->abuf->occupy_table[abuf_index].req_size = 0;
    ncp->abuf->occupy_table[abuf_index].is_used  = 0;
    ncp->abuf->tail--;

    return NC_NOERR;
}

/*----< add_record_requests() >----------------------------------------------*/
/* check if this is a record variable. if yes, add a new request for each
 * record into the list. Hereinafter, treat each request as a non-record
 * variable request
 */
static int
add_record_requests(NC_req           *reqs,
                    const MPI_Offset *stride)
{
    int    i, j, ndims;
    size_t dims_chunk;
    MPI_Offset rec_bufsize, *req0_start, *req0_count, *req0_stride;

    ndims = reqs[0].varp->ndims;
    if (stride != NULL)
        dims_chunk = (size_t)ndims * 3 * SIZEOF_MPI_OFFSET;
    else
        dims_chunk = (size_t)ndims * 2 * SIZEOF_MPI_OFFSET;

    /* pointers to start[], count[] and stride[] of reqs[0] */
    req0_start  = reqs[0].start;
    req0_count  = req0_start + ndims;
    req0_stride = (stride == NULL) ? NULL : req0_count + ndims;

    rec_bufsize = reqs[0].varp->xsz;
    for (i=1; i<ndims; i++) rec_bufsize *= req0_count[i];

    /* append each record to the end of list */
    for (i=1; i<req0_count[0]; i++) {
        MPI_Offset *reqi_start, *reqi_count, *reqi_stride; /* of reqs[i] */

        reqs[i] = reqs[0]; /* inherit most attributes from reqs[0], except
                            * below ones, including the ones need malloc
                            */

        reqs[i].start = (MPI_Offset*) NCI_Malloc(dims_chunk);
        reqi_start = reqs[i].start;
        reqi_count = reqi_start + ndims;

        if (stride != NULL) {
            reqi_stride    = reqi_count + ndims;
            reqi_start[0]  = req0_start[0] + stride[0] * i;
            reqi_stride[0] = req0_stride[0]; /* this does not matter */
        }
        else {
            reqi_stride   = NULL;
            reqi_start[0] = req0_start[0] + i;
            fSet(reqs[i].flag, NC_REQ_STRIDE_NULL);
        }

        reqi_count[0] = 1;  /* one record in each sub request */
        for (j=1; j<ndims; j++) {
            reqi_start[j]  = req0_start[j];
            reqi_count[j]  = req0_count[j];
            if (stride != NULL) reqi_stride[j] = req0_stride[j];
        }

        /* xbuf cannot be NULL    assert(reqs[0].xbuf != NULL); */

        reqs[i].buf      = (char*)(reqs[i-1].buf)  + rec_bufsize;
        reqs[i].xbuf     = (char*)(reqs[i-1].xbuf) + rec_bufsize;

        /* reqs[i].bufcount and reqs[i].buftype will not be used in
         * wait call, only the lead request's matters */
    }

    return NC_NOERR;
}

/*----< ncmpio_igetput_varm() >-----------------------------------------------*/
int
ncmpio_igetput_varm(NC               *ncp,
                    NC_var           *varp,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    const MPI_Offset  imap[],
                    void             *buf,      /* user buffer */
                    MPI_Offset        bufcount,
                    MPI_Datatype      buftype,
                    int              *reqid,    /* out, can be NULL */
                    int               reqMode,
                    int               isSameGroup) /* if part of a varn group */
{
    void *xbuf=NULL;
    int i, j, err=NC_NOERR, abuf_index=-1, el_size, buftype_is_contig;
    int need_convert, need_swap, need_swap_back_buf=0, free_xbuf=0;
    MPI_Offset bnelems=0, nbytes;
    MPI_Datatype ptype, imaptype;
    NC_req *req;

    /* decode buftype to obtain the followings:
     * ptype:    element data type (MPI primitive type) in buftype
     * bufcount: If it is -1, then this is called from a high-level API and in
     *           this case buftype will be an MPI primitive data type.
     *           If it is >=0, then this is called from a flexible API.
     * bnelems:  number of ptypes in user buffer, buf
     * nbytes:   number of bytes (in external data representation) to read from
     *           or write to the file
     * el_size:  byte size of ptype
     * buftype_is_contig: whether buftype is contiguous
     */
    err = ncmpii_buftype_decode(varp->ndims, varp->xtype, count, bufcount,
                                buftype, &ptype, &el_size, &bnelems,
                                &nbytes, &buftype_is_contig);
    if (err != NC_NOERR) return err;

#ifndef ENABLE_LARGE_REQ
    if (nbytes > INT_MAX) DEBUG_RETURN_ERROR(NC_EMAX_REQ)
#endif

    if (bnelems == 0) {
        /* zero-length request, mark this as a NULL request */
        if (!isSameGroup && reqid != NULL)
            /* only if this is not part of a group request */
            *reqid = NC_REQ_NULL;
        return NC_NOERR;
    }

    /* check if type conversion and Endianness byte swap is needed */
    need_convert = ncmpii_need_convert(ncp->format, varp->xtype, ptype);
    need_swap    = ncmpii_need_swap(varp->xtype, ptype);

    /* check whether this is a true varm call, if yes, imaptype will be a
     * newly created MPI derived data type, otherwise MPI_DATATYPE_NULL
     */
    err = ncmpii_create_imaptype(varp->ndims, count, imap, ptype, &imaptype);
    if (err != NC_NOERR) return err;

    if (fIsSet(reqMode, NC_REQ_WR)) { /* pack request to xbuf */
#if 1
        /* when user buf is used as xbuf, we need to byte-swap buf
         * back to its original contents */
        xbuf = buf;
        need_swap_back_buf = 1;

        if (fIsSet(reqMode, NC_REQ_NBB)) {
            /* for bput call, check if the remaining buffer space is sufficient
             * to accommodate this request and obtain a space for xbuf
             */
            if (ncp->abuf->size_allocated - ncp->abuf->size_used < nbytes)
                DEBUG_RETURN_ERROR(NC_EINSUFFBUF)
            err = abuf_malloc(ncp, nbytes, &xbuf, &abuf_index);
            if (err != NC_NOERR) return err;
            need_swap_back_buf = 0;
        }
        else {
            if (!buftype_is_contig || imaptype != MPI_DATATYPE_NULL ||
                need_convert ||
#ifdef DISABLE_IN_PLACE_SWAP
                need_swap
#else
                nbytes <= NC_BYTE_SWAP_BUFFER_SIZE
#endif
            ) {
                xbuf = NCI_Malloc((size_t)nbytes);
                free_xbuf = 1;
                if (xbuf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
                need_swap_back_buf = 0;
            }
        }

        /* pack user buffer, buf, to xbuf which will be used to write to file */
        err = ncmpio_pack_xbuf(ncp->format, varp, bufcount, buftype,
                               buftype_is_contig, bnelems, ptype, imaptype,
                               need_convert, need_swap, nbytes, buf, xbuf);
        if (err != NC_NOERR && err != NC_ERANGE) {
            if (fIsSet(reqMode, NC_REQ_NBB)) abuf_dealloc(ncp, abuf_index);
            else                             NCI_Free(xbuf);
            return err;
        }
#else
        void *cbuf=NULL, *lbuf=NULL;
        int position;

        /* attached buffer allocation logic
         * if (fIsSet(reqMode, NC_REQ_NBB))
         *     if contig && no imap && no convert
         *         buf   ==   lbuf   ==   cbuf    ==     xbuf memcpy-> abuf
         *                                               abuf
         *     if contig && no imap &&    convert
         *         buf   ==   lbuf   ==   cbuf convert-> xbuf == abuf
         *                                               abuf
         *     if contig &&    imap && no convert
         *         buf   ==   lbuf pack-> cbuf    ==     xbuf == abuf
         *                                abuf
         *     if contig &&    imap &&    convert
         *         buf   ==   lbuf pack-> cbuf convert-> xbuf == abuf
         *                                               abuf
         *  if noncontig && no imap && no convert
         *         buf pack-> lbuf   ==   cbuf    ==     xbuf == abuf
         *                    abuf
         *  if noncontig && no imap &&    convert
         *         buf pack-> lbuf   ==   cbuf convert-> xbuf == abuf
         *                                               abuf
         *  if noncontig &&    imap && no convert
         *         buf pack-> lbuf pack-> cbuf    ==     xbuf == abuf
         *                                abuf
         *  if noncontig &&    imap &&    convert
         *         buf pack-> lbuf pack-> cbuf convert-> xbuf == abuf
         *                                               abuf
         */

        MPI_Offset ibufsize = bnelems * el_size;
        if (ibufsize != (int)ibufsize) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

        /* Step 1: if buftype is not contiguous, i.e. a noncontiguous MPI
         * derived datatype, pack buf into a contiguous buffer, lbuf,
         */
        if (!buftype_is_contig) { /* buftype is not contiguous */
            if (imaptype == MPI_DATATYPE_NULL && !need_convert)
                /* in this case, lbuf will later become xbuf */
                lbuf = xbuf;
            else {
                /* in this case, allocate lbuf and it will be freed before
                 * constructing xbuf */
                lbuf = NCI_Malloc((size_t)ibufsize);
                if (lbuf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
            }

            /* pack buf into lbuf based on buftype */
            if (bufcount > INT_MAX) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            position = 0;
            MPI_Pack(buf, (int)bufcount, buftype, lbuf, (int)ibufsize,
                     &position, MPI_COMM_SELF);
        }
        else /* for contiguous case, we reuse buf */
            lbuf = buf;

        /* Step 2: if imap is non-contiguous, pack lbuf to cbuf */
        if (imaptype != MPI_DATATYPE_NULL) { /* true varm */
            if (!need_convert)
                /* in this case, cbuf will later become xbuf */
                cbuf = xbuf;
            else {
                /* in this case, allocate cbuf and cbuf will be freed before
                 * constructing xbuf */
                cbuf = NCI_Malloc((size_t)ibufsize);
                if (cbuf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
            }

            /* pack lbuf to cbuf based on imaptype */
            position = 0;
            MPI_Pack(lbuf, 1, imaptype, cbuf, (int)ibufsize, &position,
                     MPI_COMM_SELF);
            MPI_Type_free(&imaptype);

            /* lbuf is no longer needed */
            if (lbuf != buf) NCI_Free(lbuf);
        }
        else /* not a true varm call: reuse lbuf */
            cbuf = lbuf;

        /* Step 3: type-convert and byte-swap cbuf to xbuf, and xbuf will be
         * used in MPI write function to write to file
         */
        if (need_convert) {
            /* user buf type does not match nc var type defined in file */
            void *fillp; /* fill value in internal representation */

            /* find the fill value */
            fillp = NCI_Malloc((size_t)varp->xsz);
            ncmpio_inq_var_fill(varp, fillp);

            /* datatype conversion + byte-swap from cbuf to xbuf */
            DATATYPE_PUT_CONVERT(ncp->format, varp->xtype, xbuf, cbuf, bnelems,
                                 ptype, fillp, err)
            NCI_Free(fillp);

            /* The only error codes returned from DATATYPE_PUT_CONVERT are
             * NC_EBADTYPE or NC_ERANGE. Bad varp->xtype and itype have been
             * sanity checked at the dispatchers, so NC_EBADTYPE is not
             * possible. Thus, the only possible error is NC_ERANGE.
             * NC_ERANGE can be caused by one or more elements of buf that is
             * out of range representable by the external data type, it is not
             * considered a fatal error. The request must continue to finish.
             */
            if (cbuf != buf) NCI_Free(cbuf);
#if 0
            if (err != NC_NOERR && err != NC_ERANGE) {
                if (fIsSet(reqMode, NC_REQ_NBB)) abuf_dealloc(ncp, abuf_index);
                else                             NCI_Free(xbuf);
                return err;
            }
#endif
        }
        else {
/*
            if (xbuf == NULL) xbuf = cbuf;
            else if (cbuf == buf) memcpy(xbuf, cbuf, (size_t)nbytes);
*/
            if (cbuf == buf && xbuf != buf) memcpy(xbuf, cbuf, (size_t)nbytes);

            if (need_swap) {
                /* perform array in-place byte swap on xbuf */
                ncmpii_in_swapn(xbuf, bnelems, varp->xsz);

                if (xbuf == buf) need_swap_back_buf = 1;
                /* when user buf is used as xbuf, we need to byte-swap buf
                 * back to its original contents */
            }
        }
#endif
    }
    else { /* read request */
        /* Type conversion and byte swap for read are done at wait call, we
         * need bnelems to reverse the steps as done in write case
         */
        if (buftype_is_contig && imaptype == MPI_DATATYPE_NULL && !need_convert)
            xbuf = buf;  /* there is no buffered read (bget_var, etc.) */
        else {
            xbuf = NCI_Malloc((size_t)nbytes);
            free_xbuf = 1;
        }
    }

    if (fIsSet(reqMode, NC_REQ_WR)) {
        /* allocate or expand write request array */
        int add_reqs = IS_RECVAR(varp) ? (int)count[0] : 1;
        int rem = ncp->numPutReqs % NC_REQUEST_CHUNK;
        size_t req_alloc = ncp->numPutReqs;

        if (add_reqs > NC_REQUEST_CHUNK)
            req_alloc += add_reqs;
        else if (rem == 0)
            req_alloc += NC_REQUEST_CHUNK;
        else if (rem + add_reqs > NC_REQUEST_CHUNK)
            req_alloc += NC_REQUEST_CHUNK - rem + NC_REQUEST_CHUNK;
        else
            req_alloc = 0;

        if (req_alloc > 0)
            ncp->put_list = (NC_req*) NCI_Realloc(ncp->put_list,
                                                  req_alloc * sizeof(NC_req));
        req = ncp->put_list + ncp->numPutReqs;

        /* the new request ID will be an even number (max of write ID + 2) */
        req->id = 0;
        if (ncp->numPutReqs > 0)
            req->id = ncp->put_list[ncp->numPutReqs-1].id + 2;

        ncp->numPutReqs++;
    }
    else {  /* read request */
        /* allocate or expand read request array */
        int add_reqs = IS_RECVAR(varp) ? (int)count[0] : 1;
        int rem = ncp->numGetReqs % NC_REQUEST_CHUNK;
        size_t req_alloc = ncp->numGetReqs;

        if (add_reqs > NC_REQUEST_CHUNK)
            req_alloc += add_reqs;
        else if (rem == 0)
            req_alloc += NC_REQUEST_CHUNK;
        else if (rem + add_reqs > NC_REQUEST_CHUNK)
            req_alloc += NC_REQUEST_CHUNK - rem + NC_REQUEST_CHUNK;
        else
            req_alloc = 0;

        if (req_alloc > 0)
            ncp->get_list = (NC_req*) NCI_Realloc(ncp->get_list,
                                                  req_alloc * sizeof(NC_req));
        req = ncp->get_list + ncp->numGetReqs;

        /* the new request ID will be an odd number (max of read ID + 2) */
        req->id = 1;
        if (ncp->numGetReqs > 0)
            req->id = ncp->get_list[ncp->numGetReqs-1].id + 2;

        ncp->numGetReqs++;
    }

    /* if isSameGroup, then this request is from i_varn API */
    if (isSameGroup && reqid != NULL)
        req->id = *reqid;

    req->flag = 0;
    if (buftype_is_contig)  fSet(req->flag, NC_REQ_BUF_TYPE_IS_CONTIG);
    if (need_convert)       fSet(req->flag, NC_REQ_BUF_TYPE_CONVERT);
    if (fIsSet(reqMode, NC_REQ_WR)) {
        if (need_swap_back_buf) fSet(req->flag, NC_REQ_BUF_BYTE_SWAP);
    }
    else {
        if (need_swap) fSet(req->flag, NC_REQ_BUF_BYTE_SWAP);
    }

    req->varp        = varp;
    req->buf         = buf;
    req->xbuf        = xbuf;
    req->bufcount    = bufcount;
    req->ptype       = ptype;
    req->imaptype    = imaptype;
    req->abuf_index  = abuf_index;
    req->userBuf     = NULL;
    req->status      = NULL;

    /* for read requst and buftype is not contiguous, we duplicate buftype for
     * later in the wait call to unpack buffer based on buftype
     */
    if (fIsSet(reqMode, NC_REQ_RD) && !buftype_is_contig)
        MPI_Type_dup(buftype, &req->buftype);
    else
        req->buftype = MPI_DATATYPE_NULL;

    /* allocate a single array to store start/count/stride */
    if (stride != NULL)
        req->start = (MPI_Offset*) NCI_Malloc(varp->ndims*3*SIZEOF_MPI_OFFSET);
    else {
        req->start = (MPI_Offset*) NCI_Malloc(varp->ndims*2*SIZEOF_MPI_OFFSET);
        fSet(req->flag, NC_REQ_STRIDE_NULL);
    }

    /* copy over start/count/stride */
    j = 0;
    for (i=0; i<varp->ndims; i++) req->start[j++] = start[i];
    for (i=0; i<varp->ndims; i++) req->start[j++] = count[i];
    if (stride != NULL) {
        for (i=0; i<varp->ndims; i++) req->start[j++] = stride[i];
    }

    if (IS_RECVAR(varp) && count[0] > 1) {
        /* If this is a record variable and the number of requesting records is
         * more than 1, we split this request into multiple (sub)requests, one
         * for each record. The number of records is only preserved in the lead
         * request count[0]. The rest sub-requests count[0] are set to 1. The
         * one preserved in lead request will be used to byte-swap the I/O
         * buffer once the file read/write is complete.
         */

        /* add (count[0]-1) number of (sub)requests */
        add_record_requests(req, stride);

        if (fIsSet(reqMode, NC_REQ_WR)) ncp->numPutReqs += count[0] - 1;
        else                            ncp->numGetReqs += count[0] - 1;
    }

    /* mark req as the lead request */
    fSet(req->flag, NC_REQ_LEAD);

    /* only lead request may free xbuf */
    if (free_xbuf) fSet(req->flag, NC_REQ_XBUF_TO_BE_FREED);

    /* return the request ID */
    if (reqid != NULL) *reqid = req->id;

    return err;
}

include(`utils.m4')dnl
dnl
dnl IGETPUT_API(get/put)
dnl
define(`IGETPUT_API',dnl
`dnl
/*----< ncmpio_i$1_var() >---------------------------------------------------*/
/* start  can be NULL only when api is NC_VAR
 * count  can be NULL only when api is NC_VAR or NC_VAR1
 * stride can be NULL only when api is NC_VAR, NC_VAR1, or NC_VARA
 * imap   can be NULL only when api is NC_VAR, NC_VAR1, NC_VARA, or NC_VARS
 * bufcount is >= 0 when called from a flexible API, is -1 when called from a
 *         high-level API and in this case buftype is an MPI primitive
 *         datatype.
 * buftype is an MPI primitive data type (corresponding to the internal data
 *         type of buf, e.g. short in ncmpi_put_short is mapped to MPI_SHORT)
 *         if called from a high-level APIs. When called from a flexible API
 *         it can be an MPI derived data type or MPI_DATATYPE_NULL. If it is
 *         MPI_DATATYPE_NULL, then it means the data type of buf in memory
 *         matches the variable external data type. In this case, bufcount is
 *         ignored.
 * reqMode indicates modes (NC_REQ_COLL/NC_REQ_INDEP/NC_REQ_WR etc.)
 */
int
ncmpio_i$1_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               ifelse(`$1',`put',`const') void *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    NC *ncp=(NC*)ncdp;

    /* Note sanity check for ncdp and varid has been done in dispatchers */

    return ncmpio_igetput_varm(ncp, ncp->vars.value[varid], start, count,
                               stride, imap, (void*)buf, bufcount, buftype,
                               reqid, reqMode, 0);
}
')dnl
dnl

IGETPUT_API(put)
IGETPUT_API(get)

