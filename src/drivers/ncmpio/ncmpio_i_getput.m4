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
#include <assert.h>

#include <string.h> /* memcpy() */
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< ncmpio_abuf_malloc() >-----------------------------------------------*/
/* allocate memory space from the attached buffer pool */
int
ncmpio_abuf_malloc(NC *ncp, MPI_Offset nbytes, void **buf, int *abuf_index)
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

/*----< ncmpio_abuf_dealloc() >----------------------------------------------*/
/* deallocate (actually un-register) memory space from the attached buffer
 * pool
 */
int
ncmpio_abuf_dealloc(NC *ncp, int abuf_index)
{
    assert(abuf_index == ncp->abuf->tail - 1);

    /* mark the tail entry un-used */
    ncp->abuf->size_used -= ncp->abuf->occupy_table[abuf_index].req_size;
    ncp->abuf->occupy_table[abuf_index].req_size = 0;
    ncp->abuf->occupy_table[abuf_index].is_used  = 0;
    ncp->abuf->tail--;

    return NC_NOERR;
}

/*----< ncmpio_add_record_requests() >---------------------------------------*/
/* check if this is a record variable. if yes, add a new request for each
 * record into the list. Hereinafter, treat each request as a non-record
 * variable request
 */
int
ncmpio_add_record_requests(NC_lead_req      *lead_list,
                           NC_req           *reqs,
                           MPI_Offset        num_recs,
                           const MPI_Offset *stride)
{
    char       *xbuf;
    int         i, ndims;
    size_t      dims_chunk;
    MPI_Offset  rec_bufsize, *count;
    NC_var     *varp;

    varp = lead_list[reqs[0].lead_off].varp;
    ndims = varp->ndims;
    dims_chunk = (stride == NULL) ? ndims * 2 : ndims * 3;

    /* start[]/count[]/stride[] have been copied to reqs[0] */

    count = reqs[0].start + ndims;
    count[0] = 1; /* each non-lead request accesses one record only */

    /* calculate request size in bytes */
    rec_bufsize = reqs[0].nelems * varp->xsz;

    /* add new requests, one per record */
    xbuf = (char*)reqs[0].xbuf + rec_bufsize;
    for (i=1; i<num_recs; i++) {
        /* copy start/count/stride */
        reqs[i].start = reqs[i-1].start + dims_chunk;
        memcpy(reqs[i].start, reqs[i-1].start, dims_chunk * SIZEOF_MPI_OFFSET);

        /* jump to next stride */
        reqs[i].start[0] += (stride == NULL) ? 1 : stride[0];

        reqs[i].nelems    = reqs[0].nelems;
        reqs[i].lead_off  = reqs[0].lead_off;
        reqs[i].xbuf      = xbuf;
        xbuf             += rec_bufsize;
    }

    return NC_NOERR;
}

/*----< ncmpio_igetput_varm() >-----------------------------------------------*/
int
ncmpio_igetput_varm(NC               *ncp,
                    NC_var           *varp,
                    const MPI_Offset  start[],  /* not NULL, if not scalar */
                    const MPI_Offset  count[],  /* not NULL */
                    const MPI_Offset  stride[], /* may be NULL */
                    const MPI_Offset  imap[],   /* may be NULL */
                    void             *buf,      /* user buffer */
                    MPI_Offset        bufcount,
                    MPI_Datatype      buftype,
                    int              *reqid,    /* out, can be NULL */
                    int               reqMode)
{
    void *xbuf=NULL;
    int i, err=NC_NOERR, abuf_index=-1, isize, xsize, new_nreqs, rem;
    int mpireturn, buftype_is_contig=1, need_convert, free_xbuf=0;
    int need_swap, can_swap_in_place, need_swap_back_buf=0;
    MPI_Offset nelems=0, nbytes, *ptr;
    MPI_Datatype itype, xtype, imaptype;
    NC_lead_req *lead_req;
    NC_req *req;

    /* validity of start and count has been checked at dispatcher layer */

    nelems = 1; /* total number of array elements of this vara request */
    for (i=0; i<varp->ndims; i++)
        nelems *= count[i];

    /* xtype is the MPI element type in external representation, xsize is its
     * size in bytes. Similarly, itype and isize for internal representation.
     */
    xtype = ncmpii_nc2mpitype(varp->xtype);
    mpireturn = MPI_Type_size(xtype, &xsize);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size");

    if (buftype == MPI_DATATYPE_NULL) {
        /* In this case, bufcount is ignored and the internal buffer data type
         * match the external variable data type. No data conversion will be
         * done. In addition, it means buf is contiguous. Hereinafter, buftype
         * is ignored.
         */
        itype = xtype;
        isize = xsize;
    }
    else if (bufcount == NC_COUNT_IGNORE) {
        /* In this case, this subroutine is called from a high-level API and
         * buftype must be one of the MPI predefined datatype. We set itype to
         * buftype. itype is the MPI element type in internal representation.
         * In addition, it means the user buf is contiguous.
         */
        itype = buftype;
        mpireturn = MPI_Type_size(itype, &isize); /* buffer element size */
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size");
    }
    else { /* (bufcount > 0) */
        /* When bufcount > 0, this subroutine is called from a flexible API. If
         * buftype is noncontiguous, we pack buf into xbuf, a contiguous buffer.
         */
        MPI_Offset bnelems=0;

        /* When bufcount > NC_MAX_INT, we construct a datatype to bypass the
         * limitation of MPI file read/write APIs on the argument "count" of
         * type int.  See ncmpio_read_write() in ncmpio_file_io.c
         *
         * Note not all MPI-IO libraries support single requests larger than
         * NC_MAX_INT. In this case, MPI-IO should report an error.
         */

        /* itype (primitive MPI data type) from buftype
         * isize is the size of itype in bytes
         * bnelems is the number of itype elements in one buftype
         */
        err = ncmpii_dtype_decode(buftype, &itype, &isize, &bnelems,
                                  NULL, &buftype_is_contig);
        if (err != NC_NOERR) goto fn_exit;

        /* size in bufcount * buftype must match with counts[] */
        if (bnelems * bufcount != nelems) {
            DEBUG_ASSIGN_ERROR(err, NC_EIOMISMATCH)
            goto fn_exit;
        }
    }

    /* nbytes is the amount of this vara request in bytes */
    nbytes = nelems * xsize;

    /* Skip checking nbytes against NC_MAX_INT. Note not all MPI-IO libraries
     * support single requests larger than NC_MAX_INT. In this case, MPI-IO
     * should report an error.
     */

    /* for nonblocking API, return now if request size is zero */
    if (nbytes == 0) {
        if (reqid != NULL)
            *reqid = NC_REQ_NULL; /* mark this as a NULL request */
        goto fn_exit;
    }

    /* check if type conversion and Endianness byte swap is needed */
    need_convert = ncmpii_need_convert(ncp->format, varp->xtype, itype);
    need_swap    = NEED_BYTE_SWAP(varp->xtype, itype);

    /* check if in-place byte swap can be enabled */
    can_swap_in_place = 1;
    if (need_swap) {
        if (! fIsSet(ncp->flags, NC_MODE_SWAP_OFF)) /* hint set by user */
            can_swap_in_place = 0;
        else if (! fIsSet(ncp->flags, NC_MODE_SWAP_ON)) { /* auto mode */
            if (nbytes > NC_BYTE_SWAP_BUFFER_SIZE)
                can_swap_in_place = 0;
        }
    }

    /* check whether this is a true varm call, if yes, imaptype will be a
     * newly created MPI derived data type from imap[] and itype, otherwise
     * it is set to MPI_DATATYPE_NULL
     */
    err = ncmpii_create_imaptype(varp->ndims, count, imap, itype, &imaptype);
    if (err != NC_NOERR) goto fn_exit;

    if (fIsSet(reqMode, NC_REQ_WR)) { /* pack request to xbuf */
        if (fIsSet(reqMode, NC_REQ_NBB)) {
            /* for bput call, check if the remaining buffer space is sufficient
             * to accommodate this request and allocate a space for xbuf
             */
            if (ncp->abuf->size_allocated - ncp->abuf->size_used < nbytes) {
                DEBUG_ASSIGN_ERROR(err, NC_EINSUFFBUF)
                goto fn_exit;
            }
            err = ncmpio_abuf_malloc(ncp, nbytes, &xbuf, &abuf_index);
            if (err != NC_NOERR) goto fn_exit;
        }
        else {
            if (!buftype_is_contig || imaptype != MPI_DATATYPE_NULL ||
                need_convert || (need_swap && can_swap_in_place == 0)) {
                /* cannot use buf for I/O, must allocate xbuf */
                xbuf = NCI_Malloc((size_t)nbytes);
                free_xbuf = 1;
                if (xbuf == NULL) {
                    DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
                    goto fn_exit;
                }
            }
            else { /* when user buf is used as xbuf, we need to byte-swap buf
                    * back to its original contents */
                xbuf = buf;
                if (need_swap) need_swap_back_buf = 1;
            }
        }

        /* pack user buffer, buf, to xbuf which will be used to write to file.
         * In the meanwhile, perform byte-swap and type-conversion if required.
         */
        err = ncmpio_pack_xbuf(ncp->format, varp, bufcount, buftype,
                               buftype_is_contig, nelems, itype, isize,
                               imaptype, need_convert, need_swap, nbytes, buf,
                               xbuf);
        if (err != NC_NOERR && err != NC_ERANGE) {
            if (fIsSet(reqMode, NC_REQ_NBB))
                ncmpio_abuf_dealloc(ncp, abuf_index);
            else if (free_xbuf)
                NCI_Free(xbuf);
            goto fn_exit;
        }
    }
    else { /* read request */
        /* Type conversion and byte swap for read will be done at wait call. */
        if (buftype_is_contig && imaptype == MPI_DATATYPE_NULL && !need_convert)
            xbuf = buf;  /* there is no buffered read APIs (bget_var, etc.) */
        else {
            xbuf = NCI_Malloc((size_t)nbytes);
            free_xbuf = 1;
        }
    }

    /* add a new nonblocking request to the request queue */

    if (fIsSet(reqMode, NC_REQ_WR)) {
        /* allocate or expand the lead write request queue */
        if (ncp->numLeadPutReqs % NC_REQUEST_CHUNK == 0)
            ncp->put_lead_list = (NC_lead_req*) NCI_Realloc(ncp->put_lead_list,
                                 (ncp->numLeadPutReqs + NC_REQUEST_CHUNK) *
                                 sizeof(NC_lead_req));

        /* allocate or expand the non-lead write request queue */
        new_nreqs = IS_RECVAR(varp) ? (int)count[0] : 1;
        rem = ncp->numPutReqs % NC_REQUEST_CHUNK;
        if (rem) rem = NC_REQUEST_CHUNK - rem;

        if (ncp->put_list == NULL || new_nreqs > rem) {
            size_t req_alloc, nChunks;
            req_alloc = ncp->numPutReqs + new_nreqs;
            nChunks = req_alloc / NC_REQUEST_CHUNK;
            if (req_alloc % NC_REQUEST_CHUNK) nChunks++;
            req_alloc = nChunks * NC_REQUEST_CHUNK * sizeof(NC_req);
            ncp->put_list = (NC_req*) NCI_Realloc(ncp->put_list, req_alloc);
        }

#define SORT_LEAD_LIST_BASED_ON_VAR_BEGIN
#ifdef SORT_LEAD_LIST_BASED_ON_VAR_BEGIN
        /* add the new request to put_lead_list and keep put_lead_list sorted,
         * in an increasing order of variable begin offsets. The best is when
         * user makes varn API calls in an increasing order of variable begin
         * offsets, i.e. fixed-size variables first followed by record
         * variables and in an increasing order of variables IDs. If the
         * pending requests consist of multiple records, then keep the list
         * sorted based on the starting offsets of inidvidual records.
         */
        MPI_Offset req_off = varp->begin;
        if (IS_RECVAR(varp)) req_off += ncp->recsize * start[0];

        for (i=ncp->numLeadPutReqs-1; i>=0; i--) {
            if (ncp->put_lead_list[i].varp->begin <= req_off)
                break;
            /* make space for new lead request */
            ncp->put_lead_list[i+1] = ncp->put_lead_list[i];
            ncp->put_lead_list[i+1].nonlead_off += new_nreqs;
        }
        int lead_off = i + 1;
        lead_req = ncp->put_lead_list + lead_off;

        if (lead_off < ncp->numLeadPutReqs) {
            /* req is starting location to insert new non-lead requests */
            req = ncp->put_list + lead_req->nonlead_off;
            /* make space for new non-lead requests */
            for (i=ncp->numPutReqs-1; i>=lead_req->nonlead_off; i--) {
                ncp->put_list[i+new_nreqs] = ncp->put_list[i];
                ncp->put_list[i+new_nreqs].lead_off++;
            }
        }
        else {
            /* append new lead request at the end of ncp->put_lead_list */
            lead_req->nonlead_off = ncp->numPutReqs;
            /* append new non-lead requests at the end of ncp->put_list */
            req = ncp->put_list + ncp->numPutReqs;
        }
        req->lead_off = lead_off;
#else
        /* append new lead request at the end of ncp->put_lead_list */
        lead_req = ncp->put_lead_list + ncp->numLeadPutReqs;
        lead_req->nonlead_off = ncp->numPutReqs;
        /* append new non-lead requests at the end of ncp->put_list */
        req = ncp->put_list + ncp->numPutReqs;
        req->lead_off = ncp->numLeadPutReqs;
#endif

        lead_req->flag = 0;
        if (need_swap_back_buf) fSet(lead_req->flag, NC_REQ_BUF_BYTE_SWAP);

        /* the new request ID will be an even number (max of write ID + 2) */
        if (ncp->numLeadPutReqs == 0) {
            lead_req->id = 0;
            ncp->maxPutReqID = 0;
        } else {
            ncp->maxPutReqID += 2;
            lead_req->id = ncp->maxPutReqID;
        }

        ncp->numLeadPutReqs++;
        ncp->numPutReqs += new_nreqs;
    }
    else {  /* read request */
        /* allocate or expand the lead read request queue */
        if (ncp->numLeadGetReqs % NC_REQUEST_CHUNK == 0)
            ncp->get_lead_list = (NC_lead_req*) NCI_Realloc(ncp->get_lead_list,
                                 (ncp->numLeadGetReqs + NC_REQUEST_CHUNK) *
                                 sizeof(NC_lead_req));

        /* allocate or expand the non-lead read request queue */
        new_nreqs = IS_RECVAR(varp) ? (int)count[0] : 1;
        rem = ncp->numGetReqs % NC_REQUEST_CHUNK;
        if (rem) rem = NC_REQUEST_CHUNK - rem;

        if (ncp->get_list == NULL || new_nreqs > rem) {
            size_t req_alloc, nChunks;
            req_alloc = ncp->numGetReqs + new_nreqs;
            nChunks = req_alloc / NC_REQUEST_CHUNK;
            if (req_alloc % NC_REQUEST_CHUNK) nChunks++;
            req_alloc = nChunks * NC_REQUEST_CHUNK * sizeof(NC_req);
            ncp->get_list = (NC_req*) NCI_Realloc(ncp->get_list, req_alloc);
        }

#ifdef SORT_LEAD_LIST_BASED_ON_VARID
        /* add the new request to get_lead_list and keep get_lead_list sorted,
         * in an increasing order of variable begin offsets. The best is when
         * user makes varn API calls in an increasing order of variable begin
         * offsets, i.e. fixed-size variables first followed by record
         * variables and in an increasing order of variables IDs.
         */
        for (i=ncp->numLeadGetReqs-1; i>=0; i--) {
            if (ncp->get_lead_list[i].varp->begin <= varp->begin)
                break;
            /* make space for new lead request */
            ncp->get_lead_list[i+1] = ncp->get_lead_list[i];
            ncp->get_lead_list[i+1].nonlead_off += new_nreqs;
        }
        int lead_off = i + 1;
        lead_req = ncp->get_lead_list + lead_off;

        if (lead_off < ncp->numLeadGetReqs) {
            /* req is starting location to insert new non-lead requests */
            req = ncp->get_list + lead_req->nonlead_off;
            /* make space for new non-lead requests */
            for (i=ncp->numGetReqs-1; i>=lead_req->nonlead_off; i--) {
                ncp->get_list[i+new_nreqs] = ncp->get_list[i];
                ncp->get_list[i+new_nreqs].lead_off++;
            }
        }
        else {
            /* append new lead request at the end of ncp->get_lead_list */
            req = ncp->get_list + ncp->numGetReqs;
            /* append new non-lead requests at the end of ncp->get_list */
            lead_req->nonlead_off = ncp->numGetReqs;
        }
        req->lead_off = lead_off;
#else
        /* append new lead request at the end of ncp->get_lead_list */
        lead_req = ncp->get_lead_list + ncp->numLeadGetReqs;
        lead_req->nonlead_off = ncp->numGetReqs;
        /* append new non-lead requests at the end of ncp->get_list */
        req = ncp->get_list + ncp->numGetReqs;
        req->lead_off = ncp->numLeadGetReqs;
#endif

        lead_req->flag = 0;
        if (need_convert) fSet(lead_req->flag, NC_REQ_BUF_TYPE_CONVERT);
        if (need_swap)    fSet(lead_req->flag, NC_REQ_BUF_BYTE_SWAP);

        /* the new request ID will be an odd number (max of read ID + 2) */
        if (ncp->numLeadGetReqs == 0) {
            lead_req->id = 1;
            ncp->maxGetReqID = 1;
        } else {
            ncp->maxGetReqID += 2;
            lead_req->id = ncp->maxGetReqID;
        }

        ncp->numLeadGetReqs++;
        ncp->numGetReqs += new_nreqs;
    }

    /* set other properties for the lead request */

    lead_req->varp        = varp;
    lead_req->buf         = buf;
    lead_req->bufcount    = bufcount;
    lead_req->itype       = itype;
    lead_req->imaptype    = imaptype;
    lead_req->abuf_index  = abuf_index;
    lead_req->status      = NULL;
    lead_req->nelems      = nelems;
    lead_req->xbuf        = xbuf;
    lead_req->buftype     = MPI_DATATYPE_NULL;
    lead_req->nonlead_num = new_nreqs;

    /* only lead request free xbuf (when xbuf != buf) */
    if (free_xbuf) fSet(lead_req->flag, NC_REQ_XBUF_TO_BE_FREED);

    if (stride == NULL) fSet(lead_req->flag, NC_REQ_STRIDE_NULL);
    else {
        for (i=0; i<varp->ndims; i++)
            if (stride[i] > 1) break;
        if (i == varp->ndims) { /* all 1s */
            fSet(lead_req->flag, NC_REQ_STRIDE_NULL);
            stride = NULL;
        }
    }

    /* for read requst and buftype is not contiguous, we duplicate buftype for
     * later in the wait call to unpack buffer based on buftype
     */
    if (buftype_is_contig)
        fSet(lead_req->flag, NC_REQ_BUF_TYPE_IS_CONTIG);
    else if (fIsSet(reqMode, NC_REQ_RD)) {
        mpireturn = MPI_Type_dup(buftype, &lead_req->buftype);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_dup");
            goto fn_exit;
        }
    }

    /* allocate a single array for non-leads to store start/count/stride */
    if (varp->ndims == 0) { /* scalar variable, start may be NULL */
        lead_req->start = NULL;
        req->start = NULL;
    }
    else if (stride == NULL) {
        size_t memChunk = varp->ndims * SIZEOF_MPI_OFFSET;
        if (IS_RECVAR(varp) && count[0] > 1)
            lead_req->start = (MPI_Offset*) NCI_Malloc(memChunk * 2 * count[0]);
        else
            lead_req->start = (MPI_Offset*) NCI_Malloc(memChunk * 2);
        /* copy over start/count/stride */
        req->start = lead_req->start;
        ptr = req->start;
        memcpy(ptr, start, memChunk);
        ptr += varp->ndims;
        memcpy(ptr, count, memChunk);
    }
    else {
        size_t memChunk = varp->ndims * SIZEOF_MPI_OFFSET;
        if (IS_RECVAR(varp) && count[0] > 1)
            lead_req->start = (MPI_Offset*) NCI_Malloc(memChunk * 3 * count[0]);
        else
            lead_req->start = (MPI_Offset*) NCI_Malloc(memChunk * 3);
        /* copy over start/count/stride */
        req->start = lead_req->start;
        ptr = req->start;
        memcpy(ptr, start, memChunk);
        ptr += varp->ndims;
        memcpy(ptr, count, memChunk);
        ptr += varp->ndims;
        memcpy(ptr, stride, memChunk);
    }

    /* set the properties of non-lead request */
    req->xbuf   = xbuf;
    req->nelems = nelems;

    if (IS_RECVAR(varp)) {
        /* save the last record number accessed */
        if (stride == NULL)
            lead_req->max_rec = start[0] + count[0];
        else
            lead_req->max_rec = start[0] + stride[0] * (count[0] - 1) + 1;

        if (count[0] > 1) {
            /* If this is a record variable and the number of requesting
             * records is more than 1, we split this lead request into multiple
             * non-lead requests, one for each record. count[0] in all non-lead
             * requests are set to 1.
             */
            NC_lead_req *lead_list;

            lead_list = (fIsSet(reqMode, NC_REQ_WR)) ? ncp->put_lead_list
                                                     : ncp->get_lead_list;

            req->nelems /= count[0];

            /* add (count[0]-1) number of (sub)requests */
            ncmpio_add_record_requests(lead_list, req, count[0], stride);
        }
    }
    else { /* fixed-size variable */
        lead_req->max_rec = -1;
    }

    /* return the request ID */
    if (reqid != NULL) *reqid = lead_req->id;

fn_exit:
    return err;
}

include(`utils.m4')dnl
dnl
dnl IGETPUT_API(get/put)
dnl
define(`IGETPUT_API',dnl
`dnl
/*----< ncmpio_i$1_var() >---------------------------------------------------*/
/* start  NOT NULL (allocated in dispatcher when api is NC_VAR)
 * count  NOT NULL (allocated in dispatcher when api is NC_VAR or NC_VAR1)
 * stride can be NULL only when api is NC_VAR, NC_VAR1, or NC_VARA
 * imap   can be NULL only when api is NC_VAR, NC_VAR1, NC_VARA, or NC_VARS
 * bufcount If NC_COUNT_IGNORE, then this is called from a high-level API
 *          and buftype must be an MPI primitive data type. Otherwise,
 *          this is called from a flexible API.
 * buftype  is an MPI primitive data type (corresponding to the internal data
 *          type of buf, e.g. short in ncmpi_put_short is mapped to MPI_SHORT)
 *          if called from a high-level APIs. When called from a flexible API
 *          it can be an MPI derived data type or MPI_DATATYPE_NULL. If it is
 *          MPI_DATATYPE_NULL, then it means the data type of buf in memory
 *          matches the variable external data type. In this case, bufcount is
 *          ignored.
 * reqMode  indicates modes (NC_REQ_COLL/NC_REQ_INDEP/NC_REQ_WR etc.)
 */
int
ncmpio_i$1_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start, /* cannot be NULL */
               const MPI_Offset *count, /* cannot be NULL */
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
                               reqid, reqMode);
}
')dnl
dnl

IGETPUT_API(put)
IGETPUT_API(get)

