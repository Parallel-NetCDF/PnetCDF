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

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include <string.h> /* memcpy() */
#include <mpi.h>

#include "nc.h"
#include "ncx.h"
#include "ncmpidtype.h"
#include "macro.h"


/*----< ncmpii_abuf_malloc() >------------------------------------------------*/
/* allocate memory space from the attached buffer pool */
static int
ncmpii_abuf_malloc(NC *ncp, MPI_Offset nbytes, void **buf, int *abuf_index)
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

/*----< add_record_requests() >----------------------------------------------*/
/* check if this is a record variable. if yes, add new requests for each
 * record into the list. Hereinafter, treat each request as a non-record
 * variable request
 */
static int
add_record_requests(NC               *ncp,
                    NC_var           *varp,
                    NC_req           *reqs,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[])
{
    int    i, j;
    size_t dims_chunk = (size_t)varp->ndims * SIZEOF_MPI_OFFSET;
    MPI_Offset record_bufcount, rec_bufsize;

    record_bufcount = 1;
    for (i=1; i<varp->ndims; i++)
        record_bufcount *= reqs[0].count[i];
    rec_bufsize = varp->xsz * record_bufcount;

    /* append each record to the end of list */
    for (i=1; i<reqs[0].count[0]; i++) {

        reqs[i] = reqs[0]; /* inherit most attributes from reqs[0]
                            * except below ones, including the ones need
                            * malloc
                            */

        if (stride != NULL)
            reqs[i].start = (MPI_Offset*) NCI_Malloc(dims_chunk*3);
        else
            reqs[i].start = (MPI_Offset*) NCI_Malloc(dims_chunk*2);

        reqs[i].count = reqs[i].start + varp->ndims;

        if (stride != NULL) {
            reqs[i].stride    = reqs[i].count + varp->ndims;
            reqs[i].start[0]  = reqs[0].start[0] + stride[0] * i;
            reqs[i].stride[0] = reqs[0].stride[0];
        } else {
            reqs[i].stride   = NULL;
            reqs[i].start[0] = reqs[0].start[0] + i;
        }

        reqs[i].count[0] = 1;
        for (j=1; j<varp->ndims; j++) {
            reqs[i].start[j]  = reqs[0].start[j];
            reqs[i].count[j]  = reqs[0].count[j];
            if (stride != NULL)
                reqs[i].stride[j] = reqs[0].stride[j];
        }

        /* xbuf cannot be NULL    assert(reqs[0].xbuf != NULL); */

        reqs[i].bnelems  = record_bufcount;
        reqs[i].buf      = (char*)(reqs[i-1].buf)  + rec_bufsize;
        reqs[i].xbuf     = (char*)(reqs[i-1].xbuf) + rec_bufsize;
        reqs[i].num_recs = 0;  /* not the lead request */

        /* reqs[i].bufcount and reqs[i].buftype will not be used in
         * wait call, only the lead request's matters */
    }

    /* reset the lead request to one record at a time */
    reqs[0].bnelems  = record_bufcount;
    reqs[0].count[0] = 1;

    return NC_NOERR;
}

/*----< ncmpii_igetput_varm() >-----------------------------------------------*/
int
ncmpii_igetput_varm(NC               *ncp,
                    NC_var           *varp,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    const MPI_Offset  imap[],
                    void             *buf,      /* user buffer */
                    MPI_Offset        bufcount,
                    MPI_Datatype      buftype,
                    int              *reqid,    /* out, can be NULL */
                    int               rw_flag,
                    int               use_abuf,    /* if use attached buffer */
                    int               isSameGroup) /* if part of a varn group */
{
    void *xbuf=NULL, *cbuf=NULL, *lbuf=NULL;
    int err=NC_NOERR, status=NC_NOERR, warning=NC_NOERR;
    int i, abuf_index=-1, el_size, buftype_is_contig;
    int need_convert, need_swap, need_swap_back_buf=0, cdf_ver;
    size_t  dims_chunk;
    MPI_Offset bnelems=0, nbytes;
    MPI_Datatype ptype, imaptype=MPI_DATATYPE_NULL;
    NC_req *req;

    /* check NC_ECHAR error and calculate the followings:
     * ptype: element data type (MPI primitive type) in buftype
     * bufcount: If it is -1, then this is called from a high-level API and in
     * this case buftype will be an MPI primitive data type. If not, then this
     * is called from a flexible API. In that case, we recalculate bufcount to
     * match with count[].
     * bnelems: number of ptypes in user buffer
     * nbytes: number of bytes (in external data representation) to read/write
     * from/to the file
     * el_size: size of ptype
     * buftype_is_contig: whether buftype is contiguous
     */
    err = ncmpii_calc_datatype_elems(ncp, varp, start, count, stride, rw_flag,
                                     buftype, &ptype, &bufcount, &bnelems,
                                     &nbytes, &el_size, &buftype_is_contig);
    if (err == NC_EIOMISMATCH) DEBUG_ASSIGN_ERROR(warning, err)
    else if (err != NC_NOERR) return err;

    if (bnelems == 0) {
        /* zero-length request, mark this as a NULL request */
        if (!isSameGroup && reqid != NULL)
            /* only if this is not part of a group request */
            *reqid = NC_REQ_NULL;
        return ((warning != NC_NOERR) ? warning : NC_NOERR);
    }

    /* for bput call, check if the remaining buffer space is sufficient
     * to accommodate this request
     */
    if (rw_flag == WRITE_REQ && use_abuf &&
        ncp->abuf->size_allocated - ncp->abuf->size_used < nbytes)
        DEBUG_RETURN_ERROR(NC_EINSUFFBUF)

    if (fIsSet(ncp->flags, NC_64BIT_DATA))        cdf_ver = 5;  /* CDF-5 */
    else if (fIsSet(ncp->flags, NC_64BIT_OFFSET)) cdf_ver = 2;  /* CDF-2 */
    else                                          cdf_ver = 1;  /* CDF-1 */

    /* check if type conversion and Endianness byte swap is needed */
    need_convert = ncmpii_need_convert(cdf_ver, varp->type, ptype);
    need_swap    = ncmpii_need_swap(varp->type, ptype);

    if (imap != NULL) {
        /* check whether this is a true varm call, if yes, imaptype will be a
         * newly created MPI derived data type, otherwise MPI_DATATYPE_NULL
         */
        err = ncmpii_create_imaptype(varp, count, imap, bnelems, el_size,
                                     ptype, &imaptype);
        if (err != NC_NOERR) return err;
    }

    if (rw_flag == WRITE_REQ) { /* pack request to xbuf */
        int position, abuf_allocated=0;
        MPI_Offset outsize=bnelems*el_size;
        /* assert(bnelems > 0); */
        if (outsize != (int)outsize) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

        /* attached buffer allocation logic
         * if (use_abuf)
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

        /* Step 1: pack buf into a contiguous buffer, lbuf, if buftype is
         * not contiguous, i.e. a noncontiguous MPI derived datatype
         */
        if (!buftype_is_contig) { /* buftype is not contiguous */
            if (bufcount != (int)bufcount) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

            /* allocate lbuf */
            if (use_abuf && imaptype == MPI_DATATYPE_NULL && !need_convert) {
                status = ncmpii_abuf_malloc(ncp, nbytes, &lbuf, &abuf_index);
                if (status != NC_NOERR) return status;
                abuf_allocated = 1;
            }
            else lbuf = NCI_Malloc((size_t)outsize);

            /* pack buf into lbuf based on buftype */
            position = 0;
            MPI_Pack(buf, (int)bufcount, buftype, lbuf, (int)outsize,
                     &position, MPI_COMM_SELF);
        }
        else /* for contiguous case, we reuse buf */
            lbuf = buf;

        /* Step 2: pack lbuf to cbuf if imap is non-contiguous */
        if (imaptype != MPI_DATATYPE_NULL) { /* true varm */
            /* allocate cbuf */
            if (use_abuf && !need_convert) {
                assert(abuf_allocated == 0);
                status = ncmpii_abuf_malloc(ncp, nbytes, &cbuf, &abuf_index);
                if (status != NC_NOERR) {
                    if (lbuf != buf) NCI_Free(lbuf);
                    return status;
                }
                abuf_allocated = 1;
            }
            else cbuf = NCI_Malloc((size_t)outsize);

            /* pack lbuf to cbuf based on imaptype */
            position = 0;
            MPI_Pack(lbuf, 1, imaptype, cbuf, (int)outsize, &position,
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

        /* when user buf type != nc var type defined in file */
        if (need_convert) {
            if (use_abuf) { /* use attached buffer to allocate xbuf */
                assert(abuf_allocated == 0);
                status = ncmpii_abuf_malloc(ncp, nbytes, &xbuf, &abuf_index);
                if (status != NC_NOERR) {
                    if (cbuf != buf) NCI_Free(cbuf);
                    return status;
                }
                abuf_allocated = 1;
            }
            else xbuf = NCI_Malloc((size_t)nbytes);

            /* datatype conversion + byte-swap from cbuf to xbuf */
            DATATYPE_PUT_CONVERT(cdf_ver, varp->type, xbuf, cbuf, bnelems, ptype, status)
            /* NC_ERANGE can be caused by a subset of buf that is out of range
             * of the external data type, it is not considered a fatal error.
             * The request must continue to finish.
             */
            if (status != NC_NOERR && status != NC_ERANGE) {
                if (cbuf != buf)  NCI_Free(cbuf);
                if (xbuf != NULL) NCI_Free(xbuf);
                return status;
            }
        }
        else {
            if (use_abuf && buftype_is_contig && imaptype == MPI_DATATYPE_NULL){
                assert(abuf_allocated == 0);
                status = ncmpii_abuf_malloc(ncp, nbytes, &xbuf, &abuf_index);
                if (status != NC_NOERR) {
                    if (cbuf != buf) NCI_Free(cbuf);
                    return status;
                }
                memcpy(xbuf, cbuf, (size_t)nbytes);
            }
            else xbuf = cbuf;

            if (need_swap) {
#ifdef DISABLE_IN_PLACE_SWAP
                if (xbuf == buf)
#else
                if (xbuf == buf && nbytes <= NC_BYTE_SWAP_BUFFER_SIZE)
#endif
                {
                    /* allocate xbuf and copy buf to xbuf, before byte-swap */
                    xbuf = NCI_Malloc((size_t)nbytes);
                    memcpy(xbuf, buf, (size_t)nbytes);
                }
                /* perform array in-place byte swap on xbuf */
                ncmpii_in_swapn(xbuf, bnelems, ncmpix_len_nctype(varp->type));

                if (xbuf == buf) need_swap_back_buf = 1;
                /* user buf needs to be swapped back to its original contents */
            }
        }
        /* cbuf is no longer needed */
        if (cbuf != buf && cbuf != xbuf) NCI_Free(cbuf);
    }
    else { /* rw_flag == READ_REQ */
        /* Type conversion and byte swap for read are done at wait call, we
         * need bnelems to reverse the steps as done in write case
         */
        if (buftype_is_contig && imaptype == MPI_DATATYPE_NULL && !need_convert)
            xbuf = buf;  /* there is no buffered read (bget_var, etc.) */
        else
            xbuf = NCI_Malloc((size_t)nbytes);
    }

    if (rw_flag == WRITE_REQ) {
        /* allocate write/read request array */
        if (ncp->numPutReqs % NC_REQUEST_CHUNK == 0)
            ncp->put_list = (NC_req*) NCI_Realloc(ncp->put_list,
                                      (ncp->numPutReqs + NC_REQUEST_CHUNK) *
                                      sizeof(NC_req));
        req = ncp->put_list + ncp->numPutReqs;

        /* the new request ID will be an even number (max of write ID + 2) */
        req->id = 0;
        if (ncp->numPutReqs > 0)
            req->id = ncp->put_list[ncp->numPutReqs-1].id + 2;

        ncp->numPutReqs++;
    }
    else {  /* READ_REQ */
        /* allocate write/read request array */
        if (ncp->numGetReqs % NC_REQUEST_CHUNK == 0)
            ncp->get_list = (NC_req*) NCI_Realloc(ncp->get_list,
                                      (ncp->numGetReqs + NC_REQUEST_CHUNK) *
                                      sizeof(NC_req));
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

    req->varp               = varp;
    req->buf                = buf;
    req->xbuf               = xbuf;
    req->bnelems            = bnelems;
    req->bufcount           = bufcount;
    req->ptype              = ptype;
    req->buftype_is_contig  = buftype_is_contig;
    req->need_swap_back_buf = need_swap_back_buf;
    req->imaptype           = imaptype;
    req->abuf_index         = abuf_index;
    req->tmpBuf             = NULL;
    req->userBuf            = NULL;
    req->status             = NULL;
    req->num_recs           = 1;   /* For record variable, this will be set to
                                    * the number of records requested. For
                                    * fixed-size variable, this will be 1.
                                    */

    /* only when read and buftype is not contiguous, we duplicate buftype for
     * later in the wait call to unpack buffer based on buftype
     */
    if (rw_flag == READ_REQ && !buftype_is_contig)
        MPI_Type_dup(buftype, &req->buftype);
    else
        req->buftype = MPI_DATATYPE_NULL;

    /* allocate start/count/stride arrays */
    dims_chunk = (size_t)varp->ndims * SIZEOF_MPI_OFFSET;
    if (stride != NULL)
        req->start = (MPI_Offset*) NCI_Malloc(dims_chunk*3);
    else
        req->start = (MPI_Offset*) NCI_Malloc(dims_chunk*2);

    req->count = req->start + varp->ndims;

    if (stride != NULL)
        req->stride = req->count + varp->ndims;
    else
        req->stride = NULL;

    /* set the values for start/count/stride */
    for (i=0; i<varp->ndims; i++) {
        req->start[i] = start[i];
        req->count[i] = count[i];
        if (stride != NULL)
            req->stride[i] = stride[i];
    }

    /* if this is a record variable and number of requesting records is > 1,
     * we split the request, one for each record
     */
    if (IS_RECVAR(varp) && req->count[0] > 1) {
        req->num_recs = req->count[0];

        add_record_requests(ncp, varp, req, start, count, stride);
        /* req->count[0] has been changed to 1 */

        if (rw_flag == WRITE_REQ) ncp->numPutReqs += req->num_recs - 1;
        else                      ncp->numGetReqs += req->num_recs - 1;
    }

    /* return the request ID */
    if (reqid != NULL) *reqid = req->id;

    return ((warning != NC_NOERR) ? warning : status);
}

include(`foreach.m4')dnl
include(`utils.m4')dnl
dnl
define(`APINAME',`ifelse(`$3',`',`ncmpi_i$1_var$2',`ncmpi_i$1_var$2_$3')')dnl
dnl
dnl IGETPUT_API(get/put, kind, itype)
dnl
define(`IGETPUT_API',dnl
`dnl
/*----< APINAME($1,$2,$3)() >------------------------------------------------*/
int
APINAME($1,$2,$3)(int ncid, int varid, ArgKind($2) BufArgs($1,$3), int *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp=NULL;
    ifelse(`$2', `',  `MPI_Offset *start, *count;',
           `$2', `1', `MPI_Offset *count;')

    if (reqid != NULL) *reqid = NC_REQ_NULL;
    status = ncmpii_sanity_check(ncid, varid, ArgStartCount($2),
                                 ifelse(`$3', `', `bufcount', `0'),
                                 API_KIND($2), 0, ReadWrite($1),
                                 NONBLOCKING_IO, &ncp, &varp);
    if (status != NC_NOERR) return status;

    ifelse(`$2', `',  `GET_FULL_DIMENSIONS(start, count)',
           `$2', `1', `GET_ONE_COUNT(count)')

    /* APINAME($1,$2,$3) is a special case of APINAME($1,m,$3) */
    status = ncmpii_igetput_varm(ncp, varp, start, count, ArgStrideMap($2),
                                 (void*)buf,
                                 ifelse(`$3', `', `bufcount, buftype',
                                                  `-1, ITYPE2MPI($3)'),
                                 reqid, ReadWrite($1), 0, 0);
    ifelse(`$2', `', `NCI_Free(start);', `$2', `1', `NCI_Free(count);')
    return status;
}
')dnl
dnl
/*---- PnetCDF flexible APIs ------------------------------------------------*/
foreach(`kind', (, 1, a, s, m),
        `foreach(`putget', (put, get),
                 `IGETPUT_API(putget,kind,)'
)')

/*---- PnetCDF high-level APIs ----------------------------------------------*/
foreach(`kind', (, 1, a, s, m),
        `foreach(`putget', (put, get),
                 `foreach(`itype', (ITYPE_LIST),
                          `IGETPUT_API(putget,kind,itype)'
)')')
