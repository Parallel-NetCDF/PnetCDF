/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include "nc.h"
#include "ncx.h"
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include "ncmpidtype.h"
#include "macro.h"


/* buffer layers:       
        
        User Level              buf     (user defined buffer of MPI_Datatype)
        MPI Datatype Level      cbuf    (contiguous buffer of ptype)
        NetCDF XDR Level        xbuf    (XDR I/O buffer)
*/

/* ftype is the variable's nc_type defined in file, eg. int64
 * btype is the I/O buffer's C data type, eg. long long
 * buftype is I/O bufer's MPI data type, eg. MPI_UNSIGNED_LONG_LONG
 * apitype is data type appeared in the API names, eg. ncmpi_get_vara_longlong
 */

static int
pack_request(NC *ncp, NC_var *varp, NC_req *req, int do_vars, void *buf,
             void *xbuf, const MPI_Offset start[], const MPI_Offset count[],
             const MPI_Offset stride[], MPI_Offset fnelems, MPI_Offset bnelems,
             MPI_Offset lnelems, MPI_Offset bufcount, MPI_Datatype buftype,
             MPI_Datatype ptype, int iscontig_of_ptypes,
             int need_swap_back_buf, int use_abuf, int *reqid);

/*----< ncmpi_iput_varm() >--------------------------------------------------*/
int
ncmpi_iput_varm(int               ncid,
                int               varid,
                const MPI_Offset  start[],
                const MPI_Offset  count[],
                const MPI_Offset  stride[],
                const MPI_Offset  imap[],
                const void       *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      buftype,
                int              *reqid)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    *reqid = NC_REQ_NULL;
    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_COLL_IO, status)

    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap,
                               (void*)buf, bufcount, buftype, reqid,
                               WRITE_REQ, 0);
}

#define IPUT_VARM_TYPE(apitype, btype, buftype)                          \
int                                                                      \
ncmpi_iput_varm_##apitype(int               ncid,                        \
                          int               varid,                       \
                          const MPI_Offset  start[],                     \
                          const MPI_Offset  count[],                     \
                          const MPI_Offset  stride[],                    \
                          const MPI_Offset  imap[],                      \
                          const btype      *op,                          \
                          int              *reqid)                       \
{                                                                        \
    int         status;                                                  \
    NC         *ncp;                                                     \
    NC_var     *varp;                                                    \
    MPI_Offset  nelems;                                                  \
                                                                         \
    *reqid = NC_REQ_NULL;                                                \
    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_COLL_IO, status)      \
                                                                         \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    GET_NUM_ELEMENTS(nelems)                                             \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap,    \
                               (void*)op, nelems, buftype, reqid,        \
                               WRITE_REQ, 0);                            \
}

/*----< ncmpi_iput_varm_text() >----------------------------------------------*/
/*----< ncmpi_iput_varm_schar() >---------------------------------------------*/
/*----< ncmpi_iput_varm_uchar() >---------------------------------------------*/
/*----< ncmpi_iput_varm_short() >---------------------------------------------*/
/*----< ncmpi_iput_varm_ushort() >--------------------------------------------*/
/*----< ncmpi_iput_varm_int() >-----------------------------------------------*/
/*----< ncmpi_iput_varm_uint() >----------------------------------------------*/
/*----< ncmpi_iput_varm_long() >----------------------------------------------*/
/*----< ncmpi_iput_varm_float() >---------------------------------------------*/
/*----< ncmpi_iput_varm_double() >--------------------------------------------*/
/*----< ncmpi_iput_varm_longlong() >------------------------------------------*/
/*----< ncmpi_iput_varm_ulonglong() >-----------------------------------------*/
IPUT_VARM_TYPE(text,      char,               MPI_CHAR)
IPUT_VARM_TYPE(schar,     schar,              MPI_BYTE)
IPUT_VARM_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
IPUT_VARM_TYPE(short,     short,              MPI_SHORT)
IPUT_VARM_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
IPUT_VARM_TYPE(int,       int,                MPI_INT)
IPUT_VARM_TYPE(uint,      uint,               MPI_UNSIGNED)
IPUT_VARM_TYPE(long,      long,               MPI_LONG)
IPUT_VARM_TYPE(float,     float,              MPI_FLOAT)
IPUT_VARM_TYPE(double,    double,             MPI_DOUBLE)
IPUT_VARM_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
IPUT_VARM_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// IPUT_VARM_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

/*----< ncmpi_iget_varm() >--------------------------------------------------*/
int
ncmpi_iget_varm(int               ncid,
                int               varid,
                const MPI_Offset  start[],
                const MPI_Offset  count[],
                const MPI_Offset  stride[],
                const MPI_Offset  imap[],
                void             *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      buftype,
                int              *reqid)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    *reqid = NC_REQ_NULL;
    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_COLL_IO, status)

    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap, buf,
                               bufcount, buftype, reqid, READ_REQ, 0);
}

#define IGET_VARM_TYPE(apitype, btype, buftype)                          \
int                                                                      \
ncmpi_iget_varm_##apitype(int               ncid,                        \
                          int               varid,                       \
                          const MPI_Offset  start[],                     \
                          const MPI_Offset  count[],                     \
                          const MPI_Offset  stride[],                    \
                          const MPI_Offset  imap[],                      \
                          btype            *ip,                          \
                          int              *reqid)                       \
{                                                                        \
    int         status;                                                  \
    NC         *ncp;                                                     \
    NC_var     *varp;                                                    \
    MPI_Offset  nelems;                                                  \
                                                                         \
    *reqid = NC_REQ_NULL;                                                \
    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_COLL_IO, status)       \
                                                                         \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    GET_NUM_ELEMENTS(nelems)                                             \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap,    \
                               ip, nelems, buftype, reqid, READ_REQ, 0); \
}

/*----< ncmpi_iget_varm_text() >----------------------------------------------*/
/*----< ncmpi_iget_varm_schar() >---------------------------------------------*/
/*----< ncmpi_iget_varm_uchar() >---------------------------------------------*/
/*----< ncmpi_iget_varm_short() >---------------------------------------------*/
/*----< ncmpi_iget_varm_ushort() >--------------------------------------------*/
/*----< ncmpi_iget_varm_int() >-----------------------------------------------*/
/*----< ncmpi_iget_varm_uint() >----------------------------------------------*/
/*----< ncmpi_iget_varm_long() >----------------------------------------------*/
/*----< ncmpi_iget_varm_float() >---------------------------------------------*/
/*----< ncmpi_iget_varm_double() >--------------------------------------------*/
/*----< ncmpi_iget_varm_longlong() >------------------------------------------*/
/*----< ncmpi_iget_varm_ulonglong() >-----------------------------------------*/
IGET_VARM_TYPE(text,      char,               MPI_CHAR)
IGET_VARM_TYPE(schar,     schar,              MPI_BYTE)
IGET_VARM_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
IGET_VARM_TYPE(short,     short,              MPI_SHORT)
IGET_VARM_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
IGET_VARM_TYPE(int,       int,                MPI_INT)
IGET_VARM_TYPE(uint,      uint,               MPI_UNSIGNED)
IGET_VARM_TYPE(long,      long,               MPI_LONG)
IGET_VARM_TYPE(float,     float,              MPI_FLOAT)
IGET_VARM_TYPE(double,    double,             MPI_DOUBLE)
IGET_VARM_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
IGET_VARM_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// IGET_VARM_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

/*----< ncmpi_buffer_attach() >----------------------------------------------*/
int
ncmpi_buffer_attach(int        ncid, 
                    MPI_Offset bufsize)
{
    int status;
    NC *ncp;

    if (bufsize <= 0) return NC_ENULLBUF;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    /* check if the buffer has been prviously attached
     * note that in nc.c, the NC object is allocated with calloc, so
     * abuf should be initialized to NULL then
     */
    if (ncp->abuf != NULL) return NC_EPREVATTACHBUF;

    ncp->abuf = (NC_buf*) NCI_Malloc(sizeof(NC_buf));

    ncp->abuf->size_allocated = bufsize;
    ncp->abuf->size_used = 0;
    ncp->abuf->table_size = NC_ABUF_DEFAULT_TABLE_SIZE;
    ncp->abuf->occupy_table = (NC_buf_status*)
               NCI_Calloc(NC_ABUF_DEFAULT_TABLE_SIZE, sizeof(NC_buf_status));
    ncp->abuf->tail = 0;
    ncp->abuf->buf = NCI_Malloc(bufsize);
    return NC_NOERR;
}

/*----< ncmpi_buffer_detach() >----------------------------------------------*/
int
ncmpi_buffer_detach(int ncid)
{
    int     status;
    NC     *ncp;
    NC_req *cur_req;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    /* check if the buffer has been prviously attached */
    if (ncp->abuf == NULL) return NC_ENULLABUF;

    /* this API assumes users are responsible for no pending bput */
    cur_req = ncp->head;
    while (cur_req != NULL) { /* check if there is a pending bput */
        if (cur_req->use_abuf)
            return NC_EPENDINGBPUT;
            /* return now, so users can call wait and try detach again */
        cur_req = cur_req->next;
    }

    NCI_Free(ncp->abuf->buf);
    NCI_Free(ncp->abuf->occupy_table);
    NCI_Free(ncp->abuf);
    ncp->abuf = NULL;

    return NC_NOERR;
}
#ifdef THIS_SEEMS_OVER_DONE_IT
/*----< ncmpi_buffer_detach() >----------------------------------------------*/
/* mimic MPI_Buffer_detach()
 * Note from MPI: Even though the 'bufferptr' argument is declared as
 * 'void *', it is really the address of a void pointer.
 */
int
ncmpi_buffer_detach(int         ncid, 
                    void       *bufptr,
                    MPI_Offset *bufsize)
{
    int     status;
    NC     *ncp;
    NC_req *cur_req;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    /* check if the buffer has been prviously attached */
    if (ncp->abuf == NULL) return NC_ENULLABUF;

    /* check MPICH2 src/mpi/pt2pt/bsendutil.c for why the bufptr is void* */
    *(void **)bufptr = ncp->abuf->buf;
    *bufsize         = ncp->abuf->size_allocated;

    /* this API assumes users are respobsible for no pending bput when called */
    cur_req = ncp->head;
    while (cur_req != NULL) { /* check if there is a pending bput */
        if (cur_req->use_abuf)
            return NC_EPENDINGBPUT;
        cur_req = cur_req->next;
    }

    NCI_Free(ncp->abuf->occupy_table);
    NCI_Free(ncp->abuf);
    ncp->abuf = NULL;

    return NC_NOERR;
}
#endif

/*----< ncmpii_abuf_malloc() >-----------------------------------------------*/
static int
ncmpii_abuf_malloc(NC *ncp, MPI_Offset nbytes, void **buf)
{
    /* extend the table size if more entries are needed */
    if (ncp->abuf->tail + 1 == ncp->abuf->table_size) {
        ncp->abuf->table_size += NC_ABUF_DEFAULT_TABLE_SIZE;
        ncp->abuf->occupy_table = (NC_buf_status*)
                   NCI_Realloc(ncp->abuf->occupy_table,
                               ncp->abuf->table_size * sizeof(NC_buf_status));
    }
    /* mark the new entry is used and store the requested buffer size */
    ncp->abuf->occupy_table[ncp->abuf->tail].is_used  = 1;
    ncp->abuf->occupy_table[ncp->abuf->tail].req_size = nbytes;

    *buf = (char*)ncp->abuf->buf + ncp->abuf->size_used;
    ncp->abuf->size_used += nbytes;
    ncp->abuf->tail++;

    return NC_NOERR;
}

/*----< ncmpii_igetput_varm() >----------------------------------------------*/
int
ncmpii_igetput_varm(NC               *ncp,
                    NC_var           *varp,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    const MPI_Offset  imap[],
                    void             *buf,
                    MPI_Offset        bufcount,
                    MPI_Datatype      buftype,
                    int              *reqid,
                    int               rw_flag,
                    int               use_abuf) /* if use attached buffer */
{
    void *xbuf=NULL, *cbuf=NULL, *lbuf=NULL;
    int err, status, warning; /* err is for API abort and status is not */
    int el_size, iscontig_of_ptypes, do_vars, isderived;
    int need_convert, need_swap, need_swap_back_buf;
    int i, dim=0, imap_contig_blocklen;
    MPI_Offset fnelems, bnelems, lnelems, nbytes;
    MPI_Datatype ptype, imaptype=MPI_DATATYPE_NULL;
    NC_req *req;

    /* "API error" will abort this API call, but not the entire program */
    err = status = warning = NC_NOERR;

    if (varp->ndims > 0) {
        assert(start != NULL);
        assert(count != NULL);
    }

    do_vars = 0;

    if (varp->ndims == 0)
        /* reduced to scalar var, only one value at one fixed place */
        do_vars = 1;

    if (imap == NULL) /* no mapping, same as vars */
        do_vars = 1;
    else {
        imap_contig_blocklen = 1;
        dim = varp->ndims;
        /* test each dim's contiguity until the 1st non-contiguous dim is
           reached */
        while ( --dim >= 0 && imap_contig_blocklen == imap[dim] ) {
            if (count[dim] < 0)
                return NC_ENEGATIVECNT;
            imap_contig_blocklen *= count[dim];
        }
        if (dim == -1) /* imap is a contiguous layout */
            do_vars = 1;
    }
    /* dim is the first dimension (C order, eg. ZYX) that has non-contiguous
       imap and it will be used only when do_vars == 1
     */

    /* find the ptype (primitive MPI data type) from buftype
     * el_size is the element size of ptype
     * bnelems is the total number of ptype elements in the I/O buffer, buf
     * fnelems is the number of nc variable elements in nc_type
     * nbytes is the amount of read/write in bytes
     */
    err = ncmpii_dtype_decode(buftype, &ptype, &el_size, &bnelems,
                              &isderived, &iscontig_of_ptypes);
    /* bnelems now is the number of ptype in a buftype */
    if (err != NC_NOERR) goto err_check;

    err = NCMPII_ECHAR(varp->type, ptype);
    if (err != NC_NOERR) goto err_check;

    CHECK_NELEMS(varp, fnelems, count, bnelems, bufcount, nbytes, err)
    /* bnelems now is the number of ptype in the whole buf */
    /* warning is set in CHECK_NELEMS() */

    /* for bput call, check if the remaining buffer space is sufficient
       to accommodate this request */
    if (rw_flag == WRITE_REQ && use_abuf &&
        ncp->abuf->size_allocated - ncp->abuf->size_used < nbytes)
        return NC_EINSUFFBUF;

    need_convert  = ncmpii_need_convert(varp->type, ptype);
    need_swap     = ncmpii_need_swap(varp->type, ptype);
    need_swap_back_buf = 0;

err_check:
    if (err != NC_NOERR) return err;

    if (bnelems == 0) {
        /* zero-length request, mark this as a NULL request */
        *reqid = NC_REQ_NULL;
        return NCcoordck(ncp, varp, start);
    }

/* Here is the pseudo code description on buffer packing
    if (iscontig_of_ptypes)
        lbuf = buf
    else
        lbuf = malloc
        pack buf -> lbuf

    if do_vars
        build imaptype
        cbuf = malloc
        pack lbuf -> cbuf
        if lbuf != buf, free lbuf
    else
        cbuf = lbuf
    lbuf = NULL

    if need convert
        if use_abuf
            xbuf = attach_buf_malloc
        else
            xbuf = malloc
        convert cbuf -> xbuf
        if cbuf != buf, free cbuf
    else
        if use_abuf
            xbuf = attach_buf_malloc
            memcpy(xbuf, cbuf)
            if cbuf != buf, free cbuf
        else
            xbuf = cbuf
        if need swap
            swap xbuf
    cbuf = NULL
*/

    if (!do_vars) {
        /* construct a derived data type, imaptype, based on imap[], and use
         * it to pack lbuf to cbuf.
         */
        MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
                        ptype, &imaptype);
        MPI_Type_commit(&imaptype);
        for (i=dim, i--; i>=0; i--) {
            MPI_Datatype tmptype;
            if (count[i] < 0)
                return ((warning != NC_NOERR) ? warning : NC_ENEGATIVECNT);
#if (MPI_VERSION < 2)
            MPI_Type_hvector(count[i], 1, imap[i]*el_size, imaptype, &tmptype);
#else
            MPI_Type_create_hvector(count[i], 1, imap[i]*el_size, imaptype,
                                    &tmptype);
#endif
            MPI_Type_free(&imaptype);
            MPI_Type_commit(&tmptype);
            imaptype = tmptype;
        }
    }

    lbuf = NULL;
    cbuf = NULL;
    xbuf = NULL;

    if (rw_flag == WRITE_REQ) {

        /* Step 1: pack buf into a contiguous buffer, lbuf */
        if (iscontig_of_ptypes) { /* buf is contiguous */
            lnelems = bnelems / bufcount;
            lbuf = buf;
        }
        else { /* pack buf into lbuf, a contiguous buffer, based on buftype */
            lnelems = bnelems;
            lbuf = NCI_Malloc(lnelems*el_size);
            status = ncmpii_data_repack(buf, bufcount, buftype,
                                        lbuf, lnelems, ptype);
            if (status != NC_NOERR) {
                NCI_Free(lbuf);
                return ((warning != NC_NOERR) ? warning : status);
            }
        }

        /* Step 2: pack lbuf to cbuf if imap is non-contiguos */
        if (do_vars) { /* reuse lbuf */
            cbuf = lbuf;
        }
        else { /* a true varm case, pack lbuf to cbuf based on imap */
            bnelems = imap_contig_blocklen * count[dim];
            for (dim--; dim>=0; dim--)
                bnelems *= count[dim];

            cbuf = NCI_Malloc(bnelems*el_size);

            /* pack lbuf to cbuf based on imaptype */
            status = ncmpii_data_repack(lbuf, 1, imaptype,
                                        cbuf, bnelems, ptype);
            MPI_Type_free(&imaptype);
            imaptype = MPI_DATATYPE_NULL;

            /* for write case, lbuf is no longer needed */
            if (lbuf != buf) NCI_Free(lbuf);

            if (status != NC_NOERR) {
                NCI_Free(cbuf);
                return ((warning != NC_NOERR) ? warning : status);
            }
        }
        lbuf = NULL; /* no longer need lbuf */

        /* Step 3: pack cbuf to xbuf and xbuf will be used to write to file */
        if (need_convert) { /* user buf type != nc var type defined in file */
            if (use_abuf) { /* use attached buffer */
                status = ncmpii_abuf_malloc(ncp, nbytes, &xbuf);
                if (status != NC_NOERR) {
                    if (cbuf != NULL && cbuf != buf) NCI_Free(cbuf);
                    return ((warning != NC_NOERR) ? warning : status);
                }
            }
            else
                xbuf = NCI_Malloc(nbytes);

            /* datatype conversion + byte-swap from cbuf to xbuf */
            DATATYPE_PUT_CONVERT(varp->type, xbuf, cbuf, bnelems, ptype, err)
            /* retain the first error status */
            if (status == NC_NOERR) status = err;
        }
        else {  /* cbuf == xbuf */
            if (use_abuf) { /* use attached buffer */
                status = ncmpii_abuf_malloc(ncp, nbytes, &xbuf);
                if (status != NC_NOERR) {
                    if (cbuf != NULL && cbuf != buf) NCI_Free(cbuf);
                    return ((warning != NC_NOERR) ? warning : status);
                }
                memcpy(xbuf, cbuf, nbytes);
            } else {
                xbuf = cbuf;
            }
            if (need_swap) {
                /* perform array in-place byte swap on xbuf */
                ncmpii_in_swapn(xbuf, fnelems, ncmpix_len_nctype(varp->type));
                if (xbuf == buf)
                    need_swap_back_buf = 1;
                    /* user buf needs to be swapped back to its original
                     * contents as now buf == cbuf == xbuf */
            }
        }
        /* cbuf is no longer needed */
        if (cbuf != buf && cbuf != xbuf) NCI_Free(cbuf);
        cbuf = NULL;
    }
    else { /* rw_flag == READ_REQ */
        /* Read is done at wait call, need lnelems and bnelems to reverse the
         * steps as done in write case */
        if (iscontig_of_ptypes)
            lnelems = bnelems / bufcount;
        else
            lnelems = bnelems;

        if (!do_vars) {
            bnelems = imap_contig_blocklen * count[dim];
            for (dim--; dim>=0; dim--)
                bnelems *= count[dim];
        }
        if (iscontig_of_ptypes && do_vars && !need_convert)
            xbuf = buf;  /* there is no buffered read (bget_var, etc.) */
        else
            xbuf = NCI_Malloc(nbytes);
    }

    /* allocate a new request object to store the write info */
    req = (NC_req*) NCI_Malloc(sizeof(NC_req));

    req->is_imap  = 0;
    req->imaptype = imaptype;
    req->rw_flag  = rw_flag;

    if (!do_vars)
        req->is_imap = 1;

    pack_request(ncp, varp, req, do_vars, buf, xbuf, start, count, stride,
                 fnelems, bnelems, lnelems, bufcount, buftype, ptype,
                 iscontig_of_ptypes, need_swap_back_buf, use_abuf,
                 reqid);

    return ((warning != NC_NOERR) ? warning : status);
}

/*----< pack_request() >------------------------------------------------------*/
static int
pack_request(NC               *ncp,
             NC_var           *varp,
             NC_req           *req,
             int               do_vars,
             void             *buf,
             void             *xbuf,
             const MPI_Offset  start[],
             const MPI_Offset  count[],
             const MPI_Offset  stride[],
             MPI_Offset        fnelems,
             MPI_Offset        bnelems,
             MPI_Offset        lnelems,
             MPI_Offset        bufcount,
             MPI_Datatype      buftype,
             MPI_Datatype      ptype,
             int               iscontig_of_ptypes,
             int               need_swap_back_buf,
             int               use_abuf,
             int              *reqid)
{
    int     i, j;
    NC_req *subreqs;

    req->varp     = varp;
    req->ndims    = varp->ndims;
    req->start    = (MPI_Offset*) NCI_Malloc(2*varp->ndims*sizeof(MPI_Offset));
    req->count    = req->start + varp->ndims;
    req->buf      = buf;
    req->xbuf     = xbuf;
    req->fnelems  = fnelems;
    req->bnelems  = bnelems;
    req->lnelems  = lnelems; /* used only for iget_varm case */
    req->buftype  = buftype;
    req->bufcount = bufcount;
    req->ptype    = ptype;   /* MPI element datatype for the I/O buffer */
    req->next     = NULL;
    req->subreqs     = NULL;
    req->num_subreqs = 0;
    req->iscontig_of_ptypes = iscontig_of_ptypes;
    req->need_swap_back_buf = need_swap_back_buf;
    req->use_abuf           = use_abuf;

    if (stride != NULL)
        req->stride = (MPI_Offset*) NCI_Malloc(varp->ndims*sizeof(MPI_Offset));
    else
        req->stride = NULL;

    for (i=0; i<varp->ndims; i++) {
        req->start[i] = start[i];
        req->count[i] = count[i];
        if (stride != NULL)
            req->stride[i] = stride[i];
    }
    /* get the starting file offset for this request */
    ncmpii_get_offset(ncp, varp, start, NULL, NULL, &req->offset_start);

    /* get the ending file offset for this request */
    ncmpii_get_offset(ncp, varp, start, count, stride, &req->offset_end);
    req->offset_end += varp->xsz - 1;

    /* check if this is a record varaible. if yes, split the request into
       subrequests, one iput request for a record access. Hereandafter,
       treat each request as a non-record variable request */

    /* check if this access is within one record, if yes, no need to create
       subrequests */
    if (IS_RECVAR(varp) && req->count[0] > 1) {
        MPI_Offset rec_bufcount = 1;
        for (i=1; i<varp->ndims; i++)
            rec_bufcount *= req->count[i];

        subreqs = (NC_req*) NCI_Malloc(req->count[0]*sizeof(NC_req));
        for (i=0; i<req->count[0]; i++) {
            MPI_Offset span;
            subreqs[i] = *req; /* inherit most attributes from req */

            /* each sub-request contains <= one record size */
            subreqs[i].start = (MPI_Offset*) NCI_Malloc(2*varp->ndims*sizeof(MPI_Offset));
            subreqs[i].count = subreqs[i].start + varp->ndims;
            if (stride != NULL) {
                subreqs[i].stride = (MPI_Offset*) NCI_Malloc(varp->ndims*sizeof(MPI_Offset));
                subreqs[i].start[0] = req->start[0] + stride[0] * i;
                subreqs[i].stride[0] = req->stride[0];
            } else {
                subreqs[i].stride = NULL;
                subreqs[i].start[0] = req->start[0] + i;
            }

            subreqs[i].count[0] = 1;
            subreqs[i].fnelems = 1;
            for (j=1; j<varp->ndims; j++) {
                subreqs[i].start[j]  = req->start[j];
                subreqs[i].count[j]  = req->count[j];
                subreqs[i].fnelems  *= subreqs[i].count[j];
                if (stride != NULL)
                    subreqs[i].stride[j] = req->stride[j];
            }
            ncmpii_get_offset(ncp, varp, subreqs[i].start, NULL, NULL,
                              &subreqs[i].offset_start);
            ncmpii_get_offset(ncp, varp, subreqs[i].start,
                              subreqs[i].count, subreqs[i].stride,
                              &subreqs[i].offset_end);
            subreqs[i].offset_end += varp->xsz - 1;

            span                = i*rec_bufcount*varp->xsz;
            subreqs[i].buf      = (char*)(req->buf)  + span;
            /* xbuf cannot be NULL    assert(req->xbuf != NULL); */
            subreqs[i].xbuf     = (char*)(req->xbuf) + span;
            subreqs[i].bufcount = rec_bufcount;
        }
        req->num_subreqs = req->count[0];
        req->subreqs     = subreqs;
    }

    /* add the new request to the internal request array (or linked list) */
    if (ncp->head == NULL) {
        req->id   = 0;
        ncp->head = req;
        ncp->tail = ncp->head;
    }
    else { /* add to the tail */
        req->id = ncp->tail->id + 1;
        ncp->tail->next = req;
        ncp->tail = ncp->tail->next;
    }
    ncp->tail->next = NULL;

    /* return the request ID */
    *reqid = ncp->tail->id;

    return NC_NOERR;
}

/*----< ncmpi_inq_buffer_usage() >--------------------------------------------*/
int
ncmpi_inq_buffer_usage(int         ncid, 
                       MPI_Offset *usage) /* in bytes */
{
    int  status;
    NC  *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    /* check if the buffer has been prviously attached */
    if (ncp->abuf == NULL) return NC_ENULLABUF;

    /* return the current usage in bytes */
    *usage = ncp->abuf->size_used;

    return NC_NOERR;
}

