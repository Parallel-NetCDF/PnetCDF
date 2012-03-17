/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

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
pack_request(NC *ncp, NC_var *varp, NC_req *req, int do_vars,
             int iscontig_of_ptypes, void *buf, void *cbuf, void *xbuf,
             const MPI_Offset start[], const MPI_Offset count[],
             const MPI_Offset stride[], MPI_Offset fnelems, MPI_Offset bnelems,
             MPI_Offset lnelems, MPI_Offset bufcount, MPI_Datatype buftype,
             MPI_Datatype ptype, int *reqid);

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
    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap,
                               (void*)buf, bufcount, buftype, reqid,
                               WRITE_REQ);
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
    CHECK_NCID                                                           \
    CHECK_WRITE_PERMISSION                                               \
    if (NC_indef(ncp)) return NC_EINDEFINE;                              \
    CHECK_VARID(varid, varp)                                             \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    GET_NUM_ELEMENTS                                                     \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap,    \
                               (void*)op, nelems, buftype, reqid,        \
                               WRITE_REQ);                               \
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
    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap, buf,
                               bufcount, buftype, reqid, READ_REQ);
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
    CHECK_NCID                                                           \
    if (NC_indef(ncp)) return NC_EINDEFINE;                              \
    CHECK_VARID(varid, varp)                                             \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    GET_NUM_ELEMENTS                                                     \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap,    \
                               ip, nelems, buftype, reqid, READ_REQ);    \
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
                    int               rw_flag)
{
    void *xbuf=NULL, *cbuf=NULL, *lbuf=NULL;
    int err, status, warning; /* err is for API abort and status is not */
    int el_size, iscontig_of_ptypes, do_vars, isderived;
    int dim=0, imap_contig_blocklen;
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

    CHECK_NELEMS(varp, fnelems, count, bnelems, bufcount, nbytes)
    /* bnelems now is the number of ptype in the whole buf */
    /* warning is set in CHECK_NELEMS() */

err_check:
    if (err != NC_NOERR) return err;

    if (bnelems == 0)
        return NCcoordck(ncp, varp, start);

    if (do_vars) {
        if (!iscontig_of_ptypes) {
            cbuf = NCI_Malloc(bnelems * el_size);

            if (rw_flag == WRITE_REQ) {
                /* pack buf into cbuf, a contiguous buffer, based on buftype */
                status = ncmpii_data_repack(buf, bufcount, buftype,
                                            cbuf, bnelems, ptype);
                if (status != NC_NOERR) {
                    NCI_Free(cbuf);
                    return ((warning != NC_NOERR) ? warning : status);
                }
            }
        } else {
            cbuf = buf;
        }
    } else { /* of if (do_vars) */
        MPI_Datatype tmptype;
        if (!iscontig_of_ptypes) {
            lnelems = bnelems;
            lbuf = NCI_Malloc(lnelems*el_size);

            if (rw_flag == WRITE_REQ) {
                /* pack buf into lbuf, a contiguous buffer, based on buftype */
                status = ncmpii_data_repack(buf, bufcount, buftype,
                                            lbuf, lnelems, ptype);
                if (status != NC_NOERR) {
                    NCI_Free(lbuf);
                    return ((warning != NC_NOERR) ? warning : status);
                }
            }
        } else {
            lnelems = bnelems / bufcount;
            lbuf = buf;
        }
        /* now lbuf is a contiguous buffer containing the data from buf */

        /* construct a derived data type based on imap[], and use it to
         * pack lbuf to another contiguous buffer, cbuf. Then we can use cbuf
         * to call getput_vars
         */
        MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
                        ptype, &imaptype);
        MPI_Type_commit(&imaptype);
        bnelems = imap_contig_blocklen * count[dim];
        for (dim--; dim>=0; dim--) {
            if (count[dim] < 0) {
                if (!iscontig_of_ptypes)
                    NCI_Free(lbuf);
                return ((warning != NC_NOERR) ? warning : NC_ENEGATIVECNT);
            }
#if (MPI_VERSION < 2)
            MPI_Type_hvector(count[dim], 1, imap[dim]*el_size, imaptype,
                             &tmptype);
#else
            MPI_Type_create_hvector(count[dim], 1, imap[dim]*el_size,
                                    imaptype, &tmptype);
#endif
            MPI_Type_free(&imaptype);
            MPI_Type_commit(&tmptype);
            imaptype = tmptype;
            bnelems *= count[dim];
        }
        cbuf = NCI_Malloc(bnelems*el_size);

        if (rw_flag == WRITE_REQ) {
            /* layout lbuf to cbuf based on imap */
            status = ncmpii_data_repack(lbuf, 1, imaptype,
                                        cbuf, bnelems, ptype);
            /* for write case, lbuf is no longer needed */
            if (!iscontig_of_ptypes)
                NCI_Free(lbuf);

            if (status != NC_NOERR) {
                NCI_Free(cbuf);
                return ((warning != NC_NOERR) ? warning : status);
            }
            lbuf = NULL;
            MPI_Type_free(&imaptype);

            /* now it is the same as using cbuf to call
               ncmpi_put_vars(ncid, varid, start, count, stride,
                              cbuf, bnelems, ptype);
             */
        }
    }
    /* up to this point, cbuf contains all the write data in a contiguous
     * space and it can be == buf or created from malloc
     */

    if ( ncmpii_need_convert(varp->type, ptype) ) {
        /* allocate new buffer for data type conversion */

        xbuf = NCI_Malloc(nbytes);

        if (rw_flag == WRITE_REQ)
            /* automatic numeric datatype conversion + swap if necessary */
            DATATYPE_PUT_CONVERT(varp->type, xbuf, cbuf, bnelems, ptype)
    } else {
        xbuf = cbuf;

        if (rw_flag == WRITE_REQ && ncmpii_need_swap(varp->type, ptype))
            /* perform array in-place byte swap */
            ncmpii_in_swapn(xbuf, fnelems, ncmpix_len_nctype(varp->type));
    }

    /* allocate a new request object to store the write info */
    req = (NC_req*) NCI_Malloc(sizeof(NC_req));

    req->is_imap  = 0;
    req->imaptype = imaptype;
    req->rw_flag  = rw_flag;
    if (rw_flag == WRITE_REQ)
        req->lbuf = NULL;  /* for write, lbuf has been freed */
    else
        req->lbuf = lbuf;  /* for read, lbuf will be used to imap back to buf */
    if (!do_vars)
        req->is_imap = 1;

    pack_request(ncp, varp, req, do_vars, iscontig_of_ptypes,
                 buf, cbuf, xbuf, start, count, stride, fnelems,
                 bnelems, lnelems, bufcount, buftype, ptype, reqid);

    return ((warning != NC_NOERR) ? warning : status);
}

/*----< pack_request() >------------------------------------------------------*/
static int
pack_request(NC               *ncp,
             NC_var           *varp,
             NC_req           *req,
             int               do_vars,
             int               iscontig_of_ptypes,
             void             *buf,
             void             *cbuf,
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
             int              *reqid)
{
    int     i, j;
    NC_req *subreqs;

    req->varp     = varp;
    req->ndims    = varp->ndims;
    req->start    = (MPI_Offset*) NCI_Malloc(2*varp->ndims*sizeof(MPI_Offset));
    req->count    = req->start + varp->ndims;
    req->buf      = buf;
    req->cbuf     = cbuf;
    req->xbuf     = xbuf;
    req->fnelems  = fnelems;
    req->bnelems  = bnelems;
    req->lnelems  = lnelems; /* used only for iget_varm case */
    req->buftype  = buftype;
    req->bufcount = bufcount;
    req->ptype    = ptype;
    req->next     = NULL;
    req->subreqs     = NULL;
    req->num_subreqs = 0;
    req->iscontig_of_ptypes = iscontig_of_ptypes;

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
            subreqs[i].cbuf     = (char*)(req->cbuf) + span;
            subreqs[i].xbuf     = (char*)(req->xbuf) + span;
            subreqs[i].lbuf     = NULL;
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

