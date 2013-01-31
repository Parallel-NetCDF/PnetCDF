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

/* Prototypes for functions used only in this file */
static int ncmpii_wait_getput(NC *ncp, int num_reqs, NC_req *req_head,
                              int rw_flag, int io_method);

static int ncmpii_mgetput(NC *ncp, int num_reqs, NC_var *varps[],
                          MPI_Offset *starts[], MPI_Offset *counts[],
                          MPI_Offset *strides[], void *bufs[],
                          MPI_Offset nbytes[], int statuses[], int rw_flag,
                          int io_method);

/*----< ncmpii_abuf_free() >--------------------------------------------------*/
static int
ncmpii_abuf_free(NC *ncp, void *buf)
{
    int i;
    /* find the index in the occupy table for this buffer
     * note that abuf->tail points to the first free entry */
    MPI_Aint buf_addr, abuf_addr, dist, buf_offset = 0;
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(buf, &buf_addr);
    MPI_Get_address(ncp->abuf->buf, &abuf_addr);
#else 
    MPI_Address(buf, &buf_addr);
    MPI_Address(ncp->abuf->buf, &abuf_addr);
#endif
    dist = buf_addr - abuf_addr;
    for (i=0; i<ncp->abuf->tail; i++) {
        if (dist <= buf_offset)
            break;
        buf_offset += ncp->abuf->occupy_table[i].req_size;
    }
    assert(i < ncp->abuf->tail);
    ncp->abuf->occupy_table[i].is_used = 0; /* set to free */

    /* coalesce the freed entries backwardly from the tail */
    while (i >= 0 && i+1 == ncp->abuf->tail) {
        ncp->abuf->size_used -= ncp->abuf->occupy_table[i].req_size;
        ncp->abuf->tail--;
        if (i > 0 && ncp->abuf->occupy_table[i-1].is_used == 0)
            i--;
    }
    return NC_NOERR;
}

#define FREE_REQUEST(req) {                                                   \
    if (req->use_abuf)                                                        \
        ncmpii_abuf_free(ncp, req->xbuf);                                     \
    else if (req->xbuf != NULL && req->xbuf != req->buf)                      \
        NCI_Free(req->xbuf);                                                  \
    req->xbuf = NULL;                                                         \
                                                                              \
    for (j=0; j<req->num_subreqs; j++)                                        \
        NCI_Free(req->subreqs[j].start);                                      \
    if (req->num_subreqs > 0)                                                 \
        NCI_Free(req->subreqs);                                               \
    NCI_Free(req->start);                                                     \
    if (req->stride != NULL)                                                  \
	NCI_Free(req->stride);                                                \
}

/*----< ncmpi_cancel() >------------------------------------------------------*/
int
ncmpi_cancel(int  ncid,
             int  num_req,
             int *req_ids,  /* [num_req] */
             int *statuses) /* [num_req], can be NULL if users don't care */
{
    int status;
    NC *ncp;

    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;

    return ncmpii_cancel(ncp, num_req, req_ids, statuses);
}

/*----< ncmpii_cancel() >-----------------------------------------------------*/
int
ncmpii_cancel(NC  *ncp,
              int  num_req,
              int *req_ids,  /* [num_req] */
              int *statuses) /* [num_req], can be NULL if users don't care */
{
    int i, j;
    NC_req *pre_req, *cur_req;

    if (num_req == 0)
        return NC_NOERR;

    /* collect the requests from the linked list */
    for (i=0; i<num_req; i++) {
        if (statuses != NULL)
            statuses[i] = NC_NOERR;
        if (req_ids[i] == NC_REQ_NULL)
            continue;

        if (ncp->head == NULL) {
            printf("Error: ncmpii_cancel() NULL request queue NULL at line %d of %s\n", __LINE__, __FILE__);
            if (statuses != NULL)
                statuses[i] = NC_EINVAL_REQUEST;
            return NC_EINVAL_REQUEST;
        }

        pre_req = NULL;
        cur_req = ncp->head;
        while (cur_req->id != req_ids[i]) { /* find from the linked list */
            pre_req = cur_req;
            cur_req = cur_req->next;
            if (cur_req == NULL) {
                printf("Error: ncmpii_cancel() no such request ID = %d\n", req_ids[i]);
                if (statuses != NULL)
                    statuses[i] = NC_EINVAL_REQUEST;
                return NC_EINVAL_REQUEST;
            }
        }
        /* found it */
        if (cur_req == ncp->head) /* move cur_req to next */
            ncp->head = cur_req->next;
        else /* move pre_req and cur_req to next */
            pre_req->next = cur_req->next;

        FREE_REQUEST(cur_req)
        NCI_Free(cur_req);
    }

    /* make sure ncp->tail pointing to the tail */
    ncp->tail = ncp->head;
    while (ncp->tail != NULL && ncp->tail->next != NULL)
        ncp->tail = ncp->tail->next;

    return NC_NOERR;
}

/*----< ncmpi_wait() >--------------------------------------------------------*/
/* ncmpi_wait() is an independent call */
int
ncmpi_wait(int ncid,
           int num_reqs,
           int *req_ids,  /* [num_reqs] */
           int *statuses) /* [num_reqs] */
{
#ifndef ENABLE_NONBLOCKING
    int i;
#endif
    int err, status=NC_NOERR;
    NC  *ncp;

    if (num_reqs == 0) return NC_NOERR;

    CHECK_NCID
#ifdef ENABLE_NONBLOCKING
    return ncmpii_wait(ncp, INDEP_IO, num_reqs, req_ids, statuses);
#else
    for (i=0; i<num_reqs; i++) { /* serve one request at a time */
        err = ncmpii_wait(ncp, INDEP_IO, 1, &req_ids[i], &statuses[i]);
        if (status == NC_NOERR) status = err;
    }
    return status; /* return the first error status, if there is any */
#endif
}

/*----< ncmpi_wait_all() >----------------------------------------------------*/
/* ncmpi_wait_all() is a collective call */
int
ncmpi_wait_all(int  ncid,
               int  num_reqs, /* number of requests */
               int *req_ids,  /* [num_reqs] */
               int *statuses) /* [num_reqs] */
{
#ifndef ENABLE_NONBLOCKING
    int  i;
#endif
    int err, status=NC_NOERR;
    NC  *ncp;

    /* the following line CANNOT be added, because ncmpi_wait_all() is a
     * collective call, all processes must participate some MPI collective
     * operations used later on.
     */
    /* if (num_reqs == 0) return NC_NOERR; */

    CHECK_NCID
#ifdef ENABLE_NONBLOCKING
    return ncmpii_wait(ncp, COLL_IO, num_reqs, req_ids, statuses);
#else
    /* This API is collective, so make it illegal in indep mode. This also
       ensures the program will returns back to collective mode. */
    if (NC_indep(ncp)) return NC_EINDEP;

    /* must enter independent mode, as num_reqs may be different among
       processes */
    err = ncmpi_begin_indep_data(ncid);
    if (status == NC_NOERR) status = err;

    for (i=0; i<num_reqs; i++) { /* serve one request at a time */
        err = ncmpii_wait(ncp, INDEP_IO, 1, &req_ids[i], &statuses[i]);
        if (status == NC_NOERR) status = err;
    }

    /* return to collective data mode */
    err = ncmpi_end_indep_data(ncid);
    if (status == NC_NOERR) status = err;

    return status; /* return the first error status, if there is any */
#endif
}

/*----< ncmpii_mset_fileview() >----------------------------------------------*/
static int
ncmpii_mset_fileview(MPI_File    fh,
                     NC         *ncp, 
                     int         ntimes, 
                     NC_var     *varp[],     /* [ntimes] */
                     MPI_Offset *starts[],   /* [ntimes] */
                     MPI_Offset *counts[],   /* [ntimes] */
                     MPI_Offset *strides[],  /* [ntimes] */
                     int         statuses[], /* [ntimes] */
                     int         rw_flag) 
{
    int i, status=NC_NOERR, mpireturn, mpi_err=NC_NOERR, *blocklens;
    MPI_Datatype full_filetype, *filetypes;
    MPI_Offset *offsets;

    if (ntimes <= 0) { /* participate collective call */
        mpireturn = MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                      MPI_INFO_NULL);
        CHECK_MPI_ERROR(mpireturn, "MPI_File_set_view", NC_EFILE);
        return mpi_err;
    }

    blocklens = (int*)          NCI_Malloc(ntimes * sizeof(int));
    offsets   = (MPI_Offset*)   NCI_Malloc(ntimes * sizeof(MPI_Offset));
    filetypes = (MPI_Datatype*) NCI_Malloc(ntimes * sizeof(MPI_Datatype));

    /* create a filetype for each variable */
    for (i=0; i<ntimes; i++) {
        filetypes[i] = MPI_BYTE;
        statuses[i] = ncmpii_vars_create_filetype(ncp,
                                                  varp[i],
                                                  starts[i],
                                                  counts[i],
                                                  strides[i],
                                                  rw_flag,
                                                  &offsets[i],
                                                  &filetypes[i]);

        if (status == NC_NOERR)
            status = statuses[i]; /* report the first error */

        blocklens[i] = 1;
        if (filetypes[i] == MPI_BYTE) { /* file type is contiguous in file */
            int j;
            blocklens[i] = varp[i]->xsz;
            for (j=0; j<varp[i]->ndims; j++)
                blocklens[i] *= counts[i][j];
        }
    }

    if (status != NC_NOERR) {
        /* even if error occurs, we still msut participate the collective
           call to MPI_File_set_view() */
        MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    }
    else if (ntimes == 1) { /* no multiple filetypes to combine */
        full_filetype = filetypes[0];
        /* filetypes[0] has been committed already */

        mpireturn = MPI_File_set_view(fh, offsets[0], MPI_BYTE, full_filetype,
                                      "native", MPI_INFO_NULL);
        CHECK_MPI_ERROR(mpireturn, "MPI_File_set_view", NC_EFILE)
    }
    else {
        /* on most 32 bit systems, MPI_Aint and MPI_Offset are different sizes.
         * Possible that on those platforms some of the beginning offsets of
         * these variables in the dataset won't fit into the aint used by
         * MPI_Type_create_struct.  Minor optimization: we don't need to do any
         * of this if MPI_Aint and MPI_Offset are the same size  */
        MPI_Aint *addrs;

        if (sizeof(MPI_Offset) != sizeof(MPI_Aint)) {
            addrs = (MPI_Aint *) NCI_Malloc(ntimes * sizeof(MPI_Aint));
            for (i=0; i< ntimes; i++) {
                addrs[i] = offsets[i];
                if (addrs[i] != offsets[i]) return NC_EAINT_TOO_SMALL;
            }
        } else {
            addrs = (MPI_Aint*) offsets; /* cast ok: types same size */
        }
#if (MPI_VERSION < 2)
        MPI_Type_struct(ntimes, blocklens, addrs, filetypes, &full_filetype);
#else
        MPI_Type_create_struct(ntimes, blocklens, addrs, filetypes, &full_filetype);
#endif
        MPI_Type_commit(&full_filetype);

        mpireturn = MPI_File_set_view(fh, 0, MPI_BYTE, full_filetype, "native",
                                      MPI_INFO_NULL);
        CHECK_MPI_ERROR(mpireturn, "MPI_File_set_view", NC_EFILE)
        MPI_Type_free(&full_filetype);

        if (sizeof(MPI_Offset) != sizeof(MPI_Aint))
            NCI_Free(addrs);
    }

    for (i=0; i<ntimes; i++) {
        if (filetypes[i] != MPI_BYTE)
            MPI_Type_free(&filetypes[i]);
    }

    NCI_Free(filetypes);
    NCI_Free(offsets);
    NCI_Free(blocklens);

    /* make NC error higher priority than MPI error */
    return (status != NC_NOERR) ? status : mpi_err;
}

/*----< ncmpii_wait() >-------------------------------------------------------*/
/* The buffer management flow is described below. The wait side starts from
   the I/O step, i.e. step 5

   for put_varm:
     1. pack buf to lbuf based on buftype
     2. create imap_type based on imap
     3. pack lbuf to cbuf based on imap_type
     4. type convert and byte swap cbuf to xbuf
     5. write from xbuf
     6. byte swap the buf to its original, if it is swapped
     7. free up temp buffers, cbuf and xbuf if they were allocated separately

   for get_varm:
     1. allocate lbuf
     2. create imap_type based on imap
     3. allocate cbuf
     4. allocate xbuf
     5. read to xbuf
     6. type convert and byte swap xbuf to cbuf
     7. unpack cbuf to lbuf based on imap_type
     8. unpack lbuf to buf based on buftype
     9. free up temp buffers, cbuf and xbuf if they were allocated separately
 */
int
ncmpii_wait(NC  *ncp,
            int  io_method, /* COLL_IO or INDEP_IO */
            int  num_reqs,  /* number of requests */
            int *req_ids,   /* [num_reqs] */
            int *statuses)  /* [num_reqs] */
{
    int i, j, err, status=NC_NOERR;
    int do_read, do_write, num_w_reqs, num_r_reqs;
    NC_req *pre_req, *cur_req;
    NC_req *w_req_head, *w_req_tail, *r_req_head, *r_req_tail;

    if (NC_indef(ncp)) return NC_EINDEFINE;
 
    /* check to see that the desired MPI file handle is opened */
    if (io_method == COLL_IO)
        CHECK_COLLECTIVE_FH
    else
        CHECK_INDEP_FH

    /* Note: 1) it is illegal num_reqs is larger than the linked list size
             2) request ids must be distinct
     */
    j = 0;
    for (i=0; i<num_reqs; i++) {
        statuses[i] = NC_NOERR;
        if (req_ids[i] == NC_REQ_NULL) /* skip zero-size request */
            continue;
        j++;
    }
    /* j now is the true number of non-zero length requests */

    if (io_method == INDEP_IO && j == 0)
        return NC_NOERR;
    /* For collective APIs, even though some processes may have zero-length
       requests, they must still participate the collective call. Hence,
       only independent APIs stop here if request is of zero length.
     */

    w_req_head = w_req_tail = NULL;
    r_req_head = r_req_tail = NULL;
    num_w_reqs = num_r_reqs = 0;

    if (j > 0) /* j is the number of valid, non-zero-length, requests */
        assert(ncp->head != NULL);

    /* extract the requests from the linked list into a new linked list.
       In the meantime coalesce the linked list */

    for (i=0; i<num_reqs; i++) {
        if (req_ids[i] == NC_REQ_NULL) /* skip zero-size request */
            continue;

        assert(ncp->head != NULL);

        pre_req = NULL;
        cur_req = ncp->head;
        while (cur_req->id != req_ids[i]) { /* find from the linked list */
            pre_req = cur_req;
            cur_req = cur_req->next;
            if (cur_req == NULL) {
                printf("Error: no such request ID = %d\n", req_ids[i]);
                if (statuses != NULL)
                    statuses[i] = NC_EINVAL_REQUEST;
                /* retain the first error status */
                if (status == NC_NOERR)
                    status = NC_EINVAL_REQUEST;
            }
        }
        /* skip this invalid nonblocking request i */
        if (cur_req == NULL) continue;

        /* found it: cur_req */
        cur_req->status = statuses + i;
        for (j=0; j<cur_req->num_subreqs; j++)
            cur_req->subreqs[j].status = cur_req->status;

        /* remove cur_req from the ncp->head linked list */
        if (cur_req == ncp->head)  /* move cur_req to next */
            ncp->head = cur_req->next;
        else /* move pre_req and cur_req to next */
            pre_req->next = cur_req->next;

        if (cur_req->rw_flag == READ_REQ) { /* add cur_req to r_req_tail */
            if (r_req_head == NULL) {
                r_req_head = cur_req;
                r_req_tail = cur_req;
            }
            else {
                r_req_tail->next = cur_req;
                r_req_tail = r_req_tail->next;
            }
            r_req_tail->next = NULL;
            num_r_reqs += (cur_req->num_subreqs == 0) ? 1 : cur_req->num_subreqs;
            /* if this request is for record variable, then count only its
               subrequests (one for each individual record) */
        }
        else { /* add cur_req to w_req_tail */
            if (w_req_head == NULL) {
                w_req_head = cur_req;
                w_req_tail = cur_req;
            }
            else {
                w_req_tail->next = cur_req;
                w_req_tail = w_req_tail->next;
            }
            w_req_tail->next = NULL;
            num_w_reqs += (cur_req->num_subreqs == 0) ? 1 : cur_req->num_subreqs;
            /* if this request is for record variable, then count only its
               subrequests (one for each individual record) */
        }
    }
    /* make sure ncp->tail pointing to the tail */
    ncp->tail = ncp->head;
    while (ncp->tail != NULL && ncp->tail->next != NULL)
        ncp->tail = ncp->tail->next;

    if (io_method == COLL_IO) {
        int io_req[2], do_io[2];  /* [0]: do read, [1]: do write */
        io_req[0] = num_r_reqs;
        io_req[1] = num_w_reqs;
        MPI_Allreduce(&io_req, &do_io, 2, MPI_INT, MPI_MAX, ncp->nciop->comm);

        /* make sure if at least one process has a non-zero request, all
           processes participate the collective read/write */
        do_read  = do_io[0];
        do_write = do_io[1];
    }
    else {
        do_read  = num_r_reqs;
        do_write = num_w_reqs;
    }

    /* carry out reads and writes separately */
    if (do_read > 0)
        err = ncmpii_wait_getput(ncp, num_r_reqs, r_req_head,
                                 READ_REQ, io_method);

    if (do_write > 0)
        err = ncmpii_wait_getput(ncp, num_w_reqs, w_req_head,
                                 WRITE_REQ, io_method);

    /* retain the first error status */
    if (status == NC_NOERR) status = err;

    /* post-IO data processing: may need byte-swap user write buf, or
                                byte-swap and type convert user read buf */

    /* connect read and write request lists into a single list */
    if (r_req_head != NULL) {
        cur_req = r_req_head;
        r_req_tail->next = w_req_head;
    }
    else
        cur_req = w_req_head;

    pre_req = NULL;
    while (cur_req != NULL) {
        MPI_Offset fnelems, bnelems;
        MPI_Datatype ptype;
        NC_var *varp = cur_req->varp;

        fnelems = cur_req->fnelems;
        ptype   = cur_req->ptype;

        if (cur_req->rw_flag == WRITE_REQ) {
            /* must byte-swap the user buffer back to its original Endianess
               only when the buffer itself has been byte-swapped before,
               i.e. NOT iscontig_of_ptypes && NOT ncmpii_need_convert() &&
               ncmpii_need_swap()
            */
            if (cur_req->need_swap_back_buf)
                ncmpii_in_swapn(cur_req->buf, fnelems,
                                ncmpix_len_nctype(varp->type));
        } else { /* for read */
            /* now, xbuf contains the data read from the file.
             * It needs to be type-converted + byte-swapped to cbuf
             */
            void *cbuf, *lbuf;

            if (ncmpii_need_convert(varp->type, ptype) ) {
                if (cur_req->is_imap || !cur_req->iscontig_of_ptypes)
                    cbuf = NCI_Malloc(fnelems * varp->xsz);
                else
                    cbuf = cur_req->buf;

                /* type convert + byte swap from xbuf to cbuf */
                DATATYPE_GET_CONVERT(varp->type, cur_req->xbuf,
                                     cbuf, cur_req->bnelems, ptype)
                /* err is set in DATATYPE_GET_CONVERT() */
                /* keep the first error */
                if (*cur_req->status == NC_NOERR) *cur_req->status = err;
                if (status == NC_NOERR) status = err;
            } else {
                if (ncmpii_need_swap(varp->type, ptype))
                    ncmpii_in_swapn(cur_req->xbuf, fnelems,
                                    ncmpix_len_nctype(varp->type));
                cbuf = cur_req->xbuf;
            }

            if (cur_req->is_imap) { /* this request was made by get_varm() */
                if (cur_req->iscontig_of_ptypes)
                    lbuf = cur_req->buf;
                else {
                    int el_size;
                    MPI_Type_size(ptype, &el_size);
                    lbuf = NCI_Malloc(cur_req->lnelems*el_size);
                }

                /* unpack cbuf to lbuf based on imaptype */
                err = ncmpii_data_repack(cbuf, cur_req->bnelems, ptype,
                                         lbuf, 1, cur_req->imaptype);
                if (*cur_req->status == NC_NOERR)
                    /* keep the first error */
                    *cur_req->status = err;
                MPI_Type_free(&cur_req->imaptype);

                if (err != NC_NOERR) {
                    FREE_REQUEST(cur_req)
                    NCI_Free(cur_req);
                    return ((status != NC_NOERR) ? status : err);
                }
                /* cbuf is no longer needed
                 * if (cur_req->is_imap) cbuf cannot be == cur_req->buf */
                if (cbuf != cur_req->xbuf) NCI_Free(cbuf);
                cbuf = NULL;

                bnelems = cur_req->lnelems;
            } else { /* get_vars */
                lbuf = cbuf;
                bnelems = cur_req->bnelems;
            }

            if (!cur_req->iscontig_of_ptypes) {
                /* unpack lbuf to buf based on buftype */
                err = ncmpii_data_repack(lbuf, bnelems,
                                         ptype, cur_req->buf,
                                         cur_req->bufcount,
                                         cur_req->buftype);
                if (*cur_req->status == NC_NOERR) /* keep the first error */
                    *cur_req->status = err;
                if (err != NC_NOERR) {
                    FREE_REQUEST(cur_req)
                    NCI_Free(cur_req);
                    return ((status != NC_NOERR) ? status : err);
                }
            }
            /* lbuf is no longer needed */
            if (lbuf != cur_req->buf && lbuf != cur_req->xbuf)
                NCI_Free(lbuf);
        }
        /* free space allocated for the request objects */
        FREE_REQUEST(cur_req)
        pre_req = cur_req;
        cur_req = cur_req->next;
        NCI_Free(pre_req);
    }

    return status;
}

/* C struct for breaking down a request to a list of offset-length segments */
typedef struct {
    MPI_Offset off;      /* starting file offset of the request */
    MPI_Offset len;      /* requested length in bytes strating from off */
    MPI_Aint   buf_addr; /* distance of this request's I/O buffer to the first
                            request to be merged */
} off_len;

/*----< off_compare() >-------------------------------------------------------*/
/* used for sorting the offsets of the off_len array */
static int
off_compare(const void *a, const void *b)
{
    if (((off_len*)a)->off > ((off_len*)b)->off) return  1;
    if (((off_len*)a)->off < ((off_len*)b)->off) return -1;
    return 0;
}

/*----< ncmpii_flatten() >----------------------------------------------------*/
/* flatten a subarray request into a list of offset-length pairs */
static int
ncmpii_flatten(int          ndim,    /* number of dimensions */
               int          el_size, /* array element size */
               MPI_Offset  *dimlen,  /* dimension lengths */
               MPI_Offset   offset,  /* starting file offset of variable */
               MPI_Aint     buf_addr,/* starting buffer address */
               MPI_Offset  *start,   /* starts of subarray */
               MPI_Offset  *count,   /* counts of subarray */
               MPI_Offset  *stride,  /* strides of subarray */
               MPI_Offset  *nseg,    /* OUT: number of segments */
               off_len     *seg)     /* OUT: array of segments */
{
    int i, j, to_free_stride=0;
    MPI_Offset seg_len, nstride, array_len, off, subarray_len;
    off_len *ptr=seg, *seg0;

    *nseg = 0;
    if (ndim < 0) return *nseg;

    if (ndim == 0) {  /* 1D record variable */
        *nseg = 1;
        seg->off      = offset;
        seg->len      = el_size;
        seg->buf_addr = buf_addr;
        return *nseg;
    }

    if (stride == NULL) { /* equivalent to {1, 1, ..., 1} */
        stride = (MPI_Offset*) NCI_Malloc(ndim * sizeof(MPI_Offset));
        for (i=0; i<ndim; i++) stride[i] = 1;
        to_free_stride = 1;
    }

    /* TODO: check if all stride[] >= 1
       Q: Is it legal if any stride[] <= 0 ? */

    /* caclulate the number of offset-length pairs */
    *nseg = (stride[ndim-1] == 1) ? 1 : count[ndim-1];
    for (i=0; i<ndim-1; i++)
        *nseg *= count[i];
    if (*nseg == 0) {  /* not reachable, an error if count[] == 0 */
        if (to_free_stride) NCI_Free(stride);
        return *nseg;
    }

    /* the length of all segments are of the same size */
    seg_len  = (stride[ndim-1] == 1) ? count[ndim-1] : 1;
    seg_len *= el_size;
    nstride  = (stride[ndim-1] == 1) ? 1 : count[ndim-1];

    /* set the offset-length pairs for the lowest dimension */
    off = offset + start[ndim-1] * el_size;
    for (i=0; i<nstride; i++) {
        ptr->off       = off;
        ptr->len       = seg_len;
        ptr->buf_addr  = buf_addr;
        buf_addr      += seg_len;
        off           += stride[ndim-1] * el_size;
        ptr++;
    }
    ndim--;

    subarray_len = nstride;
    array_len = 1;
    /* for higher dimensions */
    while (ndim > 0) {
        /* array_len is global array size from lowest up to ndim */
        array_len *= dimlen[ndim];

        /* off is the global array offset for this dimension, ndim-1 */
        off = start[ndim-1] * array_len * el_size;

        /* update all offsets from lowest up to dimension ndim-1 */
        seg0 = seg;
        for (j=0; j<subarray_len; j++) {
            seg0->off += off;
            seg0++;
        }

        /* update each plan subarray of dimension ndim-1 */
        off = array_len * stride[ndim-1] * el_size;
        for (i=1; i<count[ndim-1]; i++) {
            seg0 = seg;
            for (j=0; j<subarray_len; j++) {
                ptr->off       = seg0->off + off;
                ptr->len       = seg_len;
                ptr->buf_addr  = buf_addr;
                buf_addr      += seg_len;
                ptr++;
                seg0++;
            }
            off += array_len * stride[ndim-1] * el_size;
        }
        ndim--;  /* move to next higher dimension */
        subarray_len *= count[ndim];
    }
    if (to_free_stride) NCI_Free(stride);

    return *nseg;
}

/*----< ncmpii_merge_requests() >---------------------------------------------*/
static int
ncmpii_merge_requests(NC          *ncp,
                      int          num_reqs,
                      NC_req      *reqs,    /* [num_reqs] */
                      int          rw_flag, /* WRITE_REQ or READ_REQ */
                      void       **buf,     /* OUT: 1st I/O buf addr */
                      MPI_Offset  *nsegs,   /* OUT: no. off-len pairs */
                      off_len    **segs)    /* OUT: [*nsegs] */
{
    int i, j, err, status=NC_NOERR, num_valid_reqs, ndims, is_recvar;
    MPI_Offset  nseg, *start, *count, *shape, *stride;
    MPI_Aint addr, buf_addr;

    /* first check and remove invalid reqs[], for example out-of-boundary
       requests, and set their status to NC error codes. Because some
       requests are subrequests of the same record-variable request, their
       status pointed to the same memory address of the parent request's
       status. Hece, if one of the subrequests is invalid, we mark the
       parent request invalid and skip it.
     */
    for (i=0; i<num_reqs; i++) {
        /* note that all reqs[].status have been initialized to NC_NOERR
           in ncmpii_wait()
         */
        if (*reqs[i].status != NC_NOERR) /* skip this invalid reqs[i] */
            continue;

        is_recvar = IS_RECVAR(reqs[i].varp);

        /* Check for out-of-boundary requests */
        err = NCedgeck(ncp, reqs[i].varp, reqs[i].start, reqs[i].count);
        if (err != NC_NOERR ||
            (rw_flag == READ_REQ && is_recvar &&
             reqs[i].start[0] + reqs[i].count[0] > NC_get_numrecs(ncp))) {
            err = NCcoordck(ncp, reqs[i].varp, reqs[i].start);
            if (err != NC_NOERR) { /* status is the 1st encounter errror */
                *reqs[i].status = err;
                status = (status == NC_NOERR) ? err : status;
            }
            else {
                *reqs[i].status = NC_EEDGE;
                status = (status == NC_NOERR) ? NC_EEDGE : status;
            }
        }
    }

    /* find the number of valid requests in reqs[] */
    num_valid_reqs = 0;
    for (i=0; i<num_reqs; i++) {
        if (*reqs[i].status != NC_NOERR) continue;
        num_valid_reqs++; /* incease the number of valid requests */
    }

    *nsegs = 0;    /* total number of offset-length pairs */
    *buf   = NULL; /* I/O buffer of first valid request */
    *segs  = NULL; /* I/O buffer of first valid request */
    if (num_valid_reqs == 0)
        return status;

    /* Count the number off-len pairs from the remaining valid reqs[], so we
       can malloc a contiguous memory space for storing off-len pairs
     */
    for (i=0; i<num_reqs; i++) {
        if (*reqs[i].status != NC_NOERR) continue; /* skip invalid one */

        /* buf_addr is the buffer address of the first valid request */
        if (*buf == NULL) {
#ifdef HAVE_MPI_GET_ADDRESS
            MPI_Get_address(reqs[i].xbuf, &buf_addr);
#else
            MPI_Address(reqs[i].xbuf, &buf_addr);
#endif
            *buf = reqs[i].xbuf;
        }

        is_recvar = IS_RECVAR(reqs[i].varp);

        /* for record variable, each reqs[] is within a record */
        ndims  = (is_recvar) ? reqs[i].ndims  - 1 : reqs[i].ndims;
        count  = (is_recvar) ? reqs[i].count  + 1 : reqs[i].count;
        stride = NULL;
        if (reqs[i].stride != NULL)
            stride = (is_recvar) ? reqs[i].stride + 1 : reqs[i].stride;

        if (ndims < 0) continue;
        if (ndims == 0) {  /* 1D record variable */
            (*nsegs)++;
            continue;
        }
        nseg = 1;
        if (stride != NULL && stride[ndims-1] > 1)
            nseg = count[ndims-1];  /* count of last dimension */
        for (j=0; j<ndims-1; j++)
            nseg *= count[j];  /* all count[] except the last dimension */

        *nsegs += nseg;
    }

    /* now we can allocate a contiguous memory space for the off-len pairs */
    off_len *seg_ptr = (off_len*) NCI_Malloc(*nsegs * sizeof(off_len));
    *segs = seg_ptr;

    /* now re-run the loop to fill in the off-len pairs */
    for (i=0; i<num_reqs; i++) {
        if (*reqs[i].status != NC_NOERR) continue; /* skip invalid one */

        /* buf_addr is the buffer address of the first valid request */
#ifdef HAVE_MPI_GET_ADDRESS
        MPI_Get_address(reqs[i].xbuf, &addr);
#else
        MPI_Address(reqs[i].xbuf, &addr);
#endif
        addr -= buf_addr,  /* distance to the buf of first req */

        is_recvar = IS_RECVAR(reqs[i].varp);

        /* for record variable, each reqs[] is within a record */
        ndims  = (is_recvar) ? reqs[i].ndims  - 1 : reqs[i].ndims;
        start  = (is_recvar) ? reqs[i].start  + 1 : reqs[i].start;
        count  = (is_recvar) ? reqs[i].count  + 1 : reqs[i].count;
        shape  = (is_recvar) ? reqs[i].varp->shape  + 1 :
                               reqs[i].varp->shape;
        stride = NULL;
        if (reqs[i].stride != NULL)
            stride = (is_recvar) ? reqs[i].stride + 1 : reqs[i].stride;

        /* find the starting file offset for this record */
        MPI_Offset var_begin = reqs[i].varp->begin;
        if (is_recvar) var_begin += reqs[i].start[0] * ncp->recsize;

        /* flatten each request to a list of offset-length pairs */
        ncmpii_flatten(ndims, reqs[i].varp->xsz, shape, var_begin,
                       addr, start, count, stride,
                       &nseg,    /* OUT: number of offset-length pairs */
                       seg_ptr); /* OUT: array of offset-length pairs */
        seg_ptr += nseg; /* append the list to the end of segs array */
    }

    /* sort the off-len array, segs[], in an increasing order */
    qsort(*segs, *nsegs, sizeof(off_len), off_compare);

#ifdef USE_MULTIPLE_IO_METHOD
    /* check for any overlapping offset-length pairs. If found, we will fall
       back to the group method, i.e. perform multiple MPI-IO calls */
    for (i=0; i<*nsegs-1; i++)
        if ((*segs)[i].off + (*segs)[i].len > (*segs)[i+1].off)
            break;

    if (i < *nsegs-1) {  /* overlap found */
        NCI_Free(*segs); /* free up allocated space */
        *segs    = NULL;
        *nsegs   = 0;    /* indicating merging failed */
    }
#else
    /* merge the overlapped requests, skip the overlapped regions for those
       requests with higher j indices (i.e. requests with lower j indices
       win the writes to the overlapped regions)
     */
    for (i=0, j=1; j<*nsegs; j++) {
        if ((*segs)[i].off + (*segs)[i].len >= (*segs)[j].off + (*segs)[j].len)
            /* segment i completely covers segment j, skip j */
            continue;

        MPI_Offset gap = (*segs)[i].off + (*segs)[i].len - (*segs)[j].off;
        if (gap >= 0) { /* segments i and j overlaps */
            if ((*segs)[i].buf_addr + (*segs)[i].len ==
                (*segs)[j].buf_addr + gap) {
                /* buffers i and j are contiguous, merge j to i */
                (*segs)[i].len += (*segs)[j].len - gap;
            }
            else { /* buffers are not contiguous, reduce j's len */
                (*segs)[i+1].off      = (*segs)[j].off + gap;
                (*segs)[i+1].len      = (*segs)[j].len - gap;
                (*segs)[i+1].buf_addr = (*segs)[j].buf_addr + gap;
                i++;
            }
        }
        else { /* i and j do not overlap */
            i++;
            if (i < j) (*segs)[i] = (*segs)[j];
        }
    }

    /* update number of segments, now all off-len pairs are not overlapped */
    *nsegs = i+1;
#endif
    return status;
}

/*----< ncmpii_getput_merged_requests() >-------------------------------------*/
static int
ncmpii_getput_merged_requests(NC         *ncp,
                              int         rw_flag,  /* WRITE_REQ or READ_REQ */
                              int         io_method,/* COLL_IO or INDEP_IO */
                              MPI_Offset  nsegs,    /* no. off-len pairs */
                              off_len    *segs,     /* [nsegs] off-en pairs */
                              void       *buf)
{
    int i, j, mpireturn, status=NC_NOERR, mpi_err=NC_NOERR;
    int         *blocklengths;
    MPI_Aint    *displacements;
    MPI_Offset   next_off, next_len;
    MPI_File     fh;
    MPI_Status   mpistatus;
    MPI_Datatype filetype, buftype;

    assert(nsegs > 0);

    if (io_method == COLL_IO)
        fh = ncp->nciop->collective_fh;
    else
        fh = ncp->nciop->independent_fh;

    /* create the file view MPI derived data type by concatenating the sorted
       offset-length pairs */

    /* For filetype, the segs[].off can be further coalesced. For example,
       when writing a consecutive columns of a 2D array, even though the I/O
       buffer addresses may not be able to coalesced, the file offsets on
       the same row can be coalesced. Thus, first calculate the length of
       coalesced off-len pairs (the memory space needed for malloc) */
    next_off = segs[0].off;
    next_len = segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (next_off + next_len == segs[i].off) /* j and i are contiguous */
            next_len += segs[i].len;
        else {
            j++;
            next_off = segs[i].off;
            next_len = segs[i].len;
        }
    }
    /* j+1 is the coalesced length */
    blocklengths  = (int*)      NCI_Malloc((j+1) * sizeof(int));;
    displacements = (MPI_Aint*) NCI_Malloc((j+1) * sizeof(MPI_Aint));

    /* coalesce segs[].off and len to dispalcements[] and blocklengths[] */
    displacements[0] = segs[0].off;
    blocklengths[0]  = segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (displacements[j] + blocklengths[j] == segs[i].off)
            /* j and i are contiguous */
            blocklengths[j] += segs[i].len;
        else {
            j++;
            displacements[j] = segs[i].off;
            blocklengths[j]  = segs[i].len;
        }
    }
    /* j+1 is the coalesced length */
// printf("filetype new j+1=%d blocklengths[0]=%d\n",j+1,blocklengths[0]);

    MPI_Type_create_hindexed(j+1, blocklengths, displacements, MPI_BYTE,
                             &filetype);
    MPI_Type_commit(&filetype);

    /* now we are ready to call MPI_File_set_view */
    mpireturn = MPI_File_set_view(fh, 0, MPI_BYTE, filetype, "native",
                                  MPI_INFO_NULL);
    CHECK_MPI_ERROR(mpireturn, "MPI_File_set_view", NC_EFILE);
    MPI_Type_free(&filetype);
    NCI_Free(displacements);
    NCI_Free(blocklengths);

    /* create the I/O buffer derived data type from the I/O buffer's
       offset-length pairs */

    /* Although it is unlikely buffers can be coalesced, it will not harm to
       check it out */
    next_off = segs[0].buf_addr;
    next_len = segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (next_off + next_len == segs[i].buf_addr)
            /* j and i are contiguous */
            next_len += segs[i].len;
        else {
            j++;
            next_off = segs[i].buf_addr;
            next_len = segs[i].len;
        }
    }
    /* j+1 is the coalesced length */
    blocklengths  = (int*)      NCI_Malloc((j+1) * sizeof(int));;
    displacements = (MPI_Aint*) NCI_Malloc((j+1) * sizeof(MPI_Aint));

    /* coalesce segs[].off and len to dispalcements[] and blocklengths[] */
    displacements[0] = segs[0].buf_addr;
    blocklengths[0]  = segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (displacements[j] + blocklengths[j] == segs[i].buf_addr)
            /* j and i are contiguous */
            blocklengths[j] += segs[i].len;
        else {
            j++;
            displacements[j] = segs[i].buf_addr;
            blocklengths[j]  = segs[i].len;
        }
    }
    /* j+1 is the coalesced length */
// printf("buftype new j+1=%d blocklengths[0]=%d\n",j+1,blocklengths[0]);
    MPI_Type_create_hindexed(j+1, blocklengths, displacements, MPI_BYTE,
                             &buftype);
    MPI_Type_commit(&buftype);
    NCI_Free(displacements);
    NCI_Free(blocklengths);

    /* now we are ready to call MPI read/write APIs */
    if (rw_flag == READ_REQ) {
        if (io_method == COLL_IO) {
            mpireturn = MPI_File_read_all(fh, buf, 1, buftype, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_read_all", NC_EREAD)
        } else {
            mpireturn = MPI_File_read(fh, buf, 1, buftype, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_read", NC_EREAD)
        }
        int get_size;
        MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
        ncp->nciop->get_size += get_size;
    } else { /* WRITE_REQ */
        if (io_method == COLL_IO) {
            mpireturn = MPI_File_write_all(fh, buf, 1, buftype, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_write_all", NC_EWRITE)
        } else {
            mpireturn = MPI_File_write(fh, buf, 1, buftype, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_write", NC_EWRITE)
        }
        int put_size;
        MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
        ncp->nciop->put_size += put_size;
    }
    MPI_Type_free(&buftype);

    /* reset fileview so the entire file is visible again */
    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    /* make NC error higher priority than MPI error */
    return ((status != NC_NOERR) ? status : mpi_err);
}

/*----< req_compare() >-------------------------------------------------------*/
/* used to sort the the string file offsets of reqs[] */
static int
req_compare(const NC_req *a, const NC_req *b)
{
    if (a->offset_start > b->offset_start) return (1);
    if (a->offset_start < b->offset_start) return (-1);
    return (0);
}

/*----< ncmpii_wait_getput() >------------------------------------------------*/
static int
ncmpii_wait_getput(NC     *ncp,
                   int     num_reqs,  /* # requests including subrequests */
                   NC_req *req_head,  /* linked list not include subrequests */
                   int     rw_flag,   /* WRITE_REQ or READ_REQ */
                   int     io_method) /* COLL_IO or INDEP_IO */
{
    int i, j, k, err, status=NC_NOERR, ngroups, max_ngroups, *group_index;
    void *merged_buf;
    NC_var **varps;
    NC_req *reqs, *cur_req;
    int (*fcnt)(const void*, const void*);

    if (num_reqs == 0) {
        ngroups = 0;
    }
    else {
        /* change the linked list into an array to be sorted */
        reqs = (NC_req*) NCI_Malloc(num_reqs * sizeof(NC_req));
        i = 0;
        cur_req = req_head;
        while (cur_req != NULL) {
            reqs[i++] = *cur_req;
            if (cur_req->num_subreqs > 0) i--;
            for (j=0; j<cur_req->num_subreqs; j++)
                reqs[i++] = cur_req->subreqs[j];
            cur_req = cur_req->next;
        }
        assert(i == num_reqs);

        /* divide the requests into groups of non-interleaved requests
         * This is because MPI collective I/O requires each process's fileview
         * contains monotonic non-decreasing file offsets
         */

        /* first sort reqs[] based on reqs[].offset_start */
        fcnt = (int (*)(const void *, const void *))req_compare;
        qsort(reqs, num_reqs, sizeof(NC_req), fcnt);

        /* find the non-interleaving groups of requests */
        ngroups = 1;
        group_index = (int*) NCI_Malloc(num_reqs * sizeof(int));

        j = 1;
        group_index[0] = 0;  /* first reqs[] index of the group */
        while (j < num_reqs) {
            for (i=j; i<num_reqs; i++) {
                if (reqs[j-1].offset_end < reqs[i].offset_start) {
                    /* shift reqs[j ... i-1] to reqs[j+1 ... i] */
                    if (i > j) {
                        NC_req tmp = reqs[i];
                        for (k=i; k>j; k--)
                            reqs[k] = reqs[k-1];
                        /* move reqs[i] to reqs[j] */
                        reqs[j] = tmp;
                    }
                    j++;
                }
                /* else do nothing */
            }
            if (j < num_reqs) /* next group starts with reqs[j] */
                group_index[ngroups++] = j++;
        }
    }

    /* Each group contains requests whose starting file offsets have been
       sorted in a monotonically non-descreasing order. An MPI-IO call will
       be made for each group. The file view of requests in a group are
       appended one after another into a single one. Similarly, an I/O
       buffer derived data type is created by concatenating their I/O buffers.

       However, requests across groups can be merged in order to reduce the
       number of groups, and hence the number of  MPI-IO calls to further
       improve the performance. For example, multiple nonblocking requests
       each accessing a single column of a 2D array.

       The pitfall of the merging is that it will have to break down each
       request into a list of offset-length pairs, and merge all lists into
       a sorted list based on their offsets in an increasing order.

       Be warned! The additional memory requirement for this merging can be
       more than the I/O data itself. For example, each nonblocking request
       access a single column of a 2D array of 4-byte integer type. Each
       off-len pair represents only a 4-byte integer, but the off-len pair
       itself takes 24 bytes.
     */
    MPI_Offset  nsegs=0;   /* number of merged offset-length pairs */
    off_len    *segs=NULL; /* array of the offset-length pairs */
    if (ngroups > 1) {
        /* try merge all requests into sorted offset-length pairs */
        err = ncmpii_merge_requests(ncp, num_reqs, reqs, rw_flag,
                                    &merged_buf, &nsegs, &segs);
        if (status == NC_NOERR) status = err;
        /* this returned status is the first error encountered */

        if (nsegs > 0) ngroups = 1;
        /* when nsegs > 0, the merge is successful, the returned "segs"
           points to an array of off-len pairs containing only the valid
           requests. Invalid requests in reqs[], such as out-of-boundary
           ones, are detected and removed, setting their status, reqs[].status,
           to proper NC error codes.

           If macro USE_MULTIPLE_IO_METHOD is defined,
               when nsegs == 0, it indicates merge failed becaue the fillowing
               2 causes. One, all requests are invalid. The other is because
               an overlapped between two offset-length pairs is detected. In
               both case, we fall back to the multi-group method which makes
               an MPI-IO call for each group.
           else if macro USE_MULTIPLE_IO_METHOD is NOT defined,
               merge is successful and overlaps have been resolved. I/O will
               be performed using one MPI-IO call.
         */
    }

#ifdef USE_MULTIPLE_IO_METHOD
    if (io_method == COLL_IO)
        /* I/O will be completed in max_ngroups number of collective I/O */
        MPI_Allreduce(&ngroups, &max_ngroups, 1, MPI_INT, MPI_MAX,
                      ncp->nciop->comm);
    else /* INDEP_IO */
        max_ngroups = ngroups;
#else
    max_ngroups = 1;
#endif
// printf("ngroups=%d max_ngroups=%d nsegs=%lld\n",ngroups,max_ngroups,nsegs);

    if (nsegs > 0) { /* requests are merged */
        /* sges[] will be used to construct fileview and buffer type */
        err = ncmpii_getput_merged_requests(ncp, rw_flag, io_method, nsegs,
                                            segs, merged_buf);
        /* preserve the previous error if there is any */
        if (status == NC_NOERR) status = err;
        NCI_Free(segs);
    }
    else if (num_reqs > 0) {  /* also means ngroups > 0 */
        int *statuses;
        void **bufs;
        MPI_Offset **starts, **counts, **strides, *nbytes;

        starts   = (MPI_Offset**) NCI_Malloc(3 * num_reqs *sizeof(MPI_Offset*));
        counts   = starts + num_reqs;
        strides  = counts + num_reqs;
        nbytes   = (MPI_Offset*) NCI_Malloc(num_reqs * sizeof(MPI_Offset));
        varps    = (NC_var**)    NCI_Malloc(num_reqs * sizeof(NC_var*));
        bufs     = (void**)      NCI_Malloc(num_reqs * sizeof(void*));
        statuses = (int*)        NCI_Malloc(num_reqs * sizeof(int));

        for (i=0; i<ngroups; i++) {
            int i_num_reqs = (i == ngroups-1) ? num_reqs : group_index[i+1];
            i_num_reqs -= group_index[i]; /* number of requests in group i */

            k = group_index[i];
            for (j=0; j<i_num_reqs; j++, k++) {
                varps[j]    = reqs[k].varp;
                starts[j]   = reqs[k].start;
                counts[j]   = reqs[k].count;
                strides[j]  = reqs[k].stride;
                nbytes[j]   = reqs[k].fnelems * varps[j]->xsz;
                bufs[j]     = reqs[k].xbuf;
                statuses[j] = NC_NOERR;
            }

            err = ncmpii_mgetput(ncp, i_num_reqs, varps, starts, counts,
                                 strides, bufs, nbytes, statuses, rw_flag,
                                 io_method);
            if (status == NC_NOERR) status = err;
            /* retain the first error if there is any */

            /* update status */
            for (j=0; j<i_num_reqs; j++)
                if (*reqs[group_index[i] + j].status == NC_NOERR)
                    *reqs[group_index[i] + j].status = statuses[j];
        }
        NCI_Free(statuses);
        NCI_Free(bufs);
        NCI_Free(varps);
        NCI_Free(nbytes);
        NCI_Free(starts);
    }

    /* All processes must participate all max_ngroups collective calls,
       for example, a subset of processes have zero-size data to write */
    for (i=ngroups; i<max_ngroups; i++)
        ncmpii_mgetput(ncp, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       rw_flag, io_method);

    /* update the number of records if new records are added */
    if (rw_flag == WRITE_REQ) {
        /* Because netCDF allows only one unlimited dimension, find the
         * maximum number of records from all nonblocking requests and
         * update newnumrecs once
         */
        MPI_Offset max_newnumrecs = ncp->numrecs;
        for (i=0; i<num_reqs; i++) {
            if (*reqs[i].status == NC_NOERR && IS_RECVAR(reqs[i].varp)) {
                /* update the number of records in NC */
                MPI_Offset newnumrecs = reqs[i].start[0] + reqs[i].count[0];
                max_newnumrecs = MAX(max_newnumrecs, newnumrecs);
            }
        }
        err = ncmpii_update_numrecs(ncp, max_newnumrecs);
        /* ncmpii_update_numrecs() is collective, so even if this process
         * finds its max_newnumrecs being zero, it still needs to participate
         * the call
         */
        if (status == NC_NOERR) status = err;
        /* retain the first error if there is any */
    }

    if (num_reqs > 0) {
        NCI_Free(group_index);
        NCI_Free(reqs);
    }

    return status;
}

/*----< ncmpii_mgetput() >----------------------------------------------------*/
static int
ncmpii_mgetput(NC           *ncp,
               int           num_reqs,
               NC_var       *varps[],     /* [num_reqs] */
               MPI_Offset   *starts[],    /* [num_reqs] */
               MPI_Offset   *counts[],    /* [num_reqs] */
               MPI_Offset   *strides[],   /* [num_reqs] */
               void         *bufs[],      /* [num_reqs] */
               MPI_Offset    nbytes[],    /* [num_reqs] */
               int           statuses[],  /* [num_reqs] */
               int           rw_flag,     /* WRITE_REQ or READ_REQ */
               int           io_method)   /* COLL_IO or INDEP_IO */
{
    int i, len, status=NC_NOERR, mpireturn, mpi_err=NC_NOERR;
    void *buf;
    MPI_Status mpistatus;
    MPI_Datatype buf_type;
    MPI_File fh;

    buf_type = MPI_BYTE;
    if (io_method == COLL_IO)
        fh = ncp->nciop->collective_fh;
    else
        fh = ncp->nciop->independent_fh;

    /* set the MPI file view */
    status = ncmpii_mset_fileview(fh, ncp, num_reqs, varps, starts, counts,
                                  strides, statuses, rw_flag);
    if (status != NC_NOERR) { /* skip this request */
        if (io_method == COLL_IO)
            num_reqs = 0; /* still need to participate the successive
                             collective calls, read/write/setview */
        else
            return status;
    }

    /* create the I/O buffer derived data type */
    if (num_reqs > 0) {
        int *blocklengths = (int*) NCI_Malloc(num_reqs * sizeof(int));
        MPI_Aint *disps = (MPI_Aint*) NCI_Malloc(num_reqs*sizeof(MPI_Aint));
        MPI_Aint a0, ai;

        disps[0] = 0;
#ifdef HAVE_MPI_GET_ADDRESS
        MPI_Get_address(bufs[0], &a0);
#else 
        MPI_Address(bufs[0], &a0);
#endif
        blocklengths[0] = nbytes[0];
        for (i=1; i<num_reqs; i++) {/*loop for multi-variables*/
/* wkliao: type convert from MPI_Offset nbytes[i] to int blocklengths[i]
 *         Can we do someting smarter here ? */
            blocklengths[i] = nbytes[i];
#ifdef HAVE_MPI_GET_ADDRESS
            MPI_Get_address(bufs[i], &ai);
#else
	    MPI_Address(bufs[i], &ai);
#endif
            disps[i] = ai - a0;
        }
        /* concatenate buffer addresses into a single buffer type */
        MPI_Type_hindexed(num_reqs, blocklengths, disps, MPI_BYTE, &buf_type);
        MPI_Type_commit(&buf_type);

        NCI_Free(disps);
        NCI_Free(blocklengths);

        buf = bufs[0];
        len = 1;
    }
    else { /* num_reqs == 0, simply participate the collective call */
        buf      = NULL;
        len      = 0;
        buf_type = MPI_BYTE;
    }

    if (rw_flag == READ_REQ) {
        if (io_method == COLL_IO) {
            mpireturn = MPI_File_read_all(fh, buf, len, buf_type, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_read_all", NC_EREAD)
        } else {
            mpireturn = MPI_File_read(fh, buf, len, buf_type, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_read", NC_EREAD)
        }
        int get_size;
        MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
        ncp->nciop->get_size += get_size;

    } else { /* WRITE_REQ */
        if (io_method == COLL_IO) {
            mpireturn = MPI_File_write_all(fh, buf, len, buf_type, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_write_all", NC_EWRITE)
        } else {
            mpireturn = MPI_File_write(fh, buf, len, buf_type, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_write", NC_EWRITE)
        }
        int put_size;
        MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
        ncp->nciop->put_size += put_size;
    }

    if (buf_type != MPI_BYTE)
        MPI_Type_free(&buf_type);

    /* reset fileview so the entire file is visible again */
    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    return ((status != NC_NOERR) ? status : mpi_err);
}
