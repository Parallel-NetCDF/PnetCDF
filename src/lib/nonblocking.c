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

/* Prototypes for functions used only in this file */
static int ncmpii_wait_getput(NC *ncp, int num_reqs, NC_req *req_head,
                              int rw_flag, int io_method);

static int ncmpii_mgetput(NC *ncp, int num_reqs, NC_var *varps[],
                          MPI_Offset *starts[], MPI_Offset *counts[],
                          MPI_Offset *strides[], void *bufs[],
                          MPI_Offset nbytes[], int statuses[], int rw_flag,
                          int io_method);

static int ncmpii_mset_fileview(MPI_File fh, NC* ncp, int ntimes, NC_var **varp,
                                MPI_Offset *starts[], MPI_Offset *counts[],
                                MPI_Offset *strides[], int statuses[],
                                int rw_flag);

static int req_compare(const NC_req *a, const NC_req *b);

#define FREE_REQUEST(req) {                                          \
    if (req->xbuf == req->buf) req->xbuf = NULL;                     \
    if (req->cbuf == req->buf) req->cbuf = NULL;                     \
    if (req->xbuf != req->cbuf && req->xbuf != NULL) {               \
        NCI_Free(req->xbuf);                                         \
        req->xbuf = NULL;                                            \
    }                                                                \
    if (req->cbuf != req->buf && req->cbuf != NULL) {                \
        NCI_Free(req->cbuf);                                         \
        req->cbuf = NULL;                                            \
    }                                                                \
    if (req->is_imap && req->cbuf != NULL) { /* imap cbuf != buf */  \
        NCI_Free(req->cbuf);                                         \
        req->cbuf = NULL;                                            \
    }                                                                \
    for (j=0; j<req->num_subreqs; j++)                               \
        NCI_Free(req->subreqs[j].start);                             \
    if (req->num_subreqs > 0)                                        \
        NCI_Free(req->subreqs);                                      \
    NCI_Free(req->start);                                            \
    if (req->stride != NULL)                                         \
	NCI_Free(req->stride);                                       \
}

/*----< ncmpi_cancel() >-----------------------------------------------------*/
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

/*----< ncmpii_cancel() >----------------------------------------------------*/
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

/*----< ncmpi_wait() >-------------------------------------------------------*/
/* ncmpi_wait() is an independent call */
int
ncmpi_wait(int ncid,
           int num_reqs,
           int *req_ids,  /* [num_reqs] */
           int *statuses) /* [num_reqs] */
{
#ifndef ENABLE_NONBLOCKING
    int i, ret_st=NC_NOERR;
#endif
    int status;
    NC  *ncp;

    CHECK_NCID
#ifdef ENABLE_NONBLOCKING
    return ncmpii_wait(ncp, INDEP_IO, num_reqs, req_ids, statuses);
#else
    for (i=0; i<num_reqs; i++) { /* serve one request at a time */
        status = ncmpii_wait(ncp, INDEP_IO, 1, &req_ids[i], &statuses[i]);
        if (status != NC_NOERR)
            ret_st = status;
    }
    return ret_st; /* return the last fail status, if there is any */
#endif
}

/*----< ncmpi_wait_all() >---------------------------------------------------*/
/* ncmpi_wait_all() is a collective call */
int
ncmpi_wait_all(int  ncid,
               int  num_reqs, /* number of requests */
               int *req_ids,  /* [num_reqs] */
               int *statuses) /* [num_reqs] */
{
#ifndef ENABLE_NONBLOCKING
    int  i, ret_st=NC_NOERR;
#endif
    int status;
    NC  *ncp;

    CHECK_NCID
#ifdef ENABLE_NONBLOCKING
    return ncmpii_wait(ncp, COLL_IO, num_reqs, req_ids, statuses);
#else
    /* This API is collective, so make it illegal in indep mode. This also
       ensures the program will returns back to collective mode. */
    if (NC_indep(ncp))
        return NC_EINDEP;

    /* must enter independent mode, as num_reqs may be different among
       processes */
    ncmpi_begin_indep_data(ncid);

    for (i=0; i<num_reqs; i++) { /* serve one request at a time */
        status = ncmpii_wait(ncp, INDEP_IO, 1, &req_ids[i], &statuses[i]);
        if (status != NC_NOERR)
            ret_st = status;
    }

    /* return to collective data mode */
    ncmpi_end_indep_data(ncid);

    return ret_st; /* return the last fail status, if there is any */
#endif
}

/*----< ncmpii_mset_fileview() >---------------------------------------------*/
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
    int i, status, mpireturn, *blocklens;
    MPI_Datatype full_filetype, *filetypes;
    MPI_Offset *offsets;

    if (ntimes <= 0) { /* participate collective call */
        mpireturn = MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                      MPI_INFO_NULL);
        CHECK_MPI_ERROR("MPI_File_set_view", NC_EFILE)
        return NC_NOERR;
    }

    blocklens = (int*)          NCI_Malloc(ntimes * sizeof(int));
    offsets   = (MPI_Offset*)   NCI_Malloc(ntimes * sizeof(MPI_Offset));
    filetypes = (MPI_Datatype*) NCI_Malloc(ntimes * sizeof(MPI_Datatype));

    /* create a filetype for each variable */
    for (i=0; i<ntimes; i++) {
        status = ncmpii_vars_create_filetype(ncp,
                                             varp[i],
                                             starts[i],
                                             counts[i],
                                             strides[i],
                                             rw_flag,
                                             &offsets[i],
                                             &filetypes[i]);
        statuses[i] = status;
        if (status != NC_NOERR)
            return status;

        blocklens[i] = 1;
        if (filetypes[i] == MPI_BYTE) { /* file type is contiguous in file */
            int j;
            blocklens[i] = varp[i]->xsz;
            for (j=0; j<varp[i]->ndims; j++)
                blocklens[i] *= counts[i][j];
        }
    }

    if (ntimes == 1) {  /* no two or more filetypes to combine */
        full_filetype = filetypes[0];
        /* filetypes[0] has been committed already */

        mpireturn = MPI_File_set_view(fh, offsets[0], MPI_BYTE, full_filetype,
                                      "native", MPI_INFO_NULL);
        CHECK_MPI_ERROR("MPI_File_set_view", NC_EFILE)
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
        CHECK_MPI_ERROR("MPI_File_set_view", NC_EFILE)
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

    return NC_NOERR;
}

/*----< req_compare() >------------------------------------------------------*/
static int
req_compare(const NC_req *a, const NC_req *b)
{
    if (a->offset_start > b->offset_start) return (1);
    if (a->offset_start < b->offset_start) return (-1);
    return (0);
}

/*----< ncmpii_wait() >------------------------------------------------------*/
int
ncmpii_wait(NC  *ncp,
            int  io_method, /* COLL_IO or INDEP_IO */
            int  num_reqs,  /* number of requests */
            int *req_ids,   /* [num_reqs] */
            int *statuses)  /* [num_reqs] */
{
    int i, j, warning=NC_NOERR, status, num_w_reqs, num_r_reqs;
    int do_read, do_write;
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
    /* j now is the true number of non-zeor length request */

    if (io_method == INDEP_IO && j == 0)
        return NC_NOERR;

    /* For collective APIs, even thorugh some processes may have no
       requests at all, but they are still required to participate the call */

    w_req_head = w_req_tail = NULL;
    r_req_head = r_req_tail = NULL;
    num_w_reqs = num_r_reqs = 0;

    if (j > 0) /* j is the number of valid, non-zero-length, requests */
        assert(ncp->head != NULL);

    /* extract the requests from the linked list into a new linked list.
       In the meantime coalesce the linked list */

    for (i=0; i<num_reqs; i++) {
        if (req_ids[i] == NC_REQ_NULL) /* skip zeor-size request */
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
                return NC_EINVAL_REQUEST;
            }
        }
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
            /* if this request is for record variable, then count only
               the subrequests (for individual record) */
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
            /* if this request is for record variable, then count only
               the subrequests (for individual record) */
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

        do_read  = do_io[0];
        do_write = do_io[1];
    }
    else {
        do_read  = num_r_reqs;
        do_write = num_w_reqs;
    }

    /* call ncmpii_wait for read and write separately */
    if (do_read > 0)
        ncmpii_wait_getput(ncp, num_r_reqs, r_req_head, READ_REQ, io_method);

    if (do_write > 0)
        ncmpii_wait_getput(ncp, num_w_reqs, w_req_head, WRITE_REQ, io_method);

    /* post-IO data processing: may need byte-swap user write buf, or
                                byte-swap and type convert for read buf */

    /* connect read and write request lists into a single list */
    if (r_req_head != NULL) {
        cur_req = r_req_head;
        r_req_tail->next = w_req_head;
    }
    else
        cur_req = w_req_head;

    pre_req = NULL;
    while (cur_req != NULL) {
        MPI_Offset nelems, cnelems;
        MPI_Datatype ptype;
        NC_var *varp = cur_req->varp;

        nelems  = cur_req->nelems;
        cnelems = cur_req->cnelems;
        ptype   = cur_req->ptype;

        if (cur_req->rw_flag == WRITE_REQ) {
            /* must byte-swap the user buffer back to its original endianess
               only when the buffer itself has been byte-swapped before.
               I.e. ! iscontig_of_ptypes && not ncmpii_need_convert() &&
               ncmpii_need_swap()
            */
            if (cur_req->iscontig_of_ptypes &&
                !ncmpii_need_convert(varp->type, ptype) &&
                ncmpii_need_swap(varp->type, ptype))
                ncmpii_in_swapn(cur_req->buf, nelems, ncmpix_len_nctype(varp->type));
        } else { /* for read */
            if (cur_req->is_imap) {
                /* layout cbuf to lbuf based on imap */
                status = ncmpii_data_repack(cur_req->cbuf, cnelems,
                                            ptype, cur_req->lbuf, 1,
                                            cur_req->imaptype);
                if (*(cur_req->status) == NC_NOERR) /* keep the first error */
                    *(cur_req->status) = status;
                if (status != NC_NOERR) {
                    FREE_REQUEST(cur_req)
                    NCI_Free(cur_req);
                    return status;
                }

                MPI_Type_free(&cur_req->imaptype);

                if (!cur_req->iscontig_of_ptypes) {
                    /* repack it back, like a read-modify-write operation */
                    status = ncmpii_data_repack(cur_req->lbuf, cur_req->lnelems,
                                                ptype, cur_req->buf,
                                                cur_req->bufcount,
                                                cur_req->datatype);
                    if (cur_req->lbuf != NULL) {
                        if (cur_req->cbuf == cur_req->lbuf) cur_req->cbuf = NULL;
                        if (cur_req->xbuf == cur_req->lbuf) cur_req->xbuf = NULL;
                        NCI_Free(cur_req->lbuf);
                        cur_req->lbuf = NULL;
                        /* cur_req->lbuf will never == cur_req->buf */
                    }
                    if (*(cur_req->status) == NC_NOERR)
                        /* keep the first error */
                        *(cur_req->status) = status;
                    if (status != NC_NOERR) {
                        FREE_REQUEST(cur_req)
                        NCI_Free(cur_req);
                        return status;
                    }
                }
            }
            else {
                if ( ncmpii_need_convert(varp->type, ptype) ) {
                    /* automatic numeric datatype conversion */
                    DATATYPE_GET_CONVERT(varp->type, cur_req->xbuf,
                                         cur_req->cbuf, cnelems, ptype)
                    if (*(cur_req->status) == NC_NOERR)
                        /* keep the first error */
                        *(cur_req->status) = status;
                } else if ( ncmpii_need_swap(varp->type, ptype) ) {
                    ncmpii_in_swapn(cur_req->cbuf, nelems,
                                    ncmpix_len_nctype(varp->type));
                }
                if (!cur_req->iscontig_of_ptypes) {
                    /* derived datatype: unpack from contiguous buffer */
                    status = ncmpii_data_repack(cur_req->cbuf, cnelems, ptype,
                                                cur_req->buf, cur_req->bufcount,
                                                cur_req->datatype);
                    if (*(cur_req->status) == NC_NOERR)
                        /* keep the first error */
                        *(cur_req->status) = status;
                    if (status != NC_NOERR) {
                        FREE_REQUEST(cur_req)
                        NCI_Free(cur_req);
                        return status;
                    }
                }
            }
        }

        /* free space allocated for the request objects */
        FREE_REQUEST(cur_req)
        pre_req = cur_req;
        cur_req = cur_req->next;
        NCI_Free(pre_req);
    }

    return ((warning != NC_NOERR) ? warning : status);
}

/*----< ncmpii_wait_getput() >-----------------------------------------------*/
static int
ncmpii_wait_getput(NC     *ncp,
                   int     num_reqs,  /* # requests incuding subrequests */
                   NC_req *req_head,  /* linked list not include subrequests */
                   int     rw_flag,   /* WRITE_REQ or READ_REQ */
                   int     io_method) /* COLL_IO or INDEP_IO */
{
    int i, j, k, status, ngroups, max_ngroups, *group_index, *statuses;
    void **bufs;
    MPI_Offset **starts, **counts, **strides, *nbytes;
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

    if (io_method == COLL_IO)
        /* I/O will be completed in max_ngroups number of collective I/O */
        MPI_Allreduce(&ngroups, &max_ngroups, 1, MPI_INT, MPI_MAX,
                      ncp->nciop->comm);
    else /* INDEP_IO */
        max_ngroups = ngroups;

    if (num_reqs > 0) {
        starts   = (MPI_Offset**) NCI_Malloc(3 * num_reqs * sizeof(MPI_Offset*));
        counts   = starts + num_reqs;
        strides  = counts + num_reqs;
        nbytes   = (MPI_Offset*) NCI_Malloc(num_reqs * sizeof(MPI_Offset));
        varps    = (NC_var**)    NCI_Malloc(num_reqs * sizeof(NC_var*));
        bufs     = (void**)      NCI_Malloc(num_reqs * sizeof(void*));
        statuses = (int*)        NCI_Malloc(num_reqs * sizeof(int));
    }

    for (i=0; i<ngroups; i++) {
        int i_num_reqs = (i == ngroups-1) ? num_reqs : group_index[i+1];
        i_num_reqs -= group_index[i]; /* number of requests in group i */

        for (j=0; j<i_num_reqs; j++) {
            k = group_index[i] + j;
            varps[j]    = reqs[k].varp;
            starts[j]   = reqs[k].start;
            counts[j]   = reqs[k].count;
            strides[j]  = reqs[k].stride;
            bufs[j]     = reqs[k].xbuf;
            nbytes[j]   = reqs[k].nelems * varps[j]->xsz;
            statuses[j] = NC_NOERR;
        }

        status = ncmpii_mgetput(ncp, i_num_reqs, varps, starts, counts,
                                strides, bufs, nbytes, statuses, rw_flag,
                                io_method);
	/* XXX: hey, we set status here and then never use it... */

        /* update status */
        for (j=0; j<i_num_reqs; j++)
            if (*(reqs[group_index[i] + j].status) == NC_NOERR)
                *(reqs[group_index[i] + j].status) = statuses[j];
    }

    /* All processes must participate all max_ngroups collective calls */
    for (i=ngroups; i<max_ngroups; i++)
        ncmpii_mgetput(ncp, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       rw_flag, io_method);

    if (num_reqs > 0) {
        NCI_Free(statuses);
        NCI_Free(bufs);
        NCI_Free(varps);
        NCI_Free(nbytes);
        NCI_Free(starts);
        NCI_Free(group_index);
        NCI_Free(reqs);
    }

    return NC_NOERR;
}

/*----< ncmpii_mgetput() >-----------------------------------------------*/
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
    int i, len, status=NC_NOERR, warning=NC_NOERR, mpireturn;
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
    if (status != NC_NOERR)
        return status;

    /* create the I/O buffer derived data type */
    if (num_reqs > 0) {
        int *blocklengths = (int*) NCI_Malloc(num_reqs * sizeof(int));
        MPI_Aint *disps = (MPI_Aint*) NCI_Malloc(num_reqs*sizeof(MPI_Aint));
        MPI_Aint a0, ai;

        disps[0] = 0;
        MPI_Get_address(bufs[0], &a0);
        blocklengths[0] = nbytes[0];
        for (i=1; i<num_reqs; i++) {/*loop for multi-variables*/
/* wkliao: type convert from MPI_Offset nbytes[i] to int blocklengths[i]
 *         Can we do someting smarter here ? */
            blocklengths[i] = nbytes[i];
            MPI_Get_address(bufs[i], &ai);
            disps[i] = ai - a0;
        }
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
            CHECK_MPI_ERROR("MPI_File_read_all", NC_EREAD)
        } else {
            mpireturn = MPI_File_read(fh, buf, len, buf_type, &mpistatus);
            CHECK_MPI_ERROR("MPI_File_read", NC_EREAD)
        }
    } else { /* WRITE_REQ */
        if (io_method == COLL_IO) {
            mpireturn = MPI_File_write_all(fh, buf, len, buf_type, &mpistatus);
            CHECK_MPI_ERROR("MPI_File_write_all", NC_EWRITE)
        } else {
            mpireturn = MPI_File_write(fh, buf, len, buf_type, &mpistatus);
            CHECK_MPI_ERROR("MPI_File_write", NC_EWRITE)
        }
    }

    if (buf_type != MPI_BYTE)
        MPI_Type_free(&buf_type);

    /* reset fileview so the entire file is visible again */
    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    if (rw_flag == WRITE_REQ) {
        for (i=0; i<num_reqs; i++) {
            if (IS_RECVAR(varps[i])) {
                /* update the number of records in NC */
                MPI_Offset newnumrecs = starts[i][0] + counts[i][0];
                ncmpii_update_numrecs(ncp, newnumrecs); /* this is collective */
            }
        }
    }

    return ((warning != NC_NOERR) ? warning : status);
}
