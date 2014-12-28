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

#include <mpi.h>

#include "nc.h"
#include "ncx.h"
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

static int ncmpii_mgetput(NC *ncp, int num_reqs, NC_req *reqs, int rw_flag,
                          int io_method);

/*----< ncmpii_getput_zero_req() >--------------------------------------------*/
/* This function is called when this process has zero-length I/O request and
 * must participate all the MPI collective calls involved in the collective
 * APIs and wait_all(), which include setting fileview, collective read/write,
 * another setting fileview.
 *
 * This function is collective.
 */
static int
ncmpii_getput_zero_req(NC  *ncp,
                       int  rw_flag)      /* WRITE_REQ or READ_REQ */
{
    int err, mpireturn, status=NC_NOERR;
    MPI_Status mpistatus;
    MPI_File fh;

    fh = ncp->nciop->collective_fh;

    TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                MPI_INFO_NULL);

    if (rw_flag == READ_REQ) {
        TRACE_IO(MPI_File_read_all)(fh, NULL, 0, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_File_read_all");
            status = (err == NC_EFILE) ? NC_EREAD : err;
        }
    } else { /* WRITE_REQ */
        TRACE_IO(MPI_File_write_all)(fh, NULL, 0, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_File_write_all");
            status = (err == NC_EFILE) ? NC_EWRITE : err;
        }
    }

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                 MPI_INFO_NULL);
     */

    return status;
}

/*----< ncmpii_abuf_coalesce() >---------------------------------------------*/
/* this function should be called after all bput requests have been served */
static int
ncmpii_abuf_coalesce(NC *ncp)
{
    int i;

    if (ncp->abuf == NULL) return NC_NOERR;

    i = ncp->abuf->tail - 1;
    /* tail is always pointing to the last (empty) entry of occupy_table[] */

    /* coalesce the freed entries backwardly from the tail */
    while (i >= 0) {
        if (ncp->abuf->occupy_table[i].is_used == 0) {
            ncp->abuf->size_used -= ncp->abuf->occupy_table[i].req_size;
            i--;
        }
        else break;
    }
    ncp->abuf->tail = i + 1;
    /* This may not be ideal, as we stop at the last one that is yet to be
     * freed. There may be some freed entries before this yet-to-be-freed
     * one, but we would like to keep the available space as contiguous as
     * possible. Maybe some smart approach can be considered here.
     */

    return NC_NOERR;
}

#define FREE_REQUEST(req) {                                                   \
    if (req->abuf_index >= 0)                                                 \
        ncp->abuf->occupy_table[req->abuf_index].is_used = 0; /* mark free */ \
    else if (req->xbuf != NULL && req->xbuf != req->buf)                      \
        NCI_Free(req->xbuf);                                                  \
    req->xbuf = NULL;                                                         \
                                                                              \
    for (j=0; j<req->num_subreqs; j++)                                        \
        NCI_Free(req->subreqs[j].start);                                      \
                                                                              \
    if (req->num_subreqs > 0)                                                 \
        NCI_Free(req->subreqs);                                               \
    NCI_Free(req->start);                                                     \
}

/*----< ncmpi_cancel() >------------------------------------------------------*/
int
ncmpi_cancel(int  ncid,
             int  num_req,
             int *req_ids,  /* [num_req] */
             int *statuses) /* [num_req] */
{
    int status;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    if (NC_indef(ncp)) return NC_EINDEFINE;

    return ncmpii_cancel(ncp, num_req, req_ids, statuses);
}

/*----< ncmpii_cancel() >-----------------------------------------------------*/
int
ncmpii_cancel(NC  *ncp,
              int  num_req,
              int *req_ids,  /* [num_req] */
              int *statuses) /* [num_req] can be NULL (ignore status) */
{
    int i, j, status=NC_NOERR;
    NC_req *pre_req, *cur_req;

    if (num_req == 0)
        return NC_NOERR;

    /* collect the requests from the linked list */
    for (i=0; i<num_req; i++) {
        int found_id=-1;
        if (statuses != NULL) statuses[i] = NC_NOERR;

        if (req_ids[i] == NC_REQ_NULL)
            continue;

        if (ncp->head == NULL) {
            if (statuses != NULL) statuses[i] = NC_EINVAL_REQUEST;
            if (status == NC_NOERR) status = NC_EINVAL_REQUEST;
            continue;
        }

        pre_req = cur_req = ncp->head;
        while (cur_req != NULL) {
            /* there may be more than one node with the same ID */
            if (cur_req->id == req_ids[i]) {
                if (cur_req == ncp->head) { /* move ncp->head ahead */
                    ncp->head = cur_req->next;
                    pre_req = ncp->head;
                }
                else /* move pre_req ahead */
                    pre_req->next = cur_req->next;

                FREE_REQUEST(cur_req)
                NCI_Free(cur_req);
                found_id = req_ids[i];
                /* move cur_req ahead */
                cur_req = (pre_req == ncp->head) ? pre_req : pre_req->next;
            }
            else if (found_id >= 0) break;
            else {
                if (cur_req == pre_req) cur_req = pre_req->next;
                else {
                    pre_req = cur_req;
                    cur_req = pre_req->next;
                }
            }
        }

        if (found_id == -1) { /* no such request ID */
            if (statuses != NULL) statuses[i] = NC_EINVAL_REQUEST;
            /* retain the first error status */
            if (status == NC_NOERR) status = NC_EINVAL_REQUEST;
        }
    }

    /* make sure ncp->tail pointing to the tail */
    ncp->tail = ncp->head;
    while (ncp->tail != NULL && ncp->tail->next != NULL)
        ncp->tail = ncp->tail->next;

    return status;
}

/*----< ncmpi_inq_nreqs() >---------------------------------------------------*/
int
ncmpi_inq_nreqs(int  ncid,
                int *nreqs) /* OUT: number of pending requests */
{
    int     status, prev_id;
    NC     *ncp;
    NC_req *req_ptr;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    req_ptr = ncp->head;
    *nreqs = 0;
    prev_id = NC_REQ_NULL;
    while (req_ptr != NULL) {
        /* some requests in the linked list may share the same ID, they were
         * corresponding to the subrequests made from iput/iget/bput varn APIs.
         */
        if (req_ptr->id != prev_id) {
            prev_id = req_ptr->id;
            (*nreqs)++;
        }
        req_ptr = req_ptr->next;
    }

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
    int  status=NC_NOERR;
    NC  *ncp;

    if (num_reqs == 0) return NC_NOERR;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

#ifdef ENABLE_NONBLOCKING
    return ncmpii_wait(ncp, INDEP_IO, num_reqs, req_ids, statuses);
#else
    int i, err;
    for (i=0; i<num_reqs; i++) { /* serve one request at a time */
        err = ncmpii_wait(ncp, INDEP_IO, 1, &req_ids[i], &statuses[i]);
        if (status == NC_NOERR) status = err;
    }
    return status; /* return the first error encountered */
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
    int  status=NC_NOERR;
    NC  *ncp;
    /* the following line CANNOT be added, because ncmpi_wait_all() is a
     * collective call, all processes must participate some MPI collective
     * operations used later on.
     */
    /* if (num_reqs == 0) return NC_NOERR; */

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        /* must return the error now, parallel program might hang */
        return status;

#ifdef ENABLE_NONBLOCKING
    return ncmpii_wait(ncp, COLL_IO, num_reqs, req_ids, statuses);
#else
    int  i, err;

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

    return status; /* return the first error encountered, if there is any */
#endif
}

/*----< ncmpii_concatenate_datatypes() >--------------------------------------*/
static int
ncmpii_concatenate_datatypes(int           num,
                             int          *blocklens,     /* IN: [num] */
                             MPI_Offset   *displacements, /* IN: [num] */
                             MPI_Datatype *dtypes,        /* IN: [num] */
                             MPI_Datatype *datatype)      /* OUT: */
{
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
    int i;
#endif
    int mpireturn, status=NC_NOERR;
    MPI_Aint *addrs;

    *datatype = MPI_BYTE;

    if (num <= 0) return NC_NOERR;

    /* on most 32 bit systems, MPI_Aint and MPI_Offset are different sizes.
     * Possible that on those platforms some of the beginning offsets of
     * these variables in the dataset won't fit into the aint used by
     * MPI_Type_create_struct.  Minor optimization: we don't need to do any
     * of this if MPI_Aint and MPI_Offset are the same size  */

    /* at the configure time, size of MPI_Offset and MPI_Aint are checked */
#if SIZEOF_MPI_AINT == SIZEOF_MPI_OFFSET
    addrs = (MPI_Aint*) displacements; /* cast ok: types same size */
#else
    /* if (sizeof(MPI_Offset) != sizeof(MPI_Aint)) */
    addrs = (MPI_Aint *) NCI_Malloc(num * sizeof(MPI_Aint));
    for (i=0; i<num; i++) {
        addrs[i] = displacements[i];
        if (displacements[i] != addrs[i]) {
            NCI_Free(addrs);
            return NC_EAINT_TOO_SMALL;
        }
    }
#endif

#ifdef HAVE_MPI_TYPE_CREATE_STRUCT
    mpireturn = MPI_Type_create_struct(num, blocklens, addrs, dtypes, datatype);
#else
    mpireturn = MPI_Type_struct(num, blocklens, addrs, dtypes, datatype);
#endif
    if (mpireturn != MPI_SUCCESS)
        status = ncmpii_handle_error(mpireturn, "MPI_Type_create_struct");
    else
        MPI_Type_commit(datatype);

#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
    NCI_Free(addrs);
#endif

    return status;
}

/*----< ncmpii_construct_filetypes() >----------------------------------------*/
/* concatenate the requests into a single MPI derived filetype */
static int
ncmpii_construct_filetypes(NC           *ncp,
                           int           num_reqs,
                           NC_req       *reqs,      /* [num_reqs] */
                           int           rw_flag,
                           MPI_Datatype *filetype)  /* OUT */
{
    int i, j, err, status=NC_NOERR, *blocklens;
    MPI_Datatype *ftypes;
    MPI_Offset *displacements;

    if (num_reqs <= 0) { /* for participating collective call */
        *filetype = MPI_BYTE;
        return NC_NOERR;;
    }

    /* hereinafter, num_reqs > 0 */
    blocklens     = (int*)          NCI_Malloc(num_reqs * sizeof(int));
    displacements = (MPI_Offset*)   NCI_Malloc(num_reqs * sizeof(MPI_Offset));
    ftypes        = (MPI_Datatype*) NCI_Malloc(num_reqs * sizeof(MPI_Datatype));

    /* create a filetype for each request */
    int last_contig_req = -1; /* index of the last contiguous request */
    j = 0;                    /* index of last valid ftypes */
    for (i=0; i<num_reqs; i++, j++) {
        int is_filetype_contig;
        ftypes[j] = MPI_BYTE; /* in case the call below failed */
        err = ncmpii_vars_create_filetype(ncp,
                                          reqs[i].varp,
                                          reqs[i].start,
                                          reqs[i].count,
                                          reqs[i].stride,
                                          rw_flag,
                                          &blocklens[j],
                                          &displacements[j], /* to offset 0 */
                                          &ftypes[j],
                                          &is_filetype_contig);
        if (*reqs[i].status == NC_NOERR) *reqs[i].status = err;
        if (status == NC_NOERR) status = err; /* report the first error */

        if (is_filetype_contig) {
            if (last_contig_req >= 0 &&
                displacements[j] - displacements[last_contig_req] ==
                blocklens[last_contig_req]) {
                blocklens[last_contig_req] += blocklens[j];
                j--;
            }
            else last_contig_req = j;
        }
        else last_contig_req = -1;
    }
    /* j is the new num_reqs */
    num_reqs = j;

    if (status != NC_NOERR) {
        /* even if error occurs, we still must participate the collective
           call to MPI_File_set_view() */
        *filetype = MPI_BYTE;
    }
    else {
        /* all ftypes[] created fine, now concatenate all ftypes[] */
        err = ncmpii_concatenate_datatypes(num_reqs,
                                           blocklens,
                                           displacements,
                                           ftypes,
                                           filetype);
        if (err != NC_NOERR) *filetype = MPI_BYTE;
        if (status == NC_NOERR) status = err; /* report the first error */
    }

    for (i=0; i<num_reqs; i++) {
        if (ftypes[i] != MPI_BYTE)
            MPI_Type_free(&ftypes[i]);
    }
    NCI_Free(ftypes);
    NCI_Free(displacements);
    NCI_Free(blocklens);

    return status;
}

/*----< ncmpii_construct_buffertypes() >--------------------------------------*/
/* the input requests, reqs[], are non-interleaving requests */
static int
ncmpii_construct_buffertypes(int           num_reqs,
                             NC_req       *reqs,         /* [num_reqs] */
                             MPI_Datatype *buffer_type)  /* OUT */
{
    int i, j, status=NC_NOERR, mpireturn;

    *buffer_type = MPI_BYTE;
    if (num_reqs == 0) return NC_NOERR;

    /* create the I/O buffer derived data type */
    int *blocklengths = (int*) NCI_Malloc(num_reqs * sizeof(int));
    MPI_Aint *disps = (MPI_Aint*) NCI_Malloc(num_reqs*sizeof(MPI_Aint));
    MPI_Aint a0, ai;
    MPI_Offset int8;

    /* process only valid requests */
    for (i=0, j=0; i<num_reqs; i++) {
        /* check int overflow */
        int8  = reqs[i].bnelems * reqs[i].varp->xsz;
        blocklengths[j] = int8;
        if (int8 != blocklengths[j]) { /* skip this request */
            status = NC_EINTOVERFLOW;
            continue;
        }
#ifdef HAVE_MPI_GET_ADDRESS
        MPI_Get_address(reqs[i].xbuf, &ai);
#else
        MPI_Address(reqs[i].xbuf, &ai);
#endif
        if (j == 0) a0 = ai;
        disps[j] = ai - a0;
        j++;
    }
    /* update num_reqs to number of valid requests */
    num_reqs = j;

    if (num_reqs > 0) {
        /* concatenate buffer addresses into a single buffer type */
#ifdef HAVE_MPI_TYPE_CREATE_HINDEXED
        mpireturn = MPI_Type_create_hindexed(num_reqs, blocklengths, disps,
                                             MPI_BYTE, buffer_type);
#else
        mpireturn = MPI_Type_hindexed(num_reqs, blocklengths, disps, MPI_BYTE,
                                      buffer_type);
#endif
        if (mpireturn != MPI_SUCCESS) {
            int err = ncmpii_handle_error(mpireturn, "MPI_Type_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else
            MPI_Type_commit(buffer_type);
    }
    NCI_Free(disps);
    NCI_Free(blocklengths);

    return status;
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
    int i, j, err=NC_NOERR, status=NC_NOERR, fatal_err=NC_NOERR;
    int do_read, do_write, num_w_reqs, num_r_reqs;
    NC_req *pre_req, *cur_req;
    NC_req *w_req_head, *w_req_tail, *r_req_head, *r_req_tail;

    /* check if it is in define mode */
    if (NC_indef(ncp)) fatal_err = NC_EINDEFINE;

    /* check whether collective or independent mode */
    if (fatal_err == NC_NOERR) {
        if (io_method == INDEP_IO)
            fatal_err = ncmpii_check_mpifh(ncp, 0);
        else if (io_method == COLL_IO)
            fatal_err = ncmpii_check_mpifh(ncp, 1);
    }
    if (fatal_err != NC_NOERR) {
        if (io_method == COLL_IO)
            num_reqs = 0; /* to skip the following 2 num_reqs loops and jump to
                             MPI_Allreduce for collectively return */
        else /* INDEP_IO */
            return fatal_err;
    }

    /* Note: 1) it is illegal num_reqs is larger than the linked list size
             2) request ids must be distinct
     */
    j = 0;
    for (i=0; i<num_reqs; i++) {
        statuses[i] = NC_NOERR;
        if (req_ids[i] == NC_REQ_NULL) /* skip NULL request */
            continue;
        j++;
    }
    /* j now is the true number of non-NULL requests, but some requests
       may still be invalid, i.e. no such request IDs */

    if (io_method == INDEP_IO && j == 0)
        return NC_NOERR;
    /* For collective APIs, even though some processes may have zero-length
       requests, they must still participate the collective call. Hence,
       only independent APIs stop here if request is of zero length.
     */

    w_req_head = w_req_tail = NULL;
    r_req_head = r_req_tail = NULL;
    num_w_reqs = num_r_reqs = 0;

    /* extract the requests from the linked list into a new linked list.
       In the meantime coalesce the linked list */

    for (i=0; i<num_reqs; i++) {
        int found_id=-1;
        if (req_ids[i] == NC_REQ_NULL) /* skip zero-size request */
            continue;

        if (ncp->head == NULL) { /* this reqeust is invalid */
            statuses[i] = NC_EINVAL_REQUEST;
            /* retain the first error status */
            if (status == NC_NOERR)
                status = NC_EINVAL_REQUEST;
            continue;
            /* cannot break loop i, must continue to set the error status */
        }

        /* find req_ids[i] from the request linked list */
        pre_req = cur_req = ncp->head;
        while (cur_req != NULL) {
            /* there may be more than one linked node with the same ID */
            if (cur_req->id == req_ids[i]) { /* found */
                /* point all subrequests' statuses to status */
                cur_req->status = statuses + i;
                for (j=0; j<cur_req->num_subreqs; j++)
                    cur_req->subreqs[j].status = cur_req->status;

                /* remove cur_req from the ncp->head linked list */
                if (cur_req == ncp->head) { /* move cur_req to next */
                    ncp->head = cur_req->next;
                    pre_req = ncp->head;
                }
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
                       its subrequests (one for each individual record) */
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
                       its subrequests (one for each individual record) */
                }
                found_id = req_ids[i]; /* indicating previous found ID */
                cur_req = (pre_req == ncp->head) ? pre_req : pre_req->next;
            }
            else if (found_id >= 0) {
                /* same request IDs are stored contiguously in the linked list
                 * so if it is already found and this next request has a
                 * different ID, no need to contibue check the rest list */
                break; /* while loop */
            }
            else {
            /* move on to next request ID */
                if (cur_req == pre_req) cur_req = pre_req->next;
                else {
                    pre_req = cur_req;
                    cur_req = pre_req->next;
                }
            }
        } /* done with searching the entire linked list for this request ID */

        if (found_id == -1) { /* no such request ID */
            statuses[i] = NC_EINVAL_REQUEST;
            /* retain the first error status */
            if (status == NC_NOERR)
                status = NC_EINVAL_REQUEST;
        }
    }
    /* make sure ncp->tail pointing to the tail */
    ncp->tail = ncp->head;
    while (ncp->tail != NULL && ncp->tail->next != NULL)
        ncp->tail = ncp->tail->next;

    if (io_method == COLL_IO) {
        int mpireturn;
        int io_req[3], do_io[3];  /* [0]: read [1]: write [2]: error */
        io_req[0] = num_r_reqs;
        io_req[1] = num_w_reqs;
        io_req[2] = -fatal_err;   /* all NC errors are negative */
        TRACE_COMM(MPI_Allreduce)(io_req, do_io, 3, MPI_INT, MPI_MAX,
                                  ncp->nciop->comm);

        /* if fatal_err occurs, return the API collectively */
        if (do_io[2] != -NC_NOERR) return fatal_err;

        /* make sure if at least one process has a non-zero request, all
           processes participate the collective read/write */
        do_read  = do_io[0];
        do_write = do_io[1];
    }
    else {
        do_read  = num_r_reqs;
        do_write = num_w_reqs;
    }

    /* carry out writes and reads separately (writes first) */
    if (do_write > 0)
        err = ncmpii_wait_getput(ncp, num_w_reqs, w_req_head,
                                 WRITE_REQ, io_method);

    if (do_read > 0)
        err = ncmpii_wait_getput(ncp, num_r_reqs, r_req_head,
                                 READ_REQ, io_method);

    /* retain the first error status */
    if (status == NC_NOERR) status = err;

    /* post-IO data processing: may need byte-swap user write buf, or
                                byte-swap and type convert user read buf */

    if (w_req_head != NULL) {
        cur_req = w_req_head;
        while (cur_req != NULL) {
            /* must byte-swap the user buffer back to its original Endianess
               only when the buffer itself has been byte-swapped before,
               i.e. NOT buftype_is_contig && NOT ncmpii_need_convert() &&
               ncmpii_need_swap()
             */
            if (cur_req->need_swap_back_buf)
                ncmpii_in_swapn(cur_req->buf, cur_req->bnelems,
                                ncmpix_len_nctype(cur_req->varp->type));
            cur_req = cur_req->next;
        }
        cur_req = w_req_head;
        while (cur_req != NULL) {
            /* free temp space allocated for iput/bput varn requests
             * tmpBuf is used only by nonblocking varn APIs */
            if (cur_req->tmpBuf != NULL && cur_req->abuf_index == -1)
                NCI_Free(cur_req->tmpBuf);

            FREE_REQUEST(cur_req)
            pre_req = cur_req;
            cur_req = cur_req->next;
            NCI_Free(pre_req);
        }
        /* once the bput requests are served, we reclaim the space and try
         * coalesce the freed space */
        ncmpii_abuf_coalesce(ncp);
    }

    if (r_req_head != NULL) {
        cur_req = r_req_head;
        while (cur_req != NULL) {
            NC_var *varp = cur_req->varp;

            /* now, xbuf contains the data read from the file.
             * It needs to be type-converted + byte-swapped to cbuf
             */
            void *cbuf, *lbuf;
            int el_size, position, insize;

            MPI_Type_size(cur_req->ptype, &el_size);
            insize = cur_req->bnelems * el_size;

            if (ncmpii_need_convert(varp->type, cur_req->ptype)) {
                if (cur_req->imaptype != MPI_DATATYPE_NULL || !cur_req->buftype_is_contig)
                    cbuf = NCI_Malloc(insize);
                else
                    cbuf = cur_req->buf;

                /* type convert + byte swap from xbuf to cbuf */
                DATATYPE_GET_CONVERT(varp->type, cur_req->xbuf, cbuf,
                                     cur_req->bnelems, cur_req->ptype, err)
                /* keep the first error */
                if (*cur_req->status == NC_NOERR) *cur_req->status = err;
                if (status == NC_NOERR) status = err;
            } else {
                if (ncmpii_need_swap(varp->type, cur_req->ptype))
                    ncmpii_in_swapn(cur_req->xbuf, cur_req->bnelems,
                                    ncmpix_len_nctype(varp->type));
                cbuf = cur_req->xbuf;
            }

            if (cur_req->imaptype != MPI_DATATYPE_NULL) { /* a true varm */
                if (cur_req->buftype_is_contig)
                    lbuf = cur_req->buf;
                else
                    lbuf = NCI_Malloc(insize);

                /* unpack cbuf to lbuf based on imaptype */
                position = 0;
                MPI_Unpack(cbuf, insize, &position, lbuf, 1, cur_req->imaptype,
                           MPI_COMM_SELF);
                MPI_Type_free(&cur_req->imaptype);

                /* cbuf is no longer needed
                 * for a true varm call, cbuf cannot be == cur_req->buf */
                if (cbuf != cur_req->xbuf) NCI_Free(cbuf);
                cbuf = NULL;
            } else { /* get_vars */
                lbuf = cbuf;
            }

            if (!cur_req->buftype_is_contig) {
                /* unpack lbuf to buf based on buftype */
                position = 0;
                MPI_Unpack(lbuf, insize, &position, cur_req->buf,
                           cur_req->bufcount, cur_req->buftype, MPI_COMM_SELF);
            }
            /* lbuf is no longer needed */
            if (lbuf != cur_req->buf && lbuf != cur_req->xbuf)
                NCI_Free(lbuf);

            cur_req = cur_req->next;
        }

        cur_req = r_req_head;
        while (cur_req != NULL) {
            /* free space allocated for the request objects
             * tmpBuf is used only by nonblocking varn APIs */
            if (cur_req->tmpBuf != NULL) {
                int position=0, insize;

                /* unpack tmpBuf to userBuf and free tmpBuf
                 * Note this unpack must wait for all above unpacks are done
                 * because cur_req->buf may be part of cur_req->userBuf
                 */
                MPI_Type_size(cur_req->buftype, &insize);
                insize *= cur_req->bufcount;
                MPI_Unpack(cur_req->tmpBuf, insize, &position,
                           cur_req->userBuf, cur_req->bufcount,
                           cur_req->buftype, MPI_COMM_SELF);
                NCI_Free(cur_req->tmpBuf);
                MPI_Type_free(&cur_req->buftype);
            }
            FREE_REQUEST(cur_req)
            pre_req = cur_req;
            cur_req = cur_req->next;
            NCI_Free(pre_req);
        }
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
                      void       **buf,     /* OUT: 1st I/O buf addr */
                      MPI_Offset  *nsegs,   /* OUT: no. off-len pairs */
                      off_len    **segs)    /* OUT: [*nsegs] */
{
    int i, j, status=NC_NOERR, ndims, is_recvar;
    MPI_Offset  nseg, *start, *count, *shape, *stride;
    MPI_Aint addr, buf_addr;

    *nsegs = 0;    /* total number of offset-length pairs */
    *segs  = NULL; /* array of offset-length pairs */

    /* note invalid requests have been removed in ncmpii_wait_getput() */
    *buf = reqs[0].xbuf; /* I/O buffer of first request */

    /* buf_addr is the buffer address of the first request */
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(reqs[0].xbuf, &buf_addr);
#else
    MPI_Address(reqs[0].xbuf, &buf_addr);
#endif

    /* Count the number off-len pairs from reqs[], so we can malloc a
     * contiguous memory space for storing off-len pairs
     */
    for (i=0; i<num_reqs; i++) {
        is_recvar = IS_RECVAR(reqs[i].varp);

        /* for record variable, each reqs[] is within a record */
        ndims  = (is_recvar) ? reqs[i].varp->ndims - 1 : reqs[i].varp->ndims;
        count  = (is_recvar) ? reqs[i].count + 1       : reqs[i].count;
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
        /* buf_addr is the buffer address of the first valid request */
#ifdef HAVE_MPI_GET_ADDRESS
        MPI_Get_address(reqs[i].xbuf, &addr);
#else
        MPI_Address(reqs[i].xbuf, &addr);
#endif
        addr -= buf_addr,  /* distance to the buf of first req */

        is_recvar = IS_RECVAR(reqs[i].varp);

        /* for record variable, each reqs[] is within a record */
        ndims  = (is_recvar) ? reqs[i].varp->ndims  - 1 : reqs[i].varp->ndims;
        start  = (is_recvar) ? reqs[i].start  + 1       : reqs[i].start;
        count  = (is_recvar) ? reqs[i].count  + 1       : reqs[i].count;
        shape  = (is_recvar) ? reqs[i].varp->shape  + 1 : reqs[i].varp->shape;
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

    /* check if (*segs)[].off are in an increasing order */
    for (i=1; i<*nsegs; i++) {
        if ((*segs)[i-1].off > (*segs)[i].off)
            break;
    }
    if (i < *nsegs) /* not in an increasing order */
        /* sort the off-len array, segs[], in an increasing order */
        qsort(*segs, *nsegs, sizeof(off_len), off_compare);

    /* merge the overlapped requests, skip the overlapped regions for those
     * requests with higher j indices (i.e. requests with lower j indices
     * win the writes to the overlapped regions)
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

    return status;
}

/*----< ncmpii_construct_off_len_type() >-------------------------------------*/
static int
ncmpii_construct_off_len_type(MPI_Offset    nsegs,    /* no. off-len pairs */
                              off_len      *segs,     /* [nsegs] off-en pairs */
                              MPI_Datatype *filetype,
                              MPI_Datatype *buf_type)
{
    int i, j, mpireturn, *blocklengths;
    MPI_Aint   *displacements;
    MPI_Offset  next_off, next_len;

    assert(nsegs > 0);

    /* create the file view MPI derived data type by concatenating the sorted
       offset-length pairs */

    /* For filetype, the segs[].off can be further coalesced. For example,
     * when writing a consecutive columns of a 2D array, even though the I/O
     * buffer addresses may not be able to coalesced, the file offsets on
     * the same row can be coalesced. Thus, first calculate the length of
     * coalesced off-len pairs (the memory space needed for malloc)
     */
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
            /* TODO: take care of 4-byte int overflow problem */
        else {
            j++;
            displacements[j] = segs[i].off;
            blocklengths[j]  = segs[i].len;
        }
    }
    /* j+1 is the coalesced length */

#ifdef HAVE_MPI_TYPE_CREATE_HINDEXED
    mpireturn = MPI_Type_create_hindexed(j+1, blocklengths, displacements,
                                         MPI_BYTE, filetype);
#else
    mpireturn = MPI_Type_hindexed(j+1, blocklengths, displacements, MPI_BYTE,
                                  filetype);
#endif
    if (mpireturn != MPI_SUCCESS) {
        *filetype = MPI_BYTE;
        *buf_type = MPI_BYTE;
        NCI_Free(displacements);
        NCI_Free(blocklengths);
        return ncmpii_handle_error(mpireturn, "MPI_Type_create_hindexed");
    }
    MPI_Type_commit(filetype);
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
            /* TODO: take care of 4-byte int overflow problem */
        else {
            j++;
            displacements[j] = segs[i].buf_addr;
            blocklengths[j]  = segs[i].len;
        }
    }
    /* j+1 is the coalesced length */
#ifdef HAVE_MPI_TYPE_CREATE_HINDEXED
    mpireturn = MPI_Type_create_hindexed(j+1, blocklengths, displacements,
                                         MPI_BYTE, buf_type);
#else
    mpireturn = MPI_Type_hindexed(j+1, blocklengths, displacements, MPI_BYTE,
                                  buf_type);
#endif
    if (mpireturn != MPI_SUCCESS) {
        if (*filetype != MPI_BYTE) MPI_Type_free(filetype);
        *filetype = MPI_BYTE;
        *buf_type = MPI_BYTE;
        NCI_Free(displacements);
        NCI_Free(blocklengths);
        return ncmpii_handle_error(mpireturn, "MPI_Type_create_hindexed");
    }
    MPI_Type_commit(buf_type);
    NCI_Free(displacements);
    NCI_Free(blocklengths);

    return NC_NOERR;
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

/*----< ncmpii_req_aggregation() >--------------------------------------------*/
static int
ncmpii_req_aggregation(NC     *ncp,
                       int     num_reqs,    /* # requests */
                       NC_req *reqs,        /* sorted requests */
                       int     rw_flag,     /* WRITE_REQ or READ_REQ */
                       int     io_method,   /* COLL_IO or INDEP_IO */
                       int     interleaved) /* interleaved in reqs[] */
{
    int i, type, err, status=NC_NOERR, ngroups;
    int *group_index, *group_type;

    if (num_reqs == 0) { /* only COLL_IO can reach here for 0 request */
        assert(io_method == COLL_IO);
        /* simply participate the collective call */
        /* for nonblocking APIs, there is no way for a process to know whether
         * others write to a record variable or not. Hence, we sync number of
         * records for all WRITE_REQ (set the 3rd argument to 1 below)
         */
        return ncmpii_getput_zero_req(ncp, rw_flag);
    }
    if (! interleaved) {
        /* concatenate all filetypes into a single one and do I/O */
        return ncmpii_mgetput(ncp, num_reqs, reqs, rw_flag, io_method);
    }
    /* now some request's aggregate access region is interleaved with other's */

    /* divide the requests into groups.
     * Two types of groups: one contains requests that all are not interleaved
     * and the other contains requests that any 2 consecutive requests are
     * interleaved. All requests will be aggregated into one and carried out
     * by a single MPI-IO call.
     * This approach is because MPI collective I/O requires each process's
     * fileview must contain only monotonic non-decreasing file offsets. Thus
     * if the nonblocking requests interleave with each other (although not
     * overlapp), then we cannot simply concatenate the filetypes of individual
     * requests. This approach flattens the requests of "interleaved" groups
     * into offset-length pairs, sorts, and merges them into an aggregated
     * filetype. Similar for building an aggregated I/O buffer type.
     */

    /* first calculate the number of groups, so group_index and group_type can
       be malloc-ed. Group type: 0 for non-interleaved group and 1 for
       interleaved group.
     */
    ngroups = 1;
    type    = (reqs[0].offset_end > reqs[1].offset_start) ? 1 : 0;
    for (i=1; i<num_reqs-1; i++) {
        if (type == 0 && reqs[i].offset_end > reqs[i+1].offset_start) {
            ngroups++;
            type = 1;
        }
        else if (type == 1 && reqs[i].offset_end <= reqs[i+1].offset_start) {
            type = 0;
            if (i+2 < num_reqs && reqs[i+1].offset_end > reqs[i+2].offset_start)
                type = 1; /* next group is also interleaved */
            ngroups++;
        }
    }

    group_index = (int*) NCI_Malloc((ngroups+1) * sizeof(int));
    group_type  = (int*) NCI_Malloc((ngroups+1) * sizeof(int));

    /* calculate the starting index of each group and determine group type */
    ngroups        = 1;
    type           = (reqs[0].offset_end > reqs[1].offset_start) ? 1 : 0;
    group_index[0] = 0;
    group_type[0]  = type;
    for (i=1; i<num_reqs-1; i++) {
        if (type == 0 &&
            reqs[i].offset_end > reqs[i+1].offset_start) {
            /* reqs[i] starts an interleaved group */
            group_index[ngroups] = i;
            type = 1;
            group_type[ngroups] = type;
            ngroups++;
        }
        else if (type == 1 &&
                 reqs[i].offset_end <= reqs[i+1].offset_start) {
            /* the interleaved group ends with reqs[i] */
            group_index[ngroups] = i+1;
            type = 0;
            if (i+2 < num_reqs && reqs[i+1].offset_end > reqs[i+2].offset_start)
                type = 1; /* next group is also interleaved */
            group_type[ngroups] = type;
            ngroups++;
        }
    }
    group_index[ngroups] = num_reqs; /* to indicate end of groups */

    /* for each group, construct one filetype by concatenating if the group
     * is non-interleaved and by flatten/sort/merge if the group is
     * interleaved. At the end, all ngroups filetypes are concatenated into
     * a single filetype. Similar for constructing buffer types.
     * Then use one collective I/O to commit.
     */

    void         *buf; /* point to starting buffer, used by MPI-IO call */
    int          *f_blocklengths, *b_blocklengths;
    MPI_Aint      b_begin, b_addr, *f_disps, *b_disps;
    MPI_Datatype  filetype, buf_type, *ftypes, *btypes;

    ftypes = (MPI_Datatype*) NCI_Malloc(2 * ngroups * sizeof(MPI_Datatype));
    btypes = ftypes + ngroups;
    f_blocklengths = (int*) NCI_Malloc(2 * ngroups * sizeof(int));
    b_blocklengths = f_blocklengths + ngroups;
    f_disps = (MPI_Aint*) NCI_Malloc(2 * ngroups*sizeof(MPI_Aint));
    b_disps = f_disps + ngroups;

    buf = reqs[0].xbuf; /* the buffer of 1st request */
    b_disps[0] = 0;     /* relative to address of 1st buf */
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(buf, &b_begin);
#else
    MPI_Address(buf, &b_begin);
#endif

    /* for each group, build a filetype and a buffer type in ftypes[i] and
       btypes[i] */
    for (i=0; i<ngroups; i++) {
        NC_req *g_reqs = reqs + group_index[i];
        int     g_num_reqs = group_index[i+1] - group_index[i];
        f_disps[i] = 0;  /* file displacements always to the file offset 0 */

        if (group_type[i] == 0) {
            /* This group contains no interleaved filetypes, so we can
             * simply concatenate filetypes of this group into a single one
             */
            err = ncmpii_construct_filetypes(ncp, g_num_reqs, g_reqs, rw_flag,
                                             &ftypes[i]);
            if (status == NC_NOERR) status = err;
            if (err != NC_NOERR) { /* skip this group */
                ftypes[i] = btypes[i] = MPI_BYTE;
                f_blocklengths[i] = 0;
                continue;
            }
            f_blocklengths[i] = 1;

            /* concatenate buffer types of this group into a single one */
            err = ncmpii_construct_buffertypes(g_num_reqs, g_reqs, &btypes[i]);
            if (status == NC_NOERR) status = err;
            if (err != NC_NOERR) { /* skip this group */
                ftypes[i] = btypes[i] = MPI_BYTE;
                b_blocklengths[i] = 0;
                f_blocklengths[i] = 0;
                continue;
            }
        }
        else { /* this group is interleaved */
            /* flatten the interleaved requests in this group, so interleaved
             * requests can be sorted and merged into a monotonically
             * non-decreasing filetype. For example, multiple nonblocking
             * requests each accessing a single column of a 2D array, that each
             * produces a filetype interleaving with others'.
             *
             * The pitfall of this flattening is the additional memorys
             * requirement, as it will have to break down each request into a
             * list of offset-length pairs, and merge all lists into a sorted
             * list based on their offsets into an increasing order.
             *
             * Be warned! The additional memory requirement for this merging can
             * be more than the I/O data itself. For example, each nonblocking
             * request access a single column of a 2D array of 4-byte integer
             * type. Each off-len pair represents only a 4-byte integer, but the
             * off-len pair itself takes 24 bytes. Additional memory is also
             * required for MPI arguments of displacements and blocklengths.
             */
            MPI_Offset  nsegs=0;   /* number of merged offset-length pairs */
            off_len    *segs=NULL; /* array of the offset-length pairs */
            void       *merged_buf;

            /* merge all requests into sorted offset-length pairs */
            err = ncmpii_merge_requests(ncp, g_num_reqs, g_reqs, &merged_buf,
                                        &nsegs, &segs);
            if (status == NC_NOERR) status = err;
            if (err != NC_NOERR) { /* skip this group */
                ftypes[i] = btypes[i] = MPI_BYTE;
                b_blocklengths[i] = 0;
                f_blocklengths[i] = 0;
                if (segs != NULL) NCI_Free(segs);
                continue;
            }
            assert(nsegs > 0);

            /* sges[] will be used to construct fileview and buffer type */
            err = ncmpii_construct_off_len_type(nsegs, segs, &ftypes[i],
                                                &btypes[i]);
            /* preserve the previous error if there is any */
            if (status == NC_NOERR) status = err;
            NCI_Free(segs);
            if (err != NC_NOERR) { /* skip this group */
                ftypes[i] = btypes[i] = MPI_BYTE;
                b_blocklengths[i] = 0;
                f_blocklengths[i] = 0;
                continue;
            }
            f_blocklengths[i] = 1;
        }

        if (i > 0) {
            /* get the buffer address of the first request in this group */
#ifdef HAVE_MPI_GET_ADDRESS
            MPI_Get_address(g_reqs[0].xbuf, &b_addr);
#else
            MPI_Address(g_reqs[0].xbuf, &b_addr);
#endif
            b_disps[i] = b_addr - b_begin; /* to 1st buffer of 1st group*/
        }
        b_blocklengths[i] = 1;
    }
    NCI_Free(group_index);
    NCI_Free(group_type);

    int mpireturn, buf_len=1;

    if (ngroups == 1) {
        /* use ftypes[0] and btypes[0] directly */
        filetype = ftypes[0];
        buf_type = btypes[0];
    }
    else {
        /* concatenate all ftypes[] to filetype */
#ifdef HAVE_MPI_TYPE_CREATE_STRUCT
        mpireturn = MPI_Type_create_struct(ngroups, f_blocklengths, f_disps,
                                           ftypes, &filetype);
#else
        mpireturn = MPI_Type_struct(ngroups, f_blocklengths, f_disps, ftypes,
                                    &filetype);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_Type_create_struct");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;

            buf_len  = 0; /* skip this request */
            filetype = MPI_BYTE;
        }
        else
            MPI_Type_commit(&filetype);

        for (i=0; i<ngroups; i++) {
            if (ftypes[i] != MPI_BYTE) MPI_Type_free(&ftypes[i]);
        }

        /* concatenate all btypes[] to buf_type */
#ifdef HAVE_MPI_TYPE_CREATE_STRUCT
        mpireturn = MPI_Type_create_struct(ngroups, b_blocklengths, b_disps,
                                           btypes, &buf_type);
#else
        mpireturn = MPI_Type_struct(ngroups, b_blocklengths, b_disps, btypes,
                                    &buf_type);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_handle_error(mpireturn, "MPI_Type_create_struct");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;

            buf_len  = 0; /* skip this request */
            buf_type = MPI_BYTE;
        }
        else
            MPI_Type_commit(&buf_type);

        for (i=0; i<ngroups; i++) {
            if (btypes[i] != MPI_BYTE) MPI_Type_free(&btypes[i]);
        }
    }

    MPI_File fh;
    MPI_Status mpistatus;

    if (io_method == COLL_IO)
        fh = ncp->nciop->collective_fh;
    else
        fh = ncp->nciop->independent_fh;

    /* set the file view */
    MPI_Offset offset=0;
    err = ncmpii_file_set_view(ncp, fh, &offset, filetype);
    if (err != NC_NOERR) {
        buf_len = 0; /* skip this request */
        if (status == NC_NOERR) status = err;
    }

    if (rw_flag == READ_REQ) {
        if (io_method == COLL_IO) {
            TRACE_IO(MPI_File_read_at_all)(fh, offset, buf, buf_len, buf_type,
                                           &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn, "MPI_File_read_at_all");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR)
                    status = (err == NC_EFILE) ? NC_EREAD : err;

            }
        } else {
            TRACE_IO(MPI_File_read_at)(fh, offset, buf, buf_len, buf_type,
                                       &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn, "MPI_File_read_at");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR)
                    status = (err == NC_EFILE) ? NC_EREAD : err;
            }
        }
        if (mpireturn == MPI_SUCCESS) {
            int get_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
            ncp->nciop->get_size += get_size;
        }
    } else { /* WRITE_REQ */
        if (io_method == COLL_IO) {
            TRACE_IO(MPI_File_write_at_all)(fh, offset, buf, buf_len, buf_type,
                                            &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn, "MPI_File_write_at_all");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR)
                    status = (err == NC_EFILE) ? NC_EWRITE : err;
            }
        } else {
            TRACE_IO(MPI_File_write_at)(fh, offset, buf, buf_len, buf_type,
                                        &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn, "MPI_File_write_at");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR)
                    status = (err == NC_EFILE) ? NC_EWRITE : err;
            }
        }
        if (mpireturn == MPI_SUCCESS) {
            int put_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
            ncp->nciop->put_size += put_size;
        }
    }

    if (filetype != MPI_BYTE) MPI_Type_free(&filetype);
    if (buf_type != MPI_BYTE) MPI_Type_free(&buf_type);

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                 MPI_INFO_NULL);
     */

    NCI_Free(ftypes);
    NCI_Free(f_blocklengths);
    NCI_Free(f_disps);

    return status;
}

/*----< ncmpii_check_edge() >-------------------------------------------------*/
static
int ncmpii_check_edge(NC     *ncp,
                      NC_req *req,
                      int     rw_flag)
{
    int is_recvar, err;

    /* check for example out-of-boundary request, and set its status to NC
       error codes. Because some requests are subrequests of the same
       record-variable request, their status pointed to the same memory
       address of the parent request's status. Hece, if one of the
       subrequests is invalid, we mark the parent request invalid and skip it.

       Note that all reqs[].status have been initialized to NC_NOERR
       in ncmpii_wait()
     */
    if (*req->status != NC_NOERR) /* skip this invalid req */
        return *req->status;

    is_recvar = IS_RECVAR(req->varp);

    /* Check for out-of-boundary */
    err = NCedgeck(ncp, req->varp, req->start, req->count);
    if (err != NC_NOERR ||
        (rw_flag == READ_REQ && is_recvar &&
         req->start[0] + req->count[0] > NC_get_numrecs(ncp))) {
        err = NCcoordck(ncp, req->varp, req->start, rw_flag);
        if (err == NC_NOERR) err = NC_EEDGE;
        *req->status = err;
    }
    return err;
}

/*----< ncmpii_wait_getput() >------------------------------------------------*/
static int
ncmpii_wait_getput(NC     *ncp,
                   int     num_reqs,  /* # requests including subrequests */
                   NC_req *req_head,  /* linked list not include subrequests */
                   int     rw_flag,   /* WRITE_REQ or READ_REQ */
                   int     io_method) /* COLL_IO or INDEP_IO */
{
    int i, j, err, status=NC_NOERR, access_interleaved=0;
    NC_req *reqs, *cur_req;
    int (*fcnt)(const void*, const void*);

    /* pack the linked list into an array to be sorted and remove invalid
       requests, such as out-of-boundary ones
     */
    reqs = (NC_req*) NCI_Malloc(num_reqs * sizeof(NC_req));
    i = 0;
    cur_req = req_head;
    while (cur_req != NULL) {
        if (cur_req->num_subreqs == 0) {
            err = ncmpii_check_edge(ncp, cur_req, rw_flag);
            if (status == NC_NOERR) status = err;
            if (err == NC_NOERR) /* add this valid request to reqs[] */
                reqs[i++] = *cur_req;
        }
        else {
            for (j=0; j<cur_req->num_subreqs; j++) {
                err = ncmpii_check_edge(ncp, cur_req->subreqs+j, rw_flag);
                if (status == NC_NOERR) status = err;
                if (err != NC_NOERR) /* skip all subreqs[] */
                    break;
                /* add this valid request to reqs[] */
                reqs[i++] = cur_req->subreqs[j];
            }
        }
        cur_req = cur_req->next;
    }
    /* i is the number of valid requests */
    num_reqs = i;

    /* for independent wait and there is no valid request, we can return now */
    if (num_reqs == 0 && io_method == INDEP_IO) {
        if (reqs != NULL) NCI_Free(reqs);
        return status;
    }

    /* check if reqs[].offset_start are in an increasing order */
    for (i=1; i<num_reqs; i++) {
        if (reqs[i-1].offset_start > reqs[i].offset_start) {
            break;
        }
    }
    if (i < num_reqs) { /* not in an increasing order */
        /* sort reqs[] based on reqs[].offset_start */
        fcnt = (int (*)(const void *, const void *))req_compare;
        qsort(reqs, num_reqs, sizeof(NC_req), fcnt);
    }

    /* check for any interleaved rquests */
    for (i=1; i<num_reqs; i++) {
        if (reqs[i-1].offset_end > reqs[i].offset_start) {
            access_interleaved = 1;
            break;
        }
    }

    /* aggregate requests and carry out the I/O */
    err = ncmpii_req_aggregation(ncp, num_reqs, reqs, rw_flag, io_method,
                                 access_interleaved);
    if (status == NC_NOERR) status = err;

    /* update the number of records if new records are added */
    if (rw_flag == WRITE_REQ) {
        /* Because netCDF allows only one unlimited dimension, find the
         * maximum number of records from all nonblocking requests and
         * update newnumrecs once
         */
        MPI_Offset newnumrecs = ncp->numrecs;
        for (i=0; i<num_reqs; i++) {
            if (*reqs[i].status == NC_NOERR && IS_RECVAR(reqs[i].varp)) {
                newnumrecs = MAX(newnumrecs,
                                 reqs[i].start[0] + reqs[i].count[0]);
            }
        }

        if (io_method == COLL_IO) {
            /* sync numrecs in memory and file. Note that even this process
             * does not write to record variable, others might. Note in
             * ncmpii_sync_numrecs(), new_numrecs is checked against
             * ncp->numrecs and if NC_SHARE is set, MPI_File_sync() will
             * be called.
             */
            err = ncmpii_sync_numrecs(ncp, newnumrecs);
            if (status == NC_NOERR) status = err;
            /* retain the first error if there is any */
        }
        else { /* INDEP_IO */
            if (ncp->numrecs < newnumrecs) {
                ncp->numrecs = newnumrecs;
                set_NC_ndirty(ncp);
                /* delay numrecs sync until end_indep, redef or close */
            }
        }

        if (NC_doFsync(ncp)) { /* NC_SHARE is set */
            int mpireturn;
            if (io_method == INDEP_IO) {
                TRACE_IO(MPI_File_sync)(ncp->nciop->independent_fh);
            }
            else {
                TRACE_IO(MPI_File_sync)(ncp->nciop->collective_fh);
                TRACE_COMM(MPI_Barrier)(ncp->nciop->comm);
            }
        }
    }

    if (reqs != NULL) NCI_Free(reqs);

    return status;
}

/*----< ncmpii_mgetput() >----------------------------------------------------*/
/* all the fileviews in the requests can be concatenated one after another */
static int
ncmpii_mgetput(NC           *ncp,
               int           num_reqs,
               NC_req       *reqs,        /* [num_reqs] */
               int           rw_flag,     /* WRITE_REQ or READ_REQ */
               int           io_method)   /* COLL_IO or INDEP_IO */
{
    int i, j, len=0, status=NC_NOERR, mpireturn, err;
    void *buf=NULL;
    MPI_Status mpistatus;
    MPI_Datatype filetype, buf_type=MPI_BYTE;
    MPI_File fh;

    if (io_method == COLL_IO)
        fh = ncp->nciop->collective_fh;
    else
        fh = ncp->nciop->independent_fh;

    /* construct a MPI file type by concatenating all requests */
    status = ncmpii_construct_filetypes(ncp,num_reqs, reqs, rw_flag, &filetype);
    if (status != NC_NOERR) { /* if failed, skip this request */
        if (io_method == INDEP_IO) return status;

        /* For collective I/O, we still need to participate the successive
           collective calls: read/write/setview */
        num_reqs = 0;
        filetype = MPI_BYTE;
    }

    /* set the MPI-IO fileview */
    MPI_Offset offset=0;
    err = ncmpii_file_set_view(ncp, fh, &offset, filetype);
    if (err != NC_NOERR) {
        num_reqs = 0; /* skip this request */
        if (status == NC_NOERR) status = err;
    }

    if (filetype != MPI_BYTE) MPI_Type_free(&filetype);

    /* now construct buffer datatype */
    if (num_reqs == 0) {
        /* num_reqs == 0, simply participate the collective call */
        buf = NULL;
        len = 0;
    }
    else if (num_reqs == 1) {
        MPI_Offset int8 = reqs[0].bnelems;
        int8 *= reqs[0].varp->xsz;
        len = int8;
        if (len != int8) {
            if (status == NC_NOERR) status = NC_EINTOVERFLOW;
            len = 0; /* skip this request */
        }
        buf = reqs[0].xbuf;
    }
    else if (num_reqs > 1) { /* create the I/O buffer derived data type */
        int *blocklengths = (int*) NCI_Malloc(num_reqs * sizeof(int));
        MPI_Aint *disps = (MPI_Aint*) NCI_Malloc(num_reqs*sizeof(MPI_Aint));
        MPI_Aint a0=0, ai, a_last_contig;
        MPI_Offset int8;

        int last_contig_req = 0; /* index of the last contiguous request */
        buf = NULL;
        /* process only valid requests */
        for (i=0, j=0; i<num_reqs; i++) {
            /* check int overflow */
            int8  = reqs[i].bnelems * reqs[i].varp->xsz;
            blocklengths[j] = int8;
            if (int8 != blocklengths[j]) { /* skip this request */
                if (status == NC_NOERR) status = NC_EINTOVERFLOW;
                blocklengths[j] = 0;
                continue;
            }
#ifdef HAVE_MPI_GET_ADDRESS
            MPI_Get_address(reqs[i].xbuf, &ai);
#else
            MPI_Address(reqs[i].xbuf, &ai);
#endif
            if (j == 0) { /* first valid request */
                a_last_contig = a0 = ai;
                buf = reqs[i].xbuf;
            }
            disps[j] = ai - a0;

            if (ai - a_last_contig == blocklengths[last_contig_req])
                /* user buffer of request j is contiguous from j-1
                 * we concatenate j to j-1 */
                blocklengths[last_contig_req] += blocklengths[j];
            else if (j > 0) {
                /* not contiguous from request last_contig_req */
                last_contig_req++;
                a_last_contig = ai;
                disps[last_contig_req] = ai - a0;
                blocklengths[last_contig_req] = blocklengths[i];
            }
            j++;
        }

        /* last_contig_req is the index of last contiguous request */
        if (last_contig_req == 0) {
            /* user buffers in all requests are acutally contiguous */
            buf_type = MPI_BYTE;
            len = blocklengths[0];
        }
        else {
            /* after possible concatenating the user buffers, the true number
             * of non-contiguous buffers is last_contig_req+1 */
            num_reqs = last_contig_req+1;

            /* concatenate buffer addresses into a single buffer type */
#ifdef HAVE_MPI_TYPE_CREATE_HINDEXED
            mpireturn = MPI_Type_create_hindexed(num_reqs, blocklengths, disps,
                                                 MPI_BYTE, &buf_type);
#else
            mpireturn = MPI_Type_hindexed(num_reqs, blocklengths, disps,
                                          MPI_BYTE, &buf_type);
#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else
                mpireturn = MPI_Type_commit(&buf_type);

            len = 1;
        }
        NCI_Free(disps);
        NCI_Free(blocklengths);
    }
    /* if (buf_type == MPI_BYTE) then the whole buf is contiguous */

    if (rw_flag == READ_REQ) {
        if (io_method == COLL_IO) {
            TRACE_IO(MPI_File_read_at_all)(fh, offset, buf, len, buf_type,
                                           &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn, "MPI_File_read_at_all");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR)
                    status = (err == NC_EFILE) ? NC_EREAD : err;
            }
        } else {
            TRACE_IO(MPI_File_read_at)(fh, offset, buf, len, buf_type,
                                       &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn, "MPI_File_read_at");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR)
                    status = (err == NC_EFILE) ? NC_EREAD : err;
            }
        }
        int get_size;
        MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
        ncp->nciop->get_size += get_size;

    } else { /* WRITE_REQ */
        if (io_method == COLL_IO) {
            TRACE_IO(MPI_File_write_at_all)(fh, offset, buf, len, buf_type,
                                            &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn, "MPI_File_write_at_all");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR)
                    status = (err == NC_EFILE) ? NC_EWRITE : err;
            }
        } else {
            TRACE_IO(MPI_File_write_at)(fh, offset, buf, len, buf_type,
                                        &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_handle_error(mpireturn, "MPI_File_write_at");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR)
                    status = (err == NC_EFILE) ? NC_EWRITE : err;
            }
        }
        int put_size;
        MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
        ncp->nciop->put_size += put_size;
    }

    if (buf_type != MPI_BYTE)
        mpireturn = MPI_Type_free(&buf_type);

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                 MPI_INFO_NULL);
     */

    return status;
}

/*----< ncmpii_set_iget_callback() >-----------------------------------------*/
/* this subroutine is only used by iget_varn API family and when its buftype
 * argument indicated a noncontiguous data type. In this case, tmpBuf is never
 * == userBuf
 *
 * iget_varn() divides userBuf into pieces, each is used in an iget_varm() and
 * stored in a cur_req->buf. We cannot unpack each req->buf individually, as
 * req->buf may need to type converted. So at the end of wait(), we unpack the
 * entire tmpBuf into userBuf.
 */
int ncmpii_set_iget_callback(NC           *ncp,
                             int           reqid,
                             void         *tmpBuf,
                             void         *userBuf,
                             int           userBufCount,
                             MPI_Datatype  userBufType)
{
    NC_req *cur_req;

    if (reqid == NC_REQ_NULL) return NC_NOERR;

    if (ncp->head == NULL) return NC_EINVAL_REQUEST;

    cur_req = ncp->head;
    while (cur_req != NULL) {
        /* find the first linked-list node with reqid, there may be more than
         * one node with the same request ID
         */
        if (cur_req->id == reqid) {
            MPI_Datatype dup_buftype;
            MPI_Type_dup(userBufType, &dup_buftype);

            cur_req->tmpBuf = tmpBuf;
                     /* When not NULL, tmpBuf is an internal buffer allocated
                      * in iget_varn(). At the end of nonblocking wait call,
                      * tmpBuf will be unpacked to userBuf and freed.
                      */
            cur_req->userBuf = userBuf;
                     /* User's buffer. iget_varn() divides userBuf into pieces,
                      * each is used in an iget_varm() and stored in a
                      * cur_req->buf.
                      */
            cur_req->bufcount = userBufCount;
            cur_req->buftype  = dup_buftype;
            break;
        }
        cur_req = cur_req->next;
    }
    if (cur_req == NULL) return NC_EINVAL_REQUEST;

    return NC_NOERR;
}

/*----< ncmpii_set_iput_callback() >-----------------------------------------*/
/* this subroutine is only used by iput_varn API family, to tell wait() to
 * free the temporary buffer at the end.
 */
int ncmpii_set_iput_callback(NC   *ncp,
                             int   reqid,
                             void *tmpPutBuf)
{
    int found_req_id=NC_REQ_NULL;
    NC_req *cur_req;

    if (reqid == NC_REQ_NULL) return NC_NOERR;

    if (ncp->head == NULL) return NC_EINVAL_REQUEST;

    cur_req = ncp->head;
    while (cur_req != NULL) {
        /* there may be more than one linked-list node with the same ID */
        if (cur_req->id == reqid) {
            /* set tmpBuf only on the first node, so it is freed just once */
            if (found_req_id == NC_REQ_NULL) cur_req->tmpBuf = tmpPutBuf;
            found_req_id = reqid;
        }
        if (found_req_id != NC_REQ_NULL && cur_req->id != found_req_id)
            break; /* requests with same ID are in consecutive linked nodes */
        cur_req = cur_req->next;
    }
    if (found_req_id == NC_REQ_NULL) return NC_EINVAL_REQUEST;

    return NC_NOERR;
}
