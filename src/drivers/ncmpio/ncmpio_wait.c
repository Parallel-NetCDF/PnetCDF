/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_wait()     : dispatcher->wait()
 * ncmpi_wait_all() : dispatcher->wait()
 * ncmpi_cancel()   : dispatcher->cancel()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset() */
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/* buffer layers:

        User Level              buf     (user defined buffer of MPI_Datatype)
        MPI Datatype Level      cbuf    (contiguous buffer of itype)
        NetCDF XDR Level        xbuf    (XDR I/O buffer)
*/

/*----< abuf_coalesce() >----------------------------------------------------*/
/* this function should be called after all bput requests have been served */
static int
abuf_coalesce(NC *ncp)
{
    int i;

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

/*----< ncmpio_cancel() >-----------------------------------------------------*/
/* argument num_req can be NC_REQ_ALL, NC_GET_REQ_ALL, NC_PUT_REQ_ALL, or
 * non-negative value. This is an independent subroutine.
 */
int
ncmpio_cancel(void *ncdp,
              int   num_req,
              int  *req_ids,  /* [num_req]: IN/OUT */
              int  *statuses) /* [num_req] can be NULL (ignore status) */
{
    int status=NC_NOERR;
    size_t i, j;
    NC *ncp=(NC*)ncdp;

    if (num_req == 0) return NC_NOERR;
    if (num_req < NC_PUT_REQ_ALL) DEBUG_RETURN_ERROR(NC_EINVAL);

    /* 1.7.0 and after nonblocking APIs can be called in define mode.
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)
    */

    if (num_req == NC_GET_REQ_ALL || num_req == NC_REQ_ALL) {
        /* cancel all pending read requests, ignore req_ids and statuses */
        NC_lead_req *req = ncp->get_lead_list;
        for (i=0; i<ncp->numLeadGetReqs; i++) {
            /* free resource allocated at lead request */
            if (req->imaptype != MPI_DATATYPE_NULL)
                MPI_Type_free(&req->imaptype);
            if (!fIsSet(req->flag, NC_REQ_BUF_TYPE_IS_CONTIG))
                MPI_Type_free(&req->buftype);
            if (req->abuf_index < 0) {
                if (fIsSet(req->flag, NC_REQ_XBUF_TO_BE_FREED))
                    NCI_Free(req->xbuf); /* free xbuf */
            }
            else  /* this is bput request */
                ncp->abuf->occupy_table[req->abuf_index].is_used = 0;
            NCI_Free(req->start);
            req++;
        }
        NCI_Free(ncp->get_list);
        NCI_Free(ncp->get_lead_list);
        ncp->get_list = NULL;
        ncp->get_lead_list = NULL;
        ncp->numGetReqs = 0;
        ncp->numLeadGetReqs = 0;
    }

    if (num_req == NC_PUT_REQ_ALL || num_req == NC_REQ_ALL) {
        /* cancel all pending write requests, ignore req_ids and statuses */
        NC_lead_req *req = ncp->put_lead_list;
        for (i=0; i<ncp->numLeadPutReqs; i++) {
            if (fIsSet(req->flag, NC_REQ_BUF_BYTE_SWAP))
                /* if user buffer is in-place byte-swapped, swap it back */
                ncmpii_in_swapn(req->buf, req->nelems, req->varp->xsz);
            /* free resource allocated at lead request */
            if (req->abuf_index < 0) {
                if (fIsSet(req->flag, NC_REQ_XBUF_TO_BE_FREED))
                    NCI_Free(req->xbuf); /* free xbuf */
            }
            else  /* this is bput request */
                ncp->abuf->occupy_table[req->abuf_index].is_used = 0;
            NCI_Free(req->start);
            req++;
        }
        NCI_Free(ncp->put_list);
        NCI_Free(ncp->put_lead_list);
        ncp->put_list = NULL;
        ncp->put_lead_list = NULL;
        ncp->numPutReqs = 0;
        ncp->numLeadPutReqs = 0;
        if (ncp->abuf != NULL) { /* clear out the attached buffer usage */
            ncp->abuf->tail = 0;
            ncp->abuf->size_used = 0;
        }
    }
    if (num_req < 0) return NC_NOERR; /* done with "ALL" types cancel */

    /* check each request ID from the lead read/write request queues */
    for (i=0; i<num_req; i++) {
        if (statuses != NULL) statuses[i] = NC_NOERR;

        if (req_ids[i] == NC_REQ_NULL) continue; /* skip NULL request */

        if (req_ids[i] & 1) { /* read request (id is an odd number) */
            int found=0;
            NC_lead_req *lead_req=ncp->get_lead_list;
            for (j=0; j<ncp->numLeadGetReqs; j++) {
                if (lead_req->id == NC_REQ_NULL || lead_req->id != req_ids[i]) {
                    lead_req++;
                    continue;
                }
                /* found in the lead queue */
                found = 1;
                /* free resource allocated at lead request */
                if (lead_req->imaptype != MPI_DATATYPE_NULL)
                    MPI_Type_free(&lead_req->imaptype);
                if (!fIsSet(lead_req->flag, NC_REQ_BUF_TYPE_IS_CONTIG))
                    MPI_Type_free(&lead_req->buftype);
                if (fIsSet(lead_req->flag, NC_REQ_XBUF_TO_BE_FREED))
                    NCI_Free(lead_req->xbuf); /* free xbuf */
                NCI_Free(lead_req->start);
                lead_req->id = NC_REQ_NULL; /* marked as freed */
                break;
            }
            if (found) { /* free requests in non-lead queue */
                int k, nonlead_num=lead_req->nonlead_num;
                NC_req *req = ncp->get_list;
                /* remove non-lead requests from get_list */
                k = lead_req->nonlead_off;
                j = lead_req->nonlead_off + lead_req->nonlead_num;
                for (; j<ncp->numGetReqs; j++) {
                    req[k] = req[j]; /* coalesce get_list */
                    req[k].lead_off--;
                    k++;
                }
                ncp->numGetReqs = k;
                req_ids[i] = NC_REQ_NULL;

                /* coalesce get_lead_list */
                j = lead_req - ncp->get_lead_list + 1;
                for (; j<ncp->numLeadGetReqs; j++) {
                    ncp->get_lead_list[j-1] = ncp->get_lead_list[j];
                    ncp->get_lead_list[j-1].nonlead_off -= nonlead_num;
                }
                ncp->numLeadGetReqs--;

                continue; /* loop i, go to next request ID */
            }
            /* else means req_ids[i] is not found in get_list[] */
        }
        else { /* write request (id is an even number) */
            int found=0;
            NC_lead_req *lead_req=ncp->put_lead_list;
            for (j=0; j<ncp->numLeadPutReqs; j++) {
                if (lead_req->id == NC_REQ_NULL || lead_req->id != req_ids[i]) {
                    lead_req++;
                    continue;
                }
                /* found in the lead queue */
                found = 1;
                if (fIsSet(lead_req->flag, NC_REQ_BUF_BYTE_SWAP))
                    /* user buffer has been in-place byte-swapped */
                    ncmpii_in_swapn(lead_req->buf, lead_req->nelems,
                                    lead_req->varp->xsz);
                /* free resource allocated at lead request */
                if (lead_req->abuf_index < 0) {
                    if (fIsSet(lead_req->flag, NC_REQ_XBUF_TO_BE_FREED))
                        NCI_Free(lead_req->xbuf); /* free xbuf */
                }
                else  /* this is bput request */
                    ncp->abuf->occupy_table[lead_req->abuf_index].is_used = 0;
                NCI_Free(lead_req->start);
                lead_req->id = NC_REQ_NULL; /* marked as freed */
                break;
            }
            if (found) {
                int k, nonlead_num=lead_req->nonlead_num;
                NC_req *req = ncp->put_list;
                /* remove non-lead requests from get_list */
                k = lead_req->nonlead_off;
                j = lead_req->nonlead_off + lead_req->nonlead_num;
                for (; j<ncp->numPutReqs; j++) {
                    req[k] = req[j]; /* coalesce put_list */
                    req[k].lead_off--;
                    k++;
                }
                ncp->numPutReqs = k;
                req_ids[i] = NC_REQ_NULL;

                /* coalesce put_lead_list */
                j = lead_req - ncp->put_lead_list + 1;
                for (; j<ncp->numLeadPutReqs; j++) {
                    ncp->put_lead_list[j-1] = ncp->put_lead_list[j];
                    ncp->put_lead_list[j-1].nonlead_off -= nonlead_num;
                }
                ncp->numLeadPutReqs--;

                continue; /* loop i, go to next request ID */
            }
            /* else means req_ids[i] is not found in put_list[] */
        }
        /* no such request ID, if the program reached here */
        if (statuses != NULL) DEBUG_ASSIGN_ERROR(statuses[i], NC_EINVAL_REQUEST)
        /* retain the first error status */
        if (status == NC_NOERR) DEBUG_ASSIGN_ERROR(status, NC_EINVAL_REQUEST)
    }
    if (ncp->abuf != NULL) abuf_coalesce(ncp);

    if (ncp->numLeadGetReqs == 0) {
        NCI_Free(ncp->get_lead_list);
        NCI_Free(ncp->get_list);
        ncp->get_lead_list = NULL;
        ncp->get_list = NULL;
    }

    if (ncp->numLeadPutReqs == 0) {
        NCI_Free(ncp->put_lead_list);
        NCI_Free(ncp->put_list);
        ncp->put_lead_list = NULL;
        ncp->put_list = NULL;
    }

    return status;
}

/*----< extract_reqs() >-----------------------------------------------------*/
/* extract requests from the queues into new queues to be committed.
 * Input value of num_reqs can be NC_REQ_ALL, NC_GET_REQ_ALL, or NC_PUT_REQ_ALL
 */
static int
extract_reqs(NC      *ncp,
             int      num_reqs,
             int     *req_ids,         /* IN: [num_reqs] or NULL */
             int     *statuses,        /* IN: [num_reqs] or NULL */
             int     *num_r_lead_reqs, /* OUT: no. lead get requests */
             int     *num_r_reqs,      /* OUT: no. non-lead get requests */
             NC_req **get_list,        /* OUT: extracted get requests */
             int     *num_w_lead_reqs, /* OUT: no. lead put requests */
             int     *num_w_reqs,      /* OUT: no. non-lead put requests */
             NC_req **put_list)        /* OUT: extracted put requests */
{
    int i, j, status=NC_NOERR;
    NC_req *put_list_ptr, *get_list_ptr;

    *num_r_lead_reqs = 0;
    *num_w_lead_reqs = 0;
    *num_r_reqs      = 0;
    *num_w_reqs      = 0;

    if (num_reqs == NC_PUT_REQ_ALL || num_reqs == NC_REQ_ALL) {
        /* the entire put requests */
        for (i=0; i<ncp->numLeadPutReqs; i++)
            fSet(ncp->put_lead_list[i].flag, NC_REQ_TO_FREE);

        *num_w_lead_reqs = ncp->numLeadPutReqs;
        *num_w_reqs      = ncp->numPutReqs;
        *put_list        = ncp->put_list;
        ncp->numPutReqs  = 0;
        ncp->put_list    = NULL;
    }
    if (num_reqs == NC_GET_REQ_ALL || num_reqs == NC_REQ_ALL) {
        /* the entire get requests */
        for (i=0; i<ncp->numLeadGetReqs; i++)
            fSet(ncp->get_lead_list[i].flag, NC_REQ_TO_FREE);

        *num_r_lead_reqs = ncp->numLeadGetReqs;
        *num_r_reqs      = ncp->numGetReqs;
        *get_list        = ncp->get_list;
        ncp->numGetReqs  = 0;
        ncp->get_list    = NULL;
    }
    if (num_reqs == NC_REQ_ALL || num_reqs == NC_GET_REQ_ALL ||
                                  num_reqs == NC_PUT_REQ_ALL)
        return NC_NOERR;

    if (ncp->numGetReqs == 0 && num_reqs == ncp->numLeadPutReqs) {
        /* this is the same as NC_PUT_REQ_ALL */
        for (i=0; i<num_reqs; i++) req_ids[i] = NC_REQ_NULL;
        if (statuses != NULL) {
            for (i=0; i<ncp->numLeadPutReqs; i++) {
                ncp->put_lead_list[i].status = statuses + i;
                statuses[i] = NC_NOERR;
            }
        }
        for (i=0; i<ncp->numLeadPutReqs; i++)
            fSet(ncp->put_lead_list[i].flag, NC_REQ_TO_FREE);

        *num_w_lead_reqs = ncp->numLeadPutReqs;
        *num_w_reqs      = ncp->numPutReqs;
        *put_list        = ncp->put_list;
        ncp->numPutReqs  = 0;
        ncp->put_list    = NULL;
        return NC_NOERR;
    }
    if (ncp->numPutReqs == 0 && num_reqs == ncp->numLeadGetReqs) {
        /* this is the same as NC_GET_REQ_ALL */
        for (i=0; i<num_reqs; i++) req_ids[i] = NC_REQ_NULL;
        if (statuses != NULL) {
            for (i=0; i<ncp->numLeadGetReqs; i++) {
                ncp->get_lead_list[i].status = statuses + i;
                statuses[i] = NC_NOERR;
            }
        }
        for (i=0; i<ncp->numLeadGetReqs; i++)
            fSet(ncp->get_lead_list[i].flag, NC_REQ_TO_FREE);

        *num_r_lead_reqs = ncp->numLeadGetReqs;
        *num_r_reqs      = ncp->numGetReqs;
        *get_list        = ncp->get_list;
        ncp->numGetReqs  = 0;
        ncp->get_list    = NULL;
        return NC_NOERR;
    }
    if (num_reqs == ncp->numLeadPutReqs + ncp->numLeadGetReqs &&
        statuses == NULL) {
        /* this is the same as NC_REQ_ALL */
        for (i=0; i<num_reqs; i++) req_ids[i] = NC_REQ_NULL;

        for (i=0; i<ncp->numLeadGetReqs; i++)
            fSet(ncp->get_lead_list[i].flag, NC_REQ_TO_FREE);
        *num_w_lead_reqs = ncp->numLeadPutReqs;
        *num_w_reqs      = ncp->numPutReqs;
        *put_list        = ncp->put_list;
        ncp->numPutReqs  = 0;
        ncp->put_list    = NULL;

        for (i=0; i<ncp->numLeadPutReqs; i++)
            fSet(ncp->put_lead_list[i].flag, NC_REQ_TO_FREE);
        *num_r_lead_reqs = ncp->numLeadGetReqs;
        *num_r_reqs      = ncp->numGetReqs;
        *get_list        = ncp->get_list;
        ncp->numGetReqs  = 0;
        ncp->get_list    = NULL;
        return NC_NOERR;
    }

    /* requests are a subset of pending requests */
    for (i=0; i<num_reqs; i++) {
        int found=0;

        if (req_ids[i] == NC_REQ_NULL) { /* skip NULL request */
            if (statuses != NULL) statuses[i] = NC_NOERR;
            continue;
        }

        if (req_ids[i] % 2 == 0) { /* write requests are even numbers */
            for (j=0; j<ncp->numLeadPutReqs; j++) {
                if (fIsSet(ncp->put_lead_list[j].flag, NC_REQ_TO_FREE))
                    continue; /* this request has been processed */
                if (ncp->put_lead_list[j].id == req_ids[i]) {
                    /* make this request to be freed */
                    fSet(ncp->put_lead_list[j].flag, NC_REQ_TO_FREE);
                    if (statuses != NULL) {
                        statuses[i] = NC_NOERR;
                        ncp->put_lead_list[j].status = statuses+i;
                    }
                    (*num_w_lead_reqs)++;
                    *num_w_reqs += ncp->put_lead_list[j].nonlead_num;
                    found = 1;
                    break;
                }
            }
        }
        else { /* get requests are odd numbers */
            for (j=0; j<ncp->numLeadGetReqs; j++) {
                if (fIsSet(ncp->get_lead_list[j].flag, NC_REQ_TO_FREE))
                    continue; /* this request has been processed */
                if (ncp->get_lead_list[j].id == req_ids[i]) {
                    /* make this request to be freed */
                    fSet(ncp->get_lead_list[j].flag, NC_REQ_TO_FREE);
                    if (statuses != NULL) {
                        statuses[i] = NC_NOERR;
                        ncp->get_lead_list[j].status = statuses+i;
                    }
                    (*num_r_lead_reqs)++;
                    *num_r_reqs += ncp->get_lead_list[j].nonlead_num;
                    found = 1;
                    break;
                }
            }
        }
        if (found == 0) {
            /* no such request ID, if the program reached here */
            if (statuses != NULL)
                DEBUG_ASSIGN_ERROR(statuses[i], NC_EINVAL_REQUEST)
            /* retain the first error status */
            if (status == NC_NOERR)
                DEBUG_ASSIGN_ERROR(status, NC_EINVAL_REQUEST)
        }
    }
    if (status != NC_NOERR) return status;

    /* extract the requests from the pending queues (get_list and put_list
     * containing all nonblocking requests posted by far) into two separate
     * lists, get_list and put_list. Afterward, coalesce the pending lead and
     * non-lead request lists.
     */

    /* allocate put_list and get_list */
    if (*num_w_reqs)
        *put_list = (NC_req*) NCI_Malloc(sizeof(NC_req) * (*num_w_reqs));
    if (*num_r_reqs)
        *get_list = (NC_req*) NCI_Malloc(sizeof(NC_req) * (*num_r_reqs));

    /* copy over ncp->put_list and ncp->get_list to *put_list and *get_list */
    put_list_ptr = *put_list;
    get_list_ptr = *get_list;
    for (i=0; i<num_reqs; i++) {
        if (req_ids[i] == NC_REQ_NULL) continue; /* skip NULL request */

        if (req_ids[i] % 2 == 0) { /* write requests are even numbers */
            for (j=0; j<ncp->numLeadPutReqs; j++) {
                if (fIsSet(ncp->put_lead_list[j].flag, NC_REQ_TO_FREE) &&
                    req_ids[i] == ncp->put_lead_list[j].id) {
                    memcpy(put_list_ptr,
                           ncp->put_list + ncp->put_lead_list[j].nonlead_off,
                           sizeof(NC_req) * ncp->put_lead_list[j].nonlead_num);
                    put_list_ptr += ncp->put_lead_list[j].nonlead_num;
                    req_ids[i] = NC_REQ_NULL;
                    break;
                }
            }
        }
        else { /* read requests are odd numbers */
            for (j=0; j<ncp->numLeadGetReqs; j++) {
                if (fIsSet(ncp->get_lead_list[j].flag, NC_REQ_TO_FREE) &&
                    req_ids[i] == ncp->get_lead_list[j].id) {
                    memcpy(get_list_ptr,
                           ncp->get_list + ncp->get_lead_list[j].nonlead_off,
                           sizeof(NC_req) * ncp->get_lead_list[j].nonlead_num);
                    get_list_ptr += ncp->get_lead_list[j].nonlead_num;
                    req_ids[i] = NC_REQ_NULL;
                    break;
                }
            }
        }
    }

    /* coalesce non-lead put request queue */
    if (*num_w_reqs) {
        int k = 0;
        for (i=0; i<ncp->numLeadPutReqs; i++) {
            int off;

            /* skip the lead requests to be freed */
            if (fIsSet(ncp->put_lead_list[i].flag, NC_REQ_TO_FREE)) continue;
            off = ncp->put_lead_list[i].nonlead_off;
            if (off == k) { /* no need to coalesce this one */
                k += ncp->put_lead_list[i].nonlead_num;
                continue;
            }

            /* update offset index to the non-lead queue */
            ncp->put_lead_list[i].nonlead_off = k;

            /* src and dest may overlap */
            for (j=0; j<ncp->put_lead_list[i].nonlead_num; j++)
                ncp->put_list[k++] = ncp->put_list[off++];
        }
        ncp->numPutReqs -= *num_w_reqs;

        /* realloc ncp->put_list */
        if (ncp->numPutReqs == 0) {
            NCI_Free(ncp->put_list);
            ncp->put_list = NULL;
        }
        else {
            /* calculate the ceiling based on NC_REQUEST_CHUNK */
            int rem = NC_REQUEST_CHUNK - (ncp->numPutReqs % NC_REQUEST_CHUNK);
            size_t nChunks = ncp->numPutReqs + rem;
            ncp->put_list = (NC_req*) NCI_Realloc(ncp->put_list,
                                      nChunks * sizeof(NC_req));
        }
    }

    /* coalesce non-lead get request queue */
    if (*num_r_reqs) {
        int k = 0;
        for (i=0; i<ncp->numLeadGetReqs; i++) {
            int off;

            /* skip the lead requests to be freed */
            if (fIsSet(ncp->get_lead_list[i].flag, NC_REQ_TO_FREE)) continue;
            off = ncp->get_lead_list[i].nonlead_off;
            if (off == k) { /* no need to coalesce this one */
                k += ncp->get_lead_list[i].nonlead_num;
                continue;
            }

            /* update offset index to the non-lead queue */
            ncp->get_lead_list[i].nonlead_off = k;

            /* src and dest may overlap */
            for (j=0; j<ncp->get_lead_list[i].nonlead_num; j++)
                ncp->get_list[k++] = ncp->get_list[off++];
        }
        ncp->numGetReqs -= *num_r_reqs;

        /* realloc ncp->get_list */
        if (ncp->numGetReqs == 0) {
            NCI_Free(ncp->get_list);
            ncp->get_list = NULL;
        }
        else {
            /* calculate the ceiling based on NC_REQUEST_CHUNK */
            int rem = NC_REQUEST_CHUNK - (ncp->numGetReqs % NC_REQUEST_CHUNK);
            size_t nChunks = ncp->numGetReqs + rem;
            ncp->get_list = (NC_req*) NCI_Realloc(ncp->get_list,
                                      nChunks * sizeof(NC_req));
        }
    }

    /* coalescing lead queues has to wait after MPI-IO */

    return NC_NOERR;
}

/*----< req_commit() >-------------------------------------------------------*/
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
static int
req_commit(NC  *ncp,
           int  num_reqs,   /* number of requests */
           int *req_ids,    /* [num_reqs] */
           int *statuses,   /* [num_reqs] */
           int  coll_indep) /* NC_REQ_COLL or NC_REQ_INDEP */
{
    int i, j, err, status=NC_NOERR, do_read, do_write;
    int num_w_reqs, num_r_reqs, num_r_lead_reqs, num_w_lead_reqs;
    MPI_Offset newnumrecs=0;
    NC_req *put_list=NULL, *get_list=NULL;

    /* extract requests from put and get queues into put_list and get_list */
    err = extract_reqs(ncp, num_reqs, req_ids, statuses,
                       &num_r_lead_reqs, &num_r_reqs, &get_list,
                       &num_w_lead_reqs, &num_w_reqs, &put_list);

    /* calculate new number of records:
     * Need to update the number of records if new records have been created.
     * For nonblocking APIs, there is no way for a process to know whether
     * other processes write to a record variable or not. Hence, we must sync
     * the number of records for write requests.
     * Because netCDF classic files allow only one unlimited dimension, we can
     * scan all requests to find the maximum number of records.
     */
    newnumrecs = ncp->numrecs;
    for (i=0; i<num_w_lead_reqs; i++) {
        if (!IS_RECVAR(ncp->put_lead_list[i].varp) ||
            /* skip fixed-size variables */
            !fIsSet(ncp->put_lead_list[i].flag, NC_REQ_TO_FREE) ||
            /* skip requests not marked to be freed */
             fIsSet(ncp->put_lead_list[i].flag, NC_REQ_SKIP))
            /* skip invalid request */
            continue;
        newnumrecs = MAX(newnumrecs, ncp->put_lead_list[i].max_rec);
    }

    /* synchronize request metadata across processes if collective I/O */
    if (coll_indep == NC_REQ_COLL && ncp->nprocs > 1) {
        int mpireturn;
        MPI_Offset do_io[4];  /* [0]: read [1]: write [2]: error */
        do_io[0] = num_r_reqs;
        do_io[1] = num_w_reqs;
        do_io[2] = -err;   /* all NC errors are negative */
        do_io[3] = newnumrecs;
        TRACE_COMM(MPI_Allreduce)(MPI_IN_PLACE, do_io, 4, MPI_OFFSET,
                                  MPI_MAX, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");

        /* if error occurs, return the API collectively */
        if (do_io[2] != -NC_NOERR) return err;

        /* if at least one process has a non-zero request, all processes must
         * participate the collective read/write */
        do_read    = (do_io[0] > 0);
        do_write   = (do_io[1] > 0);
        newnumrecs = do_io[3];
    }
    else {
        if (err != NC_NOERR) return err;
        do_read  = (num_r_reqs > 0);
        do_write = (num_w_reqs > 0);
    }

    /* carry out writes and reads separately (writes first) */
    if (do_write > 0) {
        err = ncmpio_ina_nreqs(ncp, NC_REQ_WR, num_w_reqs, put_list,
                               newnumrecs);
        put_list = NULL; /* has been freed in the above call */

        /* Update the number of records if new records have been created.
         * For nonblocking APIs, there is no way for a process to know whether
         * others write to a record variable or not. Note newnumrecs has been
         * sync-ed and always >= ncp->numrecs.
         */
        if (coll_indep == NC_REQ_COLL) {
            if (newnumrecs > ncp->numrecs) {
                /* update new record number in file. Note newnumrecs is already
                 * sync-ed among all processes and in collective mode
                 * ncp->numrecs is always sync-ed in memory among processes,
                 * thus no need another MPI_Allreduce to sync it. */
                err = ncmpio_write_numrecs(ncp, newnumrecs);
                if (status == NC_NOERR) status = err;
                /* retain the first error if there is any */
                if (ncp->numrecs < newnumrecs) ncp->numrecs = newnumrecs;
            }
        }
        else { /* NC_REQ_INDEP */
            if (ncp->numrecs < newnumrecs) {
                ncp->numrecs = newnumrecs;
                set_NC_ndirty(ncp);
                /* delay numrecs sync until end_indep, redef or close */
            }
        }
    }
    if (do_read > 0) {
        err = ncmpio_ina_nreqs(ncp, NC_REQ_RD, num_r_reqs, get_list,
                               newnumrecs);
        get_list = NULL; /* has been freed in the above call */
    }

    /* retain the first error status */
    if (status == NC_NOERR) status = err;

    /* post-IO data processing: In write case, we may need to byte-swap user
     * write buf if it is used as the write buffer in MPI write call and the
     * target machine is little Endian. For read case, we may need to
     * unpack/byte-swap/type-convert a temp buffer to the user read buf
     */

    if (num_w_lead_reqs > 0) {
        j = 0;
        for (i=0; i<ncp->numLeadPutReqs; i++) {
            NC_lead_req *lead_req=ncp->put_lead_list+i;

            if (fIsSet(lead_req->flag, NC_REQ_TO_FREE)) {
                /* byte-swap the user buffer back to its original Endianness
                 * if it has been byte-swapped.
                 */
                if (fIsSet(lead_req->flag, NC_REQ_BUF_BYTE_SWAP))
                    ncmpii_in_swapn(lead_req->buf, lead_req->nelems,
                                    lead_req->varp->xsz);

                /* free resource allocated at lead request */
                if (lead_req->abuf_index < 0) {
                    if (fIsSet(lead_req->flag, NC_REQ_XBUF_TO_BE_FREED))
                        NCI_Free(lead_req->xbuf); /* free xbuf */
                }
                else if (ncp->abuf != NULL)  /* from bput API */
                    ncp->abuf->occupy_table[lead_req->abuf_index].is_used = 0;
            }
            else {
                if (j < i) {
                    /* coalesce put_lead_list[] and update put_list[].lead */
                    int k, off = ncp->put_lead_list[i].nonlead_off;
                    ncp->put_lead_list[j] = ncp->put_lead_list[i];
                    for (k=0; k<ncp->put_lead_list[i].nonlead_num; k++)
                        ncp->put_list[off++].lead_off = j;
                }
                j++;
            }
        }
        ncp->numLeadPutReqs = j;
        if (ncp->numLeadPutReqs == 0) {
            NCI_Free(ncp->put_list);
            NCI_Free(ncp->put_lead_list);
            ncp->put_list = NULL;
            ncp->put_lead_list = NULL;
        }

        /* once the bput requests are served, we reclaim the space and try
         * coalesce the freed space for the attached buffer
         */
        if (ncp->abuf != NULL) abuf_coalesce(ncp);
    }

    if (num_r_lead_reqs > 0) {
        j = 0;
        for (i=0; i<ncp->numLeadGetReqs; i++) {
            NC_lead_req *lead_req=ncp->get_lead_list+i;

            if (fIsSet(lead_req->flag, NC_REQ_TO_FREE)) {
                /* now, xbuf contains the data read from the file. It may need
                 * to be type-converted, byte-swapped, imap-unpacked, and
                 * buftype- unpacked from xbuf to buf. This is done in
                 * ncmpio_unpack_xbuf().
                 */
                int isContig, need_convert, need_swap;
                isContig = fIsSet(lead_req->flag, NC_REQ_BUF_TYPE_IS_CONTIG);
                need_convert = fIsSet(lead_req->flag, NC_REQ_BUF_TYPE_CONVERT);
                need_swap = fIsSet(lead_req->flag, NC_REQ_BUF_BYTE_SWAP);

                err = ncmpio_unpack_xbuf(ncp->format, lead_req->varp,
                                         lead_req->bufcount,
                                         lead_req->buftype,
                                         isContig,
                                         lead_req->nelems,
                                         lead_req->itype,
                                         lead_req->imaptype,
                                         need_convert,
                                         need_swap,
                                         lead_req->buf,
                                         lead_req->xbuf);
                if (err != NC_NOERR) {
                    if (lead_req->status != NULL &&
                        *lead_req->status == NC_NOERR)
                        *lead_req->status = err;
                    if (status == NC_NOERR) status = err;
                }

                if (fIsSet(lead_req->flag, NC_REQ_XBUF_TO_BE_FREED))
                    NCI_Free(lead_req->xbuf); /* free xbuf */

                if (!isContig && lead_req->buftype != MPI_DATATYPE_NULL)
                    MPI_Type_free(&lead_req->buftype);
            }
            else {
                if (j < i) {
                    /* coalesce get_lead_list[] and update get_list[].lead */
                    int k, off = ncp->get_lead_list[i].nonlead_off;
                    ncp->get_lead_list[j] = ncp->get_lead_list[i];
                    for (k=0; k<ncp->get_lead_list[i].nonlead_num; k++)
                        ncp->get_list[off++].lead_off = j;
                }
                j++;
            }
        }
        ncp->numLeadGetReqs = j;
        if (ncp->numLeadGetReqs == 0) {
            NCI_Free(ncp->get_list);
            NCI_Free(ncp->get_lead_list);
            ncp->get_list = NULL;
            ncp->get_lead_list = NULL;
        }
    }

    return status;
}

/*----< ncmpio_wait() >-------------------------------------------------------*/
int
ncmpio_wait(void *ncdp,
            int   num_reqs,
            int  *req_ids,   /* [num_reqs]: IN/OUT */
            int  *statuses,  /* [num_reqs] */
            int   reqMode)   /* only check if NC_REQ_COLL or NC_REQ_INDEP */
{
    NC *ncp = (NC*)ncdp;
    int coll_indep;

    if (NC_indef(ncp)) /* wait must be called in data mode */
        DEBUG_RETURN_ERROR(NC_EINDEFINE)

    coll_indep = (fIsSet(reqMode, NC_REQ_INDEP)) ? NC_REQ_INDEP : NC_REQ_COLL;

    /* check collective or independent mode */
    if (coll_indep == NC_REQ_INDEP && !NC_indep(ncp))
        DEBUG_RETURN_ERROR(NC_ENOTINDEP)
    else if (coll_indep == NC_REQ_COLL && NC_indep(ncp))
        DEBUG_RETURN_ERROR(NC_EINDEP)

    if (coll_indep == NC_REQ_INDEP && num_reqs == 0) return NC_NOERR;

    return req_commit(ncp, num_reqs, req_ids, statuses, coll_indep);
}

