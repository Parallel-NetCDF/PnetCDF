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
#include <limits.h> /* INT_MAX */
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

/* Prototypes for functions used only in this file */
static int wait_getput(NC *ncp, int num_reqs, NC_req *reqs, int rw_flag,
                       int coll_indep, MPI_Offset newnumrecs);

static int mgetput(NC *ncp, int num_reqs, NC_req *reqs, int rw_flag,
                   int coll_indep);

/*----< ncmpio_getput_zero_req() >-------------------------------------------*/
/* This function is called when this process has zero-length I/O request and
 * must participate all the MPI collective calls involved in the collective
 * APIs and wait_all(), which include setting fileview, collective read/write,
 * another setting fileview.
 *
 * This function is collective.
 */
int
ncmpio_getput_zero_req(NC *ncp, int reqMode)
{
    int err, mpireturn, status=NC_NOERR;
    MPI_Status mpistatus;
    MPI_File fh;

    /* do nothing if this came from an independent API */
    if (fIsSet(reqMode, NC_REQ_INDEP)) return NC_NOERR;

    fh = ncp->collective_fh;

    TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                MPI_INFO_NULL);

    if (fIsSet(reqMode, NC_REQ_RD)) {
        TRACE_IO(MPI_File_read_all)(fh, NULL, 0, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_read_all");
            err = (err == NC_EFILE) ? NC_EREAD : err;
            DEBUG_ASSIGN_ERROR(status, err)
        }
    } else { /* write request */
        TRACE_IO(MPI_File_write_all)(fh, NULL, 0, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_write_all");
            err = (err == NC_EFILE) ? NC_EWRITE : err;
            DEBUG_ASSIGN_ERROR(status, err)
        }
    }

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                 MPI_INFO_NULL);
     */

    return status;
}

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
    int i, j, status=NC_NOERR;
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

/*----< construct_filetypes() >----------------------------------------------*/
/* concatenate the requests into a single MPI derived filetype */
static int
construct_filetypes(NC           *ncp,
                    NC_lead_req  *lead_list, /* NC_REQ_WR or NC_REQ_RD */
                    int           num_reqs,
                    int          *blocklens, /* [num_reqs] temp buffer */
                    MPI_Aint     *disps,     /* [num_reqs] temp buffer */
                    NC_req       *reqs,      /* [num_reqs] */
                    MPI_Datatype *filetype)  /* OUT */
{
    int i, j, err, status=NC_NOERR, all_ftype_contig=1, last_contig_req;
    MPI_Datatype *ftypes;

    if (num_reqs <= 0) { /* for participating collective call */
        *filetype = MPI_BYTE;
        return NC_NOERR;;
    }

    /* hereinafter, num_reqs > 0 */
    ftypes = (MPI_Datatype*) NCI_Malloc((size_t)num_reqs * sizeof(MPI_Datatype));

    /* create a filetype for each request */
    last_contig_req = -1; /* index of the last contiguous request */
    j = 0;                /* index of last valid ftypes */
    for (i=0; i<num_reqs; i++, j++) {
        int is_ftype_contig, ndims;
        NC_lead_req *lead;

        lead = lead_list + reqs[i].lead_off;
        ndims = lead->varp->ndims;

        ftypes[j] = MPI_BYTE; /* in case the call below failed */

        if (ndims == 0) { /* scalar variable */
#if SIZEOF_MPI_AINT < SIZEOF_MPI_OFFSET
            if (lead->varp->begin > INT_MAX) {
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
                fSet(lead->flag, NC_REQ_SKIP); /* skip this request */
                if ( lead->status != NULL &&
                    *lead->status == NC_NOERR)
                    *lead->status = err;
                status = err; /* report first error */
            }
#endif
            disps[j]        = lead->varp->begin;
            is_ftype_contig = 1;
        }
        else { /* non-scalar variable */
            MPI_Offset offset, *count, *stride;
            count  = reqs[i].start + ndims;
            stride = fIsSet(lead->flag, NC_REQ_STRIDE_NULL) ?
                     NULL : count + ndims;

            err = ncmpio_filetype_create_vars(ncp,
                                              lead->varp,
                                              reqs[i].start,
                                              count,
                                              stride,
                                              &offset,
                                              &ftypes[j],
                                              &is_ftype_contig);

#if SIZEOF_MPI_AINT < SIZEOF_MPI_OFFSET
            if (err == NC_NOERR && offset > INT_MAX)
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
#endif
            disps[j] = (MPI_Aint)offset;

            if (err != NC_NOERR) {
                fSet(lead->flag, NC_REQ_SKIP); /* skip this request */
                if ( lead->status != NULL &&
                    *lead->status == NC_NOERR)
                    *lead->status = err;
                if (status == NC_NOERR) status = err; /* report first error */
                continue;
            }
        }

        if (is_ftype_contig) {
            MPI_Offset coalesced_len;

            /* No need to construct a filetype */
            blocklens[j] = lead->varp->xsz * reqs[i].nelems;
            coalesced_len = blocklens[j];
            if (last_contig_req >= 0)
                coalesced_len += blocklens[last_contig_req];
            /* if coalesced_len overflows 4-byte int, then skip coalescing */
            if (coalesced_len < INT_MAX && last_contig_req >= 0 &&
                disps[j] - disps[last_contig_req] ==
                blocklens[last_contig_req]) {
                blocklens[last_contig_req] = coalesced_len;
                j--;
            }
            else last_contig_req = j;
        }
        else {
            /* we will construct a filetype, set blocklen to 1 */
            blocklens[j] = 1;
            last_contig_req = -1;
            all_ftype_contig = 0;
        }
    }
    /* j is the new num_reqs */
    num_reqs = j;

    if (status != NC_NOERR) {
        /* even if error occurs, we still must participate the collective
           call to MPI_File_set_view() */
        *filetype = MPI_BYTE;
    }
    else if (num_reqs == 1 && disps[0] == 0) {
        if (ftypes[0] == MPI_BYTE)
            *filetype = MPI_BYTE;
        else
            MPI_Type_dup(ftypes[0], filetype);
    }
    else { /* if (num_reqs > 1 || (num_reqs == 1 && disps[0] > 0)) */
        /* all ftypes[] created fine, now concatenate all ftypes[] */
        if (all_ftype_contig) {
            err = MPI_Type_create_hindexed(num_reqs, blocklens, disps,
                                           MPI_BYTE, filetype);
            MPI_Type_commit(filetype);
        }
        else {
            int mpireturn;
            mpireturn = MPI_Type_create_struct(num_reqs, blocklens, disps,
                                               ftypes, filetype);
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_struct");
            else {
                MPI_Type_commit(filetype);
                err = NC_NOERR;
            }
        }

        if (err != NC_NOERR) *filetype = MPI_BYTE;
        if (status == NC_NOERR) status = err; /* report the first error */
    }

    for (i=0; i<num_reqs; i++) {
        if (ftypes[i] != MPI_BYTE)
            MPI_Type_free(&ftypes[i]);
    }
    NCI_Free(ftypes);

    return status;
}

/*----< construct_buffertypes() >--------------------------------------------*/
/* the input requests, reqs[], are non-interleaving requests */
static int
construct_buffertypes(NC_lead_req  *lead_list,
                      int           num_reqs,
                      int          *blocklens, /* [num_reqs] temp buffer */
                      MPI_Aint     *disps,     /* [num_reqs] temp buffer */
                      NC_req       *reqs,      /* [num_reqs] */
                      MPI_Datatype *buf_type)  /* OUT */
{
    int i, j, k, status=NC_NOERR, mpireturn;
    MPI_Aint a0, ai;

    *buf_type = MPI_BYTE;
    if (num_reqs == 0) return NC_NOERR;

    /* create the I/O buffer derived data type */

    /* calculate blocklens[], and disps[] */
    for (i=0, j=0; i<num_reqs; i++) {
        MPI_Offset req_size;
        NC_lead_req *lead;

        lead = lead_list + reqs[i].lead_off;

        if (fIsSet(lead->flag, NC_REQ_SKIP)) continue;

        req_size = lead->varp->xsz;
        if (lead->varp->ndims > 0) { /* non-scalar variable */
            MPI_Offset *count = reqs[i].start + lead->varp->ndims;
            if (!IS_RECVAR(lead->varp)) req_size *= count[0];
            for (k=1; k<lead->varp->ndims; k++) req_size *= count[k];
        }

        /* check int overflow */
        if (req_size > INT_MAX) { /* skip this request */
            fSet(lead->flag, NC_REQ_SKIP);
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
            continue;
        }
        blocklens[j] = (int)req_size;

        MPI_Get_address(reqs[i].xbuf, &ai);
        if (j == 0) a0 = ai;
        disps[j] = ai - a0;
        j++;
    }
    /* update num_reqs to number of valid requests */
    num_reqs = j;

    if (num_reqs > 0) {
        /* concatenate buffer addresses into a single buffer type */
        mpireturn = MPI_Type_create_hindexed(num_reqs, blocklens, disps,
                                             MPI_BYTE, buf_type);
        if (mpireturn != MPI_SUCCESS) {
            int err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else
            MPI_Type_commit(buf_type);
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
        *put_list = (NC_req*) NCI_Malloc((*num_w_reqs)*sizeof(NC_req));
    if (*num_r_reqs)
        *get_list = (NC_req*) NCI_Malloc((*num_r_reqs)*sizeof(NC_req));

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
                           ncp->put_lead_list[j].nonlead_num * sizeof(NC_req));
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
                           ncp->get_lead_list[j].nonlead_num * sizeof(NC_req));
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
    if (coll_indep == NC_REQ_COLL) {
        int mpireturn;
        MPI_Offset io_req[4], do_io[4];  /* [0]: read [1]: write [2]: error */
        io_req[0] = num_r_reqs;
        io_req[1] = num_w_reqs;
        io_req[2] = -err;   /* all NC errors are negative */
        io_req[3] = newnumrecs;
        TRACE_COMM(MPI_Allreduce)(io_req, do_io, 4, MPI_OFFSET, MPI_MAX,
                                  ncp->comm);
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
        err = wait_getput(ncp, num_w_reqs, put_list, NC_REQ_WR, coll_indep,
                          newnumrecs);
        put_list = NULL; /* has been freed in wait_getput() */
    }

    if (do_read > 0) {
        err = wait_getput(ncp, num_r_reqs, get_list, NC_REQ_RD, coll_indep,
                          newnumrecs);
        get_list = NULL; /* has been freed in wait_getput() */
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

#ifdef ENABLE_REQ_AGGREGATION
    /* check collective or independent mode */
    if (coll_indep == NC_REQ_INDEP && !NC_indep(ncp))
        DEBUG_RETURN_ERROR(NC_ENOTINDEP)
    else if (coll_indep == NC_REQ_COLL && NC_indep(ncp))
        DEBUG_RETURN_ERROR(NC_EINDEP)

    if (coll_indep == NC_REQ_INDEP && num_reqs == 0) return NC_NOERR;

    return req_commit(ncp, num_reqs, req_ids, statuses, coll_indep);
#else
    /* If request aggregation is disabled, we call an independent wait() for
     * each request
     */
    int i, status=NC_NOERR, err;

    if (coll_indep == NC_REQ_INDEP) {
        /* This is called from ncmpi_wait(), which is an independent call
         * Argument num_reqs can be NC_REQ_ALL which means to flush all pending
         * nonblocking requests. In this case, arguments req_ids and statuses
         * will be ignored.
         * Argument num_reqs must either be NC_REQ_ALL, NC_GET_REQ_ALL,
         * NC_PUT_REQ_ALL, or a non-negative value.
         * Argument statuses can be NULL, meaning the caller only cares about
         * the error code returned by this call, but not the statuses of
         * individual nonblocking requests.
         */
        if (num_reqs == 0) return NC_NOERR;

        /* This is called from ncmpi_wait which must be called in independent
         * data mode, illegal in collective mode.
         */
        if (!NC_indep(ncp)) DEBUG_RETURN_ERROR(NC_ENOTINDEP);

        if (coll_indep == NC_REQ_INDEP && num_reqs == 0) return NC_NOERR;
    }
    else {
        /* This is called from ncmpi_wait_all(), which is a collective call
         * Argument num_reqs can be NC_REQ_ALL which means to flush all pending
         * nonblocking requests. In this case, arguments req_ids and statuses
         * will be ignored.
         * Argument num_reqs must either be NC_REQ_ALL, NC_GET_REQ_ALL,
         * NC_PUT_REQ_ALL, or a non-negative value.
         * Argument statuses can be NULL, meaning the caller only cares about
         * the error code returned by this call, but not the statuses of
         * individual nonblocking requests.
         */
        /* the following line CANNOT be added, because ncmpi_wait_all() is a
         * collective call, all processes must participate some MPI collective
         * operations used later on.
         */
        /* if (num_reqs == 0) return NC_NOERR; */

        /* This is called from ncmpi_wait_all which must be called in
         * collective data mode, illegal in independent mode. This also
         * ensures the program will returns back to collective mode.
         */
        if (NC_indep(ncp)) DEBUG_RETURN_ERROR(NC_EINDEP);

        /* must enter independent mode, as num_reqs may be different among
           processes */
        err = ncmpio_begin_indep_data(ncp);
        if (status == NC_NOERR) status = err;
    }

    if (num_reqs <= NC_REQ_ALL) { /* flush all get or put pending requests */
        if (num_reqs == NC_REQ_ALL || num_reqs == NC_GET_REQ_ALL) {
            while (ncp->numLeadGetReqs) {
                /* commit one request at a time. Note ncp->numLeadGetReqs
                 * will be descreased in req_commit()
                 */
                err = req_commit(ncp, 1, &ncp->get_lead_list[0].id, NULL,
                                 NC_REQ_INDEP);
                if (status == NC_NOERR) status = err;
            }
        }
        if (num_reqs == NC_REQ_ALL || num_reqs == NC_PUT_REQ_ALL) {
            while (ncp->numLeadPutReqs) {
                /* commit one request at a time. Note ncp->numLeadPutReqs
                 * will be descreased in req_commit()
                 */
                err = req_commit(ncp, 1, &ncp->put_lead_list[0].id, NULL,
                                 NC_REQ_INDEP);
                if (status == NC_NOERR) status = err;
            }
        }
    }
    else {
        for (i=0; i<num_reqs; i++) { /* commit one request at a time */
            err = req_commit(ncp, 1, &req_ids[i],
                  (statuses == NULL) ? NULL : &statuses[i], NC_REQ_INDEP);
            if (status == NC_NOERR) status = err;
        }
    }

    if (coll_indep == NC_REQ_COLL) {
        /* return to collective data mode */
        err = ncmpio_end_indep_data(ncp);
        if (status == NC_NOERR) status = err;
    }

    return status; /* return the first error encountered, if there is any */
#endif
}

/* C struct for breaking down a request to a list of offset-length segments */
typedef struct {
    MPI_Offset off;      /* starting file offset of the request */
    MPI_Offset len;      /* requested length in bytes starting from off */
    MPI_Aint   buf_addr; /* distance of this request's I/O buffer to the first
                            request to be merged */
} off_len;

/*----< off_compare() >------------------------------------------------------*/
/* used for sorting the offsets of the off_len array */
static int
off_compare(const void *a, const void *b)
{
    if (((off_len*)a)->off > ((off_len*)b)->off) return  1;
    if (((off_len*)a)->off < ((off_len*)b)->off) return -1;
    return 0;
}

/*----< vars_flatten() >-----------------------------------------------------*/
/* flatten a subarray request into a list of offset-length pairs */
static MPI_Offset
vars_flatten(int          ndim,    /* number of dimensions */
             int          el_size, /* array element size */
             MPI_Offset  *dimlen,  /* [ndim] dimension lengths */
             MPI_Offset   offset,  /* starting file offset of variable */
             MPI_Aint     buf_addr,/* starting buffer address */
             MPI_Offset  *start,   /* [ndim] starts of subarray */
             MPI_Offset  *count,   /* [ndim] counts of subarray */
             MPI_Offset  *stride,  /* [ndim] strides of subarray */
             MPI_Offset  *nseg,    /* OUT: number of segments */
             off_len     *seg)     /* OUT: array of segments */
{
    int i, j, to_free_stride=0;
    MPI_Offset seg_len, nstride, array_len, off, subarray_len;
    off_len *ptr=seg, *seg0;

    *nseg = 0;
    if (ndim < 0) return *nseg;

    if (ndim == 0) {  /* scalar record variable */
        *nseg = 1;
        seg->off      = offset;
        seg->len      = el_size;
        seg->buf_addr = buf_addr;
        return *nseg;
    }

    if (stride == NULL) { /* equivalent to {1, 1, ..., 1} */
        stride = (MPI_Offset*) NCI_Malloc((size_t)ndim * SIZEOF_MPI_OFFSET);
        for (i=0; i<ndim; i++) stride[i] = 1;
        to_free_stride = 1;
    }

    /* TODO: check if all stride[] >= 1
       Q: Is it legal if any stride[] <= 0 ? */

    /* calculate the number of offset-length pairs */
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

/*----< merge_requests() >---------------------------------------------------*/
static int
merge_requests(NC          *ncp,
               NC_lead_req *lead_list,
               int          num_reqs,
               NC_req      *reqs,    /* [num_reqs] */
               void       **buf,     /* OUT: 1st I/O buf addr */
               MPI_Offset  *nsegs,   /* OUT: no. off-len pairs */
               off_len    **segs)    /* OUT: [*nsegs] */
{
    int i, j, status=NC_NOERR, ndims;
    MPI_Offset  nseg, *start, *count, *shape, *stride;
    MPI_Aint addr, buf_addr;

    *nsegs = 0;    /* total number of offset-length pairs */
    *segs  = NULL; /* array of offset-length pairs */

    /* note invalid requests have been removed in wait_getput() */
    *buf = reqs[0].xbuf; /* I/O buffer of first request */

    /* buf_addr is the buffer address of the first request */
    MPI_Get_address(reqs[0].xbuf, &buf_addr);

    /* Count the number off-len pairs from reqs[], so we can malloc a
     * contiguous memory space for storing off-len pairs
     */
    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead = lead_list + reqs[i].lead_off;
        ndims  = lead->varp->ndims;
        start  = reqs[i].start;
        count  = start + ndims;
        stride = count + ndims;

        /* for record variable, each reqs[] is within a record */
        if (IS_RECVAR(lead->varp)) {
            ndims--;
            start++;
            count++;
            stride++;
        }
        if (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) stride = NULL;

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
    off_len *seg_ptr = (off_len*)NCI_Malloc((size_t)(*nsegs) * sizeof(off_len));
    *segs = seg_ptr;

    /* now re-run the loop to fill in the off-len pairs */
    for (i=0; i<num_reqs; i++) {
        MPI_Offset var_begin;
        NC_lead_req *lead = lead_list + reqs[i].lead_off;

        /* buf_addr is the buffer address of the first valid request */
        MPI_Get_address(reqs[i].xbuf, &addr);
        addr -= buf_addr,  /* distance to the buf of first req */

        ndims  = lead->varp->ndims;
        start  = reqs[i].start;
        count  = start + ndims;
        stride = count + ndims;
        shape  = lead->varp->shape;

        /* find the starting file offset for this variable */
        var_begin = lead->varp->begin;

        /* for record variable, each reqs[] is within a record */
        if (IS_RECVAR(lead->varp)) {
            ndims--;
            start++;
            count++;
            stride++;
            shape++;
            /* find the starting file offset for this record */
            var_begin += reqs[i].start[0] * ncp->recsize;
        }

        if (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) stride = NULL;

        /* flatten each request to a list of offset-length pairs */
        vars_flatten(ndims, lead->varp->xsz, shape, var_begin,
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
        qsort(*segs, (size_t)(*nsegs), sizeof(off_len), off_compare);

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

/*----< type_create_off_len() >----------------------------------------------*/
static int
type_create_off_len(MPI_Offset    nsegs,    /* no. off-len pairs */
                    off_len      *segs,     /* [nsegs] off-len pairs (sorted) */
                    MPI_Datatype *filetype, /* OUT */
                    MPI_Datatype *buf_type) /* OUT */
{
    int i, j, *blocklens, mpireturn;
    MPI_Offset next_off, next_len, true_nsegs;
    MPI_Aint *disps;

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
    true_nsegs = j + 1;
    blocklens = (int*)      NCI_Malloc(true_nsegs * SIZEOF_INT);
    disps     = (MPI_Aint*) NCI_Malloc(true_nsegs * SIZEOF_MPI_AINT);

    /* coalesce segs[].off and len to disps[] and blocklens[] */
    if (segs[0].len > INT_MAX) {
        NCI_Free(disps);
        NCI_Free(blocklens);
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    }
    disps[0]     =      segs[0].off;
    blocklens[0] = (int)segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (segs[i].len > INT_MAX) {
            NCI_Free(disps);
            NCI_Free(blocklens);
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
        if (disps[j] + blocklens[j] == segs[i].off)
            /* j and i are contiguous */
            blocklens[j] += (int)segs[i].len;
            /* TODO: take care of 4-byte int overflow problem */
        else {
            j++;
            disps[j]     =      segs[i].off;
            blocklens[j] = (int)segs[i].len;
        }
    }
    /* j+1 is the coalesced length */

    mpireturn = MPI_Type_create_hindexed(j+1, blocklens, disps, MPI_BYTE,
                                         filetype);
    if (mpireturn != MPI_SUCCESS) {
        *filetype = MPI_BYTE;
        *buf_type = MPI_BYTE;
        NCI_Free(disps);
        NCI_Free(blocklens);
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");
    }
    MPI_Type_commit(filetype);

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
    if (true_nsegs < j + 1) {
        blocklens = (int*)      NCI_Realloc(blocklens, (j+1) * SIZEOF_INT);
        disps     = (MPI_Aint*) NCI_Realloc(disps,     (j+1) * SIZEOF_MPI_AINT);
    }

    /* coalesce segs[].off and len to disps[] and blocklens[] */
    if (segs[0].len > INT_MAX) {
        NCI_Free(disps);
        NCI_Free(blocklens);
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    }
    disps[0]     =      segs[0].buf_addr;
    blocklens[0] = (int)segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (segs[i].len > INT_MAX) {
            NCI_Free(disps);
            NCI_Free(blocklens);
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
        if (disps[j] + blocklens[j] == segs[i].buf_addr)
            /* j and i are contiguous */
            blocklens[j] += (int)segs[i].len;
        else {
            j++;
            disps[j]     =      segs[i].buf_addr;
            blocklens[j] = (int)segs[i].len;
        }
    }
    /* j+1 is the coalesced length */
    mpireturn = MPI_Type_create_hindexed(j+1, blocklens, disps, MPI_BYTE,
                                         buf_type);
    NCI_Free(disps);
    NCI_Free(blocklens);
    if (mpireturn != MPI_SUCCESS) {
        if (*filetype != MPI_BYTE) MPI_Type_free(filetype);
        *filetype = MPI_BYTE;
        *buf_type = MPI_BYTE;
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");
    }
    MPI_Type_commit(buf_type);

    return NC_NOERR;
}

/*----< req_compare() >------------------------------------------------------*/
/* used to sort the the string file offsets of reqs[] */
static int
req_compare(const void *a, const void *b)
{
    if (((NC_req*)a)->offset_start > ((NC_req*)b)->offset_start) return (1);
    if (((NC_req*)a)->offset_start < ((NC_req*)b)->offset_start) return (-1);
    return (0);
}

/*----< req_aggregation() >--------------------------------------------------*/
/* aggregate multiple read/write (non-contiguous) requests and call MPI-IO
 */
static int
req_aggregation(NC     *ncp,
                int     num_reqs,    /* IN/OUT: # requests */
                NC_req *reqs,        /* [num_reqs] sorted requests, to be freed
                                        in this subroutine */
                int     rw_flag,     /* NC_REQ_WR or NC_REQ_RD */
                int     coll_indep,  /* NC_REQ_COLL or NC_REQ_INDEP */
                int     interleaved) /* interleaved in reqs[] */
{
    int i, gtype, err, status=NC_NOERR, ngroups, mpireturn, buf_len;
    int *group_index, *group_type, *f_blocklens, *b_blocklens;
    int  numLeadReqs;
    NC_lead_req *lead_list;
    void *buf; /* point to starting buffer, used by MPI-IO call */
    MPI_Aint      b_begin, b_addr, *f_disps, *b_disps;
    MPI_Datatype  filetype, buf_type, *ftypes, *btypes;
    MPI_File fh;
    MPI_Offset max_end, offset;

    if (num_reqs == 0) { /* only NC_REQ_COLL can reach here for 0 request */
        assert(coll_indep == NC_REQ_COLL);
        /* simply participate the collective call */
        return ncmpio_getput_zero_req(ncp, rw_flag);
    }
    if (! interleaved) {
        /* concatenate all filetypes into a single one and do I/O */
        return mgetput(ncp, num_reqs, reqs, rw_flag, coll_indep);
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
     * overlap), then we cannot simply concatenate the filetypes of individual
     * requests. This approach flattens the requests of "interleaved" groups
     * into offset-length pairs, sorts, and merges them into an aggregated
     * filetype. Similar for building an aggregated I/O buffer type.
     */

    /* first calculate the number of groups, so group_index and group_type can
       be malloc-ed. Group type: 0 for non-interleaved group and 1 for
       interleaved group.
     */
#define INTERLEAVED    1
#define NONINTERLEAVED 0

    assert(num_reqs > 1);
    ngroups = 1;
    gtype   = (reqs[0].offset_end > reqs[1].offset_start) ?
              INTERLEAVED : NONINTERLEAVED;
    max_end = MAX(reqs[0].offset_end, reqs[1].offset_end);
    for (i=1; i<num_reqs-1; i++) {
        if (gtype == NONINTERLEAVED &&
            reqs[i].offset_end > reqs[i+1].offset_start) {
            /* Done with this NONINTERLEAVED group. Continue to construct
             * next group, starting from reqs[i], which will be INTERLEAVED. */
            ngroups++;
            gtype = INTERLEAVED;
            max_end = MAX(reqs[i].offset_end, reqs[i+1].offset_end);
        }
        else if (gtype == INTERLEAVED) {
            if (max_end <= reqs[i+1].offset_start) {
                /* Done with this INTERLEAVED group. Continue to construct
                 * next group. First check whether the next group is
                 * INTERLEAVED or NONINTERLEAVED. */
                gtype = NONINTERLEAVED;
                if (i+2 < num_reqs &&
                    reqs[i+1].offset_end > reqs[i+2].offset_start)
                    gtype = INTERLEAVED; /* next group is also interleaved */
                ngroups++;
                max_end = reqs[i+1].offset_end;
            }
            else
                max_end = MAX(max_end, reqs[i+1].offset_end);
        }
    }

    group_index = (int*) NCI_Malloc((size_t)(ngroups+1) * 2 * SIZEOF_INT);
    group_type  = group_index + (ngroups+1);

    /* calculate the starting index of each group and determine group type */
    ngroups        = 1;
    gtype          = (reqs[0].offset_end > reqs[1].offset_start) ?
                     INTERLEAVED : NONINTERLEAVED;
    max_end        = MAX(reqs[0].offset_end, reqs[1].offset_end);
    group_index[0] = 0;
    group_type[0]  = gtype;
    for (i=1; i<num_reqs-1; i++) {
        if (gtype == NONINTERLEAVED &&
            reqs[i].offset_end > reqs[i+1].offset_start) {
            /* Done with this NONINTERLEAVED group. Continue to construct
             * next group, which will be INTERLEAVED. */
            /* reqs[i] starts a new interleaved group */
            group_index[ngroups] = i;
            gtype = INTERLEAVED;
            group_type[ngroups] = gtype;
            ngroups++;
            max_end = MAX(reqs[i].offset_end, reqs[i+1].offset_end);
        }
        else if (gtype == INTERLEAVED) {
            if (max_end <= reqs[i+1].offset_start) {
                /* Done with this INTERLEAVED group. Continue to construct
                 * next group. First check whether the next group is
                 * INTERLEAVED or NONINTERLEAVED. */
                gtype = NONINTERLEAVED;
                if (i+2 < num_reqs &&
                    reqs[i+1].offset_end > reqs[i+2].offset_start)
                    gtype = INTERLEAVED; /* next group is also interleaved */
                /* the interleaved group ends with reqs[i] */
                group_index[ngroups] = i+1;
                group_type[ngroups] = gtype;
                ngroups++;
                max_end = reqs[i+1].offset_end;
            }
            else
                max_end = MAX(max_end, reqs[i+1].offset_end);
        }
    }
    group_index[ngroups] = num_reqs; /* to indicate end of groups */

    /* for each group, construct one filetype by concatenating if the group
     * is non-interleaved and by flatten/sort/merge if the group is
     * interleaved. At the end, all ngroups filetypes are concatenated into
     * a single filetype. Similar for constructing buffer types.
     * Then use one collective I/O to commit.
     */

    ftypes = (MPI_Datatype*) NCI_Malloc((size_t)ngroups*2*sizeof(MPI_Datatype));
    btypes = ftypes + ngroups;
    f_blocklens = (int*) NCI_Malloc((size_t)ngroups*2*SIZEOF_INT);
    b_blocklens = f_blocklens + ngroups;
    f_disps = (MPI_Aint*) NCI_Malloc((size_t)ngroups*2*SIZEOF_MPI_AINT);
    b_disps = f_disps + ngroups;

    buf = reqs[0].xbuf; /* the buffer of 1st request */
    b_disps[0] = 0;     /* relative to address of 1st buf */
    MPI_Get_address(buf, &b_begin);

    /* temp buffers, used by multiple calls to construct_filetypes()  */
    int *blocklens = (int*) NCI_Malloc((size_t)num_reqs*SIZEOF_INT);
    MPI_Aint *disps = (MPI_Aint*) NCI_Malloc((size_t)num_reqs*SIZEOF_MPI_AINT);

    lead_list = (rw_flag == NC_REQ_RD) ? ncp->get_lead_list
                                       : ncp->put_lead_list;
    /* for each group, build a filetype and a buffer type in ftypes[i] and
       btypes[i] */
    for (i=0; i<ngroups; i++) {
        NC_req *g_reqs = reqs + group_index[i];
        int     g_num_reqs = group_index[i+1] - group_index[i];
        f_disps[i] = 0;  /* file displacements always to the file offset 0 */

        if (group_type[i] == NONINTERLEAVED) {
            /* This group contains no interleaved filetypes, so we can
             * simply concatenate filetypes of this group into a single one
             */
            err = construct_filetypes(ncp, lead_list, g_num_reqs, blocklens,
                                      disps, g_reqs, &ftypes[i]);
            if (status == NC_NOERR) status = err;
            if (err != NC_NOERR) { /* skip this group */
                ftypes[i] = btypes[i] = MPI_BYTE;
                f_blocklens[i] = 0;
                continue;
            }
            f_blocklens[i] = 1;

            /* concatenate buffer types of this group into a single one */
            err = construct_buffertypes(lead_list, g_num_reqs, blocklens, disps,
                                        g_reqs, &btypes[i]);
            if (status == NC_NOERR) status = err;
            if (err != NC_NOERR) { /* skip this group */
                ftypes[i] = btypes[i] = MPI_BYTE;
                b_blocklens[i] = 0;
                f_blocklens[i] = 0;
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
             * The pitfall of this flattening is the additional memory
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

            /* merge all requests into sorted offset-length pairs. Note
             * g_reqs[].offset_start and offset_end are relative to the
             * beginning of file */
            err = merge_requests(ncp, lead_list, g_num_reqs, g_reqs,
                                 &merged_buf, &nsegs, &segs);
            if (status == NC_NOERR) status = err;
            if (err != NC_NOERR) { /* skip this group */
                ftypes[i] = btypes[i] = MPI_BYTE;
                b_blocklens[i] = 0;
                f_blocklens[i] = 0;
                if (segs != NULL) NCI_Free(segs);
                continue;
            }
            assert(nsegs > 0);

            /* sges[] will be used to construct fileview and buffer type */
            err = type_create_off_len(nsegs, segs, &ftypes[i], &btypes[i]);
            /* preserve the previous error if there is any */
            if (status == NC_NOERR) status = err;
            NCI_Free(segs);
            if (err != NC_NOERR) { /* skip this group */
                ftypes[i] = btypes[i] = MPI_BYTE;
                b_blocklens[i] = 0;
                f_blocklens[i] = 0;
                continue;
            }
            f_blocklens[i] = 1;
        }

        if (i > 0) {
            /* get the buffer address of the first request in this group */
            MPI_Get_address(g_reqs[0].xbuf, &b_addr);
            b_disps[i] = b_addr - b_begin; /* to 1st buffer of 1st group*/
        }
        b_blocklens[i] = 1;
    }
    NCI_Free(disps);
    NCI_Free(blocklens);
    NCI_Free(group_index);

    buf_len=1;

    if (ngroups == 1) {
        /* use ftypes[0] and btypes[0] directly */
        filetype = ftypes[0];
        buf_type = btypes[0];
    }
    else {
        /* concatenate all ftypes[] to filetype */
        mpireturn = MPI_Type_create_struct(ngroups, f_blocklens, f_disps,
                                           ftypes, &filetype);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_struct");
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
        mpireturn = MPI_Type_create_struct(ngroups, b_blocklens, b_disps,
                                           btypes, &buf_type);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_struct");
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
    NCI_Free(ftypes);
    NCI_Free(f_blocklens);
    NCI_Free(f_disps);

    /* non-lead request list is no longer used once fileview and buftype have
     * been constructed. Free the start arrays allocated at lead requests.
     */
    numLeadReqs = (rw_flag == NC_REQ_RD) ? ncp->numLeadGetReqs
                                         : ncp->numLeadPutReqs;
    for (i=0; i<numLeadReqs; i++) {
        if (!fIsSet(lead_list[i].flag, NC_REQ_TO_FREE)) continue;
        NCI_Free(lead_list[i].start);
        lead_list[i].start = NULL;
    }
    NCI_Free(reqs);

    if (coll_indep == NC_REQ_COLL)
        fh = ncp->collective_fh;
    else
        fh = ncp->independent_fh;

    /* set the MPI-IO fileview, this is a collective call */
    offset = 0;
    err = ncmpio_file_set_view(ncp, fh, &offset, filetype);
    if (filetype != MPI_BYTE) MPI_Type_free(&filetype);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        if (coll_indep == NC_REQ_INDEP) return status;
        buf_len = 0;
    }

    /* call MPI_File_read/MPI_File_write */
    err = ncmpio_read_write(ncp, rw_flag, coll_indep, offset, buf_len, buf_type,
                            buf, ((buf_type == MPI_BYTE) ? 1 : 0));
    if (status == NC_NOERR) status = err;

    if (buf_type != MPI_BYTE) MPI_Type_free(&buf_type);

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                 MPI_INFO_NULL);
     */

    return status;
}

/*----< calculate_access_range() >-------------------------------------------*/
/* Returns the file offsets of access range of this request: starting file
 * offset and end offset (exclusive).
 * Note zero-length request should never call this subroutine.
 */
static int
calculate_access_range(const NC         *ncp,
                       const NC_var     *varp,
                       const MPI_Offset *start,     /* [varp->ndims] */
                       const MPI_Offset *count,     /* [varp->ndims] */
                       const MPI_Offset *stride,    /* [varp->ndims] */
                       MPI_Offset       *start_off, /* OUT: start offset */
                       MPI_Offset       *end_off)   /* OUT: end   offset */
{
    int i, ndims = varp->ndims; /* number of dimensions of this variable */

    /*
     * varp->dsizes[] is computed from right to left product of shape
     * For example, a 3D array of size 5x4x3 in C order,
     * For fixed-size variable: dsizes[0]=60 dsizes[1]=12 dsizes[2]=3
     * For record     variable: dsizes[0]=12 dsizes[1]=12 dsizes[2]=3
     */
    if (IS_RECVAR(varp)) {
        *start_off = 0;
        *end_off   = 0;
        if (stride == NULL) {
            if (ndims > 1) {
                /* least significant dimension */
                *start_off = start[ndims-1];
                *end_off   = start[ndims-1]+(count[ndims-1]-1);
                /* the remaining dimensions */
                for (i=ndims-2; i>0; i--) {
                    *start_off += start[i]*varp->dsizes[i+1];
                    *end_off += (start[i]+(count[i]-1))*varp->dsizes[i+1];
                }
            }
            *start_off *= varp->xsz;  /* offset in bytes */
            *end_off   *= varp->xsz;
            /* handle the unlimited, most significant dimension */
            *start_off += start[0] * ncp->recsize;
            *end_off   += (start[0]+(count[0]-1)) * ncp->recsize;
        }
        else {
            if (ndims > 1) {
                /* least significant dimension */
                *start_off = start[ndims-1];
                *end_off   = start[ndims-1]+(count[ndims-1]-1)*stride[ndims-1];
                /* the remaining dimensions */
                for (i=ndims-2; i>0; i--) {
                    *start_off += start[i]*varp->dsizes[i+1];
                    *end_off += (start[i]+(count[i]-1)*stride[i]) *
                                varp->dsizes[i+1];
                }
            }
            *start_off *= varp->xsz;  /* offset in bytes */
            *end_off   *= varp->xsz;
            /* handle the unlimited, most significant dimension */
            *start_off += start[0] * ncp->recsize;
            *end_off   += (start[0]+(count[0]-1)*stride[0]) * ncp->recsize;
        }
    }
    else {
        if (stride == NULL) {
            /* first handle the least significant dimension */
            *start_off = start[ndims-1];
            *end_off   = start[ndims-1] + (count[ndims-1]-1);
            /* remaining dimensions till the most significant dimension */
            for (i=ndims-2; i>=0; i--) {
                *start_off += start[i] * varp->dsizes[i+1];
                *end_off += (start[i]+(count[i]-1)) * varp->dsizes[i+1];
            }
        }
        else {
            /* first handle the least significant dimension */
            *start_off = start[ndims-1];
            *end_off   = start[ndims-1]+(count[ndims-1]-1)*stride[ndims-1];
            /* remaining dimensions till the most significant dimension */
            for (i=ndims-2; i>=0; i--) {
                *start_off += start[i] * varp->dsizes[i+1];
                *end_off += (start[i]+(count[i]-1)*stride[i])*varp->dsizes[i+1];
            }
        }
        *start_off *= varp->xsz;  /* offset in bytes */
        *end_off   *= varp->xsz;
    }
    *start_off += varp->begin; /* beginning file offset of this variable */
    *end_off   += varp->begin + varp->xsz;

    return NC_NOERR;
}

/*----< wait_getput() >------------------------------------------------------*/
static int
wait_getput(NC         *ncp,
            int         num_reqs,   /* number of non-lead requests */
            NC_req     *reqs,       /* array of non-lead requests */
            int         rw_flag,    /* NC_REQ_WR or NC_REQ_RD */
            int         coll_indep, /* NC_REQ_COLL or NC_REQ_INDEP */
            MPI_Offset  newnumrecs) /* new number of records */
{
    int i, err, status=NC_NOERR, interleaved=0, descreasing=0;
    NC_lead_req *lead_list;

    /* move the offset calculation from request posting API calls to wait call,
     * such that posting a nonblocking request can be made in define mode
     */
    lead_list = (rw_flag == NC_REQ_RD) ? ncp->get_lead_list
                                       : ncp->put_lead_list;
    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead;
        NC_var *varp;

        lead = lead_list + reqs[i].lead_off;
        varp = lead->varp;

        if (varp->ndims == 0) { /* scalar variable */
            reqs[i].offset_start = varp->begin;
            reqs[i].offset_end   = varp->begin + varp->xsz;
        }
        else {
            /* start/count/stride have been allocated in a contiguous array */
            MPI_Offset *count, *stride;
            count  = reqs[i].start + varp->ndims;
            stride = (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) ? NULL :
                     count + varp->ndims;

            /* calculate access range of this request */
            calculate_access_range(ncp, varp, reqs[i].start, count, stride,
                                   &reqs[i].offset_start, &reqs[i].offset_end);
        }
        if (i > 0) {
            /* check if offset_start are in a monotonic nondecreasing order */
            if (reqs[i].offset_start < reqs[i-1].offset_start)
                descreasing = interleaved = 1; /* possible interleaved */
            else if (reqs[i].offset_start < reqs[i-1].offset_end)
                interleaved = 1; /* possible interleaved */
        }
    }

    /* If a decreasing order is found, sort reqs[] based on reqs[].offset_start
     * into an increasing order */
    if (descreasing)
        qsort(reqs, (size_t)num_reqs, sizeof(NC_req), req_compare);

    /* check sorted requests for true interleaved */
    if (interleaved) { /* only if possible interleaved */
        interleaved = 0;
        for (i=1; i<num_reqs; i++) { /* reqs[].offset_start has been sorted */
            if (reqs[i].offset_start < reqs[i-1].offset_end) {
                interleaved = 1;
                break;
            }
        }
    }

    /* aggregate requests and carry out the I/O (reqs will be freed in
     * req_aggregation() */
    err = req_aggregation(ncp, num_reqs, reqs, rw_flag, coll_indep,
                          interleaved);
    if (status == NC_NOERR) status = err;

    /* Update the number of records if new records have been created.
     * For nonblocking APIs, there is no way for a process to know whether
     * others write to a record variable or not. Note newnumrecs has been
     * sync-ed and always >= ncp->numrecs.
     */
    if (rw_flag == NC_REQ_WR) {
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

        if (NC_doFsync(ncp)) { /* NC_SHARE is set */
            int mpireturn;
            if (coll_indep == NC_REQ_INDEP) {
                TRACE_IO(MPI_File_sync)(ncp->independent_fh);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_sync");
                    if (status == NC_NOERR) status = err;
                }
            }
            else {
                TRACE_IO(MPI_File_sync)(ncp->collective_fh);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_sync");
                    if (status == NC_NOERR) status = err;
                }
                TRACE_COMM(MPI_Barrier)(ncp->comm);
            }
        }
    }

    return status;
}

/*----< mgetput() >----------------------------------------------------------*/
/* Before reaching to this subroutine, all the filetypes in the request array
 * are sorted in a non-decreasing order and not interleaved. This subroutine
 * concatenates the filetypes into a single fileview and calls MPI-IO function.
 * This subroutine also concatenates user buffertypes into a single derived
 * data type to be used in the MPI read/write function call.
 */
static int
mgetput(NC     *ncp,
        int     num_reqs,    /* IN: number of requests */
        NC_req *reqs,        /* [num_reqs] non-lead request list, to be freed
                                in this subroutine */
        int     rw_flag,     /* NC_REQ_WR or NC_REQ_RD */
        int     coll_indep)  /* NC_REQ_COLL or NC_REQ_INDEP */
{
    int i, j, len=0, numLeadReqs, status=NC_NOERR, mpireturn, err, *blocklens;
    void *buf=NULL;
    NC_lead_req *lead_list;
    MPI_Datatype filetype, buf_type=MPI_BYTE;
    MPI_Offset offset=0;
    MPI_File fh;
    MPI_Aint *disps;

    blocklens = (int*) NCI_Malloc((size_t)num_reqs * SIZEOF_INT);
    disps = (MPI_Aint*) NCI_Malloc((size_t)num_reqs * SIZEOF_MPI_AINT);

    lead_list = (rw_flag == NC_REQ_RD) ? ncp->get_lead_list
                                       : ncp->put_lead_list;

    /* construct an MPI file type by concatenating fileviews of all requests */
    status = construct_filetypes(ncp, lead_list, num_reqs, blocklens, disps,
                                 reqs, &filetype);
    if (status != NC_NOERR) { /* if failed, skip this request */
        if (coll_indep == NC_REQ_INDEP) {
            NCI_Free(blocklens);
            NCI_Free(disps);
            NCI_Free(reqs);
            return status;
        }

        /* For collective I/O, we still need to participate the successive
           collective calls: setview/read/write */
        filetype = MPI_BYTE;
        buf = NULL;
        len = 0;
        NCI_Free(disps);
        NCI_Free(blocklens);
        goto mpi_io;
    }

    /* now construct buffer datatype */
    if (num_reqs == 1) {
        NC_lead_req *lead = lead_list + reqs[0].lead_off;
        if (fIsSet(lead->flag, NC_REQ_SKIP))
            len = 0;
        else {
            MPI_Offset req_size = reqs[0].nelems * lead->varp->xsz;
            if (req_size > INT_MAX) { /* skip this request */
                if (status == NC_NOERR)
                    DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
                fSet(lead->flag, NC_REQ_SKIP);
                len = 0; /* skip this request */
            }
            else
                len = (int)req_size;
        }
        buf = reqs[0].xbuf;
    }
    else if (num_reqs > 1) { /* create the I/O buffer derived data type */
        int last_contig_req;
        MPI_Aint a0=0, ai, a_last_contig;

        last_contig_req = 0; /* index of the last contiguous request */
        buf = NULL;
        /* process only valid requests */
        for (i=0, j=0; i<num_reqs; i++) {
            MPI_Offset req_size;
            NC_lead_req *lead = lead_list + reqs[i].lead_off;

            if (fIsSet(lead->flag, NC_REQ_SKIP)) continue;

            req_size = reqs[i].nelems * lead->varp->xsz;

            /* check int overflow */
            if (req_size > INT_MAX) { /* int overflows, skip this request */
                if (status == NC_NOERR) /* keep the 1st encountered error */
                    DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
                fSet(lead->flag, NC_REQ_SKIP);
                continue; /* skip this request */
            }
            blocklens[j] = (int)req_size;

            MPI_Get_address(reqs[i].xbuf, &ai);
            if (j == 0) { /* first valid request */
                a_last_contig = a0 = ai;
                buf = reqs[i].xbuf;
            }
            disps[j] = ai - a0;

            req_size = blocklens[last_contig_req];
            req_size += blocklens[j];
            /* if req_size overflows 4-byte int, then skip coalescing */
            if (req_size <= INT_MAX &&
                ai - a_last_contig == blocklens[last_contig_req]) {
                /* user buffer of request j is contiguous from j-1
                 * we coalesce j to j-1 */
                blocklens[last_contig_req] += blocklens[j];
            }
            else if (j > 0) {
                /* not contiguous from request last_contig_req */
                last_contig_req++;
                a_last_contig = ai;
                disps[last_contig_req] = ai - a0;
                blocklens[last_contig_req] = blocklens[i];
            }
            j++;
        }

        /* last_contig_req is the index of last contiguous request */
        if (last_contig_req == 0) {
            /* user buffers can be concatenated into a contiguous buffer */
            buf_type = MPI_BYTE;
            len = blocklens[0];
        }
        else {
            /* after possible concatenating the user buffers, the true number
             * of non-contiguous buffers is last_contig_req+1 */
            int num_contig_reqs = last_contig_req+1;

            /* concatenate buffer addresses into a single buffer type */
            mpireturn = MPI_Type_create_hindexed(num_contig_reqs, blocklens,
                                                 disps, MPI_BYTE, &buf_type);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&buf_type);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }

            len = 1;
        }
    }
    /* if (buf_type == MPI_BYTE) then the whole buf is contiguous */

    NCI_Free(disps);
    NCI_Free(blocklens);

    /* Free up memory space allocated before entering MPI-IO calls, as MPI-IO
     * flattens the fileview and buftype which can take some space. The
     * non-lead request list, reqs, is no longer used after fileview and buftype
     * have been constructed. In addition, the start arrays allocated at lead
     * requests are no longer used.
     */
    numLeadReqs = (rw_flag == NC_REQ_RD) ? ncp->numLeadGetReqs
                                         : ncp->numLeadPutReqs;
    for (i=0; i<numLeadReqs; i++) {
        if (!fIsSet(lead_list[i].flag, NC_REQ_TO_FREE)) continue;
        NCI_Free(lead_list[i].start);
        lead_list[i].start = NULL;
    }

mpi_io:
    NCI_Free(reqs);

    if (coll_indep == NC_REQ_COLL)
        fh = ncp->collective_fh;
    else
        fh = ncp->independent_fh;

    /* set the MPI-IO fileview, this is a collective call */
    err = ncmpio_file_set_view(ncp, fh, &offset, filetype);
    if (filetype != MPI_BYTE) MPI_Type_free(&filetype);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        if (coll_indep == NC_REQ_INDEP) return status;
        len = 0;
    }

    /* call MPI_File_read/MPI_File_write */
    err = ncmpio_read_write(ncp, rw_flag, coll_indep, offset, len, buf_type,
                            buf, ((buf_type == MPI_BYTE) ? 1 : 0));
    if (status == NC_NOERR) status = err;

    if (buf_type != MPI_BYTE) MPI_Type_free(&buf_type);

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                 MPI_INFO_NULL);
     */

    return status;
}
