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
        if (ncp->nprocs > 1)
            TRACE_IO(MPI_File_read_all)(fh, NULL, 0, MPI_BYTE, &mpistatus);
        else
            TRACE_IO(MPI_File_read)(fh, NULL, 0, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_read_all");
            err = (err == NC_EFILE) ? NC_EREAD : err;
            DEBUG_ASSIGN_ERROR(status, err)
        }
    } else { /* write request */
        if (ncp->nprocs > 1)
            TRACE_IO(MPI_File_write_all)(fh, NULL, 0, MPI_BYTE, &mpistatus);
        else
            TRACE_IO(MPI_File_write)(fh, NULL, 0, MPI_BYTE, &mpistatus);
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

/*----< construct_filetypes() >----------------------------------------------*/
/* concatenate the requests into a single MPI derived filetype */
static int
construct_filetypes(NC           *ncp,
                    NC_lead_req  *lead_list, /* NC_REQ_WR or NC_REQ_RD */
                    int           num_reqs,
#ifdef HAVE_MPI_LARGE_COUNT
                    MPI_Count    *blocklens, /* [num_reqs] temp buffer */
                    MPI_Count    *disps,     /* [num_reqs] temp buffer */
#else
                    int          *blocklens, /* [num_reqs] temp buffer */
                    MPI_Aint     *disps,     /* [num_reqs] temp buffer */
#endif
                    NC_req       *reqs,      /* [num_reqs] */
                    MPI_Datatype *filetype)  /* OUT */
{
    int i, j, err, status=NC_NOERR, all_ftype_contig=1, last_contig_req;
    int mpireturn;
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
            if (lead->varp->begin > NC_MAX_INT) {
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
                fSet(lead->flag, NC_REQ_SKIP); /* skip this request */
                if ( lead->status != NULL &&
                    *lead->status == NC_NOERR)
                    *lead->status = err;
                if (status == NC_NOERR)
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
            if (err == NC_NOERR && offset > NC_MAX_INT)
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
            coalesced_len = lead->varp->xsz * reqs[i].nelems;

#ifdef HAVE_MPI_LARGE_COUNT
            blocklens[j] = coalesced_len;
#else
            if (coalesced_len > NC_MAX_INT) {
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
                if (status == NC_NOERR)
                    status = err; /* report first error */
                coalesced_len = 0;
            }
            blocklens[j] = (int)coalesced_len;
#endif
            if (last_contig_req >= 0)
                coalesced_len += blocklens[last_contig_req];
#ifdef HAVE_MPI_LARGE_COUNT
            if (last_contig_req >= 0 &&
                disps[j] - disps[last_contig_req] ==
                blocklens[last_contig_req]) {
                blocklens[last_contig_req] = coalesced_len;
                j--;
            }
            else last_contig_req = j;
#else
            /* if coalesced_len overflows 4-byte int, then skip coalescing */
            if (coalesced_len < NC_MAX_INT && last_contig_req >= 0 &&
                disps[j] - disps[last_contig_req] ==
                blocklens[last_contig_req]) {
                blocklens[last_contig_req] = (int)coalesced_len;
                j--;
            }
            else last_contig_req = j;
#endif
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
        else {
            mpireturn = MPI_Type_dup(ftypes[0], filetype);
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_dup");
        }
    }
    else { /* if (num_reqs > 1 || (num_reqs == 1 && disps[0] > 0)) */
        /* all ftypes[] created fine, now concatenate all ftypes[] */
        if (all_ftype_contig) {
#ifdef HAVE_MPI_LARGE_COUNT
            mpireturn = MPI_Type_create_hindexed_c(num_reqs, blocklens, disps,
                                                   MPI_BYTE, filetype);
#else
            mpireturn = MPI_Type_create_hindexed(num_reqs, blocklens, disps,
                                                 MPI_BYTE, filetype);
#endif
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");
            else {
                MPI_Type_commit(filetype);
                err = NC_NOERR;
            }
        }
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            mpireturn = MPI_Type_create_struct_c(num_reqs, blocklens, disps,
                                                 ftypes, filetype);
#else
            mpireturn = MPI_Type_create_struct(num_reqs, blocklens, disps,
                                               ftypes, filetype);
#endif
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
#ifdef HAVE_MPI_LARGE_COUNT
                      MPI_Count    *blocklens, /* [num_reqs] temp buffer */
                      MPI_Count    *disps,     /* [num_reqs] temp buffer */
#else
                      int          *blocklens, /* [num_reqs] temp buffer */
                      MPI_Aint     *disps,     /* [num_reqs] temp buffer */
#endif
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

#ifdef HAVE_MPI_LARGE_COUNT
        blocklens[j] = req_size;
#else
        /* check int overflow */
        if (req_size > NC_MAX_INT) { /* skip this request */
            fSet(lead->flag, NC_REQ_SKIP);
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
            continue;
        }
        blocklens[j] = (int)req_size;
#endif

        MPI_Get_address(reqs[i].xbuf, &ai);
        if (j == 0) a0 = ai;
        disps[j] = MPI_Aint_diff(ai, a0);
        j++;
    }
    /* update num_reqs to number of valid requests */
    num_reqs = j;

    if (num_reqs > 0) {
        /* concatenate buffer addresses into a single buffer type */
#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_create_hindexed_c(num_reqs, blocklens, disps,
                                               MPI_BYTE, buf_type);
#else
        mpireturn = MPI_Type_create_hindexed(num_reqs, blocklens, disps,
                                             MPI_BYTE, buf_type);
#endif
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

        if (ncp->aggregation && coll_indep == NC_REQ_COLL && ncp->nprocs > 1)
            /* intra-node write aggregation must be in collective mode */
            err = ncmpio_intra_node_aggregation(ncp, num_w_reqs, put_list, newnumrecs);
        else
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
             MPI_Offset   offset,  /* starting file offset of variable */
             MPI_Offset  *dimlen,  /* [ndim] dimension lengths */
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
/* flatten all requests into offset-length pairs, sort them into an increasing
 * order, and resolve the overlapped offset-length pairs.
 */
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
        ndims = lead->varp->ndims;
        if (ndims > 0) {
            start  = reqs[i].start;
            count  = start + ndims;
            stride = count + ndims;
        }
        else
            start = count = stride = NULL;

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
        addr = MPI_Aint_diff(addr, buf_addr);  /* distance to the buf of first req */

        ndims = lead->varp->ndims;
        if (ndims > 0) {
            start  = reqs[i].start;
            count  = start + ndims;
            stride = count + ndims;
        }
        else
            start = count = stride = NULL;

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
        vars_flatten(ndims, lead->varp->xsz, var_begin, shape,
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
            if (MPI_Aint_add((*segs)[i].buf_addr, (*segs)[i].len) ==
                MPI_Aint_add((*segs)[j].buf_addr, gap)) {
                /* buffers i and j are contiguous, merge j to i */
                (*segs)[i].len = MPI_Aint_add((*segs)[i].len, (*segs)[j].len - gap);
            }
            else { /* buffers are not contiguous, reduce j's len */
                (*segs)[i+1].off      = (*segs)[j].off + gap;
                (*segs)[i+1].len      = (*segs)[j].len - gap;
                (*segs)[i+1].buf_addr = MPI_Aint_add((*segs)[j].buf_addr, gap);
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
    int i, j, mpireturn;
    MPI_Offset next_off, next_len, true_nsegs;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *blocklens;
    MPI_Count *disps;
#else
    int *blocklens;
    MPI_Aint *disps;
#endif

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
#ifdef HAVE_MPI_LARGE_COUNT
    blocklens = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * true_nsegs);
    disps     = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * true_nsegs);

    /* coalesce segs[].off and len to disps[] and blocklens[] */
    disps[0]     = segs[0].off;
    blocklens[0] = segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (disps[j] + blocklens[j] == segs[i].off)
            /* j and i are contiguous */
            blocklens[j] += segs[i].len;
        else {
            j++;
            disps[j]     = segs[i].off;
            blocklens[j] = segs[i].len;
        }
    }
    /* Now j+1 is the coalesced length */

    mpireturn = MPI_Type_create_hindexed_c(j+1, blocklens, disps, MPI_BYTE,
                                           filetype);
#else
    blocklens = (int*)      NCI_Malloc(true_nsegs * SIZEOF_INT);
    disps     = (MPI_Aint*) NCI_Malloc(true_nsegs * SIZEOF_MPI_AINT);

    /* coalesce segs[].off and len to disps[] and blocklens[] */
    if (segs[0].len > NC_MAX_INT) {
        NCI_Free(disps);
        NCI_Free(blocklens);
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    }
    disps[0]     =      segs[0].off;
    blocklens[0] = (int)segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (segs[i].len > NC_MAX_INT) {
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
    /* Now j+1 is the coalesced length */

    mpireturn = MPI_Type_create_hindexed(j+1, blocklens, disps, MPI_BYTE,
                                         filetype);
#endif
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

#ifdef HAVE_MPI_LARGE_COUNT
    if (true_nsegs < j + 1) {
        blocklens = (MPI_Count*) NCI_Realloc(blocklens, (j+1) * sizeof(MPI_Count));
        disps     = (MPI_Count*) NCI_Realloc(disps,     (j+1) * sizeof(MPI_Count));
    }

    /* coalesce segs[].off and len to disps[] and blocklens[] */
    disps[0]     = segs[0].buf_addr;
    blocklens[0] = segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (disps[j] + blocklens[j] == segs[i].buf_addr)
            /* j and i are contiguous */
            blocklens[j] += segs[i].len;
        else {
            j++;
            disps[j]     = segs[i].buf_addr;
            blocklens[j] = segs[i].len;
        }
    }
    /* j+1 is the coalesced length */

    mpireturn = MPI_Type_create_hindexed_c(j+1, blocklens, disps, MPI_BYTE,
                                           buf_type);
#else
    if (true_nsegs < j + 1) {
        blocklens = (int*)      NCI_Realloc(blocklens, (j+1) * SIZEOF_INT);
        disps     = (MPI_Aint*) NCI_Realloc(disps,     (j+1) * SIZEOF_MPI_AINT);
    }

    /* coalesce segs[].off and len to disps[] and blocklens[] */
    if (segs[0].len > NC_MAX_INT) {
        NCI_Free(disps);
        NCI_Free(blocklens);
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    }
    disps[0]     =      segs[0].buf_addr;
    blocklens[0] = (int)segs[0].len;
    for (j=0,i=1; i<nsegs; i++) {
        if (segs[i].len > NC_MAX_INT) {
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
#endif
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
    int *group_index, *group_type, numLeadReqs;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *blocklens, *f_blocklens, *b_blocklens;
    MPI_Count *disps, *f_disps, *b_disps;
#else
    int *blocklens, *f_blocklens, *b_blocklens;
    MPI_Aint *disps, *f_disps, *b_disps;
#endif
    NC_lead_req *lead_list;
    void *buf; /* point to starting buffer, used by MPI-IO call */
    MPI_Aint      b_begin, b_addr;
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
     * requests. Codes below flatten the requests of "interleaved" groups
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

    /* temp buffers, used by multiple calls to construct_filetypes()  */
#ifdef HAVE_MPI_LARGE_COUNT
    blocklens   = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * num_reqs);
    disps       = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * num_reqs);
    f_blocklens = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * ngroups);
    f_disps     = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * ngroups);
    b_blocklens = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * ngroups);
    b_disps     = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * ngroups);
#else
    blocklens   = (int*)      NCI_Malloc(sizeof(int)      * num_reqs);
    disps       = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * num_reqs);
    f_blocklens = (int*)      NCI_Malloc(sizeof(int)      * ngroups);
    f_disps     = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ngroups);
    b_blocklens = (int*)      NCI_Malloc(sizeof(int)      * ngroups);
    b_disps     = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ngroups);
#endif

    buf = reqs[0].xbuf; /* the buffer of 1st request */
    b_disps[0] = 0;     /* relative to address of 1st buf */
    MPI_Get_address(buf, &b_begin);

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
             * non-decreasing order for constructing the filetype. For example,
             * multiple nonblocking requests each writing/reading a single
             * column of a 2D array will produces an interleaving filetype
             * among all processes.
             *
             * The pitfall of this flattening is the additional memory
             * requirement, as it will have to break down each request into a
             * list of offset-length pairs, and then merge all lists into a
             * sorted list based on their offsets into an increasing order.
             *
             * Be warned! The additional memory requirement for this merging can
             * be more than the I/O data itself. For example, in the
             * column-wise data partitioning pattern, a process makes a
             * nonblocking request for accessing a single column of a 2D array
             * of 4-byte integer type. However, each element of the column is
             * flattened into an off-len pair and this off-len pair itself
             * takes 24 bytes, sizeof(struct off_len). Additional memory is
             * also required for MPI arguments of displacements and
             * blocklengths when constructing the filetype.
             */
            MPI_Offset  nsegs=0;   /* number of merged offset-length pairs */
            off_len    *segs=NULL; /* array of the offset-length pairs */
            void       *merged_buf;

            /* flatten and merge all requests into sorted offset-length pairs.
             * Note g_reqs[].offset_start and offset_end are relative to the
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
            b_disps[i] = MPI_Aint_diff(b_addr, b_begin); /* to 1st buffer of 1st group*/
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
#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_create_struct_c(ngroups, f_blocklens, f_disps,
                                             ftypes, &filetype);
#else
        mpireturn = MPI_Type_create_struct(ngroups, f_blocklens, f_disps,
                                           ftypes, &filetype);
#endif
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
#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_create_struct_c(ngroups, b_blocklens, b_disps,
                                             btypes, &buf_type);
#else
        mpireturn = MPI_Type_create_struct(ngroups, b_blocklens, b_disps,
                                           btypes, &buf_type);
#endif
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
    NCI_Free(b_blocklens);
    NCI_Free(b_disps);

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

    fh = ncp->independent_fh;
    if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
        fh = ncp->collective_fh;

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

#if 0
/* C struct for MPI messages storing a list of offset-length pairs */
typedef struct {
    MPI_Offset off;      /* starting file offset of the request */
    MPI_Offset len;      /* requested length in bytes starting from off */
} off_len_msg;

/*----< flatten_subarray() >-------------------------------------------------*/
/* flatten a subarray request into a list of offset-length pairs */
static int
flatten_subarray(int          ndim,    /* number of dimensions */
                 int          el_size, /* array element size */
                 MPI_Offset   offset,  /* starting file offset of variable */
                 MPI_Offset  *dimlen,  /* [ndim] dimension lengths */
                 MPI_Offset  *start,   /* [ndim] starts of subarray */
                 MPI_Offset  *count,   /* [ndim] counts of subarray */
                 MPI_Offset  *stride,  /* [ndim] strides of subarray */
                 MPI_Offset  *npairs,  /* OUT: number of off-len pairs */
                 MPI_Offset  *offsets, /* OUT: array of offsets */
                 MPI_Offset  *lengths) /* OUT: array of lengths */
{
    int i, j, to_free_stride=0;
    MPI_Offset length, nstride, array_len, off, subarray_len;
    size_t idx=0, idx0;

    *npairs = 0;
    if (ndim < 0) return NC_NOERR;

    if (ndim == 0) {  /* scalar record variable */
        *npairs = 1;
        offsets[0] = offset;
        lengths[0] = el_size;
        return NC_NOERR;
    }

    if (stride == NULL) { /* equivalent to {1, 1, ..., 1} */
        stride = (MPI_Offset*) NCI_Malloc((size_t)ndim * SIZEOF_MPI_OFFSET);
        for (i=0; i<ndim; i++) stride[i] = 1;
        to_free_stride = 1;
    }

    /* TODO: check if all stride[] >= 1
       Q: Is it legal if any stride[] <= 0 ? */

    /* calculate the number of offset-length pairs */
    *npairs = (stride[ndim-1] == 1) ? 1 : count[ndim-1];
for (i=0; i<ndim-1; i++) if (count[i] == 0) printf("++++++++++++++++++++++++ ERROR %d count[%d]=0\n",__LINE__,i);
    for (i=0; i<ndim-1; i++)
        *npairs *= count[i];
    if (*npairs == 0) {  /* not reachable, an error if count[] == 0 */
        if (to_free_stride) NCI_Free(stride);
        return NC_NOERR;
    }

    /* length of each row of the subarray are of the same size */
    length  = (stride[ndim-1] == 1) ? count[ndim-1] : 1;
    length *= el_size;
    nstride  = (stride[ndim-1] == 1) ? 1 : count[ndim-1];

    /* set the offset-length pairs for the lowest dimension */
    off = offset + start[ndim-1] * el_size;
    for (i=0; i<nstride; i++) {
        offsets[idx]  = off;
        lengths[idx]  = length;
        off          += stride[ndim-1] * el_size;
        idx++;
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
        idx0 = 0;
        for (j=0; j<subarray_len; j++) {
            offsets[idx0] += off;
            idx0++;
        }

        /* update each plan subarray of dimension ndim-1 */
        off = array_len * stride[ndim-1] * el_size;
        for (i=1; i<count[ndim-1]; i++) {
            idx0 = 0;
            for (j=0; j<subarray_len; j++) {
                offsets[idx] = offsets[idx0] + off;
                lengths[idx] = length;
                idx++;
                idx0++;
            }
            off += array_len * stride[ndim-1] * el_size;
        }
        ndim--;  /* move to next higher dimension */
        subarray_len *= count[ndim];
    }
    if (to_free_stride) NCI_Free(stride);

    return NC_NOERR;
}

/*----< flatten_reqs() >-----------------------------------------------------*/
/* flatten all write requests into offset-length pairs.
 * Variable pairs is allocated here and need to be freed by the caller
 */
static int
flatten_reqs(NC            *ncp,
             int            num_reqs,  /* IN: # requests */
             const NC_req  *reqs,      /* [num_reqs] requests */
             MPI_Aint      *num_pairs, /* OUT: total number of off-len pairs */
             MPI_Offset   **offsets,   /* OUT: array of flattened offsets */
             MPI_Offset   **lengths)   /* OUT: array of flattened lengths */
{
    int i, j, status=NC_NOERR, ndims;
    MPI_Offset idx, num, *start, *count, *shape, *stride;

    *num_pairs = 0;    /* total number of offset-length pairs */

    /* Count the number off-len pairs from reqs[], so we can malloc a
     * contiguous memory space for storing off-len pairs
     */
    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;
        ndims = lead->varp->ndims;
        if (ndims > 0) {
            start  = reqs[i].start;
            count  = start + ndims;
            stride = count + ndims;
        }
        else
            start = count = stride = NULL;

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
            (*num_pairs)++;
            continue;
        }
        num = 1;
        if (stride != NULL && stride[ndims-1] > 1)
            num = count[ndims-1];  /* count of last dimension */
        for (j=0; j<ndims-1; j++)
            num *= count[j];  /* all count[] except the last dimension */

        (*num_pairs) += num;
    }

    /* now we can allocate a contiguous memory space for the off-len pairs */
    *offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (*num_pairs) * 2);
    *lengths = *offsets + (*num_pairs);
    idx = 0;

    /* now re-run the loop to fill in the off-len pairs */
    for (i=0; i<num_reqs; i++) {
        MPI_Offset var_begin;
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;

        ndims = lead->varp->ndims;
        if (ndims > 0) {
            start  = reqs[i].start;
            count  = start + ndims;
            stride = count + ndims;
        }
        else
            start = count = stride = NULL;

        shape = lead->varp->shape;

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

// printf("---------------- %d num_reqs=%d i=%d ndims=%d count=%lld\n",ncp->rank,num_reqs, i,ndims, count[0]);
        /* flatten each request into a list of offset-length pairs and
         * append to the end of offsets and lengths
         */
        flatten_subarray(ndims, lead->varp->xsz, var_begin, shape,
                         start, count, stride,
                         &num,            /* OUT: number of off-len pairs */
                         *offsets + idx,  /* OUT: array of offsets */
                         *lengths + idx); /* OUT: array of lengths */
        idx += num;
// printf("---------------- %d num_reqs=%d i=%d ndims=%d num=%lld\n",ncp->rank,num_reqs, i,ndims, num);

    }

    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;
        if (fIsSet(lead->flag, NC_REQ_TO_FREE)) {
            NCI_Free(lead->start);
            lead->start = NULL;
        }
    }

    return status;
}

/*----< construct_buf_type() >-----------------------------------------------*/
/* construct an MPI derived datatype for I/O buffers from the request list, by
 * concatenate all buffers.
 */
static int
construct_buf_type(const NC     *ncp,
                   int           num_reqs,  /* IN: # requests */
                   const NC_req *reqs,      /* [num_reqs] requests */
                   MPI_Aint     *bufLen,    /* OUT: buffer size in bytes */
                   MPI_Datatype *bufType)   /* OUT: buffer datatype */
{
    int i, err, mpireturn, status=NC_NOERR;
    NC_lead_req *lead;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *blocklens = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * num_reqs);
    MPI_Count *disps     = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * num_reqs);
#else
    int       *blocklens = (int*)       NCI_Malloc(sizeof(int)       * num_reqs);
    MPI_Aint  *disps     = (MPI_Aint*)  NCI_Malloc(sizeof(MPI_Aint)  * num_reqs);
#endif

    *bufLen = 0;
    for (i=0; i<num_reqs; i++) {
        MPI_Aint addr;

        /* displacement uses MPI_BOTTOM */
        MPI_Get_address(reqs[i].xbuf, &addr);
        disps[i] = addr;

        /* blocklens[] in bytes */
        lead = ncp->put_lead_list + reqs[i].lead_off;
        blocklens[i] = reqs[i].nelems * lead->varp->xsz;

        *bufLen += blocklens[i];
    }

    /* construct buffer derived datatype */
#ifdef HAVE_MPI_LARGE_COUNT
    mpireturn = MPI_Type_create_hindexed_c(num_reqs, blocklens, disps,
                                           MPI_BYTE, bufType);
#else
    mpireturn = MPI_Type_create_hindexed(num_reqs, blocklens, disps,
                                         MPI_BYTE, bufType);
#endif
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");
        /* return the first encountered error if there is any */
        if (status == NC_NOERR) status = err;

        *bufType = MPI_DATATYPE_NULL;
    }
    else
        MPI_Type_commit(bufType);

    NCI_Free(blocklens);
    NCI_Free(disps);

    return status;
}

/*----< intra_node_aggregation() >------------------------------------------*/
/* This is a collective call */
static int
intra_node_aggregation(NC     *ncp,
                       int     num_reqs,
                       NC_req *put_list,
                       MPI_Offset newnumrecs)
{
    int i, j, err, mpireturn, status=NC_NOERR, nreqs;
    char *recv_buf=NULL, *wr_buf = NULL;
    MPI_Aint npairs, *msg, bufLen, num_pairs;
    MPI_Offset offset=0, buf_count=0, *offsets, *lengths;
    MPI_Datatype fileType=MPI_BYTE, bufType=MPI_BYTE;
    MPI_File fh;
    MPI_Request *req;

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     */
    flatten_reqs(ncp, num_reqs, put_list, &num_pairs, &offsets, &lengths);

// #define WKL_DEBUG
#ifdef WKL_DEBUG
for (j=1; j<num_pairs; j++) {
    if (offsets[j-1] + lengths[j-1] > offsets[j]) {
        printf("%d: DECREASING num_pairs=%ld offsets[%d]=%lld len=%lld = %lld > offsets[%d]=%lld\n",ncp->rank,num_pairs,j-1,offsets[j-1], lengths[j-1], offsets[j-1] + lengths[j-1], j, offsets[j]);
        break;
    }
}
for (j=0; j<num_pairs; j++) {
if (offsets[j] == 868) printf("++++++++++++++ %d: offsets[%d]=868 len=%lld\n",ncp->rank,j,lengths[j]);
}
for (j=0; j<num_pairs; j++) {
if (offsets[j] == 233588) printf("++++++++++++++ %d: offsets[%d]=233588 len=%lld\n",ncp->rank,j,lengths[j]);
if (offsets[j] == 233628) printf("++++++++++++++ %d: offsets[%d]=233628\n",ncp->rank,j);
}

printf("%d: %s:%d num_reqs=%d num_pairs=%ld offsets=%lld %lld %lld %lld lengths=%lld %lld %lld %lld\n",ncp->rank,__func__,__LINE__,num_reqs,num_pairs,
offsets[0], offsets[1], offsets[2], offsets[3],
lengths[0], lengths[1], lengths[2], lengths[3]);
#endif

    /* construct write buffer datatype, bufType.
     * bufLen is the buffer size in bytes
     */
    construct_buf_type(ncp, num_reqs, put_list, &bufLen, &bufType);

    NCI_Free(put_list);

    /* First, tell aggregator how much to receive by sending:
     * (num_pairs and bufLen). The message size to be sent by this rank
     * is num_pairs * 2 * sizeof(MPI_Offset) + bufLen
     */
    nreqs = (ncp->rank == ncp->my_aggr) ? (ncp->num_non_aggrs - 1) : 1;

    msg = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_non_aggrs * 2);
    msg[0] = num_pairs;
    msg[1] = bufLen;
#ifdef WKL_DEBUG
printf("%s:%d rank %d num_pairs=%ld bufLen=%ld\n",__func__,__LINE__,ncp->rank,num_pairs,bufLen);
#endif

    req = (MPI_Request*) NCI_Malloc(sizeof(MPI_Request) * nreqs);
    if (ncp->rank == ncp->my_aggr) {
        for (i=1; i<ncp->num_non_aggrs; i++)
            MPI_Irecv(msg + i*2, 2, MPI_AINT, ncp->nonaggr_ranks[i], 0, ncp->comm, &req[i-1]);
    }
    else
        MPI_Isend(msg, 2, MPI_AINT, ncp->my_aggr, 0, ncp->comm, &req[0]);

    MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);

    if (ncp->rank == ncp->my_aggr) {
#ifdef WKL_DEBUG
for (i=0; i<ncp->num_non_aggrs; i++) printf("%s:%d from %d num_pairs=%ld bufLen=%ld\n",__func__,__LINE__,i,msg[2*i],msg[2*i+1]);
#endif
        /* calculate the total number of offset-length pairs */
        npairs = num_pairs;
        for (i=1; i<ncp->num_non_aggrs; i++) npairs += msg[i*2];

        /* realloc to store all pairs in a contiguous buffer */
        offsets = (MPI_Offset*) NCI_Realloc(offsets, sizeof(MPI_Offset) * npairs * 2);
        lengths = offsets + num_pairs;

        /* post requests to receive offset-length pairs from non-aggregators */
        MPI_Offset *ptr = offsets + num_pairs * 2;
        for (i=1; i<ncp->num_non_aggrs; i++) {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Irecv_c(ptr, msg[i*2] * 2, MPI_OFFSET, ncp->nonaggr_ranks[i], 0, ncp->comm, &req[i-1]);
#else
            MPI_Irecv(ptr, msg[i*2] * 2, MPI_OFFSET, ncp->nonaggr_ranks[i], 0, ncp->comm, &req[i-1]);
#endif
            ptr += msg[i*2] * 2;
        }
#ifdef WKL_DEBUG
printf("%s:%d offsets=%lld %lld %lld %lld lengths=%lld %lld %lld %lld\n",__func__,__LINE__,
offsets[0], offsets[1], offsets[2], offsets[3],
lengths[0], lengths[1], lengths[2], lengths[3]);
#endif
    }
    else {
        /* send offset-length pairs data to the aggregator */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Isend_c(offsets, num_pairs * 2, MPI_OFFSET, ncp->my_aggr, 0, ncp->comm, &req[0]);
#else
        MPI_Isend(offsets, num_pairs * 2, MPI_OFFSET, ncp->my_aggr, 0, ncp->comm, &req[0]);
#endif
        NCI_Free(msg);
    }

    MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);

    /*
     * TODO, define a datatype to combine sends of offset-length pairs with the
     * write data into a single send call.
     */
    if (ncp->rank == ncp->my_aggr) {
        /* calculate the total write account */
        buf_count = bufLen;
        for (i=1; i<ncp->num_non_aggrs; i++) buf_count += msg[i*2 + 1];

        /* Allocate receive buffer, which will be sorted into an increasing
         * order based on the file offsets. Thus, after sorting pack recv_buf
         * to wr_buf to avoid creating another buffer datatype.
         */
        recv_buf = (char*) NCI_Malloc(buf_count);
        wr_buf = (char*) NCI_Malloc(buf_count);

        /* First, pack self write data into front of the recv_buf */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count position=0;
        MPI_Pack_c(MPI_BOTTOM, 1, bufType, recv_buf, bufLen, &position, MPI_COMM_SELF);
#else
        int position=0;
        MPI_Pack(MPI_BOTTOM, 1, bufType, recv_buf, bufLen, &position, MPI_COMM_SELF);
#endif
        /* post requests to receive write data from non-aggregators */
        char *ptr = recv_buf + bufLen;
        for (i=1; i<ncp->num_non_aggrs; i++) {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Irecv_c(ptr, msg[i*2 + 1], MPI_BYTE, ncp->nonaggr_ranks[i], 0, ncp->comm, &req[i-1]);
#else
            MPI_Irecv(ptr, msg[i*2 + 1], MPI_BYTE, ncp->nonaggr_ranks[i], 0, ncp->comm, &req[i-1]);
#endif
            ptr += msg[i*2 + 1];
        }
    }
    else {
        /* send write data to the aggregator */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Isend_c(MPI_BOTTOM, 1, bufType, ncp->my_aggr, 0, ncp->comm, &req[0]);
#else
        MPI_Isend(MPI_BOTTOM, 1, bufType, ncp->my_aggr, 0, ncp->comm, &req[0]);
#endif
        NCI_Free(offsets);
    }

    MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
    NCI_Free(req);

    /* free and reset bufType */
    if (bufType != MPI_BYTE && bufType != MPI_DATATYPE_NULL)
        MPI_Type_free(&bufType);
    bufType = MPI_BYTE;

    /* aggregator sorts the offset-length pairs, along with the buffer */
    if (ncp->rank == ncp->my_aggr) {
        off_len *segs = (off_len*)NCI_Malloc(sizeof(off_len) * npairs);
        MPI_Offset *off_ptr = offsets;
        int k=0;
        for (i=0; i<ncp->num_non_aggrs; i++) {
            MPI_Offset *len_ptr = off_ptr + msg[i*2];
            for (j=0; j<msg[i*2]; j++) {
                segs[k].off = off_ptr[j];
                segs[k].len = len_ptr[j];
                if (k == 0)
                    segs[k].buf_addr = 0;
                else
                    segs[k].buf_addr = segs[k-1].buf_addr + segs[k-1].len;
                k++;
            }
            off_ptr += msg[i*2] * 2;
        }
        NCI_Free(msg);
        qsort(segs, npairs, sizeof(off_len), off_compare);

#ifdef WKL_DEBUG
for (j=1; j<npairs; j++) {
    if (segs[j-1].off >= segs[j].off) {
        printf("zzzzzzzzzzz npairs=%ld segs[%d].off=%lld >= segs[%d].off=%lld\n",npairs,j-1,segs[j-1].off, j, segs[j].off);
        break;
    }
}
for (j=1; j<npairs; j++) {
    if (segs[j-1].off + segs[j-1].len > segs[j].off) {
        printf("ooooooooverlap npairs=%ld segs[%d].off=%lld len=%lld =%lld> segs[%d].off=%lld len=%lld\n",npairs,j-1,segs[j-1].off, segs[j-1].len, segs[j-1].off + segs[j-1].len,j, segs[j].off,segs[j].len);
        break;
    }
}

printf("%d: %s:%d npairs=%ld offsets=%lld %lld %lld %lld lengths=%lld %lld %lld %lld buf_addr=%ld %ld %ld %ld\n",ncp->rank,__func__,__LINE__,npairs,
segs[0].off, segs[1].off, segs[2].off, segs[3].off,
segs[0].len, segs[1].len, segs[2].len, segs[3].len,
segs[0].buf_addr, segs[1].buf_addr, segs[2].buf_addr, segs[3].buf_addr);
printf("%d: %s:%d npairs=%ld offsets[40]=%lld %lld %lld %lld lengths=%lld %lld %lld %lld buf_addr=%ld %ld %ld %ld\n",ncp->rank,__func__,__LINE__,npairs,
segs[40].off, segs[41].off, segs[42].off, segs[43].off,
segs[40].len, segs[41].len, segs[42].len, segs[43].len,
segs[40].buf_addr, segs[41].buf_addr, segs[42].buf_addr, segs[43].buf_addr);
#endif

        NCI_Free(offsets);

#if 1
        /* merge the overlapped buffer segments, skip the overlapped regions
         * for those with higher j indices (i.e. requests with lower j indices
         * win the writes to the overlapped regions)
         */
#ifdef WKL_DEBUG
MPI_Aint wkl = segs[0].len;
for (i=0, j=1; j<npairs; j++) { wkl += segs[j].len; }
#endif
        for (i=0, j=1; j<npairs; j++) {
            if (segs[i].off + segs[i].len >= segs[j].off + segs[j].len)
                /* segment i completely covers segment j, skip j */
                continue;

            MPI_Offset gap = segs[i].off + segs[i].len - segs[j].off;
// if (gap > 0) printf("xxxxxxxxx segments %d off=%lld len=%lld = %lld and %d off=%lld overlaps\n",i,segs[i].off, segs[i].len, segs[i].off + segs[i].len, j, segs[j].off);
            if (gap >= 0) { /* segments i and j overlaps */
                if (MPI_Aint_add(segs[i].buf_addr, segs[i].len) ==
                    MPI_Aint_add(segs[j].buf_addr, gap)) {
                    /* buffers i and j are contiguous, merge j to i */
                    segs[i].len = MPI_Aint_add(segs[i].len, segs[j].len - gap);
                }
                else { /* buffers are not contiguous, reduce j's len */
                    segs[i+1].off      = segs[j].off + gap;
                    segs[i+1].len      = segs[j].len - gap;
                    segs[i+1].buf_addr = MPI_Aint_add(segs[j].buf_addr, gap);
                    i++;
                }
            }
            else { /* i and j do not overlap */
                i++;
                if (i < j) segs[i] = segs[j];
            }
        }

#ifdef WKL_DEBUG
printf("------------- %s: %d npairs=%ld i+1=%d wkl=%ld buf_count=%lld\n",__func__,__LINE__,npairs, i+1,wkl,buf_count);
wkl = 0;
for (j=0; j<i+1; j++) {
    wkl += segs[j].len;
    if (segs[j].buf_addr + segs[j].len > buf_count) {
        printf("Error: segs[%d].buf_addr=%ld + segs[i].len=%lld > buf_count=%lld\n",j,segs[j].buf_addr,segs[j].len,buf_count);
        break;
   }
}
printf("------------- %s: %d npairs=%ld i+1=%d wkl=%ld buf_count=%lld\n",__func__,__LINE__,npairs, i+1,wkl,buf_count);
#endif

        /* update number of segments, now all off-len pairs are not overlapped */
        npairs = i+1;
#endif

        /* pack recv_buf into wr_buf, a contiguous buffer */
        char *ptr = wr_buf;
        for (i=0; i<npairs; i++) {
            memcpy(ptr, recv_buf + segs[i].buf_addr, segs[i].len);
            ptr += segs[i].len;
        }
        NCI_Free(recv_buf);
        bufType = MPI_BYTE;

        /* coalesce the offset-length pairs */
        for (i=0, j=1; j<npairs; j++) {
            if (segs[i].off + segs[i].len == segs[j].off) {
                /* coalesce j into i */
                segs[i].len += segs[j].len;
            }
            else {
                i++;
                if (i < j) segs[i] = segs[j];
            }
        }

#ifdef WKL_DEBUG
printf("%s: %d npairs=%ld i+1=%d\n",__func__,__LINE__,npairs, i+1);
#endif
        /* update number of segments, now all off-len pairs are not overlapped */
        npairs = i+1;

        if (npairs == 1) {
            offset = segs[0].off;
#ifdef WKL_DEBUG
printf("%s: %d npairs=%ld i+1=%d offset=%lld\n",__func__,__LINE__,npairs, i+1,offset);
#endif
        }
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            /* construct fileview */
            MPI_Count *blocklens = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * npairs);
            MPI_Count *disps     = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * npairs);
            for (i=0; i<npairs; i++) {
                disps[i]     = segs[i].off;
                blocklens[i] = segs[i].len;
            }

            mpireturn = MPI_Type_create_hindexed_c(npairs, blocklens, disps, MPI_BYTE,
                                                   &fileType);

            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&fileType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }

#if 0
            /* construct buftype */
            for (i=0; i<npairs; i++)
                disps[i] = segs[i].buf_addr;

            mpireturn = MPI_Type_create_hindexed_c(npairs, blocklens, disps, MPI_BYTE,
                                                   &bufType);

            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&bufType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
#endif
            NCI_Free(blocklens);
            NCI_Free(disps);
#else
            /* construct fileview */
            int      *blocklens = (int*)      NCI_Malloc(SIZEOF_INT * npairs);
            MPI_Aint *disps     = (MPI_Aint*) NCI_Malloc(SIZEOF_MPI_AINT * npairs);
            for (i=0; i<npairs; i++) {
                if (segs[i].len > NC_MAX_INT) {
                    NCI_Free(disps);
                    NCI_Free(blocklens);
                    DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
                }
                disps[i]     = segs[i].off;
                blocklens[i] = segs[i].len;
            }

            mpireturn = MPI_Type_create_hindexed(npairs, blocklens, disps, MPI_BYTE,
                                                 &fileType);

            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&fileType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }

#if 0
            /* construct buftype */
            for (i=0; i<npairs; i++)
                disps[i] = segs[i].buf_addr;

            mpireturn = MPI_Type_create_hindexed(npairs, blocklens, disps, MPI_BYTE,
                                                 &bufType);

            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&bufType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
#endif
            NCI_Free(blocklens);
            NCI_Free(disps);
#endif
        }
        NCI_Free(segs);
    }

    /* aggregator writes to the file */
    /* non-aggregators participate collective write call with zero-lenghth write */
    fh = ncp->collective_fh;

    /* set the MPI-IO fileview, this is a collective call */
    err = ncmpio_file_set_view(ncp, fh, &offset, fileType);
    if (fileType != MPI_BYTE) MPI_Type_free(&fileType);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        buf_count = 0;
    }

    /* call MPI_File_write_at_all */
    err = ncmpio_read_write(ncp, NC_REQ_WR, NC_REQ_COLL, offset, buf_count,
                            bufType, wr_buf, ((bufType == MPI_BYTE) ? 1 : 0));
    if (bufType != MPI_BYTE) MPI_Type_free(&bufType);
    if (status == NC_NOERR) status = err;

    if (wr_buf != NULL) NCI_Free(wr_buf);

    /* Update the number of records if new records have been created.
     * For nonblocking APIs, there is no way for a process to know whether
     * others write to a record variable or not. Note newnumrecs has been
     * sync-ed and always >= ncp->numrecs.
     */
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

    return status;
}
#endif

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
    int i, j, numLeadReqs, status=NC_NOERR, mpireturn, err;
    void *buf=NULL;
    NC_lead_req *lead_list;
    MPI_Datatype filetype, buf_type=MPI_BYTE;
    MPI_Offset offset=0, buf_count=0;
    MPI_File fh;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *blocklens;
    MPI_Count *disps;
    blocklens = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * num_reqs);
    disps = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * num_reqs);
#else
    int *blocklens;
    MPI_Aint *disps;
    blocklens = (int*) NCI_Malloc((size_t)num_reqs * SIZEOF_INT);
    disps = (MPI_Aint*) NCI_Malloc((size_t)num_reqs * SIZEOF_MPI_AINT);
#endif

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
        buf_count = 0;
        NCI_Free(disps);
        NCI_Free(blocklens);
        goto mpi_io;
    }

    /* now construct buffer datatype */
    if (num_reqs == 1) {
        NC_lead_req *lead = lead_list + reqs[0].lead_off;
        if (fIsSet(lead->flag, NC_REQ_SKIP))
            buf_count = 0;
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            buf_count = reqs[0].nelems * lead->varp->xsz;
#else
            MPI_Offset req_size = reqs[0].nelems * lead->varp->xsz;
            if (req_size > NC_MAX_INT) { /* skip this request */
                if (status == NC_NOERR)
                    DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
                fSet(lead->flag, NC_REQ_SKIP);
                buf_count = 0; /* skip this request */
            }
            else
                buf_count = req_size;
#endif
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

#ifdef HAVE_MPI_LARGE_COUNT
            blocklens[j] = req_size;
#else
            /* check int overflow */
            if (req_size > NC_MAX_INT) { /* int overflows, skip this request */
                if (status == NC_NOERR) /* keep the 1st encountered error */
                    DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
                fSet(lead->flag, NC_REQ_SKIP);
                continue; /* skip this request */
            }
            blocklens[j] = (int)req_size;
#endif

            MPI_Get_address(reqs[i].xbuf, &ai);
            if (j == 0) { /* first valid request */
                a_last_contig = a0 = ai;
                buf = reqs[i].xbuf;
            }
            disps[j] = MPI_Aint_diff(ai, a0);

            req_size = blocklens[last_contig_req];
            req_size += blocklens[j];
#ifdef HAVE_MPI_LARGE_COUNT
            if (MPI_Aint_diff(ai, a_last_contig) == blocklens[last_contig_req]) {
                /* user buffer of request j is contiguous from j-1
                 * we coalesce j to j-1 */
                blocklens[last_contig_req] += blocklens[j];
            }
#else
            /* if req_size overflows 4-byte int, then skip coalescing */
            if (req_size <= NC_MAX_INT &&
                MPI_Aint_diff(ai,- a_last_contig) == blocklens[last_contig_req]) {
                /* user buffer of request j is contiguous from j-1
                 * we coalesce j to j-1 */
                blocklens[last_contig_req] += blocklens[j];
            }
#endif
            else if (j > 0) {
                /* not contiguous from request last_contig_req */
                last_contig_req++;
                a_last_contig = ai;
                disps[last_contig_req] = MPI_Aint_diff(ai, a0);
                blocklens[last_contig_req] = blocklens[i];
            }
            j++;
        }

        /* last_contig_req is the index of last contiguous request */
        if (last_contig_req == 0) {
            /* user buffers can be concatenated into a contiguous buffer */
            buf_type = MPI_BYTE;
            buf_count = blocklens[0];
        }
        else {
            /* after possible concatenating the user buffers, the true number
             * of non-contiguous buffers is last_contig_req+1 */
            int num_contig_reqs = last_contig_req+1;

            /* concatenate buffer addresses into a single buffer type */
#ifdef HAVE_MPI_LARGE_COUNT
            mpireturn = MPI_Type_create_hindexed_c(num_contig_reqs, blocklens,
                                                   disps, MPI_BYTE, &buf_type);
#else
            mpireturn = MPI_Type_create_hindexed(num_contig_reqs, blocklens,
                                                 disps, MPI_BYTE, &buf_type);
#endif
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

            buf_count = 1;
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

    fh = ncp->independent_fh;
    if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
        fh = ncp->collective_fh;

    /* set the MPI-IO fileview, this is a collective call */
    err = ncmpio_file_set_view(ncp, fh, &offset, filetype);
    if (filetype != MPI_BYTE) MPI_Type_free(&filetype);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        if (coll_indep == NC_REQ_INDEP) return status;
        buf_count = 0;
    }

    /* call MPI_File_read/MPI_File_write */
    err = ncmpio_read_write(ncp, rw_flag, coll_indep, offset, buf_count,
                            buf_type, buf, ((buf_type == MPI_BYTE) ? 1 : 0));
    if (status == NC_NOERR) status = err;

    if (buf_type != MPI_BYTE) MPI_Type_free(&buf_type);

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                 MPI_INFO_NULL);
     */

    return status;
}
