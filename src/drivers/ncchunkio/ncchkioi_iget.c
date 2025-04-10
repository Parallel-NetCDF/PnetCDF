/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncchkio_driver.h>
#include "ncchkio_internal.h"

static inline int
ncchkioi_init_get_req( NC_chk *ncchkp,
                        NC_chk_req *req,
                        int        varid,
                        const MPI_Offset *start,
                        const MPI_Offset *count,
                        const MPI_Offset *stride, 
                        const MPI_Offset  *imap,
                        void              *buf,
                        MPI_Offset        bufcount,
                        MPI_Datatype      buftype) {
    int err;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int overlapsize, packoff;
    MPI_Datatype ptype; // Pack datatype
    NC_chk_var *varp = ncchkp->vars.data + varid;

    // Zero out the request
    memset(req, 0, sizeof(NC_chk_req));

    // Record request
    req->starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*));
    req->start = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
    req->starts[0] = req->start;
    memcpy(req->start, start, sizeof(MPI_Offset) * varp->ndim);
    req->counts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*));
    req->count = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
    req->counts[0] = req->count;
    memcpy(req->count, count, sizeof(MPI_Offset) * varp->ndim);
    if (stride != NULL){
        req->stride = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
        memcpy(req->stride, stride, sizeof(MPI_Offset) * varp->ndim);
    }

    req->varid = varid;
    req->buf = (void*)buf;
    req->nreq = 1;
    req->buftype = buftype;
    if (varp->etype != buftype){
        if (bufcount > 0){
            req->bufcount = bufcount;
        }
        else{
            int i;

            req->bufcount = 1;
            for(i = 0; i < varp->ndim; i++){
                req->bufcount *= count[i];
            }
        }

        req->xbuf = (char*)NCI_Malloc(req->bufcount * varp->esize);
    }
    else{
        req->xbuf = req->buf;
    }

    req->xbufs = (char**)NCI_Malloc(sizeof(char*));
    req->xbufs[0] = req->xbuf;

    return NC_NOERR;
}

int
ncchkioi_iget_var(NC_chk        *ncchkp,
              int               varid,
              const MPI_Offset        *start,
              const MPI_Offset        *count,
              const MPI_Offset        *stride,
              const MPI_Offset  *imap,
              void              *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               *reqid)
{
    int err;
    int req_id;
    NC_chk_req req;

    // Init request
    err = ncchkioi_init_get_req(ncchkp, &req, varid, start, count, stride, imap, buf, bufcount, buftype);

    // Add to req list
    ncchkioi_req_list_add(&(ncchkp->getlist), &req_id);
    ncchkp->getlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id * 2;
    }

    return NC_NOERR;
}

static inline int
ncchkioi_init_get_varn_req( NC_chk *ncchkp,
                        NC_chk_req *req,
                        int        varid,
                        int        nreq,
                        MPI_Offset *const*starts,
                        MPI_Offset *const*counts, 
                        void              *buf,
                        MPI_Offset        bufcount,
                        MPI_Datatype      buftype) {
    int i, j;
    MPI_Offset rsize, boff;
    NC_chk_var *varp = ncchkp->vars.data + varid;

    // Zero out the request
    memset(req, 0, sizeof(NC_chk_req));

    // Record request
    req->starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * nreq);
    req->start = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim * nreq);
    for(i = 0; i < nreq; i++){
        req->starts[i] = req->start + i * varp->ndim;
        memcpy(req->starts[i], starts[i], sizeof(MPI_Offset) * varp->ndim);
    }
    req->counts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * nreq);
    req->count = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim * nreq);
    for(i = 0; i < nreq; i++){
        req->counts[i] = req->count + i * varp->ndim;
        memcpy(req->counts[i], counts[i], sizeof(MPI_Offset) * varp->ndim);
    }
    
    req->varid = varid;
    req->buf = (void*)buf;
    req->xbuf = (void*)buf;
    req->nreq = nreq;
    req->buftype = buftype;
    if (varp->etype != buftype){
        if (bufcount > 0){
            req->bufcount = bufcount;
        }
        else{
            req->bufcount = 0;
            for(i = 0; i < nreq; i++){
                rsize = 1;
                for(j = 0; j < varp->ndim; j++){
                    rsize *= counts[i][j];
                }
                req->bufcount += rsize;
            }
        }

        req->xbuf = (char*)NCI_Malloc(req->bufcount * varp->esize);
    }
    else{
        req->xbuf = req->buf;
    }

    // Calculate buffer for each individual request
    req->xbufs = (char**)NCI_Malloc(sizeof(char*) * nreq);
    boff = 0;
    for(i = 0; i < nreq; i++){
        req->xbufs[i] = (req->xbuf + boff);

        // Advance pointer by size of the request
        rsize = varp->esize;
        for(j = 0; j < varp->ndim; j++){
            rsize *= counts[i][j];
        }
        boff += rsize;
    }

    return NC_NOERR;
}

int
ncchkioi_iget_varn(NC_chk        *ncchkp,
              int               varid,
              int               nreq,
              MPI_Offset * const*starts,
              MPI_Offset * const*counts,
              void              *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               *reqid)
{
    int err;
    int req_id;
    NC_chk_req req;

    if (nreq > 1){
        err = ncchkioi_init_get_varn_req(ncchkp, &req, varid, nreq, starts, counts, buf, bufcount, buftype);
    }
    else{
        err = ncchkioi_init_get_req(ncchkp, &req, varid, starts[0], counts[0], NULL, NULL, buf, bufcount, buftype);
    }

    // Add to req list
    ncchkioi_req_list_add(&(ncchkp->getlist), &req_id);
    ncchkp->getlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id * 2;
    }

    return NC_NOERR;
}
