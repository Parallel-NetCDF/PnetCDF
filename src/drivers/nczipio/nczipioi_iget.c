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
#include <nczipio_driver.h>
#include "nczipio_internal.h"

static inline int
nczipioi_init_get_req( NC_zip *nczipp,
                        NC_zip_req *req,
                        int        varid,
                        MPI_Offset *start,
                        MPI_Offset *count,
                        MPI_Offset *stride, 
                        const MPI_Offset  *imap,
                        void              *buf,
                        MPI_Offset        bufcount,
                        MPI_Datatype      buftype) {
    int err;
    int i, j, k, l;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int *cstart, *cend, *citr; // Bounding box for chunks overlapping my own write region
    int overlapsize, packoff;
    MPI_Datatype ptype; // Pack datatype
    NC_zip_var *varp = nczipp->vars.data + varid;

    // Zero out the request
    memset(req, 0, sizeof(NC_zip_req));

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
    req->xbuf = (void*)buf;
    req->nreq = 1;

    return NC_NOERR;
}

int
nczipioi_iget_var(NC_zip        *nczipp,
              int               varid,
              MPI_Offset        *start,
              MPI_Offset        *count,
              MPI_Offset        *stride,
              const MPI_Offset  *imap,
              void              *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               *reqid)
{
    int err;
    int req_id;
    NC_zip_req req;

    // Init request
    err = nczipioi_init_get_req(nczipp, &req, varid, start, count, stride, imap, buf, bufcount, buftype);

    // Add to req list
    nczipioi_req_list_add(&(nczipp->getlist), &req_id);
    nczipp->getlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id;
    }

    return NC_NOERR;
}

static inline int
nczipioi_init_get_varn_req( NC_zip *nczipp,
                        NC_zip_req *req,
                        int        varid,
                        int        nreq,
                        MPI_Offset *const*starts,
                        MPI_Offset *const*counts, 
                        void              *buf,
                        MPI_Offset        bufcount,
                        MPI_Datatype      buftype) {
    int i, j;
    MPI_Offset rsize, boff;
    NC_zip_var *varp = nczipp->vars.data + varid;

    // Zero out the request
    memset(req, 0, sizeof(NC_zip_req));

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

    // Calculate buffer for each individual request
    req->xbufs = (char**)NCI_Malloc(sizeof(char*) * nreq);
    boff = 0;
    for(i = 0; i < nreq; i++){
        req->xbufs[i] = (((char*)buf) + boff);

        // Advance pointer by size of the request
        rsize = varp->esize;
        for(j = 0; j < varp->ndim; j++){
            rsize *= counts[i][j];
        }
        boff += rsize;
    }

    req->varid = varid;
    req->buf = (void*)buf;
    req->xbuf = (void*)buf;
    req->nreq = nreq;

    return NC_NOERR;
}

int
nczipioi_iget_varn(NC_zip        *nczipp,
              int               varid,
              int               nreq,
              MPI_Offset        **starts,
              MPI_Offset        **counts,
              void              *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               *reqid)
{
    int err;
    int req_id;
    NC_zip_req req;

    if (nreq > 1){
        err = nczipioi_init_get_varn_req(nczipp, &req, varid, nreq, starts, counts, buf, bufcount, buftype);
    }
    else{
        err = nczipioi_init_get_req(nczipp, &req, varid, starts[0], counts[0], NULL, NULL, buf, bufcount, buftype);
    }

    // Add to req list
    nczipioi_req_list_add(&(nczipp->getlist), &req_id);
    nczipp->getlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id;
    }

    return NC_NOERR;
}
