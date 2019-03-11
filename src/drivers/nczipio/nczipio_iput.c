/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
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

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

static int get_chunk_idx(NC_zip_var *varp, int* cord){
    int i, ret;
    
    ret = cord[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->chunkdim[i - 1] + cord[i];
    }

    return ret;
}

static int get_chunk_cord(NC_zip_var *varp, int idx, int* cord){
    int i, ret;
    
    ret = cord[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->chunkdim[i - 1] + cord[i];
    }

    for(i = varp->ndim - 1; i > 0; i--){
        cord[i] = idx % varp->chunkdim[i - 1];
        idx /= varp->chunkdim[i - 1];
    }
    cord[0] = idx;

    return 0;
}

static int get_chunk_overlap(NC_zip_var *varp, int* cord, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, int *ostart, int *ocount){
    int i, ret;

    for(i = 0; i < varp->ndim; i++){
        ostart[i] = max(start[i], cord[i] * varp->chunkdim[i]);
        ocount[i] = (min(start[i] + count[i] * stride[i], (cord[i] + 1) * varp->chunkdim[i]) - ostart[i]) / stride[i];
        if (ocount[i] < 0){
            ocount[i] = 0;
        }
    }

    return 0;
}

int
nczipioi_init_put_req_by_chunk( NC_zip     *nczipp,
                            NC_zip_req *req,
                            NC_zip_var *varp,
                            MPI_Offset *starts,
                            MPI_Offset *counts,
                            MPI_Offset *stride 
                            const void *buf) {
    int err;
    int i, j, k, l;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int *cstart, *cend, *citr; // Bounding box for chunks overlapping my own write region
    int overlapsize, packoff;
    MPI_Datatype ptype; // Pack datatype

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    // Starting, ending, current chunk position
    cstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    citr = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    cend = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    // Record request
    req.start = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
    memcpy(req.start, start, sizeof(MPI_Offset) * varp->ndim);
    req.count = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
    memcpy(req.count, count, sizeof(MPI_Offset) * varp->ndim);
    if (stride != NULL){
        req.stride = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
        memcpy(req.stride, stride, sizeof(MPI_Offset) * varp->ndim);
    }

    req.nsend = nczipioi_chunk_itr_init(varp, start, count, stride, c start, cend, citr);

    req.widx = (int*)NCI_Malloc(sizeof(int) * req.nsend);
    req.sbuf = (char**)NCI_Malloc(sizeof(char*) * req.nsend);
    req.sreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * req.nsend);
    req.sstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * req.nsend);

    req.buf = buf;
    req.xbuf = xbuf;
    req.nreq = 1;

    // Iterate through chunk
    req.nsend = 0;  // Previous estimate contains our own chunks. Now, we count real chunk
    do{
        // Chunk index
        i = get_chunk_idx(varp, citr);

        if (varp->chunk_owner[i] != nczipp->rank){
            // Overlapping size
            get_chunk_overlap(varp, citr, start, count, tstart, tssize);
            overlapsize = varp->esize;
            for(j = 0; j < varp->ndim; j++){
                overlapsize *= tssize[j];                     
            }
            printf("overlapsize = %d\n", overlapsize); fflush(stdout);
            
            // Pack type
            for(j = 0; j < varp->ndim; j++){
                tstart[j] -= start[j];
                tsize[j] = (int)count[j];
            }
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
            MPI_Type_commit(&ptype);
            
            // Allocate buffer
            req.sbuf[req.nsend] = (char*)NCI_Malloc(overlapsize + sizeof(int) * (varp->ndim * 2 + 1));

            // Pack data
            *((int*)req.sbuf[req.nsend]) = 1;
            packoff = sizeof(int);
            MPI_Pack(starts[i], varp->ndim, MPI_INT, req.sbuf[req.nsend], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
            MPI_Pack(counts[i], varp->ndim, MPI_INT, req.sbuf[req.nsend], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
            MPI_Pack(bufs[i], 1, ptype, req.sbuf[req.nsend], packoff + overlapsize, &packoff, nczipp->comm);

            // Free packtype
            MPI_Type_free(&ptype);

            // Send data to owner
            MPI_Isend(req.sbuf[req.nsend], packoff, MPI_BYTE, varp->chunk_owner[i], i, nczipp->comm, req.sreqs + req.nsend);

            // Record chunk write
            req.widx[req.nsend] = i;
            
            //Count send request
            req.nsend++;
        }
    // Move to next chunk        
    } while (nczipioi_chunk_itr_next(varp, start, count, stride, c start, cend, citr));

    // Free buffers
    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);

    NCI_Free(cstart);
    NCI_Free(ccord);
    NCI_Free(cend);

    return NC_NOERR;
}

int
nczipioi_iput_var(NC_zip        *nczipp,
              NC_zip_var        *varp,
              MPI_Offset        *starts,
              MPI_Offset        *counts,
              const void        *xbuf,
              const void        *buf,
              int               *reqid)
{
    int err;
    int req_id;
    NC_zip_req req;

    err = nczipioi_init_put_req(nczipp, &req, varp, start, count, stride, xbuf, buf);

    // Release var info
    adios_free_varinfo (v);

    // Add to req list
    nczipioi_list_add(&(nczipp->putlist), &req_id);
    ncadp->putlist.reqs[req_id] = req;
    
    if (reqid != NULL){
        *reqid = req_id;
    }

    return NC_NOERR;
}

int
nczipioi_init_put_varn_req( NC_zip *nczipp,
                        NC_zip_req *req,
                        NC_zip_var *varp,
                        int        nreq,
                        MPI_Offset *start,
                        MPI_Offset *count,
                        MPI_Offset *stride, 
                        const void *xbuf,
                        const void *buf) {
    int err;
    int i, j, k, l;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int *cstart, *cend, *citr; // Bounding box for chunks overlapping my own write region
    int *wcnt_local, *wcnt_all;   // Number of processes that writes to each chunk
    char *sbuf_base, *sbuf_cur; // Send buffer, exactly the same size as buf
    int put_size_total, put_size;  // Total size of buf and size of data of a single req
    char **bufs; // Location of data in buf for every req
    int overlapsize, packoff;
    MPI_Datatype ptype; // Pack datatype
    MPI_Request *sreqs, *rreqs;    // Send and recv req
    char **sbufs, **rbufs;   // Send and recv buffer
    MPI_Status *sstats, *rstats;    // Send and recv status
    int nsend, nrecv;
    int *rsizes;
    MPI_Message rmsg;   // Receive message
    int *zsizes, *zsizes_all, *zoffs;
    MPI_Offset zstart, zcount;
    char **zbufs;
    int zdimid;
    char name[128]; // Name of objects

    // Allocate buffering for location of data in buf for every req
    req.bufs = (char**)NCI_Malloc(sizeof(char*) * nreq);

    // Allocate buffering for write count
    wcnt_local = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    wcnt_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    // Starting, ending, current chunk position
    cstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    citr = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    cend = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    // Record request
    req.starts = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset*) * nreq);
    req.starts[0] = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
    memcpy(req.start, start, sizeof(MPI_Offset) * varp->ndim);
    req.counts = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
    memcpy(req.count, count, sizeof(MPI_Offset) * varp->ndim);
    req.nsend = nczipioi_chunk_itr_init(varp, start, count, stride, c start, cend, citr);

    req.widx = (int*)NCI_Malloc(sizeof(int) * req.nsend);
    req.sbuf = (char**)NCI_Malloc(sizeof(char*) * req.nsend);
    req.sreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * req.nsend);
    req.sstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * req.nsend);

    req.buf = buf;
    req.xbuf = xbuf;
    req.nreq = nreq;

    //Calculate local write count, we calculate offset and size of each req by the way
    memset(wcnt_local, 0, sizeof(int) * nczipp->np);

    // Iterate through chunk
    req.nsend = 0;  // Previous estimate contains our own chunks. Now, we count real chunk
    do{
        // Chunk index
        i = get_chunk_idx(varp, citr);

        if (varp->chunk_owner[i] != nczipp->rank){
            // Overlapping size
            get_chunk_overlap(varp, citr, start, count, tstart, tssize);
            overlapsize = varp->esize;
            for(j = 0; j < varp->ndim; j++){
                overlapsize *= tssize[j];                     
            }
            printf("overlapsize = %d\n", overlapsize); fflush(stdout);
            
            // Pack type
            for(j = 0; j < varp->ndim; j++){
                tstart[j] -= start[j];
                tsize[j] = (int)count[j];
            }
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
            MPI_Type_commit(&ptype);
            
            // Pack data
            MPI_Pack(starts[i], varp->ndim, MPI_INT, sbuf_cur, packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
            MPI_Pack(counts[i], varp->ndim, MPI_INT, sbuf_cur, packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
            MPI_Pack(bufs[i], 1, ptype, sbuf_cur, packoff + overlapsize, &packoff, nczipp->comm);

            // Free packtype
            MPI_Type_free(&ptype);
            req.widx[req.nsend] = i;

        }
        
    } while (nczipioi_chunk_itr_next(varp, start, count, stride, c start, cend, citr));

    // All reduce
    MPI_Allreduce(zsizes, zsizes_all, varp->nchunks, MPI_INT, MPI_MAX, nczipp->comm);
    zoffs[0] = 0;
    for(i = 1; i < varp->nchunks; i++){
        zoffs[i] = zoffs[i - 1] + zsizes_all[i - 1];
    }

    // Enter redefine mode
    nczipp->driver->redef(nczipp->ncp);

    // Define dimension  for data variable
    sprintf(name, "_compressed_data_dim_%d", varp->varid);
    err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[varp->nchunks - 1] + zsizes_all[varp->nchunks - 1], &zdimid);
    if (err != NC_NOERR) return err;

    // Define variable
    sprintf(name, "_compressed_data_%d", varp->varid);
    err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));
    if (err != NC_NOERR) return err;

    // Record offset of chunks in data variable
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkoffset", NC_INT, varp->nchunks, zoffs, MPI_INT); // Original datatype
    if (err != NC_NOERR) return err;

    // Record size of chunks
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunklen", NC_INT, varp->nchunks, zsizes_all, MPI_INT); // Original datatype
    if (err != NC_NOERR) return err;

    // Switch to data mode
    err = nczipp->driver->enddef(nczipp->ncp);
    if (err != NC_NOERR) return err;

    // Do I/O
    for(l = 0; l < varp->nmychunk; l++){
        k = varp->mychunks[l];
        zstart = (MPI_Offset)zoffs[k];
        zcount = (MPI_Offset)zsizes[k];
        printf("cache[0] = %d, cache[1] = %d, off = %lld, cnt = %lld\n", ((int*)(varp->chunk_cache[k]))[0], ((int*)(varp->chunk_cache[k]))[1], zstart, zcount); fflush(stdout);
        nczipp->driver->iput_var(nczipp->ncp, varp->datavarid, &zstart, &zcount, NULL, NULL, zbufs[l], (MPI_Offset)(zsizes_all[k]), MPI_UNSIGNED_CHAR, NULL, NC_REQ_WR | NC_REQ_NBI | NC_REQ_FLEX);
    }
    nczipp->driver->wait(nczipp->ncp, NC_REQ_ALL, NULL, NULL, NC_REQ_COLL);

    // Free buffers

    NCI_Free(bufs);

    NCI_Free(wcnt_local);
    NCI_Free(wcnt_all);

    NCI_Free(zsizes);
    NCI_Free(zsizes_all);
    NCI_Free(zoffs);
    for(l = 0; l < varp->nmychunk; l++){
        NCI_Free(zbufs[l]);
    }
    NCI_Free(zbufs);

    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);

    NCI_Free(cstart);
    NCI_Free(ccord);
    NCI_Free(cend);

    NCI_Free(sbuf_base);

    NCI_Free(sreqs);
    NCI_Free(sstats);

    NCI_Free(rreqs);
    NCI_Free(rstats);
    for(i = 0; i < nrecv; i++){
        NCI_Free(rbufs[i]);
    }
    NCI_Free(rbufs);
    NCI_Free(rsizes);

    return NC_NOERR;
}

int
nczipioi_iput_varn(NC_zip        *nczipp,
              NC_zip_var        *varp,
              int               nreq,
              MPI_Offset        *starts,
              MPI_Offset        *counts,
              const void        *xbuf,
              const void        *buf,
              int               *reqid)
{
    int err;
    int i, j, k, l;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int *cstart, *cend, *ccord; // Bounding box for chunks overlapping my own write region
    int *wcnt_local, *wcnt_all;   // Number of processes that writes to each chunk
    char *sbuf_base, *sbuf_cur; // Send buffer, exactly the same size as buf
    int put_size_total, put_size;  // Total size of buf and size of data of a single req
    char **bufs; // Location of data in buf for every req
    int overlapsize, packoff;
    MPI_Datatype ptype; // Pack datatype
    MPI_Request *sreqs, *rreqs;    // Send and recv req
    char **sbufs, **rbufs;   // Send and recv buffer
    MPI_Status *sstats, *rstats;    // Send and recv status
    int nsend, nrecv;
    int *rsizes;
    MPI_Message rmsg;   // Receive message
    int *zsizes, *zsizes_all, *zoffs;
    MPI_Offset zstart, zcount;
    char **zbufs;
    int zdimid;
    char name[128]; // Name of objects

    // Allocate buffering for location of data in buf for every req
    bufs = (char**)NCI_Malloc(sizeof(char*) * nreq);

    // Allocate buffering for write count
    wcnt_local = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    wcnt_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    // Starting, ending, current chunk position
    cstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    ccord = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    ccord_raw = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    cend = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    //Calculate local write count, we calculate offset and size of each req by the way
    memset(wcnt_local, 0, sizeof(int) * nczipp->np);
    
    // Chunk boundary
    req.nsend = 1;
    for(i = 0; i < varp->ndim; i++){
        cstart[i] = starts[i] / varp->chunkdim[i];
        cend[i] = (starts[i] + (counts[i] - 1) * stride[i]) / varp->chunkdim[i] + 1;
        req.nsend *= cend[i] - cstart[i] + 1;
    }

    err = nczipioi_init_put_varn_req(nczipp, &req, varp, nreq, start, count, stride, xbuf, buf);

    // Release var info
    adios_free_varinfo (v);

        // move on to next chunk           
        ccord[varp->ndim - 1]++;
        for(j = varp->ndim - 1; j > 0; j--){
            if (ccord[j] >= cend[j]){
                ccord[j - 1]++;
                ccord[j] = cstart[j];
            }
            else{
                break;
            }
        }
    }
    
    // Allocate send buffer
    // TODO: more efficient estimation
    sbuf_base = (char*)NCI_Malloc(put_size + varp->ndim * 2 * nreq * varp->nchunks);
    sbuf_cur = sbuf_base;

    // Calculate number of send request
    nsend = 0;
    for(i = 0; i < varp->nchunks; i++){
        if (i != nczipp->rank){
            nsend += wcnt_local[i];
        }
    }
    sreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nsend);
    sstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nsend);

    // Sync write count
    MPI_Allreduce(wcnt_local, wcnt_all, nczipp->np, MPI_INT, MPI_SUM, nczipp->comm);

    // Calculate number of recv request
    nrecv = 0;
    for(i = 0; i < varp->nchunks; i++){
        if (varp->chunk_owner[i] == nczipp->rank){
            nrecv += wcnt_all[i] - wcnt_local[i];
        }
    }
    rreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nrecv);
    rstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nrecv);
    rbufs = (char**)NCI_Malloc(sizeof(char*) * nrecv);
    rsizes = (int*)NCI_Malloc(sizeof(int) * nrecv);

    // Post send and recv
    nrecv = 0;
    nsend = 0;
    memcpy(tsize, varp->chunkdim, sizeof(int) * varp->ndim);
    for(k = 0; k < varp->nchunks; k++){
        if (varp->chunk_owner[k] == nczipp->rank){
            // Receive data from other process
            for(i = 0; i < wcnt_all[k] - wcnt_local[k]; i++){
                // Get message size, including metadata
                MPI_Mprobe(MPI_ANY_SOURCE, k, nczipp->comm, &rmsg, rstats);
                MPI_Get_count(rstats, MPI_BYTE, rsizes + nrecv);

                printf("rsize = %d\n", rsizes[i]); fflush(stdout);

                // Allocate buffer
                rbufs[nrecv] = (char*)NCI_Malloc(rsizes[nrecv]);

                // Post irecv
                MPI_Imrecv(rbufs[nrecv], rsizes[nrecv], MPI_BYTE, &rmsg, rreqs + nrecv);
                nrecv++;
            }
        }
        else{
            if (wcnt_local[k] > 0){
                packoff = 0;
                for(i = 0; i < nreq; i++){
                    get_chunk_cord(varp, k, citr);                    
                    get_chunk_overlap(varp, citr, starts[i], counts[i], tstart, tssize);
                    printf("cord = %d, start = %lld, count = %lld, tstart = %d, tssize = %d, esize = %d, ndim = %d\n", citr[0], starts[i][0], counts[i][0], tstart[0], tssize[0], esize, varp->ndim); fflush(stdout);
                    overlapsize = esize;
                    for(j = 0; j < varp->ndim; j++){
                        overlapsize *= tssize[j];                     
                    }
                    printf("overlapsize = %d\n", overlapsize); fflush(stdout);
                    if (overlapsize > 0){
                        // Pack type
                        for(j = 0; j < varp->ndim; j++){
                            tstart[j] -= starts[i][j];
                            tsize[j] = (int)counts[i][j];
                        }
                        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                        MPI_Type_commit(&ptype);
                        
                        // Pack data
                        MPI_Pack(starts[i], varp->ndim, MPI_INT, sbuf_cur, packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        MPI_Pack(counts[i], varp->ndim, MPI_INT, sbuf_cur, packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        MPI_Pack(bufs[i], 1, ptype, sbuf_cur, packoff + overlapsize, &packoff, nczipp->comm);

                        MPI_Type_free(&ptype);
                    }
                }

                printf("packoff = %d\n", packoff); fflush(stdout);
                MPI_Isend(sbuf_cur, packoff, MPI_BYTE, varp->chunk_owner[k], k, nczipp->comm, sreqs + nsend);
                sbuf_cur += packoff;
                nsend++;
            }
        }
    }

    // Wait for all send
    MPI_Waitall(nsend, sreqs, sstats);

    // handle each chunk we own
    nrecv = 0;
    memset(zsizes, 0, sizeof(int) * varp->nchunks);
    for(l = 0; l < varp->nmychunk; l++){
        k = varp->mychunks[l];
        if (varp->chunk_owner[k] == nczipp->rank){
            // TODO: bring chunk into cache from disk if needed
            if (varp->chunk_cache[k] == NULL){
                varp->chunk_cache[k] = (char*)NCI_Malloc(varp->chunksize);
            }

            // Handle our own data
            if (wcnt_local[k] > 0){
                for(i = 0; i < nreq; i++){
                    get_chunk_cord(varp, k, citr);
                    get_chunk_overlap(varp, citr, starts[i], counts[i], tstart, tssize);

                    overlapsize = esize;
                    for(j = 0; j < varp->ndim; j++){
                        overlapsize *= tssize[j];
                    }
                    
                    if (overlapsize > 0){
                        // Pack into contiguous buffer
                        // Pack type
                        for(j = 0; j < varp->ndim; j++){
                            tstart[j] -= citr[j] * varp->chunkdim[j];
                            tsize[j] = (int)counts[i][j];
                        }
                        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                        MPI_Type_commit(&ptype);

                        packoff = 0;
                        MPI_Pack(bufs[i], 1, ptype, sbuf_cur, overlapsize, &packoff, nczipp->comm);

                        MPI_Type_free(&ptype);

                        // Unpack into chunk cache
                        // Pack type
                        for(j = 0; j < varp->ndim; j++){
                            tsize[j] = (int)varp->chunkdim[j];
                        }
                        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                        MPI_Type_commit(&ptype);
                        
                        // Pack data
                        packoff = 0;
                        MPI_Unpack(sbuf_cur, overlapsize, &packoff, varp->chunk_cache[k], 1, ptype, nczipp->comm);

                        MPI_Type_free(&ptype);
                    }
                }
            }

            // Receive all data from other processes
            
            // Wait for all send requests for the chunk
            MPI_Waitall(wcnt_all[k] - wcnt_local[k], rreqs + nrecv, rstats + nrecv);

            // Process data received
            printf("nrecv = %d, wcnt_all = %d, wcnt_local = %d\n", nrecv, wcnt_all[k], wcnt_local[k]); fflush(stdout);
            for(i = nrecv; i < nrecv + wcnt_all[k] - wcnt_local[k]; i++){
                packoff = 0;
                printf("rsize_2 = %d\n", rsizes[i]); fflush(stdout);
                while(packoff < rsizes[i]){
                    MPI_Unpack(rbufs[i], rsizes[i], &packoff, tstart, varp->ndim, MPI_INT, nczipp->comm);
                    MPI_Unpack(rbufs[i], rsizes[i], &packoff, tssize, varp->ndim, MPI_INT, nczipp->comm);

                    for(j = 0; j < varp->ndim; j++){
                       tsize[j] = (int)varp->chunkdim[j];
                    }
                    
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                    MPI_Type_commit(&ptype);

                    printf("tsize = %d, tssize = %d, tstart = %d, buf = %d\n", tsize[0], tssize[0], tstart[0], *((int*)(rbufs[i] + packoff))); fflush(stdout);
                    MPI_Unpack(rbufs[i], rsizes[i], &packoff, varp->chunk_cache[k], 1, ptype, nczipp->comm);
                    printf("cache[0] = %d, cache[1] = %d\n", ((int*)(varp->chunk_cache[k]))[0], ((int*)(varp->chunk_cache[k]))[1]); fflush(stdout);
                    MPI_Type_free(&ptype);
                }
            }
            nrecv += wcnt_all[k] - wcnt_local[k];

            // Apply compression

            // Test comrpessed size
            nczipp->zip->compress(varp->chunk_cache[k], varp->chunksize, NULL, zsizes + l, varp->ndim, varp->chunkdim, etype);

            // Buffer for comrpessed data
            zbufs[l] = (char*)NCI_Malloc(zsizes[l]);

            // Perform real compression
            nczipp->zip->compress(varp->chunk_cache[k], varp->chunksize, zbufs[l], zsizes + l, varp->ndim, varp->chunkdim, etype);
        }
    }

    // All reduce
    MPI_Allreduce(zsizes, zsizes_all, varp->nchunks, MPI_INT, MPI_MAX, nczipp->comm);
    zoffs[0] = 0;
    for(i = 1; i < varp->nchunks; i++){
        zoffs[i] = zoffs[i - 1] + zsizes_all[i - 1];
    }

    // Enter redefine mode
    nczipp->driver->redef(nczipp->ncp);

    // Define dimension  for data variable
    sprintf(name, "_compressed_data_dim_%d", varp->varid);
    err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[varp->nchunks - 1] + zsizes_all[varp->nchunks - 1], &zdimid);
    if (err != NC_NOERR) return err;

    // Define variable
    sprintf(name, "_compressed_data_%d", varp->varid);
    err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));
    if (err != NC_NOERR) return err;

    // Record offset of chunks in data variable
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkoffset", NC_INT, varp->nchunks, zoffs, MPI_INT); // Original datatype
    if (err != NC_NOERR) return err;

    // Record size of chunks
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunklen", NC_INT, varp->nchunks, zsizes_all, MPI_INT); // Original datatype
    if (err != NC_NOERR) return err;

    // Switch to data mode
    err = nczipp->driver->enddef(nczipp->ncp);
    if (err != NC_NOERR) return err;

    // Do I/O
    for(l = 0; l < varp->nmychunk; l++){
        k = varp->mychunks[l];
        zstart = (MPI_Offset)zoffs[k];
        zcount = (MPI_Offset)zsizes[k];
        printf("cache[0] = %d, cache[1] = %d, off = %lld, cnt = %lld\n", ((int*)(varp->chunk_cache[k]))[0], ((int*)(varp->chunk_cache[k]))[1], zstart, zcount); fflush(stdout);
        nczipp->driver->iput_var(nczipp->ncp, varp->datavarid, &zstart, &zcount, NULL, NULL, zbufs[l], (MPI_Offset)(zsizes_all[k]), MPI_UNSIGNED_CHAR, NULL, NC_REQ_WR | NC_REQ_NBI | NC_REQ_FLEX);
    }
    nczipp->driver->wait(nczipp->ncp, NC_REQ_ALL, NULL, NULL, NC_REQ_COLL);

    // Free buffers

    NCI_Free(bufs);

    NCI_Free(wcnt_local);
    NCI_Free(wcnt_all);

    NCI_Free(zsizes);
    NCI_Free(zsizes_all);
    NCI_Free(zoffs);
    for(l = 0; l < varp->nmychunk; l++){
        NCI_Free(zbufs[l]);
    }
    NCI_Free(zbufs);

    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);

    NCI_Free(cstart);
    NCI_Free(citr);
    NCI_Free(cend);

    NCI_Free(sbuf_base);

    NCI_Free(sreqs);
    NCI_Free(sstats);

    NCI_Free(rreqs);
    NCI_Free(rstats);
    for(i = 0; i < nrecv; i++){
        NCI_Free(rbufs[i]);
    }
    NCI_Free(rbufs);
    NCI_Free(rsizes);

    return NC_NOERR;
}


int
nczipioi_iput_var(NC_zip        *nczipp,
              NC_zip_var       *varp,
              MPI_Offset* const *starts,
              MPI_Offset* const *counts,
              const void       *buf)
{
    int i, j, k, l, err;
    MPI_Datatype etype; // Variable element type in MPI
    int esize;  // Variable element size
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int *cstart, *cend, *citr; // Bounding box for chunks overlapping my own write region
    int *wcnt_local, *wcnt_all;   // Number of processes that writes to each chunk
    char *sbuf_base, *sbuf_cur; // Send buffer, exactly the same size as buf
    int put_size_total, put_size;  // Total size of buf and size of data of a single req
    char **bufs; // Location of data in buf for every req
    int overlapsize, packoff;
    MPI_Datatype ptype; // Pack datatype
    MPI_Request *sreqs, *rreqs;    // Send and recv req
    char **sbufs, **rbufs;   // Send and recv buffer
    MPI_Status *sstats, *rstats;    // Send and recv status
    int nsend, nrecv;
    int *rsizes;
    MPI_Message rmsg;   // Receive message
    int *zsizes, *zsizes_all, *zoffs;
    MPI_Offset zstart, zcount;
    char **zbufs;
    int zdimid;
    char name[128]; // Name of objects

    // Original datatype and size
    esize = NC_Type_size(varp->xtype);
    etype = ncmpii_nc2mpitype(varp->xtype);

    // Allocate buffering for location of data in buf for every req
    bufs = (char**)NCI_Malloc(sizeof(char*) * nreq);

    // Allocate buffering for write count
    wcnt_local = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    wcnt_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);

    // Allocate buffering for compression
    zsizes = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    zsizes_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    zoffs = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    zbufs = (char**)NCI_Malloc(sizeof(char*) * varp->nmychunk);

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    // Starting, ending, current chunk position
    cstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    citr = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    cend = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    //Calculate local write count, we caluculate offset and size of each req by the way
    memset(wcnt_local, 0, sizeof(int) * nczipp->np);
    put_size_total = 0;
    for(k = 0; k < nreq; k++){
        // Offset and size of each req
        bufs[k] = (char*)buf + put_size_total;
        put_size = esize;
        for(i = 0; i < varp->ndim; i++){
            put_size * counts[k][i];
        }
        put_size_total += put_size;

        // Chunk boundary
        for(i = 0; i < varp->ndim; i++){
            cstart[i] = starts[k][i] / varp->chunkdim[i];
            cend[i] = (starts[k][i] + counts[k][i] - 1) / varp->chunkdim[i] + 1;
        }

        // calculate local write count, at most one per chunk
        memcpy(citr, cstart, sizeof(int) * varp->ndim);
        while(citr[0] < cend[0]){
            j = get_chunk_idx(varp, citr);    
            wcnt_local[j] = 1;

            // move on to next chunk
            citr[varp->ndim - 1]++;
            for(j = varp->ndim - 1; j > 0; j--){
                if (citr[j] >= cend[j]){
                    citr[j - 1]++;
                    citr[j] = cstart[j];
                }
                else{
                    break;
                }
            }
        }
    }

    // Allocate send buffer
    // TODO: more efficient estimation
    sbuf_base = (char*)NCI_Malloc(put_size + varp->ndim * 2 * nreq * varp->nchunks);
    sbuf_cur = sbuf_base;

    // Calculate number of send request
    nsend = 0;
    for(i = 0; i < varp->nchunks; i++){
        if (i != nczipp->rank){
            nsend += wcnt_local[i];
        }
    }
    sreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nsend);
    sstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nsend);

    // Sync write count
    MPI_Allreduce(wcnt_local, wcnt_all, nczipp->np, MPI_INT, MPI_SUM, nczipp->comm);

    // Calculate number of recv request
    nrecv = 0;
    for(i = 0; i < varp->nchunks; i++){
        if (varp->chunk_owner[i] == nczipp->rank){
            nrecv += wcnt_all[i] - wcnt_local[i];
        }
    }
    rreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nrecv);
    rstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nrecv);
    rbufs = (char**)NCI_Malloc(sizeof(char*) * nrecv);
    rsizes = (int*)NCI_Malloc(sizeof(int) * nrecv);

    // Post send and recv
    nrecv = 0;
    nsend = 0;
    memcpy(tsize, varp->chunkdim, sizeof(int) * varp->ndim);
    for(k = 0; k < varp->nchunks; k++){
        if (varp->chunk_owner[k] == nczipp->rank){
            // Receive data from other process
            for(i = 0; i < wcnt_all[k] - wcnt_local[k]; i++){
                // Get message size, including metadata
                MPI_Mprobe(MPI_ANY_SOURCE, k, nczipp->comm, &rmsg, rstats);
                MPI_Get_count(rstats, MPI_BYTE, rsizes + nrecv);

                printf("rsize = %d\n", rsizes[i]); fflush(stdout);

                // Allocate buffer
                rbufs[nrecv] = (char*)NCI_Malloc(rsizes[nrecv]);

                // Post irecv
                MPI_Imrecv(rbufs[nrecv], rsizes[nrecv], MPI_BYTE, &rmsg, rreqs + nrecv);
                nrecv++;
            }
        }
        else{
            if (wcnt_local[k] > 0){
                packoff = 0;
                for(i = 0; i < nreq; i++){
                    get_chunk_cord(varp, k, citr);                    
                    get_chunk_overlap(varp, citr, starts[i], counts[i], tstart, tssize);
                    printf("cord = %d, start = %lld, count = %lld, tstart = %d, tssize = %d, esize = %d, ndim = %d\n", citr[0], starts[i][0], counts[i][0], tstart[0], tssize[0], esize, varp->ndim); fflush(stdout);
                    overlapsize = esize;
                    for(j = 0; j < varp->ndim; j++){
                        overlapsize *= tssize[j];                     
                    }
                    printf("overlapsize = %d\n", overlapsize); fflush(stdout);
                    if (overlapsize > 0){
                        // Pack type
                        for(j = 0; j < varp->ndim; j++){
                            tstart[j] -= starts[i][j];
                            tsize[j] = (int)counts[i][j];
                        }
                        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                        MPI_Type_commit(&ptype);
                        
                        // Pack data
                        MPI_Pack(starts[i], varp->ndim, MPI_INT, sbuf_cur, packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        MPI_Pack(counts[i], varp->ndim, MPI_INT, sbuf_cur, packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        MPI_Pack(bufs[i], 1, ptype, sbuf_cur, packoff + overlapsize, &packoff, nczipp->comm);

                        MPI_Type_free(&ptype);
                    }
                }

                printf("packoff = %d\n", packoff); fflush(stdout);
                MPI_Isend(sbuf_cur, packoff, MPI_BYTE, varp->chunk_owner[k], k, nczipp->comm, sreqs + nsend);
                sbuf_cur += packoff;
                nsend++;
            }
        }
    }

    // Wait for all send
    MPI_Waitall(nsend, sreqs, sstats);

    // handle each chunk we own
    nrecv = 0;
    memset(zsizes, 0, sizeof(int) * varp->nchunks);
    for(l = 0; l < varp->nmychunk; l++){
        k = varp->mychunks[l];
        if (varp->chunk_owner[k] == nczipp->rank){
            // TODO: bring chunk into cache from disk if needed
            if (varp->chunk_cache[k] == NULL){
                varp->chunk_cache[k] = (char*)NCI_Malloc(varp->chunksize);
            }

            // Handle our own data
            if (wcnt_local[k] > 0){
                for(i = 0; i < nreq; i++){
                    get_chunk_cord(varp, k, citr);
                    get_chunk_overlap(varp, citr, starts[i], counts[i], tstart, tssize);

                    overlapsize = esize;
                    for(j = 0; j < varp->ndim; j++){
                        overlapsize *= tssize[j];
                    }
                    
                    if (overlapsize > 0){
                        // Pack into contiguous buffer
                        // Pack type
                        for(j = 0; j < varp->ndim; j++){
                            tstart[j] -= citr[j] * varp->chunkdim[j];
                            tsize[j] = (int)counts[i][j];
                        }
                        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                        MPI_Type_commit(&ptype);

                        packoff = 0;
                        MPI_Pack(bufs[i], 1, ptype, sbuf_cur, overlapsize, &packoff, nczipp->comm);

                        MPI_Type_free(&ptype);

                        // Unpack into chunk cache
                        // Pack type
                        for(j = 0; j < varp->ndim; j++){
                            tsize[j] = (int)varp->chunkdim[j];
                        }
                        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                        MPI_Type_commit(&ptype);
                        
                        // Pack data
                        packoff = 0;
                        MPI_Unpack(sbuf_cur, overlapsize, &packoff, varp->chunk_cache[k], 1, ptype, nczipp->comm);

                        MPI_Type_free(&ptype);
                    }
                }
            }

            // Receive all data from other processes
            
            // Wait for all send requests for the chunk
            MPI_Waitall(wcnt_all[k] - wcnt_local[k], rreqs + nrecv, rstats + nrecv);

            // Process data received
            printf("nrecv = %d, wcnt_all = %d, wcnt_local = %d\n", nrecv, wcnt_all[k], wcnt_local[k]); fflush(stdout);
            for(i = nrecv; i < nrecv + wcnt_all[k] - wcnt_local[k]; i++){
                packoff = 0;
                printf("rsize_2 = %d\n", rsizes[i]); fflush(stdout);
                while(packoff < rsizes[i]){
                    MPI_Unpack(rbufs[i], rsizes[i], &packoff, tstart, varp->ndim, MPI_INT, nczipp->comm);
                    MPI_Unpack(rbufs[i], rsizes[i], &packoff, tssize, varp->ndim, MPI_INT, nczipp->comm);

                    for(j = 0; j < varp->ndim; j++){
                       tsize[j] = (int)varp->chunkdim[j];
                    }
                    
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                    MPI_Type_commit(&ptype);

                    printf("tsize = %d, tssize = %d, tstart = %d, buf = %d\n", tsize[0], tssize[0], tstart[0], *((int*)(rbufs[i] + packoff))); fflush(stdout);
                    MPI_Unpack(rbufs[i], rsizes[i], &packoff, varp->chunk_cache[k], 1, ptype, nczipp->comm);
                    printf("cache[0] = %d, cache[1] = %d\n", ((int*)(varp->chunk_cache[k]))[0], ((int*)(varp->chunk_cache[k]))[1]); fflush(stdout);
                    MPI_Type_free(&ptype);
                }
            }
            nrecv += wcnt_all[k] - wcnt_local[k];

            // Apply compression

            // Test comrpessed size
            nczipp->zip->compress(varp->chunk_cache[k], varp->chunksize, NULL, zsizes + l, varp->ndim, varp->chunkdim, etype);

            // Buffer for comrpessed data
            zbufs[l] = (char*)NCI_Malloc(zsizes[l]);

            // Perform real compression
            nczipp->zip->compress(varp->chunk_cache[k], varp->chunksize, zbufs[l], zsizes + l, varp->ndim, varp->chunkdim, etype);
        }
    }

    // All reduce
    MPI_Allreduce(zsizes, zsizes_all, varp->nchunks, MPI_INT, MPI_MAX, nczipp->comm);
    zoffs[0] = 0;
    for(i = 1; i < varp->nchunks; i++){
        zoffs[i] = zoffs[i - 1] + zsizes_all[i - 1];
    }

    // Enter redefine mode
    nczipp->driver->redef(nczipp->ncp);

    // Define dimension  for data variable
    sprintf(name, "_compressed_data_dim_%d", varp->varid);
    err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[varp->nchunks - 1] + zsizes_all[varp->nchunks - 1], &zdimid);
    if (err != NC_NOERR) return err;

    // Define variable
    sprintf(name, "_compressed_data_%d", varp->varid);
    err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));
    if (err != NC_NOERR) return err;

    // Record offset of chunks in data variable
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkoffset", NC_INT, varp->nchunks, zoffs, MPI_INT); // Original datatype
    if (err != NC_NOERR) return err;

    // Record size of chunks
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunklen", NC_INT, varp->nchunks, zsizes_all, MPI_INT); // Original datatype
    if (err != NC_NOERR) return err;

    // Switch to data mode
    err = nczipp->driver->enddef(nczipp->ncp);
    if (err != NC_NOERR) return err;

    // Do I/O
    for(l = 0; l < varp->nmychunk; l++){
        k = varp->mychunks[l];
        zstart = (MPI_Offset)zoffs[k];
        zcount = (MPI_Offset)zsizes[k];
        printf("cache[0] = %d, cache[1] = %d, off = %lld, cnt = %lld\n", ((int*)(varp->chunk_cache[k]))[0], ((int*)(varp->chunk_cache[k]))[1], zstart, zcount); fflush(stdout);
        nczipp->driver->iput_var(nczipp->ncp, varp->datavarid, &zstart, &zcount, NULL, NULL, zbufs[l], (MPI_Offset)(zsizes_all[k]), MPI_UNSIGNED_CHAR, NULL, NC_REQ_WR | NC_REQ_NBI | NC_REQ_FLEX);
    }
    nczipp->driver->wait(nczipp->ncp, NC_REQ_ALL, NULL, NULL, NC_REQ_COLL);

    // Free buffers

    NCI_Free(bufs);

    NCI_Free(wcnt_local);
    NCI_Free(wcnt_all);

    NCI_Free(zsizes);
    NCI_Free(zsizes_all);
    NCI_Free(zoffs);
    for(l = 0; l < varp->nmychunk; l++){
        NCI_Free(zbufs[l]);
    }
    NCI_Free(zbufs);

    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);

    NCI_Free(cstart);
    NCI_Free(citr);
    NCI_Free(cend);

    NCI_Free(sbuf_base);

    NCI_Free(sreqs);
    NCI_Free(sstats);

    NCI_Free(rreqs);
    NCI_Free(rstats);
    for(i = 0; i < nrecv; i++){
        NCI_Free(rbufs[i]);
    }
    NCI_Free(rbufs);
    NCI_Free(rsizes);

    return NC_NOERR;
}
