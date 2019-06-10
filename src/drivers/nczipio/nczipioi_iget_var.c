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


/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int nczipioi_iget_cb_chunk(NC_zip *nczipp, int nreq, int *reqids, int *stats){
    int err;
    int i, j;
    int nvar;
    int vid;   // Iterators for variable id
    int *varids;
    int *nreqs;  // Number of reqids in each variable
    int *nums;  // Number of reqs in each varn
    int **vreqids;
    int num, maxnum = 0;
    MPI_Offset **starts, **counts, **strides;
    MPI_Offset rsize;
    char **bufs;
    NC_zip_req *req;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_INIT)
    
    // Count total number of request in per variable for packed varn request
    nums = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt * 2);
    nreqs = nums + nczipp->vars.cnt;
    memset(nums, 0, sizeof(int) * nczipp->vars.cnt);
    memset(nreqs, 0, sizeof(int) * nczipp->vars.cnt);
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        nreqs[req->varid]++;
        nums[req->varid] += req->nreq;
    }

    /* Allocate a skip list of reqids for each vriable
     * At the same time, we find out the number of starts and counts we need to allocate
     */
    vreqids = (int**)NCI_Malloc(sizeof(int*) * nczipp->vars.cnt);
    vreqids[0] = (int*)NCI_Malloc(sizeof(int) * nreq);
    maxnum = 0;
    i = 0;
    nvar = 0;
    for(vid = 0; vid < nczipp->vars.cnt; vid++){
        if (nreqs[vid] > 0){
            // Assign buffer to reqid skip list
            vreqids[vid] = vreqids[0] + i;
            i += nreqs[vid];

            // maximum number of starts and counts we need across all variables
            if (maxnum < nums[vid]){
                maxnum = nums[vid];
            }

            // Number of variable that has request to write
            nvar++;
        }
    }

    varids = (int*)NCI_Malloc(sizeof(int) * nvar);

    // Fill up the skip list
    memset(nreqs, 0, sizeof(int) * nczipp->vars.cnt);
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        vreqids[req->varid][nreqs[req->varid]++] = reqids[i];
    }
    
    // Allocate parameters
    starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * maxnum * 2);
    counts = starts + maxnum;
    bufs =  (char**)NCI_Malloc(sizeof(char*) * maxnum);

    /* Pack requests variable by variable
     */
    nvar = 0;
    for(vid = 0; vid < nczipp->vars.cnt; vid++){
        if (nreqs[vid] > 0){
            // Fill varid in the skip list
            varids[nvar++] = vid;

            // Collect parameters
            num = 0;
            for(j = 0; j < nreqs[vid]; j++){
                req = nczipp->getlist.reqs + vreqids[vid][j];

                if (req->nreq > 1){
                    for(i = 0; i < req->nreq; i++){
                        starts[num] = req->starts[i];
                        counts[num] = req->counts[i];
                        bufs[num++] = req->xbufs[i];
                    }
                }
                else{
                    starts[num] = req->start;
                    counts[num] = req->count;
                    bufs[num++] = req->xbuf;
                }
            }
            
            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB)
            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_INIT)

            // Perform collective buffering
            nczipioi_get_varn_cb_chunk(nczipp, nczipp->vars.data + vid, num, starts, counts, NULL, (void**)bufs);

            NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_INIT)
        }
    }

    // Free buffers
    NCI_Free(nums);

    NCI_Free(vreqids[0]);
    NCI_Free(vreqids);

    NCI_Free(varids);
    
    NCI_Free(starts);
    NCI_Free(bufs);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB)
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_INIT)

    return NC_NOERR;
}

int nczipioi_iget_cb_proc(NC_zip *nczipp, int nreq, int *reqids, int *stats){
    int err;
    int i, j, k;
    int cid, cown;   // Chunk iterator
    int vid;
    int r, **reqs;
    
    MPI_Offset *ostart, *osize;
    int *tsize, *tssize, *tstart, *tssizep, *tstartp;   // Size for sub-array type
    MPI_Offset *citr; // Bounding box for chunks overlapping my own write region
    
    int *rcnt_local, *rcnt_all;   // Number of processes that writes to each proc
    int *rcnt_local_chunk, *rcnt_all_chunk;   // Number of processes that writes to each chunk

    int overlapsize;    // Size of overlaping region of request and chunk
    char *tbuf = NULL;     // Intermediate buffer
    
    int packoff; // Pack offset
    MPI_Datatype ptype; // Pack datatype

    int nsend, nrecv;   // Number of send and receive
    MPI_Request *sreq, *rreq, *sreq_re, *rreq_re;    // Send and recv req
    MPI_Status *sstat, rstat, *sstat_re;    // Send and recv status
    char **sbuf, **sbufp, **rbuf, **rbufp, **sbuf_re, **rbuf_re;   // Send and recv buffer
    int *rsize, *ssize, *rsize_re, *ssize_re;    // recv size of each message
    int *roff, *roff_re;    // recv size of each message
    int *sdst;    // recv size of each message
    int *smap;
    MPI_Message rmsg;   // Receive message
    NC_zip_var *varp;
    NC_zip_req *req;
    
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_INIT)

    // Allocate buffering for write count
    rcnt_local = (int*)NCI_Malloc(sizeof(int) * nczipp->np * 3);
    rcnt_all = rcnt_local + nczipp->np;
    smap = rcnt_all + nczipp->np;

    // Intermediate buffer for our own data
    tbuf = (char*)NCI_Malloc(nczipp->max_chunk_size);

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * nczipp->max_ndim * 3);
    tssize = tsize + nczipp->max_ndim;
    tstart = tssize + nczipp->max_ndim;
    ostart = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * nczipp->max_ndim * 3);
    osize = ostart + nczipp->max_ndim;

    // Chunk iterator
    citr = osize + nczipp->max_ndim;

    // We need to calculate the size of message of each chunk
    // This is just for allocating send buffer
    // We do so by iterating through all request and all chunks they cover
    // If we are not the owner of a chunk, we need to send message
    memset(rcnt_local, 0, sizeof(int) * nczipp->np);
    nsend = 0;

    // req->counts[r] total number of messages and build a map of accessed chunk to list of comm datastructure
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        varp = nczipp->vars.data + req->varid;
        for(r = 0; r < req->nreq; r++){
            nczipioi_chunk_itr_init(varp, req->starts[r], req->counts[r], citr, &cid); // Initialize chunk iterator
            do{
                // Chunk owner
                cown = varp->chunk_owner[cid];

                // Mapping to skip list of send requests 
                if (rcnt_local[cown] == 0 && cown != nczipp->rank){
                    smap[cown] = nsend++;
                }
                rcnt_local[cown] = 1;   // Need to send message if not owner     
            } while (nczipioi_chunk_itr_next(varp, req->starts[r], req->counts[r], citr, &cid));
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_INIT)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_SYNC)

    // Sync number of messages of each chunk
    MPI_Allreduce(rcnt_local, rcnt_all, nczipp->np, MPI_INT, MPI_SUM, nczipp->comm);
    nrecv = rcnt_all[nczipp->rank] - rcnt_local[nczipp->rank];  // We don't need to receive request form self

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_SYNC)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_PACK_REQ)

    // Allocate data structure for messaging
    sbuf = (char**)NCI_Malloc(sizeof(char*) * (nsend * 2 + nrecv));
    sbufp = sbuf + (nsend + nrecv);
    ssize = (int*)NCI_Malloc(sizeof(int) * (nsend * 2 + nrecv * 1));
    sdst = ssize + (nsend + nrecv);
    sreq = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * (nsend + nrecv));
    sstat = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * (nsend + nrecv));
    reqs = (int**)NCI_Malloc(sizeof(int*) * nsend);

    rbuf = (char**)NCI_Malloc(sizeof(char*) * (nsend + nrecv * 2));
    rbufp = rbuf + (nsend + nrecv);
    rsize = (int*)NCI_Malloc(sizeof(int) * (nsend + nrecv));
    rreq = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * (nsend + nrecv));

    sbuf_re = sbuf + nsend;
    ssize_re = ssize + nsend;
    sreq_re = sreq + nsend;
    sstat_re = sstat + nsend;

    rbuf_re = rbuf + nrecv;
    rsize_re = rsize + nrecv;
    rreq_re = rreq + nrecv;

    // req->counts[r] size of each request
    memset(ssize, 0, sizeof(int) * nsend);
    memset(rsize_re, 0, sizeof(int) * nsend);
    memset(rcnt_local, 0, sizeof(int) * nsend);
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        varp = nczipp->vars.data + req->varid;
        for(r = 0; r < req->nreq; r++){
            nczipioi_chunk_itr_init_ex(varp, req->starts[r], req->counts[r], citr, &cid, ostart, osize); // Initialize chunk iterator
            do{
                // Chunk owner
                cown = varp->chunk_owner[cid];
                if (cown != nczipp->rank){
                    j = smap[cown];
                    sdst[j] = cown; // Record a reverse map by the way

                    // req->counts[r] overlap
                    //get_chunk_overlap(varp, citr, req->starts[r], req->counts[r], ostart, osize);
                    overlapsize = varp->esize;
                    for(k = 0; k < varp->ndim; k++){
                        overlapsize *= osize[k];                     
                    }
                    ssize[j] += sizeof(int) * (varp->ndim * 2 + 2);
                    rsize_re[j] += overlapsize;
                    rcnt_local[j]++;
                }
            } while (nczipioi_chunk_itr_next_ex(varp, req->starts[r], req->counts[r], citr, &cid, ostart, osize));
        }
    }

    // Allocate buffer for send
    //memset(soff, 0, sizeof(int) * (nsend + nrecv));
    for(i = 0; i < nsend; i++){
        ssize[i] += sizeof(int);
        sbuf[i] = sbufp[i] = (char*)NCI_Malloc(ssize[i]);
        
        //CHK_ERR_PACK(rsize_re + i, 1, MPI_INT, sbuf[i], ssize[i], soff + i, nczipp->comm);
        *((int*)sbufp[i]) = rsize_re[i];    sbufp[i] += sizeof(int);

        rbuf_re[i] = (char*)NCI_Malloc(rsize_re[i]);
        reqs[i] = (int*)NCI_Malloc(sizeof(int) * rcnt_local[i] * 2);
    }

    // Pack requests
    memset(rcnt_local, 0, sizeof(int) * nsend);
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        varp = nczipp->vars.data + req->varid;
        for(r = 0; r < req->nreq; r++){
            nczipioi_chunk_itr_init_ex(varp, req->starts[r], req->counts[r], citr, &cid, ostart, osize); // Initialize chunk iterator
            do{
                // Chunk owner
                cown = varp->chunk_owner[cid];
                if (cown != nczipp->rank){
                    j = smap[cown];

                    // Get overlap region
                    //get_chunk_overlap(varp, citr, req->starts[r], req->counts[r], ostart, osize);

                    // Pack metadata
                    for(k = 0; k < varp->ndim; k++){
                        tstart[k] = (int)(ostart[k] - citr[k]);
                        tssize[k] = (int)osize[k];
                    }
                    /*
                    CHK_ERR_PACK(&(varp->varid), 1, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
                    CHK_ERR_PACK(&cid, 1, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
                    CHK_ERR_PACK(tstart, varp->ndim, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
                    CHK_ERR_PACK(tssize, varp->ndim, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
                    */
                    *((int*)sbufp[j]) = varp->varid;    sbufp[j] += sizeof(int);
                    *((int*)sbufp[j]) = cid;    sbufp[j] += sizeof(int);
                    memcpy(sbufp[j], tstart, sizeof(int) * varp->ndim); sbufp[j] += sizeof(int) * varp->ndim;
                    memcpy(sbufp[j], tssize, sizeof(int) * varp->ndim); sbufp[j] += sizeof(int) * varp->ndim;

                    // Record source of the request
                    reqs[j][rcnt_local[j]++] = i;
                    reqs[j][rcnt_local[j]++] = r;
                }
            } while (nczipioi_chunk_itr_next_ex(varp, req->starts[r], req->counts[r], citr, &cid, ostart, osize));
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_PACK_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_SEND_REQ)

    // Post send 
    for(i = 0; i < nsend; i++){
        CHK_ERR_ISEND(sbuf[i], ssize[i], MPI_BYTE, sdst[i], 0, nczipp->comm, sreq + i);
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_SEND_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_RECV_REP)

    // Post receive  
    for(i = 0; i < nsend; i++){
        MPI_Irecv(rbuf_re[i], rsize_re[i], MPI_BYTE, sdst[i], 1, nczipp->comm, rreq_re + i);
    }   

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_RECV_REP)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_RECV_REQ)

    // Post recv
    for(i = 0; i < nrecv; i++){
        // Get message size, including metadata
        MPI_Mprobe(MPI_ANY_SOURCE, 0, nczipp->comm, &rmsg, &rstat);
        MPI_Get_count(&rstat, MPI_BYTE, rsize + i);

        // Allocate buffer
        rbuf[i] = rbufp[i] = (char*)NCI_Malloc(rsize[i]);

        // Post irecv
        MPI_Imrecv(rbuf[i], rsize[i], MPI_BYTE, &rmsg, rreq + i);
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_RECV_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_SELF)

    // Handle our own data
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        varp = nczipp->vars.data + req->varid;
        for(r = 0; r < req->nreq; r++){
            nczipioi_chunk_itr_init_ex(varp, req->starts[r], req->counts[r], citr, &cid, ostart, osize); // Initialize chunk iterator
            do{
                if (varp->chunk_owner[cid] == nczipp->rank){  
                    // Pack type from chunk cache to (contiguous) intermediate buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - citr[j]);
                        tsize[j] = varp->chunkdim[j];
                        tssize[j] = (int)osize[j];
                    }
                    CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    CHK_ERR_TYPE_COMMIT(&ptype);

                    // Pack data into intermediate buffer
                    packoff = 0;
                    CHK_ERR_PACK(varp->chunk_cache[cid], 1, ptype, tbuf, varp->chunksize, &packoff, nczipp->comm);
                    MPI_Type_free(&ptype);
                    overlapsize = packoff;

                    // Pack type from (contiguous) intermediate buffer to chunk buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - req->starts[r][j]);
                        tsize[j] = (int)req->counts[r][j];
                    }
                    CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    CHK_ERR_TYPE_COMMIT(&ptype);
                    
                    // Unpack data into chunk buffer
                    packoff = 0;
                    CHK_ERR_UNPACK(tbuf, overlapsize, &packoff, req->xbufs[r], 1, ptype, nczipp->comm);
                    MPI_Type_free(&ptype);    
                }
            } while (nczipioi_chunk_itr_next_ex(varp, req->starts[r], req->counts[r], citr, &cid, ostart, osize));
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_SELF)

    //Handle incoming requests
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = varp->chunkdim[i];
    }
    for(i = 0; i < nrecv; i++){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_RECV_REQ)

        // Will wait any provide any benefit?
        MPI_Waitany(nrecv, rreq, &j, &rstat);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_RECV_REQ)

        packoff = 0;
        ssize_re[j] = *((int*)rbufp[j]);    rbufp[j] += sizeof(int);
        sbuf_re[j] = (char*)NCI_Malloc(ssize_re[j]);
        while(rbufp[j] < rbuf[j] + rsize[j]){
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_UNPACK_REQ)

            // Retrieve metadata
            vid = *((int*)rbufp[j]);    rbufp[j] += sizeof(int);
            cid = *((int*)rbufp[j]);    rbufp[j] += sizeof(int);
            varp = nczipp->vars.data + vid;
            tstartp = (int*)rbufp[j];    rbufp[j] += sizeof(int) * varp->ndim;
            tssizep = (int*)rbufp[j];    rbufp[j] += sizeof(int) * varp->ndim;

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_UNPACK_REQ)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_PACK_REP)

            // Pack type
            CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssizep, tstartp, MPI_ORDER_C, varp->etype, &ptype);
            CHK_ERR_TYPE_COMMIT(&ptype);

            // Pack data
            CHK_ERR_PACK(varp->chunk_cache[cid], 1, ptype, sbuf_re[j], ssize_re[j], &packoff, nczipp->comm);
            MPI_Type_free(&ptype);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_PACK_REP)
        }

        NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_SEND_REQ)

        // Send Response
        CHK_ERR_ISEND(sbuf_re[j], packoff, MPI_BYTE, rstat.MPI_SOURCE, 1, nczipp->comm, sreq_re + j);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_SEND_REQ)
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_SEND_REQ)

    // Wait for all request
    MPI_Waitall(nsend, sreq, sstat);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_SEND_REQ)

    //Handle reply
    memset(rcnt_local, 0, sizeof(int) * nsend);
    for(i = 0; i < nsend; i++){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_RECV_REP)

        // Will wait any provide any benefit?
        MPI_Waitany(nsend, rreq_re, &j, &rstat);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_RECV_REP)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_UNPACK_REP)

        //soff[j] = sizeof(int);  // Skip reply size
        sbufp[j] = sbuf[j] + sizeof(int);
        packoff = 0;
        while(packoff < rsize_re[j]){
            // Retrieve metadata from the request we sent
            vid = *((int*)sbufp[j]);    sbufp[j] += sizeof(int);
            cid = *((int*)sbufp[j]);    sbufp[j] += sizeof(int);
            varp = nczipp->vars.data + vid;
            tstartp = (int*)sbufp[j];    sbufp[j] += sizeof(int) * varp->ndim;
            tssizep = (int*)sbufp[j];    sbufp[j] += sizeof(int) * varp->ndim;

            k = reqs[j][rcnt_local[j]++];
            r = reqs[j][rcnt_local[j]++];
            req = nczipp->getlist.reqs + reqids[k];
            get_chunk_itr(varp, cid, citr);
            for(k = 0; k < varp->ndim; k++){
                tstartp[k] += (int)(citr[k] - req->starts[r][k]);
                tsize[k] = req->counts[r][k];
            }

            // Pack type
            CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssizep, tstartp, MPI_ORDER_C, varp->etype, &ptype);
            CHK_ERR_TYPE_COMMIT(&ptype);

            // Pack data
            //printf("Rank: %d, cid = %d, MPI_Unpack(%d, %d, %d, %d)\n", nczipp->rank, cid, j, rsize_re[j], packoff, req); fflush(stdout);
            MPI_Unpack(rbuf_re[j], rsize_re[j], &packoff, req->xbufs[r], 1, ptype, nczipp->comm);
            //printf("cache[0] = %d, cache[1] = %d\n", ((int*)(varp->chunk_cache[cid]))[0], ((int*)(varp->chunk_cache[cid]))[1]); fflush(stdout);
            MPI_Type_free(&ptype);
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_UNPACK_REP)
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET_CB_SEND_REP)

    // Wait for all Response
    MPI_Waitall(nrecv, sreq_re, sstat_re);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB_SEND_REP)

    // Free buffers
    NCI_Free(rcnt_local);

    NCI_Free(tsize);

    NCI_Free(ostart);

    NCI_Free(tbuf);

    NCI_Free(sreq);
    NCI_Free(sstat);
    NCI_Free(ssize);
    for(i = 0; i < nsend; i++){
        NCI_Free(reqs[i]);
    }
    for(i = 0; i < nsend + nrecv; i++){
        NCI_Free(sbuf[i]);
        NCI_Free(rbuf[i]);
    }
    NCI_Free(sbuf);
    NCI_Free(reqs);

    NCI_Free(rreq);
    NCI_Free(rbuf);
    NCI_Free(rsize);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET_CB)

    return NC_NOERR;
}
