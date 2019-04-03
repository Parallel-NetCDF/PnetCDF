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

int
nczipioi_get_varn_cb_chunk(NC_zip          *nczipp,
                    NC_zip_var      *varp,
                    int              nreq,
                    MPI_Offset* const *starts,
                    MPI_Offset* const *counts,
                    MPI_Offset* const *strides,
                    void            **bufs)
{
    int err;
    int i, j, k, l;
    int cid, req;   // Chunk iterator

    MPI_Offset *ostart, *osize;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int *cstart, *cend, *citr; // Bounding box for chunks overlapping my own write region
    
    int *rcnt_local, *rcnt_all;   // Number of processes that writes to each chunk

    int overlapsize;    // Size of overlaping region of request and chunk
    int overlapsize_total, overlapcnt;
    char *cbuf = NULL;     // Intermediate continuous buffer
    
    int packoff, unpackoff; // Pack offset
    MPI_Datatype ptype; // Pack datatype

    int nread;  // # chunks to read form file
    int *rids;  // Id of chunks to read from file

    int nsend, nrecv;   // Number of send and receive
    MPI_Request *sreqs, *rreqs;    // Send and recv req
    MPI_Status *sstats, *rstats;    // Send and recv status
    char **sbufs, **rbufs;   // Send and recv buffer
    int *rsizes;    // recv size of each message
    MPI_Message rmsg;   // Receive message

    // Allocate buffering for write count
    rcnt_local = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    rcnt_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    ostart = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
    osize = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);

    // Starting, ending, current chunk position
    cstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    citr = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    cend = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    // We need to calculate the size of message of each chunk
    // This is just for allocating send buffer
    // We do so by iterating through all request and all chunks they cover
    // If we are not the owner of a chunk, we need to send message
    memset(rcnt_local, 0, sizeof(int) * varp->nchunks);
    nsend = 0;
    for(req = 0; req < nreq; req++){
        // Initialize chunk iterator
        nczipioi_chunk_itr_init(varp, starts[req], counts[req], cstart, cend, citr);

        // Iterate through chunks
        do{
            // Chunk index
            cid = get_chunk_idx(varp, citr);
            
            if (varp->chunk_owner[cid] != nczipp->rank && rcnt_local[cid] == 0){
                // Count number of mnessage we need to send
                nsend++;    
            }

            rcnt_local[cid] = 1;
        } while (nczipioi_chunk_itr_next(varp, starts[req], counts[req], cstart, cend, citr));
    }

    // Sync number of messages of each chunk
    MPI_Allreduce(rcnt_local, rcnt_all, varp->nchunks, MPI_INT, MPI_SUM, nczipp->comm);

    // We need to prepare chunk in the chunk cache
    // For chunks not yet allocated, we need to read them form file collectively
    // We collect chunk id of those chunks
    // Calculate number of recv request
    // This is for all the chunks
    rids = (int*)NCI_Malloc(sizeof(int) * varp->nmychunks);
    nread = 0;
    nrecv = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];
        // We don't need message for our own data
        nrecv += rcnt_all[cid] - rcnt_local[cid];
        // Count number of chunks we need to prepare
        // We read only chunks that is required
        if (rcnt_all[cid] > 0 && varp->chunk_cache[cid] == NULL){
            rids[nread++] = cid;
        }
    }

    // Decompress chunks into chunk cache
    nczipioi_load_var(nczipp, varp, nread, rids);

    // Allocate buffer for send and recv
    // We need to accept nrecv requests and receive nsend of replies
    rreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * (nrecv + nsend));
    rstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * (nrecv + nsend));
    rbufs = (char**)NCI_Malloc(sizeof(char*) * (nrecv + nsend));
    rsizes = (int*)NCI_Malloc(sizeof(int) * (nrecv + nsend));
    // We need to send nsend requests and reply nrecv of requests
    sbufs = (char**)NCI_Malloc(sizeof(char*) * (nrecv + nsend));
    sreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * (nrecv + nsend));
    sstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * (nrecv + nsend));

    // Post send
    k = l = 0;
    for(cid = 0; cid < varp->nchunks; cid++){
        if (varp->chunk_owner[cid] == nczipp->rank){
            // We are the owner of the chunk
            // Receive data from other process
            for(j = 0; j < rcnt_all[cid] - rcnt_local[cid]; j++){
                // Get message size, including metadata
                MPI_Mprobe(MPI_ANY_SOURCE, cid, nczipp->comm, &rmsg, rstats);
                MPI_Get_count(rstats, MPI_BYTE, rsizes + k);

                //printf("rsize = %d\n", rsizes[i]); fflush(stdout);

                // Allocate buffer
                rbufs[k] = (char*)NCI_Malloc(rsizes[k]);

                // Post irecv
                //printf("Rank: %d, MPI_Imrecv(%d, %d, %d, %d)\n", nczipp->rank, rsizes[k], -1, cid, k); fflush(stdout);
                MPI_Imrecv(rbufs[k], rsizes[k], MPI_BYTE, &rmsg, rreqs + k);
                k++;
            }
        }
        else{
            // We have some request to send
            if (rcnt_local[cid] > 0){
                get_chunk_cord(varp, cid, citr);
                rsizes[nrecv + l] = overlapcnt = 0;

                // Calculate send buffer size
                for(req = 0; req < nreq; req++){
                    // Calculate chunk overlap
                    get_chunk_overlap(varp, citr, starts[req], counts[req], ostart, osize);
                    
                    overlapsize = varp->esize;
                    for(j = 0; j < varp->ndim; j++){
                        overlapsize *= osize[j];                     
                    }
                    rsizes[nrecv + l]  += overlapsize;

                    if (overlapsize > 0){
                        overlapcnt++;
                    }
                }

                // Allocate buffer
                // Faster to request the entire chunk
                if (rsizes[nrecv + l]  >= varp->chunksize){
                    rsizes[nrecv + l]  = varp->chunksize;
                    overlapcnt = 1;
                }
                sbufs[l] = (char*)NCI_Malloc(sizeof(int) * (overlapcnt * varp->ndim * 2) + 1);
                rbufs[nrecv + l] = (char*)NCI_Malloc(rsizes[nrecv + l]);

                // Pack metadata
                packoff = 0;
                MPI_Pack(rsizes + nrecv + l, 1, MPI_INT, sbufs[l], packoff + sizeof(int), &packoff, nczipp->comm);  // Include total buffer size needed
                if (rsizes[nrecv + l] == varp->chunksize){
                    // Metadata
                    memset(tstart, 0, sizeof(int) * varp->ndim);
                    MPI_Pack(tstart, varp->ndim, MPI_INT, sbufs[l], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                    MPI_Pack(varp->chunkdim, varp->ndim, MPI_INT, sbufs[l], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                }
                else{
                    for(req = 0; req < nreq; req++){
                        // Calculate chunk overlap
                        get_chunk_overlap(varp, citr, starts[req], counts[req], ostart, osize);
                        //printf("cord = %d, start = %lld, count = %lld, tstart = %d, tssize = %d, esize = %d, ndim = %d\n", citr[0], starts[req][0], counts[req][0], tstart[0], tssize[0], varp->esize, varp->ndim); fflush(stdout);
                        overlapsize = varp->esize;
                        for(j = 0; j < varp->ndim; j++){
                            overlapsize *= osize[j];                     
                        }

                        if (overlapsize > 0){
                            // Metadata
                            for(j = 0; j < varp->ndim; j++){
                                tstart[j] = (int)(ostart[j] - citr[j] * varp->chunkdim[j]);
                                tsize[j] = (int)osize[j];
                            }
                                    
                            // Pack metadata   
                            MPI_Pack(tstart, varp->ndim, MPI_INT, sbufs[l], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                            MPI_Pack(tsize, varp->ndim, MPI_INT, sbufs[l], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        }
                    }
                }

                // Send request
                MPI_Isend(sbufs[l], packoff, MPI_BYTE, varp->chunk_owner[cid], cid, nczipp->comm, sreqs + l);
                //printf("Rank: %d, MPI_Irecv(%d, %d, %d, %d)\n", nczipp->rank, overlapsize, varp->chunk_owner[cid], cid + 1024, nrecv + k); fflush(stdout);
                MPI_Irecv(rbufs[l + nrecv], rsizes[nrecv + l] , MPI_BYTE, varp->chunk_owner[cid], cid + 1024, nczipp->comm, rreqs + nrecv + l);

                l++;
            } 
        }
    }

    // Allocate intermediate buffer
    cbuf = (char*)NCI_Malloc(varp->chunksize);

    // For each chunk we own, we need to reply to incoming reqeust
    k = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];
            
        // Handle our own data first if we have any
        if (rcnt_local[cid] > 0){
            // Convert chunk id to iterator
            get_chunk_cord(varp, cid, citr);

            for(req = 0; req < nreq; req++){
                // Calculate overlapping region
                get_chunk_overlap(varp, citr, starts[req], counts[req], ostart, osize);
                overlapsize = varp->esize;
                for(j = 0; j < varp->ndim; j++){
                    overlapsize *= osize[j];                     
                }

                if (overlapsize > 0){
                    // Pack type from chunk buffer to (contiguous) intermediate buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - citr[j] * varp->chunkdim[j]);
                        tsize[j] = varp->chunkdim[j];
                        tssize[j] = (int)osize[j];
                    }
                    //printf("Rank: %d, ostart=[%lld, %lld], osize=[%lld, %lld]\n", nczipp->rank, ostart[0], ostart[1], osize[0], osize[1]); fflush(stdout);
                    //printf("Rank: %d, MPI_Type_create_subarray1([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    MPI_Type_commit(&ptype);

                    // Pack data into intermediate buffer
                    packoff = 0;
                    MPI_Pack(varp->chunk_cache[cid], 1, ptype, cbuf, varp->chunksize, &packoff, nczipp->comm);
                    overlapsize = packoff;
                    MPI_Type_free(&ptype);

                    // Pack type from (contiguous) intermediate buffer to user buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - starts[req][j]);
                        tsize[j] = (int)counts[req][j];
                    }
                    //printf("Rank: %d, MPI_Type_create_subarray([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    MPI_Type_commit(&ptype);

                    // Pack data into user buffer
                    packoff = 0;
                    //printf("Rank: %d, cid = %d, MPI_Unpack_self(%d, %d)\n", nczipp->rank, cid, overlapsize, packoff); fflush(stdout);
                    MPI_Unpack(cbuf, overlapsize, &packoff, bufs[req], 1, ptype, nczipp->comm);
                    MPI_Type_free(&ptype);
                }
            }
        }

        // Now, it is time to process data from other processes
        for(j = 0; j < varp->ndim; j++){
            tsize[j] = varp->chunkdim[j];
        }

        // Wait for all send requests related to this chunk
        // We remove the impact of -1 mark in rcnt_local[cid]
        //printf("Rank: %d, MPI_Waitall_recv(%d, %d)\n", nczipp->rank, rcnt_all[cid] - rcnt_local[cid], k); fflush(stdout);
        MPI_Waitall(rcnt_all[cid] - rcnt_local[cid], rreqs + k, rstats + k);

        // Process data received
        //printf("nrecv = %d, rcnt_all = %d, rcnt_local = %d\n", nrecv, rcnt_all[cid], rcnt_local[cid]); fflush(stdout);
        for(j = k; j < k + rcnt_all[cid] - rcnt_local[cid]; j++){
            packoff = unpackoff = 0;
            
            // Allocate buffer 
            //printf("Rank: %d, MPI_Unpack_bufsize(%d, %d, %d)\n", nczipp->rank, j, rsizes[j], unpackoff); fflush(stdout);
            MPI_Unpack(rbufs[j], rsizes[j], &unpackoff, &overlapsize, 1, MPI_INT, nczipp->comm);
            sbufs[j + nsend] = (char*)NCI_Malloc(overlapsize); // For reply

            // Pack data
            while(unpackoff < rsizes[j]){
                // Get metadata
                //printf("Rank: %d, MPI_Unpack_meta(%d, %d, %d)\n", nczipp->rank, j, rsizes[j], unpackoff); fflush(stdout);
                MPI_Unpack(rbufs[j], rsizes[j], &unpackoff, tstart, varp->ndim, MPI_INT, nczipp->comm);
                MPI_Unpack(rbufs[j], rsizes[j], &unpackoff, tssize, varp->ndim, MPI_INT, nczipp->comm);

                // Pack type
                //printf("Rank: %d, MPI_Type_create_subarray([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                MPI_Type_commit(&ptype);

                // Pack data
                MPI_Pack(varp->chunk_cache[cid], 1, ptype, sbufs[j + nsend], overlapsize, &packoff, nczipp->comm);
                MPI_Type_free(&ptype);
            }

            // Send reply
            //printf("Rank: %d, MPI_Isend(%d, %d, %d, %d)\n", nczipp->rank, packoff, varp->chunk_owner[cid], cid + 1024, k + nsend); fflush(stdout);
            MPI_Isend(sbufs[j + nsend], packoff, MPI_BYTE, rstats[j].MPI_SOURCE, cid + 1024, nczipp->comm, sreqs + j + nsend);
        }
        k += rcnt_all[cid] - rcnt_local[cid];        

        //princbuf(nczipp->rank, varp->chunk_cache[cid], varp->chunksize);
    }

    // Wait for all request sent
    //printf("Rank: %d, MPI_Waitall_send(%d, %d)\n", nczipp->rank, nsend, 0); fflush(stdout);
    MPI_Waitall(nsend, sreqs, sstats);

    // Receive replies from the owners and update the user buffer
    k = 0;
    for(cid = 0; cid < varp->nchunks; cid++){
        if (rcnt_local[cid] > 0 && varp->chunk_owner[cid] != nczipp->rank){
            get_chunk_cord(varp, cid, citr);

            // Wait for reply
            //printf("Rank: %d, MPI_Wait_recv(%d)\n", nczipp->rank, nrecv + k); fflush(stdout);
            MPI_Wait(rreqs + nrecv + k, rstats + nrecv + k);

            packoff = 0;
            for(req = 0; req < nreq; req++){
                // Calculate chunk overlap
                get_chunk_overlap(varp, citr, starts[req], counts[req], ostart, osize);
                overlapsize = varp->esize;
                for(j = 0; j < varp->ndim; j++){
                    overlapsize *= osize[j];                     
                }

                if (overlapsize > 0){
                    // Pack type from recv buffer to user buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - starts[req][j]);
                        tsize[j] = (int)counts[req][j];
                        tssize[j] = (int)osize[j];
                    }
                    //printf("Rank: %d, ostart=[%lld, %lld], osize=[%lld, %lld]\n", nczipp->rank, ostart[0], ostart[1], osize[0], osize[1]); fflush(stdout);
                    //printf("Rank: %d, MPI_Type_create_subarray4([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    //printf("Rank: %d, commit\n", nczipp->rank); fflush(stdout);
                    MPI_Type_commit(&ptype);

                    //printf("Rank: %d, wait recv, nrecv = %d, k = %d, nsend = %d\n", nczipp->rank, nrecv, k, nsend); fflush(stdout);

                    // Pack data
                    //printf("Rank: %d, cid = %d, MPI_Unpack(%d, %d, %d, %d)\n", nczipp->rank, cid, nrecv + k, rsizes[nrecv + k], packoff, req); fflush(stdout);
                    MPI_Unpack(rbufs[nrecv + k], rsizes[nrecv + k], &packoff, bufs[req], 1, ptype, nczipp->comm);
                    //printf("Rank: %d, cid = %d, MPI_Unpack_done(%d, %d, %d, %d)\n", nczipp->rank, cid, nrecv + k, rsizes[nrecv + k], packoff, req); fflush(stdout);
                    MPI_Type_free(&ptype);
                }
            }
            k++;
        }
    }

    //printf("Rank: %d, wait_final\n", nczipp->rank); fflush(stdout);
    // Wait for all send replies
    //printf("Rank: %d, MPI_Waitall_send(%d, %d)\n", nczipp->rank, nrecv, nsend); fflush(stdout);
    MPI_Waitall(nrecv, sreqs + nsend, sstats + nsend);

    //printf("Rank: %d, exiting\n", nczipp->rank); fflush(stdout);

    // Free buffers
    NCI_Free(rcnt_local);
    NCI_Free(rcnt_all);

    NCI_Free(rids);

    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);
    NCI_Free(osize);
    NCI_Free(ostart);

    NCI_Free(cstart);
    NCI_Free(citr);
    NCI_Free(cend);

    for(i = 0; i < nsend + nrecv; i++){
        NCI_Free(sbufs[i]);
        NCI_Free(rbufs[i]);
    }
    NCI_Free(sreqs);
    NCI_Free(sstats);
    NCI_Free(sbufs);
    NCI_Free(rreqs);
    NCI_Free(rstats);
    NCI_Free(rbufs);
    NCI_Free(rsizes);

    if (cbuf != NULL){
        NCI_Free(cbuf);
    }

    return NC_NOERR;
}

int
nczipioi_get_varn_cb_proc(  NC_zip          *nczipp,
                            NC_zip_var      *varp,
                            int              nreq,
                            MPI_Offset* const *starts,
                            MPI_Offset* const *counts,
                            void            **bufs)
{
    int err;
    int i, j, k;
    int cid, cown;   // Chunk iterator
    int req, **reqs;

    MPI_Offset *ostart, *osize;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    MPI_Offset *citr; // Bounding box for chunks overlapping my own write region
    
    int *rcnt_local, *rcnt_all;   // Number of processes that writes to each proc
    int *rcnt_local_chunk, *rcnt_all_chunk;   // Number of processes that writes to each chunk

    int overlapsize;    // Size of overlaping region of request and chunk
    int max_tbuf = 0;   // Size of intermediate buffer
    char *tbuf = NULL;     // Intermediate buffer
    
    int packoff; // Pack offset
    MPI_Datatype ptype; // Pack datatype

    int nread;  // # chunks to read form file
    int *rids;  // Id of chunks to read from file

    int nsend, nrecv;   // Number of send and receive
    MPI_Request *sreq, *rreq, *sreq_re, *rreq_re;    // Send and recv req
    MPI_Status *sstat, rstat, *sstat_re;    // Send and recv status
    char **sbuf, **rbuf, **sbuf_re, **rbuf_re;   // Send and recv buffer
    int *rsize, *ssize, *rsize_re, *ssize_re;    // recv size of each message
    int *roff, *soff, *roff_re, *soff_re;    // recv size of each message
    int *sdst;    // recv size of each message
    int *smap;
    MPI_Message rmsg;   // Receive message

    // Allocate buffering for write counts[req]
    rcnt_local = (int*)NCI_Malloc(sizeof(int) * (nczipp->np + varp->nchunks));
    rcnt_local_chunk = rcnt_local + nczipp->np;
    rcnt_all = (int*)NCI_Malloc(sizeof(int) * (nczipp->np + varp->nchunks));
    rcnt_all_chunk = rcnt_all + nczipp->np;
    smap = (int*)NCI_Malloc(sizeof(int) * nczipp->np);

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    ostart = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
    osize = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);

    // Current chunk position
    citr = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);

    // We need to calculate the size of message of each chunk
    // This is just for allocating send buffer
    // We do so by iterating through all request and all chunks they cover
    // If we are not the owner of a chunk, we need to send message
    memset(rcnt_local, 0, sizeof(int) * (nczipp->np + varp->nchunks));
    nsend = 0;

    // counts[req] total number of messages and build a map of accessed chunk to list of comm datastructure
    for(req = 0; req < nreq; req++){
        nczipioi_chunk_itr_init_cord(varp, starts[req], counts[req], citr); // Initialize chunk iterator
        do{
            // Chunk index and owner
            cid = get_chunk_idx_cord(varp, citr);
            cown = varp->chunk_owner[cid];

            // Mapping to skip list of send requests 
            if (rcnt_local[cown] == 0 && cown != nczipp->rank){
                smap[cown] = nsend++;
            }
            rcnt_local[cown] = 1;   // Need to send message if not owner     
            rcnt_local_chunk[cid] = 1;  // This tells the owner to prepare the chunks  
        } while (nczipioi_chunk_itr_next_cord(varp, starts[req], counts[req], citr));
    }

    // Sync number of messages of each chunk
    MPI_Allreduce(rcnt_local, rcnt_all, nczipp->np + varp->nchunks, MPI_INT, MPI_SUM, nczipp->comm);
    nrecv = rcnt_all[nczipp->rank] - rcnt_local[nczipp->rank];  // We don't need to receive request form self

    // We need to prepare chunk in the chunk cache
    // For chunks not yet allocated, we need to read them form file collectively
    // We collect chunk id of those chunks
    // Calculate number of recv request
    // This is for all the chunks
    rids = (int*)NCI_Malloc(sizeof(int) * varp->nmychunks);
    nread = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];
        // counts[req] number of chunks we need to prepare
        // We read only chunks that is required
        if (rcnt_all_chunk[cid] > 0 && varp->chunk_cache[cid] == NULL){
            rids[nread++] = cid;
        }
    }

    // Decompress chunks into chunk cache
    nczipioi_load_var(nczipp, varp, nread, rids);

    // Allocate data structure for messaging
    sbuf = (char**)NCI_Malloc(sizeof(char*) * (nsend + nrecv));
    ssize = (int*)NCI_Malloc(sizeof(int) * (nsend + nrecv));
    soff = (int*)NCI_Malloc(sizeof(int) * (nsend + nrecv));
    sdst = (int*)NCI_Malloc(sizeof(int) * nsend);
    sreq = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * (nsend + nrecv));
    sstat = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * (nsend + nrecv));
    reqs = (int**)NCI_Malloc(sizeof(int*) * nsend);

    rbuf = (char**)NCI_Malloc(sizeof(char*) * (nsend + nrecv));
    rsize = (int*)NCI_Malloc(sizeof(int) * (nsend + nrecv));
    rreq = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * (nsend + nrecv));

    sbuf_re = sbuf + nsend;
    ssize_re = ssize + nsend;
    soff_re = soff + nsend;
    sreq_re = sreq + nsend;
    sstat_re = sstat + nsend;

    rbuf_re = rbuf + nrecv;
    rsize_re = rsize + nrecv;
    rreq_re = rreq + nrecv;

    // counts[req] size of each request
    memset(ssize, 0, sizeof(int) * nsend);
    memset(rsize_re, 0, sizeof(int) * nsend);
    memset(rcnt_local, 0, sizeof(int) * nsend);
    for(req = 0; req < nreq; req++){
        nczipioi_chunk_itr_init_cord(varp, starts[req], counts[req], citr); // Initialize chunk iterator
        do{
            // Chunk index and owner
            cid = get_chunk_idx_cord(varp, citr);
            cown = varp->chunk_owner[cid];
            if (cown != nczipp->rank){
                j = smap[cown];
                sdst[j] = cown; // Record a reverse map by the way

                // counts[req] overlap
                get_chunk_overlap_cord(varp, citr, starts[req], counts[req], ostart, osize);
                overlapsize = varp->esize;
                for(i = 0; i < varp->ndim; i++){
                    overlapsize *= osize[i];                     
                }
                ssize[j] += sizeof(int) * (varp->ndim * 2 + 1);
                rsize_re[j] += overlapsize;
                rcnt_local[j]++;
            }
        } while (nczipioi_chunk_itr_next_cord(varp, starts[req], counts[req], citr));
    }

    // Allocate buffer for send
    memset(soff, 0, sizeof(int) * (nsend + nrecv));
    for(i = 0; i < nsend; i++){
        ssize[i] += sizeof(int);
        sbuf[i] = (char*)NCI_Malloc(ssize[i]);
        MPI_Pack(rsize_re + i, 1, MPI_INT, sbuf[i], ssize[i], soff + i, nczipp->comm);
        rbuf_re[i] = (char*)NCI_Malloc(rsize_re[i]);
        reqs[i] = (int*)NCI_Malloc(sizeof(int) * rcnt_local[i]);
    }

    // Pack requests
    memset(rcnt_local, 0, sizeof(int) * nsend);
    for(req = 0; req < nreq; req++){
        nczipioi_chunk_itr_init_cord(varp, starts[req], counts[req], citr); // Initialize chunk iterator
        do{
            // Chunk index and owner
            cid = get_chunk_idx_cord(varp, citr);
            cown = varp->chunk_owner[cid];
            if (cown != nczipp->rank){
                j = smap[cown];

                // Get overlap region
                get_chunk_overlap_cord(varp, citr, starts[req], counts[req], ostart, osize);

                // Pack metadata
                for(i = 0; i < varp->ndim; i++){
                    tstart[i] = (int)(ostart[i] - citr[i]);
                    tssize[i] = (int)osize[i];
                }
                MPI_Pack(&cid, 1, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
                MPI_Pack(tstart, varp->ndim, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
                MPI_Pack(tssize, varp->ndim, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);

                // Record source of the request
                reqs[j][rcnt_local[j]++] = req;
            }
        } while (nczipioi_chunk_itr_next_cord(varp, starts[req], counts[req], citr));
    }

    // Post send and receive
    for(i = 0; i < nsend; i++){
        MPI_Isend(sbuf[i], soff[i], MPI_BYTE, sdst[i], 0, nczipp->comm, sreq + i);
        //printf("Rank: %d, MPI_Irecv(%d, %d, %d, %d)\n", nczipp->rank, rsize_re[i], sdst[i], 1, i); fflush(stdout);
        MPI_Irecv(rbuf_re[i], rsize_re[i], MPI_BYTE, sdst[i], 1, nczipp->comm, rreq_re + i);
    }   

    // Post recv
    for(i = 0; i < nrecv; i++){
        // Get message size, including metadata
        MPI_Mprobe(MPI_ANY_SOURCE, 0, nczipp->comm, &rmsg, &rstat);
        MPI_Get_count(&rstat, MPI_BYTE, rsize + i);

        // Allocate buffer
        rbuf[i] = (char*)NCI_Malloc(rsize[i]);

        // Post irecv
        MPI_Imrecv(rbuf[i], rsize[i], MPI_BYTE, &rmsg, rreq + i);
    }

    tbuf = (char*)NCI_Malloc(varp->chunksize);

    // Handle our own data
    for(req = 0; req < nreq; req++){
        nczipioi_chunk_itr_init_cord(varp, starts[req], counts[req], citr); // Initialize chunk iterator
        do{
            // Chunk index and owner
            cid = get_chunk_idx_cord(varp, citr);

            if (varp->chunk_owner[cid] == nczipp->rank){
                // Get overlap region
                get_chunk_overlap_cord(varp, citr, starts[req], counts[req], ostart, osize);
                overlapsize = varp->esize;
                for(i = 0; i < varp->ndim; i++){
                    overlapsize *= osize[i];                     
                }

                if (overlapsize > 0){
                    // Pack type from chunk cache to (contiguous) intermediate buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - citr[j]);
                        tsize[j] = varp->chunkdim[j];
                        tssize[j] = (int)osize[j];
                    }
                    //printf("Rank: %d, MPI_Type_create_subarray_self([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    MPI_Type_commit(&ptype);

                    // Pack data into intermediate buffer
                    packoff = 0;
                    MPI_Pack(varp->chunk_cache[cid], 1, ptype, tbuf, varp->chunksize, &packoff, nczipp->comm);
                    MPI_Type_free(&ptype);
                    overlapsize = packoff;

                    // Pack type from (contiguous) intermediate buffer to chunk buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - starts[req][j]);
                        tsize[j] = (int)counts[req][j];
                    }
                    //printf("Rank: %d, MPI_Type_create_subarray_self2([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    MPI_Type_commit(&ptype);
                    
                    // Unpack data into chunk buffer
                    packoff = 0;
                    //printf("Rank: %d, cid = %d, MPI_Unpack_self(%d, %d)\n", nczipp->rank, cid, overlapsize, packoff); fflush(stdout);
                    MPI_Unpack(tbuf, overlapsize, &packoff, bufs[req], 1, ptype, nczipp->comm);
                    MPI_Type_free(&ptype);    
                }
            }
        } while (nczipioi_chunk_itr_next_cord(varp, starts[req], counts[req], citr));
    }

    //Handle incoming requests
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = varp->chunkdim[i];
    }
    for(i = 0; i < nrecv; i++){
        // Will wait any provide any benefit?
        MPI_Waitany(nrecv, rreq, &j, &rstat);
        packoff = 0;
        //printf("rsize_2 = %d\n", rsizes[j]); fflush(stdout);
        //printf("Rank: %d, MPI_Unpack_bufsize(%d, %d, %d)\n", nczipp->rank, j, rsize[j], packoff); fflush(stdout);
        MPI_Unpack(rbuf[j], rsize[j], &packoff, ssize_re + j, 1, MPI_INT, nczipp->comm);
        sbuf_re[j] = (char*)NCI_Malloc(ssize_re[j]);
        while(packoff < rsize[j]){
            // Retrieve metadata
            MPI_Unpack(rbuf[j], rsize[j], &packoff, &cid, 1, MPI_INT, nczipp->comm);
            MPI_Unpack(rbuf[j], rsize[j], &packoff, tstart, varp->ndim, MPI_INT, nczipp->comm);
            MPI_Unpack(rbuf[j], rsize[j], &packoff, tssize, varp->ndim, MPI_INT, nczipp->comm);

            // Pack type
            //printf("Rank: %d, cid = %d, MPI_Type_create_subarray_rep([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, cid, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            MPI_Type_commit(&ptype);

            // Pack data
            MPI_Pack(varp->chunk_cache[cid], 1, ptype, sbuf_re[j], ssize_re[j], soff_re + j, nczipp->comm);
            //printf("cache[0] = %d, cache[1] = %d\n", ((int*)(varp->chunk_cache[cid]))[0], ((int*)(varp->chunk_cache[cid]))[1]); fflush(stdout);
            MPI_Type_free(&ptype);
        }
        // Free the request
        //MPI_Request_free(rreq + j);

        // Send Response
        //printf("Rank: %d, MPI_Isend(%d, %d, %d, %d)\n", nczipp->rank, soff_re[j], rstat.MPI_SOURCE, 1, j); fflush(stdout);
        MPI_Isend(sbuf_re[j], soff_re[j], MPI_BYTE, rstat.MPI_SOURCE, 1, nczipp->comm, sreq_re + j);
    }

    // Wait for all request
    MPI_Waitall(nsend, sreq, sstat);

    //Handle reply
    memset(soff, 0, sizeof(int) * nsend);
    memset(rcnt_local, 0, sizeof(int) * nsend);
    for(i = 0; i < nsend; i++){
        // Will wait any provide any benefit?
        MPI_Waitany(nsend, rreq_re, &j, &rstat);
        soff[j] = sizeof(int);  // Skip reply size
        packoff = 0;
        while(packoff < rsize_re[j]){
            // Retrieve metadata from the request we sent
            //printf("Rank: %d, MPI_Unpack_meta(%d, %d, %d)\n", nczipp->rank, j, ssize[j], soff[j]); fflush(stdout);
            MPI_Unpack(sbuf[j], ssize[j], soff + j, &cid, 1, MPI_INT, nczipp->comm);
            MPI_Unpack(sbuf[j], ssize[j], soff + j, tstart, varp->ndim, MPI_INT, nczipp->comm);
            MPI_Unpack(sbuf[j], ssize[j], soff + j, tssize, varp->ndim, MPI_INT, nczipp->comm);
            req = reqs[j][rcnt_local[j]++];
            get_chunk_itr(varp, cid, citr);
            for(k = 0; k < varp->ndim; k++){
                tstart[k] += (int)(citr[k] - starts[req][k]);
                tsize[k] = counts[req][k];
            }

            // Pack type
            //printf("Rank: %d, cid = %d, MPI_Type_create_subarray_resp([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, cid, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            MPI_Type_commit(&ptype);

            // Pack data
            //printf("Rank: %d, cid = %d, MPI_Unpack(%d, %d, %d, %d)\n", nczipp->rank, cid, j, rsize_re[j], packoff, req); fflush(stdout);
            MPI_Unpack(rbuf_re[j], rsize_re[j], &packoff, bufs[req], 1, ptype, nczipp->comm);
            //printf("cache[0] = %d, cache[1] = %d\n", ((int*)(varp->chunk_cache[cid]))[0], ((int*)(varp->chunk_cache[cid]))[1]); fflush(stdout);
            MPI_Type_free(&ptype);
        }
    }

    // Wait for all Response
    MPI_Waitall(nrecv, sreq_re, sstat_re);

    // Free buffers
    NCI_Free(rcnt_local);
    NCI_Free(rcnt_all);
    NCI_Free(smap);

    NCI_Free(rids);

    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);
    NCI_Free(osize);
    NCI_Free(ostart);
    NCI_Free(citr);

    NCI_Free(sreq);
    NCI_Free(sstat);
    NCI_Free(ssize);
    NCI_Free(sdst);
    NCI_Free(soff);
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

    if (tbuf != NULL){
        NCI_Free(tbuf);
    }

    return NC_NOERR;
}

int
nczipioi_get_varn(NC_zip        *nczipp,
              NC_zip_var       *varp,
              int              nreq,
              MPI_Offset* const *starts,
              MPI_Offset* const *counts,
              const void       *buf)
{
    int i, j;
    MPI_Offset rsize;
    char *bptr = (char*)buf;
    char **bufs;
    
    // Calculate buffer offset of each request
    bufs = (char**)NCI_Malloc(sizeof(char*) * nreq);
    for(i = 0; i < nreq; i++){
        bufs[i] = bptr;
        rsize = varp->esize;
        for(j = 0; j < varp->ndim; j++){
            rsize *= counts[i][j];
        }
        bptr += rsize;
    }

    // Collective buffer
    switch (nczipp->comm_unit){
        case NC_ZIP_COMM_CHUNK:
            nczipioi_get_varn_cb_chunk(nczipp, varp, nreq, starts, counts, NULL, bufs);
            break;
        case NC_ZIP_COMM_PROC:
            nczipioi_get_varn_cb_proc(nczipp, varp, nreq, starts, counts, bufs);
            break;
    }
    NCI_Free(bufs);

    return NC_NOERR;
}