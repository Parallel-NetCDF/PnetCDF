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
nczipioi_get_var_cb_chunk(NC_zip          *nczipp,
                    NC_zip_var      *varp,
                    const MPI_Offset      *start,
                    const MPI_Offset      *count,
                    const MPI_Offset      *stride,
                    void            *buf)
{
    int err;
    int i, j, k;
    int cid;   // Chunk iterator

    MPI_Offset *ostart, *osize;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    MPI_Offset *citr; // Chunk iterator
    
    int *rcnt_local, *rcnt_all;   // Number of processes that writes to each chunk

    int overlapsize;    // Size of overlaping region of request and chunk
    char *cbuf = NULL;     // Intermediate continuous buffer
    
    int packoff; // Pack offset
    MPI_Datatype ptype; // Pack datatype

    int nread;  // # chunks to read form file
    int *rids;  // Id of chunks to read from file

    int nsend, nrecv;   // Number of send and receive
    MPI_Request *sreqs, *rreqs;    // Send and recv req
    MPI_Status *sstats, *rstats;    // Send and recv status
    char **sbufs, **rbufs;   // Send and recv buffer
    int *rsizes;    // recv size of each message
    MPI_Message rmsg;   // Receive message

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_INIT)

    // Allocate buffering for write count
    rcnt_local = (int*)NCI_Malloc(sizeof(int) * varp->nchunk * 2);
    rcnt_all = rcnt_local + varp->nchunk;

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim * 3);
    tssize = tsize + varp->ndim;
    tstart = tssize + varp->ndim;
    ostart = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim * 3);
    osize = ostart + varp->ndim;

    // Chunk iterator
    citr = osize + varp->ndim;

    // We need to calculate the size of message of each chunk
    // This is just for allocating send buffer
    // We do so by iterating through all request and all chunks they cover
    // If we are not the owner of a chunk, we need to send message
    memset(rcnt_local, 0, sizeof(int) * varp->nchunk);
    nsend = 0;

    // Iterate through chunks
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid);
    do{
        rcnt_local[cid] = 1;

        if (varp->chunk_owner[cid] != nczipp->rank){
            // Count number of mnessage we need to send
            nsend++;    
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_INIT)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SYNC)

    // Sync number of messages of each chunk
    MPI_Allreduce(rcnt_local, rcnt_all, varp->nchunk, MPI_INT, MPI_SUM, nczipp->comm);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SYNC)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_INIT)

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

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB)  // I/O time count separately

    // Decompress chunks into chunk cache
    nczipioi_load_var(nczipp, varp, nread, rids);

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB)

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
    k = 0;
    // Initialize chunk iterator
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid);
    // Iterate through chunks
    do{
        // We got something to send if we are not owner
        if (varp->chunk_owner[cid] != nczipp->rank){
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_PACK_REQ)

            // Calculate chunk overlap
            get_chunk_overlap(varp, citr, start, count, ostart, osize);
            //printf("cord = %d, start = %lld, count = %lld, tstart = %d, tssize = %d, esize = %d, ndim = %d\n", citr[0], starts[req][0], counts[req][0], tstart[0], tssize[0], varp->esize, varp->ndim); fflush(stdout);
            overlapsize = varp->esize;
            for(j = 0; j < varp->ndim; j++){
                overlapsize *= osize[j];                     
            }
            //printf("overlapsize = %d\n", overlapsize); fflush(stdout);

            // Allocate buffer
            sbufs[k] = (char*)NCI_Malloc(sizeof(int) * varp->ndim * 2); // For request
            rbufs[k + nrecv] = (char*)NCI_Malloc(overlapsize);   // For reply, first nrecv are for request

            // Metadata
            for(j = 0; j < varp->ndim; j++){
                tstart[j] = (int)(ostart[j] - citr[j]);
                tsize[j] = (int)osize[j];
            }
                    
            // Pack metadata
            packoff = 0;
            MPI_Pack(tstart, varp->ndim, MPI_INT, sbufs[k], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
            MPI_Pack(tsize, varp->ndim, MPI_INT, sbufs[k], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);

            // Send and receive
            //printf("packoff = %d\n", packoff); fflush(stdout);
            //printf("Rank: %d, MPI_Isend(%d, %d, %d, %d)\n", nczipp->rank, packoff, varp->chunk_owner[cid], cid, k); fflush(stdout);
            
            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_PACK_REQ)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SEND_REQ)

            MPI_Isend(sbufs[k], packoff, MPI_BYTE, varp->chunk_owner[cid], cid, nczipp->comm, sreqs + k);
            
            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REQ)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_RECV_REP)

            //printf("Rank: %d, MPI_Irecv(%d, %d, %d, %d)\n", nczipp->rank, overlapsize, varp->chunk_owner[cid], cid + 1024, nrecv + k); fflush(stdout);
            err = MPI_Irecv(rbufs[k + nrecv], overlapsize, MPI_BYTE, varp->chunk_owner[cid], cid + 1024, nczipp->comm, rreqs + nrecv + k);
            if (err != 0){
                printf("err = %d\n", err); fflush(stdout);
            }

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_RECV_REP)

            k++;
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_RECV_REQ)

    // Post recv
    k = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];
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

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_RECV_REQ)

    // Allocate intermediate buffer
    cbuf = (char*)NCI_Malloc(varp->chunksize);

    // For each chunk we own, we need to receive incoming data
    k = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];

        NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SELF)

        // Handle our own data first if we have any
        if (rcnt_local[cid] > 0){
            // Convert chunk id to iterator
            get_chunk_itr(varp, cid, citr);

            // Calculate overlapping region
            get_chunk_overlap(varp, citr, start, count, ostart, osize);
            
            // Pack type from chunk buffer to (contiguous) intermediate buffer
            for(j = 0; j < varp->ndim; j++){
                tstart[j] = (int)(ostart[j] - citr[j]);
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
                tstart[j] = (int)(ostart[j] - start[j]);
                tsize[j] = (int)count[j];
            }
            //printf("Rank: %d, MPI_Type_create_subarray([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            MPI_Type_commit(&ptype);

            // Pack data into user buffer
            packoff = 0;
            MPI_Unpack(cbuf, overlapsize, &packoff, buf, 1, ptype, nczipp->comm);
            MPI_Type_free(&ptype);
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SELF)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_RECV_REQ)

        // Wait for all send requests related to this chunk
        // We remove the impact of -1 mark in rcnt_local[cid]
        //printf("Rank: %d, MPI_Waitall_recv(%d, %d)\n", nczipp->rank, rcnt_all[cid] - rcnt_local[cid], k); fflush(stdout);
        MPI_Waitall(rcnt_all[cid] - rcnt_local[cid], rreqs + k, rstats + k);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_RECV_REQ)

        // Now, it is time to process data from other processes
        for(j = 0; j < varp->ndim; j++){
            tsize[j] = varp->chunkdim[j];
        }

        // Process data received
        //printf("nrecv = %d, rcnt_all = %d, rcnt_local = %d\n", nrecv, rcnt_all[cid], rcnt_local[cid]); fflush(stdout);
        for(j = k; j < k + rcnt_all[cid] - rcnt_local[cid]; j++){
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_UNPACK_REQ)
            
            // Get metadata
            packoff = 0;
            MPI_Unpack(rbufs[j], rsizes[j], &packoff, tstart, varp->ndim, MPI_INT, nczipp->comm);
            MPI_Unpack(rbufs[j], rsizes[j], &packoff, tssize, varp->ndim, MPI_INT, nczipp->comm);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_UNPACK_REQ)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_PACK_REP)

            // Pack type
            //printf("Rank: %d, MPI_Type_create_subarray([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            MPI_Type_commit(&ptype);

            // Allocate buffer 
            MPI_Type_size(ptype, &overlapsize);
            sbufs[j + nsend] = (char*)NCI_Malloc(overlapsize); // For reply
            
            // Pack data
            packoff = 0;
            MPI_Pack(varp->chunk_cache[cid], 1, ptype, sbufs[j + nsend], overlapsize, &packoff, nczipp->comm);
            MPI_Type_free(&ptype);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_PACK_REP)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SEND_REP)

            // Send reply
            //printf("Rank: %d, MPI_Isend(%d, %d, %d, %d)\n", nczipp->rank, packoff, varp->chunk_owner[cid], cid + 1024, k + nsend); fflush(stdout);
            MPI_Isend(sbufs[j + nsend], packoff, MPI_BYTE, rstats[j].MPI_SOURCE, cid + 1024, nczipp->comm, sreqs + j + nsend);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REP)
        }
        k += rcnt_all[cid] - rcnt_local[cid];        

        //princbuf(nczipp->rank, varp->chunk_cache[cid], varp->chunksize);
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SEND_REQ)

    // Wait for all send request
    //printf("Rank: %d, MPI_Waitall_send(%d, %d)\n", nczipp->rank, nsend, 0); fflush(stdout);
    MPI_Waitall(nsend, sreqs, sstats);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REQ)

    // Receive replies from the owners and update the user buffer
    k = 0;
    // Initialize chunk iterator
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid);
    // Iterate through chunks
    do{
        // We got something to recv if we are not owner
        if (varp->chunk_owner[cid] != nczipp->rank){
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_UNPACK_REP)
            
            // Calculate chunk overlap
            get_chunk_overlap(varp, citr, start, count, ostart, osize);

            // Pack type from recv buffer to user buffer
            for(j = 0; j < varp->ndim; j++){
                tstart[j] = (int)(ostart[j] - start[j]);
                tsize[j] = (int)count[j];
                tssize[j] = (int)osize[j];
            }
            //printf("Rank: %d, ostart=[%lld, %lld], osize=[%lld, %lld]\n", nczipp->rank, ostart[0], ostart[1], osize[0], osize[1]); fflush(stdout);
            //printf("Rank: %d, MPI_Type_create_subarray4([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            //printf("Rank: %d, commit\n", nczipp->rank); fflush(stdout);
            MPI_Type_commit(&ptype);
            MPI_Type_size(ptype, &overlapsize);

            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_RECV_REP)

            //printf("Rank: %d, wait recv, nrecv = %d, k = %d, nsend = %d\n", nczipp->rank, nrecv, k, nsend); fflush(stdout);
            // Wait for reply
            //printf("Rank: %d, MPI_Wait_recv(%d)\n", nczipp->rank, nrecv + k); fflush(stdout);
            MPI_Wait(rreqs + nrecv + k, rstats + nrecv + k);

            NC_ZIP_TIMER_STOPEX(NC_ZIP_TIMER_CB_RECV_REP, NC_ZIP_TIMER_CB_UNPACK_REP)

            // Pack data
            //printf("Rank: %d, pack\n", nczipp->rank); fflush(stdout);
            packoff = 0;
            MPI_Unpack(rbufs[nrecv + k], overlapsize, &packoff, buf, 1, ptype, nczipp->comm);
            MPI_Type_free(&ptype);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_UNPACK_REP)

            k++;
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SEND_REP)

    //printf("Rank: %d, wait_final\n", nczipp->rank); fflush(stdout);
    // Wait for all send replies
    //printf("Rank: %d, MPI_Waitall_send(%d, %d)\n", nczipp->rank, nrecv, nsend); fflush(stdout);
    MPI_Waitall(nrecv, sreqs + nsend, sstats + nsend);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REP)

    //printf("Rank: %d, exiting\n", nczipp->rank); fflush(stdout);

    // Free buffers
    NCI_Free(rcnt_local);

    NCI_Free(rids);

    NCI_Free(tsize);

    NCI_Free(ostart);

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

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REP)

    return NC_NOERR;
}


int
nczipioi_get_var_cb_proc(      NC_zip          *nczipp,
                            NC_zip_var      *varp,
                            const MPI_Offset      *start,
                            const MPI_Offset      *count,
                            const MPI_Offset      *stride,
                            void            *buf)
{
    int err;
    int i, j, k;
    int cid, cown;   // Chunk iterator

    MPI_Offset *ostart, *osize;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    MPI_Offset *citr; // Chunk iterator
    
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

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_INIT)

    // Allocate buffering for write count
    rcnt_local = (int*)NCI_Malloc(sizeof(int) * (nczipp->np * 3 + varp->nchunk * 2));
    rcnt_local_chunk = rcnt_local + nczipp->np;
    rcnt_all = rcnt_local_chunk + varp->nchunk;
    rcnt_all_chunk = rcnt_all + nczipp->np;
    smap = rcnt_all_chunk + varp->nchunk;

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim * 3);
    tssize = tsize + varp->ndim;
    tstart = tssize + varp->ndim;
    ostart = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim * 3);
    osize = ostart + varp->ndim;

    // Chunk iterator
    citr = osize + varp->ndim;

    // We need to calculate the size of message of each chunk
    // This is just for allocating send buffer
    // We do so by iterating through all request and all chunks they cover
    // If we are not the owner of a chunk, we need to send message
    memset(rcnt_local, 0, sizeof(int) * (nczipp->np + varp->nchunk));
    nsend = 0;

    // Count total number of messages and build a map of accessed chunk to list of comm datastructure
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid); // Initialize chunk iterator
    do{
        // Chunk owner
        cown = varp->chunk_owner[cid];

        // Mapping to skip list of send requests 
        if (rcnt_local[cown] == 0 && cown != nczipp->rank){
            smap[cown] = nsend++;
        }
        rcnt_local[cown] = 1;   // Need to send message if not owner     
        rcnt_local_chunk[cid] = 1;  // This tells the owner to prepare the chunks  
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_INIT)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SYNC)

    // Sync number of messages of each chunk
    MPI_Allreduce(rcnt_local, rcnt_all, nczipp->np + varp->nchunk, MPI_INT, MPI_SUM, nczipp->comm);
    nrecv = rcnt_all[nczipp->rank] - rcnt_local[nczipp->rank];  // We don't need to receive request form self

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SYNC)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_INIT)

    // We need to prepare chunk in the chunk cache
    // For chunks not yet allocated, we need to read them form file collectively
    // We collect chunk id of those chunks
    // Calculate number of recv request
    // This is for all the chunks
    rids = (int*)NCI_Malloc(sizeof(int) * varp->nmychunks);
    nread = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];
        // Count number of chunks we need to prepare
        // We read only chunks that is required
        if (rcnt_all_chunk[cid] > 0 && varp->chunk_cache[cid] == NULL){
            rids[nread++] = cid;
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB)  // I/O time count separately

    // Decompress chunks into chunk cache
    nczipioi_load_var(nczipp, varp, nread, rids);

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_PACK_REQ)

    // Allocate data structure for messaging
    sbuf = (char**)NCI_Malloc(sizeof(char*) * (nsend + nrecv));
    ssize = (int*)NCI_Malloc(sizeof(int) * (nsend * 3 + nrecv * 2));
    soff = ssize + (nsend + nrecv);
    sdst = soff + (nsend + nrecv);
    sreq = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * (nsend + nrecv));
    sstat = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * (nsend + nrecv));

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

    // Count size of each request
    memset(ssize, 0, sizeof(int) * nsend);
    memset(rsize_re, 0, sizeof(int) * nsend);
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid); // Initialize chunk iterator
    do{
        // Chunk owner               
        cown = varp->chunk_owner[cid];
        if (cown != nczipp->rank){
            j = smap[cown];
            sdst[j] = cown; // Record a reverse map by the way

            // Count overlap
            get_chunk_overlap(varp, citr, start, count, ostart, osize);
            overlapsize = varp->esize;
            for(i = 0; i < varp->ndim; i++){
                overlapsize *= osize[i];                     
            }
            ssize[j] += sizeof(int) * (varp->ndim * 2 + 1);
            rsize_re[j] += overlapsize;
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    // Allocate buffer for send
    memset(soff, 0, sizeof(int) * (nsend + nrecv));
    for(i = 0; i < nsend; i++){
        ssize[i] += sizeof(int);
        sbuf[i] = (char*)NCI_Malloc(ssize[i]);
        MPI_Pack(rsize_re + i, 1, MPI_INT, sbuf[i], ssize[i], soff + i, nczipp->comm);
        rbuf_re[i] = (char*)NCI_Malloc(rsize_re[i]);
    }

    // Pack requests
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid); // Initialize chunk iterator
    do{
        // Chunk owner
        cown = varp->chunk_owner[cid];
        if (cown != nczipp->rank){
            j = smap[cown];

            // Get overlap region
            get_chunk_overlap(varp, citr, start, count, ostart, osize);

            // Pack metadata
            for(i = 0; i < varp->ndim; i++){
                tstart[i] = (int)(ostart[i] - citr[i]);
                tssize[i] = (int)osize[i];
            }
            MPI_Pack(&cid, 1, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
            MPI_Pack(tstart, varp->ndim, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
            MPI_Pack(tssize, varp->ndim, MPI_INT, sbuf[j], ssize[j], soff + j, nczipp->comm);
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_PACK_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SEND_REQ)

    // Post send 
    for(i = 0; i < nsend; i++){
        MPI_Isend(sbuf[i], soff[i], MPI_BYTE, sdst[i], 0, nczipp->comm, sreq + i);
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_RECV_REP)

    // Post receive  
    for(i = 0; i < nsend; i++){
        MPI_Irecv(rbuf_re[i], rsize_re[i], MPI_BYTE, sdst[i], 1, nczipp->comm, rreq_re + i);
    }   

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_RECV_REP)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_RECV_REQ)

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

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_RECV_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SELF)

    tbuf = (char*)NCI_Malloc(varp->chunksize);

    // Handle our own data
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid); // Initialize chunk iterator
    do{
        if (varp->chunk_owner[cid] == nczipp->rank){
            // Get overlap region
            overlapsize = get_chunk_overlap(varp, citr, start, count, ostart, osize);

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
                    tstart[j] = (int)(ostart[j] - start[j]);
                    tsize[j] = (int)count[j];
                }
                //printf("Rank: %d, MPI_Type_create_subarray_self2([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                MPI_Type_commit(&ptype);
                
                // Unpack data into chunk buffer
                packoff = 0;
                MPI_Unpack(tbuf, overlapsize, &packoff, buf, 1, ptype, nczipp->comm);
                MPI_Type_free(&ptype);    
            }
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SELF)

    //Handle incoming requests
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = varp->chunkdim[i];
    }
    for(i = 0; i < nrecv; i++){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_RECV_REQ)

        // Will wait any provide any benefit?
        MPI_Waitany(nrecv, rreq, &j, &rstat);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_RECV_REQ)

        packoff = 0;
        //printf("rsize_2 = %d\n", rsizes[j]); fflush(stdout);
        MPI_Unpack(rbuf[j], rsize[j], &packoff, ssize_re + j, 1, MPI_INT, nczipp->comm);
        sbuf_re[j] = (char*)NCI_Malloc(ssize_re[j]);
        while(packoff < rsize[j]){
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_UNPACK_REQ)

            // Retrieve metadata
            MPI_Unpack(rbuf[j], rsize[j], &packoff, &cid, 1, MPI_INT, nczipp->comm);
            MPI_Unpack(rbuf[j], rsize[j], &packoff, tstart, varp->ndim, MPI_INT, nczipp->comm);
            MPI_Unpack(rbuf[j], rsize[j], &packoff, tssize, varp->ndim, MPI_INT, nczipp->comm);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_UNPACK_REQ)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_PACK_REP)

            // Pack type
            //printf("Rank: %d, cid = %d, MPI_Type_create_subarray_rep([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, cid, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            MPI_Type_commit(&ptype);

            // Pack data
            MPI_Pack(varp->chunk_cache[cid], 1, ptype, sbuf_re[j], ssize_re[j], soff_re + j, nczipp->comm);
            //printf("cache[0] = %d, cache[1] = %d\n", ((int*)(varp->chunk_cache[cid]))[0], ((int*)(varp->chunk_cache[cid]))[1]); fflush(stdout);
            MPI_Type_free(&ptype);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_PACK_REP)
        }

        NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SEND_REQ)

        // Send Response
        //printf("Rank: %d, MPI_Isend(%d, %d, %d, %d)\n", nczipp->rank, soff_re[j], rstat.MPI_SOURCE, 1, j); fflush(stdout);
        MPI_Isend(sbuf_re[j], soff_re[j], MPI_BYTE, rstat.MPI_SOURCE, 1, nczipp->comm, sreq_re + j);
        
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REQ)
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SEND_REQ)

    // Wait for all request
    MPI_Waitall(nsend, sreq, sstat);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REQ)

    //Handle reply
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = count[i];
    }
    memset(soff, 0, sizeof(int) * nsend);
    for(i = 0; i < nsend; i++){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_RECV_REP)

        // Will wait any provide any benefit?
        MPI_Waitany(nsend, rreq_re, &j, &rstat);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_RECV_REP)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_UNPACK_REP)

        soff[j] = sizeof(int);  // Skip reply size
        packoff = 0;
        while(packoff < rsize_re[j]){
            // Retrieve metadata from the request we sent
            MPI_Unpack(sbuf[j], ssize[j], soff + j, &cid, 1, MPI_INT, nczipp->comm);
            MPI_Unpack(sbuf[j], ssize[j], soff + j, tstart, varp->ndim, MPI_INT, nczipp->comm);
            MPI_Unpack(sbuf[j], ssize[j], soff + j, tssize, varp->ndim, MPI_INT, nczipp->comm);
            get_chunk_itr(varp, cid, citr);
            for(k = 0; k < varp->ndim; k++){
                tstart[k] += (int)(citr[k] - start[k]);
            }

            // Pack type
            //printf("Rank: %d, cid = %d, MPI_Type_create_subarray_resp([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, cid, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
            MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            MPI_Type_commit(&ptype);

            // Pack data
            MPI_Unpack(rbuf_re[j], rsize_re[j], &packoff, buf, 1, ptype, nczipp->comm);
            //printf("cache[0] = %d, cache[1] = %d\n", ((int*)(varp->chunk_cache[cid]))[0], ((int*)(varp->chunk_cache[cid]))[1]); fflush(stdout);
            MPI_Type_free(&ptype);
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_UNPACK_REP)
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_CB_SEND_REP)

    // Wait for all Response
    MPI_Waitall(nrecv, sreq_re, sstat_re);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB_SEND_REP)

    // Free buffers
    NCI_Free(rcnt_local);


    NCI_Free(rids);

    NCI_Free(tsize);

    NCI_Free(ostart);

    NCI_Free(sreq);
    NCI_Free(sstat);
    NCI_Free(ssize);
    for(i = 0; i < nsend + nrecv; i++){
        NCI_Free(sbuf[i]);
        NCI_Free(rbuf[i]);
    }
    NCI_Free(sbuf);

    NCI_Free(rreq);
    NCI_Free(rbuf);
    NCI_Free(rsize);

    if (tbuf != NULL){
        NCI_Free(tbuf);
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_CB)

    return NC_NOERR;
}
