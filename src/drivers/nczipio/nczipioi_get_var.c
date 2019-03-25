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
nczipioi_get_var(NC_zip          *nczipp,
                    NC_zip_var      *varp,
                    MPI_Offset      *start,
                    MPI_Offset      *count,
                    MPI_Offset      *stride,
                    void            *buf)
{
    int err;
    int i, j, k;
    int cid;   // Chunk iterator

    MPI_Offset *ostart, *osize;
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int *cstart, *cend, *citr; // Bounding box for chunks overlapping my own write region
    
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
    memset(rcnt_local, 0, sizeof(int) * nczipp->np);
    nsend = 0;

    // Initialize chunk iterator
    nczipioi_chunk_itr_init(varp, start, count, cstart, cend, citr);

    // Iterate through chunks
    do{
        // Chunk index
        cid = get_chunk_idx(varp, citr);
        
        rcnt_local[cid] = 1;

        if (varp->chunk_owner[cid] != nczipp->rank){
            // Count number of mnessage we need to send
            nsend++;    
        }

    } while (nczipioi_chunk_itr_next(varp, start, count, cstart, cend, citr));


    // Sync number of messages of each chunk
    MPI_Allreduce(rcnt_local, rcnt_all, nczipp->np, MPI_INT, MPI_SUM, nczipp->comm);

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
    k = 0;
    // Initialize chunk iterator
    nczipioi_chunk_itr_init(varp, start, count, cstart, cend, citr);
    // Iterate through chunks
    do{
        // Chunk index
        cid = get_chunk_idx(varp, citr);

        // We got something to send if we are not owner
        if (varp->chunk_owner[cid] != nczipp->rank){
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
                tstart[j] = (int)(ostart[j] - citr[j] * varp->chunkdim[j]);
                tsize[j] = (int)osize[j];
            }
                    
            // Pack metadata
            packoff = 0;
            MPI_Pack(tstart, varp->ndim, MPI_INT, sbufs[k], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
            MPI_Pack(tsize, varp->ndim, MPI_INT, sbufs[k], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);

            // Send and receive
            //printf("packoff = %d\n", packoff); fflush(stdout);
            //printf("Rank: %d, MPI_Isend(%d, %d, %d, %d)\n", nczipp->rank, packoff, varp->chunk_owner[cid], cid, k); fflush(stdout);
            MPI_Isend(sbufs[k], packoff, MPI_BYTE, varp->chunk_owner[cid], cid, nczipp->comm, sreqs + k);
            //printf("Rank: %d, MPI_Irecv(%d, %d, %d, %d)\n", nczipp->rank, overlapsize, varp->chunk_owner[cid], cid + 1024, nrecv + k); fflush(stdout);
            err = MPI_Irecv(rbufs[k + nrecv], overlapsize, MPI_BYTE, varp->chunk_owner[cid], cid + 1024, nczipp->comm, rreqs + nrecv + k);
            if (err != 0){
                printf("err = %d\n", err); fflush(stdout);
            }
            k++;
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, cstart, cend, citr));

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


    // Allocate intermediate buffer
    cbuf = (char*)NCI_Malloc(varp->chunksize);

    // For each chunk we own, we need to receive incoming data
    k = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];
            
        // Handle our own data first if we have any
        if (rcnt_local[cid] > 0){
            // Convert chunk id to iterator
            get_chunk_cord(varp, cid, citr);

            // Calculate overlapping region
            get_chunk_overlap(varp, citr, start, count, ostart, osize);
            
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
            // Get metadata
            packoff = 0;
            MPI_Unpack(rbufs[j], rsizes[j], &packoff, tstart, varp->ndim, MPI_INT, nczipp->comm);
            MPI_Unpack(rbufs[j], rsizes[j], &packoff, tssize, varp->ndim, MPI_INT, nczipp->comm);

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

            // Send reply
            //printf("Rank: %d, MPI_Isend(%d, %d, %d, %d)\n", nczipp->rank, packoff, varp->chunk_owner[cid], cid + 1024, k + nsend); fflush(stdout);
            MPI_Isend(sbufs[j + nsend], packoff, MPI_BYTE, rstats[j].MPI_SOURCE, cid + 1024, nczipp->comm, sreqs + j + nsend);
        }
        k += rcnt_all[cid] - rcnt_local[cid];        

        //princbuf(nczipp->rank, varp->chunk_cache[cid], varp->chunksize);
    }

    // Wait for all send request
    //printf("Rank: %d, MPI_Waitall_send(%d, %d)\n", nczipp->rank, nsend, 0); fflush(stdout);
    MPI_Waitall(nsend, sreqs, sstats);

    // Receive replies from the owners and update the user buffer
    k = 0;
    // Initialize chunk iterator
    nczipioi_chunk_itr_init(varp, start, count, cstart, cend, citr);
    // Iterate through chunks
    do{
        // Chunk index
        cid = get_chunk_idx(varp, citr);

        // We got something to recv if we are not owner
        if (varp->chunk_owner[cid] != nczipp->rank){
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

            //printf("Rank: %d, wait recv, nrecv = %d, k = %d, nsend = %d\n", nczipp->rank, nrecv, k, nsend); fflush(stdout);
            // Wait for reply
            //printf("Rank: %d, MPI_Wait_recv(%d)\n", nczipp->rank, nrecv + k); fflush(stdout);
            MPI_Wait(rreqs + nrecv + k, rstats + nrecv + k);

            // Pack data
            //printf("Rank: %d, pack\n", nczipp->rank); fflush(stdout);
            packoff = 0;
            MPI_Unpack(rbufs[nrecv + k], overlapsize, &packoff, buf, 1, ptype, nczipp->comm);
            MPI_Type_free(&ptype);

            k++;
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, cstart, cend, citr));

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