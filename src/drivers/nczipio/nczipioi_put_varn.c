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

int printbuf(int rank, char* buf, int len){
    int i;
    printf("Rank %d: ", rank);
    for(i = 0; i < len; i++){
        printf("%02x ", buf[i]);
    }
    printf("\n");
}

int
nczipioi_put_varn_cb(  NC_zip        *nczipp,
                    NC_zip_var       *varp,
                    int              nreq,
                    MPI_Offset* const *starts,
                    MPI_Offset* const *counts,
                    MPI_Offset* const *strides,
                    void              **bufs)
{
    int err;
    int i, j, k;
    int cid, req;   // Chunk and request iterator

    int *tsize, *tssize, *tstart;   // Size for sub-array type
    MPI_Offset *ostart, *osize;
    int *cstart, *cend, *citr; // Bounding box for chunks overlapping my own write region
    
    int *wcnt_local, *wcnt_all;   // Number of processes that writes to each chunk

    int overlapsize;    // Size of overlaping region of request and chunk
    int max_tbuf;   // Size of intermediate buffer
    char *tbuf = NULL;     // Intermediate buffer
    
    int packoff; // Pack offset
    MPI_Datatype ptype; // Pack datatype

    int nsend, nrecv;   // Number of send and receive
    MPI_Request *sreqs, *rreqs;    // Send and recv req
    MPI_Status *sstats, *rstats;    // Send and recv status
    char **sbufs, **rbufs;   // Send and recv buffer
    int *rsizes;    // recv size of each message
    MPI_Message rmsg;   // Receive message

    // Allocate buffering for write count
    wcnt_local = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    wcnt_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);

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
    memset(wcnt_local, 0, sizeof(int) * nczipp->np);
    nsend = 0;
    max_tbuf = 0;
    for(req = 0; req < nreq; req++){
        // Initialize chunk iterator
        nczipioi_chunk_itr_init(varp, starts[req], counts[req], cstart, cend, citr);

        // Iterate through chunks
        do{
            // Chunk index
            cid = get_chunk_idx(varp, citr);

            // Calculate overlapping
            get_chunk_overlap(varp, citr, starts[req], counts[req], ostart, osize);
            //printf("cord = %d, start = %lld, count = %lld, tstart = %d, tssize = %d, esize = %d, ndim = %d\n", citr[0], starts[i][0], counts[i][0], tstart[0], tssize[0], esize, varp->ndim); fflush(stdout);
            overlapsize = varp->esize;
            for(j = 0; j < varp->ndim; j++){
                overlapsize *= osize[j];                     
            }

            if (varp->chunk_owner[cid] != nczipp->rank){
                // Count number of mnessage we need to send
                if (wcnt_local[cid] == 0){
                    nsend++;
                }
                wcnt_local[cid] += overlapsize + sizeof(int) * 2 * varp->ndim;
            }
            else{
                // We mark covered chunk of our own to prevent unnecessary calculation of overlap
                // -1 is purely a mark, we need to add 1 back to global message count
                wcnt_local[cid] = -1;

                // Record max overlapsize so we know how large the intermediate buffer is needed later
                if (max_tbuf < overlapsize){
                    max_tbuf = overlapsize;
                }
            }

        } while (nczipioi_chunk_itr_next(varp, starts[req], counts[req], cstart, cend, citr));
    }

    // Allocate buffer for sending
    sbufs = (char**)NCI_Malloc(sizeof(char*) * nsend);
    sreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nsend);
    sstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nsend);
    j = 0;
    // Allocate buffer for data
    for(cid = 0; cid < varp->nchunks; cid++){
        //printf("Rank: %d, wcnt_local[%d] = %d\n", nczipp->rank, cid, wcnt_local[cid]); fflush(stdout);
        // Count number of mnessage we need to send
        if (wcnt_local[cid] > 0){
            // Add space for number of reqs
            sbufs[j++] = (char*)NCI_Malloc(wcnt_local[cid]);
            // We don't need message size anymore, wcnt_local is used to track number of message from now on 
            wcnt_local[cid] = 1;
        }
    }

    // Sync number of messages of each chunk
    MPI_Allreduce(wcnt_local, wcnt_all, varp->nchunks, MPI_INT, MPI_SUM, nczipp->comm);

    // Calculate number of recv request
    // This is for all the chunks
    nrecv = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];
        // We don't need message for our own data
        nrecv += wcnt_all[cid] - wcnt_local[cid];
    }
    rreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nrecv);
    rstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nrecv);
    rbufs = (char**)NCI_Malloc(sizeof(char*) * nrecv);
    rsizes = (int*)NCI_Malloc(sizeof(int) * nrecv);

    // Post send and recv
    nrecv = 0;
    nsend = 0;
    for(cid = 0; cid < varp->nchunks; cid++){
        if (varp->chunk_owner[cid] == nczipp->rank){
            // We are the owner of the chunk
            // Receive data from other process
            for(i = 0; i < wcnt_all[cid] - wcnt_local[cid]; i++){
                // Get message size, including metadata
                MPI_Mprobe(MPI_ANY_SOURCE, cid, nczipp->comm, &rmsg, rstats);
                MPI_Get_count(rstats, MPI_BYTE, rsizes + nrecv);

                //printf("rsize = %d\n", rsizes[i]); fflush(stdout);

                // Allocate buffer
                rbufs[nrecv] = (char*)NCI_Malloc(rsizes[nrecv]);

                // Post irecv
                MPI_Imrecv(rbufs[nrecv], rsizes[nrecv], MPI_BYTE, &rmsg, rreqs + nrecv);
                nrecv++;
            }
        }
        else{
            // If we any of our request overlap with this chunk, we need to send data
            // We send only 1 message for 1 chunk
            if (wcnt_local[cid] > 0){
                packoff = 0;
                // Get chunk iterator
                get_chunk_cord(varp, cid, citr);  
                for(req = 0; req < nreq; req++){
                    // Calculate chunk overlap
                    get_chunk_overlap(varp, citr, starts[req], counts[req], ostart, osize);
                    //printf("cord = %d, start = %lld, count = %lld, tstart = %d, tssize = %d, esize = %d, ndim = %d\n", citr[0], starts[req][0], counts[req][0], tstart[0], tssize[0], varp->esize, varp->ndim); fflush(stdout);
                    overlapsize = varp->esize;
                    for(j = 0; j < varp->ndim; j++){
                        overlapsize *= osize[j];                     
                    }
                    //printf("overlapsize = %d\n", overlapsize); fflush(stdout);

                    // If current request have any overlap with the chunk, we pack the data and metadata
                    if (overlapsize > 0){
                        // Pack type
                        for(j = 0; j < varp->ndim; j++){
                            tstart[j] = (int)(ostart[j] - starts[req][j]);
                            tsize[j] = (int)counts[req][j];
                            tssize[j] = (int)osize[j];
                        }
                        //printf("Rank: %d, MPI_Type_create_subarray([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                        MPI_Type_commit(&ptype);
                        
                        // Pack data
                        //printf("Rank: %d, Pack([%d, %d], [%d, %d], ndim = %d\n", nczipp->rank, starts[req][0], starts[req][1], counts[req][0], counts[req][1], varp->ndim); fflush(stdout);

                        // Start and count related to chunk position
                        for(j = 0; j < varp->ndim; j++){
                            tstart[j] = (int)(ostart[j] - (MPI_Offset)citr[j] * (MPI_Offset)varp->chunkdim[j]);
                            tsize[j] = (int)osize[j];
                        }
                        MPI_Pack(tstart, varp->ndim, MPI_INT, sbufs[nsend], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        MPI_Pack(tsize, varp->ndim, MPI_INT, sbufs[nsend], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        MPI_Pack(bufs[req], 1, ptype, sbufs[nsend], packoff + overlapsize, &packoff, nczipp->comm);

                        MPI_Type_free(&ptype);
                    }
                }

                // Send the request
                //printf("packoff = %d\n", packoff); fflush(stdout);
                MPI_Isend(sbufs[nsend], packoff, MPI_BYTE, varp->chunk_owner[cid], cid, nczipp->comm, sreqs + nsend);
                nsend++;
            }
        }
    }

    // Allocate intermediate buffer
    if (max_tbuf > 0){
        tbuf = (char*)NCI_Malloc(max_tbuf);
    }

    // Wait for all send
    MPI_Waitall(nsend, sreqs, sstats);

    // For each chunk we own, we need to receive incoming data
    nrecv = 0;
    for(i = 0; i < varp->nmychunks; i++){
        cid = varp->mychunks[i];
            
        // Handle our own data first if we have any
        if (wcnt_local[cid] < 0){
            for(req = 0; req < nreq; req++){
                // Convert chunk id to iterator
                get_chunk_cord(varp, cid, citr);

                // Calculate overlapping region
                get_chunk_overlap(varp, citr, starts[req], counts[req], ostart, osize);
                overlapsize = varp->esize;
                for(j = 0; j < varp->ndim; j++){
                    overlapsize *= osize[j];
                }
                
                // If anything overlaps
                if (overlapsize > 0){
                    // Pack type from user buffer to (contiguous) intermediate buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - starts[req][j]);
                        tsize[j] = (int)counts[req][j];
                        tssize[j] = (int)osize[j];
                    }
                    //printf("Rank: %d, MPI_Type_create_subarray([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    MPI_Type_commit(&ptype);

                    // Pack data into intermediate buffer
                    packoff = 0;
                    MPI_Pack(bufs[req], 1, ptype, tbuf, overlapsize, &packoff, nczipp->comm);

                    MPI_Type_free(&ptype);

                    // Pack type from (contiguous) intermediate buffer to chunk buffer
                    for(j = 0; j < varp->ndim; j++){
                        tstart[j] = (int)(ostart[j] - (MPI_Offset)citr[j] * (MPI_Offset)varp->chunkdim[j]);
                        tsize[j] = varp->chunkdim[j];
                    }
                    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                    MPI_Type_commit(&ptype);
                    
                    // Unpack data into chunk buffer
                    packoff = 0;
                    MPI_Unpack(tbuf, overlapsize, &packoff, varp->chunk_cache[cid], 1, ptype, nczipp->comm);

                    MPI_Type_free(&ptype);
                }
            }
        }

        // Now, it is time to process data from other processes
        
        // Wait for all send requests related to this chunk
        // We remove the impact of -1 mark in wcnt_local[cid]
        MPI_Waitall(wcnt_all[cid] - wcnt_local[cid], rreqs + nrecv, rstats + nrecv);

        // Process data received
        //printf("nrecv = %d, wcnt_all = %d, wcnt_local = %d\n", nrecv, wcnt_all[cid], wcnt_local[cid]); fflush(stdout);
        for(j = nrecv; j < nrecv + wcnt_all[cid] - wcnt_local[cid]; j++){
            packoff = 0;
            //printf("rsize_2 = %d\n", rsizes[j]); fflush(stdout);
            while(packoff < rsizes[j]){
                MPI_Unpack(rbufs[j], rsizes[j], &packoff, tstart, varp->ndim, MPI_INT, nczipp->comm);
                MPI_Unpack(rbufs[j], rsizes[j], &packoff, tssize, varp->ndim, MPI_INT, nczipp->comm);

                for(k = 0; k < varp->ndim; k++){
                    tsize[k] = varp->chunkdim[k];
                }

                //printf("Rank: %d, MPI_Type_create_subarray([%d, %d], [%d, %d], [%d, %d]\n", nczipp->rank, tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
                MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                MPI_Type_commit(&ptype);

                //printf("tsize = %d, tssize = %d, tstart = %d, buf = %d\n", tsize[0], tssize[0], tstart[0], *((int*)(rbufs[j] + packoff))); fflush(stdout);
                MPI_Unpack(rbufs[j], rsizes[j], &packoff, varp->chunk_cache[cid], 1, ptype, nczipp->comm);
                //printf("cache[0] = %d, cache[1] = %d\n", ((int*)(varp->chunk_cache[cid]))[0], ((int*)(varp->chunk_cache[cid]))[1]); fflush(stdout);
                MPI_Type_free(&ptype);
            }
        }
        nrecv += wcnt_all[cid] - wcnt_local[cid];        

        //printbuf(nczipp->rank, varp->chunk_cache[cid], varp->chunksize);
    }

    // Free buffers
    NCI_Free(wcnt_local);
    NCI_Free(wcnt_all);

    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);
    NCI_Free(osize);
    NCI_Free(ostart);

    NCI_Free(cstart);
    NCI_Free(citr);
    NCI_Free(cend);

    NCI_Free(sreqs);
    NCI_Free(sstats);
    for(i = 0; i < nsend; i++){
        NCI_Free(sbufs[i]);
    }
    NCI_Free(sbufs);

    NCI_Free(rreqs);
    NCI_Free(rstats);
    for(i = 0; i < nrecv; i++){
        NCI_Free(rbufs[i]);
    }
    NCI_Free(rbufs);
    NCI_Free(rsizes);

    if (tbuf != NULL){
        NCI_Free(tbuf);
    }

    return NC_NOERR;
}


int
nczipioi_put_varn_new(NC_zip        *nczipp,
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
    nczipioi_put_varn_cb(nczipp, varp, nreq, starts, counts, NULL, bufs);

    // Write the compressed variable
    nczipioi_save_var(nczipp, varp);

    NCI_Free(bufs);

    return NC_NOERR;
}