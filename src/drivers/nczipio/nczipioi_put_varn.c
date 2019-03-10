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

int
nczipioi_handle_nreq(NC_zip        *nczipp,
              NC_zip_var       *varp,
              int              nreq,
              MPI_Offset* const *starts,
              MPI_Offset* const *counts,
              MPI_Offset* const *strides,
              void              **bufs)
{
    int i, j, k, l, err;
    MPI_Datatype etype; // Variable element type in MPI
    int esize;  // Variable element size
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int *cstart, *cend, *citr; // Bounding box for chunks overlapping my own write region
    int *wcnt_local, *wcnt_all;   // Number of processes that writes to each chunk
    char *sbuf_base, *sbuf_cur; // Send buffer, exactly the same size as buf
    int put_size_total, put_size;  // Total size of buf and size of data of a single req
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
    zbufs = (char**)NCI_Malloc(sizeof(char*) * varp->nmychunks);

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

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
    for(k = 0; k < nreq; k++){
        // Initialize chunk iterator
        nczipioi_chunk_itr_init(varp, starts[k], counts[k], NULL, cstart, cend, citr);

        // Iterate through chunks
        do{
            // Chunk index
            i = get_chunk_idx(varp, citr);

            if (varp->chunk_owner[i] != nczipp->rank){
                get_chunk_overlap(varp, citr, starts[i], counts[i], tstart, tssize);
                //printf("cord = %d, start = %lld, count = %lld, tstart = %d, tssize = %d, esize = %d, ndim = %d\n", citr[0], starts[i][0], counts[i][0], tstart[0], tssize[0], esize, varp->ndim); fflush(stdout);
                
                overlapsize = varp->esize;
                for(j = 0; j < varp->ndim; j++){
                    overlapsize *= tssize[j];                     
                }

                // Count number of mnessage we need to send
                if (wcnt_local[i] == 0 && varp->chunk_owner[i] != nczipp->rank){
                    nsend++;
                }

                wcnt_local[i] += overlapsize + sizeof(int) * 2 * varp->ndim;
            }
            else{
                // We mark covered chunk of our own to prevent unnecessary calculation of overlap
                // -1 is purely a mark, we need to add 1 back to global message count
                wcnt_local[i] = -1;
            }

        } while (nczipioi_chunk_itr_next(varp, start, count, stride, c start, cend, citr));
    }

    // Allocate buffer for sending
    sbufs = (char**)NCI_Malloc(sizeof(char*) * nsend);
    sreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nsend);
    sstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nsend);
    j = 0;
    // Allocate buffer for data
    for(k = 0; k < varp->nchunks; k++){
        // Count number of mnessage we need to send
        if (wcnt_local[k] > 0){
            // Add space for number of reqs
            sbufs[j++] = (char*)NCI_Malloc(wcnt_local[k]);
            // We don't need message size anymore, wcnt_local is used to track number of message from now on 
            wcnt_local[k] = 1;
        }
    }

    // Sync number of messages of each chunk
    MPI_Allreduce(wcnt_local, wcnt_all, nczipp->np, MPI_INT, MPI_SUM, nczipp->comm);

    // Calculate number of recv request
    // This is for all the chunks
    nrecv = 0;
    for(l = 0; l < varp->nmychunkss; l++){
        i = varp->mychunks[l];
        // We add 1 back form the -1 we use as mark
        nrecv += wcnt_all[i] + 1;
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
            // We are the owner of the chunk
            // Receive data from other process
            for(i = 0; i <= wcnt_all[k]; i++){
                // Get message size, including metadata
                MPI_Mprobe(MPI_ANY_SOURCE, k, nczipp->comm, &rmsg, rstats);
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
            if (wcnt_local[k] < 0){
                packoff = 0;
                // Get chunk iterator
                get_chunk_cord(varp, k, citr);  
                for(i = 0; i < nreq; i++){
                    // Calculate chunk overlap
                    get_chunk_overlap(varp, citr, starts[i], counts[i], tstart, tssize);
                    printf("cord = %d, start = %lld, count = %lld, tstart = %d, tssize = %d, esize = %d, ndim = %d\n", citr[0], starts[i][0], counts[i][0], tstart[0], tssize[0], esize, varp->ndim); fflush(stdout);
                    overlapsize = esize;
                    for(j = 0; j < varp->ndim; j++){
                        overlapsize *= tssize[j];                     
                    }
                    printf("overlapsize = %d\n", overlapsize); fflush(stdout);

                    // If current request have any overlap with the chunk, we pack the data and metadata
                    if (overlapsize > 0){
                        // Pack type
                        for(j = 0; j < varp->ndim; j++){
                            tstart[j] -= starts[i][j];
                            tsize[j] = (int)counts[i][j];
                        }
                        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, etype, &ptype);
                        MPI_Type_commit(&ptype);
                        
                        // Pack data
                        MPI_Pack(starts[i], varp->ndim, MPI_INT, sbufs[nsend], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        MPI_Pack(counts[i], varp->ndim, MPI_INT, sbufs[nsend], packoff + sizeof(int) * varp->ndim, &packoff, nczipp->comm);
                        MPI_Pack(bufs[i], 1, ptype, sbufs[nsend], packoff + overlapsize, &packoff, nczipp->comm);

                        MPI_Type_free(&ptype);
                    }
                }

                // Send the request
                printf("packoff = %d\n", packoff); fflush(stdout);
                MPI_Isend(sbufs[nsend], packoff, MPI_BYTE, varp->chunk_owner[k], k, nczipp->comm, sreqs + nsend);
                nsend++;
            }
        }
    }

    // For each chunk we own, we need to receive incoming data
    nrecv = 0;
    memset(zsizes, 0, sizeof(int) * varp->nchunks);
    for(l = 0; l < varp->nmychunks; l++){
        k = varp->mychunks[l];
        
        if (varp->chunk_owner[k] == nczipp->rank){

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

    // Wait for all send
    MPI_Waitall(nsend, sreqs, sstats);

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
    for(l = 0; l < varp->nmychunks; l++){
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
    for(l = 0; l < varp->nmychunks; l++){
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
nczipioi_put_varn(NC_zip        *nczipp,
              NC_zip_var       *varp,
              int              nreq,
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
    zbufs = (char**)NCI_Malloc(sizeof(char*) * varp->nmychunks);

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
    for(l = 0; l < varp->nmychunks; l++){
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
    for(l = 0; l < varp->nmychunks; l++){
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
    for(l = 0; l < varp->nmychunks; l++){
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
