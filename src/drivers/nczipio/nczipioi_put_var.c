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
nczipioi_put_var_cb_chunk(NC_zip          *nczipp,
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
    int *tsize, *tssize, *tstart, *tsizep, *tstartp;   // Size for sub-array type
    MPI_Offset *citr;

    int *wcnt_local, *wcnt_all;   // Number of processes that writes to each chunk

    int overlapsize;    // Size of overlaping region of request and chunk
    int max_tbuf = 0;   // Size of intermediate buffer
    char *tbuf = NULL;     // Intermediate buffer
    
    int packoff; // Pack offset
    MPI_Datatype ptype; // Pack datatype

    int nsend, nrecv;   // Number of send and receive
    MPI_Request *sreqs, *rreqs;    // Send and recv req
    MPI_Status *sstats, *rstats;    // Send and recv status
    char **sbufs, **rbufs;   // Send and recv buffer
    int *rsizes;    // recv size of each message
    MPI_Message rmsg;   // Receive message

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_INIT)

    // Allocate buffering for write count
    wcnt_local = (int*)NCI_Malloc(sizeof(int) * varp->nchunk * 2);
    wcnt_all = wcnt_local + varp->nchunk;

    // Allocate buffering for overlaping index
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim * 3);
    tsize = tstart + varp->ndim;
    tssize = tsize + varp->ndim;
    ostart = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim * 3);
    osize = ostart + varp->ndim;

    // Chunk iterator
    citr = osize + varp->ndim;
    
    // We need to calculate the size of message of each chunk
    // This is just for allocating send buffer
    // We do so by iterating through all request and all chunks they cover
    // If we are not the owner of a chunk, we need to send message
    memset(wcnt_local, 0, sizeof(int) * varp->nchunk);
    nsend = 0;

    // Iterate through chunks
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid); // Initialize chunk iterator
    do{
        if (varp->chunk_owner[cid] != nczipp->rank){
            // Count number of mnessage we need to send
            nsend++;
            wcnt_local[cid] = 1;
        }
        else{
            // We mark covered chunk of our own to prevent unnecessary calculation of overlap
            // -1 is purely a mark, we need to add 1 back to global message count
            wcnt_local[cid] = -1;
            max_tbuf = varp->chunksize;
        }
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    // Allocate buffer for sending
    sbufs = (char**)NCI_Malloc(sizeof(char*) * nsend);
    sreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nsend);
    sstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nsend);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_INIT)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_SYNC)

    // Sync number of messages of each chunk
    CHK_ERR_ALLREDUCE(wcnt_local, wcnt_all, varp->nchunk, MPI_INT, MPI_SUM, nczipp->comm); 

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_SYNC)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_PACK_REQ)

    // Calculate number of recv request
    // This is for all the chunks
    nrecv = 0;
    for(i = 0; i < varp->nmychunk; i++){
        cid = varp->mychunks[i];
        // We don't need message for our own data
        nrecv += wcnt_all[cid] - wcnt_local[cid];
    }
    rreqs = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nrecv);
    rstats = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nrecv);
    rbufs = (char**)NCI_Malloc(sizeof(char*) * nrecv);
    rsizes = (int*)NCI_Malloc(sizeof(int) * nrecv);

    // Post send
    nsend = 0;
    // Iterate through chunks
    nczipioi_chunk_itr_init_ex(varp, start, count, citr, &cid, ostart, osize); // Initialize chunk iterator
    do{
        // We got something to send if we are not owner
        if (varp->chunk_owner[cid] != nczipp->rank){
            // Calculate chunk overlap
            overlapsize = varp->esize;
            for(j = 0; j < varp->ndim; j++){
                overlapsize *= osize[j];                     
            }

            // Allocate buffer
            sbufs[nsend] = (char*)NCI_Malloc(overlapsize + sizeof(int) * varp->ndim * 2);
            
            // Metadata
            packoff = 0;
            tstartp = (int*)sbufs[nsend]; packoff += varp->ndim * sizeof(int);
            tsizep = (int*)(sbufs[nsend] + packoff); packoff += varp->ndim * sizeof(int);
            for(j = 0; j < varp->ndim; j++){
                tstartp[j] = (int)(ostart[j] - citr[j]);
                tsizep[j] = (int)osize[j];
            }
 
            // Pack type
            for(j = 0; j < varp->ndim; j++){
                tstart[j] = (int)(ostart[j] - start[j]);
                tsize[j] = (int)count[j];
                tssize[j] = (int)osize[j];
            }
            err = MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            if (err != MPI_SUCCESS){
                err = ncmpii_error_mpi2nc(err, "MPI_Type_create_subarray");
                DEBUG_RETURN_ERROR(err)
            }
            err = MPI_Type_commit(&ptype);
            if (err != MPI_SUCCESS){
                err = ncmpii_error_mpi2nc(err, "MPI_Type_commit");
                DEBUG_RETURN_ERROR(err)
            }

            // Pack data
            CHK_ERR_PACK(buf, 1, ptype, sbufs[nsend], packoff + overlapsize, &packoff, nczipp->comm);

            MPI_Type_free(&ptype);

            NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_SEND_REQ)
            
            // Send the request
            CHK_ERR_ISEND(sbufs[nsend], packoff, MPI_BYTE, varp->chunk_owner[cid], cid, nczipp->comm, sreqs + nsend);

            NC_ZIP_TIMER_STOPEX(NC_ZIP_TIMER_PUT_CB_SEND_REQ, NC_ZIP_TIMER_PUT_CB_PACK_REQ)
            nsend++;
        }
    } while (nczipioi_chunk_itr_next_ex(varp, start, count, citr, &cid, ostart, osize));

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_PACK_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_RECV_REQ)

    // Post recv
    nrecv = 0;
    for(j = 0; j < varp->nmychunk; j++){
        cid = varp->mychunks[j];
        // We are the owner of the chunk
        // Receive data from other process
        for(i = 0; i < wcnt_all[cid] - wcnt_local[cid]; i++){
            // Get message size, including metadata
            MPI_Mprobe(MPI_ANY_SOURCE, cid, nczipp->comm, &rmsg, rstats);
            MPI_Get_count(rstats, MPI_BYTE, rsizes + nrecv);

            // Allocate buffer
            rbufs[nrecv] = (char*)NCI_Malloc(rsizes[nrecv]);

            // Post irecv
            MPI_Imrecv(rbufs[nrecv], rsizes[nrecv], MPI_BYTE, &rmsg, rreqs + nrecv);
            nrecv++;
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_RECV_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_SEND_REQ)

    // Wait for all send
    MPI_Waitall(nsend, sreqs, sstats);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_SEND_REQ)

    // Allocate intermediate buffer
    if (max_tbuf > 0){
        tbuf = (char*)NCI_Malloc(max_tbuf);
    }

    // For each chunk we own, we need to receive incoming data
    nrecv = 0;
    for(i = 0; i < varp->nmychunk; i++){
        cid = varp->mychunks[i];

        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_SELF)   

        // Handle our own data first if we have any
        if (wcnt_local[cid] < 0){
            // Convert chunk id to iterator
            get_chunk_itr(varp, cid, citr);

            // Calculate overlapping region
            overlapsize = get_chunk_overlap(varp, citr, start, count, ostart, osize);

            if (overlapsize > 0){
                // Pack type from user buffer to (contiguous) intermediate buffer
                for(j = 0; j < varp->ndim; j++){
                    tstart[j] = (int)(ostart[j] - start[j]);
                    tsize[j] = (int)count[j];
                    tssize[j] = (int)osize[j];
                }
                CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                CHK_ERR_TYPE_COMMIT(&ptype);

                // Pack data into intermediate buffer
                packoff = 0;
                MPI_Pack(buf, 1, ptype, tbuf, varp->chunksize, &packoff, nczipp->comm);
                overlapsize = packoff;

                MPI_Type_free(&ptype);

                // Pack type from (contiguous) intermediate buffer to chunk buffer
                for(j = 0; j < varp->ndim; j++){
                    tstart[j] = (int)(ostart[j] - citr[j]);
                    tsize[j] = varp->chunkdim[j];
                }
                MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
                MPI_Type_commit(&ptype);
                
                // Unpack data into chunk buffer
                packoff = 0;
                MPI_Unpack(tbuf, overlapsize, &packoff, varp->chunk_cache[cid], 1, ptype, nczipp->comm);

                MPI_Type_free(&ptype); 

                varp->dirty[cid] = 1;   
            }
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_SELF)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_RECV_REQ)

        // Now, it is time to process data from other processes
        
        // Wait for all send requests related to this chunk
        // We remove the impact of -1 mark in wcnt_local[cid]
        MPI_Waitall(wcnt_all[cid] - wcnt_local[cid], rreqs + nrecv, rstats + nrecv);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_RECV_REQ)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_UNPACK_REQ)

        // Process data received
        for(j = nrecv; j < nrecv + wcnt_all[cid] - wcnt_local[cid]; j++){
            // Metadata
            packoff = 0;
            tstartp = (int*)rbufs[j]; packoff += varp->ndim * sizeof(int);
            tsizep = (int*)(rbufs[j] + packoff); packoff += varp->ndim * sizeof(int);
            
            // Pack type
            CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, varp->chunkdim, tsizep, tstartp, MPI_ORDER_C, varp->etype, &ptype);
            CHK_ERR_TYPE_COMMIT(&ptype);

            // Data
            CHK_ERR_UNPACK(rbufs[j], rsizes[j], &packoff, varp->chunk_cache[cid], 1, ptype, nczipp->comm);
            MPI_Type_free(&ptype);

            varp->dirty[cid] = 1;
        }
        nrecv += wcnt_all[cid] - wcnt_local[cid]; 

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_UNPACK_REQ)
    }

    // Free buffers
    NCI_Free(wcnt_local);

    NCI_Free(tstart);
    
    NCI_Free(ostart);

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

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB)

    return NC_NOERR;
}

int
nczipioi_put_var_cb_proc(   NC_zip          *nczipp,
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
    int *tsize, *tssize, *tstart, *tssizep, *tstartp;   // Size for sub-array type
    MPI_Offset *citr; // Bounding box for chunks overlapping my own write region
    
    int *wcnt_local, *wcnt_all;   // Number of processes that writes to each chunk

    int overlapsize;    // Size of overlaping region of request and chunk
    int max_tbuf = 0;   // Size of intermediate buffer
    char *tbuf = NULL;     // Intermediate buffer
    
    int packoff; // Pack offset
    MPI_Datatype ptype; // Pack datatype

    int nsend, nrecv;   // Number of send and receive
    MPI_Request *sreq, *rreq;    // Send and recv req
    MPI_Status *sstat, rstat;    // Send and recv status
    char **sbuf, **sbufp, **rbuf, **rbufp;   // Send and recv buffer
    int *rsize, *ssize;    // recv size of each message
    int *sdst;    // recv size of each message
    int *smap;
    MPI_Message rmsg;   // Receive message

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_INIT)

    // Allocate buffering for write count
    wcnt_local = (int*)NCI_Malloc(sizeof(int) * nczipp->np * 3);
    wcnt_all = wcnt_local + nczipp->np;
    smap = wcnt_all + nczipp->np;

    // Allocate buffering for overlaping index
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim * 3);
    tssize = tstart + varp->ndim;
    tsize = tssize + varp->ndim;
    ostart = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim * 3);
    osize = ostart + varp->ndim;

    // Chunk iterator
    citr = osize + varp->ndim;

    // We need to calculate the size of message of each chunk
    // This is just for allocating send buffer
    // We do so by iterating through all request and all chunks they cover
    // If we are not the owner of a chunk, we need to send message
    memset(wcnt_local, 0, sizeof(int) * nczipp->np);
    nsend = 0;

    // Count total number of messages and build a map of accessed chunk to list of comm datastructure
    nczipioi_chunk_itr_init(varp, start, count, citr, &cid); // Initialize chunk iterator
    do{
        // Chunk owner
        cown = varp->chunk_owner[cid];

        // Mapping to skip list of send requests 
        if (wcnt_local[cown] == 0 && cown != nczipp->rank){
            smap[cown] = nsend++;
        }
        wcnt_local[cown] = 1;   // Need to send message if not owner       
    } while (nczipioi_chunk_itr_next(varp, start, count, citr, &cid));

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_INIT)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_SYNC)

    // Sync number of messages of each chunk
    MPI_Allreduce(wcnt_local, wcnt_all, nczipp->np, MPI_INT, MPI_SUM, nczipp->comm);
    nrecv = wcnt_all[nczipp->rank] - wcnt_local[nczipp->rank];  // We don't need to receive request form self

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_SYNC)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_PACK_REQ)

    // Allocate data structure for messaging
    sbuf = (char**)NCI_Malloc(sizeof(char*) * nsend * 2);
    sbufp = sbuf + nsend;
    ssize = (int*)NCI_Malloc(sizeof(int) * nsend * 2);
    sdst = ssize + nsend;
    sreq = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nsend);
    sstat = (MPI_Status*)NCI_Malloc(sizeof(MPI_Status) * nsend);

    rbuf = (char**)NCI_Malloc(sizeof(char*) * nrecv * 2);
    rbufp = rbuf + nrecv;
    rsize = (int*)NCI_Malloc(sizeof(int) * nrecv);
    rreq = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nrecv);

    // Count size of each request
    memset(ssize, 0, sizeof(int) * nsend);
    nczipioi_chunk_itr_init_ex(varp, start, count, citr, &cid, ostart, osize); // Initialize chunk iterator
    do{
        // Chunk owner
        cown = varp->chunk_owner[cid];
        if (cown != nczipp->rank){
            j = smap[cown];
            sdst[j] = cown; // Record a reverse map by the way

            // Count overlap
            overlapsize = varp->esize;
            for(i = 0; i < varp->ndim; i++){
                overlapsize *= osize[i];                     
            }
            ssize[j] += overlapsize + sizeof(int) * (varp->ndim * 2 + 1);
        }
    } while (nczipioi_chunk_itr_next_ex(varp, start, count, citr, &cid, ostart, osize));

    // Allocate buffer for send
    for(i = 0; i < nsend; i++){
        sbuf[i] = sbufp[i] = (char*)NCI_Malloc(ssize[i]);
    }

    // Pack requests
    nczipioi_chunk_itr_init_ex(varp, start, count, citr, &cid, ostart, osize); // Initialize chunk iterator
    do{
        // Chunk owner
        cown = varp->chunk_owner[cid];
        if (cown != nczipp->rank){
            j = smap[cown];

            // Metadata
            *((int*)sbufp[j]) = cid; sbufp[j] += sizeof(int);
            tstartp = (int*)sbufp[j];  sbufp[j] += varp->ndim * sizeof(int);
            tssizep = (int*)sbufp[j];  sbufp[j] += varp->ndim * sizeof(int);
            for(i = 0; i < varp->ndim; i++){
                tstartp[i] = (int)(ostart[i] - citr[i]);
                tssizep[i] = (int)osize[i];
            }

            // Pack type 
            for(i = 0; i < varp->ndim; i++){
                tstart[i] = (int)(ostart[i] - start[i]);
                tsize[i] = (int)count[i];
            }
            CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssizep, tstart, MPI_ORDER_C, varp->etype, &ptype);
            CHK_ERR_TYPE_COMMIT(&ptype);

            // Data
            packoff = 0;
            CHK_ERR_PACK(buf, 1, ptype, sbufp[j], ssize[j], &packoff, nczipp->comm);    sbufp[j] += packoff;
            MPI_Type_free(&ptype);
        }
    } while (nczipioi_chunk_itr_next_ex(varp, start, count, citr, &cid, ostart, osize));

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_PACK_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_SEND_REQ)

    // Post send
    for(i = 0; i < nsend; i++){
        CHK_ERR_ISEND(sbuf[i], ssize[i], MPI_BYTE, sdst[i], 0, nczipp->comm, sreq + i);
    }    

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_SEND_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_RECV_REQ)

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

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_RECV_REQ)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_SELF)

    tbuf = (char*)NCI_Malloc(varp->chunksize);

    // Handle our own data
    nczipioi_chunk_itr_init_ex(varp, start, count, citr, &cid, ostart, osize); // Initialize chunk iterator
    do{
        if (varp->chunk_owner[cid] == nczipp->rank){
            // Pack type from user buffer to (contiguous) intermediate buffer
            for(j = 0; j < varp->ndim; j++){
                tstart[j] = (int)(ostart[j] - start[j]);
                tsize[j] = (int)count[j];
                tssize[j] = (int)osize[j];
            }
            CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            CHK_ERR_TYPE_COMMIT(&ptype);

            // Pack data into intermediate buffer
            packoff = 0;
            CHK_ERR_PACK(buf, 1, ptype, tbuf, varp->chunksize, &packoff, nczipp->comm);
            MPI_Type_free(&ptype);
            overlapsize = packoff;

            // Pack type from (contiguous) intermediate buffer to chunk buffer
            for(j = 0; j < varp->ndim; j++){
                tstart[j] = (int)(ostart[j] - citr[j]);
                tsize[j] = varp->chunkdim[j];
            }            
            CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, varp->etype, &ptype);
            CHK_ERR_TYPE_COMMIT(&ptype);
            
            // Unpack data into chunk buffer
            packoff = 0;
            CHK_ERR_UNPACK(tbuf, overlapsize, &packoff, varp->chunk_cache[cid], 1, ptype, nczipp->comm);
            MPI_Type_free(&ptype);    

            varp->dirty[cid] = 1;
        }
    } while (nczipioi_chunk_itr_next_ex(varp, start, count, citr, &cid, ostart, osize));

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_SELF)

    //Handle incoming requests
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = varp->chunkdim[i];
    }
    for(i = 0; i < nrecv; i++){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_RECV_REQ)

        // Will wait any provide any benefit?
        MPI_Waitany(nrecv, rreq, &j, &rstat);
        
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_RECV_REQ)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_UNPACK_REQ)
                
        while(rbufp[j] < rbuf[j] + rsize[j]){
            // Metadata
            cid = *(int*)(rbufp[j]); rbufp[j] += sizeof(int);
            tstartp = (int*)rbufp[j];  rbufp[j] += varp->ndim * sizeof(int);
            tssizep = (int*)rbufp[j];  rbufp[j] += varp->ndim * sizeof(int);
           
            // Pack type
            CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssizep, tstartp, MPI_ORDER_C, varp->etype, &ptype);
            CHK_ERR_TYPE_COMMIT(&ptype);

            // Data
            packoff = 0;
            CHK_ERR_UNPACK(rbufp[j], rsize[j], &packoff, varp->chunk_cache[cid], 1, ptype, nczipp->comm);   rbufp[j] += packoff;
            MPI_Type_free(&ptype);

            varp->dirty[cid] = 1;
        }
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_UNPACK_REQ)
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_CB_SEND_REQ)

    MPI_Waitall(nsend, sreq, sstat);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB_SEND_REQ)

    // Free buffers
    NCI_Free(wcnt_local);

    NCI_Free(tstart);

    NCI_Free(ostart);

    NCI_Free(sreq);
    NCI_Free(sstat);
    NCI_Free(ssize);

    for(i = 0; i < nsend; i++){
        NCI_Free(sbuf[i]);
    }
    NCI_Free(sbuf);

    NCI_Free(rreq);
    for(i = 0; i < nrecv; i++){
        NCI_Free(rbuf[i]);
    }
    NCI_Free(rbuf);
    NCI_Free(rsize);

    if (tbuf != NULL){
        NCI_Free(tbuf);
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_CB)

    return NC_NOERR;
}

int
nczipioi_put_var(NC_zip        *nczipp,
              NC_zip_var       *varp,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              void       *buf)
{
    if (varp->isrec){
        if (nczipp->recsize < start[0] + count[0]){
            nczipp->recsize = start[0] + count[0];
        }
        MPI_Allreduce(MPI_IN_PLACE, &(nczipp->recsize), 1, MPI_LONG_LONG, MPI_MAX, nczipp->comm);   // Sync number of recs
        if (varp->dimsize[0] < nczipp->recsize){
            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_RESIZE)

            nczipioi_var_resize(nczipp, varp);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_RESIZE)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT)
        }
    }

    // Collective buffer
    switch (nczipp->comm_unit){
        case NC_ZIP_COMM_CHUNK:
            nczipioi_put_var_cb_chunk(nczipp, varp, start, count, stride, buf);
            break;
        case NC_ZIP_COMM_PROC:
            nczipioi_put_var_cb_proc(nczipp, varp, start, count, stride, buf);
            break;
    }

    // Write the compressed variable
    nczipioi_save_var(nczipp, varp);

    return NC_NOERR;
}