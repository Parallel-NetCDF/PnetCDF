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
#include "../ncmpio/ncmpio_NC.h"

int nczipioi_save_var(NC_zip *nczipp, NC_zip_var *varp) {
    int i, j, k, l, err;
    int *zsizes, *zsizes_all;
    MPI_Datatype mtype, ftype;  // Memory and file datatype
    int wcnt;
    int *lens;
    MPI_Aint *disps;
    MPI_Status status;
    MPI_Offset *zoffs;
    MPI_Offset voff;
    void **zbufs;
    int zdimid, zvarid;
    int put_size;
    char name[128]; // Name of objects
    NC *ncp = (NC*)(nczipp->ncp);

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO)

    // Allocate buffer for compression
    zsizes = (int*)NCI_Malloc(sizeof(int) * varp->nchunk);
    zbufs = (void**)NCI_Malloc(sizeof(void*) * varp->nmychunk);
    zsizes_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunk);
    zoffs = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (varp->nchunk + 1));

    // Allocate buffer for I/O
    wcnt = 0;
    for(l = 0; l < varp->nmychunk; l++){
        k = varp->mychunks[l];
        if (varp->dirty[k]){
            wcnt++;
        }
    }
    if (nczipp->rank == varp->chunk_owner[0]){
        wcnt += 1;
    }
    lens = (int*)NCI_Malloc(sizeof(int) * wcnt);
    disps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);

    memset(zsizes, 0, sizeof(int) * varp->nchunk);

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_COM)

    // Compress each chunk we own
    if (varp->zip != NULL){
        varp->zip->init(MPI_INFO_NULL);
        for(l = 0; l < varp->nmychunk; l++){
            k = varp->mychunks[l];

            if (varp->dirty[k]){
                // Apply compression
                varp->zip->compress_alloc(varp->chunk_cache[k]->buf, varp->chunksize, zbufs + l, zsizes + k, varp->ndim, varp->chunkdim, varp->etype);
            }
        }
        varp->zip->finalize();
    }
    else{
        for(l = 0; l < varp->nmychunk; l++){
            k = varp->mychunks[l];
            if (varp->dirty[k]){
                zbufs[l] = varp->chunk_cache[k]->buf;
                zsizes[k] = varp->chunksize;
            }
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_COM)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_SYNC)

    // Sync compressed data size with other processes
    CHK_ERR_ALLREDUCE(zsizes, zsizes_all, varp->nchunk, MPI_INT, MPI_MAX, nczipp->comm);

    if (varp->metaoff < 0 || varp->expanded){ 
        zoffs[0] = varp->nchunkalloc * (sizeof(long long) + sizeof(int));
    }
    else{
        zoffs[0] = 0;
    }
    for(i = 0; i < varp->nchunk; i++){
        zoffs[i + 1] = zoffs[i] + zsizes_all[i];
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_SYNC)
    
    if (zoffs[varp->nchunk] > 0){   // No need to do I/O if no dirty chunk to write
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_INIT)

        /* Write comrpessed variable
        * We start by defining data variable and writing metadata
        * Then, we create buffer type and file type for data
        * Finally MPI collective I/O is used for writing data
        */

        // Enter redefine mode
        nczipp->driver->redef(nczipp->ncp);

        // Prepare data variable
        
        // Define dimension for data variable
        sprintf(name, "_datablock_dim_%d", nczipp->nwrite);
        err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[varp->nchunk], &zdimid);
        if (err != NC_NOERR) return err;

        // Define data variable
        sprintf(name, "_datablock_%d", nczipp->nwrite);
        err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(zvarid));
        if (err != NC_NOERR) return err;

        // Mark as data variable
        i = NC_ZIP_VAR_DATA;
        err = nczipp->driver->put_att(nczipp->ncp, zvarid, "_varkind", NC_INT, 1, &i, MPI_INT);
        if (err != NC_NOERR) return err;

        // Record serial
        nczipp->nwrite++;
        err = nczipp->driver->put_att(nczipp->ncp, NC_GLOBAL, "_nwrite", NC_INT, 1, &(nczipp->nwrite), MPI_INT);
        if (err != NC_NOERR) return err;

        // Metadata offset
        if (varp->metaoff < 0){
            err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_metaoffset", NC_INT64, 1, &(varp->metaoff), MPI_LONG_LONG);
            if (err != NC_NOERR) return err;
        }

        // Switch to data mode
        err = nczipp->driver->enddef(nczipp->ncp);
        if (err != NC_NOERR) return err;

        // Update metadata
        voff = ncp->vars.value[zvarid]->begin;
        for(i = 0; i < varp->nchunk; i++){
            if (zsizes_all[i] > 0){
                varp->chunk_index[i].len = zsizes_all[i];
                varp->chunk_index[i].off = zoffs[i] + voff - ncp->begin_var;
            }
        }

        if (varp->metaoff < 0 || varp->expanded){
            varp->metaoff = voff - ncp->begin_var;
            err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_metaoffset", NC_INT64, 1, &(varp->metaoff), MPI_LONG_LONG);
            if (err != NC_NOERR) return err;

            // unset expand flag
            varp->expanded = 0;
            
        }

        /* Carry out coll I/O
        * OpenMPI will fail when set view or do I/O on type created with MPI_Type_create_hindexed when count is 0
        * We use a dummy call inplace of type with 0 count
        */
        if (wcnt > 0){
            // Create file type
            l = 0;
            if (nczipp->rank == varp->chunk_owner[0]){  // First chunk owner writes metadata
                lens[l] = (varp->nchunk) * sizeof(long long);
                disps[l++] = (MPI_Aint)varp->metaoff + ncp->begin_var;

                lens[l] = (varp->nchunk) * sizeof(int);
                disps[l++] = (MPI_Aint)(varp->metaoff + ncp->begin_var + sizeof(long long) * varp->nchunkalloc);
            }
            for(i = 0; i < varp->nmychunk; i++){
                k = varp->mychunks[i];

                // Record compressed size
                if (varp->dirty[k]){
                    lens[l] = zsizes[k];
                    disps[l++] = (MPI_Aint)(varp->chunk_index[k].off) + ncp->begin_var;
                }
            }
            MPI_Type_create_hindexed(wcnt, lens, disps, MPI_BYTE, &ftype);
            CHK_ERR_TYPE_COMMIT(&ftype);

            // Create memory buffer type
            l = 0;
            if (nczipp->rank == varp->chunk_owner[0]){  // First chunk owner writes metadata
                lens[l] = (varp->nchunk) * sizeof(long long);
                disps[l++] = (MPI_Aint)varp->chunk_index;
            }
            for(i = 0; i < varp->nmychunk; i++){
                k = varp->mychunks[i];

                // Record compressed size
                if (varp->dirty[k]){
                    lens[l] = zsizes[k];
                    disps[l++] = (MPI_Aint)zbufs[i];
                }
            }
            err = MPI_Type_create_hindexed(wcnt, lens, disps, MPI_BYTE, &mtype);
            CHK_ERR_TYPE_COMMIT(&mtype);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_INIT)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_WR)

    #ifndef WORDS_BIGENDIAN // NetCDF data is big endian
            if (nczipp->rank == varp->chunk_owner[0]){
                //ncmpii_in_swapn(varp->chunk_index, varp->nchunk, sizeof(long long));
                //ncmpii_in_swapn(varp->data_lens, varp->nchunk, sizeof(int));
            }
    #endif

            // Perform MPI-IO
            // Set file view
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
            // Write data
            CHK_ERR_WRITE_AT_ALL(ncp->collective_fh, 0, MPI_BOTTOM, 1, mtype, &status);
            // Restore file view
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    #ifndef WORDS_BIGENDIAN // Switch back to little endian
            if (nczipp->rank == varp->chunk_owner[0]){
                //ncmpii_in_swapn(varp->chunk_index, varp->nchunk, sizeof(long long));
                //ncmpii_in_swapn(varp->data_lens, varp->nchunk, sizeof(int));
            }
    #endif

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_WR)

    #ifdef _USE_MPI_GET_COUNT
            MPI_Get_count(&status, MPI_BYTE, &put_size);
    #else
            MPI_Type_size(mtype, &put_size);
    #endif
            nczipp->putsize += put_size;

            // Free type
            MPI_Type_free(&ftype);
            MPI_Type_free(&mtype);
        }
        else{
            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_INIT)
            NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_WR)

            // Follow coll I/O with dummy call
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
            CHK_ERR_WRITE_AT_ALL(ncp->collective_fh, 0, MPI_BOTTOM, 0, MPI_BYTE, &status);
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_WR)
        }
    }

    // Free buffers
    NCI_Free(zsizes);
    NCI_Free(zsizes_all);
    NCI_Free(zoffs);
    for(l = 0; l < varp->nmychunk; l++){
        k = varp->mychunks[l];
        if (varp->dirty[k]){
            if (varp->zip != NULL){
                free(zbufs[l]);
            }
            varp->dirty[k] = 0;
        }
    }
    NCI_Free(zbufs);

    NCI_Free(lens);
    NCI_Free(disps);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO)

    return NC_NOERR;
}

int nczipioi_save_nvar(NC_zip *nczipp, int nvar, int *varids) {
    int i, j, k, l, err;
    int vid;    // Iterator for variable id
    int cid;    // Iterator for chunk id
    int max_nchunks = 0;
    int *zsizes, *zsizes_all;
    int *reqids;
    int nreq;
    MPI_Offset *zoffs, *zoffsp;
    MPI_Offset start, count, oldzoff, voff;
    MPI_Datatype mtype, ftype;  // Memory and file datatype
    int wcnt, ccnt, wcur, ccur;
    int *lens;
    MPI_Aint *mdisps, *fdisps;
    MPI_Status status;
    int put_size;
    void **zbufs;
    int *zdels;
    int zdimid, zvarid;
    char name[128]; // Name of objects
    NC_zip_var *varp;
    NC *ncp = (NC*)(nczipp->ncp);
    NC_var *ncvarp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_INIT)

    wcnt = 0;
    ccnt = 0;
    for(i = 0; i < nvar; i++){
        varp = nczipp->vars.data + varids[i];
        if (nczipp->rank == varp->chunk_owner[0]){
            wcnt += 1;
        }
        for(l = 0; l < varp->nmychunk; l++){
            k = varp->mychunks[l];
            if (varp->dirty[k]){
                ccnt++;
            }
        }
        total_nchunks += varp->nchunk + 1;
    }
    wcnt += ccnt;

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_INIT)

    // Allocate reqid for metadata
    reqids = (int*)NCI_Malloc(sizeof(int) * nvar * 2);

    // Allocate buffer for compression
    zsizes = (int*)NCI_Malloc(sizeof(int) * max_nchunks);
    zsizes_all = (int*)NCI_Malloc(sizeof(int) * max_nchunks);
    zbufs = (void**)NCI_Malloc(sizeof(void*) * ccnt);
    zdels = (int*)NCI_Malloc(sizeof(int) * ccnt);
    zoffs = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (total_nchunks + 1));

    // Allocate buffer file type
    mdisps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);
    lens = (int*)NCI_Malloc(sizeof(int) * wcnt);
    fdisps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);

    ccur = 0;
    zsizesp = zsizes + nvar;
    zsizes_allp = zsizes_all + nvar;
    for(vid = 0; vid < nvar; vid++){
        varp = nczipp->vars.data + varids[vid];

        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_COM)

        oldzoff = zoffs[varp->nchunk];

        memset(zsizes, 0, sizeof(int) * varp->nchunk);

        // Compress each chunk we own
        if (varp->zip != NULL){
            varp->zip->init(MPI_INFO_NULL);
            for(l = 0; l < varp->nmychunk; l++){
                cid = varp->mychunks[l];

                // Apply compression
                if (varp->dirty[cid]){
                    zdels[ccur] = 1;
                    varp->zip->compress_alloc(varp->chunk_cache[cid]->buf, varp->chunksize, zbufs + (ccur++), zsizesp + cid, varp->ndim, varp->chunkdim, varp->etype);
                }
            }
            varp->zip->finalize();
        }
        else{
            for(l = 0; l < varp->nmychunk; l++){
                cid = varp->mychunks[l];
                if (varp->dirty[cid]){
                    zsizesp[cid] = varp->chunksize;
                    zdels[ccur] = 0;
                    zbufs[ccur++] = varp->chunk_cache[cid]->buf;
                }
            }
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_COM)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_SYNC)

        // Sync compressed data size with other processes
        CHK_ERR_IALLREDUCE(zsizesp, zsizes_allp, varp->nchunk, MPI_INT, MPI_MAX, nczipp->comm, reqs + vid);

        if (varp->metaoff < 0 || varp->expanded){ 
            zsizes_all[vid] = varp->nchunkalloc * (sizeof(long long) + sizeof(int));
        }
        else{
            zsizes_all[vid] = 0;
        }
        
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_SYNC)

        zsizesp += varp->nchunk;
        zsizes_allp += varp->nchunk;
    }
    
    /* Write comrpessed variable
    * We start by defining data variable and writing metadata
    * Then, we create buffer type and file type for data
    * Finally MPI collective I/O is used for writing data
    */
   
    zsizes_allp = zsizes_all + nvar;
    for(vid = 0; vid < nvar; vid++){
        varp = nczipp->vars.data + varids[vid];

        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_SYNC)

        CHK_ERR_WAIT(reqs + vid, &status);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_SYNC)

        zsizes_allp += varp->nchunk;
    }

    zoffs[0] = 0;
    for(i = 0; i < total_nchunks; i++){
        zoffs[i + 1] = zoffs[i] + zsizes_all[i];
    }    

    if (zoffs[total_nchunks] > 0){   // No need to do I/O if no dirty chunk to write
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_INIT)

        // Prepare data variable

        // Enter redefine mode
        nczipp->driver->redef(nczipp->ncp);

        // Define dimension for data variable
        sprintf(name, "_datablock_dim_%d", nczipp->nwrite);
        err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[total_nchunks], &zdimid);
        if (err != NC_NOERR) return err;

        // Define data variable
        sprintf(name, "_datablock_%d", nczipp->nwrite);
        err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &zvarid);
        if (err != NC_NOERR) return err;

        // Mark as data variable
        i = NC_ZIP_VAR_DATA;
        err = nczipp->driver->put_att(nczipp->ncp, zvarid, "_varkind", NC_INT, 1, &i, MPI_INT);
        if (err != NC_NOERR) return err;

        // Record serial
        nczipp->nwrite++;
        err = nczipp->driver->put_att(nczipp->ncp, NC_GLOBAL, "_nwrite", NC_INT, 1, &(nczipp->nwrite), MPI_INT);
        if (err != NC_NOERR) return err;

        // Metadata offset
        for(vid = 0; vid < nvar; vid++){
            varp = nczipp->vars.data + varids[vid];
            if (varp->metaoff < 0){
                err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_metaoffset", NC_INT64, 1, &(varp->metaoff), MPI_LONG_LONG);
                if (err != NC_NOERR) return err;
            }
        }

        // Switch back to data mode
        err = nczipp->driver->enddef(nczipp->ncp);
        if (err != NC_NOERR) return err;

        voff = ncp->vars.value[zvarid]->begin;

        wcur = ccur = 0;
        for(vid = 0; vid < nvar; vid++){
            varp = nczipp->vars.data + varids[vid];

            if (varp->metaoff < 0 || varp->expanded){
                varp->metaoff = zoffs[vid] + voff - ncp->begin_var;
                err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_metaoffset", NC_INT64, 1, &(varp->metaoff), MPI_LONG_LONG);
                if (err != NC_NOERR) return err;

                // unset expand flag
                varp->expanded = 0;
            }

            if (nczipp->rank == varp->chunk_owner[0]){  // First chunk owner writes metadata
                lens[wcur] = varp->nchunk * sizeof(NC_zip_chunk_index_entry);
                fdisps[wcur] = (MPI_Aint)varp->metaoff + ncp->begin_var;
                mdisps[wcur++] = (MPI_Aint)(varp->chunk_index);
                
                //lens[wcur] = varp->nchunk * sizeof(int);
                //fdisps[wcur] = (MPI_Aint)(varp->metaoff + ncp->begin_var + sizeof(long long) * varp->nchunkalloc);
                //mdisps[wcur++] = (MPI_Aint)(varp->data_lens);
            }
        }

        nczipioi_sort_file_offset(wcur, fdisps, mdisps, lens);

        zsizes_allp = zsizes_all + nvar;
        zoffsp = zoffs + nvar;
        for(vid = 0; vid < nvar; vid++){
            varp = nczipp->vars.data + varids[vid];

            for(cid = 0; cid < varp->nchunk; cid++){
                if (zsizes_allp[cid] > 0){
                    varp->chunk_index[cid].len = zsizes_allp[cid];
                    varp->chunk_index[cid].off = zoffsp[cid] + voff - ncp->begin_var;
                }
            }

            /* Paramemter for file and memory type 
            * We do not know variable file offset until the end of define mode
            * We will add the displacement later
            */
            for(i = 0; i < varp->nmychunk; i++){
                cid = varp->mychunks[i];

                // Record parameter
                if (varp->dirty[cid]){
                    lens[wcur] = varp->chunk_index[cid].len;
                    fdisps[wcur] = (MPI_Aint)(varp->chunk_index[cid].off) + ncp->begin_var;
                    mdisps[wcur++] = (MPI_Aint)zbufs[ccur++];
                }
            }

            // Clear dirty flag
            memset(varp->dirty, 0, varp->nchunk * sizeof(int));

            zsizes_allp += varp->nchunk;
            zoffsp += varp->nchunk;
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_INIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT_IO_WR)

        /* Carry our coll I/O
        * OpenMPI will fail when set view or do I/O on type created with MPI_Type_create_hindexed when count is 0
        * We use a dummy call inplace of type with 0 count
        */
        if (wcnt > 0){
            // Create file type
            MPI_Type_create_hindexed(wcnt, lens, fdisps, MPI_BYTE, &ftype);
            CHK_ERR_TYPE_COMMIT(&ftype);

            // Create memmory type
            MPI_Type_create_hindexed(wcnt, lens, mdisps, MPI_BYTE, &mtype);
            CHK_ERR_TYPE_COMMIT(&mtype);

    #ifndef WORDS_BIGENDIAN // NetCDF data is big endian
            for(vid = 0; vid < nvar; vid++){
                varp = nczipp->vars.data +  varids[vid];
                if (nczipp->rank == varp->chunk_owner[0]){
                    //ncmpii_in_swapn(varp->chunk_index, varp->nchunk + 1, sizeof(long long));
                    //ncmpii_in_swapn(varp->data_lens, varp->nchunk + 1, sizeof(int));
                }
            }
    #endif     

            // Perform MPI-IO
            // Set file view
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
            // Write data
            CHK_ERR_WRITE_AT_ALL(ncp->collective_fh, 0, MPI_BOTTOM, 1, mtype, &status);
            // Restore file view
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    #ifndef WORDS_BIGENDIAN // Switch back to little endian
            for(vid = 0; vid < nvar; vid++){
                varp = nczipp->vars.data +  varids[vid];
                if (nczipp->rank == varp->chunk_owner[0]){
                    //ncmpii_in_swapn(varp->chunk_index, varp->nchunk + 1, sizeof(long long));
                    //ncmpii_in_swapn(varp->data_lens, varp->nchunk + 1, sizeof(int));
                }
            }
    #endif    

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_WR)

    #ifdef _USE_MPI_GET_COUNT
            MPI_Get_count(&status, MPI_BYTE, &put_size);
    #else
            MPI_Type_size(mtype, &put_size);
    #endif
            nczipp->putsize += put_size;

            // Free type
            MPI_Type_free(&ftype);
            MPI_Type_free(&mtype);
        }
        else{
            // Follow coll I/O with dummy call
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
            CHK_ERR_WRITE_AT_ALL(ncp->collective_fh, 0, MPI_BOTTOM, 0, MPI_BYTE, &status);
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

            NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO_WR)
        }
    }

    // Free buffers
    NCI_Free(zsizes);
    if (varp->zip != NULL){
        for(l = 0; l < ccnt; l++){
            free(zbufs[l]);
        }
    }
    NCI_Free(zbufs);

    NCI_Free(lens);
    NCI_Free(fdisps);
    NCI_Free(mdisps);

    NCI_Free(reqids);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT_IO)

    return NC_NOERR;
}
