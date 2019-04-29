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

int nczipioi_var_init(NC_zip *nczipp, NC_zip_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
    int i, j, err;
    int valid;
    MPI_Offset len;
    NC_zip_var *var;

    if (varp->varkind == NC_ZIP_VAR_COMPRESSED){
        if (varp->chunkdim == NULL){    // This is a new uninitialized variable 
            // Update dimsize on rec dim
            if (nczipp->recdim >= 0){
                if (varp->dimsize[0] < nczipp->recsize){
                    varp->dimsize[0] = nczipp->recsize;
                }
            }

            // Determine its block size
            varp->chunkdim = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
            varp->nchunks = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

            // First check attribute
            valid = 1;
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_chunkdim", NULL, &len);
            if (err == NC_NOERR && len == varp->ndim){
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_chunkdim", varp->chunkdim, MPI_INT);
                if (err != NC_NOERR){
                    valid = 0;
                }
                //chunkdim must be at leasst 1
                for(j = 0; j < varp->ndim; j++){ 
                    if (varp->chunkdim[j] <= 0){
                        valid = 0;
                        printf("Warning: block size invalid, use default");
                        break;
                    }
                }
            }
            else{
                valid = 0;
            }
            
            // Default block size is same as dim size, only 1 blocks
            if (!valid){
                err = nczipioi_calc_chunk_size(nczipp, varp, nreq, starts, counts);
                if (err != NC_NOERR){
                    return err;
                }
                if (!nczipp->delay_init){
                    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkdim", NC_INT, varp->ndim, varp->chunkdim, MPI_INT);
                    if (err != NC_NOERR){
                        return err;
                    }
                }
            }

            // Calculate total # chunks, # chunks along each dim, chunksize
            varp->nchunk = 1;
            varp->chunksize = NC_Type_size(varp->xtype);
            for(i = 0; i < varp->ndim; i++){ //chunkdim must be at leasst 1
                if (varp->dimsize[i] % varp->chunkdim[i] == 0){
                    varp->nchunks[i] = (int)varp->dimsize[i] / varp->chunkdim[i];
                }
                else{
                    varp->nchunks[i] = (int)varp->dimsize[i] / varp->chunkdim[i] + 1;
                }
                varp->nchunk *= varp->nchunks[i];
                varp->chunksize *= varp->chunkdim[i];
            }

            // Calculate number of chunks below each dimension
            varp->cidsteps = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
            varp->cidsteps[varp->ndim - 1] = 1;
            for(i = varp->ndim - 2; i >= 0; i--){
                varp->cidsteps[i] = varp->cidsteps[i + 1] * varp->nchunks[i + 1];
            }

            // Determine block ownership
            varp->chunk_owner = (int*)NCI_Malloc(sizeof(int) * varp->nchunk);
            varp->chunk_cache = (char**)NCI_Malloc(sizeof(char*) * varp->nchunk);
            memset(varp->chunk_cache, 0, sizeof(char*) * varp->nchunk);
            // We infer owners by reqs
            err = nczipioi_calc_chunk_owner(nczipp, varp, nreq, starts, counts);
            if (err != NC_NOERR){
                return err;
            }

            // Build skip list of my own chunks
            varp->nmychunks = 0;
            for(j = 0; j < varp->nchunk; j++){ 
                if (varp->chunk_owner[j] == nczipp->rank){
                    varp->nmychunks++;
                }
            }
            varp->mychunks = (int*)NCI_Malloc(sizeof(int) * varp->nmychunks);
            varp->nmychunks = 0;
            for(j = 0; j < varp->nchunk; j++){ 
                if (varp->chunk_owner[j] == nczipp->rank){
                    varp->mychunks[varp->nmychunks++] = j;
                    if (varp->isnew){   // Only apply to new var, old var will be read when it is needed
                        varp->chunk_cache[j] = (void*)NCI_Malloc(varp->chunksize);  // Allocate buffer for blocks we own
                    }
                }
            }
            
            // Update global chunk count
            nczipp->nmychunks += varp->nmychunks;

            /*
            if (nczipp->rank == 0){
                printf("Var %d, cown = [", varp->varid);
                for(i = 0; i < varp->nchunk; i++)
                    printf("%d, ", varp->chunk_owner[i]);
                printf("]\n");
            }
            */

            // Determine block offset
            varp->data_offs = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (varp->nchunk + 1));
            varp->data_lens = (int*)NCI_Malloc(sizeof(int) * varp->nchunk);
            // Try if there are offset recorded in attributes, it can happen after opening a file
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_chunkoffsets", NULL, &len);
            if (err == NC_NOERR && varp->nchunk == len - 1){
                err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_chunklens", NULL, &len);
                if (err == NC_NOERR && varp->nchunk == len){
                    err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_chunkoffsets", varp->data_offs, MPI_LONG_LONG);
                    err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_chunklens", varp->data_lens, MPI_INT);
                }
                else{
                    // If not, 0 len means no data avaiable
                    if (err != NC_NOERR){
                        memset(varp->data_offs, 0, sizeof(MPI_Offset) * varp->nchunk);
                        memset(varp->data_lens, 0, sizeof(int) * (varp->nchunk + 1));
                    }
                }
            }
            else{
                // If not, 0 len means no data avaiable
                if (err != NC_NOERR){
                    memset(varp->data_offs, 0, sizeof(MPI_Offset) * (varp->nchunk + 1));
                    memset(varp->data_lens, 0, sizeof(int) * varp->nchunk);
                }
            }

            /* Select compression driver based on attribute */
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_zipdriver", NULL, &len);
            if (err == NC_NOERR && len == 1){
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_zipdriver", &(varp->zipdriver), MPI_INT);
            }
            else{
                varp->zipdriver = 0;
            }
            switch (varp->zipdriver){
                case NC_ZIP_DRIVER_DUMMY:
                    varp->zip = nczip_dummy_inq_driver();
                break;
            }

            // Get variable id
            if (varp->isnew == 0){
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_datavarid", &(varp->datavarid), MPI_INT);
            }

            // Update max ndim and chunksize
            if (nczipp->max_ndim < varp->ndim){
                nczipp->max_ndim = varp->ndim;
            }
            if (nczipp->max_chunk_size < varp->chunksize){
                nczipp->max_chunk_size = varp->chunksize;
            }
        }   
    }

    return NC_NOERR;
}

int nczipioi_var_resize(NC_zip *nczipp, NC_zip_var *varp) {
    int i, j, err;
    int valid;
    MPI_Offset len;
    NC_zip_var *var;

    if (varp->varkind == NC_ZIP_VAR_COMPRESSED && varp->isrec){
        if (varp->dimsize[0] < nczipp->recsize){
            int oldnchunk, oldnrec;
            int chunkperrec;
            int oldnmychunk;

            oldnrec = varp->dimsize[0];
            oldnchunk = varp->nchunk;
            varp->dimsize[0] = nczipp->recsize;

            // Calculate new # chunks, # chunks along each dim, chunksize
            varp->nchunk = 1;
            for(i = 0; i < varp->ndim; i++){ //chunkdim must be at leasst 1
                if (varp->dimsize[i] % varp->chunkdim[i] == 0){
                    varp->nchunks[i] = (int)varp->dimsize[i] / varp->chunkdim[i];
                }
                else{
                    varp->nchunks[i] = (int)varp->dimsize[i] / varp->chunkdim[i]; + 1;
                }
                varp->nchunk *= varp->nchunks[i];
            }

            // Extend offset and len list
            varp->data_offs = (MPI_Offset*)NCI_Realloc(varp->data_offs, sizeof(MPI_Offset) * (varp->nchunk + 1));
            varp->data_lens = (int*)NCI_Realloc(varp->data_lens, sizeof(int) * varp->nchunk);
            memset(varp->data_offs + oldnchunk, 0, sizeof(int) * (varp->nchunk - oldnchunk));
            memset(varp->data_lens + oldnchunk, 0, sizeof(int) * (varp->nchunk - oldnchunk));

            // Extend block ownership list
            varp->chunk_owner = (int*)NCI_Realloc(varp->chunk_owner, sizeof(int) * varp->nchunk);
            varp->chunk_cache = (char**)NCI_Realloc(varp->chunk_cache, sizeof(char*) * varp->nchunk);
            if (oldnrec > 0){
                chunkperrec = oldnchunk / oldnrec;
                for(i = oldnchunk; i < varp->nchunk; i += chunkperrec){
                    // We reuse chunk mapping of other records
                    memcpy(varp->chunk_owner + i, varp->chunk_owner, sizeof(int) * chunkperrec);
                }
            }
            else{
                varp->nmychunks = 0;
                if (nczipp->blockmapping == NC_ZIP_MAPPING_STATIC){
                    for(j = 0; j < varp->nchunk; j++){ 
                        varp->chunk_owner[j] = j % nczipp->np;
                    }
                }
            }

            // Extend skip list of my own chunks
            oldnmychunk = varp->nmychunks;
            for(i = oldnchunk; i < varp->nchunk; i ++){
                if (varp->chunk_owner[i] == nczipp->rank){
                    varp->nmychunks++;
                    varp->chunk_cache[i] = (void*)NCI_Malloc(varp->chunksize);  // Allocate buffer for blocks we own
                }
            }

            if (oldnmychunk < varp->nmychunks){
                varp->mychunks = (int*)NCI_Realloc(varp->mychunks, sizeof(int) * varp->nmychunks);
                for(i = oldnchunk; i < varp->nchunk; i++){ 
                    if (varp->chunk_owner[i] == nczipp->rank){
                        varp->mychunks[oldnmychunk++] = i;
                    }
                }

                if (oldnmychunk != varp->nmychunks){
                    printf("Error\n");
                }
            }
        }
    }
    else{
        // Notify ncmpio driver
    }

    return NC_NOERR;
}

void nczipioi_var_free(NC_zip_var *varp) {
    int i;

    if (varp->chunkdim != NULL){
        NCI_Free(varp->dimsize);
        NCI_Free(varp->chunkdim);
        NCI_Free(varp->dimids);
        NCI_Free(varp->nchunks);
        NCI_Free(varp->cidsteps);
        NCI_Free(varp->data_offs);
        NCI_Free(varp->chunk_owner);
        NCI_Free(varp->data_lens);
        for(i = 0; i < varp->nmychunks; i++){
            if (varp->chunk_cache[varp->mychunks[i]] != NULL){
                NCI_Free(varp->chunk_cache[varp->mychunks[i]]);
            }
        }
        NCI_Free(varp->chunk_cache);
        NCI_Free(varp->mychunks);
    }
}

int nczipioi_load_var(NC_zip *nczipp, NC_zip_var *varp, int nchunk, int *cids) {
    int err;
    int i;
    int cid;
    int get_size;

    int dsize;
    MPI_Offset bsize;

    int *lens;
    MPI_Aint *disps;
    MPI_Status status;
    MPI_Datatype ftype;  // Memory and file datatype

    int *zsizes;
    MPI_Offset *zoffs;
    char **zbufs;

    NC *ncp = (NC*)(nczipp->ncp);
    NC_var *ncvarp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_INIT)

    // -1 means all chunks
    if (nchunk < 0){
        nchunk = varp->nmychunks;
        cids = varp->mychunks;
    }

    zsizes = varp->data_lens;
    zoffs = varp->data_offs;

    // Allocate buffer for I/O
    lens = (int*)NCI_Malloc(sizeof(int) * nchunk);
    disps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * nchunk);
    zbufs = (char**)NCI_Malloc(sizeof(char*) * nchunk);

    /* Carry our coll I/O
     * OpenMPI will fail when set view or do I/O on type created with MPI_Type_create_hindexed when count is 0
     * We use a dummy call inplace of type with 0 count
     */
    if (nchunk > 0){
        // Create file type
        ncvarp = ncp->vars.value[varp->datavarid];
        bsize = 0;
        for(i = 0; i < nchunk; i++){
            cid = cids[i];
            // offset and length of compressed chunks
            lens[i] = zsizes[cid];
            disps[i] = (MPI_Aint)zoffs[cid] + (MPI_Aint)ncvarp->begin;
            // At the same time, we record the size of buffer we need
            bsize += (MPI_Offset)lens[i];
        }
        MPI_Type_create_hindexed(nchunk, lens, disps, MPI_BYTE, &ftype);
        MPI_Type_commit(&ftype);

        // Allocate buffer for compressed data
        // We allocate it continuously so no mem type needed
        zbufs[0] = (char*)NCI_Malloc(bsize);
        for(i = 1; i < nchunk; i++){
            zbufs[i] = zbufs[i - 1] + zsizes[cids[i - 1]];
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_INIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_RD)

        // Perform MPI-IO
        // Set file view
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
        // Write data
        MPI_File_read_at_all(ncp->collective_fh, 0, zbufs[0], bsize, MPI_BYTE, &status);
        // Restore file view
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_RD)

#ifdef _USE_MPI_GET_COUNT
        MPI_Get_count(&status, MPI_BYTE, &get_size);
#else
        MPI_Type_size(ftype, &get_size);
#endif
        nczipp->getsize += get_size;

        // Free type
        MPI_Type_free(&ftype);
    }
    else{
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_INIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_RD)

        // Follow coll I/O with dummy call
        CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
        CHK_ERR_READ_AT_ALL(ncp->collective_fh, 0, &i, 0, MPI_BYTE, &status);
        CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_RD)
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_DECOM)

    // Decompress each chunk
    // Allocate chunk cache if not allocated
    dsize = varp->chunksize;
    for(i = 0; i < nchunk; i++){
        cid = cids[i];
        if (varp->chunk_cache[cid] == NULL){
            varp->chunk_cache[cid] = (char*)NCI_Malloc(varp->chunksize);
        }
        
        varp->zip->decompress(zbufs[i], zsizes[i], varp->chunk_cache[cid], &dsize, varp->ndim, varp->chunkdim, varp->etype);

        if(dsize != varp->chunksize){
            printf("Decompress Error\n");
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_DECOM)

    // Free buffers
    if (nchunk > 0){
        NCI_Free(zbufs[0]);
    }
    NCI_Free(zbufs);

    NCI_Free(lens);
    NCI_Free(disps);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO)

    return NC_NOERR;
}

int nczipioi_load_nvar(NC_zip *nczipp, int nvar, int *varids) {
    int err;
    int i, j, k;
    int cid, vid;
    int get_size;

    int nchunk;

    int dsize;
    MPI_Offset bsize;

    int *lens;
    MPI_Aint *disps;
    MPI_Status status;
    MPI_Datatype ftype;  // Memory and file datatype

    char **zbufs;

    NC *ncp = (NC*)(nczipp->ncp);
    NC_zip_var *varp;
    NC_var *ncvarp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_INIT)

    // -1 means all chunks
    nchunk = 0;
    for(i = 0; i < nvar; i++){
        varp = nczipp->vars.data + varids[i];
    
        for(j = 0; j < varp->nmychunks; j++){
            cid = varp->mychunks[j];

            // We only need to read when it is not in cache
            if (varp->chunk_cache[cid] == NULL){
                nchunk++;
            }
        }
    }

    // Allocate buffer for I/O
    lens = (int*)NCI_Malloc(sizeof(int) * nchunk);
    disps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * nchunk);
    zbufs = (char**)NCI_Malloc(sizeof(char*) * nchunk);

    /* Carry our coll I/O
     * OpenMPI will fail when set view or do I/O on type created with MPI_Type_create_hindexed when count is 0
     * We use a dummy call inplace of type with 0 count
     */
    if (nchunk > 0){
        // Create file type
        bsize = 0;
        k = 0;
        for(i = 0; i < nvar; i++){
            varp = nczipp->vars.data + varids[i];
            ncvarp = ncp->vars.value[varp->datavarid];
        
            for(j = 0; j < varp->nmychunks; j++){
                cid = varp->mychunks[j];

                // We only need to read when it is not in cache
                if (varp->chunk_cache[cid] == NULL){
                    // offset and length of compressed chunks
                    lens[k] = varp->data_lens[cid];
                    disps[k] = (MPI_Aint)(varp->data_offs[cid]) + (MPI_Aint)ncvarp->begin;
                    // At the same time, we record the size of buffer we need
                    bsize += (MPI_Offset)lens[k++];
                }
            }
        }

        MPI_Type_create_hindexed(nchunk, lens, disps, MPI_BYTE, &ftype);
        MPI_Type_commit(&ftype);

        // Allocate buffer for compressed data
        // We allocate it continuously so no mem type needed
        zbufs[0] = (char*)NCI_Malloc(bsize);
        for(j = 1; j < nchunk; j++){
            zbufs[j] = zbufs[j - 1] + lens[j - 1];
        }    

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_INIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_RD)

        // Perform MPI-IO
        // Set file view
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
        // Write data
        MPI_File_read_at_all(ncp->collective_fh, 0, zbufs[0], bsize, MPI_BYTE, &status);
        // Restore file view
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_RD)

#ifdef _USE_MPI_GET_COUNT
        MPI_Get_count(&status, MPI_BYTE, &get_size);
#else
        MPI_Type_size(ftype, &get_size);
#endif
        nczipp->getsize += get_size;

        // Free type
        MPI_Type_free(&ftype);

    }
    else{
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_RD)

        // Follow coll I/O with dummy call
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
        MPI_File_read_at_all(ncp->collective_fh, 0, &i, 0, MPI_BYTE, &status);
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_RD)
    }
    
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_DECOM)

    // Decompress each chunk
    k = 0;
    for(i = 0; i < nvar; i++){
        varp = nczipp->vars.data + varids[i];
        dsize = varp->chunksize;

        for(j = 0; j < varp->nmychunks; j++){
            cid = varp->mychunks[j];

            // Allocate chunk cache if not allocated
            if (varp->chunk_cache[cid] == NULL){
                varp->chunk_cache[cid] = (char*)NCI_Malloc(varp->chunksize);
            }

            // Perform decompression
            varp->zip->decompress(zbufs[k], lens[k], varp->chunk_cache[cid], &dsize, varp->ndim, varp->chunkdim, varp->etype);
            if(dsize != varp->chunksize){
                printf("Decompress Error\n");
            }

            k++;
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_DECOM)

    // Free buffers
    if (nchunk > 0){
        NCI_Free(zbufs[0]);
    }
    NCI_Free(zbufs);

    NCI_Free(lens);
    NCI_Free(disps);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO)

    return NC_NOERR;
}


int nczipioi_save_var(NC_zip *nczipp, NC_zip_var *varp) {
    int i, j, k, l, err;
    int *zsizes, *zsizes_all;
    MPI_Datatype mtype, ftype;  // Memory and file datatype
    int wcnt;
    int *lens;
    MPI_Aint *disps;
    MPI_Status status;
    MPI_Offset *zoffs;
    void **zbufs;
    int zdimid;
    int put_size;
    char name[128]; // Name of objects
    NC *ncp = (NC*)(nczipp->ncp);
    NC_var *ncvarp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO)

    // Allocate buffer for compression
    zsizes = (int*)NCI_Malloc(sizeof(int) * varp->nchunk);
    zbufs = (void**)NCI_Malloc(sizeof(void*) * varp->nmychunks);
    zsizes_all = varp->data_lens;
    zoffs = varp->data_offs;

    // Allocate buffer for I/O
    wcnt = varp->nmychunks;
    lens = (int*)NCI_Malloc(sizeof(int) * wcnt);
    disps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);

    memset(zsizes, 0, sizeof(int) * varp->nchunk);

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_COMP)

    // Compress each chunk we own
    memset(zsizes, 0, sizeof(int) * varp->nchunk);
    for(l = 0; l < varp->nmychunks; l++){
        k = varp->mychunks[l];

        // Apply compression
        varp->zip->compress_alloc(varp->chunk_cache[k], varp->chunksize, zbufs + l, zsizes + k, varp->ndim, varp->chunkdim, varp->etype);

        // Record compressed size
        lens[l] = zsizes[k];
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_COMP)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_SYNC)

    // Sync compressed data size with other processes
    MPI_Allreduce(zsizes, zsizes_all, varp->nchunk, MPI_INT, MPI_MAX, nczipp->comm);
    zoffs[0] = 0;
    for(i = 0; i < varp->nchunk; i++){
        zoffs[i + 1] = zoffs[i] + zsizes_all[i];
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_SYNC)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_INIT)

    /* Write comrpessed variable
     * We start by defining data variable and writing metadata
     * Then, we create buffer type and file type for data
     * Finally MPI collective I/O is used for writing data
     */

    // Enter redefine mode
    nczipp->driver->redef(nczipp->ncp);

    // Define dimension for data variable
    sprintf(name, "_compressed_data_dim_%d", varp->varid);
    err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[varp->nchunk], &zdimid);
    if (err != NC_NOERR) return err;

    // Define data variable
    sprintf(name, "_compressed_data_%d", varp->varid);
    err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));
    if (err != NC_NOERR) return err;

    // Mark as data variable
    i = NC_ZIP_VAR_DATA;
    err = nczipp->driver->put_att(nczipp->ncp, varp->datavarid, "_varkind", NC_INT, 1, &i, MPI_INT);
    if (err != NC_NOERR) return err;

    // Record offset of chunks in data variable
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkoffsets", NC_INT64, varp->nchunk + 1, zoffs, MPI_LONG_LONG);
    if (err != NC_NOERR) return err;

    // Record size of chunks
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunklens", NC_INT, varp->nchunk, zsizes_all, MPI_INT);
    if (err != NC_NOERR) return err;

    // Record data variable id
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_datavarid", NC_INT, 1, &(varp->datavarid), MPI_INT);
    if (err != NC_NOERR) return err;

    // Switch to data mode
    err = nczipp->driver->enddef(nczipp->ncp);
    if (err != NC_NOERR) return err;

    /* Carry our coll I/O
     * OpenMPI will fail when set view or do I/O on type created with MPI_Type_create_hindexed when count is 0
     * We use a dummy call inplace of type with 0 count
     */
    if (wcnt > 0){
        // Create file type
        ncvarp = ncp->vars.value[varp->datavarid];
        for(l = 0; l < varp->nmychunks; l++){
            k = varp->mychunks[l];

            // Record compressed size
            lens[l] = zsizes[k];
            disps[l] = (MPI_Aint)zoffs[k] + (MPI_Aint)ncvarp->begin;
        }
        MPI_Type_create_hindexed(wcnt, lens, disps, MPI_BYTE, &ftype);
        MPI_Type_commit(&ftype);

        // Create memory buffer type
        for(l = 0; l < varp->nmychunks; l++){
            k = varp->mychunks[l];

            // Record compressed size
            lens[l] = zsizes[k];
            disps[l] = (MPI_Aint)zbufs[l];
        }
        err = MPI_Type_create_hindexed(wcnt, lens, disps, MPI_BYTE, &mtype);
        MPI_Type_commit(&mtype);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_INIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_WR)

        // Perform MPI-IO
        // Set file view
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
        // Write data
        MPI_File_write_at_all(ncp->collective_fh, 0, MPI_BOTTOM, 1, mtype, &status);
        // Restore file view
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_WR)

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
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_INIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_WR)

        // Follow coll I/O with dummy call
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
        MPI_File_write_at_all(ncp->collective_fh, 0, MPI_BOTTOM, 0, MPI_BYTE, &status);
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_WR)
    }

    // Free buffers
    NCI_Free(zsizes);
    for(l = 0; l < varp->nmychunks; l++){
        NCI_Free(zbufs[l]);
    }
    NCI_Free(zbufs);

    NCI_Free(lens);
    NCI_Free(disps);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO)

    return NC_NOERR;
}

int nczipioi_save_nvar(NC_zip *nczipp, int nvar, int *varids) {
    NC_zip_var *varp;
    int i, j, k, l, err;
    int vid;    // Iterator for variable id
    int cid;    // Iterator for chunk id
    int max_nchunks = 0;
    int *zsizes, *zsizes_all;
    MPI_Offset *zoffs;
    MPI_Datatype mtype, ftype;  // Memory and file datatype
    int wcnt, wcur;
    int *mlens, *flens;
    MPI_Aint *mdisps, *fdisps;
    MPI_Status status;
    int put_size;
    void **zbufs;
    int zdimid;
    char name[128]; // Name of objects
    NC *ncp = (NC*)(nczipp->ncp);
    NC_var *ncvarp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_INIT)

    wcnt = 0;
    for(i = 0; i < nvar; i++){
        vid = varids[i];
        wcnt += nczipp->vars.data[vid].nmychunks;
        if (max_nchunks < nczipp->vars.data[vid].nchunk){
            max_nchunks = nczipp->vars.data[vid].nchunk;
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_INIT)

    // Allocate buffer for compression
    zsizes = (int*)NCI_Malloc(sizeof(int) * max_nchunks);
    zbufs = (void**)NCI_Malloc(sizeof(void*) * wcnt);

    // Allocate buffer file type
    mlens = (int*)NCI_Malloc(sizeof(int) * wcnt);
    mdisps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);
    flens = (int*)NCI_Malloc(sizeof(int) * wcnt);
    fdisps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);

    // Enter redefine mode
    nczipp->driver->redef(nczipp->ncp);

    wcur = 0;
    for(vid = 0; vid < nvar; vid++){
        varp = nczipp->vars.data + varids[vid];

        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_COMP)

        zsizes_all = varp->data_lens;
        zoffs = varp->data_offs;

        memset(zsizes, 0, sizeof(int) * varp->nchunk);

        // Compress each chunk we own
        memset(zsizes, 0, sizeof(int) * varp->nchunk);
        for(l = 0; l < varp->nmychunks; l++){
            cid = varp->mychunks[l];

            // Apply compression
            varp->zip->compress_alloc(varp->chunk_cache[cid], varp->chunksize, zbufs + wcur + l, zsizes + cid, varp->ndim, varp->chunkdim, varp->etype);
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_COMP)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_SYNC)

        // Sync compressed data size with other processes
        MPI_Allreduce(zsizes, zsizes_all, varp->nchunk, MPI_INT, MPI_MAX, nczipp->comm);
        zoffs[0] = 0;
        for(cid = 0; cid < varp->nchunk; cid++){
            zoffs[cid + 1] = zoffs[cid] + zsizes_all[cid];
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_SYNC)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_INIT)

        /* Write comrpessed variable
        * We start by defining data variable and writing metadata
        * Then, we create buffer type and file type for data
        * Finally MPI collective I/O is used for writing data
        */

        // Define dimension for data variable
        sprintf(name, "_compressed_data_dim_%d", varp->varid);
        err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[varp->nchunk], &zdimid);
        if (err != NC_NOERR) return err;

        // Define data variable
        sprintf(name, "_compressed_data_%d", varp->varid);
        err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));
        if (err != NC_NOERR) return err;

        // Mark as data variable
        i = NC_ZIP_VAR_DATA;
        err = nczipp->driver->put_att(nczipp->ncp, varp->datavarid, "_varkind", NC_INT, 1, &i, MPI_INT);
        if (err != NC_NOERR) return err;

        // Record offset of chunks in data variable
        err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkoffsets", NC_INT64, varp->nchunk + 1, zoffs, MPI_LONG_LONG);
        if (err != NC_NOERR) return err;

        // Record size of chunks
        err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunklens", NC_INT, varp->nchunk, zsizes_all, MPI_INT);
        if (err != NC_NOERR) return err;

        // Record data variable id
        err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_datavarid", NC_INT, 1, &(varp->datavarid), MPI_INT);
        if (err != NC_NOERR) return err;

        /* Paramemter for file and memory type 
         * We do not know variable file offset until the end of define mode
         * We will add the displacement later
         */
        for(l = 0; l < varp->nmychunks; l++){
            cid = varp->mychunks[l];

            // Record parameter
            flens[wcur + l] = zsizes[cid];
            fdisps[wcur + l] = (MPI_Aint)zoffs[cid];
            mlens[wcur + l] = zsizes[cid];
            mdisps[wcur + l] = (MPI_Aint)zbufs[wcur + l];
        }

        // Move to parameters for next variable
        wcur += varp->nmychunks;

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_INIT)
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_INIT)

    // Switch back to data mode
    err = nczipp->driver->enddef(nczipp->ncp);
    if (err != NC_NOERR) return err;

    /* Now it's time to add variable file offset to displacements
     * File type offset need to be specified in non-decreasing order
     * We assume ncmpio place variable according to the order they are declared
     */
    wcur = 0;
    for(vid = 0; vid < nvar; vid++){
        varp = nczipp->vars.data +  varids[vid];
        ncvarp = ncp->vars.value[varp->datavarid];
        for(l = 0; l < varp->nmychunks; l++){
            cid = varp->mychunks[l];
            // Adjust file displacement
            fdisps[wcur++] += (MPI_Aint)ncvarp->begin;
        }
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_INIT)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_IO_WR)

    /* Carry our coll I/O
     * OpenMPI will fail when set view or do I/O on type created with MPI_Type_create_hindexed when count is 0
     * We use a dummy call inplace of type with 0 count
     */
    if (wcnt > 0){
         // Create file type
        MPI_Type_create_hindexed(wcnt, flens, fdisps, MPI_BYTE, &ftype);
        MPI_Type_commit(&ftype);

        // Create memmory type
        MPI_Type_create_hindexed(wcnt, mlens, mdisps, MPI_BYTE, &mtype);
        MPI_Type_commit(&mtype);

        // Perform MPI-IO
        // Set file view
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
        // Write data
        MPI_File_write_at_all(ncp->collective_fh, 0, MPI_BOTTOM, 1, mtype, &status);
        // Restore file view
        MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_WR)

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

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO_WR)
    }

    // Free buffers
    NCI_Free(zsizes);
    for(l = 0; l < varp->nmychunks; l++){
        NCI_Free(zbufs[l]);
    }
    NCI_Free(zbufs);

    NCI_Free(flens);
    NCI_Free(fdisps);
    NCI_Free(mlens);
    NCI_Free(mdisps);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_IO)

    return NC_NOERR;
}

int nczipioi_init_nvar(NC_zip *nczipp, int nput, int *putreqs, int nget, int *getreqs){
    int err;
    int i, j;
    int nflag;
    unsigned int *flag, *flag_all;
    int nvar;
    int *rcnt, *roff;
    int *vids, *vmap;
    MPI_Offset **starts, **counts;
    NC_zip_req *req;

    CHK_ERR_ALLREDUCE(MPI_IN_PLACE, &(nczipp->recsize), 1, MPI_LONG_LONG, MPI_MAX, nczipp->comm);   // Sync number of recs

    // Flag of touched vars
    nflag = nczipp->vars.cnt / 32 + 1;
    flag = (unsigned int*)NCI_Malloc(sizeof(int) * nflag);
    flag_all = (unsigned int*)NCI_Malloc(sizeof(int) * nflag);
    memset(flag, 0, sizeof(int) * nflag);
    for(i = 0; i < nput; i++){
        req = nczipp->putlist.reqs + putreqs[i];
        flag[req->varid >> 5] |= 1u << (req->varid % 32);
    }
    for(i = 0; i < nget; i++){
        req = nczipp->getlist.reqs + getreqs[i];
        flag[req->varid >> 5] |= 1u << (req->varid % 32);
    }

    // Sync flag
    CHK_ERR_ALLREDUCE(flag, flag_all, nflag, MPI_UNSIGNED, MPI_BOR, nczipp->comm);

    // Build a skip list of touched vars
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (flag_all[i >> 5] & (1u << (i % 32))) {
            if ((nczipp->vars.data + i)->chunkdim == NULL){   // If not yet inited
                nvar++;
            }
            else{   
                flag_all[i >> 5] ^= (1u << (i % 32));
                if ((nczipp->vars.data + i)->dimsize[0] < nczipp->recsize){
                    nczipioi_var_resize(nczipp, nczipp->vars.data + i);
                }
            }
        }
    }
    vids = (int*)NCI_Malloc(sizeof(int) * nvar);
    vmap = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt);
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (flag_all[i >> 5] & (1u << (i % 32))) {
            vids[nvar] = i;
            vmap[i] = nvar++;
        }
    }

    // Count reqs for each var
    roff = (int*)NCI_Malloc(sizeof(int) * (nvar + 1));
    rcnt = (int*)NCI_Malloc(sizeof(int) * nvar);
    memset(rcnt, 0, sizeof(int) * nvar);
    for(i = 0; i < nput; i++){
        req = nczipp->putlist.reqs + putreqs[i];
        j = req->varid;
        if (flag_all[j >> 5] & (1u << (j % 32))) {
            rcnt[vmap[j]] += req->nreq;
        }
    }
    for(i = 0; i < nget; i++){
        req = nczipp->getlist.reqs + getreqs[i];
        j = req->varid;
        if (flag_all[j >> 5] & (1u << (j % 32))) {
            rcnt[vmap[j]] += req->nreq;
        }
    }
    roff[0] = 0;
    for(i = 0; i < nvar; i++){
        roff[i + 1] = roff[i] + rcnt[i];
    }

    // Gather starts and counts
    starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * roff[nvar] * 2);
    counts = starts + roff[nvar];
    memset(rcnt, 0, sizeof(int) * nvar);
    for(i = 0; i < nput; i++){
        req = nczipp->putlist.reqs + putreqs[i];
        j = req->varid;
        if (flag_all[j >> 5] & (1u << (j % 32))) {
            j = vmap[req->varid];
            if (req->nreq > 1){
                memcpy(starts + roff[j] + rcnt[j], req->starts, sizeof(MPI_Offset*) * req->nreq);
                memcpy(counts + roff[j] + rcnt[j], req->counts, sizeof(MPI_Offset*) * req->nreq);
                rcnt[j] += req->nreq;
            }
            else{
                starts[roff[j] + rcnt[j]] = req->start;
                counts[roff[j] + (rcnt[j]++)] = req->count;     
            }
        }
    }
    for(i = 0; i < nget; i++){
        req = nczipp->getlist.reqs + getreqs[i];
        j = req->varid;
        if (flag_all[j >> 5] & (1u << (j % 32))) {
            j = vmap[req->varid];
            if (req->nreq > 1){
                memcpy(starts + roff[j] + rcnt[j], req->starts, sizeof(MPI_Offset*) * req->nreq);
                memcpy(counts + roff[j] + rcnt[j], req->counts, sizeof(MPI_Offset*) * req->nreq);
                rcnt[j] += req->nreq;
            }
            else{
                starts[roff[j] + rcnt[j]] = req->start;
                counts[roff[j] + (rcnt[j]++)] = req->count;     
            }
        }
    }

    for(i = 0; i < nvar; i++){
        err = nczipioi_var_init(nczipp, nczipp->vars.data + vids[i], rcnt[i], starts + roff[i], counts + roff[i]);
        if (err != NC_NOERR){
            return err;
        }
    }

    NCI_Free(flag);
    NCI_Free(vids);
    NCI_Free(vmap);
    NCI_Free(roff);
    NCI_Free(rcnt);
    NCI_Free(starts);

    return NC_NOERR;
}