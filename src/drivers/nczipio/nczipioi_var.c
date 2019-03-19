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

int nczipioi_var_init(NC_zip *nczipp, NC_zip_var *varp) {
    int i, j, err;
    int valid;
    MPI_Offset len;
    NC_zip_var *var;

    if (varp->varkind == NC_ZIP_VAR_COMPRESSED){
        if (varp->chunkdim == NULL){    // This is a new uninitialized variable 
            // Determine its block size
            varp->chunkdim = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
            
            // First check attribute
            valid = 1;
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_chunkdim", NULL, &len);
            if (err == NC_NOERR && len == varp->ndim){
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_chunkdim", varp->chunkdim, MPI_UNSIGNED_LONG_LONG);
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
                for(j = 0; j < varp->ndim; j++){ 
                    varp->chunkdim[j] = (int)varp->dimsize[j];
                }
                err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkdim", NC_INT, varp->ndim, varp->chunkdim, MPI_INT);
                if (err != NC_NOERR){
                    return err;
                }
            }

            // Calculate block size
            varp->chunksize = NC_Type_size(varp->xtype);
            for(i = 0; i < varp->ndim; i++){
                varp->chunksize *= varp->chunkdim[i];
            }

            // Calculate number of blocks
            varp->nchunks = 1;
            len = 1;
            for(j = 0; j < varp->ndim; j++){ //chunkdim must be at leasst 1
                if (varp->chunkdim[j] % varp->chunkdim[j] == 0){
                    varp->nchunks *= varp->chunkdim[j] / varp->chunkdim[j];
                }
                else{
                    varp->nchunks *= varp->chunkdim[j] / varp->chunkdim[j] + 1;
                }
                len *= varp->chunkdim[j];   // Block size
            }

            // Determine block ownership
            varp->chunk_owner = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
            varp->chunk_cache = (char**)NCI_Malloc(sizeof(char*) * varp->nchunks);
            memset(varp->chunk_cache, 0, sizeof(char*) * varp->nchunks);
            varp->nmychunks = 0;
            if (nczipp->blockmapping == NC_ZIP_MAPPING_STATIC){
                for(j = 0; j < varp->nchunks; j++){ 
                    varp->chunk_owner[j] = j % nczipp->np;
                    if (varp->chunk_owner[j] == nczipp->rank){
                        varp->chunk_cache[j] = (void*)NCI_Malloc(varp->chunksize);  // Allocate buffer for blocks we own
                    }
                    varp->nmychunks++;
                }
            }

            // Build skip list of my own chunks
            varp->mychunks = (int*)NCI_Malloc(sizeof(int) * varp->nmychunks);
            varp->nmychunks = 0;
            for(j = 0; j < varp->nchunks; j++){ 
                if (varp->chunk_owner[j] == nczipp->rank){
                    varp->mychunks[varp->nmychunks++] = j;
                }
            }

            // Determine block offset
            varp->data_offs = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->nchunks);
            varp->data_lens = (int*)NCI_Malloc(sizeof(int) * (varp->nchunks + 1));
            // Try if there are offset recorded in attributes, it can happen after opening a file
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_chunkoffsets", NULL, &len);
            if (err == NC_NOERR && varp->nchunks == len - 1){
                err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_chunklens", NULL, &len);
                if (err == NC_NOERR && varp->nchunks == len - 1){
                    err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_chunkoffsets", varp->data_offs, MPI_UNSIGNED_LONG_LONG);
                    err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_chunklens", varp->data_lens, MPI_INT);
                }
            }

            // If not, 0 len means no data avaiable
            if (err != NC_NOERR){
                memset(varp->data_offs, 0, sizeof(MPI_Offset) * varp->nchunks);
                memset(varp->data_lens, 0, sizeof(MPI_Offset) * varp->nchunks);
            }
        }   
    }

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
    char **zbufs;
    int zdimid;
    char name[128]; // Name of objects
    NC *ncp = (NC*)(nczipp->ncp);
    NC_var *ncvarp;

    // Allocate buffer for compression
    zsizes = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
    zbufs = (char**)NCI_Malloc(sizeof(char*) * varp->nmychunks);
    zsizes_all = varp->data_lens;
    zoffs = varp->data_offs;

    // Allocate buffer for I/O
    wcnt = varp->nmychunks;
    lens = (int*)NCI_Malloc(sizeof(int) * wcnt);
    disps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);

    memset(zsizes, 0, sizeof(int) * varp->nchunks);

    // Compress each chunk we own
    memset(zsizes, 0, sizeof(int) * varp->nchunks);
    for(l = 0; l < varp->nmychunks; l++){
        k = varp->mychunks[l];

        // Apply compression
        nczipp->zip->compress_alloc(varp->chunk_cache[k], varp->chunksize, zbufs + l, zsizes + k, varp->ndim, varp->chunkdim, varp->etype);

        // Record compressed size
        lens[l] = zsizes[k];
    }

    // Sync compressed data size with other processes
    MPI_Allreduce(zsizes, zsizes_all, varp->nchunks, MPI_INT, MPI_MAX, nczipp->comm);
    zoffs[0] = 0;
    for(i = 0; i < varp->nchunks; i++){
        zoffs[i + 1] = zoffs[i] + zsizes_all[i];
    }

    /* Write comrpessed variable
     * We start by defining data variable and writing metadata
     * Then, we create buffer type and file type for data
     * Finally MPI collective I/O is used for writing data
     */

    // Enter redefine mode
    nczipp->driver->redef(nczipp->ncp);

    // Define dimension for data variable
    sprintf(name, "_compressed_data_dim_%d", varp->varid);
    err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[varp->nchunks], &zdimid);
    if (err != NC_NOERR) return err;

    // Define data variable
    sprintf(name, "_compressed_data_%d", varp->varid);
    err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));
    if (err != NC_NOERR) return err;

    // Record offset of chunks in data variable
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkoffsets", NC_INT, varp->nchunks, zoffs, MPI_INT);
    if (err != NC_NOERR) return err;

    // Record size of chunks
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunklens", NC_INT, varp->nchunks, zsizes_all, MPI_INT);
    if (err != NC_NOERR) return err;

    // Switch to data mode
    err = nczipp->driver->enddef(nczipp->ncp);
    if (err != NC_NOERR) return err;

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
    MPI_Type_create_hindexed(wcnt, lens, disps, MPI_BYTE, &mtype);
    MPI_Type_commit(&mtype);

    // Perform MPI-IO
    // Set file view
    MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
    // Write data
    MPI_File_write_at_all(ncp->collective_fh, 0, MPI_BOTTOM, 1, mtype, &status);
    // Restore file view
    MPI_File_set_view(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    // Free type
    MPI_Type_commit(&ftype);
    MPI_Type_commit(&mtype);

    // Free buffers
    NCI_Free(zsizes);
    for(l = 0; l < varp->nmychunks; l++){
        NCI_Free(zbufs[l]);
    }
    NCI_Free(zbufs);

    NCI_Free(lens);
    NCI_Free(disps);

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
    char **zbufs;
    int zdimid;
    char name[128]; // Name of objects
    NC *ncp = (NC*)(nczipp->ncp);
    NC_var *ncvarp;

    wcnt = 0;
    for(vid = 0; vid < nvar; vid++){
        wcnt += nczipp->vars.data[vid].nmychunks;
        if (max_nchunks < nczipp->vars.data[vid].nchunks){
            max_nchunks = nczipp->vars.data[vid].nchunks;
        }
    }

    // Allocate buffer for compression
    zsizes = (int*)NCI_Malloc(sizeof(int) * max_nchunks);
    zbufs = (char**)NCI_Malloc(sizeof(char*) * wcnt);

    // Allocate buffer file type
    mlens = (int*)NCI_Malloc(sizeof(int) * wcnt);
    mdisps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);
    flens = (int*)NCI_Malloc(sizeof(int) * wcnt);
    fdisps = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * wcnt);

    // Enter redefine mode
    nczipp->driver->redef(nczipp->ncp);

    wcur = 0;
    for(vid = 0; vid < nvar; vid++){
        varp = nczipp->vars.data + vid;

        zsizes_all = varp->data_lens;
        zoffs = varp->data_offs;

        memset(zsizes, 0, sizeof(int) * varp->nchunks);

        // Compress each chunk we own
        memset(zsizes, 0, sizeof(int) * varp->nchunks);
        for(l = 0; l < varp->nmychunks; l++){
            cid = varp->mychunks[l];

            // Apply compression
            nczipp->zip->compress_alloc(varp->chunk_cache[cid], varp->chunksize, zbufs + wcur + l, zsizes + cid, varp->ndim, varp->chunkdim, varp->etype);
        }

        // Sync compressed data size with other processes
        MPI_Allreduce(zsizes, zsizes_all, varp->nchunks, MPI_INT, MPI_MAX, nczipp->comm);
        zoffs[0] = 0;
        for(cid = 0; cid < varp->nchunks; cid++){
            zoffs[cid + 1] = zoffs[cid] + zsizes_all[cid];
        }

        /* Write comrpessed variable
        * We start by defining data variable and writing metadata
        * Then, we create buffer type and file type for data
        * Finally MPI collective I/O is used for writing data
        */

        // Define dimension for data variable
        sprintf(name, "_compressed_data_dim_%d", varp->varid);
        err = nczipp->driver->def_dim(nczipp->ncp, name, zoffs[varp->nchunks], &zdimid);
        if (err != NC_NOERR) return err;

        // Define data variable
        sprintf(name, "_compressed_data_%d", varp->varid);
        err = nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));
        if (err != NC_NOERR) return err;

        // Record offset of chunks in data variable
        err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkoffsets", NC_INT, varp->nchunks, zoffs, MPI_INT);
        if (err != NC_NOERR) return err;

        // Record size of chunks
        err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunklens", NC_INT, varp->nchunks, zsizes_all, MPI_INT);
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
            mlens[l] = zsizes[cid];
            mdisps[l] = (MPI_Aint)zbufs[wcur + l];
        }

        // Move to parameters for next variable
        wcur += varp->nmychunks;
    }

    // Switch back to data mode
    err = nczipp->driver->enddef(nczipp->ncp);
    if (err != NC_NOERR) return err;

    /* Now it's time to add variable file offset to displacements
     * File type offset need to be specified in non-decreasing order
     * We assume ncmpio place variable according to the order they are declared
     */
    wcur = 0;
    for(vid = 0; vid < nvar; vid++){
        varp = nczipp->vars.data + vid;
        ncvarp = ncp->vars.value[varp->datavarid];
        for(l = 0; l < varp->nmychunks; l++){
            cid = varp->mychunks[l];
            // Adjust file displacement
            fdisps[wcur++] += (MPI_Aint)ncvarp->begin;
        }
    }

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

    // Free type
    MPI_Type_commit(&ftype);
    MPI_Type_commit(&mtype);

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

    return NC_NOERR;
}