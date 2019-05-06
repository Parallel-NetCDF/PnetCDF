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
            err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_offvarid", &(varp->offvarid), MPI_LONG_LONG);
            err |= nczipp->driver->get_att(nczipp->ncp, varp->varid, "_lenvarid", &(varp->lenvarid), MPI_INT);
            if (err == NC_NOERR){
                MPI_Offset start, count;
                
                start = 0;
                count = varp->nchunk + 1;
                err = nczipp->driver->get_var(nczipp->ncp, varp->offvarid, &start, &count, NULL, NULL, varp->data_offs, -1, MPI_LONG_LONG, NC_REQ_RD | NC_REQ_BLK | NC_REQ_HL | NC_REQ_COLL);
                if (err != NC_NOERR) return err;

                count = varp->nchunk;
                err = nczipp->driver->get_var(nczipp->ncp, varp->lenvarid, &start, &count, NULL, NULL, varp->data_lens, -1, MPI_INT, NC_REQ_RD | NC_REQ_BLK | NC_REQ_HL | NC_REQ_COLL);
                if (err != NC_NOERR) return err;
            }
            else {
                varp->offvarid = varp->lenvarid = -1;
                memset(varp->data_offs, 0, sizeof(MPI_Offset) * (varp->nchunk + 1));
                memset(varp->data_lens, 0, sizeof(int) * varp->nchunk);
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

            // Get metadata, they may not exist so we don't catch error
            if (varp->isnew == 0){
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_datavarid", &(varp->datavarid), MPI_INT);
                if (err != NC_NOERR){
                    varp->datavarid = -1;
                }
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_dataserial", &(varp->dataserial), MPI_INT);
                err |= nczipp->driver->get_att(nczipp->ncp, varp->varid, "_metaserial", &(varp->metaserial), MPI_INT);
                if (err != NC_NOERR){
                    varp->dataserial = varp->metaserial = 0;
                }
            }
            else{
                varp->datavarid = -1;
                varp->dataserial = varp->metaserial = 0;
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

            varp->expanded = 1;

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
    flag = (unsigned int*)NCI_Malloc(sizeof(int) * nflag * 2);
    flag_all = flag + nflag;
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

int nczipioi_resize_nvar(NC_zip *nczipp, int nput, int *putreqs, int nget, int *getreqs){
    int err;
    int i;
    int nflag;
    unsigned int *flag, *flag_all;
    int nvar;
    int *vids;
    NC_zip_req *req;
    NC_zip_var *varp;

    CHK_ERR_ALLREDUCE(MPI_IN_PLACE, &(nczipp->recsize), 1, MPI_LONG_LONG, MPI_MAX, nczipp->comm);   // Sync number of recs

    // Flag of touched vars
    nflag = nczipp->vars.cnt / 32 + 1;
    flag = (unsigned int*)NCI_Malloc(sizeof(int) * nflag * 2);
    flag_all = flag + nflag;
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
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (flag_all[i >> 5] & (1u << (i % 32))) {
            vids[nvar] = i;
        }
    }

    // Count reqs for each var
    for(i = 0; i < nvar; i++){
        varp = nczipp->vars.data + vids[i];
        err = nczipioi_var_resize(nczipp, varp);
        if (err != NC_NOERR){
            return err;
        }
    }

    NCI_Free(flag);
    NCI_Free(vids);

    return NC_NOERR;
}