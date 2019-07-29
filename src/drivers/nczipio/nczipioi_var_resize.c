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

int nczipioi_var_resize(NC_zip *nczipp, NC_zip_var *varp) {
    int i, j, err;
    int cid;
    int valid;
    MPI_Offset len;
    NC_zip_var *var;

    if (varp->varkind == NC_ZIP_VAR_COMPRESSED && varp->isrec){
        if (varp->dimsize[0] < nczipp->recsize){
            int oldnchunk, oldnrec;
            int chunkperrec;
            int oldnmychunk;

            oldnrec = varp->nrec;
            oldnchunk = varp->nchunk;
            oldnmychunk = varp->nmychunk;
            varp->nrec = varp->dimsize[0] = varp->nchunks[0] = nczipp->recsize;
            varp->nchunk = varp->nchunkrec * varp->nrec;

            // Extend metadata list if needed
            if (varp->nrec > varp->nrecalloc){
                while(varp->nrecalloc < varp->nrec){
                    varp->nrecalloc *= NC_ZIP_REC_MULTIPLIER;
                }
                varp->nchunkalloc = varp->nrecalloc * varp->nchunkrec;

                varp->chunk_owner = (int*)NCI_Realloc(varp->chunk_owner, sizeof(int) * varp->nchunkalloc);
                varp->dirty = (int*)NCI_Realloc(varp->dirty, sizeof(int) * varp->nchunkalloc);
                varp->chunk_cache = (char**)NCI_Realloc(varp->chunk_cache, sizeof(char*) * varp->nchunkalloc);
                for(i = 0; i < oldnmychunk; i++){
                    cid = varp->mychunks[i];
                    if (varp->chunk_cache[cid] != NULL){
                        varp->chunk_cache[cid]->ref = varp->chunk_cache + cid;
                    }
                }

                varp->data_offs = (MPI_Offset*)NCI_Realloc(varp->data_offs, sizeof(MPI_Offset) * (varp->nchunkalloc + 1));
                varp->data_lens = (int*)NCI_Realloc(varp->data_lens, sizeof(int) * (varp->nchunkalloc + 1));

                varp->mychunks = (int*)NCI_Realloc(varp->mychunks, sizeof(int) * varp->nrecalloc * varp->nmychunkrec);

                varp->expanded = 1;
            }
            memset(varp->data_offs + oldnchunk, 0, sizeof(int) * (varp->nchunk - oldnchunk));
            memset(varp->data_lens + oldnchunk, 0, sizeof(int) * (varp->nchunk - oldnchunk));
            memset(varp->dirty + oldnchunk, 0, sizeof(int) * (varp->nchunk - oldnchunk));
            memset(varp->chunk_cache + oldnchunk, 0, sizeof(char*) * (varp->nchunk - oldnchunk));

            // Extend block ownership list
            if (oldnchunk > 0){
                for(i = oldnchunk; i < varp->nchunk; i += varp->nchunkrec){
                    // We reuse chunk mapping of other records
                    memcpy(varp->chunk_owner + i, varp->chunk_owner, sizeof(int) * varp->nchunkrec);
                }
                varp->nmychunk = varp->nmychunkrec * varp->nrec;
                for(i = oldnmychunk; i < varp->nmychunk; i += varp->nmychunkrec){
                    // We reuse chunk mapping of other records
                    memcpy(varp->mychunks + i, varp->mychunks, sizeof(int) * varp->nmychunkrec);
                }
            }
            else{
                err = nczipioi_calc_chunk_owner(nczipp, varp, 0, NULL, NULL); CHK_ERR

                varp->nmychunkrec = 0;
                for(i = 0; i < varp->nchunkrec; i ++){
                    if (varp->chunk_owner[i] == nczipp->rank){
                        varp->nmychunkrec++;
                    }
                }
                varp->mychunks = (int*)NCI_Malloc(sizeof(int) * varp->nmychunkrec * varp->nrecalloc);

                if (nczipp->cache_limit_hint == -1){
                    nczipp->cache_limit += varp->nmychunkrec * varp->chunksize;
                }
            }

            varp->nmychunk = oldnmychunk;
            for(i = oldnchunk; i < varp->nchunk; i ++){
                if (varp->chunk_owner[i] == nczipp->rank){
                    varp->mychunks[varp->nmychunk++] = i;
                    //varp->chunk_cache[i] = (void*)NCI_Malloc(varp->chunksize);  // Allocate buffer for blocks we own
                    //memset(varp->chunk_cache[i], 0 , varp->chunksize);
                }
            }

            // Update global chunk count
            nczipp->nmychunks += varp->nmychunk - oldnmychunk;
        }
    }
    else{
        // Notify ncmpio driver
    }

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

    // Resize each var
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (flag_all[i >> 5] & (1u << (i % 32))) {
            flag_all[i >> 5] ^= (1u << (i % 32));
            if ((nczipp->vars.data + i)->dimsize[0] < nczipp->recsize){
                nczipioi_var_resize(nczipp, nczipp->vars.data + i);
            }
        }
    }

    NCI_Free(flag);

    return NC_NOERR;
}