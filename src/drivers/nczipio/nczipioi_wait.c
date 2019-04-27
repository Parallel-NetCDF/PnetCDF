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


/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int nczipioi_wait_put_reqs(NC_zip *nczipp, int nreq, int *reqids, int *stats){
    int err;
    int i, j, k, l;
    int nvar, needinit;
    int maxnreq = 0;
    int *varids;
    int *nreqs, *offs;  // Number of reqids in each variable
    int *smap;
    int *init;
    MPI_Offset **starts, **counts;
    NC_zip_req *req;
    NC_zip_var *varp;

    // Build a skip list of touched vars
    nreqs = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt);
    smap = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt);
    memset(nreqs, 0, sizeof(int) * nczipp->vars.cnt);
    for(i = 0; i < nreq; i++){
        req = nczipp->putlist.reqs + reqids[i];
        nreqs[req->varid] += req->nreq;
    }
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (nreqs[i]){
            nvar++;
        }
    }
    varids = (int*)NCI_Malloc(sizeof(int) * nvar);
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (nreqs[i]){
            varids[nvar] = i;
            smap[i] = nvar++;
        }
    }

    if (nczipp->delay_init){
        init = (int*)NCI_Malloc(sizeof(int) * nvar * 2);
        offs = init + nvar;
        
        offs[0] = 0;
        for(i = 0; i < nvar; i++){
            varp = nczipp->vars.data + varids[i];
            if (varp->chunkdim == NULL){
                offs[i + 1] = offs[i] + nreqs[varids[i]];
                init[i] = 1;
            }
            else{
                offs[i + 1] = offs[i];
            }
        }

        memset(nreqs, 0, sizeof(int) * nvar); //reuse it for counter

        starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * offs[nvar] * 2);
        counts = starts + offs[nvar];

        for(i = 0; i < nreq; i++){
            req = nczipp->putlist.reqs + reqids[i];
            k = smap[req->varid];
            if (init[k]){
                if (req->nreq > 1){
                    for(j = 0; j < req->nreq; j++){
                        starts[offs[k] + nreqs[k]] = req->starts[j];
                        counts[offs[k] + (nreqs[k]++)] = req->counts[j];
                    }
                }
                else{
                    starts[offs[k] + nreqs[k]] = req->start;
                    counts[offs[k] + (nreqs[k]++)] = req->count;     
                }
            }
        }

        for(i = 0; i < nvar; i++){
            varp = nczipp->vars.data + varids[i];
            if (init[i]){
                nczipioi_var_init(nczipp, varp, 1, nreqs[i], starts + offs[i], counts + offs[i]);
            }
            else if (varp->isrec && varp->dimsize[0] < nczipp->recsize){
                nczipioi_var_resize(nczipp, varp);
            }
        }

        NCI_Free(init);
        NCI_Free(starts);
    }

    // Perform collective buffer
    if (nczipp->comm_unit == NC_ZIP_COMM_CHUNK){
        nczipioi_iput_cb_chunk(nczipp, nreq, reqids, stats);
    }
    else{
        nczipioi_iput_cb_proc(nczipp, nreq, reqids, stats);
    }

    // Perform I/O for comrpessed variables
    nczipioi_save_nvar(nczipp, nvar, varids);

    // Free buffers
    NCI_Free(varids);
    NCI_Free(nreqs);
    NCI_Free(smap);

    return NC_NOERR;
}

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int nczipioi_wait_get_reqs(NC_zip *nczipp, int nreq, int *reqids, int *stats){
    int err;
    int i, j, k;
    int nvar;
    int *nreqs, *offs;  // Number of reqids in each variable
    int *smap;
    int *init;
    int *varids;
    MPI_Offset **starts, **counts;
    NC_zip_req *req;
    NC_zip_var *varp;

    // Build a skip list of touched vars
    nreqs = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt);
    memset(nreqs, 0, sizeof(int) * nczipp->vars.cnt);
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        nreqs[req->varid] = 1;
    }
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (nreqs[i]){
            nvar++;
        }
    }
    varids = (int*)NCI_Malloc(sizeof(int) * nvar);
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (nreqs[i]){
            varids[nvar++] = i;
        }
    }

    if (nczipp->delay_init){
        init = (int*)NCI_Malloc(sizeof(int) * nvar * 2);
        offs = init + nvar;
        
        offs[0] = 0;
        for(i = 0; i < nvar; i++){
            varp = nczipp->vars.data + varids[i];
            if (varp->chunkdim == NULL){
                offs[i + 1] = offs[i] + nreqs[varids[i]];
                init[i] = 1;
            }
            else{
                offs[i + 1] = offs[i];
            }
        }

        memset(nreqs, 0, sizeof(int) * nvar); //reuse it for counter

        starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * offs[nvar] * 2);
        counts = starts + offs[nvar];

        for(i = 0; i < nreq; i++){
            req = nczipp->putlist.reqs + reqids[i];
            k = smap[req->varid];
            if (init[k]){
                if (req->nreq > 1){
                    for(j = 0; j < req->nreq; j++){
                        starts[offs[k] + nreqs[k]] = req->starts[j];
                        counts[offs[k] + (nreqs[k]++)] = req->counts[j];
                    }
                }
                else{
                    starts[offs[k] + nreqs[k]] = req->start;
                    counts[offs[k] + (nreqs[k]++)] = req->count;     
                }
            }
        }

        for(i = 0; i < nvar; i++){
            varp = nczipp->vars.data + varids[i];
            if (init[i]){
                nczipioi_var_init(nczipp, varp, 0, nreqs[i], starts + offs[i], counts + offs[i]);
            }
            else if (varp->isrec && varp->dimsize[0] < nczipp->recsize){
                nczipioi_var_resize(nczipp, varp);
            }
        }

        NCI_Free(starts);
        NCI_Free(init);
    }

    // Perform I/O for comrpessed variables
    nczipioi_load_nvar(nczipp, nvar, varids);

    // Perform collective buffer
    if (nczipp->comm_unit == NC_ZIP_COMM_CHUNK){
        nczipioi_iget_cb_chunk(nczipp, nreq, reqids, stats);
    }
    else{
        nczipioi_iget_cb_proc(nczipp, nreq, reqids, stats);
        //nczipioi_iget_cb_chunk(nczipp, nreq, reqids, stats);
    }

    // Free buffers
    NCI_Free(varids);
    NCI_Free(nreqs);
}

int
nczipioi_wait(NC_zip *nczipp, int nreqs, int *reqids, int *stats, int reqMode){
    int err;
    int i;
    int nput = 0, nget = 0;
    int *putreqs = NULL, *getreqs = NULL;
    int *putstats = NULL, *getstats = NULL;

    if (nreqs == NC_REQ_ALL || nreqs == NC_PUT_REQ_ALL){
        nput = nczipp->putlist.nused;
        putreqs = (int*)NCI_Malloc(sizeof(int) * nput);
        memcpy(putreqs, nczipp->putlist.ids, nput * sizeof(int));
    }
    if(nreqs == NC_REQ_ALL || nreqs == NC_GET_REQ_ALL){
        nget = nczipp->getlist.nused;
        getreqs = (int*)NCI_Malloc(sizeof(int) * nget);
        memcpy(getreqs, nczipp->getlist.ids, nget * sizeof(int));
    }

    if (nreqs > 0){
        // Count number of get and put requests
        for(i = 0; i < nreqs; i++){
            if (reqids[i] & 1){
                nput++;
            }
        }

        // Allocate buffer
        nget = nreqs - nput;
        putreqs = (int*)NCI_Malloc(sizeof(int) * nput);
        getreqs = (int*)NCI_Malloc(sizeof(int) * nget);
        
        // Build put and get req list
        nput = nget = 0;
        for(i = 0; i < nreqs; i++){
            if (reqids[i] & 1){
                putreqs[nput++] = reqids[i] >> 1;
            }
            else{
                getreqs[nget++] = reqids[i] >> 1;
            }
        }
    }

    if (stats != NULL){
        putstats = (int*)NCI_Malloc(sizeof(int) * nput);
        getstats = (int*)NCI_Malloc(sizeof(int) * nget);
    }
    else{
        putstats = NULL;
        getstats = NULL;
    }

    if (nput > 0){
        nczipioi_wait_put_reqs(nczipp, nput, putreqs, putstats);
    }
    
    if (nget > 0){
        nczipioi_wait_get_reqs(nczipp, nget, getreqs, getstats);
    }

    // Assign stats
    if (stats != NULL){
        nput = nget = 0;
        for(i = 0; i < nreqs; i++){
            if (reqids[i] & 1){
                stats[i] = putstats[nput++];
            }
            else{
                stats[i] = getstats[nget++];
            }
        }

        NCI_Free(putstats);
        NCI_Free(getstats);
    }

    // Remove from req list
    for(i = 0; i < nput; i++){
        nczipioi_req_list_remove(&(nczipp->putlist), putreqs[i]);
    }
    for(i = 0; i < nget; i++){
        nczipioi_req_list_remove(&(nczipp->getlist), getreqs[i]);
    }

    free(putreqs);
    free(getreqs);

    return NC_NOERR;
}
