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
    int i, j;
    int nvar;
    int vid;   // Iterators for variable id
    int *varids;
    int *nreqs;  // Number of reqids in each variable
    int *nums;  // Number of reqs in each varn
    int **vreqids;
    int num, maxnum = 0;
    MPI_Offset **starts, **counts, **strides;
    MPI_Offset rsize;
    char **bufs;
    NC_zip_req *req;

    // Count total number of request in per variable for packed varn request
    nums = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt);
    nreqs = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt);
    memset(nums, 0, sizeof(int) * nczipp->vars.cnt);
    memset(nreqs, 0, sizeof(int) * nczipp->vars.cnt);
    for(i = 0; i < nreq; i++){
        req = nczipp->putlist.reqs + reqids[i];
        nreqs[req->varid]++;
        nums[req->varid] += req->nreq;
    }

    /* Allocate a skip list of reqids for each vriable
     * At the same time, we find out the number of starts and counts we need to allocate
     */
    vreqids = (int**)NCI_Malloc(sizeof(int*) * nczipp->vars.cnt);
    vreqids[0] = (int*)NCI_Malloc(sizeof(int) * nreq);
    maxnum = 0;
    i = 0;
    nvar = 0;
    for(vid = 0; vid < nczipp->vars.cnt; vid++){
        if (nreqs[vid] > 0){
            // Assign buffer to reqid skip list
            vreqids[vid] = vreqids[0] + i;
            i += nreqs[vid];

            // maximum number of starts and counts we need across all variables
            if (maxnum < nums[vid]){
                maxnum = nums[vid];
            }

            // Number of variable that has request to write
            nvar++;
        }
    }

    varids = (int*)NCI_Malloc(sizeof(int) * nvar);

    // Fill up the skip list
    memset(nreqs, 0, sizeof(int) * nczipp->vars.cnt);
    for(i = 0; i < nreq; i++){
        req = nczipp->putlist.reqs + reqids[i];
        vreqids[req->varid][nreqs[req->varid]++] = reqids[i];
    }
    
    // Allocate parameters
    starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * maxnum);
    counts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * maxnum);
    bufs =  (char**)NCI_Malloc(sizeof(char*) * maxnum);

    /* Pack requests variable by variable
     */
    nvar = 0;
    for(vid = 0; vid < nczipp->vars.cnt; vid++){
        if (nreqs[vid] > 0){
            // Fill varid in the skip list
            varids[nvar++] = vid;

            // Collect parameters
            num = 0;
            for(j = 0; j < nreqs[vid]; j++){
                req = nczipp->putlist.reqs + vreqids[vid][j];

                if (req->nreq > 1){
                    for(i = 0; i < req->nreq; i++){
                        starts[num] = req->starts[i];
                        counts[num] = req->counts[i];
                        bufs[num++] = req->xbufs[i];
                    }
                }
                else{
                    starts[num] = req->start;
                    counts[num] = req->count;
                    bufs[num++] = req->xbuf;
                }
            }

            // Perform collective buffering
            nczipioi_put_varn_cb(nczipp, nczipp->vars.data + vid, num, starts, counts, NULL, bufs);
        }
    }

    // Perform I/O for comrpessed variables
    nczipioi_save_nvar(nczipp, nvar, varids);

    // Free buffers
    NCI_Free(nums);
    NCI_Free(nreqs);

    NCI_Free(vreqids[0]);
    NCI_Free(vreqids);

    NCI_Free(varids);
    
    NCI_Free(starts);
    NCI_Free(counts);
    NCI_Free(bufs);

    return NC_NOERR;
}

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int nczipioi_wait_get_reqs(NC_zip *nczipp, int nreq, int *reqids, int *stats){
    int err;
    int i, j;
    int nvar;
    int vid;   // Iterators for variable id
    int *varids;
    int *nreqs;  // Number of reqids in each variable
    int *nums;  // Number of reqs in each varn
    int **vreqids;
    int num, maxnum = 0;
    MPI_Offset **starts, **counts, **strides;
    MPI_Offset rsize;
    char **bufs;
    NC_zip_req *req;

    // Count total number of request in per variable for packed varn request
    nums = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt);
    nreqs = (int*)NCI_Malloc(sizeof(int) * nczipp->vars.cnt);
    memset(nums, 0, sizeof(int) * nczipp->vars.cnt);
    memset(nreqs, 0, sizeof(int) * nczipp->vars.cnt);
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        nreqs[req->varid]++;
        nums[req->varid] += req->nreq;
    }

    /* Allocate a skip list of reqids for each vriable
     * At the same time, we find out the number of starts and counts we need to allocate
     */
    vreqids = (int**)NCI_Malloc(sizeof(int*) * nczipp->vars.cnt);
    vreqids[0] = (int*)NCI_Malloc(sizeof(int) * nreq);
    maxnum = 0;
    i = 0;
    nvar = 0;
    for(vid = 0; vid < nczipp->vars.cnt; vid++){
        if (nreqs[vid] > 0){
            // Assign buffer to reqid skip list
            vreqids[vid] = vreqids[0] + i;
            i += nreqs[vid];

            // maximum number of starts and counts we need across all variables
            if (maxnum < nums[vid]){
                maxnum = nums[vid];
            }

            // Number of variable that has request to write
            nvar++;
        }
    }

    varids = (int*)NCI_Malloc(sizeof(int) * nvar);

    // Fill up the skip list
    memset(nreqs, 0, sizeof(int) * nczipp->vars.cnt);
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        vreqids[req->varid][nreqs[req->varid]++] = reqids[i];
    }
    
    // Allocate parameters
    starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * maxnum);
    counts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * maxnum);
    bufs =  (char**)NCI_Malloc(sizeof(char*) * maxnum);

    /* Pack requests variable by variable
     */
    nvar = 0;
    for(vid = 0; vid < nczipp->vars.cnt; vid++){
        if (nreqs[vid] > 0){
            // Fill varid in the skip list
            varids[nvar++] = vid;

            // Collect parameters
            num = 0;
            for(j = 0; j < nreqs[vid]; j++){
                req = nczipp->getlist.reqs + vreqids[vid][j];

                if (req->nreq > 1){
                    for(i = 0; i < req->nreq; i++){
                        starts[num] = req->starts[i];
                        counts[num] = req->counts[i];
                        bufs[num++] = req->xbufs[i];
                    }
                }
                else{
                    starts[num] = req->start;
                    counts[num] = req->count;
                    bufs[num++] = req->xbuf;
                }
            }

            // Perform collective buffering
            nczipioi_get_varn_cb(nczipp, nczipp->vars.data + vid, num, starts, counts, NULL, bufs);
        }
    }

    // Free buffers
    NCI_Free(nums);
    NCI_Free(nreqs);

    NCI_Free(vreqids[0]);
    NCI_Free(vreqids);

    NCI_Free(varids);
    
    NCI_Free(starts);
    NCI_Free(counts);
    NCI_Free(bufs);

    return NC_NOERR;
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

    NCI_Free(putreqs);
    NCI_Free(getreqs);

    return NC_NOERR;
}
