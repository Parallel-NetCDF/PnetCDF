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
    int i;
    unsigned int j;
    int nvar, nflag;
    unsigned int *flag, *flag_all;
    int *vids;
    NC_zip_req *req;

    // Flag of touched vars
    nflag = nczipp->vars.cnt / 32 + 1;
    flag = (unsigned int*)NCI_Malloc(sizeof(int) * nflag * 2);
    flag_all = flag + nflag;
    memset(flag, 0, sizeof(int) * nflag);
    for(i = 0; i < nreq; i++){
        req = nczipp->putlist.reqs + reqids[i];
        flag[req->varid >> 5] |= 1u << (req->varid % 32);
    }
    
    // Sync flag
    CHK_ERR_ALLREDUCE(flag, flag_all, nflag, MPI_UNSIGNED, MPI_BOR, nczipp->comm);

    // Build a skip list of touched vars
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (flag_all[i >> 5] & (1u << (i % 32))) {
            nvar++;
        }
    }
    vids = (int*)NCI_Malloc(sizeof(int) * nvar);
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (flag_all[i >> 5] & (1u << (i % 32))) {
            vids[nvar++] = i;
        }
    }  

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB)
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB_WAIT)

    // Perform collective buffer
    if (nczipp->comm_unit == NC_ZIP_COMM_CHUNK){
        nczipioi_iput_cb_chunk(nczipp, nreq, reqids, stats);
    }
    else{
        nczipioi_iput_cb_proc(nczipp, nreq, reqids, stats);
    }

    // Perform I/O for comrpessed variables
    nczipioi_save_nvar(nczipp, nvar, vids);

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB_WAIT)

    // Free buffers
    NCI_Free(vids);
    NCI_Free(flag);

    return NC_NOERR;
}

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int nczipioi_wait_get_reqs(NC_zip *nczipp, int nreq, int *reqids, int *stats){
    int err;
    int i;
    unsigned int j;
    int nvar, nflag;
    unsigned int *flag, *flag_all;
    int *vids;
    NC_zip_req *req;

    // Flag of touched vars
    nflag = nczipp->vars.cnt / 32 + 1;
    flag = (unsigned int*)NCI_Malloc(sizeof(int) * nflag * 2);
    flag_all = flag + nflag;
    memset(flag, 0, sizeof(int) * nflag);
    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        flag[req->varid >> 5] |= 1u << (req->varid % 32);
    }
    
    // Sync flag
    CHK_ERR_ALLREDUCE(flag, flag_all, nflag, MPI_UNSIGNED, MPI_BOR, nczipp->comm);

    // Build a skip list of touched vars
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (flag_all[i >> 5] & (1u << (i % 32))) {
            nvar++;
        }
    }
    vids = (int*)NCI_Malloc(sizeof(int) * nvar);
    nvar = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        if (flag_all[i >> 5] & (1u << (i % 32))) {
            vids[nvar++] = i;
        }
    }  

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB)
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB_WAIT)

    // Perform I/O for comrpessed variables
    //nczipioi_load_nvar(nczipp, nvar, vids);

    // Perform collective buffer
    if (nczipp->comm_unit == NC_ZIP_COMM_CHUNK){
        nczipioi_iget_cb_chunk(nczipp, nreq, reqids, stats);
    }
    else{
        nczipioi_iget_cb_proc(nczipp, nreq, reqids, stats);
        //nczipioi_iget_cb_chunk(nczipp, nreq, reqids, stats);
    }

    for(i = 0; i < nreq; i++){
        req = nczipp->getlist.reqs + reqids[i];
        if (req->buf != req->xbuf){
            void *cbuf=(void*)req->buf;
            
            err = nczipioiconvert(req->xbuf, cbuf, nczipp->vars.data[req->varid].etype, req->buftype, req->bufcount);
            if (err != NC_NOERR) return err;

            if (cbuf != req->buf) NCI_Free(cbuf);
        }
    }

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB_WAIT)
    
    // Free buffers
    NCI_Free(vids);
    NCI_Free(flag);

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

    if (nczipp->delay_init){
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB)
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB_WAIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT_META)

        err = nczipioi_init_nvar(nczipp, nput, putreqs, nget, getreqs);   // nput + nget = real nreq
        if (err != NC_NOERR){
            return err;
        }
        
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT)
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT_META)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB_WAIT)
    }
    else{
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB)
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB_WAIT)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_RESIZE)

        // Sync number of rec
        err = nczipioi_resize_nvar(nczipp, nput, putreqs, nget, getreqs);   // nput + nget = real nreq
        if (err != NC_NOERR){
            return err;
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_RESIZE)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB)
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB_WAIT)
    }

    if (stats != NULL){
        putstats = (int*)NCI_Malloc(sizeof(int) * nput);
        getstats = (int*)NCI_Malloc(sizeof(int) * nget);
        memset(putstats, 0, sizeof(int) * nput);
        memset(getstats, 0, sizeof(int) * nget);
    }
    else{
        putstats = NULL;
        getstats = NULL;
    }

    if ((nczipp->mode & NC_WRITE) && nreqs != NC_GET_REQ_ALL){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_PUT)
        nczipioi_wait_put_reqs(nczipp, nput, putreqs, putstats);
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_PUT)
    }
    
    if (nreqs != NC_PUT_REQ_ALL){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_GET)
        nczipioi_wait_get_reqs(nczipp, nget, getreqs, getstats);
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_GET)
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
