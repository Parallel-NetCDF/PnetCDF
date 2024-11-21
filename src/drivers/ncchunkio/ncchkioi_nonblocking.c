/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <pnc_debug.h>
#include <common.h>
#include <stdio.h>
#include <string.h>
#include <ncchkio_driver.h>
#include <ncchkio_internal.h>

#define PUT_ARRAY_SIZE 128 /* Size of initial put list */
#define SIZE_MULTIPLIER 2    /* When metadata buffer is full, we'll NCI_Reallocate it to META_BUFFER_MULTIPLIER times the original size*/

/* getlist is a module in ADIOS driver that manage nonblocking get request object
 * It consist of a pool of request object (reqs) and request ids (ids)
 * It's implemented by 3 array of the same number of entries
 * The id i corresponds to the i-th request object
 * We issue request object by issuing the corresponding request id
 * ids is initialized with increasing id, ie. ids[i] = i
 * ids are issued from the begining of ids array
 * We keep track of location of id in ids array in pos array. Initially, pos[i] = i
 * A pointer nused keep track of issued ids, it also marks the position of next unused ready to be issued
 * ids[0:nused] => active (used) request ids
 * ids[nused:nalloc] => available (unused) request ids
 * When issuing an id, we take from ids from the position marked by nused and increase nused by 1
 * When recycling an id, we swap if with the right before position marked by nused and decrease nused by 1 so that it falls to unused pool
 * NOTE: We does not guarantee to issue id in continuous and increasing order
 * NOTE: ids is simply a pool housing reqeust ids, the position od id within ids is not fixed and has no meaning
 *
 * Eaxmple:
 * Initial:
 * ids = 0 1 2 3
 *       ^
 *   nused = 0
 * After issuing 2 ids:
 * undefined|Avaiable ids --->
 * ids = 0 1 2 3
 *           ^
 *       nused = 2
 * Recycling id 0
 *        |Avaiable ids --->
 * ids = 1 0 2 3
 *         ^
 *     nused = 1
* Recycling id 1
 *      |Avaiable ids --->
 * ids = 1 0 2 3
 *       ^
 *   nused = 0
 */

/*
 * Initialize the put list
 * ids[0:nused] => active (used) request ids
 * ids[nused:nalloc] => available (unused) request ids
 */
int ncchkioi_req_list_init(NC_chk_req_list *lp) {
    int err=NC_NOERR;
    int i;

    /* Initialize parameter and allocate the array  */
    lp->nused = 0;
    lp->nalloc = PUT_ARRAY_SIZE;
    lp->reqs = (NC_chk_req*)NCI_Malloc(lp->nalloc * sizeof(NC_chk_req));
    lp->ids = (int*)NCI_Malloc(lp->nalloc * SIZEOF_INT);
    CHK_PTR(lp->ids)
    lp->pos = (int*)NCI_Malloc(lp->nalloc * SIZEOF_INT);
    CHK_PTR(lp->pos)
    if (lp->reqs == NULL || lp->ids == NULL) {
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }

    /* Initialize values of ids and reqs
     * Assign increasing unique id
     */
    for (i=0; i<lp->nalloc; i++) {
        lp->ids[i] = i; // Unique ids
        lp->pos[i] = i; // Not in use
    }

err_out:;
    return err;
}

/*
 * Enlarge the put list
 * When there are no more unused ids to issue, we must add more ids to the pool
 * We simply enlarge ids and reqs array
 * We initialize the extended part as usual
 */
static int ncchkioi_req_list_resize(NC_chk_req_list *lp)
{
    int i;
    size_t nsize;
    void *ptr;

    /* Calculate new size */
    nsize = lp->nalloc * SIZE_MULTIPLIER;

    /* Realloc reqs and ids */
    ptr = NCI_Realloc(lp->reqs, nsize * sizeof(NC_chk_req));
    if (ptr == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM);
    lp->reqs = (NC_chk_req*)ptr;

    ptr = NCI_Realloc(lp->ids, nsize * SIZEOF_INT);
    if (ptr == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM);
    lp->ids = (int*)ptr;

    ptr = NCI_Realloc(lp->pos, nsize * SIZEOF_INT);
    if (ptr == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM);
    lp->pos = (int*)ptr;

    /* Initialize values of ids and reqs
     * Assign increasing unique id
     */
    for (i=lp->nalloc; i<nsize; i++) {
        lp->ids[i] = i; // Unique ids
        lp->pos[i] = i; // Default position
    }

    lp->nalloc = nsize;

    return NC_NOERR;
}

/*
 * Clean up the put list
 */
int ncchkioi_req_list_free(NC_chk_req_list *lp)
{
    NCI_Free(lp->reqs);
    NCI_Free(lp->ids);
    NCI_Free(lp->pos);

    return NC_NOERR;
}

/*
 * Allocate a new request object from the getlist with id
 * We first check if there are unused ids
 * We increase the size of pool, bringing in new ids if there aren't
 * Then we issue the ids at position nused and increase it by 1
 */
int ncchkioi_req_list_add(NC_chk_req_list *lp, int *id)
{
    int err;

    /* Increase size if necessary */
    if (lp->nused == lp->nalloc) {
        err = ncchkioi_req_list_resize(lp);
        if (err != NC_NOERR) return err;
    }

    /* Get the first unused id marked by nused */
    *id = lp->ids[lp->nused++];

    return NC_NOERR;
}

/*
 * Recycle a request object in the put list
 * We need to maintain the position of each request id in the ids list
 * ids[0:nused] => active (used) request ids
 * ids[nused:nalloc] => available (unused) request ids
 */
int ncchkioi_req_list_remove(NC_chk_req_list *lp, int reqid) {
    NC_chk_req * req = lp->reqs + reqid;

    /* Clean up request */
    if (req->start != NULL){
        NCI_Free(req->start);
    }
    if (req->count != NULL){
        NCI_Free(req->count);
    }
    if (req->starts != NULL){
        NCI_Free(req->starts);
    }
    if (req->counts != NULL){
        NCI_Free(req->counts);
    }
    if (req->stride != NULL){
        NCI_Free(req->stride);
    }
    if (req->xbufs != NULL){
        NCI_Free(req->xbufs);
    }
    if (req->xbuf != req->buf){
        NCI_Free(req->xbuf);
    }

    /* Return id to the list */
    lp->nused--;
    lp->ids[lp->pos[reqid]] = lp->ids[lp->nused];
    lp->pos[lp->ids[lp->nused]] = lp->pos[reqid];
    lp->ids[lp->nused] = reqid;
    lp->pos[reqid] = lp->nused;

    return NC_NOERR;
}
