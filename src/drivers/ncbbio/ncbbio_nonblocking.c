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
#include <ncbbio_driver.h>

#define PUT_ARRAY_SIZE 128 /* Size of initial put list */
#define SIZE_MULTIPLIER 2    /* When metadata buffer is full, we'll NCI_Reallocate it to META_BUFFER_MULTIPLIER times the original size*/

/* putlist is a module in ncbbio driver that manage nonblocking put request object
 * It consist of a pool of request object (reqs) and request ids (ids)
 * It's implemented by 2 array of the same number of entries
 * The id i corresponds to the i-th request object
 * We issue request object by issuing the corresponding request id
 * ids is initialized with increasing id, ie. ids[i] = i
 * ids are issued from the begining of ids array
 * A pointer nused keep track of issued ids, it also marks the position of next unused ready to be issued
 * When issuing an id, we take from ids from the position marked by nused and increase nused by 1
 * When recycle an id, we simple put it to the position right before position marked by nused and decrease nused by 1
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
 * ids = 0 0 2 3
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
 *
 */
int ncbbio_put_list_init(NC_bb *ncbbp) {
    int i;
    NC_bb_put_list *lp = &(ncbbp->putlist);

    /* Initialize parameter and allocate the array  */
    lp->nused = 0;
    lp->nalloc = PUT_ARRAY_SIZE;
    lp->reqs = (NC_bb_put_req*)NCI_Malloc(lp->nalloc * sizeof(NC_bb_put_req));
    lp->ids = (int*)NCI_Malloc(lp->nalloc * SIZEOF_INT);
    if (lp->reqs == NULL || lp->ids == NULL) {
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }

    /* Initialize values of ids and reqs
     * Assign increasing unique id
     */
    for (i=0; i<lp->nalloc; i++) {
        lp->ids[i] = i; // Unique ids
        lp->reqs[i].valid = 0;  // Not in use
    }

    return NC_NOERR;
}

/*
 * Enlarge the put list
 * When there are no more unused ids to issue, we must add more ids to the pool
 * We simply enlarge ids and reqs array
 * We initialize the extended part as usual
 */
static int ncbbio_put_list_resize(NC_bb *ncbbp)
{
    int i;
    ssize_t nsize;
    void *ptr;
    NC_bb_put_list *lp = &(ncbbp->putlist);

    /* Calculate new size */
    nsize = lp->nalloc * SIZE_MULTIPLIER;

    /* Realloc reqs and ids */
    ptr = NCI_Realloc(lp->reqs, nsize * sizeof(NC_bb_put_req));
    if (ptr == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM);

    lp->reqs = (NC_bb_put_req*)ptr;
    ptr = NCI_Realloc(lp->ids, nsize * SIZEOF_INT);
    if (ptr == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM);

    lp->ids = (int*)ptr;

    /* Initialize values of ids and reqs
     * Assign increasing unique id
     */
    for (i=lp->nalloc; i<nsize; i++) {
        lp->ids[i] = i; // Unique ids
        lp->reqs[i].valid = 0; // Not in use
    }

    lp->nalloc = nsize;

    return NC_NOERR;
}

/*
 * Clean up the put list
 */
int ncbbio_put_list_free(NC_bb *ncbbp)
{
    NCI_Free(ncbbp->putlist.reqs);
    NCI_Free(ncbbp->putlist.ids);

    return NC_NOERR;
}

/*
 * Allocate a new request object from the putlist with id
 * We first checkif there are unused ids
 * We increase the size of pool, bringing in new ids if there aren't
 * Then we issue the ids at position nused and increase it by 1
 */
int ncbbio_put_list_add(NC_bb *ncbbp, int *id)
{
    int err;
    NC_bb_put_list *lp = &(ncbbp->putlist);

    /* Increase size if necessary */
    if (lp->nused == lp->nalloc) {
        err = ncbbio_put_list_resize(ncbbp);
        if (err != NC_NOERR) return err;
    }

    /* Get the first unused id marked by nused */
    *id = lp->ids[lp->nused++];
    // Initialize new request object
    lp->reqs[*id].valid = 1;    // Id in use means request object in use
    lp->reqs[*id].ready = 0;    // Status not avaiable yet
    lp->reqs[*id].status = NC_NOERR;

    return NC_NOERR;
}

/*
 * Recycle a request object in the put list
 * We simply put the id back to ids array
 * We put it at the empty slot right before position marked by nused
 * Decrease nused by 1 to mark the recycled id as unused
 */
int ncbbio_put_list_remove(NC_bb *ncbbp, int reqid) {
    NC_bb_put_list *lp = &(ncbbp->putlist);

    /* Mark entry as invalid
     * When the id is recycled, we also take back the request object
     */
    lp->reqs[reqid].valid = 0;

    /* Return id to the list */
    lp->ids[--lp->nused] = reqid;

    return NC_NOERR;
}

/*
 * Process put request
 * We process the request in the request object and recycle the request id
 * We first check if the status is ready
 * If it is, the corresponding log entry must have been flushed and the result is already avaiable in status, we return it directly
 * If not, the operation is still pending in the log file.
 * We do a log flush. Since log entries are linked to the request object, the status should be updated when the corresponding entry is flushed
 * Log module is responsible to fill up the status of corresponding request object
 * The request should be ready after log flush
 */
int ncbbio_handle_put_req(NC_bb *ncbbp, int reqid, int *stat)
{
    int err, status = NC_NOERR;
    NC_bb_put_list *lp = &(ncbbp->putlist);
    NC_bb_put_req *req;

    /* Filter invalid reqid
     * Valid id range from 0 ~ nalloc - 1
     */
    if (reqid >= lp->nalloc || reqid < 0) {
        if (stat != NULL) *stat = NC_EINVAL_REQUEST;
        return status;
    }

    // Locate the req object, which is reqs[reqid]
    req = lp->reqs + reqid;

    /* Filter invalid reqid
     * The request object must be in used
     * If not, the id is not yet issued, and hence invalid
     */
    if (!req->valid) {
        if (stat != NULL) *stat = NC_EINVAL_REQUEST;
        return status;
    }

    /* Flush is done whenever a wait is called
     * Log module is responsible to update the request obejct when entries are flushed
     * This should never happen
     * If it do, we have mising log entry or corrupt metadata index
     */
    if (!req->ready) {
        printf("Fatal error: nonblocking request not in log file\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // Return status to the user
    if (stat != NULL) *stat = req->status;

    // Recycle req object to the pool
    err = ncbbio_put_list_remove(ncbbp, reqid);
    if (status == NC_NOERR) status = err;

    return status;
}

/*
 * Process all put request
 * We didbn't keep track of issued ids
 * To process all issued ids, we need to do a linear search on reqs and process all request object that is in use
 */
int ncbbio_handle_all_put_req(NC_bb *ncbbp) {
    int i, err, status = NC_NOERR;
    NC_bb_put_list *lp = &(ncbbp->putlist);

    // Search through req object array for object in use */
    for (i=0; i<lp->nalloc; i++) {
        if (lp->reqs[i].valid) {
            err = ncbbio_handle_put_req(ncbbp, i, NULL);
            if (status == NC_NOERR) status = err;
        }
    }

    return NC_NOERR;
}

/*
 * Process put request
 * If corresponding log entries are already flushed, we retrun error
 * If not, we mark those log entries as invalid, preventing them from being flushed
 * After processing the request, we simply recycle the id
 */
int ncbbio_cancel_put_req(NC_bb *ncbbp, int reqid, int *stat) {
    int i, err, status = NC_NOERR;
    NC_bb_put_list *lp = &(ncbbp->putlist);
    NC_bb_put_req *req;

    /* Filter invalid reqid
     * Valid id range from 0 ~ nalloc - 1
     */
    if (reqid >= lp->nalloc || reqid < 0) {
        if (stat != NULL) *stat = NC_EINVAL_REQUEST;
        return status;
    }

    // Locate the req object, which is reqs[reqid]
    req = lp->reqs + reqid;

    /* Filter invalid reqid
     * The request object must be in used
     * If not, the id is not yet issued, and hence invalid
     */
    if (!req->valid) {
        if (stat != NULL) *stat = NC_EINVAL_REQUEST;
        return status;
    }

    /* If log entry is already flushed, it's too late to cancel
     */
    if (req->ready) {
        if (stat != NULL) *stat = NC_EFLUSHED;    // Fail
    }
    else {
        if (stat != NULL) *stat = NC_NOERR;   // Success

        // Mark log entries as invalid
        for (i=req->entrystart; i<req->entryend; i++)
            ncbbp->metaidx.entries[i].valid = 0;
    }

    // Recycle req object to the pool
    err = ncbbio_put_list_remove(ncbbp, reqid);
    if (status == NC_NOERR) status = err;

    return status;
}

/*
 * Cancel all put request
 * We didbn't keep track of issued ids
 * To process all issued ids, we need to do a linear search on reqs and process all request object that is in use
 */
int ncbbio_cancel_all_put_req(NC_bb *ncbbp) {
    int i, err, status = NC_NOERR;
    NC_bb_put_list *lp = &(ncbbp->putlist);

    // Search through req object list for valid objects */
    for (i=0; i<lp->nalloc; i++) {
        if (lp->reqs[i].valid) {
            err = ncbbio_cancel_put_req(ncbbp, i, NULL);
            if (status == NC_NOERR) status = err;
        }
    }

    return NC_NOERR;
}


