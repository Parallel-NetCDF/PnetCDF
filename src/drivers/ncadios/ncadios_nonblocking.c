/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <pnc_debug.h>
#include <common.h>
#include <stdio.h>
#include <string.h>
#include <ncadios_driver.h>
#include <ncadios_internal.h>

#define PUT_ARRAY_SIZE 128 /* Size of initial put list */
/* When list is full, we will reallocate it to
 * SIZE_MULTIPLIER times the original size
 */
#define SIZE_MULTIPLIER 2

/* getlist is a module in ADIOS driver that manage nonblocking get request
 * object.
 * It consist of a pool of request object (reqs) and request ids (ids)
 * It's implemented by 3 array of the same number of entries
 * The id i corresponds to the i-th request object
 * We issue request object by issuing the corresponding request id
 * ids is initialized with increasing id, ie. ids[i] = i
 * ids are issued from the begining of ids array
 * We keep track of location of id in ids array in pos array.
 * Initially, pos[i] = i
 * A pointer nused keep track of issued ids, it also marks the position of
 * next unused ready to be issued.
 * ids[0:nused] => active (used) request ids
 * ids[nused:nalloc] => available (unused) request ids
 * When issuing an id, we take from ids from the position marked by nused and
 * increase nused by 1.
 * When recycling an id, we swap if with the right before position marked by
 * nused and decrease nused by 1 so that it falls to unused pool.
 * NOTE: We does not guarantee to issue id in continuous and increasing order
 * NOTE: ids is simply a pool housing reqeust ids, the position od id within
 * ids is not fixed and has no meaning.
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
int ncadiosi_get_list_init(NC_ad_get_list *lp) {
    int i;

    /* Initialize parameter and allocate the array  */
    lp->nused = 0;
    lp->nalloc = PUT_ARRAY_SIZE;
    lp->reqs = (NC_ad_get_req*)NCI_Malloc(lp->nalloc * sizeof(NC_ad_get_req));
    lp->ids = (int*)NCI_Malloc(lp->nalloc * SIZEOF_INT);
    lp->pos = (int*)NCI_Malloc(lp->nalloc * SIZEOF_INT);
    if (lp->reqs == NULL || lp->ids == NULL) {
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }

    /* Initialize values of ids and reqs
     * Assign increasing unique id
     */
    for (i=0; i<lp->nalloc; i++) {
        lp->ids[i] = i; /* Unique ids */
        lp->pos[i] = i; /* Position of id i in ids array */
    }

    return NC_NOERR;
}

/*
 * Enlarge the put list
 * When there are no more unused ids to issue, we must add more ids to the pool
 * We simply enlarge ids and reqs array
 * We initialize the extended part as usual
 */
static int ncadiosi_get_list_resize(NC_ad_get_list *lp)
{
    int i;
    size_t nsize;
    void *ptr;

    /* Calculate new size */
    nsize = lp->nalloc * SIZE_MULTIPLIER;

    /* Realloc reqs and ids */
    ptr = NCI_Realloc(lp->reqs, nsize * sizeof(NC_ad_get_req));
    if (ptr == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM);
    lp->reqs = (NC_ad_get_req*)ptr;

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
        lp->ids[i] = i; /* Unique ids */
        lp->pos[i] = i; /* Position of id i in ids array */
    }

    lp->nalloc = nsize;

    return NC_NOERR;
}

/*
 * Clean up the put list
 */
int ncadiosi_get_list_free(NC_ad_get_list *lp)
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
int ncadiosi_get_list_add(NC_ad_get_list *lp, int *id)
{
    int err;

    /* Increase size if necessary */
    if (lp->nused == lp->nalloc) {
        err = ncadiosi_get_list_resize(lp);
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
int ncadiosi_get_list_remove(NC_ad_get_list *lp, int reqid) {
    /* Return id to the list */
    lp->nused--;
    lp->ids[lp->pos[reqid]] = lp->ids[lp->nused];
    lp->pos[lp->ids[lp->nused]] = lp->pos[reqid];
    lp->ids[lp->nused] = reqid;
    lp->pos[reqid] = lp->nused;

    return NC_NOERR;
}

/*
 * All adios perform read and mark requests corresponds to posted operation as completed
 */
int ncadiosi_perform_read(NC_ad *ncadp) {
    int i, err;
    NC_ad_get_list *lp = &(ncadp->getlist);

    /* Read all posted operation */
    err = adios_perform_reads (ncadp->fp, 1);
    if (err != 0){
        err = ncmpii_error_adios2nc(adios_errno, "Open");
        DEBUG_RETURN_ERROR(err);
    }

    /* All current active request has been read by ADIOS */
    for(i = 0; i < lp->nused; i++){
        lp->reqs[lp->ids[i]].ready = 1;
    }

    return NC_NOERR;
}

/*
 * Process a read request and return it's status
 * We need to call ADIOS perform read if data hasn't been read
 * If data is available, we perform necessary type and shape converion
 */
int ncadiosi_handle_get_req(NC_ad *ncadp, NC_ad_get_req *req){
    int err, status = NC_NOERR;
    int cesize;

    /* Perform ADIOS read if this request haven't been read */
    if (!(req->ready)){
        ncadiosi_perform_read(ncadp);
    }

    /* If type do not match */
    if (req->vtype != req->buftype){
        err = ncadiosiconvert(req->xbuf, req->cbuf, req->vtype, req->buftype,
                                (int)(req->ecnt));
        if (status == NC_NOERR){
            status = err;
        }
        NCI_Free(req->xbuf);
    }

    /* If imap is used (memory buffer not contiguous) */
    if (req->cbuf != req->buf){
        int position = 0;

        MPI_Unpack(req->cbuf, req->cbsize, &position, req->buf, 1,
                    req->imaptype, MPI_COMM_SELF);
        MPI_Type_free(&(req->imaptype));

        NCI_Free(req->cbuf);
    }

    /* Free up structured used to post ADIOS operation */
    if (req->points != NULL){
        NCI_Free(req->points);
    }

    adios_selection_delete(req->sel);

    /* Record get size */
    MPI_Type_size(req->vtype, &cesize);
    ncadp->getsize += cesize * (MPI_Offset)req->ecnt;

    return status;
}

/*
 * Process put request
 * If the request exists in active pool, we process it and return the id
 */
int ncadiosi_wait_get_req(NC_ad *ncadp, int reqid, int *stat)
{
    int err, status = NC_NOERR;
    NC_ad_get_list *lp = &(ncadp->getlist);

    /* Filter invalid reqid
     * Valid id range from 0 ~ nalloc - 1
     */
    if (reqid >= lp->nalloc || reqid < 0 || lp->pos[reqid] >= lp->nused) {
        status = NC_EINVAL_REQUEST;
    }
    else{
        /* Locate the req object, which is reqs[reqid] */
        status = ncadiosi_handle_get_req(ncadp, lp->reqs + reqid);

        /* Recycle req object to the pool */
        err = ncadiosi_get_list_remove(lp, reqid);
        if (err != NC_NOERR){
            return err;
        }
    }

    /* Return status to the user */
    if (stat != NULL) *stat = status;

    return NC_NOERR;
}

/*
 * Process all put request
 */
int ncadiosi_wait_all_get_req(NC_ad *ncadp) {
    int err, status = NC_NOERR;
    NC_ad_get_list *lp = &(ncadp->getlist);

    /* Search through req object array for object in use */
    while(lp->nused) {
        err = ncadiosi_wait_get_req(ncadp, lp->ids[lp->nused - 1], NULL);
        if (status == NC_NOERR) status = err;
    }

    return status;
}

/*
 * Initialize a read request structure and post corresponding ADIOS oepration
 */
int
ncadiosi_init_get_req( NC_ad *ncadp,
              NC_ad_get_req *r,
              ADIOS_VARINFO *v,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              void             *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype)
{
    int err;
    int i;
    size_t esize;
    int cesize;
    int sstart, scount, sstride;

    r->ready = 0;
    r->buf = buf;
    r->buftype = buftype;

    /* Calculate number of elements in single record */
    r->ecnt = 1;
    for(i = 0; i < v->ndim; i++){
        r->ecnt *= (size_t)count[i];
    }

    /* If user buffer is contiguous */
    if (imap == NULL){
        r->cbuf = r->buf;
    }
    else{
        err = ncmpii_create_imaptype(v->ndim, count, imap, buftype,
                                        &(r->imaptype));
        if (err != NC_NOERR) {
            return err;
        }
        MPI_Type_size(buftype, &cesize);
        r->cbsize = r->ecnt * (size_t)cesize;
        r->cbuf = NCI_Malloc(r->cbsize);
    }

    /* PnetCDF allows accessing in different type
     * Check if we need to convert
     */
    esize = (size_t)adios_type_size(v->type, NULL);
    r->vtype = ncadios_to_mpi_type(v->type);
    if (r->vtype == buftype){
        r->xbuf = r->cbuf;
    }
    else{
        r->xbuf = NCI_Malloc(esize * r->ecnt);
    }

    /* Time step dimension must be treated specially */
    if (v->nsteps > 1){
        sstart = (int)start[0];
        start++;
        scount = (int)count[0];
        count++;
        if (stride != NULL){
            sstride = (int)stride[0];
            stride++;
        }
        else{
            sstride = 1;
        }
    }
    else{
        sstart = 0;
        scount = 1;
        sstride = 1;
    }

    /* ADIOS selection */
    /* If stride is not used, we can use bounding box selection
     * Otherwise, we need to specify every points
     */
    if (stride == NULL){
        r->sel = adios_selection_boundingbox (v->ndim, (uint64_t*)start,
                                                (uint64_t*)count);
        r->points = NULL;
    }
    else{
        uint64_t *p, *cur;

        /* Somehow ADIOS doe not deep copy points, we need to keep it in the
         * request structure.
         */
        r->points = (uint64_t*)NCI_Malloc(sizeof(uint64_t) * r->ecnt * v->ndim);
        p = (uint64_t*)NCI_Malloc(sizeof(uint64_t) * v->ndim);
        cur = r->points;

        memset(p, 0, sizeof(uint64_t) * v->ndim);

        /* Iterate through every cells accessed */
        while(p[0] < count[0]) {
            for(i = 0; i < v->ndim; i++){
                *cur = p[i] * (uint64_t)stride[i];
                cur++;
            }

            p[v->ndim - 1]++;
            for(i = v->ndim - 1; i > 0; i--){
                if (p[i] >= count[i]){
                    p[i - 1]++;
                    p[i] = 0;
                }
                else{
                    break;
                }
            }
        }

        r->sel = adios_selection_points(v->ndim, (uint64_t)r->ecnt, r->points);

        NCI_Free(p);
    }
    if (r->sel == NULL){
        err = ncmpii_error_adios2nc(adios_errno, "select");
        DEBUG_RETURN_ERROR(err);
    }

    /* Post read operation */
    if (sstride > 1){
        /* ADIOS does not support stripe on time steps, post one step at a
         * time
         */
        for(i = 0; i < scount; i++){
            err = adios_schedule_read_byid (ncadp->fp, r->sel, v->varid,
                            sstart + i * sstride, 1,
                            (void*)(((char*)r->xbuf) + i * esize * r->ecnt));
            if (err != 0){
                err = ncmpii_error_adios2nc(adios_errno, "schedule_read");
                DEBUG_RETURN_ERROR(err);
            }
        }
    }
    else{
        err = adios_schedule_read_byid (ncadp->fp, r->sel, v->varid, sstart,
                                        scount, r->xbuf);
        if (err != 0){
            err = ncmpii_error_adios2nc(adios_errno, "schedule_read");
            DEBUG_RETURN_ERROR(err);
        }
    }

    return NC_NOERR;
}

/*
 * Non-blocking get
 * Initialize a read request structure
 * Obtain a request id from the pool
 * Put the request into list
 */
int
ncadiosi_iget_var(NC_ad *ncadp,
              int               varid,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              void             *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int *reqid)
{
    int err;
    int req_id;
    NC_ad_get_req r;
    ADIOS_VARINFO *v;

    /* Get ADIOS variable */
    v = adios_inq_var(ncadp->fp, ncadp->vars.data[varid].name);
    if (v == NULL){
        err = ncmpii_error_adios2nc(adios_errno, "get_var");
        DEBUG_RETURN_ERROR(err);
    }

    /* Create a read request */
    err = ncadiosi_init_get_req(ncadp, &r, v, start, count, stride, imap, buf,
                                bufcount, buftype);
    if (err != NC_NOERR){
        return err;
    }

    /* Release var info */
    adios_free_varinfo (v);

    /* Add to req list */
    ncadiosi_get_list_add(&(ncadp->getlist), &req_id);
    ncadp->getlist.reqs[req_id] = r;

    if (reqid != NULL){
        *reqid = req_id;
    }

    return NC_NOERR;
}
