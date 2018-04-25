/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs
 *
 * ncmpi_create()           : dispatcher->create()
 * ncmpi_open()             : dispatcher->open()
 * ncmpi_close()            : dispatcher->close()
 * ncmpi_enddef()           : dispatcher->enddef()
 * ncmpi__enddef()          : dispatcher->_enddef()
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_inq_misc()         : dispatcher->inq_misc()
 * ncmpi_wait()             : dispatcher->wait()
 * ncmpi_wait_all()         : dispatcher->wait()
 * ncmpi_cancel()           : dispatcher->cancel()
 *
 * ncmpi_set_fill()         : dispatcher->set_fill()
 * ncmpi_fill_var_rec()     : dispatcher->fill_rec()
 * ncmpi_def_var_fill()     : dispatcher->def_var_fill()
 * ncmpi_inq_var_fill()     : dispatcher->inq()
 *
 * ncmpi_sync()             : dispatcher->sync()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 *
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strlen() */

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncdwio_driver.h>

int
ncdwio_create(MPI_Comm     comm,
             const char  *path,
             int          cmode,
             int          ncid,
             MPI_Info     info,
             void       **ncpp)  /* OUT */
{
    int err;
    void *ncp=NULL;
    NC_dw *ncdwp;
    PNC_driver *driver=NULL;

    /* TODO: use comde to determine the true driver */
    driver = ncmpio_inq_driver();
    if (driver == NULL) return NC_ENOTNC;

    err = driver->create(comm, path, cmode, ncid, info, &ncp);
    if (err != NC_NOERR) return err;

    /* Create a NC_dw object and save its driver pointer */
    ncdwp = (NC_dw*) NCI_Malloc(sizeof(NC_dw));
    if (ncdwp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    ncdwp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (ncdwp->path == NULL) {
        NCI_Free(ncdwp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(ncdwp->path, path);
    ncdwp->mode = cmode;
    ncdwp->ncmpio_driver = driver;  // ncmpio driver
    ncdwp->ncid = ncid;
    ncdwp->isindep = 0; // Start at collective mode
    ncdwp->ncp = ncp;   // NC object used by ncmpio driver
    ncdwp->recdimsize = 0;
    ncdwp->recdimid = -1;   // Id of record dimension
    ncdwp->max_ndims = 0;   // Highest dimensionality among all variables
    MPI_Comm_dup(comm, &(ncdwp->comm));
    MPI_Info_dup(info, &(ncdwp->info));
    ncdwio_extract_hint(ncdwp, info);   // Translate MPI hint into hint flags

    /* Log init delayed to enddef */
    ncdwp->inited = 0;

    *ncpp = ncdwp;

    return NC_NOERR;
}

int
ncdwio_open(MPI_Comm     comm,
           const char  *path,
           int          omode,
           int          ncid,
           MPI_Info     info,
           void       **ncpp)
{
    int err, format;
    void *ncp=NULL;
    NC_dw *ncdwp;
    PNC_driver *driver=NULL;

    err = ncmpi_inq_file_format(path, &format);
    if (err != NC_NOERR) return err;

    if (format == NC_FORMAT_CLASSIC ||
        format == NC_FORMAT_CDF2 ||
        format == NC_FORMAT_CDF5) {
        driver = ncmpio_inq_driver();
    }
    if (driver == NULL) return NC_ENOTNC;

    err = driver->open(comm, path, omode, ncid, info, &ncp);
    if (err != NC_NOERR) return err;

    /* Create a NC_dw object and save its driver pointer */
    ncdwp = (NC_dw*) NCI_Malloc(sizeof(NC_dw));
    if (ncdwp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    ncdwp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (ncdwp->path == NULL) {
        NCI_Free(ncdwp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(ncdwp->path, path);
    ncdwp->mode = omode;
    ncdwp->ncmpio_driver = driver;  // ncmpio driver
    ncdwp->ncid = ncid;
    ncdwp->isindep = 0; // Start at collective mode
    ncdwp->ncp = ncp;   // NC object used by ncmpio driver
    ncdwp->recdimsize = 0;
    ncdwp->recdimid = -1;   // Id of record dimension
    ncdwp->max_ndims = 0;   // Highest dimensionality among all variables
    MPI_Comm_dup(comm, &(ncdwp->comm));
    MPI_Info_dup(info, &(ncdwp->info));
    ncdwio_extract_hint(ncdwp, info);   // Translate MPI hint into hint flags

    /* Opened file is in data mode
     * We must initialize the log for if file is not opened for read only
     */
    if (omode != NC_NOWRITE ){
        /* Init log file */
        err = ncdwio_log_create(ncdwp, info);
        if (err != NC_NOERR) {
            NCI_Free(ncdwp);
            return err;
        }
        // Initialize put list for nonblocking put operation
        ncdwio_put_list_init(ncdwp);
        if (err != NC_NOERR) {
            return err;
        }
        // Initialize metadata index for log entries
        ncdwio_metaidx_init(ncdwp);
        if (err != NC_NOERR) {
            return err;
        }
        // Mark as initialized
        ncdwp->inited = 1;
    }
    else{
        ncdwp->inited = 0;
    }
    *ncpp = ncdwp;

    return NC_NOERR;
}

/*
 * Close the netcdf file
 * The logfile must be closed
 */
int
ncdwio_close(void *ncdp)
{
    int err, status = NC_NOERR;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    if (ncdwp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    /* If log is initialized, we must close the log file
     * Putlist and metadata index also needs to be cleaned up
     */
    if (ncdwp->inited){
        // Close log file
        err = ncdwio_log_close(ncdwp);
        if (status == NC_NOERR) {
            status = err;
        }
        // Clean up put list
        ncdwio_put_list_free(ncdwp);
        // Clean up metadata index
        ncdwio_metaidx_free(ncdwp);
    }

    // Call ncmpio driver
    err = ncdwp->ncmpio_driver->close(ncdwp->ncp);
    if (status == NC_NOERR) {
        status = err;
    }

    // Cleanup NC-dw object
    MPI_Comm_free(&(ncdwp->comm));
    MPI_Info_free(&(ncdwp->info));
    NCI_Free(ncdwp->path);
    NCI_Free(ncdwp);

    return status;
}

/*
 * Handle enddef driver request
 * If the log file is not initialized yet, we must initialize it
 * Otherwise, we update the header of logfile with new information
 */
int
ncdwio_enddef(void *ncdp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    // Call ncmpio enddef
    err = ncdwp->ncmpio_driver->enddef(ncdwp->ncp);
    if (err != NC_NOERR) return err;

    /* If logfile are not initialized, we initialize the logfile
     * The file is newly created
     * If the log is already initialized, user must have called redef to define new structure
     * We update the information in logfile header accordingly
     */
    if (ncdwp->inited){
        /* Update log with new information */
        err = ncdwio_log_enddef(ncdwp);
        if (err != NC_NOERR) return err;
    }
    else {
        /* Init log file */
        err = ncdwio_log_create(ncdwp, ncdwp->info);
        if (err != NC_NOERR) {
            return err;
        }
        // Initialize put list for nonblocking put operation
        ncdwio_put_list_init(ncdwp);
        if (err != NC_NOERR) {
            return err;
        }
        // Initialize metadata index for log entries
        ncdwio_metaidx_init(ncdwp);
        if (err != NC_NOERR) {
            return err;
        }
        // Mark as initialized
        ncdwp->inited = 1;
    }

    return NC_NOERR;
}

/*
 * This function handle the enddef driver request
 * If the log file is not initialized yet, we must initialize it
 * Otherwise, we update the header of logfile with new information
 */
int
ncdwio__enddef(void       *ncdp,
              MPI_Offset  h_minfree,
              MPI_Offset  v_align,
              MPI_Offset  v_minfree,
              MPI_Offset  r_align)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    // Call ncmpio enddef
    err = ncdwp->ncmpio_driver->_enddef(ncdwp->ncp, h_minfree, v_align, v_minfree,
                               r_align);
    if (err != NC_NOERR) return err;

    /* If logfile are not initialized, we initialize the logfile
     * The file is newly created
     * If the log is already initialized, user must have called redef to define new structure
     * We update the information in logfile header accordingly
     */
    if (ncdwp->inited){
        /* Update log with new information */
        err = ncdwio_log_enddef(ncdwp);
        if (err != NC_NOERR) return err;
    }
    else {
        /* Init log file */
        err = ncdwio_log_create(ncdwp, ncdwp->info);
        if (err != NC_NOERR) {
            return err;
        }
        // Initialize put list for nonblocking put operation
        ncdwio_put_list_init(ncdwp);
        if (err != NC_NOERR) {
            return err;
        }
        // Initialize metadata index for log entries
        ncdwio_metaidx_init(ncdwp);
        if (err != NC_NOERR) {
            return err;
        }
        // Mark as initialized
        ncdwp->inited = 1;
    }

    return NC_NOERR;
}

int
ncdwio_redef(void *ncdp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    /*
     * Flush log entries to the file system on redefine
     */
    /*
    if (ncdwp->inited) {
        err = ncdwio_log_flush(ncdwp);
        if (err != NC_NOERR) return err;
    }
    */

    err = ncdwp->ncmpio_driver->redef(ncdwp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_begin_indep_data(void *ncdp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->begin_indep_data(ncdwp->ncp);
    if (err != NC_NOERR) return err;

    /* Independent mode
     * We keep track of current IO mode so we know what mode to use when flushing the log
     */
    ncdwp->isindep = 1;

    return NC_NOERR;
}

int
ncdwio_end_indep_data(void *ncdp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->end_indep_data(ncdwp->ncp);
    if (err != NC_NOERR) return err;

    /* Collective mode
     * We keep track of current IO mode so we know what mode to use when flushing the log
     */
    ncdwp->isindep = 0;

    return NC_NOERR;
}

int
ncdwio_abort(void *ncdp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    if (ncdwp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    err = ncdwp->ncmpio_driver->abort(ncdwp->ncp);

    MPI_Comm_free(&(ncdwp->comm));
    NCI_Free(ncdwp->path);
    NCI_Free(ncdwp);

    return err;
}

int
ncdwio_inq(void *ncdp,
          int  *ndimsp,
          int  *nvarsp,
          int  *nattsp,
          int  *xtendimp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->inq(ncdwp->ncp, ndimsp, nvarsp, nattsp, xtendimp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_inq_misc(void       *ncdp,
               int        *pathlen,
               char       *path,
               int        *num_fix_varsp,
               int        *num_rec_varsp,
               int        *striping_size,
               int        *striping_count,
               MPI_Offset *header_size,
               MPI_Offset *header_extent,
               MPI_Offset *recsize,
               MPI_Offset *put_size,
               MPI_Offset *get_size,
               MPI_Info   *info_used,
               int        *nreqs,
               MPI_Offset *usage,
               MPI_Offset *buf_size)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->inq_misc(ncdwp->ncp, pathlen, path, num_fix_varsp,
                                num_rec_varsp, striping_size, striping_count,
                                header_size, header_extent, recsize, put_size,
                                get_size, info_used, nreqs, usage, buf_size);
    if (err != NC_NOERR) return err;

    /* Data that is pending in the log is not counted by the ncmpio driver, we add the size of data in the log to put size
     * ncmpio driver does not handle put requests, we add number of pending put requests to ureqs
     */
    if (ncdwp->inited) {
        /* Add the size of data log to reflect pending put in the log */
        if (put_size != NULL){
            *put_size += (MPI_Offset)ncdwp->datalogsize - 8;

        }

        /* Add number of write requests to nreqs */
        if (nreqs != NULL){
            *nreqs += ncdwp->putlist.nused;

        }
    }

    /* Export dw related settings
     * ncdwio driver has it's own hints that is not handled by ncmpio driver
     */
    if (info_used != NULL){
        ncdwio_export_hint(ncdwp, *info_used);
    }

    return NC_NOERR;
}

/*
 * Cancel non-blocking request
 * IN       ncdp:    NC_dw object
 * IN    num_req:    Number of request to be canceled
 * IN    req_ids:    Reqest ids to be canceled
 * OUT  statuses:    Result of cancelation
 *
 * We only keep track of put requests, get requests are handled by the ncmpio driver
 * Given an array of request ids, we separate them into put and get request ids
 * Put requests are handled by the ncdwio driver, get requests are forwarded to the ncmpio driver
 * Put requests have odd ids, get request have even ids
 * We first reorder req_ids so that the first part contains all put requests and the second part contains only get requests
 * We call ncmpio_wait using the first part, then we handle the second part
 * After handling the requests, we restore the original order
 * We keep track of every swap operation we used to reorder req_ids
 * We apply them in reverse order to restore the original order
 */
int
ncdwio_cancel(void *ncdp,
             int   num_req,
             int  *req_ids,
             int  *statuses)
{
    int i, j, err, status = NC_NOERR;
    int tmp, stat;
    int nput;   // How many put request
    int *swapidx;   // Swap target
    NC_dw *ncdwp = (NC_dw*)ncdp;

   /*
    * If num_req is one of all requests, we don't need to handle request ids
    */
    if (num_req == NC_REQ_ALL || num_req == NC_PUT_REQ_ALL || num_req == NC_GET_REQ_ALL){
        // Cancel all put requests
        if (num_req == NC_REQ_ALL || num_req == NC_PUT_REQ_ALL){
            err = ncdwio_cancel_all_put_req(ncdwp);
            if (status == NC_NOERR){
                status = err;
            }
        }
        // Cancel all get requests
        if (num_req == NC_REQ_ALL || num_req == NC_GET_REQ_ALL){
            err = ncdwp->ncmpio_driver->cancel(ncdwp->ncp, num_req, NULL, NULL);
            if (status == NC_NOERR){
                status = err;
            }
        }

        return status;
    }

    /* Allocate buffer for tracking swap operation
     * swapidx stores the target location that swaps with current location
     * if swapidx[i] = j, it means the i-th entry is swapped with j-th entry in req_ids
     * We do onw swap for each put request, so there are at most num_req swaps
     */
    swapidx = (int*)NCI_Malloc(SIZEOF_INT * num_req);

    /* Count the number of put requests and swap it to the first section
     * nput is number of put request known so far, it also mark the end of the first section
     * When we find one put reqeust, we swap it with the entry at nput and increase nput by 1
     */
    nput = 0;
    for(i = 0; i < num_req; i++){
        if ((req_ids[i] & 1) == 0 || req_ids[i] == NC_REQ_NULL){    // Even id means a put request or NULL request
            // We are swapping req_ids[nput] with req_ids[i]
            swapidx[nput] = i;
            // Perform swap
            tmp = req_ids[i];
            req_ids[i] = req_ids[nput];
            req_ids[nput++] = tmp;
        }
    }

    /* If we have put requests
     * The internal put list uses a continuous id, so we translate it by dividing the id by 2
     */
    if (nput > 0){
        /* Cancel the request one by one */
        for(i = 0; i < nput; i++){
            // Skip NULL requests
            if (req_ids[i] == NC_REQ_NULL){
                stat = NC_NOERR;
            }
            else{
                /* Try canceling the request
                * Cancelation can fail if the request is already flushed to the file
                */
                err = ncdwio_cancel_put_req(ncdwp, (req_ids[i] / 2), &stat);
                if (status == NC_NOERR){
                    status = err;
                }
            }

            if (statuses != NULL){
                statuses[i] = stat;
            }
        }
    }

    /* If we have get requests
     * ncmpio driver has it's own request id management, so we translate it by dividing the id by 2
     */
    if (num_req > nput){
        // Translate reqid to ncmpio reqid
        for(i = nput; i < num_req; i++){
            req_ids[i] /= 2;
        }
        // Call ncmpio cancel
        err = ncdwp->ncmpio_driver->cancel(ncdwp->ncp, num_req - nput, req_ids + nput, statuses + nput);
        if (status == NC_NOERR){
            status = err;
        }
        // Translate reqid back to ncdwio reqid
        for(i = nput; i < num_req; i++){
            req_ids[i] *= 2;
        }
    }

    /* After processing the requests, we need to resotre the original order
     * We read swapsidx in reverse order
     * There are exactly nput swaps
     * Since req_ids[i] was swapped with req_ids[swapidx[i]], we must repeat it to swap it back
     */
    for(i = nput - 1; i > -1; i--){
        j = swapidx[i];
        tmp = req_ids[i];
        req_ids[i] = req_ids[j];
        req_ids[j] = tmp;
    }
    // the order of statuses must be restored as well since they are recorded after reordering
    if (statuses != NULL){
        for(i = nput - 1; i > -1; i--){
            j = swapidx[i];
            tmp = statuses[i];
            statuses[i] = statuses[j];
            statuses[j] = tmp;
        }
    }

    /* If no error happened, set req_ids to NC_REQ_NULL */
    if (status == NC_NOERR && req_ids != NULL){
        for(i = 0; i < num_req; i++){
            req_ids[i] = NC_REQ_NULL;
        }
    }

    // Free the tracking buffer
    NCI_Free(swapidx);

    return status;
}

/*
 * Handle non-blocking request
 * IN       ncdp:    NC_dw object
 * IN    num_req:    Number of request to be canceled
 * IN    req_ids:    Reqest ids to be canceled
 * OUT  statuses:    Result of cancelation
 *
 * We only keep track of put requests, get requests are handled by the ncmpio driver
 * Given an array of request ids, we separate them into put and get request ids
 * Put requests are handled by the ncdwio driver, get requests are forwarded to the ncmpio driver
 * Put requests have odd ids, get request have even ids
 * We first reorder req_ids so that the first part contains all put requests and the second part contains only get requests
 * We call ncmpio_wait using the first part, then we handle the second part
 * After handling the requests, we restore the original order
 * We keep track of every swap operation we used to reorder req_ids
 * We apply them in reverse order to restore the original order
 *
 * Here we give an example of reordering:
 * Initial:
 * req_ids = 0 1 2 3 4
 * After reordering:
 *      put reqs|get reqs
 * req_ids = 1 3 0 2 4
 *               ^
 *              nput = 2
 * swapidx = 1 3
 */
int
ncdwio_wait(void *ncdp,
           int   num_reqs,
           int  *req_ids,
           int  *statuses,
           int   reqMode)
{
    int i, j, err, status = NC_NOERR;
    int tmp, stat;
    int nput;   // How many put request
    int *swapidx;   // Swap target
    NC_dw *ncdwp = (NC_dw*)ncdp;

   /*
    * If num_reqs is one of all requests, we don't need to handle request ids
    */
    if (num_reqs == NC_REQ_ALL || num_reqs == NC_PUT_REQ_ALL || num_reqs == NC_GET_REQ_ALL){
        // Cancel all put requests
        if (num_reqs == NC_REQ_ALL || num_reqs == NC_PUT_REQ_ALL){
            err = ncdwio_handle_all_put_req(ncdwp);
            if (status == NC_NOERR){
                status = err;
            }
        }
        // Cancel all get requests
        if (num_reqs == NC_REQ_ALL || num_reqs == NC_GET_REQ_ALL){
            err = ncdwp->ncmpio_driver->wait(ncdwp->ncp, num_reqs, NULL, NULL, reqMode);
            if (status == NC_NOERR){
                status = err;
            }
        }

        return status;
    }

    /* Allocate buffer for tracking swap operation
     * swapidx stores the target location that swaps with current location
     * if swapidx[i] = j, it means the i-th entry is swapped with j-th entry in req_ids
     * We do onw swap for each put request, so there are at most num_reqs swaps
     */
    swapidx = (int*)NCI_Malloc(SIZEOF_INT * num_reqs);

    /* Count the number of put requests and swap it to the first section
     * nput is number of put request known so far, it also mark the end of the first section
     * When we find one put reqeust, we swap it with the entry at nput and increase nput by 1
     */
    nput = 0;
    for(i = 0; i < num_reqs; i++){
        if ((req_ids[i] & 1) == 0 || req_ids[i] == NC_REQ_NULL){    // Even id means a put request or NULL request
            // We are swapping req_ids[nput] with req_ids[i]
            swapidx[nput] = i;
            // Perform swap
            tmp = req_ids[i];
            req_ids[i] = req_ids[nput];
            req_ids[nput++] = tmp;
        }
    }

    /* If we have put requests
     * The internal put list uses a continuous id, so we translate it by dividing the id by 2
     */
    if (nput > 0){
        /* Handle the request one by one */
        for(i = 0; i < nput; i++){
            // Handle request, skiping NULL requests
            if (req_ids[i] == NC_REQ_NULL){
                stat = NC_NOERR;
            }
            else{
                // Waiting can fail if there's problem writing request to file
                err = ncdwio_handle_put_req(ncdwp, (req_ids[i] / 2), &stat);
                if (status == NC_NOERR){
                    status = err;
                }
            }

            if (statuses != NULL){
                statuses[i] = stat;
            }
        }
    }

    /* If we have get requests
     * ncmpio driver has it's own request id management, so we translate it by dividing the id by 2
     * We need to flush the log so new data can be read
     */
    if (num_reqs > nput){
        // Flush the log if log is initialized
        if (ncdwp->inited){
            err = ncdwio_log_flush(ncdwp);
            if (status == NC_NOERR){
                status = err;
            }
        }
        // Translate reqid to ncmpio reqid
        for(i = nput; i < num_reqs; i++){
            req_ids[i] /= 2;
        }
        // Call ncmpio wait
        err = ncdwp->ncmpio_driver->wait(ncdwp->ncp, num_reqs - nput, req_ids + nput, statuses + nput, reqMode);
        if (status == NC_NOERR){
            status = err;
        }
        // Translate reqid back to ncdwio reqid
        for(i = nput; i < num_reqs; i++){
            req_ids[i] *= 2;
        }
    }

    /* After processing the requests, we need to resotre the original order
     * We read swapsidx in reverse order
     * There are exactly nput swaps
     * Since req_ids[i] was swapped with req_ids[swapidx[i]], we must repeat it to swap it back
     */
    for(i = nput - 1; i > -1; i--){
        j = swapidx[i];
        tmp = req_ids[i];
        req_ids[i] = req_ids[j];
        req_ids[j] = tmp;
    }
    // the order of statuses must be restored as well since they are recorded after reordering
    if (statuses != NULL){
        for(i = nput - 1; i > -1; i--){
            j = swapidx[i];
            tmp = statuses[i];
            statuses[i] = statuses[j];
            statuses[j] = tmp;
        }
    }

    /* If no error happened, set req_ids to NC_REQ_NULL */
    if (status == NC_NOERR && req_ids != NULL){
        for(i = 0; i < num_reqs; i++){
            req_ids[i] = NC_REQ_NULL;
        }
    }

    // Free the tracking buffer
    NCI_Free(swapidx);

    return status;
}

int
ncdwio_set_fill(void *ncdp,
               int   fill_mode,
               int  *old_fill_mode)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->set_fill(ncdwp->ncp, fill_mode, old_fill_mode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_fill_var_rec(void      *ncdp,
                   int        varid,
                   MPI_Offset recno)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->fill_var_rec(ncdwp->ncp, varid, recno);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_def_var_fill(void       *ncdp,
                   int         varid,
                   int         no_fill,
                   const void *fill_value)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->def_var_fill(ncdwp->ncp, varid, no_fill, fill_value);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_sync_numrecs(void *ncdp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->sync_numrecs(ncdwp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_sync(void *ncdp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    /* Flush on sync */
    if (ncdwp->inited) {
        err = ncdwio_log_flush(ncdwp);
        if (err != NC_NOERR) return err;
    }

    err = ncdwp->ncmpio_driver->sync(ncdwp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

