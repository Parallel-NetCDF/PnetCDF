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
 * ncmpi_flush()            : dispatcher->flush()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 *
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <string.h> /* strlen() */
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncbbio_driver.h>

int
ncbbio_create(MPI_Comm     comm,
              const char  *path,
              int          cmode,
              int          ncid,
              MPI_Info     info,
              void       **ncpp)  /* OUT */
{
    int err;
    void *ncp=NULL;
    NC_bb *ncbbp;
    PNC_driver *driver=NULL;

    /* TODO: use cmode to determine the true driver */
    driver = ncmpio_inq_driver();
    if (driver == NULL) DEBUG_RETURN_ERROR(NC_ENOTNC)

    err = driver->create(comm, path, cmode, ncid, info, &ncp);
    if (err != NC_NOERR) return err;

    /* Create a NC_bb object and save its driver pointer */
    ncbbp = (NC_bb*) NCI_Malloc(sizeof(NC_bb));
    if (ncbbp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    ncbbp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (ncbbp->path == NULL) {
        NCI_Free(ncbbp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(ncbbp->path, path);
    ncbbp->mode          = cmode;
    ncbbp->ncmpio_driver = driver;  /* ncmpio driver */
    ncbbp->ncid          = ncid;
    ncbbp->ncp           = ncp;   /* NC object used by ncmpio driver */
    ncbbp->recdimsize    = 0;
    ncbbp->recdimid      = -1;   /* Id of record dimension */
    ncbbp->max_ndims     = 0;   /* Highest dimensionality among all variables */
    ncbbp->datalog_fd    = NULL;
    ncbbp->metalog_fd    = NULL;
    ncbbp->flag          = NC_MODE_CREATE | NC_MODE_DEF;
    ncbbp->logcomm       = MPI_COMM_SELF;
    ncbbp->comm          = comm;  /* reuse comm duplicated in dispatch layer */

    driver->inq_misc(ncp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, &ncbbp->info, NULL, NULL, NULL);
    ncbbio_extract_hint(ncbbp, info); /* extract bb related hints */

    /* Log init delayed to enddef */
    ncbbp->inited = 0;

    *ncpp = ncbbp;

    return NC_NOERR;
}

int
ncbbio_open(MPI_Comm     comm,
            const char  *path,
            int          omode,
            int          ncid,
            MPI_Info     info,
            void       **ncpp)
{
    int err;
    void *ncp=NULL;
    NC_bb *ncbbp;
    PNC_driver *driver=NULL;

    driver = ncmpio_inq_driver();
    if (driver == NULL) DEBUG_RETURN_ERROR(NC_ENOTNC)

    err = driver->open(comm, path, omode, ncid, info, &ncp);
    if (err != NC_NOERR) return err;

    /* Create a NC_bb object and save its driver pointer */
    ncbbp = (NC_bb*) NCI_Malloc(sizeof(NC_bb));
    if (ncbbp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    ncbbp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (ncbbp->path == NULL) {
        NCI_Free(ncbbp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(ncbbp->path, path);
    ncbbp->mode          = omode;
    ncbbp->ncmpio_driver = driver;  /* ncmpio driver */
    ncbbp->ncid          = ncid;
    ncbbp->ncp           = ncp;   /* NC object used by ncmpio driver */
    ncbbp->recdimsize    = 0;
    ncbbp->recdimid      = -1;   /* Id of record dimension */
    ncbbp->max_ndims     = 0;   /* Highest dimensionality among all variables */
    ncbbp->flag          = 0;
    ncbbp->logcomm       = MPI_COMM_SELF;
    ncbbp->comm          = comm;  /* reuse comm duplicated in dispatch layer */

    driver->inq_misc(ncp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, &ncbbp->info, NULL, NULL, NULL);
    ncbbio_extract_hint(ncbbp, info); /* extract bb related hints */

    /* Opened file is in data mode
     * We must initialize the log for if file is not opened for read only
     */
    if (omode != NC_NOWRITE) {
        /* Init log file */
        err = ncbbio_log_create(ncbbp, ncbbp->info);
        if (err != NC_NOERR) {
            NCI_Free(ncbbp);
            return err;
        }
        /* Initialize put list for nonblocking put operation */
        err = ncbbio_put_list_init(ncbbp);
        if (err != NC_NOERR) {
            NCI_Free(ncbbp);
            return err;
        }
        /* Initialize metadata index for log entries */
        err = ncbbio_metaidx_init(ncbbp);
        if (err != NC_NOERR) {
            NCI_Free(ncbbp);
            return err;
        }
        /* Mark as initialized */
        ncbbp->inited = 1;
    }
    else {
        fSet(ncbbp->flag, NC_MODE_RDONLY);
        ncbbp->inited = 0;
    }
    *ncpp = ncbbp;

    return NC_NOERR;
}

/*
 * Close the netcdf file
 * The logfile must be closed
 */
int
ncbbio_close(void *ncdp)
{
    int status=NC_NOERR, err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    if (ncbbp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    /* If log is initialized, we must close the log file
     * Putlist and metadata index also needs to be cleaned up
     */
    if (ncbbp->inited) {
        /* Close log file */
        err = ncbbio_log_close(ncbbp, 1);
        if (status == NC_NOERR) status = err;
        /* Clean up put list */
        err = ncbbio_put_list_free(ncbbp);
        if (status == NC_NOERR) status = err;
        /* Clean up metadata index */
        err = ncbbio_metaidx_free(ncbbp);
        if (status == NC_NOERR) status = err;
    }

    /* Call ncmpio driver */
    err = ncbbp->ncmpio_driver->close(ncbbp->ncp);
    if (status == NC_NOERR) status = err;

    /* Cleanup NC-bb object */
    if (ncbbp->info != MPI_INFO_NULL)
        MPI_Info_free(&ncbbp->info);
    NCI_Free(ncbbp->path);
    NCI_Free(ncbbp);

    return status;
}

/*
 * Handle enddef driver request
 * If the log file is not initialized yet, we must initialize it
 * Otherwise, we update the header of logfile with new information
 */
int
ncbbio_enddef(void *ncdp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Call ncmpio enddef */
    err = ncbbp->ncmpio_driver->enddef(ncbbp->ncp);
    if (err != NC_NOERR) return err;

    /* If logfile are not initialized, we initialize the logfile. The file is
     * newly created.  If the log is already initialized, user must have called
     * redef to define new structure.  We update the information in logfile
     * header accordingly.
     */
    if (ncbbp->inited)
        /* Update log with new information */
        err = ncbbio_log_enddef(ncbbp);
    else
        err = ncbbio_init(ncbbp);

    if (err != NC_NOERR) return err;

    fClr(ncbbp->flag, NC_MODE_DEF);

    return NC_NOERR;
}

/*
 * Initialize the driver for data handling
 * Must be called on before first data operation, including nonblocking request
 * Initialize the log
 * Initialize metadata index
 * Initialize nonblocking request pool
 */
int
ncbbio_init(NC_bb *ncbbp)
{
    int err;

    /* If logfile are not initialized, we initialize the logfile */
    if (!ncbbp->inited) {
        /* Init log file */
        err = ncbbio_log_create(ncbbp, ncbbp->info);
        if (err != NC_NOERR) return err;

        /* Initialize put list for nonblocking put operation */
        err = ncbbio_put_list_init(ncbbp);
        if (err != NC_NOERR) return err;

        /* Initialize metadata index for log entries */
        err = ncbbio_metaidx_init(ncbbp);
        if (err != NC_NOERR) return err;

        /* Mark as initialized */
        ncbbp->inited = 1;
    }

    return NC_NOERR;
}

/*
 * This function handle the enddef driver request
 * If the log file is not initialized yet, we must initialize it
 * Otherwise, we update the header of logfile with new information
 */
int
ncbbio__enddef(void       *ncdp,
               MPI_Offset  h_minfree,
               MPI_Offset  v_align,
               MPI_Offset  v_minfree,
               MPI_Offset  r_align)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Call ncmpio enddef */
    err = ncbbp->ncmpio_driver->_enddef(ncbbp->ncp, h_minfree, v_align,
                                        v_minfree, r_align);
    if (err != NC_NOERR) return err;

    /* If logfile are not initialized, we initialize the logfile.  The file is
     * newly created.  If the log is already initialized, user must have called
     * redef to define new structure.  We update the information in logfile
     * header accordingly.
     */
    if (ncbbp->inited) {
        /* Update log with new information */
        err = ncbbio_log_enddef(ncbbp);
        if (err != NC_NOERR) return err;
    }
    else {
        /* Init log file */
        err = ncbbio_log_create(ncbbp, ncbbp->info);
        if (err != NC_NOERR) return err;

        /* Initialize put list for nonblocking put operation */
        err = ncbbio_put_list_init(ncbbp);
        if (err != NC_NOERR) return err;

        /* Initialize metadata index for log entries */
        err = ncbbio_metaidx_init(ncbbp);
        if (err != NC_NOERR) return err;

        /* Mark as initialized */
        ncbbp->inited = 1;
    }

    fClr(ncbbp->flag, NC_MODE_DEF);

    return NC_NOERR;
}

int
ncbbio_redef(void *ncdp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /*
     * Flush log entries to the file system on redefine.  After redef, our
     * record in the log can become out dated due to change in varid, dimsize
     * ... etc.  Flush the log to ensure we have a fresh start.
     */

    if (ncbbp->inited) {
        err = ncbbio_log_flush(ncbbp);
        if (err != NC_NOERR) return err;
    }

    err = ncbbp->ncmpio_driver->redef(ncbbp->ncp);
    if (err != NC_NOERR) return err;

    fSet(ncbbp->flag, NC_MODE_DEF);

    return NC_NOERR;
}

int
ncbbio_begin_indep_data(void *ncdp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->begin_indep_data(ncbbp->ncp);
    if (err != NC_NOERR) return err;

    /* Independent mode
     * We keep track of current IO mode so we know what mode to use when
     * flushing the log
     */
    fSet(ncbbp->flag, NC_MODE_INDEP);

    return NC_NOERR;
}

int
ncbbio_end_indep_data(void *ncdp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->end_indep_data(ncbbp->ncp);
    if (err != NC_NOERR) return err;

    /* Collective mode
     * We keep track of current IO mode so we know what mode to use when
     * flushing the log
     */
    fClr(ncbbp->flag, NC_MODE_INDEP);

    return NC_NOERR;
}

int
ncbbio_abort(void *ncdp)
{
    int status=NC_NOERR, err;
    int replay = 1;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    if (ncbbp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    /* If log is initialized, we must close the log file
     * Putlist and metadata index also needs to be cleaned up
     */
    if (ncbbp->inited) {
        /* No flushing on abort if in define mode */
        if (fIsSet(ncbbp->flag, NC_MODE_DEF))
            replay = 0;

        /* Close log file */
        err = ncbbio_log_close(ncbbp, replay);
        if (status == NC_NOERR) status = err;
        /* Clean up put list */
        err = ncbbio_put_list_free(ncbbp);
        if (status == NC_NOERR) status = err;
        /* Clean up metadata index */
        err = ncbbio_metaidx_free(ncbbp);
        if (status == NC_NOERR) status = err;
    }

    /* Call ncmpio driver */
    err = ncbbp->ncmpio_driver->abort(ncbbp->ncp);
    if (status == NC_NOERR) status = err;

    if (ncbbp->info != MPI_INFO_NULL)
        MPI_Info_free(&ncbbp->info);
    NCI_Free(ncbbp->path);
    NCI_Free(ncbbp);

    return status;
}

int
ncbbio_inq(void *ncdp,
           int  *ndimsp,
           int  *nvarsp,
           int  *nattsp,
           int  *xtendimp)
{
    NC_bb *ncbbp = (NC_bb*)ncdp;

    return ncbbp->ncmpio_driver->inq(ncbbp->ncp, ndimsp, nvarsp, nattsp,
                                     xtendimp);
}

int
ncbbio_inq_misc(void       *ncdp,
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
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->inq_misc(ncbbp->ncp, pathlen, path,
                                         num_fix_varsp, num_rec_varsp,
                                         striping_size, striping_count,
                                         header_size, header_extent, recsize,
                                         put_size, get_size, info_used, nreqs,
                                         usage, buf_size);
    if (err != NC_NOERR) return err;

    /* Data that is pending in the log is not counted by the ncmpio driver, we
     * add the size of data in the log to put size.  ncmpio driver does not
     * handle put requests, we add number of pending put requests to nreqs
     */
    if (ncbbp->inited) {
        /* Add the size of data log to reflect pending put in the log */
        if (put_size != NULL)
            *put_size += (MPI_Offset)ncbbp->datalogsize - 8;

        /* Add number of write requests to nreqs */
        if (nreqs != NULL)
            *nreqs += ncbbp->putlist.nused;
    }

    /* Export bb related settings
     * ncbbio driver has it's own hints that is not handled by ncmpio driver
     */
    if (info_used != NULL)
        ncbbio_export_hint(ncbbp, info_used);

    return NC_NOERR;
}

/*
 * Cancel non-blocking request
 * IN       ncdp:    NC_bb object
 * IN    num_req:    Number of request to be canceled
 * IN    req_ids:    Reqest ids to be canceled
 * OUT  statuses:    Result of cancellation
 *
 * We only keep track of put requests, get requests are handled by the ncmpio
 * driver.  Given an array of request ids, we separate them into put and get
 * request ids.  Put requests are handled by the ncbbio driver, get requests
 * are forwarded to the ncmpio driver.  Put requests have odd ids, get request
 * have even ids.  We first reorder req_ids so that the first part contains all
 * put requests and the second part contains only get requests.  We call
 * ncmpio_wait using the first part, then we handle the second part.  After
 * handling the requests, we restore the original order.  We keep track of
 * every swap operation we used to reorder req_ids.  We apply them in reverse
 * order to restore the original order.
 */
int
ncbbio_cancel(void *ncdp,
              int   num_req,
              int  *req_ids,
              int  *statuses)
{
    int i, j, err, status = NC_NOERR;
    int tmp, stat;
    int nput;   /* How many put request */
    int *swapidx;   /* Swap target */
    NC_bb *ncbbp = (NC_bb*)ncdp;

   /*
    * If num_req is one of all requests, we don't need to handle request ids
    */
    if (num_req == NC_REQ_ALL || num_req == NC_PUT_REQ_ALL ||
                                 num_req == NC_GET_REQ_ALL) {
        /* Cancel all put requests */
        if (num_req == NC_REQ_ALL || num_req == NC_PUT_REQ_ALL) {
            /* Only handle put reqs in wr mode */
            if (ncbbp->inited) {
                err = ncbbio_cancel_all_put_req(ncbbp);
                if (status == NC_NOERR) status = err;
            }
        }
        /* Cancel all get requests */
        if (num_req == NC_REQ_ALL || num_req == NC_GET_REQ_ALL) {
            err = ncbbp->ncmpio_driver->cancel(ncbbp->ncp, num_req, NULL, NULL);
            if (status == NC_NOERR) status = err;
        }
        return status;
    }

    /* Allocate buffer for tracking swap operation.
     * swapidx stores the target location that swaps with current location.  if
     * swapidx[i] = j, it means the i-th entry is swapped with j-th entry in
     * req_ids.  We do onw swap for each put request, so there are at most
     * num_req swaps.
     */
    swapidx = (int*)NCI_Malloc(SIZEOF_INT * num_req);

    /* Count the number of put requests and swap it to the first section.  nput
     * is number of put request known so far, it also mark the end of the first
     * section.  When we find one put request, we swap it with the entry at
     * nput and increase nput by 1.
     */
    nput = 0;
    for (i=0; i<num_req; i++) {
        if ((req_ids[i] & 1) == 0 || req_ids[i] == NC_REQ_NULL) {
            /* Even id means a put request or NULL request */
            /* We are swapping req_ids[nput] with req_ids[i] */
            swapidx[nput] = i;
            /* Perform swap */
            tmp = req_ids[i];
            req_ids[i] = req_ids[nput];
            req_ids[nput++] = tmp;
        }
    }

    /* If we have put requests
     * The internal put list uses a continuous id, so we translate it by
     * dividing the id by 2
     */
    if (nput > 0) {
        /* Only handle put reqs in wr mode */
        if (ncbbp->inited) {
            for (i=0; i<nput; i++) {
                /* Skip NULL requests */
                if (req_ids[i] == NC_REQ_NULL)
		    stat = NC_NOERR;
                else {
                    /* Try canceling the request
		    * Cancellation can fail if the request is already flushed to
		    * the file
                    */
                    err = ncbbio_cancel_put_req(ncbbp, (req_ids[i] / 2), &stat);
                    if (status == NC_NOERR) status = err;
                }

                if (statuses != NULL) statuses[i] = stat;
            }
        }
        else {
            for (i=0; i<nput; i++) {
                /* Any put req other than NULL req is invalid in rd only mode */
                if (req_ids[i] == NC_REQ_NULL)
                    stat = NC_NOERR;
                else
		    /* Waiting can fail if there's problem writing request to */
		    /* file */
		    stat = NC_EINVAL_REQUEST;

                if (statuses != NULL) statuses[i] = stat;
            }
        }
    }

    /* If we have get requests
     * ncmpio driver has it's own request id management, so we translate it by
     * dividing the id by 2
     */
    if (num_req > nput) {
        /* Translate reqid to ncmpio reqid */
        for (i=nput; i<num_req; i++)
            req_ids[i] /= 2;

        /* Call ncmpio cancel */
        err = ncbbp->ncmpio_driver->cancel(ncbbp->ncp, num_req - nput,
	                                   req_ids + nput, statuses + nput);
        if (status == NC_NOERR) status = err;

        /* Translate reqid back to ncbbio reqid */
        for (i=nput; i<num_req; i++)
            req_ids[i] *= 2;
    }

    /* After processing the requests, we need to restore the original order
     * We read swapsidx in reverse order
     * There are exactly nput swaps
     * Since req_ids[i] was swapped with req_ids[swapidx[i]], we must repeat it
     * to swap it back
     */
    for (i=nput-1; i>-1; i--) {
        j = swapidx[i];
        tmp = req_ids[i];
        req_ids[i] = req_ids[j];
        req_ids[j] = tmp;
    }
    /* the order of statuses must be restored as well since they are recorded */
    /* after reordering */
    if (statuses != NULL) {
        for (i=nput-1; i>-1; i--) {
            j = swapidx[i];
            tmp = statuses[i];
            statuses[i] = statuses[j];
            statuses[j] = tmp;
        }
    }

    /* If no error happened, set req_ids to NC_REQ_NULL */
    if (status == NC_NOERR && req_ids != NULL) {
        for (i=0; i<num_req; i++)
            req_ids[i] = NC_REQ_NULL;
    }

    /* Free the tracking buffer */
    NCI_Free(swapidx);

    return status;
}

/*
 * Handle non-blocking request
 * IN       ncdp:    NC_bb object
 * IN    num_req:    Number of request to be canceled
 * IN    req_ids:    Reqest ids to be canceled
 * OUT  statuses:    Result of cancellation
 *
 * We only keep track of put requests, get requests are handled by the ncmpio
 * driver.  Given an array of request ids, we separate them into put and get
 * request ids.  Put requests are handled by the ncbbio driver, get requests
 * are forwarded to the ncmpio driver.  Put requests have odd ids, get request
 * have even ids.  We first reorder req_ids so that the first part contains all
 * put requests and the second part contains only get requests.  We call
 * ncmpio_wait using the first part, then we handle the second part.  After
 * handling the requests, we restore the original order.  We keep track of
 * every swap operation we used to reorder req_ids.  We apply them in reverse
 * order to restore the original order.
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
ncbbio_wait(void *ncdp,
            int   num_reqs,
            int  *req_ids,
            int  *statuses,
            int   reqMode)
{
    int i, j, err, status = NC_NOERR;
    int tmp, stat;
    int nput;   /* How many put request */
    int *swapidx;   /* Swap target */
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Flush the log if log is initialized */
    if (ncbbp->inited) {
        err = ncbbio_log_flush(ncbbp);
        if (status == NC_NOERR) status = err;
    }

   /*
    * If num_reqs is one of all requests, we don't need to handle request ids
    */
    if (num_reqs == NC_REQ_ALL || num_reqs == NC_PUT_REQ_ALL ||
                                  num_reqs == NC_GET_REQ_ALL) {
        /* Wait all put requests */
        if (num_reqs == NC_REQ_ALL || num_reqs == NC_PUT_REQ_ALL) {
            /* Only handle put reqs in wr mode */
            if (ncbbp->inited) {
                err = ncbbio_handle_all_put_req(ncbbp);
                if (status == NC_NOERR) status = err;
            }
        }
        /* Wait all get requests */
        if (num_reqs == NC_REQ_ALL || num_reqs == NC_GET_REQ_ALL) {
            err = ncbbp->ncmpio_driver->wait(ncbbp->ncp, num_reqs, NULL, NULL,
	                                     reqMode);
            if (status == NC_NOERR) status = err;
        }

        return status;
    }

    /* Allocate buffer for tracking swap operation
     * swapidx stores the target location that swaps with current location
     * if swapidx[i] = j, it means the i-th entry is swapped with j-th entry in
     * req_ids.  We do onw swap for each put request, so there are at most
     * num_reqs swaps
     */
    swapidx = (int*)NCI_Malloc(SIZEOF_INT * num_reqs);

    /* Count the number of put requests and swap it to the first section
     * nput is number of put request known so far, it also mark the end of the
     * first section When we find one put request, we swap it with the entry at
     * nput and increase nput by 1
     */
    nput = 0;
    for (i=0; i<num_reqs; i++) {
        if ((req_ids[i] & 1) == 0 || req_ids[i] == NC_REQ_NULL) {
	    /* Even id means a put request or NULL request */
            /* We are swapping req_ids[nput] with req_ids[i] */
            swapidx[nput] = i;
            /* Perform swap */
            tmp = req_ids[i];
            req_ids[i] = req_ids[nput];
            req_ids[nput++] = tmp;
        }
    }

    /* If we have put requests
     * The internal put list uses a continuous id, so we translate it by
     * dividing the id by 2
     */
    if (nput > 0) {
        /* Handle requests if in wr mode */
        if (ncbbp->inited) {
            /* Handle the request one by one */
            for (i=0; i<nput; i++) {
                /* Handle request, skipping NULL requests */
                if (req_ids[i] == NC_REQ_NULL)
                    stat = NC_NOERR;
                else {
                    /* Wait can fail if there's problem writing request to file */
                    err = ncbbio_handle_put_req(ncbbp, (req_ids[i] / 2), &stat);
                    if (status == NC_NOERR) status = err;
                }

                if (statuses != NULL) statuses[i] = stat;
            }
        }
        else {
            for (i=0; i<nput; i++) {
                /* Any put req other than NULL req is invalid in rd only mode */
                if (req_ids[i] == NC_REQ_NULL)
                    stat = NC_NOERR;
                else
                    /* Wait can fail if there's problem writing request to file */
                    stat = NC_EINVAL_REQUEST;

                if (statuses != NULL)
                    statuses[i] = stat;
            }
        }
    }

    /* If we have get requests
     * ncmpio driver has it's own request id management, so we translate it by
     * dividing the id by 2.  We need to flush the log so new data can be read
     */
    if (num_reqs > nput || !(fIsSet(ncbbp->flag, NC_MODE_INDEP))) {
        /* Translate reqid to ncmpio reqid */
        for (i=nput; i<num_reqs; i++)
            req_ids[i] /= 2;

        /* Call ncmpio wait */
        err = ncbbp->ncmpio_driver->wait(ncbbp->ncp, num_reqs - nput,
	                                 req_ids + nput, statuses + nput,
					 reqMode);
        if (status == NC_NOERR) status = err;

        /* Translate reqid back to ncbbio reqid */
        for (i=nput; i<num_reqs; i++) req_ids[i] *= 2;
    }

    /* After processing the requests, we need to restore the original order
     * We read swapsidx in reverse order
     * There are exactly nput swaps
     * Since req_ids[i] was swapped with req_ids[swapidx[i]], we must repeat it
     * to swap it back
     */
    for (i=nput-1; i>-1; i--) {
        j = swapidx[i];
        tmp = req_ids[i];
        req_ids[i] = req_ids[j];
        req_ids[j] = tmp;
    }
    /* the order of statuses must be restored as well since they are recorded */
    /* after reordering */
    if (statuses != NULL) {
        for (i=nput-1; i>-1; i--) {
            j = swapidx[i];
            tmp = statuses[i];
            statuses[i] = statuses[j];
            statuses[j] = tmp;
        }
    }

    /* If no error happened, set req_ids to NC_REQ_NULL */
    if (status == NC_NOERR && req_ids != NULL) {
        for (i=0; i<num_reqs; i++)
            req_ids[i] = NC_REQ_NULL;
    }

    /* Free the tracking buffer */
    NCI_Free(swapidx);

    return status;
}

int
ncbbio_set_fill(void *ncdp,
                int   fill_mode,
                int  *old_fill_mode)
{
    NC_bb *ncbbp = (NC_bb*)ncdp;

    return ncbbp->ncmpio_driver->set_fill(ncbbp->ncp, fill_mode, old_fill_mode);
}

int
ncbbio_fill_var_rec(void      *ncdp,
                    int        varid,
                    MPI_Offset recno)
{
    NC_bb *ncbbp = (NC_bb*)ncdp;

    return ncbbp->ncmpio_driver->fill_var_rec(ncbbp->ncp, varid, recno);
}

int
ncbbio_def_var_fill(void       *ncdp,
                    int         varid,
                    int         no_fill,
                    const void *fill_value)
{
    NC_bb *ncbbp = (NC_bb*)ncdp;

    return ncbbp->ncmpio_driver->def_var_fill(ncbbp->ncp, varid, no_fill,
                                              fill_value);
}

int
ncbbio_sync_numrecs(void *ncdp)
{
    NC_bb *ncbbp = (NC_bb*)ncdp;

    return ncbbp->ncmpio_driver->sync_numrecs(ncbbp->ncp);
}

int
ncbbio_sync(void *ncdp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Flush on sync */
    if (ncbbp->inited) {
        err = ncbbio_log_flush(ncbbp);
        if (err != NC_NOERR) return err;
    }

    err = ncbbp->ncmpio_driver->sync(ncbbp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_flush(void *ncdp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Flush the burst buffer */
    if (ncbbp->inited) {
        err = ncbbio_log_flush(ncbbp);
        if (err != NC_NOERR) return err;
    }

    err = ncbbp->ncmpio_driver->flush(ncbbp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

