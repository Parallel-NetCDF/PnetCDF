/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_def_var()                  : dispatcher->def_var()
 * ncmpi_inq_varid()                : dispatcher->inq_varid()
 * ncmpi_inq_var()                  : dispatcher->inq_var()
 * ncmpi_rename_var()               : dispatcher->rename_var()
 *
 * ncmpi_get_var<kind>()            : dispatcher->get_var()
 * ncmpi_put_var<kind>()            : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>()     : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>()     : dispatcher->put_var()
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 *
 * ncmpi_iget_var<kind>()           : dispatcher->iget_var()
 * ncmpi_iput_var<kind>()           : dispatcher->iput_var()
 * ncmpi_iget_var<kind>_<type>()    : dispatcher->iget_var()
 * ncmpi_iput_var<kind>_<type>()    : dispatcher->iput_var()
 *
 * ncmpi_buffer_attach()            : dispatcher->buffer_attach()
 * ncmpi_buffer_detach()            : dispatcher->buffer_detach()
 * ncmpi_bput_var<kind>_<type>()    : dispatcher->bput_var()
 *
 * ncmpi_get_varn_<type>()          : dispatcher->get_varn()
 * ncmpi_put_varn_<type>()          : dispatcher->put_varn()
 *
 * ncmpi_iget_varn_<type>()         : dispatcher->iget_varn()
 * ncmpi_iput_varn_<type>()         : dispatcher->iput_varn()
 * ncmpi_bput_varn_<type>()         : dispatcher->bput_varn()
 *
 * ncmpi_get_vard()                 : dispatcher->get_vard()
 * ncmpi_put_vard()                 : dispatcher->put_vard()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncbbio_driver.h>

int
ncbbio_def_var(void       *ncdp,
              const char *name,
              nc_type     xtype,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->def_var(ncbbp->ncp, name, xtype, ndims, dimids, varidp);
    if (err != NC_NOERR) return err;

    /* Update max_ndims */
    if (ndims > ncbbp->max_ndims){
        ncbbp->max_ndims = ndims;
    }

    return NC_NOERR;
}

int
ncbbio_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->inq_varid(ncbbp->ncp, name, varid);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_inq_var(void       *ncdp,
              int         varid,
              char       *name,
              nc_type    *xtypep,
              int        *ndimsp,
              int        *dimids,
              int        *nattsp,
              MPI_Offset *offsetp,
              int        *no_fillp,
              void       *fill_valuep)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->inq_var(ncbbp->ncp, varid, name, xtypep, ndimsp, dimids,
                               nattsp, offsetp, no_fillp, fill_valuep);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->rename_var(ncbbp->ncp, varid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_get_var(void             *ncdp,
              int               varid,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              void             *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               reqMode)
{
    int err, status = NC_NOERR;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Flush on read */
    if(ncbbp->inited){
        err = ncbbio_log_flush(ncbbp);
        if (status == NC_NOERR){
            status = err;
        }
    }

    err = ncbbp->ncmpio_driver->get_var(ncbbp->ncp, varid, start, count, stride, imap,
                               buf, bufcount, buftype, reqMode);
    if (status == NC_NOERR){
        status = err;
    }

    return status;
}

int
ncbbio_put_var(void             *ncdp,
              int               varid,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              const void       *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               reqMode)
{
    int err;
    void *cbuf=(void*)buf;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Resolve imap */
    if (imap != NULL || bufcount != -1) {
        /* pack buf to cbuf -------------------------------------------------*/
        /* If called from a true varm API or a flexible API, ncmpii_pack()
         * packs user buf into a contiguous cbuf (need to be freed later).
         * Otherwise, cbuf is simply set to buf. ncmpii_pack() also returns
         * etype (MPI primitive datatype in buftype), and nelems (number of
         * etypes in buftype * bufcount)
         */
        int ndims;
        MPI_Offset nelems;
        MPI_Datatype etype;

        err = ncbbp->ncmpio_driver->inq_var(ncbbp->ncp, varid, NULL, NULL, &ndims, NULL,
                                   NULL, NULL, NULL, NULL);
        if (err != NC_NOERR) return err;;

        if (buftype != MPI_DATATYPE_NULL){
            err = ncmpii_pack(ndims, count, imap, (void*)buf, bufcount, buftype,
                            &nelems, &etype, &cbuf);
            if (err != NC_NOERR) return err;
        }
        else{
            etype = MPI_DATATYPE_NULL;
            nelems = 0;
        }

        imap     = NULL;
        bufcount = (nelems == 0) ? 0 : -1;  /* make it a high-level API */
        buftype  = etype;                   /* an MPI primitive type */
    }

    /* Add log entry */
    err = ncbbio_log_put_var(ncbbp, varid, start, count, stride, cbuf, buftype, NULL);

    if (cbuf != buf) NCI_Free(cbuf);

    return err;
}

int
ncbbio_iget_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               void             *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->iget_var(ncbbp->ncp, varid, start, count, stride, imap,
                                buf, bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    // Translate from ncmpio id to ncbbio id
    *reqid = *reqid * 2 + 1;

    return NC_NOERR;
}

/*
 * Perform non-blocking put operation
 * We handle it the same way as blocking put
 * We allocate a request object that links to corresponding entries in the log
 * The id identifying the request object is then given to the user
 */
int
ncbbio_iput_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int i, err, id;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Initialize driver if not yet initialized 
     * Since PnetCDF allow nonblocking requests in define mode, we must initialize the driver if it hasn't been initialized
     */
    if (!ncbbp->inited){
        err = ncbbio_init(ncbbp);
        if (err != NC_NOERR){
            return err;
        }
    }

    // Create a new put request with id
    err = ncbbio_put_list_add(ncbbp, &id);
    if (err != NC_NOERR){
        return err;
    }

    // Translate putlist id to ncbbio id
    if (reqid != NULL){
        *reqid = id * 2;
    }

    /* We must link the request object to corresponding log entries
     * For varn operation, we will create multiple log entries, so it's a 1 to many mapping
     * Assuming the program runs under single thread, those entries are a continuous region within the metadata log
     * We record the first and last metadata entries associated with this request
     * This is done by checking number of lof entries before and after handling this operation
     * We also need to link those log entries to the request object so statues can be sset when log is flushed
     * Since one log entry can be associate with at most one request, it is at most 1 to 1
     */

    // Number of log entries before recording current operation to log
    ncbbp->putlist.reqs[id].entrystart = ncbbp->metaidx.nused;
    ncbbp->putlist.reqs[id].entrystart = ncbbp->metaidx.nused;

    err = ncbbio_put_var(ncdp, varid, start, count, stride, imap, buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR){
        ncbbio_put_list_remove(ncbbp, id);
        return err;
    }

    // Number of log entries after recording current operation to log
    ncbbp->putlist.reqs[id].entryend = ncbbp->metaidx.nused;

    /*
     * If new entry is created in the log, link those entries to the request
     * The entry may go directly to the ncmpio driver if it is too large
     * If there are no entry created, we mark this request as completed
     */
    if (ncbbp->putlist.reqs[id].entryend > ncbbp->putlist.reqs[id].entrystart) {
        for (i = ncbbp->putlist.reqs[id].entrystart; i < ncbbp->putlist.reqs[id].entryend; i++) {
            ncbbp->metaidx.entries[i].reqid = id;
        }
    }
    else{
        ncbbp->putlist.reqs[id].ready = 1;
        ncbbp->putlist.reqs[id].status = NC_NOERR;
    }

    return NC_NOERR;
}

int
ncbbio_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->buffer_attach(ncbbp->ncp, bufsize);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_buffer_detach(void *ncdp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->buffer_detach(ncbbp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_bput_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int err;

    /* bput same as iput in bb driver */
    err = ncbbio_iput_var(ncdp, varid, start, count, stride, imap, buf, bufcount, buftype, reqid, reqMode);

    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
int
ncbbio_get_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               void              *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int err, status = NC_NOERR;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Flush on read */
    if(ncbbp->inited){
        err = ncbbio_log_flush(ncbbp);
        if (status == NC_NOERR){
            status = err;
        }
    }

    err = ncbbp->ncmpio_driver->get_varn(ncbbp->ncp, varid, num, starts, counts, buf,
                                bufcount, buftype, reqMode);
    if (status == NC_NOERR){
        status = err;
    }

    return status;
}

/*
 * We treate varn operation as n of var operation
 */
int
ncbbio_put_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               const void        *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int i, err, status = NC_NOERR;
    MPI_Offset size;
    void *cbuf = (void*)buf;
    void *bufp;
    NC_bb *ncbbp = (NC_bb*)ncdp;
    MPI_Datatype ptype = buftype;

    /* It is illegal for starts to be NULL unless num is 0*/
    if (num > 0 && starts == NULL){
        DEBUG_RETURN_ERROR(NC_ENULLSTART)
    }

    /* Resolve flexible api so we can calculate size of each put_var */
    if (bufcount != -1){
        int isderived, iscontig_of_ptypes, elsize, position = 0;
        MPI_Offset bnelems;

        err = ncmpii_dtype_decode(buftype, &ptype, &elsize, &bnelems,
                                  &isderived, &iscontig_of_ptypes);
        if (err != NC_NOERR) return err;

        if (!iscontig_of_ptypes) { /* pack only if non-contiguous */
            bnelems *= elsize;
            if (bnelems != (int)bnelems) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

            cbuf = NCI_Malloc(bnelems);
            MPI_Pack(buf, (int)bufcount, buftype, cbuf, (int)bnelems,
                     &position, MPI_COMM_SELF);
        }
    }

    /* Decompose it into num put_vara calls */
    bufp = cbuf;
    for(i = 0; i < num; i++){
        err = ncbbio_log_put_var(ncbbp, varid, starts[i], counts[i], NULL, bufp, ptype, &size);
        if (status == NC_NOERR){
            status = err;
        }
        bufp = (void*)(((char*)bufp) + size);
    }

    if (cbuf != buf){
        NCI_Free(cbuf);
    }

    return status;
}

int
ncbbio_iget_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                void               *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->iget_varn(ncbbp->ncp, varid, num, starts, counts, buf,
                                 bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    // Translate form ncmpio id to ncbbio id
    *reqid = *reqid * 2 + 1;

    return NC_NOERR;
}

/*
 * Perform non-blocking put operation
 * We handle it the same way as blocking put
 * We allocate a request object that links to corresponding entries in the log
 * The id identifying the request object is then given to the user
 */
int
ncbbio_iput_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                const void         *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int i, err, id;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* It is illegal for starts to be NULL unless num is 0*/
    if (num > 0 && starts == NULL){
        DEBUG_RETURN_ERROR(NC_ENULLSTART)
    }

    /* Initialize driver if not yet initialized 
     * Since PnetCDF allow nonblocking requests in define mode, we must initialize the driver if it hasn't been initialized
     */
    if (!ncbbp->inited){
        err = ncbbio_init(ncbbp);
        if (err != NC_NOERR){
            return err;
        }
    }

    // Create a new put request with id
    err = ncbbio_put_list_add(ncbbp, &id);
    if (err != NC_NOERR){
        return err;
    }

    // Translate putlist id to ncbbio id
    if (reqid != NULL){
        *reqid = id * 2;
    }

    /* We must link the request object to corresponding log entries
     * For varn operation, we will create multiple log entries, so it's a 1 to many mapping
     * Assuming the program runs under single thread, those entries are a continuous region within the metadata log
     * We record the first and last metadata entries associated with this request
     * This is done by checking number of lof entries before and after handling this operation
     * We also need to link those log entries to the request object so statues can be sset when log is flushed
     * Since one log entry can be associate with at most one request, it is at most 1 to 1
     */

    // Number of log entries before recording current operation to log
    ncbbp->putlist.reqs[id].entrystart = ncbbp->metaidx.nused;

    // Handle the IO operation same as blocking one
    err = ncbbio_put_varn(ncdp, varid, num, starts, counts, buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR){
        ncbbio_put_list_remove(ncbbp, id);
        return err;
    }

    // Number of log entries after recording current operation to log
    ncbbp->putlist.reqs[id].entryend = ncbbp->metaidx.nused;

    /* If new entry is created in the log, link those entries to the request
     * The entry may go directly to the ncmpio driver if it is too large
     * If there are no entry created, we mark this request as completed
     */
    if (ncbbp->putlist.reqs[id].entryend > ncbbp->putlist.reqs[id].entrystart) {
        for (i = ncbbp->putlist.reqs[id].entrystart; i < ncbbp->putlist.reqs[id].entryend; i++) {
            ncbbp->metaidx.entries[i].reqid = id;
        }
    }
    else{
        ncbbp->putlist.reqs[id].ready = 1;
        ncbbp->putlist.reqs[id].status = NC_NOERR;
    }

    return NC_NOERR;
}

int
ncbbio_bput_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                const void         *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int err;

    /* bput same as iput in bb driver */
    err = ncbbio_iput_varn(ncdp, varid, num, starts, counts, buf,
                                 bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_get_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err, status = NC_NOERR;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* Flush on read */
    if(ncbbp->inited){
        err = ncbbio_log_flush(ncbbp);
        if (status == NC_NOERR){
            status = err;
        }
    }

    err = ncbbp->ncmpio_driver->get_vard(ncbbp->ncp, varid, filetype, buf, bufcount,
                                buftype, reqMode);
    if (status == NC_NOERR){
        status = err;
    }

    return status;
}

int
ncbbio_put_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    /* BB driver does not support vard */
    err = ncbbp->ncmpio_driver->put_vard(ncbbp->ncp, varid, filetype, buf, bufcount,
                                buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

