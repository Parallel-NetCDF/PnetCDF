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
#include <ncdwio_driver.h>

int
ncdwio_def_var(void       *ncdp,
              const char *name,
              nc_type     xtype,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->def_var(ncdwp->ncp, name, xtype, ndims, dimids, varidp);
    if (err != NC_NOERR) return err;

    /* Update max_ndims */
    if (ndims > ncdwp->max_ndims){
        ncdwp->max_ndims = ndims;
    }

    return NC_NOERR;
}

int
ncdwio_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->inq_varid(ncdwp->ncp, name, varid);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_inq_var(void       *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->inq_var(ncdwp->ncp, varid, name, xtypep, ndimsp, dimids,
                               nattsp, offsetp, no_fillp, fill_valuep);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->rename_var(ncdwp->ncp, varid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_get_var(void             *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;

    /* Flush on read */
    if(ncdwp->inited){
        err = ncdwio_log_flush(ncdwp);
        if (status == NC_NOERR){
            status = err;
        }
    }

    err = ncdwp->ncmpio_driver->get_var(ncdwp->ncp, varid, start, count, stride, imap,
                               buf, bufcount, buftype, reqMode);
    if (status == NC_NOERR){
        status = err;
    }

    return status;
}

int
ncdwio_put_var(void             *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;

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

        err = ncdwp->ncmpio_driver->inq_var(ncdwp->ncp, varid, NULL, NULL, &ndims, NULL,
                                   NULL, NULL, NULL, NULL);
        if (err != NC_NOERR) return err;;

        err = ncmpii_pack(ndims, count, imap, (void*)buf, bufcount, buftype,
                          &nelems, &etype, &cbuf);
        if (err != NC_NOERR) return err;

        imap     = NULL;
        bufcount = (nelems == 0) ? 0 : -1;  /* make it a high-level API */
        buftype  = etype;                   /* an MPI primitive type */
    }

    /* Add log entry */
    err = ncdwio_log_put_var(ncdwp, varid, start, count, stride, cbuf, buftype, NULL);

    if (cbuf != buf) NCI_Free(cbuf);

    return err;
}

int
ncdwio_iget_var(void             *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->iget_var(ncdwp->ncp, varid, start, count, stride, imap,
                                buf, bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    // Translate from ncmpio id to ncdwio id
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
ncdwio_iput_var(void             *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;

    // Create a new put request with id
    err = ncdwio_put_list_add(ncdwp, &id);
    if (err != NC_NOERR){
        return err;
    }

    // Translate putlist id to ncdwio id
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
    ncdwp->putlist.reqs[id].entrystart = ncdwp->metaidx.nused;
    ncdwp->putlist.reqs[id].entrystart = ncdwp->metaidx.nused;

    err = ncdwio_put_var(ncdp, varid, start, count, stride, imap, buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR){
        ncdwio_put_list_remove(ncdwp, id);
        return err;
    }

    // Number of log entries after recording current operation to log
    ncdwp->putlist.reqs[id].entryend = ncdwp->metaidx.nused;

    /*
     * If new entry is created in the log, link those entries to the request
     * The entry may go directly to the ncmpio driver if it is too large
     * If there are no entry created, we mark this request as completed
     */
    if (ncdwp->putlist.reqs[id].entryend > ncdwp->putlist.reqs[id].entrystart) {
        for (i = ncdwp->putlist.reqs[id].entrystart; i < ncdwp->putlist.reqs[id].entryend; i++) {
            ncdwp->metaidx.entries[i].reqid = id;
        }
    }
    else{
        ncdwp->putlist.reqs[id].ready = 1;
        ncdwp->putlist.reqs[id].status = NC_NOERR;
    }

    return NC_NOERR;
}

int
ncdwio_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->buffer_attach(ncdwp->ncp, bufsize);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_buffer_detach(void *ncdp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->buffer_detach(ncdwp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_bput_var(void             *ncdp,
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
    err = ncdwio_iput_var(ncdp, varid, start, count, stride, imap, buf, bufcount, buftype, reqid, reqMode);

    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
int
ncdwio_get_varn(void              *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;

    /* Flush on read */
    if(ncdwp->inited){
        err = ncdwio_log_flush(ncdwp);
        if (status == NC_NOERR){
            status = err;
        }
    }

    err = ncdwp->ncmpio_driver->get_varn(ncdwp->ncp, varid, num, starts, counts, buf,
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
ncdwio_put_varn(void              *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;
    MPI_Datatype ptype = buftype;

    /* It is illegal for starts to be NULL unless num is 0*/
    if (num > 0 && starts == NULL){
        DEBUG_RETURN_ERROR(NC_ENULLSTART)
    }

    /* Resolve flexible api so we can calculate size of each put_var */
    if (bufcount != -1){
        int isderived, iscontig_of_ptypes;
        int elsize, position;
        MPI_Offset bnelems = 0;

        err = ncmpii_dtype_decode(buftype, &ptype, &elsize, &bnelems, &isderived, &iscontig_of_ptypes);
        if (err != NC_NOERR){
            return err;
        }

        cbuf = NCI_Malloc(elsize * bnelems);
        MPI_Pack(buf, (int)bufcount, buftype, cbuf, (int)(elsize * bnelems), &position, MPI_COMM_SELF);
    }

    /* Decompose it into num put_vara calls */
    bufp = cbuf;
    for(i = 0; i < num; i++){
        err = ncdwio_log_put_var(ncdwp, varid, starts[i], counts[i], NULL, bufp, ptype, &size);
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
ncdwio_iget_varn(void               *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;

    err = ncdwp->ncmpio_driver->iget_varn(ncdwp->ncp, varid, num, starts, counts, buf,
                                 bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    // Translate form ncmpio id to ncdwio id
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
ncdwio_iput_varn(void               *ncdp,
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
    NC_dw *ncdwp = (NC_dw*)ncdp;

    /* It is illegal for starts to be NULL unless num is 0*/
    if (num > 0 && starts == NULL){
        DEBUG_RETURN_ERROR(NC_ENULLSTART)
    }

    // Create a new put request with id
    err = ncdwio_put_list_add(ncdwp, &id);
    if (err != NC_NOERR){
        return err;
    }

    // Translate putlist id to ncdwio id
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
    ncdwp->putlist.reqs[id].entrystart = ncdwp->metaidx.nused;

    // Handle the IO operation same as blocking one
    err = ncdwio_put_varn(ncdp, varid, num, starts, counts, buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR){
        ncdwio_put_list_remove(ncdwp, id);
        return err;
    }

    // Number of log entries after recording current operation to log
    ncdwp->putlist.reqs[id].entryend = ncdwp->metaidx.nused;

    /* If new entry is created in the log, link those entries to the request
     * The entry may go directly to the ncmpio driver if it is too large
     * If there are no entry created, we mark this request as completed
     */
    if (ncdwp->putlist.reqs[id].entryend > ncdwp->putlist.reqs[id].entrystart) {
        for (i = ncdwp->putlist.reqs[id].entrystart; i < ncdwp->putlist.reqs[id].entryend; i++) {
            ncdwp->metaidx.entries[i].reqid = id;
        }
    }
    else{
        ncdwp->putlist.reqs[id].ready = 1;
        ncdwp->putlist.reqs[id].status = NC_NOERR;
    }

    return NC_NOERR;
}

int
ncdwio_bput_varn(void               *ncdp,
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
    err = ncdwio_iput_varn(ncdp, varid, num, starts, counts, buf,
                                 bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_get_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err, status = NC_NOERR;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    /* Flush on read */
    if(ncdwp->inited){
        err = ncdwio_log_flush(ncdwp);
        if (status == NC_NOERR){
            status = err;
        }
    }

    err = ncdwp->ncmpio_driver->get_vard(ncdwp->ncp, varid, filetype, buf, bufcount,
                                buftype, reqMode);
    if (status == NC_NOERR){
        status = err;
    }

    return status;
}

int
ncdwio_put_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;

    /* BB driver does not support vard */
    err = ncdwp->ncmpio_driver->put_vard(ncdwp->ncp, varid, filetype, buf, bufcount,
                                buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

