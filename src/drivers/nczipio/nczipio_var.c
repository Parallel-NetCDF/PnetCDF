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
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>
#include "nczipio_internal.h"

int
nczipio_def_var(void       *ncdp,
              const char *name,
              nc_type     xtype,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int i, err;
    NC_zip *nczipp = (NC_zip*)ncdp;
    NC_zip_var var;

    var.ndim = ndims;
    var.chunkdim = NULL;
    var.data_offs = NULL;
    var.chunk_owner = NULL;
    var.xtype = xtype;
    var.esize = NC_Type_size(xtype);
    var.etype = ncmpii_nc2mpitype(xtype);

    if (ndims > 3 || ndims < 1) { // Does not support higher dimensional vars
        var.varkind = NC_ZIP_VAR_RAW;
        var.dimsize = NULL;

        err = nczipp->driver->def_var(nczipp->ncp, name, xtype, ndims, dimids, &var.varid);  // We use it to save the id of data variable
        if (err != NC_NOERR) return err;
        
        err = nczipp->driver->put_att(nczipp->ncp, var.varid, "_varkind", NC_INT, 1, &(var.varkind), MPI_INT);   // Comressed var?
        if (err != NC_NOERR) return err;
    }
    else{
        err = nczipp->driver->def_var(nczipp->ncp, name, NC_INT, 0, NULL, &var.varid);  // We use it to save the id of data variable
        if (err != NC_NOERR) return err;
        
        var.varkind = NC_ZIP_VAR_COMPRESSED;
        var.dimids = (int*)NCI_Malloc(sizeof(int) * ndims);
        memcpy(var.dimids, dimids, sizeof(int) * ndims);
        var.dimsize = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * ndims);
        for(i = 0; i < ndims; i++){
            nczipp->driver->inq_dim(nczipp->ncp, dimids[i], NULL, var.dimsize + i);
        }
        if (var.dimids[0] == nczipp->recdim){
            var.isrec = 1;
        }

        err = nczipp->driver->put_att(nczipp->ncp, var.varid, "_ndim", NC_INT, 1, &ndims, MPI_INT); // Original dimensions
        if (err != NC_NOERR) return err;
        err = nczipp->driver->put_att(nczipp->ncp, var.varid, "_dimids", NC_INT, ndims, dimids, MPI_INT);   // Dimensiona IDs
        if (err != NC_NOERR) return err;
        err = nczipp->driver->put_att(nczipp->ncp, var.varid, "_datatype", NC_INT, 1, &xtype, MPI_INT); // Original datatype
        if (err != NC_NOERR) return err;
        err = nczipp->driver->put_att(nczipp->ncp, var.varid, "_varkind", NC_INT, 1, &(var.varkind), MPI_INT);   // Comressed var?
        if (err != NC_NOERR) return err;
    }

    err = nczipioi_var_list_add(&(nczipp->vars), var);
    if (err != NC_NOERR) return err;

    *varidp = nczipp->vars.cnt - 1;

    return NC_NOERR;
}

int
nczipio_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
    int i, vid, err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->inq_varid(nczipp->ncp, name, &vid);
    if (err != NC_NOERR) return err;

    if (varid != NULL){
        for(i = 0; i < nczipp->vars.cnt; i++){
            if (nczipp->vars.data[i].varid == vid){
                *varid = i;
                break;
            }
        }
        if (i >= nczipp->vars.cnt){
            DEBUG_RETURN_ERROR(NC_ENOTVAR)
        }
    }

    return NC_NOERR;
}

int
nczipio_inq_var(void       *ncdp,
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
    NC_zip *nczipp = (NC_zip*)ncdp;
    NC_zip_var *varp;

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }

    varp = nczipp->vars.data + varid;

    err = nczipp->driver->inq_var(nczipp->ncp, varp->varid, name, xtypep, NULL, NULL,
                               nattsp, offsetp, no_fillp, fill_valuep);
    if (err != NC_NOERR) return err;

    if (ndimsp != NULL){
        *ndimsp = varp->ndim;
    }

    if (dimids != NULL){
        memcpy(dimids, varp->dimids, sizeof(int) * varp->ndim);
    }

    return NC_NOERR;
}

int
nczipio_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;
    NC_zip_var *varp;

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    err = nczipp->driver->rename_var(nczipp->ncp, varp->varid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_get_var(void             *ncdp,
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
    int err=NC_NOERR, status;
    void *cbuf=(void*)buf;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    if (varp->varkind == NC_ZIP_VAR_RAW){
        return nczipp->driver->get_var(nczipp->ncp, varp->varid, start, count, stride, imap, buf, bufcount, buftype, reqMode);
    }

    if (varp->isrec && (varp->dimsize[0] < nczipp->recsize) && (start[0] + count[0] >= varp->dimsize[0])){
        nczipioi_var_resize(nczipp, varp);
    }

    // Collective buffer
    switch (nczipp->comm_unit){
        case NC_ZIP_COMM_CHUNK:
            status = nczipioi_get_var_cb_chunk(nczipp, varp, start, count, stride, buf);
            break;
        case NC_ZIP_COMM_PROC:
            status = nczipioi_get_var_cb_proc(nczipp, varp, start, count, stride, buf);
            break;
    }

    return (err == NC_NOERR) ? status : err; /* first error encountered */
}

int
nczipio_put_var(void             *ncdp,
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
    int err=NC_NOERR, status;
    void *cbuf=(void*)buf;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (reqMode == NC_REQ_INDEP){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;
    
    if (varp->varkind == NC_ZIP_VAR_RAW){
        return nczipp->driver->put_var(nczipp->ncp, varp->varid, start, count, stride, imap, buf, bufcount, buftype, reqMode);
    }

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

        err = nczipp->driver->inq_var(nczipp->ncp, varid, NULL, NULL, &ndims, NULL,
                                   NULL, NULL, NULL, NULL);
        if (err != NC_NOERR) goto err_check;

        err = ncmpii_pack(ndims, count, imap, (void*)buf, bufcount, buftype,
                          &nelems, &etype, &cbuf);
        if (err != NC_NOERR) goto err_check;

        imap     = NULL;
        bufcount = (nelems == 0) ? 0 : -1;  /* make it a high-level API */
        buftype  = etype;                   /* an MPI primitive type */
    }

err_check:
    if (err != NC_NOERR) {
        if (reqMode & NC_REQ_INDEP) return err;
        reqMode |= NC_REQ_ZERO; /* participate collective call */
    }

    status = nczipioi_put_var(nczipp, varp, start, count, stride, cbuf);
    if (cbuf != buf) NCI_Free(cbuf);

    return (err == NC_NOERR) ? status : err; /* first error encountered */
}


int
nczipio_iget_var(void             *ncdp,
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
    void *cbuf=(void*)buf;
    void *xbuf=(void*)buf;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (reqMode == NC_REQ_INDEP){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    if (varp->varkind == NC_ZIP_VAR_RAW){
        err = nczipp->driver->iget_var(nczipp->ncp, varp->varid, start, count, stride, imap, buf, bufcount, buftype, reqid, reqMode);
        if (err != NC_NOERR){
            return err;
        }
        *reqid = *reqid * 2 + 1;
        return NC_NOERR;
    }

    nczipioi_iget_var(nczipp, varid, start, count, stride, imap, buf, bufcount, buftype, reqid);
    *reqid *= 2;

    return NC_NOERR;
}

int
nczipio_iput_var(void             *ncdp,
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
    int err=NC_NOERR, status;
    void *cbuf=(void*)buf;
    void *xbuf=(void*)buf;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (reqMode == NC_REQ_INDEP){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    if (varp->varkind == NC_ZIP_VAR_RAW){
        err = nczipp->driver->iput_var(nczipp->ncp, varp->varid, start, count, stride, imap, buf, bufcount, buftype, reqid, reqMode);
        if (err != NC_NOERR){
            return err;
        }
        *reqid = *reqid * 2 + 1;
        return NC_NOERR;
    }

    if (varp->isrec){
        if (nczipp->recsize < start[0] + count[0]){
            nczipp->recsize = start[0] + count[0];
        }
    }

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

        err = nczipp->driver->inq_var(nczipp->ncp, varid, NULL, NULL, &ndims, NULL,
                                   NULL, NULL, NULL, NULL);
        if (err != NC_NOERR) goto err_check;

        err = ncmpii_pack(ndims, count, imap, (void*)buf, bufcount, buftype,
                          &nelems, &etype, &cbuf);
        if (err != NC_NOERR) goto err_check;

        imap     = NULL;
        bufcount = (nelems == 0) ? 0 : -1;  /* make it a high-level API */
        buftype  = etype;                   /* an MPI primitive type */
    }

err_check:
    if (err != NC_NOERR) {
        if (reqMode & NC_REQ_INDEP) return err;
        reqMode |= NC_REQ_ZERO; /* participate collective call */
    }

    xbuf = cbuf;
    status = nczipioi_iput_var(nczipp, varid, start, count, stride, xbuf, buf, reqid);
    (*reqid) *= 2;

    //if (cbuf != buf) NCI_Free(cbuf);

    return (err == NC_NOERR) ? status : err; /* first error encountered */

    return NC_NOERR;
}

int
nczipio_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    err = nczipp->driver->buffer_attach(nczipp->ncp, bufsize);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_buffer_detach(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    err = nczipp->driver->buffer_detach(nczipp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_bput_var(void             *ncdp,
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
    NC_zip *nczipp = (NC_zip*)ncdp;

    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    err = nczipp->driver->bput_var(nczipp->ncp, varid, start, count, stride, imap,
                                buf, bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
int
nczipio_get_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               void              *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int err;
    int i;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (reqMode == NC_REQ_INDEP){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    if (varp->varkind == NC_ZIP_VAR_RAW){
        return nczipp->driver->get_varn(nczipp->ncp, varp->varid, num, starts, counts, buf, bufcount, buftype, reqMode);
    }

    if (varp->isrec && (varp->dimsize[0] < nczipp->recsize)){
        for(i = 0; i < num; i++){
            if (starts[i][0] + counts[i][0] >= varp->dimsize[0]){
                nczipioi_var_resize(nczipp, varp);
                break;
            }
        }
    }

    err = nczipioi_get_varn(nczipp, varp, num, starts, counts, buf);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_put_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               const void        *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int err;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (reqMode == NC_REQ_INDEP){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    if (varp->varkind == NC_ZIP_VAR_RAW){
        return nczipp->driver->put_varn(nczipp->ncp, varp->varid, num, starts, counts, buf, bufcount, buftype, reqMode);
    }

    err = nczipioi_put_varn(nczipp, varp, num, starts, counts, buf);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_iget_varn(void               *ncdp,
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
    void *cbuf=(void*)buf;
    void *xbuf=(void*)buf;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (reqMode == NC_REQ_INDEP){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    if (varp->varkind == NC_ZIP_VAR_RAW){
        err = nczipp->driver->iget_varn(nczipp->ncp, varp->varid, num, starts, counts, buf, bufcount, buftype, reqid, reqMode);
        if (err != NC_NOERR){
            return err;
        }
        *reqid = *reqid * 2 + 1;
        return NC_NOERR;
    }

    xbuf = cbuf;
    nczipioi_iget_varn(nczipp, varid, num, starts, counts, buf, bufcount, buftype, reqid);
    *reqid *= 2;

    return NC_NOERR;
}

int
nczipio_iput_varn(void               *ncdp,
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
    int i;
    void *cbuf=(void*)buf;
    void *xbuf=(void*)buf;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (reqMode == NC_REQ_INDEP){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    if (varp->isrec){
        for(i = 0; i < num; i++){
            if (nczipp->recsize < starts[i][0] + counts[i][0]){
                nczipp->recsize = starts[i][0] + counts[i][0];
            }
        }
    }

    if (varp->varkind == NC_ZIP_VAR_RAW){
        err = nczipp->driver->iput_varn(nczipp->ncp, varp->varid, num, starts, counts, buf, bufcount, buftype, reqid, reqMode);
        if (err != NC_NOERR){
            return err;
        }
        *reqid = *reqid * 2 + 1;
        return NC_NOERR;
    }

    xbuf = cbuf;
    nczipioi_iput_varn(nczipp, varid, num, starts, counts, xbuf, buf, reqid);
    *reqid *= 2;

    return NC_NOERR;
}

int
nczipio_bput_varn(void               *ncdp,
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
    NC_zip *nczipp = (NC_zip*)ncdp;

    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    err = nczipp->driver->bput_varn(nczipp->ncp, varid, num, starts, counts, buf,
                                 bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_get_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    err = nczipp->driver->get_vard(nczipp->ncp, varid, filetype, buf, bufcount,
                                buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_put_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    err = nczipp->driver->put_vard(nczipp->ncp, varid, filetype, buf, bufcount,
                                buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

