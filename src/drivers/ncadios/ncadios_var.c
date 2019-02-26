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
#include <ncadios_driver.h>
#include <ncadios_internal.h>
#include <string.h>

int
ncadios_def_var(void       *ncdp,
              const char *name,
              nc_type     xtype,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    return ncadiosi_inq_varid(ncadp, (char*)name, varid);
}

int
ncadios_inq_var(void       *ncdp,
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
    int i, j;
    NC_ad *ncadp = (NC_ad*)ncdp;
    NC_ad_var var;

    var = ncadp->vars.data[varid];
    
    if (xtypep != NULL){
        *xtypep = var.type;
    }

    if (ndimsp != NULL){
        *ndimsp = var.ndim;
    }

    if (dimids != NULL){
        memcpy(dimids, var.dimids, var.ndim * SIZEOF_INT);
    }

    if (nattsp != NULL){
        *nattsp = var.atts.cnt;
    }

    if (name != NULL){
        strcpy(name, var.name);
    }

    /* Not supported by adios, set to 0 */
    if (offsetp != NULL){
        *offsetp = 0;
    }
    if (no_fillp != NULL){
        *no_fillp = 0;
    }
    if (fill_valuep != NULL){
        //*fill_valuep = 0;
    }

    return NC_NOERR;
}

int
ncadios_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_get_var(void             *ncdp,
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
    int err;
    int i;
    NC_ad *ncadp = (NC_ad*)ncdp;
    ADIOS_VARINFO * v;
    ADIOS_SELECTION *sel;
    MPI_Datatype    vtype;
    size_t esize, ecnt = 0;
    int cesize;
    void *xbuf, *cbuf;
    uint64_t *points;
    int sstart, scount, sstride;

    v = adios_inq_var(ncadp->fp, ncadp->vars.data[varid].name);
    if (v == NULL){
        err = ncmpii_error_adios2nc(adios_errno, "get_var");
        DEBUG_RETURN_ERROR(err);
    }

    // Determine if it is record variable
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

    // Calculate number of elements in single record
    ecnt = 1;
    for(i = 0; i < v->ndim; i++){
        ecnt *= (size_t)count[i];
    }

    // If user buffer is contiguous
    if (imap == NULL){
        cbuf = buf;
    }
    else{
        MPI_Type_size(buftype, &cesize);
        cbuf = NCI_Malloc((size_t)cesize * ecnt * scount);
    }

    // PnetCDF allows accessing in different type
    // Check if we need to convert
    vtype = ncadios_to_mpi_type(v->type);
    if (vtype == buftype){
        xbuf = cbuf;
    }
    else{
        esize = (size_t)adios_type_size(v->type, NULL);
        xbuf = NCI_Malloc(esize * ecnt * scount);
    }

    // If stride is not used, we can use bounding box selection
    // Otherwise, we need to specify every points
    if (stride == NULL){
        sel = adios_selection_boundingbox (v->ndim, (uint64_t*)start, (uint64_t*)count);
    }
    else{
        uint64_t *p, *cur;

        points = (uint64_t*)NCI_Malloc(sizeof(uint64_t) * ecnt * v->ndim);
        p = (uint64_t*)NCI_Malloc(sizeof(uint64_t) * v->ndim);
        cur = points;

        memset(p, 0, sizeof(uint64_t) * v->ndim);

        // Iterate through every cells accessed
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
            }
        }

        sel = adios_selection_points(v->ndim, (uint64_t)ecnt, points);
    }
    if (sel == NULL){
        err = ncmpii_error_adios2nc(adios_errno, "select");
        DEBUG_RETURN_ERROR(err);
    }
    if (sstride > 1){
        for(i = 0; i < scount; i++){
            err = adios_schedule_read_byid (ncadp->fp, sel, v->varid, sstart + i * sstride, 1, (void*)(((char*)xbuf) + i * esize * ecnt));
            if (err != 0){
                err = ncmpii_error_adios2nc(adios_errno, "Open");
                DEBUG_RETURN_ERROR(err);
            }
        }
    }
    else{
        err = adios_schedule_read_byid (ncadp->fp, sel, v->varid, sstart, scount, xbuf);
        if (err != 0){
            err = ncmpii_error_adios2nc(adios_errno, "Open");
            DEBUG_RETURN_ERROR(err);
        }
    }
    err = adios_perform_reads (ncadp->fp, 1);
    if (err != 0){
        err = ncmpii_error_adios2nc(adios_errno, "Open");
        DEBUG_RETURN_ERROR(err);
    }

    if (stride != NULL){
        NCI_Free(points);
    }

    if (vtype != buftype){
        err = ncadiosiconvert(xbuf, cbuf, vtype, buftype, (int)ecnt * scount);
        if (err != NC_NOERR){
            return err;
        }
        NCI_Free(xbuf);
    }

    if (imap != NULL){
        int position;
        MPI_Datatype imaptype;

        if (scount > 1){
            count--;
        }

        err = ncmpii_create_imaptype(v->ndim, count, imap, buftype, &imaptype);
        if (err != NC_NOERR) {
            return err;
        }
        position = 0;
        MPI_Unpack(cbuf, cesize * (int)ecnt * scount, &position, buf, 1, imaptype, MPI_COMM_SELF);
        MPI_Type_free(&imaptype);

        NCI_Free(cbuf);
    }

    adios_free_varinfo (v);

    MPI_Type_size(vtype, &cesize);
    ncadp->getsize += cesize * ecnt;

    return NC_NOERR;
}

int
ncadios_put_var(void             *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return (err == NC_NOERR) ? status : err; /* first error encountered */
}

int
ncadios_iget_var(void             *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    return ncadiosi_iget_var(ncadp, varid, start, count, stride, imap, buf, bufcount, buftype, reqid);
} 

int
ncadios_iput_var(void             *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* TODO: Nonblocking support */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_buffer_detach(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* TODO: Nonblocking support */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_bput_var(void             *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}
int
ncadios_get_varn(void              *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* No support for varn at this time
     * It make varn difficult to be implemented efficiently 
     */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

}

int
ncadios_put_varn(void              *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_iget_varn(void               *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* TODO: nonblocking support */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_iput_varn(void               *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_bput_varn(void               *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_get_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* ADIOS has not vard interface */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_put_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

