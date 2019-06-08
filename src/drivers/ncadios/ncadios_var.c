/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
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
    }

    return NC_NOERR;
}

int
ncadios_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    NC_ad *ncadp = (NC_ad*)ncdp;
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

    /* Release var info */
    adios_free_varinfo (v);

    /* Handle the request */
    err = ncadiosi_handle_get_req(ncadp, &r);

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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    return ncadiosi_iget_var(ncadp, varid, start, count, stride, imap, buf,
                                bufcount, buftype, reqid);
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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_buffer_detach(void *ncdp)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    /* No support for varn at this time
     * ADIOS interface make varn difficult to be implemented efficiently
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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    /* No varn support */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    /* ADIOS has no vard interface */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

