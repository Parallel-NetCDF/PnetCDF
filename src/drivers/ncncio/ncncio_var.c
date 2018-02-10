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
#include <ncncio_driver.h>

int
ncncio_def_var(void       *ncdp,
              const char *name,
              nc_type     xtype,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

int
ncncio_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_inq_varid */
    err = nc_inq_varid(nc4p->ncid, name, varid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
ncncio_inq_var(void       *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_inq_var_all */
    err = nc_inq_var_all(nc4p->ncid, varid, name, xtypep, ndimsp, dimids, nattsp, 
                        NULL, NULL, NULL, NULL, NULL, NULL,
                        no_fillp, fill_valuep, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
ncncio_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

/* 
ncncio_get_var is implemented iin ncmpio_get_put.m4
*/

int
ncncio_put_var(void             *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)
}

int
ncncio_iget_var(void             *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* TODO: Support nonblocking IO */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

int
ncncio_iput_var(void             *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

int
ncncio_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* We don't use buffer */
    /* Todo, count buffer size to report in inq */

    return NC_NOERR;
}

int
ncncio_buffer_detach(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* We don't use buffer */
    /* Todo, count buffer size to report in inq */

    return NC_NOERR;
}

int
ncncio_bput_var(void             *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}
int
ncncio_get_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               void              *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int i, err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call ncncio_get_var for N times */
    for(i = 0; i < num; i++){
        err = ncncio_get_var(ncdp, varid, starts[i], counts[i], NULL, NULL, buf, bufcount, buftype);
        if (err != NC_NOERR){
            return err;
        }
    }

    return NC_NOERR;
}

int
ncncio_put_varn(void              *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

int
ncncio_iget_varn(void               *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call ncncio_get_varn */
    err = ncncio_get_varn(ncdp, varid, num, starts, counts, buf, bufcount, buftype,Reqmode);
    if (err != NC_NOERR) return err;

    /* TODO: Issue dummy id */
    reqid = NC_REQ_NULL;

    return NC_NOERR;
}

int
ncncio_iput_varn(void               *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

int
ncncio_bput_varn(void               *ncdp,
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
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

int
ncncio_get_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* vard not supported in NetCDF */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

int
ncncio_put_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOT_SUPPORTED)

    return NC_NOERR;
}

