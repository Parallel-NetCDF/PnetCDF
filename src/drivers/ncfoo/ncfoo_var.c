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
#include <ncfoo_driver.h>

int
ncfoo_def_var(void       *ncdp,
              const char *name,
              nc_type     xtype,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->def_var(foo->ncp, name, xtype, ndims, dimids, varidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->inq_varid(foo->ncp, name, varid);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_inq_var(void       *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->inq_var(foo->ncp, varid, name, xtypep, ndimsp, dimids,
                               nattsp, offsetp, no_fillp, fill_valuep);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->rename_var(foo->ncp, varid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_get_var(void             *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->get_var(foo->ncp, varid, start, count, stride, imap,
                               buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_put_var(void             *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

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

        err = foo->driver->inq_var(foo->ncp, varid, NULL, NULL, &ndims, NULL,
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

    status = foo->driver->put_var(foo->ncp, varid, start, count, stride, imap,
                                  cbuf, bufcount, buftype, reqMode);
    if (cbuf != buf) NCI_Free(cbuf);

    return (err == NC_NOERR) ? status : err; /* first error encountered */
}

int
ncfoo_iget_var(void             *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->iget_var(foo->ncp, varid, start, count, stride, imap,
                                buf, bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_iput_var(void             *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->iput_var(foo->ncp, varid, start, count, stride, imap,
                                buf, bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->buffer_attach(foo->ncp, bufsize);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_buffer_detach(void *ncdp)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->buffer_detach(foo->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_bput_var(void             *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->bput_var(foo->ncp, varid, start, count, stride, imap,
                                buf, bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
int
ncfoo_get_varn(void              *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->get_varn(foo->ncp, varid, num, starts, counts, buf,
                                bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_put_varn(void              *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->put_varn(foo->ncp, varid, num, starts, counts, buf,
                                bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_iget_varn(void               *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->iget_varn(foo->ncp, varid, num, starts, counts, buf,
                                 bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_iput_varn(void               *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->iput_varn(foo->ncp, varid, num, starts, counts, buf,
                                 bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_bput_varn(void               *ncdp,
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
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->bput_varn(foo->ncp, varid, num, starts, counts, buf,
                                 bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_get_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->get_vard(foo->ncp, varid, filetype, buf, bufcount,
                                buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_put_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->put_vard(foo->ncp, varid, filetype, buf, bufcount,
                                buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

