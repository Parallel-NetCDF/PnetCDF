/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_inq_attname() : dispatcher->inq_attname()
 * ncmpi_inq_attid()   : dispatcher->inq_attid()
 * ncmpi_inq_att()     : dispatcher->inq_att()
 * ncmpi_rename_att()  : dispatcher->inq_rename_att()
 * ncmpi_copy_att()    : dispatcher->inq_copy_att()
 * ncmpi_del_att()     : dispatcher->inq_del_att()
 * ncmpi_get_att()     : dispatcher->inq_get_att()
 * ncmpi_put_att()     : dispatcher->inq_put_arr()
 *
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncchkio_driver.h>
#include "ncchkio_internal.h"

int
ncchkio_inq_attname(void *ncdp,
                  int   varid,
                  int   attid,
                  char *name)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->inq_attname(ncchkp->ncp, varid, attid, name);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncchkio_inq_attid(void       *ncdp,
                int         varid,
                const char *name,
                int        *attidp)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->inq_attid(ncchkp->ncp, varid, name, attidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncchkio_inq_att(void       *ncdp,
              int         varid,
              const char *name,
              nc_type    *datatypep,
              MPI_Offset *lenp)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->inq_att(ncchkp->ncp, varid, name, datatypep, lenp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncchkio_rename_att(void       *ncdp,
                 int         varid,
                 const char *name,
                 const char *newname)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->rename_att(ncchkp->ncp, varid, name, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}


int
ncchkio_copy_att(void       *ncdp_in,
               int         varid_in,
               const char *name,
               void       *ncdp_out,
               int         varid_out)
{
    int err;
    NC_chk *foo_in  = (NC_chk*)ncdp_in;
    NC_chk *foo_out = (NC_chk*)ncdp_out;

    err = foo_in->driver->copy_att(foo_in->ncp,  varid_in, name,
                                   foo_out->ncp, varid_out);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncchkio_del_att(void       *ncdp,
              int         varid,
              const char *name)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->del_att(ncchkp->ncp, varid, name);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncchkio_get_att(void         *ncdp,
              int           varid,
              const char   *name,
              void         *buf,
              MPI_Datatype  itype)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->get_att(ncchkp->ncp, varid, name, buf, itype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncchkio_put_att(void         *ncdp,
              int           varid,
              const char   *name,
              nc_type       xtype,
              MPI_Offset    nelems,
              const void   *buf,
              MPI_Datatype  itype)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->put_att(ncchkp->ncp, varid, name, xtype, nelems, buf,
                               itype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
