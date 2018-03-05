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
#include <ncdwio_driver.h>

int
ncdwio_inq_attname(void *ncdp,
                  int   varid,
                  int   attid,
                  char *name)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->inq_attname(ncdwp->ncp, varid, attid, name);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_inq_attid(void       *ncdp,
                int         varid,
                const char *name,
                int        *attidp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->inq_attid(ncdwp->ncp, varid, name, attidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_inq_att(void       *ncdp,
              int         varid,
              const char *name,
              nc_type    *datatypep,
              MPI_Offset *lenp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->inq_att(ncdwp->ncp, varid, name, datatypep, lenp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_rename_att(void       *ncdp,
                 int         varid,
                 const char *name,
                 const char *newname)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->rename_att(ncdwp->ncp, varid, name, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}


int
ncdwio_copy_att(void       *ncdp_in,
               int         varid_in,
               const char *name,
               void       *ncdp_out,
               int         varid_out)
{
    int err;
    NC_dw *ncdwp_in  = (NC_dw*)ncdp_in;
    NC_dw *ncdwp_out = (NC_dw*)ncdp_out;
    
    err = ncdwp_in->ncmpio_driver->copy_att(ncdwp_in->ncp,  varid_in, name,
                                   ncdwp_out->ncp, varid_out);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_del_att(void       *ncdp,
              int         varid,
              const char *name)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->del_att(ncdwp->ncp, varid, name);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_get_att(void         *ncdp,
              int           varid,
              const char   *name,
              void         *buf,
              MPI_Datatype  itype)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->get_att(ncdwp->ncp, varid, name, buf, itype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_put_att(void         *ncdp,
              int           varid,
              const char   *name,
              nc_type       xtype,
              MPI_Offset    nelems,
              const void   *buf,
              MPI_Datatype  itype)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->put_att(ncdwp->ncp, varid, name, xtype, nelems, buf,
                               itype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
