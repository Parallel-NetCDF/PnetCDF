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

#include <ncbbio_driver.h>

int
ncbbio_inq_attname(void *ncdp,
                  int   varid,
                  int   attid,
                  char *name)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->inq_attname(ncbbp->ncp, varid, attid, name);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_inq_attid(void       *ncdp,
                int         varid,
                const char *name,
                int        *attidp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->inq_attid(ncbbp->ncp, varid, name, attidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_inq_att(void       *ncdp,
              int         varid,
              const char *name,
              nc_type    *datatypep,
              MPI_Offset *lenp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->inq_att(ncbbp->ncp, varid, name, datatypep, lenp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_rename_att(void       *ncdp,
                 int         varid,
                 const char *name,
                 const char *newname)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->rename_att(ncbbp->ncp, varid, name, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}


int
ncbbio_copy_att(void       *ncdp_in,
               int         varid_in,
               const char *name,
               void       *ncdp_out,
               int         varid_out)
{
    int err;
    NC_bb *ncbbp_in  = (NC_bb*)ncdp_in;
    NC_bb *ncbbp_out = (NC_bb*)ncdp_out;

    err = ncbbp_in->ncmpio_driver->copy_att(ncbbp_in->ncp,  varid_in, name,
                                   ncbbp_out->ncp, varid_out);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_del_att(void       *ncdp,
              int         varid,
              const char *name)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->del_att(ncbbp->ncp, varid, name);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_get_att(void         *ncdp,
              int           varid,
              const char   *name,
              void         *buf,
              MPI_Datatype  itype)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->get_att(ncbbp->ncp, varid, name, buf, itype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncbbio_put_att(void         *ncdp,
              int           varid,
              const char   *name,
              nc_type       xtype,
              MPI_Offset    nelems,
              const void   *buf,
              MPI_Datatype  itype)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->put_att(ncbbp->ncp, varid, name, xtype, nelems, buf,
                               itype);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
