/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_def_dim()    : dispatcher->def_dim()
 * ncmpi_inq_dimid()  : dispatcher->inq_dimid()
 * ncmpi_inq_dim()    : dispatcher->inq_dim()
 * ncmpi_rename_dim() : dispatcher->rename_dim()
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
ncfoo_def_dim(void       *ncdp,
              const char *name,
              MPI_Offset  size,
              int        *dimidp)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->def_dim(foo->ncp, name, size, dimidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_inq_dimid(void       *ncdp,
                const char *name,
                int        *dimid)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->inq_dimid(foo->ncp, name, dimid);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_inq_dim(void       *ncdp,
              int         dimid,
              char       *name,
              MPI_Offset *sizep)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->inq_dim(foo->ncp, dimid, name, sizep);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncfoo_rename_dim(void       *ncdp,
                 int         dimid,
                 const char *newname)
{
    int err;
    NC_foo *foo = (NC_foo*)ncdp;

    err = foo->driver->rename_dim(foo->ncp, dimid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
