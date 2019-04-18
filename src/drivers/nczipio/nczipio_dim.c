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
#include <nczipio_driver.h>
#include "nczipio_internal.h"

int
nczipio_def_dim(void       *ncdp,
              const char *name,
              MPI_Offset  size,
              int        *dimidp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->def_dim(nczipp->ncp, name, size, dimidp);
    if (err != NC_NOERR) return err;

    if (size == NC_UNLIMITED){
        nczipp->recdim = *dimidp;
    }

    return NC_NOERR;
}

int
nczipio_inq_dimid(void       *ncdp,
                const char *name,
                int        *dimid)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->inq_dimid(nczipp->ncp, name, dimid);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_inq_dim(void       *ncdp,
              int         dimid,
              char       *name,
              MPI_Offset *sizep)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->inq_dim(nczipp->ncp, dimid, name, sizep);
    if (err != NC_NOERR) return err;

    if (dimid == nczipp->recdim){   // update # records
        if (*sizep < nczipp->recsize){
            *sizep = nczipp->recsize;
        }
    }

    return NC_NOERR;
}

int
nczipio_rename_dim(void       *ncdp,
                 int         dimid,
                 const char *newname)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->rename_dim(nczipp->ncp, dimid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
