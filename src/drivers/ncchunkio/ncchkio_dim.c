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
#include <ncchkio_driver.h>
#include "ncchkio_internal.h"

int
ncchkio_def_dim(void       *ncdp,
              const char *name,
              MPI_Offset  size,
              int        *dimidp)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->def_dim(ncchkp->ncp, name, size, dimidp);
    if (err != NC_NOERR) return err;

    if (size == NC_UNLIMITED){
        ncchkp->recdim = *dimidp;
    }

    return NC_NOERR;
}

int
ncchkio_inq_dimid(void       *ncdp,
                const char *name,
                int        *dimid)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->inq_dimid(ncchkp->ncp, name, dimid);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncchkio_inq_dim(void       *ncdp,
              int         dimid,
              char       *name,
              MPI_Offset *sizep)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->inq_dim(ncchkp->ncp, dimid, name, sizep);
    if (err != NC_NOERR) return err;

    if (dimid == ncchkp->recdim){   // update # records
        if (*sizep < ncchkp->recsize){
            *sizep = ncchkp->recsize;
        }
    }

    return NC_NOERR;
}

int
ncchkio_rename_dim(void       *ncdp,
                 int         dimid,
                 const char *newname)
{
    int err;
    NC_chk *ncchkp = (NC_chk*)ncdp;

    err = ncchkp->driver->rename_dim(ncchkp->ncp, dimid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
