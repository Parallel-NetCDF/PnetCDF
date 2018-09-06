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

#include <ncbbio_driver.h>

int
ncbbio_def_dim(void       *ncdp,
               const char *name,
               MPI_Offset  size,
               int        *dimidp)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->def_dim(ncbbp->ncp, name, size, dimidp);
    if (err != NC_NOERR) return err;

    /*
     * Record record dimension
     * Note: Assume only 1 rec dim
     */
    if (size == NC_UNLIMITED) ncbbp->recdimid = *dimidp;

    return NC_NOERR;
}

int
ncbbio_inq_dimid(void       *ncdp,
                 const char *name,
                 int        *dimid)
{
    NC_bb *ncbbp = (NC_bb*)ncdp;

    return ncbbp->ncmpio_driver->inq_dimid(ncbbp->ncp, name, dimid);
}

int
ncbbio_inq_dim(void       *ncdp,
               int         dimid,
               char       *name,
               MPI_Offset *sizep)
{
    int err;
    NC_bb *ncbbp = (NC_bb*)ncdp;

    err = ncbbp->ncmpio_driver->inq_dim(ncbbp->ncp, dimid, name, sizep);
    if (err != NC_NOERR) return err;

    /*
     * Update size of record dimension with pending records in the log
     * Note: Assume only 1 rec dim
     */
    if (sizep != NULL) {
        if (dimid == ncbbp->recdimid && *sizep < ncbbp->recdimsize)
            *sizep = ncbbp->recdimsize;
    }

    return NC_NOERR;
}

int
ncbbio_rename_dim(void       *ncdp,
                  int         dimid,
                  const char *newname)
{
    NC_bb *ncbbp = (NC_bb*)ncdp;

    return ncbbp->ncmpio_driver->rename_dim(ncbbp->ncp, dimid, newname);
}
