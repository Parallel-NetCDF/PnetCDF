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

/* Note, netcdf header must come first due to conflicting constant definition */
#include <netcdf.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nc4io_driver.h>

int
nc4io_def_dim(void       *ncdp,
              const char *name,
              MPI_Offset  size,
              int        *dimidp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_def_dim */
    err = nc_def_dim(nc4p->ncid, name, (size_t) size, dimidp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_inq_dimid(void       *ncdp,
                const char *name,
                int        *dimid)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_inq_dimid */
    err = nc_inq_dimid(nc4p->ncid, name, dimid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_inq_dim(void       *ncdp,
              int         dimid,
              char       *name,
              MPI_Offset *sizep)
{
    int err;
    size_t len;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_inq_dim */
    err = nc_inq_dim(nc4p->ncid, dimid, name, &len);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Convert to MPI_Offset */
    if (sizep != NULL){
        *sizep = (MPI_Offset)len;
    }

    return NC_NOERR;
}

int
nc4io_rename_dim(void       *ncdp,
                 int         dimid,
                 const char *newname)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* New name can not be longer than old one in data mode */
    if (!fIsSet(nc4p->flag, NC_MODE_DEF)){
        char oldname[NC_MAX_NAME + 1];
        err = nc_inq_dimname(nc4p->ncid, dimid, oldname);
        if (err != NC_NOERR){
            DEBUG_RETURN_ERROR(err);
        }
        if (strlen(newname) > strlen(oldname)){
            DEBUG_RETURN_ERROR(NC_ENOTINDEFINE);
        }
    }

    /* Call nc_rename_dim */
    err = nc_rename_dim(nc4p->ncid, dimid, newname);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}
