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
 * ncmpi_put_att()     : dispatcher->inq_put_att()
 *
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
nc4io_inq_attname(void *ncdp,
                  int   varid,
                  int   attid,
                  char *name)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_inq_attname */
    err = nc_inq_attname(nc4p->ncid, varid, attid, name);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_inq_attid(void       *ncdp,
                int         varid,
                const char *name,
                int        *attidp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_inq_attid */
    err = nc_inq_attid(nc4p->ncid, varid, name, attidp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_inq_att(void       *ncdp,
              int         varid,
              const char *name,
              nc_type    *datatypep,
              MPI_Offset *lenp)
{
    int err;
    size_t len;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_inq_att */
    err = nc_inq_att(nc4p->ncid, varid, name, datatypep, &len);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    if (lenp != NULL) *lenp = (MPI_Offset)len;

    return NC_NOERR;
}

int
nc4io_rename_att(void       *ncdp,
                 int         varid,
                 const char *name,
                 const char *newname)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* New name can not be longer than old one in data mode */
    if (!fIsSet(nc4p->flag, NC_MODE_DEF)){
        if (strlen(newname) > strlen(name)){
            DEBUG_RETURN_ERROR(NC_ENOTINDEFINE);
        }
    }

    /* Call nc_rename_att */
    err = nc_rename_att(nc4p->ncid, varid, name, newname);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}


int
nc4io_copy_att(void       *ncdp_in,
               int         varid_in,
               const char *name,
               void       *ncdp_out,
               int         varid_out)
{
    int err;
    NC_nc4 *nc4p_in  = (NC_nc4*)ncdp_in;
    NC_nc4 *nc4p_out = (NC_nc4*)ncdp_out;

    /* Call nc_copy_att */
    err = nc_copy_att(nc4p_in->ncid, varid_in, name, nc4p_out->ncid, varid_out);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_del_att(void       *ncdp,
              int         varid,
              const char *name)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_del_att */
    err = nc_del_att(nc4p->ncid, varid, name);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}


/*
nc4io_get_att is implemented in nc4io_get_put.m4
*/

/*
nc4io_put_att is implemented in nc4io_get_put.m4
*/
