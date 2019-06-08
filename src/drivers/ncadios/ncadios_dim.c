/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

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
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncadios_driver.h>
#include <ncadios_internal.h>

int
ncadios_def_dim(void       *ncdp,
              const char *name,
              MPI_Offset  size,
              int        *dimidp)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_inq_dimid(void       *ncdp,
                const char *name,
                int        *dimid)
{
    NC_ad *ncadp = (NC_ad*)ncdp;

    return ncadiosi_inq_dimid(ncadp, (char*)name, dimid);
}

int
ncadios_inq_dim(void       *ncdp,
              int         dimid,
              char       *name,
              MPI_Offset *sizep)
{
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* ADIOS read API does not expose dimension information.
     * We rely on a modified bp2ncd utility to build up our
     * own dimensional list.
     */
    if (name != NULL){
        strcpy(name, ncadp->dims.data[dimid].name);
    }

    if (sizep != NULL){
        *sizep = ncadp->dims.data[dimid].len;
    }

    return NC_NOERR;
}

int
ncadios_rename_dim(void       *ncdp,
                 int         dimid,
                 const char *newname)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}
