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
#include <ncdwio_driver.h>

int
ncdwio_def_dim(void       *ncdp,
              const char *name,
              MPI_Offset  size,
              int        *dimidp)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->def_dim(ncdwp->ncp, name, size, dimidp);
    if (err != NC_NOERR) return err;
    
    /* 
     * Record record dimension
     * Note: Assume only 1 rec dim
     */
    if (size == NC_UNLIMITED){
        ncdwp->recdimid = *dimidp;
    }

    return NC_NOERR;
}

int
ncdwio_inq_dimid(void       *ncdp,
                const char *name,
                int        *dimid)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->inq_dimid(ncdwp->ncp, name, dimid);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncdwio_inq_dim(void       *ncdp,
              int         dimid,
              char       *name,
              MPI_Offset *sizep)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->inq_dim(ncdwp->ncp, dimid, name, sizep);
    if (err != NC_NOERR) return err;
    
    /* 
     * Update size of record dimension with pending records in the log
     * Note: Assume only 1 rec dim
     */
    if (dimid == ncdwp->recdimid){
        if (*sizep < ncdwp->recdimsize){
            *sizep = ncdwp->recdimsize;
        }
    }

    return NC_NOERR;
}

int
ncdwio_rename_dim(void       *ncdp,
                 int         dimid,
                 const char *newname)
{
    int err;
    NC_dw *ncdwp = (NC_dw*)ncdp;
    
    err = ncdwp->ncmpio_driver->rename_dim(ncdwp->ncp, dimid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
