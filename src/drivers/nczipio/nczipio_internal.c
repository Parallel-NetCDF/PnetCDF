/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>
#include "nczipio_internal.h"

int
nczipioi_init(NC_zip *nczipp){
    int err;

    /* Initialize var list */
    err = nczipioi_var_list_init(&(nczipp->vars));
    if (err != NC_NOERR) return err;

    /* Initialize nonblocking list */
    err = nczipioi_req_list_init(&(nczipp->getlist));
    if (err != NC_NOERR) return err;
    err = nczipioi_req_list_init(&(nczipp->putlist));
    if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
    memset(&(nczipp->profile), 0, sizeof(NC_zip_timers));
#endif
}

int
nczipioi_parse_var_info(NC_zip *nczipp){
    int err;
    int vid;
    int i;
    int nvar;
    NC_zip_var var;

    err = nczipp->driver->inq(nczipp->ncp, NULL, &nvar, NULL, NULL);

    for(vid = 0; vid < nvar; vid++){
        err = nczipp->driver->get_att(nczipp->ncp, vid, "_varkind", &(var.varkind), MPI_INT);   // Comressed var?
        if (err != NC_NOERR || var.varkind == NC_ZIP_VAR_DATA){
            continue;
        }

        var.varid = vid;
        
        if (var.varkind == NC_ZIP_VAR_COMPRESSED){
            err = nczipp->driver->get_att(nczipp->ncp, var.varid, "_ndim", &(var.ndim), MPI_INT); // Original dimensions
            if (err != NC_NOERR) return err;

            var.dimids = (int*)NCI_Malloc(sizeof(int) * var.ndim);
            var.dimsize = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * var.ndim);

            err = nczipp->driver->get_att(nczipp->ncp, var.varid, "_dimids", var.dimids, MPI_INT);   // Dimensiona IDs
            if (err != NC_NOERR) return err;

            for(i = 0; i < var.ndim; i++){
                nczipp->driver->inq_dim(nczipp->ncp, var.dimids[i], NULL, var.dimsize + i);
            }

            err = nczipp->driver->get_att(nczipp->ncp, var.varid, "_datatype", &(var.xtype), MPI_INT); // Original datatype
            if (err != NC_NOERR) return err;

            var.esize = NC_Type_size(var.xtype);
            var.etype = ncmpii_nc2mpitype(var.xtype);
            var.chunkdim = NULL;

            nczipioi_var_init(nczipp, &var, 0);
        }
    
        err = nczipioi_var_list_add(&(nczipp->vars), var);
        if (err != NC_NOERR) return err;
    }

    return NC_NOERR;
}