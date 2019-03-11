/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>
#include "nczipio_internal.h"


int
nczipioi_wait(NC_zip *nczipp, int nreqs, int *reqids, int *stats, int reqMode){
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->wait(nczipp->ncp, num_reqs, req_ids, statuses, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipioi_wait(NC_zip *nczipp, int nreqs, int *reqids, int *stats, int reqMode){
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->wait(nczipp->ncp, num_reqs, req_ids, statuses, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}