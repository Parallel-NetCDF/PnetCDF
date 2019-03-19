/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs
 *
 * ncmpi_create()           : dispatcher->create()
 * ncmpi_open()             : dispatcher->open()
 * ncmpi_close()            : dispatcher->close()
 * ncmpi_enddef()           : dispatcher->enddef()
 * ncmpi__enddef()          : dispatcher->_enddef()
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_inq_misc()         : dispatcher->inq_misc()
 * ncmpi_wait()             : dispatcher->wait()
 * ncmpi_wait_all()         : dispatcher->wait()
 * ncmpi_cancel()           : dispatcher->cancel()
 *
 * ncmpi_set_fill()         : dispatcher->set_fill()
 * ncmpi_fill_var_rec()     : dispatcher->fill_rec()
 * ncmpi_def_var_fill()     : dispatcher->def_var_fill()
 * ncmpi_inq_var_fill()     : dispatcher->inq()
 *
 * ncmpi_sync()             : dispatcher->sync()
 * ncmpi_flush()             : dispatcher->flush()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 *
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strlen() */

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>
#include "nczipio_internal.h"

int
nczipio_create(MPI_Comm     comm,
             const char  *path,
             int          cmode,
             int          ncid,
             MPI_Info     info,
             void       **ncpp)  /* OUT */
{
    int err;
    int one = 1;
    void *ncp=NULL;
    NC_zip *nczipp;
    PNC_driver *driver=NULL;

    /* TODO: use comde to determine the true driver */
    driver = ncmpio_inq_driver();
    if (driver == NULL) return NC_ENOTNC;

    err = driver->create(comm, path, cmode, ncid, info, &ncp);
    if (err != NC_NOERR) return err;

    /* Create a NC_zip object and save its driver pointer */
    nczipp = (NC_zip*) NCI_Malloc(sizeof(NC_zip));
    if (nczipp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nczipp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (nczipp->path == NULL) {
        NCI_Free(nczipp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(nczipp->path, path);
    nczipp->mode   = cmode;
    nczipp->driver = driver;
    nczipp->flag   = 0;
    nczipp->ncp    = ncp;
    nczipp->comm   = comm;
    MPI_Comm_rank(comm, &(nczipp->rank));
    MPI_Comm_size(comm, &(nczipp->np));

    err = nczipioi_extract_hint(nczipp, info);
    if (err != NC_NOERR) return err;

    err = driver->put_att(nczipp->ncp, NC_GLOBAL, "_comressed", NC_INT, 1, &one, MPI_INT); // Mark this file as compressed
    if (err != NC_NOERR) return err;

    nczipioi_init(nczipp);

    *ncpp = nczipp;

    return NC_NOERR;
}

int
nczipio_open(MPI_Comm     comm,
           const char  *path,
           int          omode,
           int          ncid,
           MPI_Info     info,
           void       **ncpp)
{
    int err, format;
    void *ncp=NULL;
    NC_zip *nczipp;
    PNC_driver *driver=NULL;

    err = ncmpi_inq_file_format(path, &format);
    if (err != NC_NOERR) return err;

    if (format == NC_FORMAT_CLASSIC ||
        format == NC_FORMAT_CDF2 ||
        format == NC_FORMAT_CDF5) {
        driver = ncmpio_inq_driver();
    }
    if (driver == NULL) return NC_ENOTNC;

    err = driver->open(comm, path, omode, ncid, info, &ncp);
    if (err != NC_NOERR) return err;

    /* Create a NC_zip object and save its driver pointer */
    nczipp = (NC_zip*) NCI_Malloc(sizeof(NC_zip));
    if (nczipp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nczipp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (nczipp->path == NULL) {
        NCI_Free(nczipp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(nczipp->path, path);
    nczipp->mode   = omode;
    nczipp->driver = driver;
    nczipp->flag   = 0;
    nczipp->ncp    = ncp;
    nczipp->comm   = comm;
    MPI_Comm_rank(comm, &(nczipp->rank));
    MPI_Comm_size(comm, &(nczipp->np));

    nczipioi_init(nczipp);

    *ncpp = nczipp;

    return NC_NOERR;
}

int
nczipio_close(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (nczipp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    err = nczipp->driver->close(nczipp->ncp);

    err = nczipioi_var_list_free(&(nczipp->vars));

    NCI_Free(nczipp->path);
    NCI_Free(nczipp);

    return err;
}

int
nczipio_enddef(void *ncdp)
{
    int i, err;
    NC_zip *nczipp = (NC_zip*)ncdp;
    
    for(i = 0; i < nczipp->vars.cnt; i++){
        nczipioi_var_init(nczipp, nczipp->vars.data + i);
    }

    err = nczipp->driver->enddef(nczipp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio__enddef(void       *ncdp,
              MPI_Offset  h_minfree,
              MPI_Offset  v_align,
              MPI_Offset  v_minfree,
              MPI_Offset  r_align)
{
    int i, err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    for(i = 0; i < nczipp->vars.cnt; i++){
        nczipioi_var_init(nczipp, nczipp->vars.data + i);
    }

    err = nczipp->driver->_enddef(nczipp->ncp, h_minfree, v_align, v_minfree,
                               r_align);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_redef(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->redef(nczipp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_begin_indep_data(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->begin_indep_data(nczipp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_end_indep_data(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->end_indep_data(nczipp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_abort(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (nczipp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    err = nczipp->driver->abort(nczipp->ncp);

    NCI_Free(nczipp->path);
    NCI_Free(nczipp);

    return err;
}

int
nczipio_inq(void *ncdp,
          int  *ndimsp,
          int  *nvarsp,
          int  *nattsp,
          int  *xtendimp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->inq(nczipp->ncp, ndimsp, nvarsp, nattsp, xtendimp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_inq_misc(void       *ncdp,
               int        *pathlen,
               char       *path,
               int        *num_fix_varsp,
               int        *num_rec_varsp,
               int        *striping_size,
               int        *striping_count,
               MPI_Offset *header_size,
               MPI_Offset *header_extent,
               MPI_Offset *recsize,
               MPI_Offset *put_size,
               MPI_Offset *get_size,
               MPI_Info   *info_used,
               int        *nreqs,
               MPI_Offset *usage,
               MPI_Offset *buf_size)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->inq_misc(nczipp->ncp, pathlen, path, num_fix_varsp,
                                num_rec_varsp, striping_size, striping_count,
                                header_size, header_extent, recsize, put_size,
                                get_size, info_used, nreqs, usage, buf_size);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_cancel(void *ncdp,
             int   num_req,
             int  *req_ids,
             int  *statuses)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->cancel(nczipp->ncp, num_req, req_ids, statuses);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_wait(void *ncdp,
           int   num_reqs,
           int  *req_ids,
           int  *statuses,
           int   reqMode)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    nczipioi_wait(nczipp, num_reqs, req_ids, statuses, reqMode);

    return NC_NOERR;
}

int
nczipio_set_fill(void *ncdp,
               int   fill_mode,
               int  *old_fill_mode)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->set_fill(nczipp->ncp, fill_mode, old_fill_mode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_fill_var_rec(void      *ncdp,
                   int        varid,
                   MPI_Offset recno)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->fill_var_rec(nczipp->ncp, varid, recno);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_def_var_fill(void       *ncdp,
                   int         varid,
                   int         no_fill,
                   const void *fill_value)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->def_var_fill(nczipp->ncp, varid, no_fill, fill_value);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_sync_numrecs(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->sync_numrecs(nczipp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_sync(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->sync(nczipp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_flush(void *ncdp)
{
    int err;
    NC_zip *nczipp = (NC_zip*)ncdp;

    err = nczipp->driver->flush(nczipp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipioi_init(NC_zip *nczipp){
    int err;

    /* Initialize var list */
    err = nczipioi_var_list_init(&(nczipp->vars));
    if (err != NC_NOERR) return err;

    /* Select compression driver based on hint */
    switch (nczipp->zipdriver){
        case NC_ZIP_DRIVER_DUMMY:
            nczipp->zip = nczip_dummy_inq_driver();
        break;
    }
}