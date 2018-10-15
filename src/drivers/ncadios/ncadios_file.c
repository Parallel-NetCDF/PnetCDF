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
#include <ncadios_driver.h>

int
ncadios_create(MPI_Comm     comm,
             const char  *path,
             int          cmode,
             int          ncid,
             MPI_Info     info,
             void       **ncpp)  /* OUT */
{
    int err;
    void *ncp=NULL;
    NC_ad *ncadp;
    PNC_driver *driver=NULL;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    /* TODO: use comde to determine the true driver */

    /* Create a NC_ad object and save its driver pointer */
    ncadp = (NC_ad*) NCI_Malloc(sizeof(NC_ad));
    if (ncadp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    ncadp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (ncadp->path == NULL) {
        NCI_Free(ncadp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(ncadp->path, path);
    ncadp->mode   = cmode;
    ncadp->flag   = 0;
    ncadp->comm   = comm;

    *ncpp = ncadp;

    return NC_NOERR;
}

int
ncadios_open(MPI_Comm     comm,
           const char  *path,
           int          omode,
           int          ncid,
           MPI_Info     info,
           void       **ncpp)
{
    int err, format;
    int i;
    void *ncp=NULL;
    NC_ad *ncadp;
    PNC_driver *driver=NULL;

    if (fIsSet(omode, NC_WRITE)){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }


    /* Create a NC_ad object and save its driver pointer */
    ncadp = (NC_ad*) NCI_Malloc(sizeof(NC_ad));
    if (ncadp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    
    ncadp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (ncadp->path == NULL) {
        NCI_Free(ncadp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(ncadp->path, path);
    ncadp->mode   = omode;
    ncadp->flag   = 0;
    ncadp->comm   = comm;

    *ncpp = ncadp;

    ncadiosi_var_list_init(&(ncadp->vars));
    ncadiosi_att_list_init(&(ncadp->atts));
    ncadiosi_dim_list_init(&(ncadp->dims));

    ncadiosi_parse_header(ncadp);

    /* Open with adios */
    ncadp->fp = adios_read_open_file (path, ADIOS_READ_METHOD_BP, comm);
    if (ncadp->fp == NULL) {
        DEBUG_RETURN_ERROR(NC_EADIOS);
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    /* Build dimensionality list */
    ncadp->ndims = (int*)NCI_Malloc(sizeof(int) * ncadp->fp->nvars);
    for (i = 0; i < ncadp->fp->nvars; i++) {
        ADIOS_VARINFO *v = adios_inq_var_byid (ncadp->fp, i);
        adios_inq_var_stat (ncadp->fp, v, 0, 0);
        ncadp->ndims[i] = v->ndim;
    }

    return NC_NOERR;
}

int
ncadios_close(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    ncadiosi_var_list_free(&(ncadp->vars));
    ncadiosi_att_list_free(&(ncadp->atts));
    ncadiosi_dim_list_free(&(ncadp->dims));

    if (ncadp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    err = adios_read_close(ncadp->fp);
    if (err != 0){
        DEBUG_RETURN_ERROR(NC_EADIOS)
    }

    NCI_Free(ncadp->ndims);
    NCI_Free(ncadp->path);
    NCI_Free(ncadp);

    return err;
}

int
ncadios_enddef(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios__enddef(void       *ncdp,
              MPI_Offset  h_minfree,
              MPI_Offset  v_align,
              MPI_Offset  v_minfree,
              MPI_Offset  r_align)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_redef(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_begin_indep_data(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Make sure we are in data mode */
    if (fIsSet(ncadp->flag, NC_MODE_DEF)){
        DEBUG_RETURN_ERROR(NC_EINDEFINE);
    }

    /* Set indep flag */
    fSet(ncadp->flag, NC_MODE_INDEP);

    return NC_NOERR;
}

int
ncadios_end_indep_data(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Make sure we are in data mode */
    if (fIsSet(ncadp->flag, NC_MODE_DEF)){
        DEBUG_RETURN_ERROR(NC_EINDEFINE);
    }

    /* Clear indep flag */
    fClr(ncadp->flag, NC_MODE_INDEP);

    return NC_NOERR;
}

int
ncadios_abort(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    if (ncadp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    NCI_Free(ncadp->path);
    NCI_Free(ncadp);

    return err;
}

int
ncadios_inq(void *ncdp,
          int  *ndimsp,
          int  *nvarsp,
          int  *nattsp,
          int  *xtendimp)
{
    int err;
    int i;
    NC_ad *ncadp = (NC_ad*)ncdp;

    if (ndimsp != NULL){
        *ndimsp = ncadp->dims.cnt;
    }

    if (nvarsp != NULL){
        *nvarsp = ncadp->vars.cnt;
    }

    if (nattsp != NULL){
        *nattsp = ncadp->fp->nattrs;
    }

    if (xtendimp != NULL){
        *xtendimp = -1;
    }

    return NC_NOERR;
}

int
ncadios_inq_misc(void       *ncdp,
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    if (pathlen != NULL){
        *pathlen = strlen(ncadp->path);
    }

    if (path != NULL){
        strcpy(path, ncadp->path);
    }

    if (num_fix_varsp != NULL){
        *num_fix_varsp = ncadp->vars.cnt;
    }

    if (num_rec_varsp != NULL){
        *num_rec_varsp = 0;
    }

    if (striping_size != NULL){
        *striping_size = 0;
    }

    if (striping_count != NULL){
        *striping_count = 0;
    }

    if (header_size != NULL){
        *header_size = 0;
    }

    if (header_extent != NULL){
        *header_extent = 0;
    }

    if (recsize != NULL){
        *recsize = 0;
    }

    if (put_size != NULL){
        *put_size = 0;
    }

    //TODO: Count get size
    if (get_size != NULL){
        *get_size = 0;
    }

    if (info_used != NULL){
        *info_used = MPI_INFO_NULL;
    }

    //TODO: Wire up nonblocking req
    if (nreqs != NULL){
        *nreqs = 0;
    }

    //TODO: Wire up nonblocking req
    if (usage != NULL){
        *usage = 0;
    }

    //TODO: Wire up nonblocking req
    if (buf_size != NULL){
        *buf_size = MPI_INFO_NULL;
    }
    
    return NC_NOERR;
}

int
ncadios_cancel(void *ncdp,
             int   num_req,
             int  *req_ids,
             int  *statuses)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* TODO: Nonblocking IO support */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_wait(void *ncdp,
           int   num_reqs,
           int  *req_ids,
           int  *statuses,
           int   reqMode)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* TODO: Nonblocking IO support */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_set_fill(void *ncdp,
               int   fill_mode,
               int  *old_fill_mode)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_fill_var_rec(void      *ncdp,
                   int        varid,
                   MPI_Offset recno)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_def_var_fill(void       *ncdp,
                   int         varid,
                   int         no_fill,
                   const void *fill_value)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_sync_numrecs(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_sync(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_flush(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

