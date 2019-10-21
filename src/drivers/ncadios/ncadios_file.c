/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

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
#include <ncadios_internal.h>

int
ncadios_create(MPI_Comm     comm,
             const char  *path,
             int          cmode,
             int          ncid,
             MPI_Info     info,
             void       **ncpp)  /* OUT */
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_open(MPI_Comm     comm,
           const char  *path,
           int          omode,
           int          ncid,
           MPI_Info     info,
           void       **ncpp)
{
    int err, parse_done=0;
    NC_ad *ncadp;

    if (fIsSet(omode, NC_WRITE)){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }

    /* Create a NC_ad object and save its driver pointer */
    ncadp = (NC_ad*) NCI_Malloc(sizeof(NC_ad));
    if (ncadp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    ncadp->path = (char*) NCI_Malloc(strlen(path) + 1);
    if (ncadp->path == NULL) {
        NCI_Free(ncadp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(ncadp->path, path);
    ncadp->mode   = omode;
    ncadp->flag   = 0;
    ncadp->comm   = comm;
    ncadp->getsize = 0;
    MPI_Comm_rank(ncadp->comm, &(ncadp->rank));

    *ncpp = ncadp;

    /*
     * Use a modified bp2ncd utility to parse metadata related information
     * to guarantee the driver conforms to the converted nc file
     * We do not use bp2ncd for attributes, we rely on ADIOS read
     * API for attributes
     * Rank 0 parse the header and boardcast to other ranks
     */

    ncadiosi_var_list_init(&(ncadp->vars));
    ncadiosi_dim_list_init(&(ncadp->dims));

    if (ncadp->rank == 0) {
        err = ncadiosi_parse_header_bp2ncd(ncadp);
        if (err == 0){
            parse_done = 1;
        }
        else{
            parse_done = 0;
        }
    }

    /* Open with ADIOS read API */
    ncadp->fp = adios_read_open_file (path, ADIOS_READ_METHOD_BP, comm);
    if (ncadp->fp == NULL) {
        err = ncmpii_error_adios2nc(adios_errno, "Open");
        DEBUG_RETURN_ERROR(err);
    }

    if (ncadp->rank == 0) {
        /* bp2ncd does not support all type of files
         * In case it fail, we parse the metadata using our own rule
         */
        if (!parse_done){
            /* Reset var and dim list by free and realloc */
            ncadiosi_var_list_free(&(ncadp->vars));
            ncadiosi_dim_list_free(&(ncadp->dims));
            ncadiosi_var_list_init(&(ncadp->vars));
            ncadiosi_dim_list_init(&(ncadp->dims));

            ncadiosi_parse_header_readall(ncadp);
        }

        /* This require fp be opened */
        ncadiosi_parse_attrs(ncadp);
    }
    ncadios_sync_header(ncadp);

    /* Parse information regarding record dim */
    ncadiosi_parse_rec_dim(ncadp);

    /* Init non-blocking req list */
    ncadiosi_get_list_init(&(ncadp->getlist));

    return NC_NOERR;
}

int
ncadios_close(void *ncdp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    ncadiosi_var_list_free(&(ncadp->vars));
    ncadiosi_dim_list_free(&(ncadp->dims));

    if (ncadp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    err = adios_read_close(ncadp->fp);
    if (err != 0){
        err = ncmpii_error_adios2nc(adios_errno, "open");
        DEBUG_RETURN_ERROR(err);
    }

    ncadiosi_get_list_free(&(ncadp->getlist));
    NCI_Free(ncadp->path);
    NCI_Free(ncadp);

    return err;
}

int
ncadios_enddef(void *ncdp)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios__enddef(void       *ncdp,
              MPI_Offset  h_minfree,
              MPI_Offset  v_align,
              MPI_Offset  v_minfree,
              MPI_Offset  r_align)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_redef(void *ncdp)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_begin_indep_data(void *ncdp)
{
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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_inq(void *ncdp,
          int  *ndimsp,
          int  *nvarsp,
          int  *nattsp,
          int  *xtendimp)
{
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
    NC_ad *ncadp = (NC_ad*)ncdp;

    if (pathlen != NULL){
        *pathlen = strlen(ncadp->path);
    }

    if (path != NULL){
        strcpy(path, ncadp->path);
    }

    if (num_fix_varsp != NULL){
        /* All variables - number of record variables */
        int i, j;
        *num_fix_varsp = ncadp->vars.cnt;
        for(i = 0; i < ncadp->vars.cnt; i++){
            for(j = 0; j < ncadp->vars.data[i].ndim; j++){
                if (ncadp->dims.data[ncadp->vars.data[i].dimids[j]].len
                    == NC_UNLIMITED){
                    *num_rec_varsp -= 1;
                    break;
                }
            }
        }
    }

    if (num_rec_varsp != NULL){
        /* We count those variable with unlimited dim as rec variable */
        int i, j;
        *num_rec_varsp = 0;
        for(i = 0; i < ncadp->vars.cnt; i++){
            for(j = 0; j < ncadp->vars.data[i].ndim; j++){
                if (ncadp->dims.data[ncadp->vars.data[i].dimids[j]].len
                     == NC_UNLIMITED){
                    *num_rec_varsp += 1;
                    break;
                }
            }
        }
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
        *recsize = ncadp->nrec;
    }

    if (put_size != NULL){
        *put_size = 0;
    }

    if (get_size != NULL){
        *get_size = ncadp->getsize;
    }

    if (info_used != NULL){
        *info_used = MPI_INFO_NULL;
    }

    if (nreqs != NULL){
        *nreqs = ncadp->getlist.nused;
    }

    if (usage != NULL){
        *usage = 0;
    }

    if (buf_size != NULL){
        *buf_size = 0;
    }

    return NC_NOERR;
}

int
ncadios_cancel(void *ncdp,
             int   num_req,
             int  *req_ids,
             int  *statuses)
{
    /* Nonblocking IO does not support canceling due to ADIOS limitation */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_wait(void *ncdp,
           int   num_reqs,
           int  *req_ids,
           int  *statuses,
           int   reqMode)
{
    int err, status = NC_NOERR;
    int i;
    NC_ad *ncadp = (NC_ad*)ncdp;

    if (num_reqs == NC_REQ_ALL || num_reqs == NC_GET_REQ_ALL){
        /* Handle all active requests in the pool */
        err = ncadiosi_wait_all_get_req(ncadp);
        if (status == NC_NOERR){
            status = err;
        }
    }
    else{
        if (statuses == NULL){
            for(i = 0; i < num_reqs; i++){
                /* Handle request one by one */
                err = ncadiosi_wait_get_req(ncadp, req_ids[i], NULL);
                if (status == NC_NOERR){
                    status = err;
                }
            }
        }
        else{
            for(i = 0; i < num_reqs; i++){
                /* Handle request one by one */
                err = ncadiosi_wait_get_req(ncadp, req_ids[i], statuses + i);
                if (status == NC_NOERR){
                    status = err;
                }
            }
        }
    }

    return NC_NOERR;
}

int
ncadios_set_fill(void *ncdp,
               int   fill_mode,
               int  *old_fill_mode)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_fill_var_rec(void      *ncdp,
                   int        varid,
                   MPI_Offset recno)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_def_var_fill(void       *ncdp,
                   int         varid,
                   int         no_fill,
                   const void *fill_value)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_sync_numrecs(void *ncdp)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_sync(void *ncdp)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_flush(void *ncdp)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

