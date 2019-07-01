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
#ifdef PNETCDF_PROFILING
    double t0, t1;

    t0 = MPI_Wtime();
#endif

    /* TODO: use comde to determine the true driver */
    driver = ncmpio_inq_driver();
    if (driver == NULL) return NC_ENOTNC;

    err = driver->create(comm, path, cmode | NC_64BIT_DATA, ncid, info, &ncp);
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
    nczipp->mode   = cmode | NC_WRITE;
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

    nczipioi_init(nczipp, 1);

    *ncpp = nczipp;

    // Timer array is not avaiable until init, can't use NC_ZIP_TIMER_START
#ifdef PNETCDF_PROFILING
    t0 = MPI_Wtime() - t0;
    nczipp->profile.tt[NC_ZIP_TIMER_INIT] += t0;
    nczipp->profile.tt[NC_ZIP_TIMER_TOTAL] += t0;
#endif

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
    int err;
    int one = 0;
    void *ncp=NULL;
    NC_zip *nczipp = NULL;
    PNC_driver *driver=NULL;
#ifdef PNETCDF_PROFILING
    double t0;

    t0 = MPI_Wtime();
#endif

    /* TODO: use comde to determine the true driver */
    driver = ncmpio_inq_driver();
    if (driver == NULL) { 
        DEBUG_ASSIGN_ERROR(err, NC_ENOTNC)
        goto errout;
    }

    err = driver->open(comm, path, omode, ncid, info, &ncp);
    if (err != NC_NOERR) goto errout;

    /* Create a NC_zip object and save its driver pointer */
    nczipp = (NC_zip*) NCI_Malloc(sizeof(NC_zip));
    if (nczipp == NULL){
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto errout;
    }

    nczipp->path = (char*) NCI_Malloc(strlen(path)+1);
    if (nczipp->path == NULL) {
        NCI_Free(nczipp);
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto errout;
    }
    strcpy(nczipp->path, path);
    nczipp->mode   = omode;
    nczipp->driver = driver;
    if (nczipp->mode & NC_WRITE){
        nczipp->flag = 0;
    }
    else{
        nczipp->flag |= NC_MODE_RDONLY;
    }
    nczipp->ncp    = ncp;
    nczipp->comm   = comm;
    MPI_Comm_rank(comm, &(nczipp->rank));
    MPI_Comm_size(comm, &(nczipp->np));

    err = nczipioi_extract_hint(nczipp, info);
    if (err != NC_NOERR) goto errout;

    err = driver->get_att(nczipp->ncp, NC_GLOBAL, "_comressed", &one, MPI_INT); // Mark this file as compressed
    if (err != NC_NOERR){
        if (err == NC_ENOTATT){
            err = NC_EINVAL;
        }
         goto errout;
    }
    
    // Not compressed file
    if (one != 1){
        NCI_Free(nczipp->path);
        NCI_Free(nczipp);
        DEBUG_RETURN_ERROR(NC_EINVAL)
    }

    nczipioi_init(nczipp, 0);

    nczipioi_parse_var_info(nczipp);

    *ncpp = nczipp;

    // Timer array is not avaiable until init, can't use NC_ZIP_TIMER_START
#ifdef PNETCDF_PROFILING
    t0 = MPI_Wtime() - t0;
    nczipp->profile.tt[NC_ZIP_TIMER_INIT] += t0;
    nczipp->profile.tt[NC_ZIP_TIMER_TOTAL] += t0;
#endif

    return NC_NOERR;

errout:
    if (ncp != NULL){
        driver->close(nczipp->ncp);
    }  
    if (nczipp != NULL){
        if (nczipp->path != NULL){
            NCI_Free(nczipp->path);
        }      
        NCI_Free(nczipp);
    }

    return err;
}

int
nczipio_close(void *ncdp)
{
    int err;
#ifdef PNETCDF_PROFILING
    MPI_Offset put_size, get_size;
    char *_env_str = getenv("PNETCDF_SHOW_PERFORMANCE_INFO");
#endif
    NC_zip *nczipp = (NC_zip*)ncdp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_FINALIZE)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_TOTAL)

    if (nczipp == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    if (!(nczipp->flag & NC_MODE_RDONLY)){
        int i;

        NC_ZIP_TIMER_START(NC_ZIP_TIMER_FINALIZE_META)

        err = nczipp->driver->redef(nczipp->ncp);
        if (err != NC_NOERR){
            return err;
        }

        //record chunk dim
        for(i = 0; i < nczipp->vars.cnt; i++){
            if (nczipp->vars.data[i].isnew){            
                err = nczipp->driver->put_att(nczipp->ncp, nczipp->vars.data[i].varid, "_chunkdim", NC_INT, nczipp->vars.data[i].ndim, nczipp->vars.data[i].chunkdim, MPI_INT);
                if (err != NC_NOERR){
                    return err;
                }
                err = nczipp->driver->put_att(nczipp->ncp, nczipp->vars.data[i].varid, "_zipdriver", NC_INT, 1, &(nczipp->vars.data[i].zipdriver), MPI_INT);
                if (err != NC_NOERR){
                    return err;
                }
            }
        }

        // Record recsize
        err = nczipp->driver->put_att(nczipp->ncp, NC_GLOBAL, "_recsize", NC_INT64, 1, &(nczipp->recsize), MPI_LONG_LONG); // Mark this file as compressed
        if (err != NC_NOERR) return err;

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_FINALIZE_META)
    }

#ifdef PNETCDF_PROFILING
    err = nczipp->driver->inq_misc(nczipp->ncp, NULL, NULL, NULL,
                            NULL, NULL, NULL,
                            NULL, NULL, NULL, &put_size,
                            &get_size, NULL, NULL, NULL, NULL);
    nczipp->putsize += put_size;
    nczipp->getsize += get_size;
#endif
    err = nczipp->driver->close(nczipp->ncp);

    err = nczipioi_var_list_free(&(nczipp->vars));

    nczipioi_req_list_free(&(nczipp->putlist));
    nczipioi_req_list_free(&(nczipp->getlist));

#ifdef PNETCDF_PROFILING
    if (_env_str != NULL && *_env_str != '0') {                        
        nczipioi_print_profile(nczipp);
    }            
#endif

    NCI_Free(nczipp->path);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_FINALIZE)
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_TOTAL)

    NCI_Free(nczipp);

    return err;
}

int
nczipio_enddef(void *ncdp)
{
    int i, err;
    MPI_Offset rsize;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_TOTAL)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT)

    // Reserve header space
    rsize = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        varp = nczipp->vars.data + i;
        if (varp->varkind == NC_ZIP_VAR_COMPRESSED){
            rsize += ((8 + 16) + 4 + 8 + 4) * 8 + ((8 + 16) + 4 + 8 + 4 * varp->ndim) * 2; // Atts
            rsize += ((8 + 32) + 8) * 3; // dims
            rsize += ((8 + 32) + 8 + 8 + (8 + 8 + (8 + 12 + 4 + 8 + 4)) + 4 + 8 + 8) * 3; // vars
        }
        else{
            rsize += ((8 + 16) + 4 + 8 + 4); // Atts
        }
    }
    rsize *= 4;   // 4 times for future expension

    err = nczipp->driver->_enddef(nczipp->ncp, rsize, 0, 0, 0);
    if (err != NC_NOERR) return err;

    if (!(nczipp->delay_init)){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT_META)

        for(i = 0; i < nczipp->vars.cnt; i++){
            nczipioi_var_init(nczipp, nczipp->vars.data + i, 0, NULL, NULL);
        }
        
        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT_META)
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT)
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_TOTAL)

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
    MPI_Offset rsize;
    NC_zip_var *varp;
    NC_zip *nczipp = (NC_zip*)ncdp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_TOTAL)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT)

    // Reserve header space
    rsize = 0;
    for(i = 0; i < nczipp->vars.cnt; i++){
        varp = nczipp->vars.data + i;
        if (varp->varkind == NC_ZIP_VAR_COMPRESSED){
            rsize += ((8 + 16) + 4 + 8 + 4) * 8 + ((8 + 16) + 4 + 8 + 4 * varp->ndim) * 2; // Atts
            rsize += ((8 + 32) + 8) * 3; // dims
            rsize += ((8 + 32) + 8 + 8 + (8 + 8 + (8 + 12 + 4 + 8 + 4)) + 4 + 8 + 8) * 3; // vars
        }
        else{
            rsize += ((8 + 16) + 4 + 8 + 4); // Atts
        }
    }
    rsize *= 4;   // 4 times for future expension

    err = nczipp->driver->_enddef(nczipp->ncp, h_minfree + rsize, v_align, v_minfree,
                               r_align);
    if (err != NC_NOERR) return err;
    
    if (!(nczipp->delay_init)){
        NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT_META)

        for(i = 0; i < nczipp->vars.cnt; i++){
            nczipioi_var_init(nczipp, nczipp->vars.data + i, 0, NULL, NULL);
        }

        NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT_META)
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT)
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_TOTAL)

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

    err = nczipp->driver->inq(nczipp->ncp, ndimsp, NULL, nattsp, xtendimp);
    if (err != NC_NOERR) return err;

    if (nvarsp != NULL){
        *nvarsp = nczipp->vars.cnt;
    }

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

    if (num_fix_varsp != NULL){
        *num_fix_varsp = nczipp->vars.cnt;
    }

    if (nreqs != NULL){
        *nreqs = nczipp->putlist.nused + nczipp->getlist.nused;
    }

    if (put_size != NULL){
        *put_size += nczipp->putsize;
    }

    if (get_size != NULL){
        *get_size += nczipp->getsize;
    }

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
    int err, status = NC_NOERR;
    int i;
    int ncom = 0, nraw = 0;
    int *rawreqs = NULL, *comreqs = NULL;
    int *rawstats = NULL, *comstats = NULL;
    NC_zip *nczipp = (NC_zip*)ncdp;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_TOTAL)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB)
    NC_ZIP_TIMER_START(NC_ZIP_TIMER_NB_WAIT)

    if (num_reqs < 0){  // NC_REQ_ALL || nreqs == NC_PUT_REQ_ALL || nreqs == NC_GET_REQ_ALL
        err = nczipioi_wait(nczipp, num_reqs, NULL, NULL, reqMode);
        if (status == NC_NOERR){
            status = err;
        }
        err = nczipp->driver->wait(nczipp->ncp, num_reqs, NULL, NULL, reqMode);
        if (status == NC_NOERR){
            status = err;
        }
        return status;
    }

    if (num_reqs > 0){
        // Count number of get and put requests
        for(i = 0; i < num_reqs; i++){
            if (req_ids[i] & 1){
                nraw++;
            }
        }

        // Allocate buffer
        ncom = num_reqs - nraw;
        rawreqs = (int*)NCI_Malloc(sizeof(int) * nraw);
        comreqs = (int*)NCI_Malloc(sizeof(int) * ncom);
        
        // Build put and get req list
        nraw = ncom = 0;
        for(i = 0; i < num_reqs; i++){
            if (req_ids[i] & 1){
                rawreqs[nraw++] = req_ids[i] >> 1;
            }
            else{
                comreqs[ncom++] = req_ids[i] >> 1;
            }
        }
    }

    if (statuses != NULL){
        rawstats = (int*)NCI_Malloc(sizeof(int) * nraw);
        comstats = (int*)NCI_Malloc(sizeof(int) * ncom);
    }
    else{
        rawstats = NULL;
        comstats = NULL;
    }

    if (nraw > 0){
        err = nczipp->driver->wait(nczipp->ncp, nraw, rawreqs, rawstats, reqMode);
        if (status == NC_NOERR){
            status = err;
        }
    }
    
    if (ncom > 0){
        err = nczipioi_wait(nczipp, ncom, comreqs, comstats, reqMode);
        if (status == NC_NOERR){
            status = err;
        }
    }

    // Assign stats
    if (statuses != NULL){
        nraw = ncom = 0;
        for(i = 0; i < num_reqs; i++){
            if (req_ids[i] & 1){
                statuses[i] = rawstats[nraw++];
            }
            else{
                statuses[i] = comstats[ncom++];
            }
        }

        NCI_Free(rawstats);
        NCI_Free(comstats);
    }
    
    NCI_Free(rawreqs);
    NCI_Free(comreqs);
    
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_TOTAL)
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB)
    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_NB_WAIT)

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
