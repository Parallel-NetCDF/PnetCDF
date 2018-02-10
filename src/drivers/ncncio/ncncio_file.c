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
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 *
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

/* Note, netcdf header must come first due to conflicting constant definition */
#include <netcdf.h>
#include <netcdf_par.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strlen() */

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncncio_driver.h>

int
ncncio_create(MPI_Comm     comm,
             const char  *path,
             int          cmode,
             int          ncid,
             MPI_Info     info,
             void       **ncpp)  /* OUT */
{
    int err;
    void *ncp=NULL;
    NC_nc4 *nc4p;

    /* TODO: support NC4 write */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    /* Create a NC_nc4 object and save its driver pointer */
    nc4p = (NC_nc4*) NCI_Malloc(sizeof(NC_nc4));
    if (nc4p == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nc4p->path = (char*) NCI_Malloc(strlen(path)+1);
    if (nc4p->path == NULL) {
        NCI_Free(nc4p);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(nc4p->path, path);
    nc4p->mode   = cmode;
    nc4p->flag   = 0;
    nc4p->ncid  = ncid;
    nc4p->comm   = comm;

    *ncpp = nc4p;

    return NC_NOERR;
}

int
ncncio_open(MPI_Comm     comm,
           const char  *path,
           int          omode,
           int          ncid,
           MPI_Info     info,
           void       **ncpp)
{
    int err, format;
    int ncidtmp;
    NC_nc4 *nc4p;

    err = ncmpi_inq_file_format(path, &format);
    if (err != NC_NOERR) return err;

    /* Read only driver */
    if (omode & NC_WRITE){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
    }

    if (format == NC_FORMAT_NETCDF4_CLASSIC ||
        format == NC_FORMAT_NETCDF4) {
        /* Open with netcdf */
        err = nc_open_par(path, omode | NC_MPIIO, comm, info, &ncidtmp);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
    }
    else{
        DEBUG_RETURN_ERROR(NC_ENOTNC4);
    }

    /* Create a NC_nc4 object and save its driver pointer */
    nc4p = (NC_nc4*) NCI_Malloc(sizeof(NC_nc4));
    if (nc4p == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nc4p->path = (char*) NCI_Malloc(strlen(path)+1);
    if (nc4p->path == NULL) {
        NCI_Free(nc4p);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(nc4p->path, path);
    nc4p->mode   = omode;
    nc4p->flag   = 0;
    nc4p->ncid  = ncid;
    nc4p->comm   = comm;
    nc4p->ncid = ncidtmp;

    *ncpp = nc4p;

    return NC_NOERR;
}

int
ncncio_close(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    if (nc4p == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    /* Close with netcdf */
    err = nc_close(nc4p->ncid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    NCI_Free(nc4p->path);
    NCI_Free(nc4p);

    return err;
}

int
ncncio_enddef(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    return NC_NOERR;
}

int
ncncio__enddef(void       *ncdp,
              MPI_Offset  h_minfree,
              MPI_Offset  v_align,
              MPI_Offset  v_minfree,
              MPI_Offset  r_align)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    return NC_NOERR;
}

int
ncncio_redef(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    return NC_NOERR;
}

int
ncncio_begin_indep_data(void *ncdp)
{
    int i, err, nvar;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Get number of variables */
    err = nc_inq(nc4p->ncid, NULL, &nvar, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Set all variables to indep mode */
    for(i = 0; i < nvar; i++){
        err = nc_var_par_access(nc4p->ncid, i, NC_INDEPENDENT);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
    }

    return NC_NOERR;
}

int
ncncio_end_indep_data(void *ncdp)
{
    int i, err, nvar;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Get number of variables */
    err = nc_inq(nc4p->ncid, NULL, &nvar, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Set all variables to coll mode */
    for(i = 0; i < nvar; i++){
        err = nc_var_par_access(nc4p->ncid, i, NC_COLLECTIVE);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
    }

    return NC_NOERR;
}

int
ncncio_abort(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    if (nc4p == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    /* Call nc_abort */
    err = nc_abort(nc4p->ncid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    NCI_Free(nc4p->path);
    NCI_Free(nc4p);

    return err;
}

int
ncncio_inq(void *ncdp,
          int  *ndimsp,
          int  *nvarsp,
          int  *nattsp,
          int  *xtendimp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_inq */
    err = nc_inq(nc4p->ncid, ndimsp, nvarsp, nattsp, xtendimp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
ncncio_inq_misc(void       *ncdp,
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
    int i, err, flag, mpireturn;
    char value[MPI_MAX_INFO_VAL];
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Get the file pathname which was used to open/create the ncid's file.
     * path must already be allocated. Ignored if NULL */
    if (nc4p->path == NULL) {
        if (pathlen != NULL) *pathlen = 0;
        if (path    != NULL) *path = '\0';
    } else {
        if (pathlen != NULL) *pathlen = (int)strlen(nc4p->path);
        if (path    != NULL) strcpy(path, nc4p->path);
    }

    /* Calculate number of record and fix sized variables
     * NOTE: We assume there are only 1 record dims in NetCDF 4 classic model
     * NOTE: In NetCDF4 enhenced model, there can be multiple record dims, only the first one is considered as record dim here
     */
    if (num_fix_varsp != NULL || num_rec_varsp != NULL || recsize != NULL) {
        int nvar, ndim;
        int *dims;
        int *vars;
        int unlimdim;
        int nrec, nfix;

        /* Record dimid */
        err = nc_inq_unlimdim(nc4p->ncid, &unlimdim);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

        /* Record size */
        if (recsize != NULL){
            size_t udlen;
            err = nc_inq_dimlen(nc4p->ncid, unlimdim, &udlen);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
            *recsize = (MPI_Offset)udlen;
        }

        if (num_fix_varsp != NULL || num_rec_varsp != NULL){
            int j;

            /* Get all variables */
            err = nc_inq_varids(nc4p->ncid, &nvar, NULL);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
            vars = NCI_Malloc(SIZEOF_INT * nvar);
            if (vars == NULL){
                DEBUG_RETURN_ERROR(NC_ENOMEM);
            }
            err = nc_inq_varids(nc4p->ncid, NULL, vars);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

            /* Iterate through all variables */
            for (i = 0; i < nvar; i++){
                /* Get all dimensions */
                err = nc_inq_varndims(nc4p->ncid, vars[i], &ndim);
                if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
                dims = NCI_Malloc(SIZEOF_INT * ndim);
                if (dims == NULL){
                    DEBUG_RETURN_ERROR(NC_ENOMEM);
                }
                err = nc_inq_vardimid(nc4p->ncid, vars[i], dims);
                if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

                /* Iterate through all dimensions */
                for (j = 0; j < ndim; j++){
                    if (dims[j] == unlimdim){
                        break;
                    }
                }

                NCI_Free(dims);

                /* If non of the dimension is record dim, count as fixed var */
                if (j == ndim){
                    nfix += 1;
                }
                else{
                    nrec += 1;
                }
            }

            NCI_Free(vars);

            if (num_fix_varsp != NULL){
                *num_fix_varsp = nfix;
            }

            if (num_rec_varsp != NULL){
                *num_rec_varsp = nrec;
            }
        }

    }

    /* NetCDF does not expose any MPI related info */
    if (striping_size != NULL) {
        *striping_size = 0;
    }
    if (striping_count != NULL) {
        *striping_count = 0;
    }

    /* Read only */
    if (put_size != NULL) *put_size = 0;

    /* TODO: Calculate get size */
    if (get_size != NULL) *get_size = 0;

    /* NetCDF does not expose such info */
    if (header_size != NULL) *header_size = 0;
    if (header_extent != NULL) *header_extent = 0;

    /* TODO: suporting nonblocking IO */
    if (nreqs != NULL) {
        *nreqs = 0;
    }

    /* We don't use user buffer */
    if (usage != NULL) *usage = 0;

    /* We don't use user buffer */
    if (buf_size != NULL) *buf_size = 0;

    return NC_NOERR;
}

int
ncncio_cancel(void *ncdp,
             int   num_req,
             int  *req_ids,
             int  *statuses)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* TODO: Nonblocking IO */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    return NC_NOERR;
}

int
ncncio_wait(void *ncdp,
           int   num_reqs,
           int  *req_ids,
           int  *statuses,
           int   reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* TODO: Nonblocking IO */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    return NC_NOERR;
}

int
ncncio_set_fill(void *ncdp,
               int   fill_mode,
               int  *old_fill_mode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_set_fill */
    err = nc_set_fill(nc4p->ncid, fill_mode, old_fill_mode);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
ncncio_fill_var_rec(void      *ncdp,
                   int        varid,
                   MPI_Offset recno)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* NetCDF does not support this natively */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    return NC_NOERR;
}

int
ncncio_def_var_fill(void       *ncdp,
                   int         varid,
                   int         no_fill,
                   const void *fill_value)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_def_var_fill */
    err = nc_def_var_fill(nc4p->ncid, varid, no_fill, fill_value);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
ncncio_sync_numrecs(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Not supported in NetCDF, just sync everything */
    return ncncio_sync(ncdp);

    return NC_NOERR;
}

int
ncncio_sync(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_sync */
    err = nc_sync(nc4p->ncid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

