/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_def_var()                  : dispatcher->def_var()
 * ncmpi_inq_varid()                : dispatcher->inq_varid()
 * ncmpi_inq_var()                  : dispatcher->inq_var()
 * ncmpi_rename_var()               : dispatcher->rename_var()
 *
 * ncmpi_get_var<kind>()            : dispatcher->get_var()
 * ncmpi_put_var<kind>()            : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>()     : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>()     : dispatcher->put_var()
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 *
 * ncmpi_iget_var<kind>()           : dispatcher->iget_var()
 * ncmpi_iput_var<kind>()           : dispatcher->iput_var()
 * ncmpi_iget_var<kind>_<type>()    : dispatcher->iget_var()
 * ncmpi_iput_var<kind>_<type>()    : dispatcher->iput_var()
 *
 * ncmpi_buffer_attach()            : dispatcher->buffer_attach()
 * ncmpi_buffer_detach()            : dispatcher->buffer_detach()
 * ncmpi_bput_var<kind>_<type>()    : dispatcher->bput_var()
 *
 * ncmpi_get_varn_<type>()          : dispatcher->get_varn()
 * ncmpi_put_varn_<type>()          : dispatcher->put_varn()
 *
 * ncmpi_iget_varn_<type>()         : dispatcher->iget_varn()
 * ncmpi_iput_varn_<type>()         : dispatcher->iput_varn()
 * ncmpi_bput_varn_<type>()         : dispatcher->bput_varn()
 *
 * ncmpi_get_vard()                 : dispatcher->get_vard()
 * ncmpi_put_vard()                 : dispatcher->put_vard()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

/* Note, netcdf header must come first due to conflicting constant definition */
#include <netcdf.h>
#include <netcdf_par.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nc4io_driver.h>

int
nc4io_def_var(void       *ncdp,
              const char *name,
              nc_type     xtype,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_def_var */
    err = nc_def_var(nc4p->ncid, name, xtype, ndims, dimids, varidp);

    /* Default mode in NetCDF is indep, set to coll if in coll mode */
    if (!(nc4p->flag & NC_MODE_INDEP)){
        err = nc_var_par_access(nc4p->ncid, *varidp, NC_COLLECTIVE);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
    }

    return NC_NOERR;
}

int
nc4io_inq_varid(void       *ncdp,
                const char *name,
                int        *varid)
{
    int err;
    int vid;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_inq_varid */
    err = nc_inq_varid(nc4p->ncid, name, &vid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* NetCDF does not support NULL varid
     * When varid is NULL, NC_NOERR will always return even given invalid name 
     */
    if (varid != NULL){
        *varid = vid;
    }

    return NC_NOERR;
}

int
nc4io_inq_var(void       *ncdp,
              int         varid,
              char       *name,
              nc_type    *xtypep,
              int        *ndimsp,
              int        *dimids,
              int        *nattsp,
              MPI_Offset *offsetp,
              int        *no_fillp,
              void       *fill_valuep)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call NC_inq_var_all */
    err = NC_inq_var_all(nc4p->ncid, varid, name, xtypep, ndimsp, dimids, nattsp, 
                        NULL, NULL, NULL, NULL, NULL, NULL,
                        no_fillp, fill_valuep, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_rename_var(void       *ncdp,
                 int         varid,
                 const char *newname)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* New name can not be longer than old one in data mode */
    if (!fIsSet(nc4p->flag, NC_MODE_DEF)){
        char oldname[NC_MAX_NAME + 1];
        err = nc_inq_varname(nc4p->ncid, varid, oldname);
        if (err != NC_NOERR){
            DEBUG_RETURN_ERROR(err);
        }
        if (strlen(newname) > strlen(oldname)){
            DEBUG_RETURN_ERROR(NC_ENOTINDEFINE);
        }
    }

    /* Call nc_rename_var */
    err = nc_rename_var(nc4p->ncid, varid, newname);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

/* 
nc4io_get_var is implemented iin ncmpio_get_put.m4
*/

/* 
nc4io_put_var is implemented iin ncmpio_get_put.m4
*/

int
nc4io_iget_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               void             *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc4io_get_var */
    err = nc4io_get_var(ncdp, varid, start, count, stride, imap, buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    /* TODO: Issue dummy id */
    *reqid = NC_REQ_NULL;

    return NC_NOERR;
}

int
nc4io_iput_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* We do not support nonblocking I/O so far */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    /* Call nc4io_put_var */
    err = nc4io_put_var(ncdp, varid, start, count, stride, imap, buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    /* TODO: Issue dummy id */
    *reqid = NC_REQ_NULL;

    return NC_NOERR;
}

int
nc4io_buffer_attach(void       *ncdp,
                    MPI_Offset  bufsize)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* We don't use buffer */
    /* Todo, count buffer size to report in inq */

    return NC_NOERR;
}

int
nc4io_buffer_detach(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* We don't use buffer */
    /* Todo, count buffer size to report in inq */

    return NC_NOERR;
}

int
nc4io_bput_var(void             *ncdp,
               int               varid,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               const MPI_Offset *imap,
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype,
               int              *reqid,
               int               reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc4io_iput_var */
    err = nc4io_iput_var(ncdp, varid, start, count, stride, imap, buf, bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}
int
nc4io_get_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               void              *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int i, j, ndim, err, status = NC_NOERR;
    int elsize;
    int num_max, num_min;
    int isindep;
    size_t putsize;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* No support for varn at this time
     * NetCDF hang when collective call consists of 0 count on part of process 
     * It make varn difficult to be implemented efficiently 
     */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    isindep = fIsSet(nc4p->flag, NC_MODE_INDEP);

    /* Check arguments */
    if (starts == NULL){
        if (isindep){
            DEBUG_RETURN_ERROR(NC_ENULLSTART)
        }
        else{
            // Participate coll I/O
            num = 0;
            DEBUG_ASSIGN_ERROR(status, NC_ENULLSTART);
        }
    }

    /* Get variable dimensionality */
    err = nc_inq_varndims(nc4p->ncid, varid, &ndim);
    if (err != NC_NOERR){ 
        if (isindep){
            DEBUG_RETURN_ERROR(err)
        }
        else{
            // Participate coll I/O
            num = 0;
            if (status == NC_NOERR){
                DEBUG_ASSIGN_ERROR(status, err);
            }
        }
    }

    /* 
     * Get element size 
     * Use size of buftype for high-level API
     * Flexible APi not supported
     */
    if (bufcount == -1){
        MPI_Type_size(buftype, &elsize);
    }
    else{
        nc_type type;
        
        if (isindep){
            DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
        }
        else{
            // Participate coll I/O
            num = 0;
            if (status == NC_NOERR){
                DEBUG_ASSIGN_ERROR(status, NC_ENOTSUPPORT);
            }
        }
        
        err = nc_inq_vartype(nc4p->ncid, varid, &type);
        if (err != NC_NOERR){ 
            DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
        }
        
        elsize = nc4io_nc_type_size(type);
    }
    
    if (!isindep){
        MPI_Allreduce(&num, &num_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&num, &num_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (num_max != num_min){
            // Num mismatch, need to switch to indep mode
            err = nc_var_par_access(nc4p->ncid, varid, NC_INDEPENDENT);
            if (err != NC_NOERR){ 
                num = 0;
                if (status == NC_NOERR){
                    DEBUG_ASSIGN_ERROR(status, err);
                }
            }
        }
    }

    /* Call nc4io_put_var for num times */
    for(i = 0; i < num; i++){
            err = nc4io_get_var(ncdp, varid, starts[i], counts[i], NULL, NULL, buf, bufcount, buftype, reqMode);
            if (err != NC_NOERR){ 
                if (isindep){
                    DEBUG_RETURN_ERROR(err)
                }
                else{
                    // Participate coll I/O
                    if (status == NC_NOERR){
                        DEBUG_ASSIGN_ERROR(status, err);
                    }
                }
            }
            
            /* Calculate the size of each put_var */
            putsize = (size_t)elsize;
            for(j = 0; j < ndim; j++){
                putsize *= counts[i][j];
            }

            /* Move buffer pointer */
            buf = (void*)(((char*)buf) + putsize);
    }

    if (!isindep){
        if (num_max != num_min){
            // switch back to coll mode
            err = nc_var_par_access(nc4p->ncid, varid, NC_COLLECTIVE);
            if (err != NC_NOERR){ 
                if (status == NC_NOERR){
                    DEBUG_ASSIGN_ERROR(status, err);
                }
            }
        }
    }

    return status;
}

int
nc4io_put_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               const void        *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{
    int i, j, ndim, err, status = NC_NOERR;
    int elsize;
    int num_max, num_min;
    int isindep;
    size_t putsize;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* No support for varn at this time
     * NetCDF hang when collective call consists of 0 count on part of process 
     * It make varn difficult to be implemented efficiently 
     */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    isindep = fIsSet(nc4p->flag, NC_MODE_INDEP);

    /* Check arguments */
    if (starts == NULL){
        if (isindep){
            DEBUG_RETURN_ERROR(NC_ENULLSTART)
        }
        else{
            // Participate coll I/O
            num = 0;
            DEBUG_ASSIGN_ERROR(status, NC_ENULLSTART);
        }
    }

    /* Get variable dimensionality */
    err = nc_inq_varndims(nc4p->ncid, varid, &ndim);
    if (err != NC_NOERR){ 
        if (isindep){
            DEBUG_RETURN_ERROR(err)
        }
        else{
            // Participate coll I/O
            num = 0;
            if (status == NC_NOERR){
                DEBUG_ASSIGN_ERROR(status, err);
            }
        }
    }

    /* 
     * Get element size 
     * Use size of buftype for high-level API
     * Flexible APi not supported
     */
    if (bufcount == -1){
        MPI_Type_size(buftype, &elsize);
    }
    else{
        nc_type type;
        
        if (isindep){
            DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
        }
        else{
            // Participate coll I/O
            num = 0;
            if (status == NC_NOERR){
                DEBUG_ASSIGN_ERROR(status, NC_ENOTSUPPORT);
            }
        }
        
        err = nc_inq_vartype(nc4p->ncid, varid, &type);
        if (err != NC_NOERR){ 
            DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
        }
        
        elsize = nc4io_nc_type_size(type);
    }
    
    if (!isindep){
        MPI_Allreduce(&num, &num_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&num, &num_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (num_max != num_min){
            // Num mismatch, need to switch to indep mode
            err = nc_var_par_access(nc4p->ncid, varid, NC_INDEPENDENT);
            if (err != NC_NOERR){ 
                num = 0;
                if (status == NC_NOERR){
                    DEBUG_ASSIGN_ERROR(status, err);
                }
            }
        }
    }

    /* Call nc4io_put_var for num times */
    for(i = 0; i < num; i++){
        err = nc4io_put_var(ncdp, varid, starts[i], counts[i], NULL, NULL, buf, bufcount, buftype, reqMode);
        if (err != NC_NOERR){ 
            if (isindep){
                DEBUG_RETURN_ERROR(err)
            }
            else{
                // Participate coll I/O
                if (status == NC_NOERR){
                    DEBUG_ASSIGN_ERROR(status, err);
                }
            }
        }
        
        /* Calculate the size of each put_var */
        putsize = (size_t)elsize;
        for(j = 0; j < ndim; j++){
            putsize *= counts[i][j];
        }

        /* Move buffer pointer */
        buf = (void*)(((char*)buf) + putsize);
    }

    if (!isindep){
        if (num_max != num_min){
            // switch back to coll mode
            err = nc_var_par_access(nc4p->ncid, varid, NC_COLLECTIVE);
            if (err != NC_NOERR){ 
                if (status == NC_NOERR){
                    DEBUG_ASSIGN_ERROR(status, err);
                }
            }
        }
    }

    return status;
}

int
nc4io_iget_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                void               *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc4io_get_varn */
    err = nc4io_get_varn(ncdp, varid, num, starts, counts, buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    /* TODO: Issue dummy id */
    *reqid = NC_REQ_NULL;

    return NC_NOERR;
}

int
nc4io_iput_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                const void         *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* We do not support nonblocking I/O so far */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    
    /* Call nc4io_put_varn */
    err = nc4io_put_varn(ncdp, varid, num, starts, counts, buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    /* TODO: Issue dummy id */
    *reqid = NC_REQ_NULL;

    return NC_NOERR;
}

int
nc4io_bput_varn(void               *ncdp,
                int                 varid,
                int                 num,
                MPI_Offset* const  *starts,
                MPI_Offset* const  *counts,
                const void         *buf,
                MPI_Offset          bufcount,
                MPI_Datatype        buftype,
                int                *reqid,
                int                 reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc4io_iput_varn */
    err = nc4io_iput_varn(ncdp, varid, num, starts, counts, buf, bufcount, buftype, reqid, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nc4io_get_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* vard not supported in NetCDF */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    return NC_NOERR;
}

int
nc4io_put_vard(void         *ncdp,
               int           varid,
               MPI_Datatype  filetype,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int           reqMode)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* vard not supported in NetCDF */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)

    return NC_NOERR;
}

