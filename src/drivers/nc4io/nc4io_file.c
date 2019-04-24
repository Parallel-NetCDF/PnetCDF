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
#include <nc4io_driver.h>

int
nc4io_create(MPI_Comm     comm,
             const char  *path,
             int          cmode,
             int          ncid,
             MPI_Info     info,
             void       **ncpp)  /* OUT */
{
    char *filename;
    int err, ncidtmp;
    NC_nc4 *nc4p;

    /* remove the file system type prefix name if there is any.
     * For example, path=="lustre:/home/foo/testfile.nc",
     * use "/home/foo/testfile.nc" when calling nc_create_par()
     */
    filename = strchr(path, ':');
    if (filename == NULL) filename = (char*)path; /* no prefix */
    else                  filename++;

    /* add NC_MPIIO in case NetCDF 4.6.1 and earlier is used.
     * NC_MPIIO is ignored in 4.6.2 and after.
     */
    cmode |= NC_MPIIO;
    err = nc_create_par(filename, cmode, comm, info, &ncidtmp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Set fill mdoe to NC_NOFILL
     * netcdf default fill mode is NC_FILL */
    err = nc_set_fill(ncidtmp, NC_NOFILL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Create a NC_nc4 object and save its driver pointer */
    nc4p = (NC_nc4*) NCI_Malloc(sizeof(NC_nc4));
    if (nc4p == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nc4p->path = (char*) NCI_Malloc(strlen(filename)+1);
    if (nc4p->path == NULL) {
        NCI_Free(nc4p);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(nc4p->path, filename);
    nc4p->mode = cmode | NC_WRITE;
    nc4p->flag = NC_MODE_DEF;
    nc4p->ncid = ncid;
    nc4p->comm = comm;
    nc4p->ncid = ncidtmp;
    nc4p->putsize = 0;
    nc4p->getsize = 0;
    if (info == MPI_INFO_NULL)
        nc4p->mpiinfo = MPI_INFO_NULL;
    else
        MPI_Info_dup(info, &nc4p->mpiinfo);

    *ncpp = nc4p;

    return NC_NOERR;
}

int
nc4io_open(MPI_Comm     comm,
           const char  *path,
           int          omode,
           int          ncid,
           MPI_Info     info,
           void       **ncpp)
{
    char *filename;
    int err, ncidtmp;
    NC_nc4 *nc4p;

    /* remove the file system type prefix name if there is any.
     * For example, path=="lustre:/home/foo/testfile.nc",
     * use "/home/foo/testfile.nc" when calling nc_open_par()
     */
    filename = strchr(path, ':');
    if (filename == NULL) filename = (char*)path; /* no prefix */
    else                  filename++;

    /* add NC_MPIIO in case NetCDF 4.6.1 and earlier is used.
     * NC_MPIIO is ignored in 4.6.2 and after.
     */
    omode |= NC_MPIIO;
    err = nc_open_par(path, omode, comm, info, &ncidtmp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Create a NC_nc4 object and save its driver pointer */
    nc4p = (NC_nc4*) NCI_Malloc(sizeof(NC_nc4));
    if (nc4p == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nc4p->path = (char*) NCI_Malloc(strlen(filename)+1);
    if (nc4p->path == NULL) {
        NCI_Free(nc4p);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    strcpy(nc4p->path, filename);
    nc4p->mode = omode;
    nc4p->flag = 0;
    nc4p->ncid = ncid;
    nc4p->comm = comm;
    nc4p->ncid = ncidtmp;
    nc4p->putsize = 0;
    nc4p->getsize = 0;
    if (info == MPI_INFO_NULL)
        nc4p->mpiinfo = MPI_INFO_NULL;
    else
        MPI_Info_dup(info, &nc4p->mpiinfo);

    if (!fIsSet(omode, NC_WRITE)) fSet(nc4p->flag, NC_MODE_RDONLY);

    *ncpp = nc4p;

    return NC_NOERR;
}

int
nc4io_close(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    if (nc4p == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    /* Close with netcdf */
    err = nc_close(nc4p->ncid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    if (nc4p->mpiinfo != MPI_INFO_NULL)
        MPI_Info_free(&nc4p->mpiinfo);

    NCI_Free(nc4p->path);
    NCI_Free(nc4p);

    return err;
}

int
nc4io_enddef(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_enddef */
    err = nc_enddef(nc4p->ncid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Clear def flag */
    fClr(nc4p->flag, NC_MODE_DEF);

    return NC_NOERR;
}

int
nc4io__enddef(void       *ncdp,
              MPI_Offset  h_minfree,
              MPI_Offset  v_align,
              MPI_Offset  v_minfree,
              MPI_Offset  r_align)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    if (h_minfree != (size_t)h_minfree || v_align != (size_t)v_align ||
        v_minfree != (size_t)v_minfree || r_align != (size_t)r_align)
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    /* Call nc__enddef */
    err = nc__enddef(nc4p->ncid, (size_t)h_minfree, (size_t)v_align, (size_t)v_minfree, (size_t)r_align);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Clear def flag */
    fClr(nc4p->flag, NC_MODE_DEF);

    return NC_NOERR;
}

int
nc4io_redef(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Call nc_redef */
    err = nc_redef(nc4p->ncid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Set def flag */
    fSet(nc4p->flag, NC_MODE_DEF);

    return NC_NOERR;
}

int
nc4io_begin_indep_data(void *ncdp)
{
    int i, err, nvar;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Make sure we are in data mode */
    if (fIsSet(nc4p->flag, NC_MODE_DEF))
        DEBUG_RETURN_ERROR(NC_EINDEFINE);

    /* Get number of variables */
    err = nc_inq(nc4p->ncid, NULL, &nvar, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Set all variables to indep mode */
    for (i=0; i<nvar; i++) {
        err = nc_var_par_access(nc4p->ncid, i, NC_INDEPENDENT);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
    }

    /* Set indep flag */
    fSet(nc4p->flag, NC_MODE_INDEP);

    return NC_NOERR;
}

int
nc4io_end_indep_data(void *ncdp)
{
    int i, err, nvar;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Make sure we are in data mode */
    if (fIsSet(nc4p->flag, NC_MODE_DEF))
        DEBUG_RETURN_ERROR(NC_EINDEFINE);

    /* Get number of variables */
    err = nc_inq(nc4p->ncid, NULL, &nvar, NULL, NULL);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* Set all variables to coll mode */
    for (i=0; i<nvar; i++) {
        err = nc_var_par_access(nc4p->ncid, i, NC_COLLECTIVE);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
    }

    /* Clear indep flag */
    fClr(nc4p->flag, NC_MODE_INDEP);

    return NC_NOERR;
}

int
nc4io_abort(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    if (nc4p == NULL) DEBUG_RETURN_ERROR(NC_EBADID)

    /* Call nc_abort */
    err = nc_abort(nc4p->ncid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    if (nc4p->mpiinfo != MPI_INFO_NULL)
        MPI_Info_free(&nc4p->mpiinfo);

    NCI_Free(nc4p->path);
    NCI_Free(nc4p);

    return err;
}

int
nc4io_inq(void *ncdp,
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
nc4io_inq_misc(void       *ncdp,
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
    int i, j, err;
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

    if (info_used != NULL) MPI_Info_dup(nc4p->mpiinfo, info_used);

    /* Calculate number of record and fix sized variables
     * NOTE: We assume there are only 1 record dims in NetCDF 4 classic model
     * NOTE: In NetCDF4 enhenced model, there can be multiple record dims, only
     * the first one is considered as record dim here
     */
    if (num_fix_varsp != NULL || num_rec_varsp != NULL || recsize != NULL) {
        int nvar, ndim, *dims, *vars, nrecdim, *isrec, *recdims, nrec=0, nfix=0;

        /* Record dimid (TODO: non-classic model NetCDF4 file may have more
         * than one unlimited dimension) 
         * We flag each dimension that is reported as unlimited dimension
         */
        err = nc_inq(nc4p->ncid, &ndim, NULL, NULL, NULL);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);
        isrec = (int*)NCI_Malloc(SIZEOF_INT * ndim);
        memset(isrec, 0, SIZEOF_INT * ndim);

        /* Get list of unlimited dim */
        err = nc_inq_unlimdims(nc4p->ncid, &nrecdim, NULL);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
        recdims = (int*)NCI_Malloc(SIZEOF_INT * nrecdim);
        err = nc_inq_unlimdims(nc4p->ncid, NULL, recdims);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

        /* Set up the flag */
        for(i = 0; i < nrecdim; i++){
            isrec[recdims[i]] = 1;
        }

        /* number of variables */
        err = nc_inq_varids(nc4p->ncid, &nvar, NULL);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

        vars = NCI_Malloc(SIZEOF_INT * nvar);
        if (vars == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

        /* IDs of all variables */
        err = nc_inq_varids(nc4p->ncid, NULL, vars);
        if (err != NC_NOERR) {
            NCI_Free(vars);
            NCI_Free(recdims);
            NCI_Free(isrec);
            DEBUG_RETURN_ERROR(err)
        }

        if (recsize != NULL) *recsize = 0;

        /* Iterate through all variables */
        for (i=0; i<nvar; i++) {
            /* number of dimensions */
            err = nc_inq_varndims(nc4p->ncid, vars[i], &ndim);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

            dims = NCI_Malloc(SIZEOF_INT * ndim);
            if (dims == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM);

            /* IDs of dimensions */
            err = nc_inq_vardimid(nc4p->ncid, vars[i], dims);
            if (err != NC_NOERR) {
                NCI_Free(vars);
                NCI_Free(dims);
                NCI_Free(recdims);
                NCI_Free(isrec);
                DEBUG_RETURN_ERROR(err)
            }

            /* Iterate through all dimensions */
            for (j=0; j<ndim; j++)
                if (isrec[dims[j]])
                    break;
            /* If none of the dimension is record dim, count as fixed var */
            if (j == ndim) nfix += 1;
            else           nrec += 1;

            if (recsize == NULL || j == ndim) {
                NCI_Free(dims);
                continue;
            }

            /* calculate sum of all individual record sizes */
            nc_type xtype;
            int xtype_size;
            MPI_Offset var_size;

            /* variable external data type */
            err = nc_inq_vartype(nc4p->ncid, vars[i], &xtype);
            if (err != NC_NOERR) {
                NCI_Free(vars);
                NCI_Free(dims);
                NCI_Free(recdims);
                NCI_Free(isrec);
                DEBUG_RETURN_ERROR(err)
            }
            /* size of variable external data type */
            err = ncmpii_xlen_nc_type(xtype, &xtype_size);
            if (err != NC_NOERR) {
                NCI_Free(vars);
                NCI_Free(dims);
                NCI_Free(recdims);
                NCI_Free(isrec);
                DEBUG_RETURN_ERROR(err)
            }
            var_size = xtype_size;
            for (j=0; j<ndim; j++) {
                size_t dim_size;
                if (isrec[dims[j]]) continue;
                err = nc_inq_dimlen(nc4p->ncid, dims[j], &dim_size);
                if (err != NC_NOERR) {
                    NCI_Free(vars);
                    NCI_Free(dims);
                    NCI_Free(recdims);
                    NCI_Free(isrec);
                    DEBUG_RETURN_ERROR(err)
                }
                var_size *= dim_size;
            }
            *recsize += var_size;

            NCI_Free(dims);
        }
        NCI_Free(vars);
        NCI_Free(recdims);
        NCI_Free(isrec);

        if (num_fix_varsp != NULL) *num_fix_varsp = nfix;
        if (num_rec_varsp != NULL) *num_rec_varsp = nrec;
    }

    /* NetCDF does not expose any MPI related info */
    if (striping_size  != NULL) DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    if (striping_count != NULL) DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    /* Put size counted by the driver */
    if (put_size != NULL){
        *put_size = nc4p->putsize;
    }

    /* Get size counted by the driver */
    if (get_size != NULL){
        *get_size = nc4p->getsize;
    }

    /* NetCDF does not expose such info */
    if (header_size   != NULL) DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    if (header_extent != NULL) DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    /* TODO: suporting nonblocking IO */
    if (nreqs != NULL) DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    /* bput APIs are not supported yet */
    if (usage != NULL) DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    /* bput APIs are not supported yet */
    if (buf_size != NULL) DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
nc4io_cancel(void *ncdp,
             int   num_req,
             int  *req_ids,
             int  *statuses)
{
    /* We do not support nonblocking I/O yet */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
nc4io_wait(void *ncdp,
           int   num_reqs,
           int  *req_ids,
           int  *statuses,
           int   reqMode)
{
    /* We do not support nonblocking I/O yet */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
nc4io_set_fill(void *ncdp,
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
nc4io_fill_var_rec(void      *ncdp,
                   int        varid,
                   MPI_Offset recno)
{
    /* NetCDF does not support this natively */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
}

int
nc4io_def_var_fill(void       *ncdp,
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
nc4io_sync_numrecs(void *ncdp)
{
    /* Will NetCDF take care of this internally? */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
}

int
nc4io_sync(void *ncdp)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;

    /* Make sure we are in data mode */
    if (fIsSet(nc4p->flag, NC_MODE_DEF))
        DEBUG_RETURN_ERROR(NC_EINDEFINE);

    /* Call nc_sync */
    err = nc_sync(nc4p->ncid);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_flush(void *ncdp)
{
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT)
}

