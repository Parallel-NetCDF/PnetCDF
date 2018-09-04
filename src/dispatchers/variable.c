/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h> /* INT_MAX */

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

/*----< ncmpi_def_var() >----------------------------------------------------*/
/* this API is collective, and must be called in define mode */
int
ncmpi_def_var(int         ncid,    /* IN:  file ID */
              const char *name,    /* IN:  name of variable */
              nc_type     type,
              int         ndims,
              const int  *dimids,
              int        *varidp)
{
    int i, err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (!(pncp->flag & NC_MODE_DEF)) { /* must be called in define mode */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

    if (name == NULL || *name == 0) { /* name cannot be NULL or NULL string */
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    if (strlen(name) > NC_MAX_NAME) { /* name length */
        DEBUG_ASSIGN_ERROR(err, NC_EMAXNAME)
        goto err_check;
    }

    /* check if the name string is legal for netcdf format */
    err = ncmpii_check_name(name, pncp->format);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR(err)
        goto err_check;
    }

    /* the max data type supported by CDF-5 is NC_UINT64 */
    if (type <= 0 || type > NC_UINT64) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADTYPE)
        goto err_check;
    }

    /* For NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, and
     * NC_FORMAT_NETCDF4_CLASSIC files, only classic types are allowed. */
    if (type > NC_DOUBLE) {
        if (pncp->format <= NC_FORMAT_64BIT_OFFSET) {
            DEBUG_ASSIGN_ERROR(err, NC_ESTRICTCDF2)
            goto err_check;
        }
        else if (pncp->format == NC_FORMAT_NETCDF4_CLASSIC) {
            DEBUG_ASSIGN_ERROR(err, NC_ESTRICTNC3)
            goto err_check;
        }
    }

    /* Argument ndims is of type "int". Its max value will be less than
     * INT_MAX. Thus if NC_MAX_VAR_DIMS == INT_MAX, then there is no need to
     * check whether or not ndims > NC_MAX_VAR_DIMS.
     *
     * When checking against NC_MAX_VAR_DIMS, because there is no error code
     * corresponding to this, we use NC_EMAXDIMS
     */
#if NC_MAX_VAR_DIMS < INT_MAX
    if (ndims > NC_MAX_VAR_DIMS) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXDIMS)
        goto err_check;
    }
#endif
    if (ndims < 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto err_check;
    }

    /* Note we no longer limit the number of variables, as CDF file formats
     * impose no such limit. Thus, the value of NC_MAX_VARS has been changed
     * to NC_MAX_INT, as argument nvars is of type signed int in API
     * ncmpi_inq_nvars()
     */
    if (pncp->nvars == NC_MAX_VARS) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXVARS)
        goto err_check;
    }

    /* check whether new name is already in use, for this API (def_var) the
     * name should NOT already exist */
    err = pncp->driver->inq_varid(pncp->ncp, name, NULL);
    if (err != NC_ENOTVAR) {
        DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
        goto err_check;
    }
    else
        err = NC_NOERR;

    /* check dimids[] */
    if (ndims > 0 && dimids == NULL) { /* for non-scalar variable */
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto err_check;
    }
    for (i=0; i<ndims; i++) {
        if (dimids[i] < 0 || pncp->ndims == 0 || dimids[i] >= pncp->ndims) {
            DEBUG_ASSIGN_ERROR(err, NC_EBADDIM)
            goto err_check;
        }
    }

err_check:
    if (pncp->flag & NC_MODE_SAFE) {
        int root_name_len, root_ndims, minE, rank, mpireturn;
        char *root_name=NULL;

        /* first check the error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;

        MPI_Comm_rank(pncp->comm, &rank);

        /* check if name is consistent among all processes */
        assert(name != NULL);
        root_name_len = strlen(name) + 1;
        TRACE_COMM(MPI_Bcast)(&root_name_len, 1, MPI_INT, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast root_name_len");

        root_name = (char*) NCI_Malloc((size_t)root_name_len);
        if (rank == 0) strcpy(root_name, name);
        TRACE_COMM(MPI_Bcast)(root_name, root_name_len, MPI_CHAR, 0,pncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(root_name);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && strcmp(root_name, name))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_NAME)
        NCI_Free(root_name);

        /* check if type is consistent among all processes */
        nc_type root_type=type;
        TRACE_COMM(MPI_Bcast)(&root_type, 1, MPI_INT, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_type != type)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_TYPE)

        /* check if ndims is consistent among all processes */
        root_ndims=ndims;
        TRACE_COMM(MPI_Bcast)(&root_ndims, 1, MPI_INT, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_ndims != ndims)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_NDIMS)

        /* check if dimids is consistent among all processes */
        if (root_ndims > 0) {
            int *root_dimids = (int*)NCI_Malloc((size_t)root_ndims *SIZEOF_INT);
            if (dimids != NULL)
                memcpy(root_dimids, dimids, (size_t)root_ndims*SIZEOF_INT);
            else
                memset(root_dimids, 0, (size_t)root_ndims*SIZEOF_INT);
            TRACE_COMM(MPI_Bcast)(root_dimids, root_ndims,MPI_INT,0,pncp->comm);
            if (mpireturn != MPI_SUCCESS) {
                NCI_Free(root_dimids);
                return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
            }
            if (err == NC_NOERR && dimids != NULL &&
                memcmp(root_dimids, dimids, (size_t)root_ndims*SIZEOF_INT))
                DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_DIMIDS)
            NCI_Free(root_dimids);
        }

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }

    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_def_var() */
    err = pncp->driver->def_var(pncp->ncp, name, type, ndims, dimids, varidp);
    if (err != NC_NOERR) return err;

    assert(*varidp == pncp->nvars);

    /* add new variable into pnc-vars[] */
    if (pncp->nvars % PNC_VARS_CHUNK == 0)
        pncp->vars = NCI_Realloc(pncp->vars,
                                 (pncp->nvars+PNC_VARS_CHUNK)*sizeof(PNC_var));

    pncp->vars[*varidp].ndims  = ndims;
    pncp->vars[*varidp].xtype  = type;
    pncp->vars[*varidp].recdim = -1;   /* if fixed-size variable */
    pncp->vars[*varidp].shape  = NULL;
    if (ndims > 0) {
        if (dimids[0] == pncp->unlimdimid) { /* record variable */
            pncp->vars[*varidp].recdim = pncp->unlimdimid;
            pncp->nrec_vars++;
        }

        pncp->vars[*varidp].shape = (MPI_Offset*)
                                    NCI_Malloc(ndims * SIZEOF_MPI_OFFSET);
        for (i=0; i<ndims; i++) {
            /* obtain size of dimension i */
            err = pncp->driver->inq_dim(pncp->ncp, dimids[i], NULL,
                                        pncp->vars[*varidp].shape+i);
            if (err != NC_NOERR) return err;
        }
    }
    pncp->nvars++;

    return NC_NOERR;
}

/*----< ncmpi_def_var_fill() >-----------------------------------------------*/
/* this API is collective, and must be called in define mode */
int
ncmpi_def_var_fill(int         ncid,    /* IN:  file ID */
                   int         varid,
                   int         nofill,
                   const void *fill_value)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (!(pncp->flag & NC_MODE_DEF)) { /* must be called in define mode */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

    if (varid == NC_GLOBAL) {
        /* setting global _FillValue through this API is not allowed */
        DEBUG_ASSIGN_ERROR(err, NC_EGLOBAL)
        goto err_check;
    }

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

err_check:
    if (pncp->flag & NC_MODE_SAFE) {
        int minE, mpireturn;
        /* check error code so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }

    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_def_var_fill() */
    return pncp->driver->def_var_fill(pncp->ncp, varid, nofill, fill_value);
}

/*----< ncmpi_inq_varid() >--------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_inq_varid(int         ncid,    /* IN:  file ID */
                const char *name,    /* IN:  name of variable */
                int        *varidp)  /* OUT: variable ID */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (name == NULL || *name == 0) DEBUG_RETURN_ERROR(NC_EBADNAME)

    if (strlen(name) > NC_MAX_NAME) DEBUG_RETURN_ERROR(NC_EMAXNAME)

    /* calling the subroutine that implements ncmpi_inq_varid() */
    return pncp->driver->inq_varid(pncp->ncp, name, varidp);
}

/*----< ncmpi_inq_var() >----------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_inq_var(int      ncid,    /* IN:  file ID */
              int      varid,   /* IN:  variable ID */
              char    *name,    /* OUT: name of variable */
              nc_type *xtypep,
              int     *ndimsp,
              int     *dimids,
              int     *nattsp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL &&
        (name != NULL || xtypep != NULL || ndimsp != NULL || dimids != NULL))
        DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars))
        DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* calling the subroutine that implements ncmpi_inq_var() */
    return pncp->driver->inq_var(pncp->ncp, varid, name, xtypep, ndimsp,
                                 dimids, nattsp, NULL, NULL, NULL);
}

/*----< ncmpi_inq_varname() >------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_inq_varname(int   ncid,    /* IN:  file ID */
                  int   varid,   /* IN:  variable ID */
                  char *name)    /* OUT: name of variable */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* calling the subroutine that implements ncmpi_inq_varname() */
    return pncp->driver->inq_var(pncp->ncp, varid, name, NULL, NULL,
                                 NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_vartype() >------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_inq_vartype(int      ncid,    /* IN:  file ID */
                  int      varid,   /* IN:  variable ID */
                  nc_type *xtypep)  /* OUT: external type of variable */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    *xtypep = pncp->vars[varid].xtype;
    return NC_NOERR;

#if 0
    /* calling the subroutine that implements ncmpi_inq_vartype() */
    return pncp->driver->inq_var(pncp->ncp, varid, NULL, xtypep, NULL,
                                 NULL, NULL, NULL, NULL, NULL);
#endif
}

/*----< ncmpi_inq_varndims() >-----------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_inq_varndims(int  ncid,    /* IN:  file ID */
                   int  varid,   /* IN:  variable ID */
                   int *ndimsp)  /* OUT: number of dimensions of variable */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    *ndimsp = pncp->vars[varid].ndims;
    return NC_NOERR;

#if 0
    /* calling the subroutine that implements ncmpi_inq_varndims() */
    return pncp->driver->inq_var(pncp->ncp, varid, NULL, NULL, ndimsp,
                                 NULL, NULL, NULL, NULL, NULL);
#endif
}

/*----< ncmpi_inq_vardimid() >-----------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_inq_vardimid(int  ncid,    /* IN:  file ID */
                   int  varid,   /* IN:  variable ID */
                   int *dimids)  /* OUT: dimension IDs of variable */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* calling the subroutine that implements ncmpi_inq_vardimid() */
    return pncp->driver->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
                                 dimids, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_varnatts() >-----------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_inq_varnatts(int  ncid,    /* IN:  file ID */
                   int  varid,   /* IN:  variable ID */
                   int *nattsp)  /* OUT: number of attributes of variable */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars))
         DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* calling the subroutine that implements ncmpi_inq_varnatts() */
    return pncp->driver->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
                                 NULL, nattsp, NULL, NULL, NULL);
}

/*----< ncmpi_inq_varoffset() >----------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_inq_varoffset(int         ncid,   /* IN: file ID */
                    int         varid,  /* IN: variable ID */
                    MPI_Offset *offset) /* OUT: starting file offset */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* calling the subroutine that implements ncmpi_inq_varoffset() */
    return pncp->driver->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
                                 NULL, NULL, offset, NULL, NULL);
}

/*----< ncmpi_inq_var_fill() >-----------------------------------------------*/
/* this API can be called independently and in both data and define mode */
int
ncmpi_inq_var_fill(int   ncid,
                   int   varid,
                   int  *no_fill,    /* OUT: 1 not fill mode, 0 fill mode */
                   void *fill_value) /* OUT: user-defined or default fill value */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* calling the subroutine that implements ncmpi_inq_var_fill() */
    return pncp->driver->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
                                 NULL, NULL, NULL, no_fill, fill_value);
}

/*----< ncmpi_fill_var_rec() >-----------------------------------------------*/
/* this API is collective and can only be called in collective data mode */
int
ncmpi_fill_var_rec(int        ncid,
                   int        varid,
                   MPI_Offset recno)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (pncp->flag & NC_MODE_RDONLY) { /* cannot be read-only */
        DEBUG_ASSIGN_ERROR(err, NC_EPERM)
        goto err_check;
    }

    if (pncp->flag & NC_MODE_DEF) { /* must be called in data mode */
        DEBUG_ASSIGN_ERROR(err, NC_EINDEFINE)
        goto err_check;
    }

    /* using NC_GLOBAL in varid is illegal for this API */
    if (varid == NC_GLOBAL) {
         DEBUG_ASSIGN_ERROR(err, NC_EGLOBAL)
        goto err_check;
    }

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    if (pncp->vars[varid].recdim == -1) { /* not a record variable */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTRECVAR)
        goto err_check;
    }

    if (pncp->flag & NC_MODE_INDEP) { /* must be called in collective mode */
        DEBUG_ASSIGN_ERROR(err, NC_EINDEP)
        goto err_check;
    }

err_check:
    if (pncp->flag & NC_MODE_SAFE) {
        int minE, mpireturn;
        /* check error code so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }

    /* calling the subroutine that implements ncmpi_fill_var_rec() */
    return pncp->driver->fill_var_rec(pncp->ncp, varid, recno);
}

/*----< ncmpi_rename_var() >-------------------------------------------------*/
/* This is a collective subroutine */
int
ncmpi_rename_var(int         ncid,    /* IN: file ID */
                 int         varid,   /* IN: variable ID */
                 const char *newname) /* IN: name of variable */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (pncp->flag & NC_MODE_RDONLY) { /* cannot be read-only */
        DEBUG_ASSIGN_ERROR(err, NC_EPERM)
        goto err_check;
    }

    if (varid == NC_GLOBAL) { /* Global is error in this context */
        DEBUG_ASSIGN_ERROR(err, NC_EGLOBAL)
        goto err_check;
    }

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    if (newname == NULL || *newname == 0) { /* cannot be NULL or NULL string */
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    if (strlen(newname) > NC_MAX_NAME) { /* newname length */
        DEBUG_ASSIGN_ERROR(err, NC_EMAXNAME)
        goto err_check;
    }

    /* check whether new name is legal */
    err = ncmpii_check_name(newname, pncp->format);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR(err)
        goto err_check;
    }

    /* check whether new name is already in use, for this API (rename) the
     * name should NOT already exist */
    err = pncp->driver->inq_varid(pncp->ncp, newname, NULL);
    if (err != NC_ENOTVAR) { /* expecting NC_ENOTVAR */
        DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
        goto err_check;
    }
    else err = NC_NOERR;  /* reset err */

err_check:
    if (pncp->flag & NC_MODE_SAFE) {
        int root_name_len, root_varid, minE, rank, mpireturn;
        char *root_name=NULL;

        /* First check error code so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;

        MPI_Comm_rank(pncp->comm, &rank);

        /* check if newname is consistent among all processes */
        assert(newname != NULL);
        root_name_len = strlen(newname) + 1;
        TRACE_COMM(MPI_Bcast)(&root_name_len, 1, MPI_INT, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast root_name_len");

        root_name = (char*) NCI_Malloc((size_t)root_name_len);
        if (rank == 0) strcpy(root_name, newname);
        TRACE_COMM(MPI_Bcast)(root_name, root_name_len, MPI_CHAR, 0,pncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(root_name);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && strcmp(root_name, newname))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_NAME)
        NCI_Free(root_name);

        /* check if varid is consistent across all processes */
        root_varid = varid;
        TRACE_COMM(MPI_Bcast)(&root_varid, 1, MPI_INT, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_varid != varid)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }

    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_rename_var() */
    return pncp->driver->rename_var(pncp->ncp, varid, newname);
}

