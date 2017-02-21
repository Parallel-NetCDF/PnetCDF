/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdlib.h>

#include <pnetcdf.h>
#include <dispatch.h>

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
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_def_var() */
    return pncp->dispatch->def_var(pncp->ncp, name, type, ndims, dimids, varidp);
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

    /* calling the subroutine that implements ncmpi_def_var_fill() */
    return pncp->dispatch->def_var_fill(pncp->ncp, varid, nofill, fill_value);
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

    /* calling the subroutine that implements ncmpi_inq_varid() */
    return pncp->dispatch->inq_varid(pncp->ncp, name, varidp);
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
        return NC_EGLOBAL;

    /* calling the subroutine that implements ncmpi_inq_var() */
    return pncp->dispatch->inq_var(pncp->ncp, varid, name, xtypep, ndimsp,
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
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* calling the subroutine that implements ncmpi_inq_varname() */
    return pncp->dispatch->inq_var(pncp->ncp, varid, name, NULL, NULL,
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
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* calling the subroutine that implements ncmpi_inq_vartype() */
    return pncp->dispatch->inq_var(pncp->ncp, varid, NULL, xtypep, NULL,
                                   NULL, NULL, NULL, NULL, NULL);
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
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* calling the subroutine that implements ncmpi_inq_varndims() */
    return pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, ndimsp,
                                   NULL, NULL, NULL, NULL, NULL);
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
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* calling the subroutine that implements ncmpi_inq_vardimid() */
    return pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
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

    /* calling the subroutine that implements ncmpi_inq_varnatts() */
    return pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
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
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* calling the subroutine that implements ncmpi_inq_varoffset() */
    return pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
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
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* calling the subroutine that implements ncmpi_inq_var_fill() */
    return pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
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

    /* using NC_GLOBAL in varid is illegal for this API */
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* calling the subroutine that implements ncmpi_fill_var_rec() */
    return pncp->dispatch->fill_rec(pncp->ncp, varid, recno);
}

/*----< ncmpi_rename_var() >-------------------------------------------------*/
/* This is an independent subroutine */
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

    /* calling the subroutine that implements ncmpi_rename_var() */
    return pncp->dispatch->rename_var(pncp->ncp, varid, newname);
}

/*----< ncmpi_cancel() >-----------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_cancel(int  ncid,
             int  num_reqs, /* number of requests */
             int *req_ids,  /* [num_reqs]: IN/OUT */
             int *statuses) /* [num_reqs], can be NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_cancel() */
    return pncp->dispatch->cancel(pncp->ncp, num_reqs, req_ids, statuses);
}
