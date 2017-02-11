/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include <stdlib.h>

#include <dispatch.h>
#include <pnetcdf.h>

/*----< ncmpi_def_var() >----------------------------------------------------*/
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
    err = pncp->dispatch->def_var(pncp->ncp, name, type, ndims, dimids, varidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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
    err = pncp->dispatch->inq_varid(pncp->ncp, name, varidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_var() */
    err = pncp->dispatch->inq_var(pncp->ncp, varid, name, xtypep, ndimsp,
                                  dimids, nattsp, NULL, NULL, NULL);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_varname() */
    err = pncp->dispatch->inq_var(pncp->ncp, varid, name, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_vartype() */
    err = pncp->dispatch->inq_var(pncp->ncp, varid, NULL, xtypep, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_varndims() */
    err = pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, ndimsp,
                                  NULL, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_vardimid() */
    err = pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
                                  dimids, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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
    err = pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
                                  NULL, nattsp, NULL, NULL, NULL);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_varoffset() */
    err = pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
                                  NULL, NULL, offset, NULL, NULL);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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

    /* using NC_GLOBAL in varid is illegal for this API. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     */
    if (varid == NC_GLOBAL) return NC_EGLOBAL;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_var_fill() */
    err = pncp->dispatch->inq_var(pncp->ncp, varid, NULL, NULL, NULL,
                                  NULL, NULL, NULL, no_fill, fill_value);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
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
    err = pncp->dispatch->rename_var(pncp->ncp, varid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

