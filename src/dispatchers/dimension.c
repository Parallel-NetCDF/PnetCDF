/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include <stdlib.h>

#include <dispatch.h>
#include <pnetcdf.h>

/*----< ncmpi_def_dim() >----------------------------------------------------*/
int
ncmpi_def_dim(int         ncid,    /* IN:  file ID */
              const char *name,    /* IN:  name of dimension */
              MPI_Offset  size,    /* IN:  dimension size */
              int        *dimidp)  /* OUT: dimension ID */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_def_dim() */
    err = pncp->dispatch->def_dim(pncp->ncp, name, size, dimidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_dimid() >--------------------------------------------------*/
int
ncmpi_inq_dimid(int         ncid,    /* IN:  file ID */
                const char *name,    /* IN:  name of dimension */
                int        *dimidp)  /* OUT: dimension ID */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_dimid() */
    err = pncp->dispatch->inq_dimid(pncp->ncp, name, dimidp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_dim() >----------------------------------------------------*/
int
ncmpi_inq_dim(int         ncid,    /* IN:  file ID */
              int         dimid,   /* IN:  dimension ID */
              char       *name,    /* OUT: name of dimension */
              MPI_Offset *lengthp) /* OUT: length of dimension */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_dim() */
    err = pncp->dispatch->inq_dim(pncp->ncp, dimid, name, lengthp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_dimname() >------------------------------------------------*/
int
ncmpi_inq_dimname(int   ncid,    /* IN:  file ID */
                  int   dimid,   /* IN:  dimension ID */
                  char *name)    /* OUT: name of dimension */
{
    return ncmpi_inq_dim(ncid, dimid, name, NULL);
}

/*----< ncmpi_inq_dimlen() >-------------------------------------------------*/
int
ncmpi_inq_dimlen(int         ncid,
                 int         dimid,
                 MPI_Offset *lenp)
{
    return ncmpi_inq_dim(ncid, dimid, NULL, lenp);
}

/*----< ncmpi_rename_dim() >-------------------------------------------------*/
int
ncmpi_rename_dim(int         ncid,    /* IN: file ID */
                 int         dimid,   /* IN: dimension ID */
                 const char *newname) /* IN: name of dimension */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_rename_dim() */
    err = pncp->dispatch->rename_dim(pncp->ncp, dimid, newname);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

