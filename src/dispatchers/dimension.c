/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include <stdlib.h>

#include <dispatch.h>
#include <pnetcdf.h>

/*----< ncmpi_def_dim() >----------------------------------------------------*/
/* This is a collective subroutine. */
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
    return pncp->dispatch->def_dim(pncp->ncp, name, size, dimidp);
}

/*----< ncmpi_inq_dimid() >--------------------------------------------------*/
/* This is an independent subroutine. */
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
    return pncp->dispatch->inq_dimid(pncp->ncp, name, dimidp);
}

/*----< ncmpi_inq_dim() >----------------------------------------------------*/
/* This is an independent subroutine. */
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
    return pncp->dispatch->inq_dim(pncp->ncp, dimid, name, lengthp);
}

/*----< ncmpi_inq_dimname() >------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_dimname(int   ncid,    /* IN:  file ID */
                  int   dimid,   /* IN:  dimension ID */
                  char *name)    /* OUT: name of dimension */
{
    return ncmpi_inq_dim(ncid, dimid, name, NULL);
}

/*----< ncmpi_inq_dimlen() >-------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_dimlen(int         ncid,
                 int         dimid,
                 MPI_Offset *lenp)
{
    return ncmpi_inq_dim(ncid, dimid, NULL, lenp);
}

/*----< ncmpi_rename_dim() >-------------------------------------------------*/
/* This is a collective subroutine. */
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
    return pncp->dispatch->rename_dim(pncp->ncp, dimid, newname);
}

