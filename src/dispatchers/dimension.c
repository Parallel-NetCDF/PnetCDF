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

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

/*----< ncmpi_def_dim() >----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_def_dim(int         ncid,    /* IN:  file ID */
              const char *name,    /* IN:  name of dimension */
              MPI_Offset  size,    /* IN:  dimension size */
              int        *dimidp)  /* OUT: dimension ID */
{
    int err=NC_NOERR, dimid;
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

    /* check if the name string is legal for the netcdf format */
    err = ncmpii_check_name(name, pncp->format);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR(err)
        goto err_check;
    }

    /* MPI_Offset is usually a signed value, but serial netcdf uses size_t.
     * In 1999 ISO C standard, size_t is an unsigned integer type of at least
     * 16 bit. */
    if (pncp->format == NC_FORMAT_CDF2) { /* CDF-2 format, max is 2^32-4 */
        if (size > NC_MAX_UINT - 3 || (size < 0))
            /* "-3" handles rounded-up size */
            err = NC_EDIMSIZE;
    } else if (pncp->format == NC_FORMAT_CDF5) { /* CDF-5 format */
        if (size < 0)
            err = NC_EDIMSIZE;
    } else if (pncp->format == NC_FORMAT_NETCDF4 ||
               pncp->format == NC_FORMAT_NETCDF4_CLASSIC) { /* NetCDF-4 format */
        if (size < 0)
            err = NC_EDIMSIZE;
    } else { /* CDF-1 format, max is 2^31-4 */
        if (size > NC_MAX_INT - 3 || (size < 0))
            /* "-3" handles rounded-up size */
            err = NC_EDIMSIZE;
    }
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR(err)
        goto err_check;
    }

    if (size == NC_UNLIMITED && pncp->unlimdimid != -1) {
        /* netcdf allows one unlimited dimension defined per file */
        DEBUG_ASSIGN_ERROR(err, NC_EUNLIMIT) /* already defined */
        goto err_check;
    }

    /* Note we no longer limit the number of dimensions, as CDF file formats
     * impose no such limit. Thus, the value of NC_MAX_DIMS has been changed
     * to NC_MAX_INT, as argument ndims in ncmpi_inq_varndims() is of type
     * signed int.
     */
    if (pncp->ndims == NC_MAX_DIMS) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXDIMS)
        goto err_check;
    }

    /* check if the name string is previously used */
    err = pncp->driver->inq_dimid(pncp->ncp, name, NULL);
    if (err != NC_EBADDIM) {
        DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
        goto err_check;
    }
    else err = NC_NOERR;

err_check:
    if (pncp->flag & NC_MODE_SAFE) {
        int root_name_len, minE, rank, mpireturn;
        char *root_name=NULL;
        MPI_Offset root_size;

        /* check the error so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR)
            return minE;

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
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_NAME)
        NCI_Free(root_name);

        /* check if sizes are consistent across all processes */
        root_size = size;
        TRACE_COMM(MPI_Bcast)(&root_size, 1, MPI_OFFSET, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_size != size)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_SIZE)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR)
            return minE;
    }

    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_def_dim() */
    err = pncp->driver->def_dim(pncp->ncp, name, size, &dimid);
    if (err != NC_NOERR) return err;

    if (size == NC_UNLIMITED && pncp->unlimdimid == -1)
        pncp->unlimdimid = dimid;

    pncp->ndims++;

    if (dimidp != NULL) *dimidp = dimid;

    return NC_NOERR;
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

    if (name == NULL || *name == 0) DEBUG_RETURN_ERROR(NC_EBADNAME)

    if (strlen(name) > NC_MAX_NAME) DEBUG_RETURN_ERROR(NC_EMAXNAME)

    /* calling the subroutine that implements ncmpi_inq_dimid() */
    return pncp->driver->inq_dimid(pncp->ncp, name, dimidp);
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

    if (dimid < 0 || dimid >= pncp->ndims) DEBUG_RETURN_ERROR(NC_EBADDIM)

    /* calling the subroutine that implements ncmpi_inq_dim() */
    return pncp->driver->inq_dim(pncp->ncp, dimid, name, lengthp);
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
    int err=NC_NOERR, inq_id, skip_rename=0;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (pncp->flag & NC_MODE_RDONLY) { /* cannot be read-only */
        DEBUG_ASSIGN_ERROR(err, NC_EPERM)
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

    /* check if the newname string is legal for the netcdf format */
    err = ncmpii_check_name(newname, pncp->format);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR(err)
        goto err_check;
    }

    /* check NC_EBADDIM for whether dimid is valid */
    if (dimid < 0 || dimid >= pncp->ndims) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADDIM)
        goto err_check;
    }

    /* check if the name string is previously used */
    err = pncp->driver->inq_dimid(pncp->ncp, newname, &inq_id);
    if (err == NC_NOERR) { /* name already exist */
        if (inq_id == dimid) /* same name, same dimid, skip rename */
            skip_rename = 1;
        else
            DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
    }
    else if (err == NC_EBADDIM) { /* cannot find dimid with this name */
        err = NC_NOERR;
    }

err_check:
    if (pncp->flag & NC_MODE_SAFE) {
        int root_name_len, root_dimid, minE, rank, mpireturn;
        char *root_name=NULL;

        /* check the error so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;

        MPI_Comm_rank(pncp->comm, &rank);

        /* check if name is consistent among all processes */
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
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_NAME)
        NCI_Free(root_name);

        /* check if dimid is consistent across all processes */
        root_dimid = dimid;
        TRACE_COMM(MPI_Bcast)(&root_dimid, 1, MPI_INT, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_dimid != dimid)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }

    if (err != NC_NOERR) return err;

    if (skip_rename) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_rename_dim() */
    return pncp->driver->rename_dim(pncp->ncp, dimid, newname);
}

