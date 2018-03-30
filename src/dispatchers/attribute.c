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

/*----< ncmpi_inq_att() >----------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_att(int         ncid,
              int         varid,
              const char *name, /* input, attribute name */
              nc_type    *xtypep,
              MPI_Offset *lenp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars))
        DEBUG_RETURN_ERROR(NC_ENOTVAR)

    if (name == NULL || *name == 0) DEBUG_RETURN_ERROR(NC_EBADNAME)

    if (strlen(name) > NC_MAX_NAME) DEBUG_RETURN_ERROR(NC_EMAXNAME)

    /* calling the subroutine that implements ncmpi_inq_att() */
    return pncp->driver->inq_att(pncp->ncp, varid, name, xtypep, lenp);
}

/*----< ncmpi_inq_atttype() >------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_atttype(int         ncid,
                  int         varid,
                  const char *name, /* input, attribute name */
                  nc_type    *xtypep)
{
    return ncmpi_inq_att(ncid, varid, name, xtypep, NULL);
}

/*----< ncmpi_inq_attlen() >-------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_attlen(int         ncid,
                 int         varid,
                 const char *name, /* input, attribute name */
                 MPI_Offset *lenp)
{
    return ncmpi_inq_att(ncid, varid, name, NULL, lenp);
}

/*----< ncmpi_inq_attid() >--------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_attid(int         ncid,
                int         varid,
                const char *name,
                int        *attnump)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars))
        DEBUG_RETURN_ERROR(NC_ENOTVAR)

    if (name == NULL || *name == 0) DEBUG_RETURN_ERROR(NC_EBADNAME)

    if (strlen(name) > NC_MAX_NAME) DEBUG_RETURN_ERROR(NC_EMAXNAME)

    /* calling the subroutine that implements ncmpi_inq_attid() */
    return pncp->driver->inq_attid(pncp->ncp, varid, name, attnump);
}

/*----< ncmpi_inq_attname() >------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_attname(int   ncid,
                  int   varid,
                  int   attnum,
                  char *name) /* output, attribute name */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars))
        DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* calling the subroutine that implements ncmpi_inq_attname() */
    return pncp->driver->inq_attname(pncp->ncp, varid, attnum, name);
}

/*----< ncmpi_copy_att() >---------------------------------------------------*/
/* This is a collective subroutine.
 * ncid_out must be in define mode. If varid_in's attribute name has alreay
 * existed in varid_out, it means to overwrite the attribute in varid_out.
 * In this case, if the space used by varid_in's attribute is larger than
 * varid_out's, then this API must be called when the file is in define mode.
 */
int
ncmpi_copy_att(int         ncid_in,
               int         varid_in,
               const char *name,
               int         ncid_out,
               int         varid_out)
{
    int err;
    PNC *pncp_in, *pncp_out;

    /* check if ncid_in is valid */
    err = PNC_check_id(ncid_in, &pncp_in);
    if (err != NC_NOERR) return err;

    /* check if ncid_out is valid */
    err = PNC_check_id(ncid_out, &pncp_out);
    if (err != NC_NOERR) return err;

    if (pncp_out->flag & NC_MODE_RDONLY) { /* cannot be read-only */
        DEBUG_ASSIGN_ERROR(err, NC_EPERM)
        goto err_check;
    }

    /* check whether variable ID is valid */
    if (varid_in != NC_GLOBAL &&
        (varid_in < 0 || varid_in >= pncp_in->nvars)) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    /* check whether variable ID is valid */
    if (varid_out != NC_GLOBAL &&
        (varid_out < 0 || varid_out >= pncp_out->nvars)) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    if (name == NULL || *name == 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    if (strlen(name) > NC_MAX_NAME) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXNAME)
        goto err_check;
    }

err_check:
    if (pncp_out->flag & NC_MODE_SAFE) {
        int root_ids[2], root_name_len, minE, rank, mpireturn;
        char *root_name=NULL;

        /* check the error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN,
                                  pncp_out->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;

        MPI_Comm_rank(pncp_out->comm, &rank);

        /* check if name is consistent among all processes */
        assert(name != NULL);
        root_name_len = strlen(name) + 1;
        TRACE_COMM(MPI_Bcast)(&root_name_len, 1, MPI_INT, 0, pncp_out->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast root_name_len");

        root_name = (char*) NCI_Malloc((size_t)root_name_len);
        if (rank == 0) strcpy(root_name, name);
        TRACE_COMM(MPI_Bcast)(root_name, root_name_len, MPI_CHAR, 0,
                              pncp_out->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(root_name);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && strcmp(root_name, name))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_NAME)
        NCI_Free(root_name);

        /* check if varid_in, varid_out, are consistent across all
         * processes */
        root_ids[0] = varid_in;
        root_ids[1] = varid_out;
        TRACE_COMM(MPI_Bcast)(&root_ids, 2, MPI_INT, 0, pncp_out->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && (root_ids[0] != varid_in ||
            root_ids[1] != varid_out))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN,
                                  pncp_out->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }

    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_copy_att() */
    return pncp_in->driver->copy_att(pncp_in->ncp,  varid_in, name,
                                     pncp_out->ncp, varid_out);
}

/*----< ncmpi_rename_att() >-------------------------------------------------*/
/* This is a collective subroutine. If the new name is longer than the old
 * name, this API must be called in define mode.
 */
int
ncmpi_rename_att(int         ncid,
                 int         varid,
                 const char *name,
                 const char *newname)
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

    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars)) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    if (name == NULL || *name == 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    if (strlen(name) > NC_MAX_NAME) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXNAME)
        goto err_check;
    }

    if (newname == NULL || *newname == 0) {
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
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_NAME)
        NCI_Free(root_name);

        /* check if newname is consistent among all processes */
        assert(newname != NULL);
        root_name_len = strlen(newname) + 1;
        TRACE_COMM(MPI_Bcast)(&root_name_len, 1, MPI_INT, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast root_name_len");

        root_name = (char*) NCI_Malloc((size_t)root_name_len);
        strcpy(root_name, newname);
        TRACE_COMM(MPI_Bcast)(root_name, root_name_len, MPI_CHAR, 0,pncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(root_name);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && strcmp(root_name, newname))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_NAME)
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

    /* calling the subroutine that implements ncmpi_rename_att() */
    return pncp->driver->rename_att(pncp->ncp, varid, name, newname);
}

/*----< ncmpi_del_att() >----------------------------------------------------*/
/* This is a collective subroutine.
 * This API must be called in define mode.
 */
int
ncmpi_del_att(int         ncid,
              int         varid,
              const char *name)
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

    if (!(pncp->flag & NC_MODE_DEF)) { /* must be called in define mode */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars))
        DEBUG_RETURN_ERROR(NC_ENOTVAR)

    if (name == NULL || *name == 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    if (strlen(name) > NC_MAX_NAME) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXNAME)
        goto err_check;
    }

err_check:
    if (pncp->flag & NC_MODE_SAFE) {
        int root_varid, root_name_len, minE, rank, mpireturn;
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
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_NAME)
        NCI_Free(root_name);

        /* check if varid is consistent across all processes */
        root_varid = varid;
        TRACE_COMM(MPI_Bcast)(&root_varid, 1, MPI_INT, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_varid != varid)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN,pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }

    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_del_att() */
    return pncp->driver->del_att(pncp->ncp, varid, name);
}

