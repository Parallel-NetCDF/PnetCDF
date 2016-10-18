/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "nc.h"
#include "macro.h"
#include "ncmpidtype.h"

/*----< ncmpii_sanity_check() >----------------------------------------------*/
int ncmpii_sanity_check(int               ncid,
                        int               varid,
                        const MPI_Offset *start,
                        const MPI_Offset *count,
                        const MPI_Offset  bufcount,
                        MPI_Datatype      buftype,  /* internal datatype */
                        enum API_KIND     api,
                        int               isFlexibleAPI,
                        int               mustInDataMode,
                        int               rw_flag,
                        int               io_method,
                        NC              **ncp,   /* OUT */
                        NC_var          **varp)  /* OUT */
{
    /* all errors detected here are fatal, must return immediately */
    int err;

    /* check if ncid is valid */
    err = ncmpii_NC_check_id(ncid, ncp);
    if (err != NC_NOERR) return err;
    /* For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */

    /* if this call must be made in data mode, check if currently is in define
     * mode */
    if (mustInDataMode && NC_indef(*ncp)) {
        DEBUG_ASSIGN_ERROR(err, NC_EINDEFINE)
        goto fn_exit;
    }

    /* check file write permission if this is write request */
    if (rw_flag == WRITE_REQ && NC_readonly(*ncp)) {
        DEBUG_ASSIGN_ERROR(err, NC_EPERM)
        goto fn_exit;
    }

    if (io_method != NONBLOCKING_IO) { /* for blocking APIs */
        /* check if in the right collective or independent mode and initialize
         * MPI file handlers */
        err = ncmpii_check_mpifh(*ncp, io_method);
        if (err != NC_NOERR) goto fn_exit;
    }

    /* check if varid is valid */
    err = ncmpii_NC_lookupvar(*ncp, varid, varp);
    if (err != NC_NOERR) goto fn_exit;

    /* check NC_ECHAR */
    if (isFlexibleAPI) {
        /* when buftype == MPI_DATATYPE_NULL, bufcount is ignored and this API
         * assumes argument buf's data type matches the data type of variable
         * defined in the file - no data conversion will be done.
         */
        if (buftype != MPI_DATATYPE_NULL) {
            int isderived, el_size, buftype_is_contig;
            MPI_Datatype ptype;
            MPI_Offset   bnelems=0;

            err = ncmpii_dtype_decode(buftype, &ptype, &el_size, &bnelems,
                                      &isderived, &buftype_is_contig);
            if (err != NC_NOERR) goto fn_exit;

            err = NCMPII_ECHAR((*varp)->type, ptype);
            if (err != NC_NOERR) goto fn_exit;
        }
        /* else case types are matched */
    }
    else {
        err = NCMPII_ECHAR((*varp)->type, buftype);
        if (err != NC_NOERR) goto fn_exit;
    }

    /* for API var1, vara, vars, varm, and varn, start cannot be NULL */
    if (start == NULL && api >= API_VAR1 && (*varp)->ndims > 0) {
        DEBUG_ASSIGN_ERROR(err, NC_ENULLSTART)
        goto fn_exit;
    }

    /* for API vara, vars, and varm, count cannot be NULL */
    if (count == NULL && api >= API_VARA && (*varp)->ndims > 0) {
        DEBUG_ASSIGN_ERROR(err, NC_ENULLCOUNT)
        goto fn_exit;
    }

    /* for flexible APIs, bufcount cannot be negative */
    if (bufcount < 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto fn_exit;
    }

fn_exit:
    if ((*ncp)->safe_mode == 1 && io_method == COLL_IO) {
        int min_st, mpireturn;
        TRACE_COMM(MPI_Allreduce)(&err, &min_st, 1, MPI_INT, MPI_MIN, (*ncp)->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR) err = min_st;
    }
    return err;
}

