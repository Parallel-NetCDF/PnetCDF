/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include "ncconfig.h"
#endif

#include <mpi.h>
#include "nc.h"

/*----< ncmpii_sanity_check() >----------------------------------------------*/
int ncmpii_sanity_check(int               ncid,
                        int               varid,
                        const MPI_Offset *start,
                        const MPI_Offset *count,
                        MPI_Offset        bufcount,
                        enum API_KIND     api,
                        int               isFlexAPI,
                        int               rw_flag,
                        int               io_method,
                        NC              **ncp,
                        NC_var          **varp)
{
    /* all errors detected here are fatal, must return immediately */
    int status;

    /* check if ncid is valid */
    status = ncmpii_NC_check_id(ncid, ncp);
    if (status != NC_NOERR) return status;
    /* For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */

    /* in this context, define mode is wrong, must be called in data mode */
    if (NC_indef(*ncp)) {
        status = NC_EINDEFINE;
        goto fn_exit;
    }

    /* check file write permission if this is write request */
    if (rw_flag == WRITE_REQ && NC_readonly(*ncp)) {
        status = NC_EPERM;
        goto fn_exit;
    }

    if (io_method != NONBLOCKING_IO) { /* for blocking APIs */
        /* check if in the right collective or independent mode and initialize
         * MPI file handlers */
        status = ncmpii_check_mpifh(*ncp, io_method);
        if (status != NC_NOERR) goto fn_exit;
    }

    /* check if varid is valid */
    status = ncmpii_NC_lookupvar(*ncp, varid, varp);
    if (status != NC_NOERR) goto fn_exit;

    /* for API var1, vara, vars, varm, and varn, start cannot be NULL */
    if (start == NULL && api >= API_VAR1 && (*varp)->ndims > 0) {
        status = NC_ENULLSTART;
        goto fn_exit;
    }

    /* for API vara, vars, and varm, count cannot be NULL */
    if (count == NULL && api >= API_VARA && (*varp)->ndims > 0) {
        status = NC_ENULLCOUNT;
        goto fn_exit;
    }

    /* for flexible APIs, bufcount cannot be negative */
    if (isFlexAPI && bufcount < 0) {
        status = NC_EINVAL;
        goto fn_exit;
    }

fn_exit:
    if ((*ncp)->safe_mode == 1 && io_method == COLL_IO) {
        int min_st;
        MPI_Allreduce(&status, &min_st, 1, MPI_INT, MPI_MIN, (*ncp)->nciop->comm);
        if (status == NC_NOERR) status = min_st;
    }
    return status;
}

