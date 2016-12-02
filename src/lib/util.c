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
#include "ncx.h"

/*----< ncmpii_sanity_check() >----------------------------------------------*/
/* check the following errors and in that precedence.
 * NC_EBADID, NC_EPERM, NC_EINDEFINE, NC_EINDEP/NC_ENOTINDEP, NC_ENOTVAR,
 * NC_ECHAR, NC_EINVALCOORDS, NC_EEDGE, NC_ESTRIDE, NC_EINVAL.
 */
int ncmpii_sanity_check(int               ncid,
                        int               varid,
                        const MPI_Offset *start,
                        const MPI_Offset *count,
                        const MPI_Offset *stride,
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
    int i, firstDim, err;

    /* check if ncid is valid (check NC_EBADID) */
    err = ncmpii_NC_check_id(ncid, ncp);
    if (err != NC_NOERR) return err;
    /* For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */

    /* check file write permission if this is write request */
    if (rw_flag == WRITE_REQ && NC_readonly(*ncp)) {
        DEBUG_ASSIGN_ERROR(err, NC_EPERM)
        goto fn_exit;
    }

    /* if this call must be made in data mode, check if currently is in define
     * mode */
    if (mustInDataMode && NC_indef(*ncp)) {
        DEBUG_ASSIGN_ERROR(err, NC_EINDEFINE)
        goto fn_exit;
    }

    if (io_method != NONBLOCKING_IO) { /* for blocking APIs */
        /* check if in the right collective or independent mode and initialize
         * MPI file handlers */
        err = ncmpii_check_mpifh(*ncp, io_method);
        if (err != NC_NOERR) goto fn_exit;
    }

    /* check if varid is valid (check NC_ENOTVAR) */
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

    /* for flexible APIs, bufcount cannot be negative */
    if (bufcount < 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto fn_exit;
    }

    if ((*varp)->ndims == 0) { /* scalar variable: ignore start/count/stride */
        err = NC_NOERR;
        goto fn_exit;
    }

    if (api <= API_VAR) { /* var/varn/vard APIs, start/count/stride are NULL */
        err = NC_NOERR;
        goto fn_exit;
    }

    /* Now only check var1, vara, vars, and varm APIs */

    /* Check NC_EINVALCOORDS
     * for API var1/vara/vars/varm, start cannot be NULL, except for scalars */
    if (start == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
        goto fn_exit;
    }

    firstDim = 0;
    if (IS_RECVAR(*varp)) {
        if (start[0] < 0) { /* no negative value */
            DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
            goto fn_exit;
        }

        if ((*ncp)->format < 5 &&    /* not CDF-5 */
            start[0] > X_UINT_MAX) { /* sanity check */
            DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
            goto fn_exit;
        }

        /* for record variable, [0] is the NC_UNLIMITED dimension */
        if (rw_flag == READ_REQ) {
            /* read cannot go beyond current numrecs */
#ifdef RELAX_COORD_BOUND
            if (start[0] >  (*ncp)->numrecs) {
                DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
                goto fn_exit;
            }
            if (start[0] == (*ncp)->numrecs) {
                if (api == API_VAR1) {
                    /* for var1 APIs, count[0] is considered of 1 */
                    DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
                    goto fn_exit;
                }
                else if (count != NULL && count[0] > 0) {
                    DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
                    goto fn_exit;
                }
            }
#else
            if (start[0] >= (*ncp)->numrecs) {
                DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
                goto fn_exit;
            }
#endif          
        }
        firstDim = 1; /* skip checking the record dimension */
    }

    for (i=firstDim; i<(*varp)->ndims; i++) {
#ifdef RELAX_COORD_BOUND
        if (start[i] < 0 || start[i] > (*varp)->shape[i]) {
            DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
            goto fn_exit;
        }
        if (start[i] == (*ncp)->numrecs) {
            if (api == API_VAR1) {
                /* for var1 APIs, count[0] is considered of 1 */
                DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
                goto fn_exit;
            }
            else if (count != NULL && count[i] > 0) {
                DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
                goto fn_exit;
            }
        }
#else
        if (start[i] < 0 || start[i] >= (*varp)->shape[i]) {
            DEBUG_ASSIGN_ERROR(err, NC_EINVALCOORDS)
            goto fn_exit;
        }
#endif
    }

    if (api <= API_VAR1) {
        /* var1/var APIs have no count argument */
        err = NC_NOERR;
        goto fn_exit;
    }

    /* Check NC_EEDGE
     * for API vara/vars/varm, count cannot be NULL, except for scalars */
    if (count == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_EEDGE)
        goto fn_exit;
    }
    firstDim = 0;
    if (IS_RECVAR(*varp)) {
        if (count[0] < 0) { /* no negative count[] */
            DEBUG_ASSIGN_ERROR(err, NC_ENEGATIVECNT)
            goto fn_exit;
        }
        /* for record variable, [0] is the NC_UNLIMITED dimension */
        if (rw_flag == READ_REQ) { /* read cannot go beyond current numrecs */
            if (stride == NULL) {  /* for vara APIs */
                if (start[0] + count[0] > (*ncp)->numrecs) {
                    DEBUG_ASSIGN_ERROR(err, NC_EEDGE)
                    goto fn_exit;
                }
            }
            else { /* for vars/varm APIs */
                if (count[0] > 0 &&
                    start[0] + (count[0]-1) * stride[0] >= (*ncp)->numrecs) {
                    DEBUG_ASSIGN_ERROR(err, NC_EEDGE)
                    goto fn_exit;
                }
            }
        }
        firstDim = 1; /* skip checking the record dimension */
    }

    for (i=firstDim; i<(*varp)->ndims; i++) {
        if ((*varp)->shape[i] < 0) {
            DEBUG_ASSIGN_ERROR(err, NC_EEDGE)
            goto fn_exit;
        }
        if (count[i] < 0) { /* no negative count[] */
            DEBUG_ASSIGN_ERROR(err, NC_ENEGATIVECNT)
            goto fn_exit;
        }

        if (stride == NULL) { /* for vara APIs */
            if (count[i] > (*varp)->shape[i] ||
                start[i] + count[i] > (*varp)->shape[i]) {
                DEBUG_ASSIGN_ERROR(err, NC_EEDGE)
                goto fn_exit;
            }
        }
        else { /* for vars APIs */
            if (count[i] > 0 &&
                start[i] + (count[i]-1) * stride[i] >= (*varp)->shape[i]) {
                DEBUG_ASSIGN_ERROR(err, NC_EEDGE)
                goto fn_exit;
            }
        }
    }

    if (api <= API_VARA) {
        /* vara APIs have no stride argument */
        err = NC_NOERR;
        goto fn_exit;
    }

    /* Check NC_ESTRIDE */
    for (i=0; i<(*varp)->ndims; i++) {
        if (stride != NULL && stride[i] == 0) {
            DEBUG_ASSIGN_ERROR(err, NC_ESTRIDE)
            goto fn_exit;
        }
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

