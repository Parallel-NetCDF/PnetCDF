/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <mpi.h>

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

/*----< ncmpii_start_count_stride_check() >----------------------------------*/
/*
 * Check whether start, count, stride values are valid for the variable.
 * Report error codes: NC_EINVALCOORDS, NC_EEDGE, or NC_ESTRIDE.
 * Note that even if the request size is zero, this check is enforced in both
 * netCDF and PnetCDF. Otherwise, many test cases under test directory can
 * fail. Arguments count and stride can be NULL.
 */
int
ncmpii_start_count_stride_check(int               format,
                                int               api,
                                int               ndims,
                                int               numrecs,
                                const MPI_Offset *shape,
                                const MPI_Offset *start,
                                const MPI_Offset *count,
                                const MPI_Offset *stride,
                                const int         reqMode) /* read or write */
{
    int i=0, firstDim;

    if (ndims == 0) return NC_NOERR; /* 'scalar' variable */

    if (api <= API_VAR) /* var/varn/vard APIs, start/count/stride are NULL */
        return NC_NOERR;

    /* Now only need to check var1, vara, vars, and varm APIs */

    /* Check NC_EINVALCOORDS error for argument start[]
     * for API var1/vara/vars/varm, start cannot be NULL, except for scalars
     * and negative start[] is illegal */
    if (start == NULL || start[0] < 0) DEBUG_RETURN_ERROR(NC_EINVALCOORDS)

    firstDim = 0;
    /* check NC_EINVALCOORDS for record dimension */
    if (shape[0] == NC_UNLIMITED) {
        if (format < 5 && start[0] > NC_MAX_UINT) /* CDF-1 and 2 */
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)

        /* for record variable, [0] is the NC_UNLIMITED dimension */
        if (fIsSet(reqMode, NC_REQ_RD)) {
            /* read cannot go beyond current numrecs */
#ifdef RELAX_COORD_BOUND
            if (start[0] > numrecs) DEBUG_RETURN_ERROR(NC_EINVALCOORDS)

            if (start[0] == numrecs) {
                if (api == API_VAR1) {
                    /* for var1 APIs, count[0] is considered of 1 */
                    DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
                }
                else if (count != NULL && count[0] > 0) {
                    DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
                }
            }
#else
            if (start[0] >= numrecs) DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
#endif
        }
        firstDim = 1; /* done for checking the record dimension */
    }

    /* continue to check NC_EINVALCOORDS for the rest dimensions */
    for (i=firstDim; i<ndims; i++) {
#ifdef RELAX_COORD_BOUND
        if (start[i] < 0 || start[i] > shape[i])
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)

        if (start[i] == shape[i]) {
            if (api == API_VAR1) {
                /* for var1 APIs, count[0] is considered of 1 */
                DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
            }
            else if (count != NULL && count[i] > 0) {
                DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
            }
        }
#else
        if (start[i] < 0 || start[i] >= shape[i])
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
#endif
    }

    /* Now check NC_EEDGE error for argument count[] */
    if (api <= API_VAR1) /* var1/var APIs have no count argument */
        return NC_NOERR;

    /* for API vara/vars/varm, count cannot be NULL, except for scalars */
    if (count == NULL) DEBUG_RETURN_ERROR(NC_EEDGE)

    firstDim = 0;
    /* check NC_EINVALCOORDS for record dimension */
    if (shape[0] == NC_UNLIMITED) {
        if (count[0] < 0)  /* no negative count[] */
            DEBUG_RETURN_ERROR(NC_ENEGATIVECNT)

        /* for record variable, [0] is the NC_UNLIMITED dimension */
        if (fIsSet(reqMode, NC_REQ_RD)) {
            /* read cannot go beyond current numrecs */
            if (stride == NULL) {  /* for vara APIs */
                if (start[0] + count[0] > numrecs)
                    DEBUG_RETURN_ERROR(NC_EEDGE)
            }
            else { /* for vars/varm APIs */
                if (count[0] > 0 &&
                    start[0] + (count[0]-1) * stride[0] >= numrecs)
                    DEBUG_RETURN_ERROR(NC_EEDGE)
            }
        }
        firstDim = 1; /* skip checking the record dimension */
    }

    /* continue to check NC_EEDGE for the rest dimensions */
    for (i=firstDim; i<ndims; i++) {
        if (shape[i] < 0)
            DEBUG_RETURN_ERROR(NC_EEDGE)

        if (count[i] < 0) /* no negative count[] */
            DEBUG_RETURN_ERROR(NC_ENEGATIVECNT)

        if (stride == NULL) { /* for vara APIs */
            if (count[i] > shape[i] || start[i] + count[i] > shape[i])
                DEBUG_RETURN_ERROR(NC_EEDGE)
        }
        else { /* for vars APIs */
            if (count[i] > 0 && start[i] + (count[i]-1) * stride[i] >= shape[i])
                DEBUG_RETURN_ERROR(NC_EEDGE)
        }
    }

    /* Now check NC_ESTRIDE error for argument stride[] */
    if (api <= API_VARA || stride == NULL)
        /* vara APIs have no stride argument */
        return NC_NOERR;

    /* Check NC_ESTRIDE for non-positive values. We did not check
     * stride[i] >= shape[i], as it is caught as NC_EEDGE error above */
    for (i=0; i<ndims; i++) {
        if (stride[i] <= 0)
            DEBUG_RETURN_ERROR(NC_ESTRIDE)
    }
    return NC_NOERR;
}
