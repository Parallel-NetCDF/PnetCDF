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
#include <string.h>
#include <strings.h>  /* strcasecmp() */
#include <assert.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncx.h"

/*----< ncmpio_sanity_check() >----------------------------------------------*/
/* check the following errors and in that precedence.
 * NC_EBADID, NC_EPERM, NC_EINDEFINE, NC_EINDEP/NC_ENOTINDEP, NC_ENOTVAR,
 * NC_ECHAR, NC_EINVALCOORDS, NC_EEDGE, NC_ESTRIDE, NC_EINVAL.
 */
int ncmpio_sanity_check(NC                *ncp,
                        int               varid,
                        const MPI_Offset  bufcount,
                        MPI_Datatype      buftype,  /* internal datatype */
                        int               reqMode,
                        NC_var          **varp)  /* OUT */
{
    /* all errors detected here are fatal, must return immediately */
    int err=NC_NOERR;

    /* check file write permission if this is write request */
    if (fIsSet(reqMode, NC_REQ_WR) && NC_readonly(ncp))
        DEBUG_RETURN_ERROR(NC_EPERM)

    if (fIsSet(reqMode, NC_REQ_BLK)) {
        /* blocking APIs must be called in data mode */
        if (NC_indef(ncp))
            DEBUG_RETURN_ERROR(NC_EINDEFINE)

        /* for blocking APIs, check if in the right collective or independent
         * mode, nonblocking APIs can be called in either mode */
        if (fIsSet(reqMode, NC_REQ_INDEP) && !NC_indep(ncp))
            DEBUG_RETURN_ERROR(NC_ENOTINDEP)
        else if (fIsSet(reqMode, NC_REQ_COLL) && NC_indep(ncp))
            DEBUG_RETURN_ERROR(NC_EINDEP)
    }

    if (fIsSet(reqMode, NC_REQ_ZERO)) return NC_NOERR;

    /* check if varid is valid (check NC_ENOTVAR) */
    err = ncmpio_NC_lookupvar(ncp, varid, varp);
    if (err != NC_NOERR) return err;

    /* check NC_ECHAR */
    if (fIsSet(reqMode, NC_REQ_FLEX)) {
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
            if (err != NC_NOERR) return err;

            err = NCMPII_ECHAR((*varp)->type, ptype);
            if (err != NC_NOERR) return err;
        }
        /* else case: itype matches xtype */

        /* for flexible APIs, bufcount cannot be negative */
        if (bufcount < 0) DEBUG_RETURN_ERROR(NC_EINVAL)
    }
    else { /* called from a high-level API */
        err = NCMPII_ECHAR((*varp)->type, buftype);
        if (err != NC_NOERR) return err;
    }
    return NC_NOERR;
}
