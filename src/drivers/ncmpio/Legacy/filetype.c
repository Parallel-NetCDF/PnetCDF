/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
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
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncx.h"

/*----< ncmpio_calc_datatype_elems() >---------------------------------------*/
/* Obtain the following metadata about buftype:
 * etype:    element data type (MPI primitive type) in buftype
 * bufcount: If it is -1, then this is called from a high-level API and in
 *           this case buftype will be an MPI primitive data type and bufcount
 *           is assigned to a number match with count[]. If bufcount is not -1,
 *           then this is called from a flexible API.
 * nelems:   number of etypes in user buffer
 * xnbytes:  number of bytes (in external data representation) to read/write
 *           from/to the file
 * esize:    byte size of etype
 * isContig: whether buftype is contiguous
 */
int
ncmpio_calc_datatype_elems(const NC_var     *varp,
                           const MPI_Offset *count,
                           MPI_Datatype      buftype,
                           MPI_Datatype     *etype,    /* out */
                           MPI_Offset       *bufcount, /* in/out */
                           MPI_Offset       *nelems,   /* out */
                           MPI_Offset       *xnbytes,  /* out */
                           int              *esize,    /* out */
                           int              *isContig) /* out */
{
    int i;
    MPI_Offset fnelems;

    /* fnelems is the total number of nc_type elements calculated from
     * count[]. count[] is the access count to the variable defined in
     * the netCDF file.
     */
    fnelems = 1;
    for (i=0; i<varp->ndims; i++)
        fnelems *= count[i];

    if (*bufcount == -1) { /* the subroutine is called from a high-level API */
        *bufcount = fnelems;
        *nelems   = fnelems;
        *etype    = buftype; /* buftype is an MPI primitive data type */
        *isContig = 1;
        *xnbytes  = *nelems * varp->xsz;
        MPI_Type_size(buftype, esize);
    }
    else if (buftype == MPI_DATATYPE_NULL) {
        /* This is called from a flexible API and buftype is set by user to
         * MPI_DATATYPE_NULL. In this case, bufcount is ignored and set by
         * this subroutine to a number match count[], and etype to match the
         * variable's external NC data type.
         */
        *bufcount = fnelems;
        *nelems   = fnelems;
        *etype    = ncmpii_nc2mpitype(varp->xtype);
        *esize    = varp->xsz;
        *xnbytes  = *nelems * varp->xsz;
        *isContig = 1;
    }
    else { /* This is called from a flexible API */
        int err, isderived;
        /* check some metadata of the MPI derived datatype */
        err = ncmpii_dtype_decode(buftype, etype, esize, nelems, &isderived,
                                  isContig);
        if (err != NC_NOERR) return err;

        /* make nelems the number of etype in the whole user buf */
        *nelems  *= *bufcount;
        *xnbytes  = *nelems * varp->xsz;

        /* check mismatch between nelems and fnelems */
        if (fnelems != *nelems) DEBUG_RETURN_ERROR(NC_EIOMISMATCH)
    }
    return NC_NOERR;
}

