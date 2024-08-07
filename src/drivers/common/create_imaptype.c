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

#include <mpi.h>

#include <pnetcdf.h>
#include <pnc_debug.h>
#include <common.h>

/*----< ncmpii_create_imaptype() >-------------------------------------------*/
/* Check if a request is a true varm call. If yes, create an MPI derived
 * data type, imaptype, using imap[]
 */
int
ncmpii_create_imaptype(int               ndims,
                       const MPI_Offset *count,   /* [ndims] */
                       const MPI_Offset *imap,    /* [ndims] */
                       MPI_Datatype      itype,   /* element type */
                       MPI_Datatype     *imaptype)/* out */
{
    int dim, el_size, mpireturn;
    MPI_Offset bnelems, imap_contig_blocklen;

    /* check if this is a vars call or a true varm call */
    *imaptype = MPI_DATATYPE_NULL;

    if (imap == NULL) /* no mapping, same as vars */
        return NC_NOERR;

    if (ndims == 0) /* scalar var, only one value at one fixed place */
        return NC_NOERR;

    for (bnelems=1, dim=0; dim<ndims; dim++) bnelems *= count[dim];
    if (bnelems == 1) return NC_NOERR;

    /* test each dim's contiguity in imap[] until the 1st non-contiguous
     * dim is reached */
    imap_contig_blocklen = 1;
    dim = ndims;
    while (--dim >= 0 && imap_contig_blocklen == imap[dim])
        imap_contig_blocklen *= count[dim];

    if (dim == -1) /* imap is a contiguous layout */
        return NC_NOERR;

    /* itype: element data type (MPI primitive type) */
    mpireturn = MPI_Type_size(itype, &el_size);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size");


    /* We have a true varm call, as imap gives non-contiguous layout.
     * User buffer will be packed (write case) or unpacked (read case)
     * to/from a contiguous buffer based on imap[], before/after MPI-IO.
     * First, we construct a derived data type, imaptype, based on
     * imap[], and use it to pack lbuf to cbuf (for write), or unpack
     * cbuf to lbuf (for read).
     * dim is the first dimension (C order, eg. ZYX) that has
     * non-contiguous imap.
     */
    if (imap_contig_blocklen > NC_MAX_INT || count[dim] > NC_MAX_INT ||
        imap[dim] > NC_MAX_INT) {
#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_vector_c(count[dim], imap_contig_blocklen,
                                      imap[dim], itype, imaptype);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_error_mpi2nc(mpireturn,"MPI_Type_vector_c");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
#else
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
#endif
    }
    else {
        mpireturn = MPI_Type_vector((int)count[dim], (int)imap_contig_blocklen,
                                    (int)imap[dim], itype, imaptype);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_error_mpi2nc(mpireturn,"MPI_Type_vector");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
    }

    mpireturn = MPI_Type_commit(imaptype);
    if (mpireturn != MPI_SUCCESS) {
        ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
        DEBUG_RETURN_ERROR(NC_EMPI)
    }

    for (dim--; dim>=0; dim--) {
        MPI_Datatype tmptype;

        if (count[dim] > NC_MAX_INT) {
#ifdef HAVE_MPI_LARGE_COUNT
            mpireturn = MPI_Type_create_hvector_c(count[dim], 1,
                        imap[dim]*el_size, *imaptype, &tmptype);
            if (mpireturn != MPI_SUCCESS) {
                ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hvector_c");
                DEBUG_RETURN_ERROR(NC_EMPI)
            }
#else
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
#endif
        }
        else {
            mpireturn = MPI_Type_create_hvector((int)count[dim], 1,
                        imap[dim]*el_size, *imaptype, &tmptype);
            if (mpireturn != MPI_SUCCESS) {
                ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hvector");
                DEBUG_RETURN_ERROR(NC_EMPI)
            }
        }

        mpireturn = MPI_Type_free(imaptype);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_error_mpi2nc(mpireturn,"MPI_Type_free");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
        mpireturn = MPI_Type_commit(&tmptype);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
        *imaptype = tmptype;
    }
    return NC_NOERR;
}

