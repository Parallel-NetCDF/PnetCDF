/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include <pnetcdf.h>
#include <pnc_debug.h>
#include <common.h>

/*----< ncmpii_pack() >------------------------------------------------------*/
/* This subroutine packs buf, if it is noncontiguous, into a contiguous
 * buffer, cbuf. Whether buf is contiguous or not depends on buftype and imap.
 * If both buftype and imap indicate a contiguous buffer, then cbuf will be
 * set to buf. Otherwise, cbuf will be malloc-ed and it needs to be freed
 * later.
 * buftype is decoded to find the element type, ptype, used to create the type.
 * bnelems is calculated base don bufcount and buftype.
 */
int
ncmpii_pack(int                ndims,
            const MPI_Offset  *count,
            const MPI_Offset  *imap,     /* can be NULL */
            void              *buf,      /* user buffer */
            MPI_Offset         bufcount, /* -1: from high-level API */
            MPI_Datatype       buftype,  /* MPI derived datatype */
            MPI_Offset        *bnelems,  /* OUT: no. of ptypes in buf */
            MPI_Datatype      *ptype,    /* OUT: MPI primitive datatype */
            void             **cbuf)     /* OUT: a contiguous buffer */
{
    void *lbuf=NULL;
    int i, err=NC_NOERR, mpireturn;
    MPI_Offset buf_size, nelems;
    MPI_Datatype etype, imaptype=MPI_DATATYPE_NULL;

#ifdef HAVE_MPI_TYPE_SIZE_C
    MPI_Count type_size;
    mpireturn = MPI_Type_size_c(buftype, &type_size);
#elif defined(HAVE_MPI_TYPE_SIZE_X)
    MPI_Count type_size;
    mpireturn = MPI_Type_size_x(buftype, &type_size);
#else
    int type_size;
    mpireturn = MPI_Type_size(buftype, &type_size);
#endif
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size");
        DEBUG_RETURN_ERROR(err)
    }
    else if (type_size == MPI_UNDEFINED)
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    buf_size = type_size;

    *cbuf = buf;

    for (nelems=1, i=0; i<ndims; i++) nelems *= count[i];

    lbuf  = buf;
    if (bufcount == -1) { /* called from a high-level API */
        buf_size *= nelems;
        etype     = buftype;

        if (bnelems != NULL) *bnelems = nelems;
        if (ptype   != NULL) *ptype   = buftype;

        if (buf_size == 0) /* zero-length request */
            return NC_NOERR;

        if (imap == NULL) /* not called from a true varm API */
            return NC_NOERR;
    }
    else { /* called from a flexible API */
        int el_size, isDerived, isContig;
        MPI_Offset num_ptypes=0;

        buf_size *= bufcount;

        /* check if buftype is an MPI predefined primitive datatype */
        err = ncmpii_dtype_decode(buftype, &etype, &el_size, &num_ptypes,
                                  &isDerived, &isContig);
        if (err != NC_NOERR) return err;

        if (buf_size == 0) { /* zero-length request */
            if (bnelems != NULL) *bnelems = 0;
            if (ptype   != NULL) *ptype   = etype;
            return NC_NOERR;
        }

        num_ptypes *= bufcount;

        if (bnelems != NULL) *bnelems = num_ptypes;
        if (ptype   != NULL) *ptype   = etype;

        /* check if number of elements of bufcount/buftype matches count[] */
        if (num_ptypes != nelems)
            DEBUG_RETURN_ERROR(NC_EIOMISMATCH)

        if (isDerived) { /* Not a predefined datatype */
            /* allocate lbuf and pack buf into lbuf */
            lbuf = NCI_Malloc((size_t)buf_size);
            if (lbuf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count position=0;
            MPI_Pack_c(buf, (MPI_Count)bufcount, buftype, lbuf,
                       (MPI_Count)buf_size, &position, MPI_COMM_SELF);
#else
            int position=0;
            if (buf_size > NC_MAX_INT)
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

            MPI_Pack(buf, (int)bufcount, buftype, lbuf, (int)buf_size,
                     &position, MPI_COMM_SELF);
#endif
        }
    }

    /* Construct a derived datatype, imaptype, if this is a true varm call */
    err = ncmpii_create_imaptype(ndims, count, imap, etype, &imaptype);
    if (err != NC_NOERR) return err;

    /* Step 2: pack lbuf to cbuf if imap is non-contiguous */
    if (imaptype != MPI_DATATYPE_NULL) { /* true varm */
        /* pack lbuf to cbuf, a contiguous buffer, using imaptype */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count position=0;
        *cbuf = NCI_Malloc((size_t)buf_size);
        MPI_Pack_c(lbuf, 1, imaptype, *cbuf, (MPI_Count)buf_size, &position,
                   MPI_COMM_SELF);
#else
        int position=0;
        if (buf_size > NC_MAX_INT)
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

        *cbuf = NCI_Malloc((size_t)buf_size);
        MPI_Pack(lbuf, 1, imaptype, *cbuf, (int)buf_size, &position,
                 MPI_COMM_SELF);
#endif
        MPI_Type_free(&imaptype);
    }
    else /* reuse lbuf */
        *cbuf = lbuf;

    /* lbuf is no longer needed */
    if (lbuf != buf && lbuf != *cbuf) NCI_Free(lbuf);

    return NC_NOERR;
}

