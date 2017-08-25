/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <string.h>  /* memcpy(), memcmp() */
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "nc.h"
#include "ncx.h"

typedef enum {
    NC_UNSPECIFIED =  0,
/*  NC_BITFIELD    =  7, */
/*  NC_STRING      =  8, */
    NC_DIMENSION   = 10,  /* \x00 \x00 \x00 \x0A */
    NC_VARIABLE    = 11,  /* \x00 \x00 \x00 \x0B */
    NC_ATTRIBUTE   = 12   /* \x00 \x00 \x00 \x0C */
} NC_tag;

/*
 * "magic number" at beginning of file: 0x43444601 (big endian)
 */
static const char ncmagic1[] = {'C', 'D', 'F', 0x01};
static const char ncmagic2[] = {'C', 'D', 'F', 0x02};
static const char ncmagic5[] = {'C', 'D', 'F', 0x05};


/*----< compute_var_shape() >------------------------------------------------*/
/* Recompute the shapes of all variables
 * Sets ncp->begin_var to start of first variable.
 * Sets ncp->begin_rec to start of first record variable.
 * Returns -1 on error. The only possible error is an reference to a non
 * existent dimension, which would occur for a corrupt netcdf file.
 */
static int
compute_var_shape(NC *ncp)
{
    int i, err;
    NC_var *first_var = NULL;       /* first "non-record" var */
    NC_var *first_rec = NULL;       /* first "record" var */

    if (ncp->vars.ndefined == 0) return NC_NOERR;

    ncp->begin_var = ncp->xsz;
    ncp->begin_rec = ncp->xsz;
    ncp->recsize   = 0;

    for (i=0; i<ncp->vars.ndefined; i++) {
        /* ncp->vars.value[i]->len will be recomputed from dimensions in
         * ncmpio_NC_var_shape64() */
        err = ncmpio_NC_var_shape64(ncp->vars.value[i], &ncp->dims);
        if (err != NC_NOERR) return err;

        if (IS_RECVAR(ncp->vars.value[i])) {
            if (first_rec == NULL) first_rec = ncp->vars.value[i];
            ncp->recsize += ncp->vars.value[i]->len;
        }
        else { /* fixed-size variable */
            if (first_var == NULL) first_var = ncp->vars.value[i];
            /*
             * Overwritten each time thru.
             * Usually overwritten in first_rec != NULL clause.
             */
            ncp->begin_rec = ncp->vars.value[i]->begin
                           + ncp->vars.value[i]->len;
        }
    }

    if (first_rec != NULL) {
        if (ncp->begin_rec > first_rec->begin)
            DEBUG_RETURN_ERROR(NC_ENOTNC) /* not a netCDF file or corrupted */

        ncp->begin_rec = first_rec->begin;
        /*
         * for special case of exactly one record variable, pack value
         */
        if (ncp->recsize == first_rec->len)
            ncp->recsize = *first_rec->dsizes * first_rec->xsz;
    }

    if (first_var != NULL)
        ncp->begin_var = first_var->begin;
    else
        ncp->begin_var = ncp->begin_rec;

    if (ncp->begin_var <= 0 || ncp->xsz > ncp->begin_var ||
        ncp->begin_rec <= 0 || ncp->begin_var > ncp->begin_rec)
        DEBUG_RETURN_ERROR(NC_ENOTNC) /* not a netCDF file or corrupted */

    return NC_NOERR;
}

#define X_SIZEOF_NC_TYPE X_SIZEOF_INT
#define X_SIZEOF_NC_TAG  X_SIZEOF_INT

/*----< hdr_len_NC_name() >--------------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_name(const NC_string *ncstrp,
                int              sizeof_NON_NEG)     /* NON_NEG */
{
    /* netCDF file format:
     * name       = nelems  namestring
     * nelems     = NON_NEG
     * namestring = ID1 [IDN ...] padding
     * ID1        = alphanumeric | '_'
     * IDN        = alphanumeric | special1 | special2
     * padding    = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz = sizeof_NON_NEG; /* nelems */

    assert(ncstrp != NULL);

    if (ncstrp->nchars != 0)  /* namestring */
        sz += _RNDUP(ncstrp->nchars, X_ALIGN);

    return sz;
}

/*----< hdr_len_NC_dim() >---------------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_dim(const NC_dim *dimp,
               int           sizeof_NON_NEG)     /* NON_NEG */
{
    /* netCDF file format:
     *  ...
     * dim        = name  dim_length
     * dim_length = NON_NEG
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;

    assert(dimp != NULL);

    sz = hdr_len_NC_name(dimp->name, sizeof_NON_NEG); /* name */
    sz += sizeof_NON_NEG;                             /* dim_length */

    return sz;
}

/*----< hdr_len_NC_dimarray() >----------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_dimarray(const NC_dimarray *ncap,
                    int                sizeof_NON_NEG)     /* NON_NEG */
{
    /* netCDF file format:
     *  ...
     * dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_DIMENSION = \x00 \x00 \x00 \x0A         // tag for list of dimensions
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i;
    MPI_Offset xlen;

    xlen = X_SIZEOF_NC_TAG;           /* NC_DIMENSION */
    xlen += sizeof_NON_NEG;           /* nelems */

    if (ncap == NULL) /* ABSENT: no dimension is defined */
        return xlen;

    /* [dim ...] */
    for (i=0; i<ncap->ndefined; i++)
        xlen += hdr_len_NC_dim(ncap->value[i], sizeof_NON_NEG);

    return xlen;
}

/*----< hdr_len_NC_attr() >--------------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_attr(const NC_attr *attrp,
                int            sizeof_NON_NEG)     /* NON_NEG */
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     * nc_type = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * nelems  = NON_NEG       // number of elements in following sequence
     * values  = bytes | chars | shorts | ints | floats | doubles
     * bytes   = [BYTE ...]  padding
     * chars   = [CHAR ...]  padding
     * shorts  = [SHORT ...]  padding
     * ints    = [INT ...]
     * floats  = [FLOAT ...]
     * doubles = [DOUBLE ...]
     * padding = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG = <non-negative INT> |  // CDF-1 and CDF-2
     *           <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;

    assert(attrp != NULL);

    sz  = hdr_len_NC_name(attrp->name, sizeof_NON_NEG); /* name */
    sz += X_SIZEOF_NC_TYPE;                             /* nc_type */
    sz += sizeof_NON_NEG;                               /* nelems */
    sz += attrp->xsz;                                   /* [values ...] */

    return sz;
}

/*----< hdr_len_NC_attrarray() >---------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_attrarray(const NC_attrarray *ncap,
                     int                 sizeof_NON_NEG)     /* NON_NEG */
{
    /* netCDF file format:
     *  ...
     * att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_ATTRIBUTE = \x00 \x00 \x00 \x0C         // tag for list of attributes
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i;
    MPI_Offset xlen;

    xlen = X_SIZEOF_NC_TAG;        /* NC_ATTRIBUTE */
    xlen += sizeof_NON_NEG;        /* nelems */

    if (ncap == NULL) /* ABSENT: no attribute is defined */
        return xlen;

    for (i=0; i<ncap->ndefined; i++) /* [attr ...] */
        xlen += hdr_len_NC_attr(ncap->value[i], sizeof_NON_NEG);

    return xlen;
}

/*----< hdr_len_NC_var() >---------------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_var(const NC_var *varp,
               int           sizeof_off_t, /* OFFSET */
               int           sizeof_NON_NEG)     /* NON_NEG */
{
    /* netCDF file format:
     * netcdf_file = header data
     * header      = magic numrecs dim_list gatt_list var_list
     *  ...
     * var         = name nelems [dimid ...] vatt_list nc_type vsize begin
     * nelems      = NON_NEG
     * dimid       = NON_NEG
     * vatt_list   = att_list
     * nc_type     = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * vsize       = NON_NEG
     * begin       = OFFSET        // Variable start location.
     * OFFSET      = <non-negative INT> |  // CDF-1
     *               <non-negative INT64>  // CDF-2 and CDF-5
     * NON_NEG     = <non-negative INT> |  // CDF-1 and CDF-2
     *               <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;

    assert(varp != NULL);

    /* for CDF-1, sizeof_off_t == 4 && sizeof_NON_NEG == 4
     * for CDF-2, sizeof_off_t == 8 && sizeof_NON_NEG == 4
     * for CDF-5, sizeof_off_t == 8 && sizeof_NON_NEG == 8
     */
    sz = hdr_len_NC_name(varp->name, sizeof_NON_NEG);         /* name */
    sz += sizeof_NON_NEG;                                     /* nelems */
    sz += sizeof_NON_NEG * varp->ndims;                       /* [dimid ...] */
    sz += hdr_len_NC_attrarray(&varp->attrs, sizeof_NON_NEG); /* vatt_list */
    sz += X_SIZEOF_NC_TYPE;                                   /* nc_type */
    sz += sizeof_NON_NEG;                                     /* vsize */
    sz += sizeof_off_t;                                       /* begin */

    return sz;
}

/*----< hdr_len_NC_vararray() >----------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_vararray(const NC_vararray *ncap,
                    int                sizeof_NON_NEG,     /* NON_NEG */
                    int                sizeof_off_t) /* OFFSET */
{
    /* netCDF file format:
     * netcdf_file = header  data
     * header      = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var_list    = ABSENT | NC_VARIABLE   nelems  [var ...]
     * ABSENT      = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *               ZERO  ZERO64  // for CDF-5
     * ZERO        = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64      = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_VARIABLE = \x00 \x00 \x00 \x0B         // tag for list of variables
     * nelems      = NON_NEG       // number of elements in following sequence
     * NON_NEG     = <non-negative INT> |        // CDF-1 and CDF-2
     *               <non-negative INT64>        // CDF-5
     */
    int i;
    MPI_Offset xlen;

    xlen = X_SIZEOF_NC_TAG;           /* NC_VARIABLE */
    xlen += sizeof_NON_NEG;           /* nelems */

    if (ncap == NULL) /* ABSENT: no variable is defined */
        return xlen;

    /* for CDF-1, sizeof_off_t == 4 && sizeof_NON_NEG == 4
     * for CDF-2, sizeof_off_t == 8 && sizeof_NON_NEG == 4
     * for CDF-5, sizeof_off_t == 8 && sizeof_NON_NEG == 8
     */
    for (i=0; i<ncap->ndefined; i++)  /* [var ...] */
        xlen += hdr_len_NC_var(ncap->value[i], sizeof_off_t, sizeof_NON_NEG);

    return xlen;
}

/*----< hdr_put_NC_name() >--------------------------------------------------*/
inline static int
hdr_put_NC_name(bufferinfo      *pbp,
                const NC_string *ncstrp)
{
    /* netCDF file format:
     *  ...
     * name       = nelems  namestring
     * nelems     = NON_NEG
     * namestring = ID1 [IDN ...] padding
     * ID1        = alphanumeric | '_'
     * IDN        = alphanumeric | special1 | special2
     * padding    = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    int err;

    /* copy nelems */
    if (pbp->version < 5) {
        if (ncstrp->nchars != (uint)ncstrp->nchars)
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        err = ncmpix_put_uint32((void**)(&pbp->pos), (uint)ncstrp->nchars);
    }
    else
        err = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)ncstrp->nchars);
    if (err != NC_NOERR) return err;

    /* copy namestring */
    err = ncmpix_pad_putn_text(&pbp->pos, ncstrp->nchars, ncstrp->cp);

    return err;
}

/*----< hdr_put_NC_dim() >---------------------------------------------------*/
inline static int
hdr_put_NC_dim(bufferinfo   *pbp,
               const NC_dim *dimp)
{
    /* netCDF file format:
     *  ...
     * dim        = name  dim_length
     * dim_length = NON_NEG
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    int err;

    /* copy name */
    err = hdr_put_NC_name(pbp, dimp->name);
    if (err != NC_NOERR) return err;

    /* copy dim_length */
    if (pbp->version < 5) {
        /* TODO: Isn't checking dimension size already done in def_dim()? */
        if (dimp->size != (uint)dimp->size) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        err = ncmpix_put_uint32((void**)(&pbp->pos), (uint)dimp->size);
    }
    else
        err = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)dimp->size);

    return err;
}

/*----< hdr_put_NC_dimarray() >----------------------------------------------*/
inline static int
hdr_put_NC_dimarray(bufferinfo        *pbp,
                    const NC_dimarray *ncap)
{
    /* netCDF file format:
     *  ...
     * dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_DIMENSION = \x00 \x00 \x00 \x0A         // tag for list of dimensions
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i, status;

    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        if (pbp->version < 5)
            status = ncmpix_put_uint32((void**)(&pbp->pos), 0);
        else
            status = ncmpix_put_uint64((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_DIMENSION */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_DIMENSION);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        if (pbp->version < 5)
            status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)ncap->ndefined);
        else
            status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)ncap->ndefined);
        if (status != NC_NOERR) return status;

        /* copy [dim ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = hdr_put_NC_dim(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< hdr_put_NC_attrV() >-------------------------------------------------*/
/*
 * Put the values of an attribute
 */
inline static int
hdr_put_NC_attrV(bufferinfo    *pbp,
                 const NC_attr *attrp)
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     *  ...
     * values  = bytes | chars | shorts | ints | floats | doubles
     * bytes   = [BYTE ...]  padding
     * chars   = [CHAR ...]  padding
     * shorts  = [SHORT ...]  padding
     * ints    = [INT ...]
     * floats  = [FLOAT ...]
     * doubles = [DOUBLE ...]
     * padding = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     */
    MPI_Offset padding, sz;

    /* ncmpio_xlen_nc_type() returns the element size (unaligned) of attrp->type
       attrp->xsz is the aligned total size of attribute values
     */
    sz = ncmpio_xlen_nc_type(attrp->type);
    sz *= attrp->nelems;
    padding = attrp->xsz - sz;

    if (sz != (size_t) sz) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    memcpy(pbp->pos, attrp->xvalue, (size_t)sz);
    pbp->pos = (void *)((char *)pbp->pos + sz);

    if (padding > 0) {
        /* zero-padding is per buffer, not per element */
        memset(pbp->pos, 0, (size_t)padding);
        pbp->pos = (void *)((char *)pbp->pos + padding);
    }

    return NC_NOERR;
}

/*----< hdr_put_NC_attr() >--------------------------------------------------*/
inline static int
hdr_put_NC_attr(bufferinfo    *pbp,
                const NC_attr *attrp)
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     * nc_type = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * nelems  = NON_NEG       // number of elements in following sequence
     * NON_NEG = <non-negative INT> |  // CDF-1 and CDF-2
     *           <non-negative INT64>  // CDF-5
     */
    int status;

    /* copy name */
    status = hdr_put_NC_name(pbp, attrp->name);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)attrp->type);
    if (status != NC_NOERR) return status;

    /* copy nelems */
    if (pbp->version < 5) {
        if (attrp->nelems != (uint)attrp->nelems)
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)attrp->nelems);
    }
    else
        status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)attrp->nelems);
    if (status != NC_NOERR) return status;

    /* copy [values ...] */
    status = hdr_put_NC_attrV(pbp, attrp);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

/*----< hdr_put_NC_attrarray() >---------------------------------------------*/
inline static int
hdr_put_NC_attrarray(bufferinfo         *pbp,
                     const NC_attrarray *ncap)
{
    /* netCDF file format:
     *  ...
     * att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_ATTRIBUTE = \x00 \x00 \x00 \x0C         // tag for list of attributes
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i, status;

    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        if (pbp->version < 5)
            status = ncmpix_put_uint32((void**)(&pbp->pos), 0);
        else
            status = ncmpix_put_uint64((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_ATTRIBUTE */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_ATTRIBUTE);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        if (pbp->version < 5)
            status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)ncap->ndefined);
        else
            status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)ncap->ndefined);
        if (status != NC_NOERR) return status;

        /* copy [attr ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = hdr_put_NC_attr(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< hdr_put_NC_var() >---------------------------------------------------*/
static int
hdr_put_NC_var(bufferinfo   *pbp,
               const NC_var *varp)
{
    /* netCDF file format:
     * netcdf_file = header data
     * header      = magic numrecs dim_list gatt_list var_list
     *  ...
     * var         = name nelems [dimid ...] vatt_list nc_type vsize begin
     * nelems      = NON_NEG
     * dimid       = NON_NEG
     * vatt_list   = att_list
     * nc_type     = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * vsize       = NON_NEG
     * begin       = OFFSET        // Variable start location.
     * OFFSET      = <non-negative INT> |  // CDF-1
     *               <non-negative INT64>  // CDF-2 and CDF-5
     * NON_NEG     = <non-negative INT> |  // CDF-1 and CDF-2
     *               <non-negative INT64>  // CDF-5
     */
    int i, status;

    /* copy name */
    status = hdr_put_NC_name(pbp, varp->name);
    if (status != NC_NOERR) return status;

    /* copy nelems */
    if (pbp->version < 5)
        status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)varp->ndims);
    else
        status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)varp->ndims);
    if (status != NC_NOERR) return status;

    /* copy [dimid ...] */
    for (i=0; i<varp->ndims; i++) {
        if (pbp->version < 5)
            status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)varp->dimids[i]);
        else
            status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)varp->dimids[i]);
        if (status != NC_NOERR) return status;
    }

    /* copy vatt_list */
    status = hdr_put_NC_attrarray(pbp, &varp->attrs);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)varp->type);
    if (status != NC_NOERR) return status;

    /* copy vsize */
    /* in CDF-1 and CDF-2, a variable's size in the header is a 32-bit integer
     * in CDF-5, it is a 64-bit integer
     */
    if (pbp->version < 5) {
        /* Special case, when there is no record variable, the last fixed-size
         * variable can be larger than 2 GiB if its file starting offset is
         * less than 2 GiB. This checking has already been done in the call
         * to ncmpio_NC_check_vlens() in ncmpio_NC_enddef().
         *
         * if (varp->len != (int)varp->len) DEBUG_RETURN_ERROR(NC_EVARSIZE)
         */
        uint vsize = (uint)varp->len;
        if (varp->len > 4294967292LL) { /* 2^32 - 4 bytes */
            /* CDF-2 specification: use 2^32-1 for vsize when the variable
             * size is larger than 2^32-4 bytes
             */
            vsize = 4294967295U;
        }
        status = ncmpix_put_uint32((void**)(&pbp->pos), vsize);
    }
    else {
        status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)varp->len);
    }
    if (status != NC_NOERR) return status;

    /* copy begin */
    /* in CDF-1 header, a variable's starting file offset is a 32-bit integer
     * in CDF-2 and CDF-5, it is a 64-bit integer
     */
    if (pbp->version == 1) {
        if (varp->begin != (uint)varp->begin) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)varp->begin);
    }
    else
        status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)varp->begin);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

/*----< hdr_put_NC_vararray() >----------------------------------------------*/
static int
hdr_put_NC_vararray(bufferinfo        *pbp,
                    const NC_vararray *ncap)
{
    /* netCDF file format:
     * netcdf_file = header  data
     * header      = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var_list    = ABSENT | NC_VARIABLE   nelems  [var ...]
     * ABSENT      = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *               ZERO  ZERO64  // for CDF-5
     * ZERO        = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64      = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_VARIABLE = \x00 \x00 \x00 \x0B         // tag for list of variables
     * nelems      = NON_NEG       // number of elements in following sequence
     * NON_NEG     = <non-negative INT> |        // CDF-1 and CDF-2
     *               <non-negative INT64>        // CDF-5
     */
    int i, status;

    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        if (pbp->version < 5)
            status = ncmpix_put_uint32((void**)(&pbp->pos), 0);
        else
            status = ncmpix_put_uint64((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_VARIABLE */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_VARIABLE);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        if (pbp->version < 5)
            status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)ncap->ndefined);
        else
            status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)ncap->ndefined);
        if (status != NC_NOERR) return status;

        /* copy [var ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = hdr_put_NC_var(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< ncmpio_hdr_put_NC() >------------------------------------------------*/
/* fill the file header into the I/O buffer, buf
 * this function is collective */
int
ncmpio_hdr_put_NC(NC *ncp, void *buf)
{
    int status;
    bufferinfo putbuf;
    MPI_Offset nrecs=0;

    putbuf.comm          = ncp->comm;
    putbuf.collective_fh = ncp->collective_fh;
    putbuf.offset        = 0;
    putbuf.pos           = buf;
    putbuf.base          = buf;
    putbuf.size          = ncp->xsz;
    putbuf.safe_mode     = ncp->safe_mode;

    /* netCDF file format:
     * netcdf_file  = header  data
     * header       = magic  numrecs  dim_list  gatt_list  var_list
     */

    /* copy "magic", 4 characters */
    if (ncp->format == 5) {
        putbuf.version = 5;
        status = ncmpix_putn_text(&putbuf.pos, sizeof(ncmagic5), ncmagic5);
    }
    else if (ncp->format == 2) {
        putbuf.version = 2;
        status = ncmpix_putn_text(&putbuf.pos, sizeof(ncmagic2), ncmagic2);
    }
    else {
        putbuf.version = 1;
        status = ncmpix_putn_text(&putbuf.pos, sizeof(ncmagic1), ncmagic1);
    }
    if (status != NC_NOERR) return status;

    /* copy numrecs, number of records */
    nrecs = ncp->numrecs;
    if (ncp->format < 5) {
        if (nrecs != (uint)nrecs) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        status = ncmpix_put_uint32((void**)(&putbuf.pos), (uint)nrecs);
    }
    else {
        status = ncmpix_put_uint64((void**)(&putbuf.pos), (uint64)nrecs);
    }
    if (status != NC_NOERR) return status;

    /* copy dim_list */
    status = hdr_put_NC_dimarray(&putbuf, &ncp->dims);
    if (status != NC_NOERR) return status;

    /* copy gatt_list */
    status = hdr_put_NC_attrarray(&putbuf, &ncp->attrs);
    if (status != NC_NOERR) return status;

    /* copy var_list */
    status = hdr_put_NC_vararray(&putbuf, &ncp->vars);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

/*----< hdr_fetch() >--------------------------------------------------------*/
/* Fetch the next header chunk.  the chunk is 'gbp->size' bytes big
 * Takes care to not overwrite leftover (unused) data in the buffer before
 * fetching a new chunk: the current approach is to re-read the extra data.
 *
 * NOTE: An alternate approach (which we do not do) would be to save the old
 *       data, read the next chunk and then copy the old data into the new
 *       chunk.  This alternate approach might help if it is important for
 *       reads to be aligned.
 */
static int
hdr_fetch(bufferinfo *gbp) {
    int rank, err=NC_NOERR, mpireturn;
    MPI_Offset slack;        /* any leftover data in the buffer */
    MPI_Aint pos_addr, base_addr;

    assert(gbp->base != NULL);

#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    slack = gbp->size - (pos_addr - base_addr);
    /* . if gbp->pos and gbp->base are the same, there is no leftover buffer
     *   data to worry about.
     * In the other extreme, where gbp->size == (gbp->pos - gbp->base), then
     * all data in the buffer has been consumed */
    if (slack == gbp->size) slack = 0;

    if (gbp->size != (int)gbp->size) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    memset(gbp->base, 0, (size_t)gbp->size);
    gbp->pos = gbp->base;

    MPI_Comm_rank(gbp->comm, &rank);
    if (rank == 0) {
        MPI_Status mpistatus;
        /* fileview is already entire file visible and MPI_File_read_at does
           not change the file pointer */
        TRACE_IO(MPI_File_read_at)(gbp->collective_fh,
                                   (gbp->offset)-slack, gbp->base,
                                   (int)gbp->size, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpio_handle_error(mpireturn, "MPI_File_read_at");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD)
        }
        else {
            int get_size; /* actual read amount can be smaller */
            MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
            gbp->get_size += get_size;
        }
    }
    /* we might have to backtrack */
    gbp->offset += (gbp->size - slack);

    if (gbp->safe_mode == 1) {
        TRACE_COMM(MPI_Bcast)(&err, 1, MPI_INT, 0, gbp->comm);
        if (err != NC_NOERR) return err;
    }

    /* broadcast root's read (full or partial header) to other processes */
    TRACE_COMM(MPI_Bcast)(gbp->base, (int)gbp->size, MPI_BYTE, 0, gbp->comm);

    return err;
}

/*----< hdr_check_buffer() >--------------------------------------------------*/
/* Ensure that 'nextread' bytes are available.  */
inline static int
hdr_check_buffer(bufferinfo *gbp,
                 MPI_Offset  nextread)
{
    MPI_Aint pos_addr, base_addr;

#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    if (pos_addr + nextread <= base_addr + gbp->size)
        return NC_NOERR;

    return hdr_fetch(gbp);
}

/*----< hdr_get_uint32() >---------------------------------------------------*/
inline static int
hdr_get_uint32(bufferinfo *gbp,
               uint        *xp)
{
    /* in CDF-1 format, all integers are 32-bit
     * in CDF-2 format, only variable begin (starting file offset) is 64-bit
     * in CDF-5 format, both variable's begin and size are 64-bit
     */
    int status = hdr_check_buffer(gbp, 4); /* size of int32 == 4 */
    if (status != NC_NOERR) return status;

    status = ncmpix_get_uint32((const void **)(&gbp->pos), xp);
    return status;
}

/*----< hdr_get_uint64() >---------------------------------------------------*/
inline static int
hdr_get_uint64(bufferinfo *gbp,
               uint64     *xp)
{
    /* in CDF-1 format, all integers are 32-bit
     * in CDF-2 format, only variable begin (starting file offset) is 64-bit
     * in CDF-5 format, both variable's begin and size are 64-bit
     */
    int status = hdr_check_buffer(gbp, 8); /* size of int64 == 8 */
    if (status != NC_NOERR) return status;

    status = ncmpix_get_uint64((const void **)(&gbp->pos), xp);
    return status;
}

/*----< hdr_get_NC_tag() >----------------------------------------------------*/
inline static int
hdr_get_NC_tag(bufferinfo *gbp,
               NC_tag     *tagp)
{
    /* NC_tag is 4-byte integer: NC_DIMENSION, NC_VARIABLE, NC_ATTRIBUTE */
    uint type = 0;
    int status = hdr_check_buffer(gbp, 4);
    if (status != NC_NOERR) return status;

    /* get an external unsigned 4-byte integer from the file */
    status = ncmpix_get_uint32((const void**)(&gbp->pos), &type);
    if (status != NC_NOERR) return status;

    *tagp = (NC_tag) type;
    return NC_NOERR;
}

/*----< hdr_get_nc_type() >---------------------------------------------------*/
inline static int
hdr_get_nc_type(bufferinfo *gbp,
                nc_type    *typep)
{
    /* nc_type is 4-byte integer, X_SIZEOF_INT */
    int status;
    uint type;

    status = hdr_check_buffer(gbp, X_SIZEOF_INT);
    if (status != NC_NOERR) return status;

    status = ncmpix_get_uint32((const void**)(&gbp->pos), &type);
    if (status != NC_NOERR) return status;

    if (type != NC_CHAR    &&
        type != NC_BYTE    &&
        type != NC_UBYTE   &&
        type != NC_SHORT   &&
        type != NC_USHORT  &&
        type != NC_INT     &&
        type != NC_UINT    &&
        type != NC_FLOAT   &&
        type != NC_DOUBLE  &&
        type != NC_INT64   &&
        type != NC_UINT64
       )
        DEBUG_RETURN_ERROR(NC_EBADTYPE)

    *typep = (nc_type) type;
    return NC_NOERR;
}

/*----< hdr_get_NC_name() >---------------------------------------------------*/
static int
hdr_get_NC_name(bufferinfo  *gbp,
                NC_string  **ncstrpp)
{
    /* netCDF file format:
     *  ...
     * name       = nelems  namestring
     * nelems     = NON_NEG
     * namestring = ID1 [IDN ...] padding
     * ID1        = alphanumeric | '_'
     * IDN        = alphanumeric | special1 | special2
     * padding    = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    int err;
    char *cpos;
    NC_string *ncstrp;
    MPI_Aint pos_addr, base_addr;
    MPI_Offset nchars, nbytes, padding, bufremain, strcount;

    /* get nelems */
    if (gbp->version < 5) {
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        nchars = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        nchars = (MPI_Offset)tmp;
    }
    if (err != NC_NOERR) return err;

    /* Allocate a NC_string structure large enough to hold nchars characters */
    ncstrp = ncmpio_new_NC_string((size_t)nchars, NULL);
    if (ncstrp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    nbytes = nchars * X_SIZEOF_CHAR;
    padding = _RNDUP(X_SIZEOF_CHAR * ncstrp->nchars, X_ALIGN)
            - X_SIZEOF_CHAR * ncstrp->nchars;
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    bufremain = gbp->size - (pos_addr - base_addr);
    cpos = ncstrp->cp;

    /* get namestring with padding */
    while (nbytes > 0) {
        if (bufremain > 0) {
            strcount = MIN(bufremain, nbytes);
            if (strcount != (size_t)strcount)
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            memcpy(cpos, gbp->pos, (size_t)strcount);
            nbytes -= strcount;
            gbp->pos = (void *)((char *)gbp->pos + strcount);
            cpos += strcount;
            bufremain -= strcount;
        } else {
            err = hdr_fetch(gbp);
            if (err != NC_NOERR) {
                ncmpio_free_NC_string(ncstrp);
                return err;
            }
            bufremain = gbp->size;
        }
    }

    /* handle the padding */
    if (padding > 0) {
        /* CDF specification: Header padding uses null (\x00) bytes. */
        char pad[X_ALIGN-1];
        memset(pad, 0, X_ALIGN-1);
        if (memcmp(gbp->pos, pad, (size_t)padding) != 0) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header non-zero padding found\n",__FILE__,__func__,__LINE__);
#endif
            ncmpio_free_NC_string(ncstrp);
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }
        gbp->pos = (void *)((char *)gbp->pos + padding);
    }

    *ncstrpp = ncstrp;

    return NC_NOERR;
}

/*----< hdr_get_NC_dim() >----------------------------------------------------*/
inline static int
hdr_get_NC_dim(bufferinfo  *gbp,
               NC_dim     **dimpp)
{
    /* netCDF file format:
     *  ...
     * dim        = name  dim_length
     * dim_length = NON_NEG
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    int status;
    NC_string *ncstrp;
    NC_dim *dimp;

    /* get name */
    status = hdr_get_NC_name(gbp, &ncstrp);
    if (status != NC_NOERR) return status;

    dimp = ncmpio_new_x_NC_dim(ncstrp);
    if (dimp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* get dim_length */
    if (gbp->version < 5) {
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        dimp->size = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        dimp->size = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) {
        ncmpio_free_NC_dim(dimp); /* frees dim */
        return status;
    }

    *dimpp = dimp;
    return NC_NOERR;
}

/*----< hdr_get_NC_dimarray() >-----------------------------------------------*/
static int
hdr_get_NC_dimarray(bufferinfo  *gbp,
                    NC_dimarray *ncap)
{
    /* netCDF file format:
     *  ...
     * dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_DIMENSION = \x00 \x00 \x00 \x0A         // tag for list of dimensions
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i, status;
    NC_tag tag = NC_UNSPECIFIED;
    MPI_Offset ndefined;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* get NC_tag (NC_DIMENSION) */
    status = hdr_get_NC_tag(gbp, &tag);
    if (status != NC_NOERR) return status;

    /* get nelems */
    if (gbp->version < 5) {
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        ndefined = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        ndefined = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) return status;
    if (ndefined != (int)ndefined) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    ncap->ndefined = (int)ndefined;
    /* TODO: we should allow ndefined > 2^32, considering change the data type
     * of ndefined from int to MPI_Offset */

    ncap->unlimited_id = -1;

    if (ndefined == 0) {
        if (tag != NC_DIMENSION && tag != NC_UNSPECIFIED) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header format for NC_DIMENSION\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }
    } else {
        if (tag != NC_DIMENSION) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header format for NC_DIMENSION\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }

        ncap->value = (NC_dim**) NCI_Malloc((size_t)ndefined * sizeof(NC_dim*));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        ncap->nalloc = (int)ndefined;
        /* TODO: we should allow nalloc > 2^32, considering change the data
         * type of nalloc from int to MPI_Offset */

        for (i=0; i<ndefined; i++) {
            status = hdr_get_NC_dim(gbp, ncap->value + i);
            if (status != NC_NOERR) { /* error: fail to get the next dim */
                ncap->ndefined = i;
                ncmpio_free_NC_dimarray(ncap);
                return status;
            }
            if (ncap->value[i]->size == NC_UNLIMITED)
                ncap->unlimited_id = i; /* ID of unlimited dimension */
        }
    }
    return NC_NOERR;
}

/*----< hdr_get_NC_attrV() >--------------------------------------------------*/
static int
hdr_get_NC_attrV(bufferinfo *gbp,
                 NC_attr    *attrp)
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     *  ...
     * values  = bytes | chars | shorts | ints | floats | doubles
     * bytes   = [BYTE ...]  padding
     * chars   = [CHAR ...]  padding
     * shorts  = [SHORT ...]  padding
     * ints    = [INT ...]
     * floats  = [FLOAT ...]
     * doubles = [DOUBLE ...]
     * padding = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     */
    void *value = attrp->xvalue;
    MPI_Offset nbytes, padding, bufremain, attcount;
    MPI_Aint pos_addr, base_addr;

    nbytes = attrp->nelems * ncmpio_xlen_nc_type(attrp->type);
    padding = attrp->xsz - nbytes;
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    bufremain = gbp->size - (pos_addr - base_addr);

    /* get values */
    while (nbytes > 0) {
        if (bufremain > 0) {
            attcount = MIN(bufremain, nbytes);
            if (attcount != (size_t)attcount)
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            memcpy(value, gbp->pos, (size_t)attcount);
            nbytes -= attcount;
            gbp->pos = (void *)((char *)gbp->pos + attcount);
            value = (void *)((char *)value + attcount);
            bufremain -= attcount;
        } else {
            int err;
            err = hdr_fetch(gbp);
            if (err != NC_NOERR) return err;
            bufremain = gbp->size;
        }
    }

    /* handle the padding */
    if (padding > 0) {
        /* CDF specification: Header padding uses null (\x00) bytes. */
        char pad[X_ALIGN-1];
        memset(pad, 0, X_ALIGN-1);
        if (memcmp(gbp->pos, pad, (size_t)padding) != 0) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header non-zero padding found\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }
        gbp->pos = (void *)((char *)gbp->pos + padding);
    }
    return NC_NOERR;
}

/*----< hdr_get_NC_attr() >---------------------------------------------------*/
static int
hdr_get_NC_attr(bufferinfo  *gbp,
                NC_attr    **attrpp)
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     * nc_type = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * nelems  = NON_NEG       // number of elements in following sequence
     * NON_NEG = <non-negative INT> |  // CDF-1 and CDF-2
     *           <non-negative INT64>  // CDF-5
     */
    int status;
    NC_string *strp;
    nc_type type;
    MPI_Offset nelems;
    NC_attr *attrp;

    /* get name */
    status = hdr_get_NC_name(gbp, &strp);
    if (status != NC_NOERR) return status;

    /* get nc_type */
    status = hdr_get_nc_type(gbp, &type);
    if (status != NC_NOERR) {
        ncmpio_free_NC_string(strp);
        return status;
    }

    /* get nelems */
    if (gbp->version < 5) {
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        nelems = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        nelems = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) {
        ncmpio_free_NC_string(strp);
        return status;
    }

    /* allocate space for attribute object */
    attrp = ncmpio_new_x_NC_attr(strp, type, nelems);
    if (attrp == NULL) {
        ncmpio_free_NC_string(strp);
        return status;
    }

    /* get [values ...] */
    status = hdr_get_NC_attrV(gbp, attrp);
    if (status != NC_NOERR) {
        ncmpio_free_NC_attr(attrp);
        return status;
    }

    *attrpp = attrp;
    return NC_NOERR;
}

/*----< hdr_get_NC_attrarray() >----------------------------------------------*/
static int
hdr_get_NC_attrarray(bufferinfo   *gbp,
                     NC_attrarray *ncap)
{
    /* netCDF file format:
     *  ...
     * att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_ATTRIBUTE = \x00 \x00 \x00 \x0C         // tag for list of attributes
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i, status;
    NC_tag tag = NC_UNSPECIFIED;
    MPI_Offset ndefined;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* get NC_tag (NC_ATTRIBUTE) */
    status = hdr_get_NC_tag(gbp, &tag);
    if (status != NC_NOERR) return status;

    /* get nelems */
    if (gbp->version < 5) {
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        ndefined = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        ndefined = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) return status;
    if (ndefined != (int)ndefined) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    ncap->ndefined = (int)ndefined;

    if (ndefined == 0) {
        if (tag != NC_ATTRIBUTE && tag != NC_UNSPECIFIED) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header format for NC_ATTRIBUTE\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }
    } else {
        if (tag != NC_ATTRIBUTE) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header format for NC_ATTRIBUTE\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }

        ncap->value = (NC_attr**)NCI_Malloc((size_t)ndefined *sizeof(NC_attr*));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        ncap->nalloc = (int)ndefined;

        /* get [attr ...] */
        for (i=0; i<ndefined; i++) {
            status = hdr_get_NC_attr(gbp, ncap->value + i);
            if (status != NC_NOERR) { /* Error: fail to get the next att */
                ncap->ndefined = i;
                ncmpio_free_NC_attrarray(ncap);
                return status;
            }
        }
    }
    return NC_NOERR;
}

/*----< hdr_get_NC_var() >---------------------------------------------------*/
static int
hdr_get_NC_var(bufferinfo  *gbp,
               NC_var     **varpp)
{
    /* netCDF file format:
     * netcdf_file = header data
     * header      = magic numrecs dim_list gatt_list var_list
     *  ...
     * var         = name nelems [dimid ...] vatt_list nc_type vsize begin
     * nelems      = NON_NEG
     * dimid       = NON_NEG
     * vatt_list   = att_list
     * nc_type     = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * vsize       = NON_NEG
     * begin       = OFFSET        // Variable start location.
     * OFFSET      = <non-negative INT> |  // CDF-1
     *               <non-negative INT64>  // CDF-2 and CDF-5
     * NON_NEG     = <non-negative INT> |  // CDF-1 and CDF-2
     *               <non-negative INT64>  // CDF-5
     */
    int status;
    NC_string *strp;
    MPI_Offset ndims=0, dim;
    MPI_Offset tmp_dimids=0;
    NC_var *varp;

    /* get name */
    status = hdr_get_NC_name(gbp, &strp);
    if (status != NC_NOERR) return status;

    /* nelems */
    if (gbp->version < 5) {
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        ndims = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        ndims = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) {
         ncmpio_free_NC_string(strp);
         return status;
    }
    if (ndims != (int)ndims) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    /* allocate space for var object */
    varp = ncmpio_new_x_NC_var(strp, (int)ndims);
    if (varp == NULL) {
        ncmpio_free_NC_string(strp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    /* get [dimid ...] */
    for (dim=0; dim<ndims; dim++) {
        status = hdr_check_buffer(gbp, (gbp->version < 5 ? 4 : 8));
        if (status != NC_NOERR) {
            ncmpio_free_NC_var(varp);
            return status;
        }
        if (gbp->version < 5) {
            uint tmp;
            status = hdr_get_uint32(gbp, &tmp);
            tmp_dimids = (MPI_Offset)tmp;
        }
        else {
            uint64 tmp;
            status = hdr_get_uint64(gbp, &tmp);
            tmp_dimids = (MPI_Offset)tmp;
        }
        /* TODO: consider change the data type of dimids from int to
         * MPI_Offset */
        varp->dimids[dim] = (int)tmp_dimids;
        if (status != NC_NOERR) {
           return status;
        }
    }

    /* get vatt_list */
    status = hdr_get_NC_attrarray(gbp, &varp->attrs);
    if (status != NC_NOERR) {
        ncmpio_free_NC_var(varp);
        return status;
    }

    /* get nc_type */
    status = hdr_get_nc_type(gbp, &varp->type);
    if (status != NC_NOERR) {
        ncmpio_free_NC_var(varp);
        return status;
    }

    /* get vsize */
    if (gbp->version < 5) {
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        varp->len = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        varp->len = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) {
        ncmpio_free_NC_var(varp);
        return status;
    }
    /* As described in CDF-2 format specification, vsize is redundant.
       Its value may be computed from the product of dimension lengths.
       In CDF-2, vsize is a 4-byte integer. So, if we define a variable of
       less than 2^32 elements but size > 2^32-4 bytes, then vsize in CDF-2
       will overflow. Recompute varp->len can ignore an overflowed value in
       vsize stored in the file and hence bypass the limitation of CDF-2 on
       variable size of 2^32-4 bytes.

       Later on, back to ncmpio_hdr_get_NC(), compute_var_shape() is
       called which recomputes varp->len using the dimension values and hence
       overwrites the value read from file above.

       In summary, PnetCDF now ignores the value of vsize stored in the file
       header.
     */

    /* next element is 'begin' */
    status = hdr_check_buffer(gbp, (gbp->version == 1 ? 4 : 8));
    if (status != NC_NOERR) {
        ncmpio_free_NC_var(varp);
        return status;
    }

    /* get begin */
    if (gbp->version == 1) {
        uint tmp=0;
        status = ncmpix_get_uint32((const void **)(&gbp->pos), &tmp);
        varp->begin = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp=0;
        status = ncmpix_get_uint64((const void **)(&gbp->pos), &tmp);
        varp->begin = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) {
        ncmpio_free_NC_var(varp);
        return status;
    }

    *varpp = varp;
    return NC_NOERR;
}

/*----< hdr_get_NC_vararray() >----------------------------------------------*/
static int
hdr_get_NC_vararray(bufferinfo  *gbp,
                    NC_vararray *ncap)
{
    /* netCDF file format:
     * netcdf_file = header  data
     * header      = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var_list    = ABSENT | NC_VARIABLE   nelems  [var ...]
     * ABSENT      = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *               ZERO  ZERO64  // for CDF-5
     * ZERO        = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64      = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_VARIABLE = \x00 \x00 \x00 \x0B         // tag for list of variables
     * nelems      = NON_NEG       // number of elements in following sequence
     * NON_NEG     = <non-negative INT> |        // CDF-1 and CDF-2
     *               <non-negative INT64>        // CDF-5
     */
    int i, status;
    NC_tag tag = NC_UNSPECIFIED;
    MPI_Offset ndefined;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* get NC_tag (NC_VARIABLE) from gbp buffer */
    status = hdr_get_NC_tag(gbp, &tag);
    if (status != NC_NOERR) return status;

    /* get nelems (number of variables) from gbp buffer */
    if (gbp->version < 5) {
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        ndefined = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        ndefined = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) return status;
    if (ndefined != (int)ndefined) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    ncap->ndefined = (int)ndefined;
    /* TODO: we should allow ndefined > 2^32, considering change the data type
     * of ndefined from int to MPI_Offset */

    if (ndefined == 0) { /* no variable defined */
        if (tag != NC_VARIABLE && tag != NC_UNSPECIFIED) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header format for NC_VARIABLE\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }
    } else {
        if (tag != NC_VARIABLE) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header format for NC_VARIABLE\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }

        ncap->value = (NC_var**) NCI_Malloc((size_t)ndefined * sizeof(NC_var*));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        ncap->nalloc = (int)ndefined;

        /* get [var ...] */
        for (i=0; i<ndefined; i++) {
            status = hdr_get_NC_var(gbp, ncap->value + i);
            ncap->value[i]->varid = i;
            if (status != NC_NOERR) { /* Error: fail to get the next var */
                ncap->ndefined = i;
                ncmpio_free_NC_vararray(ncap);
                return status;
            }
        }
    }

    return NC_NOERR;
}

/*----< ncmpio_hdr_len_NC() >------------------------------------------------*/
MPI_Offset
ncmpio_hdr_len_NC(const NC *ncp)
{
    /* netCDF file format:
     * netcdf_file = header  data
     * header      = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * numrecs     = NON_NEG | STREAMING   // length of record dimension
     * NON_NEG     = <non-negative INT> |  // CDF-1 and CDF-2
     *               <non-negative INT64>  // CDF-5
     */

    int sizeof_NON_NEG, sizeof_off_t;
    MPI_Offset xlen;

    assert(ncp != NULL);

    if (ncp->format == 5) {        /* CDF-5 */
        sizeof_NON_NEG = X_SIZEOF_INT64; /* 8-byte integer for all integers */
        sizeof_off_t   = X_SIZEOF_INT64; /* 8-byte integer for var begin */
    }
    else if (ncp->format == 2) { /* CDF-2 */
        sizeof_NON_NEG = X_SIZEOF_INT; /* 4-byte integer in CDF-1 */
        sizeof_off_t   = X_SIZEOF_INT64; /* 8-byte integer for var begin */
    }
    else { /* CDF-1 */
        sizeof_NON_NEG = X_SIZEOF_INT; /* 4-byte integer in CDF-1 */
        sizeof_off_t   = X_SIZEOF_INT; /* 4-byte integer in CDF-1 */
    }

    xlen  = sizeof(ncmagic1);                                                /* magic */
    xlen += sizeof_NON_NEG;                                                  /* numrecs */
    xlen += hdr_len_NC_dimarray(&ncp->dims,   sizeof_NON_NEG);               /* dim_list */
    xlen += hdr_len_NC_attrarray(&ncp->attrs, sizeof_NON_NEG);               /* gatt_list */
    xlen += hdr_len_NC_vararray(&ncp->vars,   sizeof_NON_NEG, sizeof_off_t); /* var_list */

    return xlen; /* return the header size (not yet aligned) */
}

/*----< ncmpio_hdr_get_NC() >------------------------------------------------*/
/*  CDF format specification
 *      netcdf_file  = header  data
 *      header       = magic  numrecs  dim_list  gatt_list  var_list
 *      magic        = 'C'  'D'  'F'  VERSION
 *      VERSION      = \x01 |                      // classic format
 *                     \x02 |                      // 64-bit offset format
 *                     \x05                        // 64-bit data format
 *      numrecs      = NON_NEG | STREAMING         // length of record dimension
 *      dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
 *      gatt_list    = att_list                    // global attributes
 *      att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
 *      var_list     = ABSENT | NC_VARIABLE   nelems  [var ...]
 */
int
ncmpio_hdr_get_NC(NC *ncp)
{
    int status;
    bufferinfo getbuf;
    char magic[sizeof(ncmagic1)];
    MPI_Offset nrecs = 0;
    MPI_Aint pos_addr, base_addr;

    assert(ncp != NULL);

    /* Initialize the get buffer that stores the header read from the file */
    getbuf.comm          = ncp->comm;
    getbuf.collective_fh = ncp->collective_fh;
    getbuf.get_size      = 0;
    getbuf.offset        = 0;   /* read from start of the file */
    getbuf.safe_mode     = ncp->safe_mode;

    /* CDF-5's minimum header size is 4 bytes more than CDF-1 and CDF-2's */
    getbuf.size = _RNDUP( MAX(MIN_NC_XSZ+4, ncp->chunk), X_ALIGN );

    if (getbuf.size != (size_t)getbuf.size) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    getbuf.pos = getbuf.base = (void *)NCI_Malloc((size_t)getbuf.size);

    /* Fetch the next header chunk. The chunk is 'gbp->size' bytes big */
    status = hdr_fetch(&getbuf);
    if (status != NC_NOERR) return status;

    /* processing the header from getbuf, the get buffer */

    /* First get the file format information, magic */
    memset(magic, 0, sizeof(magic));
    status = ncmpix_getn_text((const void **)(&getbuf.pos), sizeof(magic),
                              magic);
    if (status != NC_NOERR) return status;

    /* check if the first three bytes are 'C','D','F' */
    if (memcmp(magic, ncmagic1, sizeof(ncmagic1)-1) != 0) {
        /* check if is HDF5 file */
        char signature[8], *hdf5_signature="\211HDF\r\n\032\n";
        memcpy(signature, magic, 4);
        ncmpix_getn_text((const void **)(&getbuf.pos), 4, signature+4);
        if (memcmp(signature, hdf5_signature, 8) == 0) {
            DEBUG_ASSIGN_ERROR(status, NC_ENOTNC3)
            if (ncp->safe_mode)
                fprintf(stderr,"Error: file %s is HDF5 format\n",ncp->path);
        }
        else
            DEBUG_ASSIGN_ERROR(status, NC_ENOTNC)
        goto fn_exit;
    }

    /* check version number in last byte of magic */
    if (magic[sizeof(ncmagic1)-1] == 0x1) {
        getbuf.version = ncp->format = 1;
    } else if (magic[sizeof(ncmagic1)-1] == 0x2) {
        getbuf.version = ncp->format = 2;
#if SIZEOF_MPI_OFFSET < 8
        /* take the easy way out: if we can't support all CDF-2
         * files, return immediately */
        NCI_Free(getbuf.base);
        DEBUG_RETURN_ERROR(NC_ESMALL)
#endif
    } else if (magic[sizeof(ncmagic1)-1] == 0x5) {
        getbuf.version = ncp->format = 5;
#if SIZEOF_MPI_OFFSET < 8
        NCI_Free(getbuf.base);
        DEBUG_RETURN_ERROR(NC_ESMALL)
#endif
    } else {
        NCI_Free(getbuf.base);
        DEBUG_RETURN_ERROR(NC_ENOTNC) /* not a netCDF file */
    }

    /** Ensure that 'nextread' bytes (numrecs) are available. */
    status = hdr_check_buffer(&getbuf, (getbuf.version < 5) ? 4 : 8);
    if (status != NC_NOERR) goto fn_exit;

    /* get numrecs from getbuf into ncp */
    if (getbuf.version < 5) {
        uint tmp=0;
        status = ncmpix_get_uint32((const void **)(&getbuf.pos), &tmp);
        nrecs = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp=0;
        status = ncmpix_get_uint64((const void **)(&getbuf.pos), &tmp);
        nrecs = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) goto fn_exit;

    ncp->numrecs = nrecs;

#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(getbuf.pos,  &pos_addr);
    MPI_Get_address(getbuf.base, &base_addr);
#else
    MPI_Address(getbuf.pos,  &pos_addr);
    MPI_Address(getbuf.base, &base_addr);
#endif
    assert(pos_addr < base_addr + getbuf.size);

    /* get dim_list from getbuf into ncp */
    status = hdr_get_NC_dimarray(&getbuf, &ncp->dims);
    if (status != NC_NOERR) goto fn_exit;

    /* get gatt_list from getbuf into ncp */
    status = hdr_get_NC_attrarray(&getbuf, &ncp->attrs);
    if (status != NC_NOERR) goto fn_exit;

    /* get var_list from getbuf into ncp */
    status = hdr_get_NC_vararray(&getbuf, &ncp->vars);
    if (status != NC_NOERR) goto fn_exit;

    /* get the un-aligned size occupied by the file header */
    ncp->xsz = ncmpio_hdr_len_NC(ncp);

    /* Recompute the shapes of all variables
     * Sets ncp->begin_var to start of first variable.
     * Sets ncp->begin_rec to start of first record variable.
     */
    status = compute_var_shape(ncp);
    if (status != NC_NOERR) goto fn_exit;

    status = ncmpio_NC_check_vlens(ncp);
    if (status != NC_NOERR) goto fn_exit;

fn_exit:
    ncp->get_size += getbuf.get_size;
    NCI_Free(getbuf.base);

    return status;
}

/*----< ncmpio_write_header() >---------------------------------------------*/
/* This function is collective (even in independent data mode).
 * It is called only in data mode (collective or independent) and by
 * 1. ncmpi_rename_att()
 * 2. ncmpi_copy_att()
 * 3. ncmpi_put_att()
 * 4. ncmpi_rename_dim()
 * 5. ncmpi_rename_var()
 */
int ncmpio_write_header(NC *ncp)
{
    int rank, status=NC_NOERR, mpireturn, err;
    MPI_File fh;

    /* Write the entire header to the file. This function may be called from
     * a rename API. In that case, we cannot just change the variable name in
     * the file header, because if the file space occupied by the name shrinks,
     * all metadata following the new name must be moved ahead.
     */

    fh = ncp->collective_fh;
    if (NC_indep(ncp))
        fh = ncp->independent_fh;

    /* update file header size, as this subroutine may be called from a rename
     * API (var or attribute) and the new name is smaller/bigger which changes
     * the header size. We recalculate ncp->xsz by getting the un-aligned size
     * occupied by the file header */
    ncp->xsz = ncmpio_hdr_len_NC(ncp);

    MPI_Comm_rank(ncp->comm, &rank);
    if (rank == 0) { /* only root writes to file header */
        MPI_Status mpistatus;
        void *buf = NCI_Malloc((size_t)ncp->xsz); /* header's write buffer */

        /* copy header object to write buffer */
        status = ncmpio_hdr_put_NC(ncp, buf);

        if (ncp->xsz != (int)ncp->xsz) {
            NCI_Free(buf);
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
        TRACE_IO(MPI_File_write_at)(fh, 0, buf, (int)ncp->xsz, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpio_handle_error(mpireturn, "MPI_File_write_at");
            if (status == NC_NOERR) {
                err = (err == NC_EFILE) ? NC_EWRITE : err;
                DEBUG_ASSIGN_ERROR(status, err)
            }
        }
        else {
            ncp->put_size += ncp->xsz;
        }
        NCI_Free(buf);
    }

    if (ncp->safe_mode == 1) {
        /* broadcast root's status, because only root writes to the file */
        int root_status = status;
        TRACE_COMM(MPI_Bcast)(&root_status, 1, MPI_INT, 0, ncp->comm);
        /* root's write has failed, which is serious */
        if (root_status == NC_EWRITE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
        if (mpireturn != MPI_SUCCESS) {
            ncmpio_handle_error(mpireturn,"MPI_Bcast");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
    }

    if (NC_doFsync(ncp)) { /* NC_SHARE is set */
        TRACE_IO(MPI_File_sync)(fh);
        if (mpireturn != MPI_SUCCESS) {
            ncmpio_handle_error(mpireturn,"MPI_File_sync");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
        TRACE_COMM(MPI_Barrier)(ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            ncmpio_handle_error(mpireturn,"MPI_Barrier");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
    }

    return status;
}

