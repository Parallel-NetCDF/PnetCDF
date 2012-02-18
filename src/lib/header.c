/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#include "nc.h"

#include <mpi.h>
#include <assert.h>
#include <string.h>  /* memcpy() */
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>

#include "ncx.h"
#include "macro.h"

#ifdef SIZEOF_INT
# if SIZEOF_INT == 4
#  define lld(x) (x)
# elif  SIZEOF_INT == 8
#  define lld(x) (long long)(x)
# endif
#endif

/* Prototypes for functions used only in this file */
static MPI_Offset hdr_len_NC_name(const NC_string *ncstrp, int sizeof_t);
static MPI_Offset hdr_len_NC_dim(const NC_dim *dimp, int sizeof_t);
static MPI_Offset hdr_len_NC_dimarray(const NC_dimarray *ncap, int sizeof_t);
static MPI_Offset hdr_len_NC_attr(const NC_attr *attrp, int sizeof_t);
static MPI_Offset hdr_len_NC_attrarray(const NC_attrarray *ncap, int sizeof_t);
static MPI_Offset hdr_len_NC_var(const NC_var *varp, int sizeof_off_t, int sizeof_t);
static MPI_Offset hdr_len_NC_vararray(const NC_vararray *ncap, int sizeof_t, int sizeof_off_t);
static int hdr_put_NCtype(bufferinfo *pbp, NCtype type);
static int hdr_put_nc_type(bufferinfo *pbp, const nc_type *typep);
static int hdr_put_NC_name(bufferinfo *pbp, const NC_string *ncstrp);
static int hdr_put_NC_attrV(bufferinfo *pbp, const NC_attr *attrp);
static int hdr_put_NC_dim(bufferinfo *pbp, const NC_dim *dimp);
static int hdr_put_NC_attr(bufferinfo *pbp, const NC_attr *attrp);
static int hdr_put_NC_var(bufferinfo *pbp, const NC_var *varp);
static int hdr_put_NC_dimarray(bufferinfo *pbp, const NC_dimarray *ncap);
static int hdr_put_NC_attrarray(bufferinfo *pbp, const NC_attrarray *ncap);
static int hdr_put_NC_vararray(bufferinfo *pbp, const NC_vararray *ncap);
static int hdr_fetch(bufferinfo *gbp);
static int hdr_check_buffer(bufferinfo *gbp, MPI_Offset nextread);
static int hdr_get_NCtype(bufferinfo *gbp, NCtype *typep);
static int hdr_get_size_t(bufferinfo *gbp, MPI_Offset *sp);
static int hdr_get_NC_string(bufferinfo *gbp, NC_string **ncstrpp);
static int hdr_get_NC_dim(bufferinfo *gbp, NC_dim **dimpp);
static int hdr_get_NC_dimarray(bufferinfo *gbp, NC_dimarray *ncap);
static int hdr_get_nc_type(bufferinfo *gbp, nc_type *typep);
static int hdr_get_NC_attrV(bufferinfo *gbp, NC_attr *attrp);
static int hdr_get_NC_attr(bufferinfo *gbp, NC_attr **attrpp);
static int hdr_get_NC_attrarray(bufferinfo *gbp, NC_attrarray *ncap);
static int hdr_get_NC_var(bufferinfo *gbp, NC_var **varpp);
static int hdr_get_NC_vararray(bufferinfo *gbp, NC_vararray *ncap);

/*
 * "magic number" at beginning of file: 0x43444601 (big endian) 
 */
static const schar ncmagic[]  = {'C', 'D', 'F', 0x02}; 
static const schar ncmagic1[] = {'C', 'D', 'F', 0x01}; 
static const schar ncmagic2[] = {'C', 'D', 'F', 0x05}; 

/*
 * Recompute the shapes of all variables
 * Sets ncp->begin_var to start of first variable.
 * Sets ncp->begin_rec to start of first record variable.
 * Returns -1 on error. The only possible error is an reference
 * to a non existent dimension, which would occur for a corrupt
 * netcdf file.
 */ 
int
ncmpii_NC_computeshapes(NC *ncp)
{
    NC_var **vpp = (NC_var **)ncp->vars.value;
    NC_var *const *const end = &vpp[ncp->vars.ndefined];
    NC_var *first_var = NULL;       /* first "non-record" var */
    NC_var *first_rec = NULL;       /* first "record" var */
    int status;
 
    ncp->begin_var = (MPI_Offset) ncp->xsz;
    ncp->begin_rec = (MPI_Offset) ncp->xsz;
    ncp->recsize = 0;
 
    if (ncp->vars.ndefined == 0) return NC_NOERR;

    for ( /*NADA*/; vpp < end; vpp++) {    
        status = ncmpii_NC_var_shape64(ncp, *vpp, &ncp->dims);

        if (status != NC_NOERR) return status ;
 
        if (IS_RECVAR(*vpp)) {
            if (first_rec == NULL)
                first_rec = *vpp;
            ncp->recsize += (*vpp)->len;
        }
        else {
            if (first_var == NULL) 
            first_var = *vpp;
            /*
             * Overwritten each time thru.
             * Usually overwritten in first_rec != NULL clause.
             */
            ncp->begin_rec = (*vpp)->begin + (MPI_Offset)(*vpp)->len;
        }
    }
 
    if (first_rec != NULL) {
        assert(ncp->begin_rec <= first_rec->begin);
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

    assert(ncp->begin_var > 0);
    assert(ncp->xsz <= ncp->begin_var);
    assert(ncp->begin_rec > 0);
    assert(ncp->begin_var <= ncp->begin_rec);
 
    return NC_NOERR;
} 

/* 
 * To compute how much space will the xdr'd header take 
 */
 
#define X_SIZEOF_NC_TYPE X_SIZEOF_INT
#define X_SIZEOF_NCTYPE X_SIZEOF_INT

/*----< hdr_len_NC_name() >--------------------------------------------------*/
static MPI_Offset
hdr_len_NC_name(const NC_string *ncstrp,
                int              sizeof_t)     /* NON_NEG */
{
    /* netCDF file format:
     * name         = nelems  namestring
     * nelems       = NON_NEG
     * namestring   = ID1 [IDN ...] padding
     * ID1          = alphanumeric | '_'
     * IDN          = alphanumeric | special1 | special2
     * padding      = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz = sizeof_t; /* nelems */
 
    assert(ncstrp != NULL);
 
    if (ncstrp->nchars != 0)  /* namestring */
        sz += _RNDUP(ncstrp->nchars, X_ALIGN);
 
    return sz;
}
 
/*----< hdr_len_NC_dim() >---------------------------------------------------*/
static MPI_Offset
hdr_len_NC_dim(const NC_dim *dimp,
               int           sizeof_t)     /* NON_NEG */
{
    /* netCDF file format:
     *  ...
     * dim          = name  dim_length
     * dim_length   = NON_NEG
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;
 
    assert(dimp != NULL);
 
    sz = hdr_len_NC_name(dimp->name, sizeof_t); /* name */
    sz += sizeof_t;                             /* dim_length */
 
    return sz;
}
 
/*----< hdr_len_NC_dimarray() >----------------------------------------------*/
static MPI_Offset
hdr_len_NC_dimarray(const NC_dimarray *ncap,
                    int                sizeof_t)     /* NON_NEG */
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
    MPI_Offset xlen;

    xlen = X_SIZEOF_NCTYPE;           /* NC_DIMENSION */
    xlen += sizeof_t;                 /* nelems */

    if (ncap == NULL) /* ABSENT: no dimension is defined */
        return xlen;

    const NC_dim **dpp = (const NC_dim **)ncap->value;
    const NC_dim *const *const end = &dpp[ncap->ndefined];
    for (/*NADA*/; dpp < end; dpp++)  /* [dim ...] */
        xlen += hdr_len_NC_dim(*dpp, sizeof_t);

    return xlen;
} 

/*----< hdr_len_NC_attr() >--------------------------------------------------*/
static MPI_Offset
hdr_len_NC_attr(const NC_attr *attrp,
                int            sizeof_t)     /* NON_NEG */
{
    /* netCDF file format:
     *  ...
     * attr         = name  nc_type  nelems  [values ...]
     * nc_type      = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * nelems       = NON_NEG       // number of elements in following sequence
     * values       = bytes | chars | shorts | ints | floats | doubles
     * bytes        = [BYTE ...]  padding
     * chars        = [CHAR ...]  padding
     * shorts       = [SHORT ...]  padding
     * ints         = [INT ...]
     * floats       = [FLOAT ...]
     * doubles      = [DOUBLE ...]
     * padding      = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;
 
    assert(attrp != NULL);
 
    sz  = hdr_len_NC_name(attrp->name, sizeof_t); /* name */
    sz += X_SIZEOF_NC_TYPE;                       /* nc_type */
    sz += sizeof_t;                               /* nelems */
    sz += attrp->xsz;                             /* [values ...] */
 
    return sz;
}
 
/*----< hdr_len_NC_attrarray() >---------------------------------------------*/
static MPI_Offset
hdr_len_NC_attrarray(const NC_attrarray *ncap,
                     int                 sizeof_t)     /* NON_NEG */
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

    xlen = X_SIZEOF_NCTYPE;        /* NC_ATTRIBUTE */
    xlen += sizeof_t;              /* nelems */

    if (ncap == NULL) /* ABSENT: no attribute is defined */
        return xlen;

    for (i=0; i<ncap->ndefined; i++) /* [attr ...] */
        xlen += hdr_len_NC_attr(ncap->value[i], sizeof_t);

    return xlen;
}
 
/*----< hdr_len_NC_var() >---------------------------------------------------*/
static MPI_Offset
hdr_len_NC_var(const NC_var *varp,
               int           sizeof_off_t, /* OFFSET */
               int           sizeof_t)     /* NON_NEG */
{
    /* netCDF file format:
     * netcdf_file  = header  data
     * header       = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var          = name  nelems  [dimid ...]  vatt_list  nc_type  vsize  begin
     * nelems       = NON_NEG
     * dimid        = NON_NEG
     * vatt_list    = att_list
     * nc_type      = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * vsize        = NON_NEG
     * begin        = OFFSET        // Variable start location.
     * OFFSET       = <non-negative INT> |  // CDF-1
     *                <non-negative INT64>  // CDF-2 and CDF-5
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;
 
    assert(varp != NULL);
 
    /* for CDF-1, sizeof_off_t == 4 && sizeof_t == 4
     * for CDF-2, sizeof_off_t == 8 && sizeof_t == 4
     * for CDF-5, sizeof_off_t == 8 && sizeof_t == 8
     */
    sz = hdr_len_NC_name(varp->name, sizeof_t);         /* name */
    sz += sizeof_t;                                     /* nelems */
    sz += sizeof_t * varp->ndims;                       /* [dimid ...] */
    sz += hdr_len_NC_attrarray(&varp->attrs, sizeof_t); /* vatt_list */
    sz += X_SIZEOF_NC_TYPE;                             /* nc_type */
    sz += sizeof_t;                                     /* vsize */
    sz += sizeof_off_t;                                 /* begin */
 
    return sz;
} 

/*----< hdr_len_NC_vararray() >----------------------------------------------*/
static MPI_Offset
hdr_len_NC_vararray(const NC_vararray *ncap,
                    int                sizeof_t,     /* NON_NEG */
                    int                sizeof_off_t) /* OFFSET */
{
    /* netCDF file format:
     * netcdf_file  = header  data
     * header       = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var_list     = ABSENT | NC_VARIABLE   nelems  [var ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_VARIABLE  = \x00 \x00 \x00 \x0B         // tag for list of variables
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i;
    MPI_Offset xlen;

    xlen = X_SIZEOF_NCTYPE;           /* NC_VARIABLE */
    xlen += sizeof_t;                 /* nelems */

    if (ncap == NULL) /* ABSENT: no variable is defined */
        return xlen;

    /* for CDF-1, sizeof_off_t == 4 && sizeof_t == 4
     * for CDF-2, sizeof_off_t == 8 && sizeof_t == 4
     * for CDF-5, sizeof_off_t == 8 && sizeof_t == 8
     */
    for (i=0; i<ncap->ndefined; i++)  /* [var ...] */
        xlen += hdr_len_NC_var(ncap->value[i], sizeof_off_t, sizeof_t);

    return xlen;
}
 
/*----< ncmpii_hdr_len_NC() >------------------------------------------------*/
MPI_Offset
ncmpii_hdr_len_NC(const NC *ncp)
{
    /* netCDF file format:
     * netcdf_file  = header  data
     * header       = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * numrecs      = NON_NEG | STREAMING   // length of record dimension
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */

    int sizeof_t, sizeof_off_t;
    MPI_Offset xlen;
 
    assert(ncp != NULL);

    sizeof_t     = X_SIZEOF_INT; /* 4-byte integer in CDF-1 */
    sizeof_off_t = X_SIZEOF_INT; /* 4-byte integer in CDF-1 */

    if (fIsSet(ncp->flags, NC_64BIT_DATA)) {      /* CDF-5 */
        sizeof_t     = 8; /* 8-byte integer for all integers */
        sizeof_off_t = 8; /* 8-byte integer for var begin */
    }
    else if (fIsSet(ncp->flags, NC_64BIT_OFFSET)) /* CDF-2 */
        sizeof_off_t = 8; /* 8-byte integer for var begin */

    xlen  = sizeof(ncmagic);                                           /* magic */
    xlen += sizeof_t;                                                  /* numrecs */
    xlen += hdr_len_NC_dimarray(&ncp->dims,   sizeof_t);               /* dim_list */
    xlen += hdr_len_NC_attrarray(&ncp->attrs, sizeof_t);               /* gatt_list */
    xlen += hdr_len_NC_vararray(&ncp->vars,   sizeof_t, sizeof_off_t); /* var_list */

    return xlen; /* return the header size (not yet aligned) */
} 

/* Begin Of put NC */

/* Begin Of get NC */

/*----< hdr_put_NCtype() >---------------------------------------------------*/
static int
hdr_put_NCtype(bufferinfo *pbp, NCtype type) {
    int status;
    const int itype = (int)type;

    status = ncmpix_put_int_int(pbp->pos, &itype);
    if (status != NC_NOERR) return status;

    pbp->pos = (void *)((char *)pbp->pos + X_SIZEOF_INT); 

    return status;
}

/*----< hdr_put_nc_type() >--------------------------------------------------*/
static int 
hdr_put_nc_type(bufferinfo    *pbp,
                const nc_type *typep)
{
    int status;
    const int itype = (int) *typep;
  
    status =  ncmpix_put_int_int(pbp->pos, &itype);
    if (status != NC_NOERR) return status;

    pbp->pos = (void *)((char *)pbp->pos + X_SIZEOF_INT);
  
    return status; 
}

/*----< hdr_put_NC_name() >--------------------------------------------------*/
static int 
hdr_put_NC_name(bufferinfo      *pbp,
                const NC_string *ncstrp)
{
    /* netCDF file format:
     *  ...
     * name         = nelems  namestring
     * nelems       = NON_NEG
     * namestring   = ID1 [IDN ...] padding
     * ID1          = alphanumeric | '_'
     * IDN          = alphanumeric | special1 | special2
     * padding      = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */
    int status;

    /* copy nelems */
    status = ncmpix_put_size_t(&pbp->pos, ncstrp->nchars, pbp->version == 5 ? 8 : 4);
    if (status != NC_NOERR) return status;

    /* copy namestring */
    status = ncmpix_pad_putn_text(&pbp->pos, ncstrp->nchars, ncstrp->cp); 
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

/*----< hdr_put_NC_attrV() >-------------------------------------------------*/
/*
 * Put the values of an attribute 
 */
static int
hdr_put_NC_attrV(bufferinfo    *pbp,
                 const NC_attr *attrp)
{
    /* netCDF file format:
     *  ...
     * attr         = name  nc_type  nelems  [values ...]
     *  ...
     * values       = bytes | chars | shorts | ints | floats | doubles
     * bytes        = [BYTE ...]  padding
     * chars        = [CHAR ...]  padding
     * shorts       = [SHORT ...]  padding
     * ints         = [INT ...]
     * floats       = [FLOAT ...]
     * doubles      = [DOUBLE ...]
     * padding      = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     */
    void *value = attrp->xvalue;
    MPI_Offset padding, esz;

    /* esz is the element size (unaligned)
       attrp->xsz is an aligned element size (for byte and short)
     */
    esz = ncmpix_len_nctype(attrp->type);
    padding = attrp->xsz - esz * attrp->nelems;

    (void) memcpy(pbp->pos, value, esz * attrp->nelems);
    pbp->pos = (void *)((char *)pbp->pos + esz * attrp->nelems);

    /* padding is per buffer, not per element */
    memset(pbp->pos, 0, padding);
    pbp->pos = (void *)((char *)pbp->pos + padding);
    
    return NC_NOERR;
}

/*----< hdr_put_NC_dim() >---------------------------------------------------*/
static int 
hdr_put_NC_dim(bufferinfo   *pbp,
               const NC_dim *dimp)
{
    /* netCDF file format:
     *  ...
     * dim          = name  dim_length
     * dim_length   = NON_NEG
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */
    int status; 

    /* copy name */
    status = hdr_put_NC_name(pbp, dimp->name);
    if (status != NC_NOERR) return status;
 
    /* copy dim_length */
    status = ncmpix_put_size_t(&pbp->pos, dimp->size, pbp->version == 5 ? 8 : 4);
    if (status != NC_NOERR) return status;
  
    return NC_NOERR;
}

/*----< hdr_put_NC_attr() >--------------------------------------------------*/
static int
hdr_put_NC_attr(bufferinfo    *pbp,
                const NC_attr *attrp)
{
    /* netCDF file format:
     *  ...
     * attr         = name  nc_type  nelems  [values ...]
     * nc_type      = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */
    int status;

    /* copy name */
    status = hdr_put_NC_name(pbp, attrp->name);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = hdr_put_nc_type(pbp, &attrp->type);
    if (status != NC_NOERR) return status;

    /* copy nelems */
    status = ncmpix_put_size_t(&pbp->pos, attrp->nelems, pbp->version == 5 ? 8 : 4);
    if (status != NC_NOERR) return status;

    /* copy [values ...] */
    status = hdr_put_NC_attrV(pbp, attrp);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

/*----< hdr_put_NC_var() >---------------------------------------------------*/
static int
hdr_put_NC_var(bufferinfo   *pbp,
               const NC_var *varp)
{
    /* netCDF file format:
     * netcdf_file  = header  data
     * header       = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var          = name  nelems  [dimid ...]  vatt_list  nc_type  vsize  begin
     * nelems       = NON_NEG
     * dimid        = NON_NEG
     * vatt_list    = att_list
     * nc_type      = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * vsize        = NON_NEG
     * begin        = OFFSET        // Variable start location.
     * OFFSET       = <non-negative INT> |  // CDF-1
     *                <non-negative INT64>  // CDF-2 and CDF-5
     * NON_NEG      = <non-negative INT> |  // CDF-1 and CDF-2
     *                <non-negative INT64>  // CDF-5
     */
    int i, status, sizeof_t, sizeof_off_t;

    sizeof_t     = pbp->version == 5 ? 8 : 4;  /* for vsize */
    sizeof_off_t = pbp->version == 1 ? 4 : 8;  /* for begin */

    /* copy name */
    status = hdr_put_NC_name(pbp, varp->name);
    if (status != NC_NOERR) return status;

    /* copy nelems */
    status = ncmpix_put_size_t(&pbp->pos, varp->ndims, sizeof_t);
    if (status != NC_NOERR) return status;

    /* copy [dimid ...] */
    for (i=0; i<varp->ndims; i++) {
        status = ncmpix_put_size_t(&pbp->pos, varp->dimids[i], sizeof_t);
        if (status != NC_NOERR) return status;
    }

    /* copy vatt_list */
    status = hdr_put_NC_attrarray(pbp, &varp->attrs);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = hdr_put_nc_type(pbp, &varp->type); 
    if (status != NC_NOERR) return status;

    /* copy vsize */
    /* in CDF-1 and CDF-2, a variable's size in the header is a 32-bit integer
     * in CDF-5, it is a 64-bit integer
     */
    status = ncmpix_put_size_t(&pbp->pos, varp->len, sizeof_t);
    if (status != NC_NOERR) return status;

    /* copy begin */
    /* in CDF-1, a variable's starting file offset in the header is a 32-bit integer
     * in CDF-2 and CDF-5, it is a 64-bit integer
     */
    status = ncmpix_put_size_t(&pbp->pos, varp->begin, sizeof_off_t);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

/*----< hdr_put_NC_dimarray() >----------------------------------------------*/
static int 
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
    int status;
  
    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = hdr_put_NCtype(pbp, NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        status = ncmpix_put_size_t(&pbp->pos, 0, pbp->version == 5 ? 8 : 4);
        if (status != NC_NOERR) return status;
    }
    else {
        const NC_dim **dpp = (const NC_dim **)ncap->value; 
        const NC_dim *const *const end = &dpp[ncap->ndefined]; 

        /* copy NC_DIMENSION */
        status = hdr_put_NCtype(pbp, NC_DIMENSION);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        status = ncmpix_put_size_t(&pbp->pos, ncap->ndefined, pbp->version == 5 ? 8 : 4);
        if (status != NC_NOERR) return status;

        /* copy [dim ...] */
        for ( /*NADA*/; dpp < end; dpp++) {
            status = hdr_put_NC_dim(pbp, *dpp);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< hdr_put_NC_attrarray() >---------------------------------------------*/
static int
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
    int status;

    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = hdr_put_NCtype(pbp, NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        status = ncmpix_put_size_t(&pbp->pos, 0, pbp->version == 5 ? 8 : 4);
        if (status != NC_NOERR) return status;
    }
    else {
        const NC_attr **app = (const NC_attr **)ncap->value;
        const NC_attr *const *const end = &app[ncap->ndefined]; 

        /* copy NC_ATTRIBUTE */
        status = hdr_put_NCtype(pbp, NC_ATTRIBUTE);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        status = ncmpix_put_size_t(&pbp->pos, ncap->ndefined, pbp->version == 5 ? 8 : 4);
        if (status != NC_NOERR) return status;
  
        /* copy [attr ...] */
        for ( /*NADA*/; app < end; app++) {
            status = hdr_put_NC_attr(pbp, *app);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< hdr_put_NC_vararray() >----------------------------------------------*/
static int
hdr_put_NC_vararray(bufferinfo        *pbp,
                    const NC_vararray *ncap)
{
    /* netCDF file format:
     * netcdf_file  = header  data
     * header       = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var_list     = ABSENT | NC_VARIABLE   nelems  [var ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_VARIABLE  = \x00 \x00 \x00 \x0B         // tag for list of variables
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int status;

    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = hdr_put_NCtype(pbp, NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        status = ncmpix_put_size_t(&pbp->pos, 0, pbp->version == 5 ? 8 : 4);
        if (status != NC_NOERR) return status;
    }
    else { 
        const NC_var **vpp = (const NC_var **)ncap->value; 
        const NC_var *const *const end = &vpp[ncap->ndefined]; 

        /* copy NC_VARIABLE */
        status = hdr_put_NCtype(pbp, NC_VARIABLE);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        status = ncmpix_put_size_t(&pbp->pos, ncap->ndefined, pbp->version == 5 ? 8 : 4);
        if (status != NC_NOERR) return status;

        /* copy [var ...] */
        for (/*NADA*/; vpp < end; vpp++) {
            status = hdr_put_NC_var(pbp, *vpp); 
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< ncmpii_hdr_put_NC() >------------------------------------------------*/
int 
ncmpii_hdr_put_NC(NC   *ncp,
                  void *buf) {
    int status;
    bufferinfo putbuf;
    MPI_Offset nrecs=0; 

    putbuf.nciop  = NULL;
    putbuf.offset = 0;
    putbuf.pos    = buf;
    putbuf.base   = buf;
    putbuf.size   = ncp->xsz;

    /* netCDF file format:
     * netcdf_file  = header  data
     * header       = magic  numrecs  dim_list  gatt_list  var_list
     */

    /* copy magic */
    if (ncp->flags & NC_64BIT_DATA) {
        putbuf.version = 5;
        status = ncmpix_putn_schar_schar(&putbuf.pos, sizeof(ncmagic2), ncmagic2);
    }
    else if (ncp->flags & NC_64BIT_OFFSET) {
        putbuf.version = 2;
        status = ncmpix_putn_schar_schar(&putbuf.pos, sizeof(ncmagic), ncmagic);
    }
    else {
        putbuf.version = 1;
        status = ncmpix_putn_schar_schar(&putbuf.pos, sizeof(ncmagic1), ncmagic1);
    }
    if (status != NC_NOERR) return status;

    /* copy numrecs */
    nrecs = ncp->numrecs; 
    status = ncmpix_put_size_t(&putbuf.pos, nrecs, putbuf.version == 5 ? 8 : 4);
    if (status != NC_NOERR) return status;

    // assert((char *)putbuf.pos < (char *)putbuf.base + putbuf.size);

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

/* End Of put NC */

/* Begin Of get NC */

/*
 * Fetch the next header chunk.  the chunk is 'gbp->size' bytes big
 * Takes care to not overwrite leftover (unused) data in the buffer before
 * fetching a new chunk: the current aproach is to re-read the extra data.  
 *
 * NOTE: An alternate approach (which we do not do) would be to save the old
 *       data, read the next chunk and then copy the old data into the new
 *       chunk.  This alternate aproach might help if it is important for reads
 *       to be aligned.
 */
static int
hdr_fetch(bufferinfo *gbp) {
  int rank;
  MPI_Comm comm;
  int mpireturn;
  MPI_Offset slack;        /* any leftover data in the buffer */

  assert(gbp->base != NULL);
  
  comm = gbp->nciop->comm;
  MPI_Comm_rank(comm, &rank);

  /* XXX: this pointer math might not be good on 64 bit platforms */
  slack = gbp->size - ((char *)gbp->pos - (char *)gbp->base);
  /* . if gbp->pos and gbp->base are the same, there is no leftover buffer data
   *   to worry about.  
   * In the other extreme, where gbp->size == (gbp->pos - gbp->base), then all
   * data in the buffer has been consumed */

  if (slack == gbp->size) slack = 0;

  memset(gbp->base, 0, gbp->size);
  gbp->pos = gbp->base;
  gbp->index = 0;

  /* fileview is already entire file visible and MPI_File_read_at does not 
     change the file pointer */

  if (rank == 0) {
      MPI_Status mpistatus;
      mpireturn = MPI_File_read_at(gbp->nciop->collective_fh, (gbp->offset)-slack, gbp->base, 
                                   gbp->size, MPI_BYTE, &mpistatus);  
      if (mpireturn != MPI_SUCCESS) {
          ncmpii_handle_error(rank, mpireturn, "MPI_File_read_at");
          MPI_Finalize();
          return NC_EREAD;
      }
  }
  /* we might have had to backtrack */
  gbp->offset += (gbp->size - slack); 

  MPI_Bcast(gbp->base, gbp->size, MPI_BYTE, 0, comm);

  return NC_NOERR;
}

/*
 * Ensure that 'nextread' bytes are available.
 */
static int
hdr_check_buffer(bufferinfo *gbp, MPI_Offset nextread) {
  if ((char *)gbp->pos + nextread <= (char *)gbp->base + gbp->size)
    return NC_NOERR;
  return hdr_fetch(gbp);
}

static int
hdr_get_NCtype(bufferinfo *gbp, NCtype *typep) {
  int type = 0;
//  int sizeof_t = (gbp->version == 5) ? 8 : 4;
  int sizeof_t = 4;
  int status = hdr_check_buffer(gbp, sizeof_t);
  if (status != NC_NOERR)
    return status;

  status =  ncmpix_get_int_int(gbp->pos, &type);
  gbp->pos = (void *)((char *)gbp->pos + X_SIZEOF_INT);
  gbp->index += X_SIZEOF_INT;
  if (status != NC_NOERR)
    return status;
  *typep = (NCtype) type;
  return NC_NOERR;
}

static int
hdr_get_size_t(bufferinfo *gbp, MPI_Offset *sp) {
  /* in CDF-1 format, all integers are 32-bit
   * in CDF-2 format, only the variable begin (starting file offset) is 64-bit
   * in CDF-5 format, both variable's begin and size are 64-bit
   */
  MPI_Offset sizeof_t = (gbp->version == 5) ? 8 : 4;
  int status = hdr_check_buffer(gbp, sizeof_t);
  if (status != NC_NOERR)
    return status; 
  gbp->index += sizeof_t;
  return ncmpix_get_size_t((const void **)(&gbp->pos), sp, sizeof_t);
}

static int
hdr_get_NC_string(bufferinfo *gbp, NC_string **ncstrpp) {
  int status;
  MPI_Offset  nchars = 0, nbytes, padding, bufremain, strcount; 
  NC_string *ncstrp;
  char *cpos;
  char pad[X_ALIGN-1];

  status = hdr_get_size_t(gbp, &nchars);
  if (status != NC_NOERR)
    return status;

  ncstrp = ncmpii_new_NC_string(nchars, NULL);

  if (ncstrp == NULL)
    return NC_ENOMEM;

  nbytes = nchars * X_SIZEOF_CHAR;
  padding = _RNDUP(X_SIZEOF_CHAR * ncstrp->nchars, X_ALIGN)
            - X_SIZEOF_CHAR * ncstrp->nchars;
  bufremain = gbp->size - (MPI_Offset)((char *)gbp->pos - (char *)gbp->base);
  cpos = ncstrp->cp;

  while (nbytes > 0) {
    if (bufremain > 0) {
      strcount = MIN(bufremain, nbytes); 
      (void) memcpy(cpos, gbp->pos, strcount);
      nbytes -= strcount;
      gbp->pos = (void *)((char *)gbp->pos + strcount);
      gbp->index += strcount;
      cpos += strcount; 
      bufremain -= strcount;
    } else {
      status = hdr_fetch(gbp);
      if(status != NC_NOERR) {
        ncmpii_free_NC_string(ncstrp);
        return status;
      } 
      bufremain = gbp->size;
    }
  }

  if (padding > 0) {
    memset(pad, 0, X_ALIGN-1);
    if (memcmp(gbp->pos, pad, padding) != 0) {
      ncmpii_free_NC_string(ncstrp);
      return EINVAL;
    }
    gbp->pos = (void *)((char *)gbp->pos + padding);
    gbp->index += padding;
  }
  
  *ncstrpp = ncstrp;
  
  return NC_NOERR;  
}

static int
hdr_get_NC_dim(bufferinfo *gbp, NC_dim **dimpp) {
  int status;
  NC_string *ncstrp;
  NC_dim *dimp;

  status = hdr_get_NC_string(gbp, &ncstrp);
  if (status != NC_NOERR)
    return status;

  dimp = ncmpii_new_x_NC_dim(ncstrp);
  if(dimp == NULL)
    return NC_ENOMEM;

  status = hdr_get_size_t(gbp, &dimp->size);
  if(status != NC_NOERR) {
    ncmpii_free_NC_dim(dimp); /* frees name */
    return status;
  }

  *dimpp = dimp;
  return NC_NOERR;
}

static int
hdr_get_NC_dimarray(bufferinfo *gbp, NC_dimarray *ncap) {
  int status;
  NCtype type = NC_UNSPECIFIED; 
  NC_dim **dpp, **end;
  MPI_Offset tmp;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL);

  status = hdr_get_NCtype(gbp, &type);
  if(status != NC_NOERR)
    return status; 

  status = hdr_get_size_t(gbp, &tmp);
  if (status != NC_NOERR) return status;
  ncap->ndefined = tmp;

  if(ncap->ndefined == 0) {
    if (type != NC_DIMENSION && type != NC_UNSPECIFIED)
      return EINVAL;
  } else {
    if(type != NC_DIMENSION)
      return EINVAL;

    ncap->value = (NC_dim **) NCI_Malloc(ncap->ndefined * sizeof(NC_dim *));
    if(ncap->value == NULL)
      return NC_ENOMEM;
    ncap->nalloc = ncap->ndefined;

    dpp = ncap->value;
    end = &dpp[ncap->ndefined];

    for( /*NADA*/; dpp < end; dpp++) {
      status = hdr_get_NC_dim(gbp, dpp);
      if (status != NC_NOERR) {
        ncap->ndefined = dpp - ncap->value;
        ncmpii_free_NC_dimarrayV(ncap);
        return status;
      }
    }
  }

  return NC_NOERR;
}

static int
hdr_get_nc_type(bufferinfo *gbp, nc_type *typep) {
  int type = 0;
  int status = hdr_check_buffer(gbp, X_SIZEOF_INT);
  if(status != NC_NOERR)
    return status;

  status =  ncmpix_get_int_int(gbp->pos, &type);
  gbp->pos = (void *)((char *)gbp->pos + X_SIZEOF_INT); 
  gbp->index += X_SIZEOF_INT;
  if(status != NC_NOERR)
    return status;

  if (   type != NC_BYTE
      && type != NC_CHAR
      && type != NC_SHORT
      && type != NC_INT
      && type != NC_FLOAT
      && type != NC_DOUBLE) 
    return EINVAL; 
 
  *typep = (nc_type) type;

  return NC_NOERR;
}

MPI_Offset 
ncmpix_len_nctype(nc_type type) {
  switch(type) {
    case NC_BYTE:
    case NC_CHAR:
        return X_SIZEOF_CHAR;
    case NC_SHORT:
        return X_SIZEOF_SHORT;
    case NC_INT:
        return X_SIZEOF_INT;
    case NC_FLOAT:
        return X_SIZEOF_FLOAT;
    case NC_DOUBLE:
        return X_SIZEOF_DOUBLE;
    default: 
    assert("ncmpix_len_nctype bad type" == 0);
  }
  return 0;                
}

/*
 * Get the values of an attribute  
 */
static int
hdr_get_NC_attrV(bufferinfo *gbp, NC_attr *attrp) {
  int status;
  void *value = attrp->xvalue;
  char pad[X_ALIGN-1]; 
  MPI_Offset nbytes, esz, padding, bufremain, attcount;

  esz = ncmpix_len_nctype(attrp->type);
  padding = attrp->xsz - esz * attrp->nelems;
  bufremain = gbp->size - (MPI_Offset)((char *)gbp->pos - (char *)gbp->base);
  nbytes = esz * attrp->nelems;

  while (nbytes > 0) {
    if (bufremain > 0) {
      attcount = MIN(bufremain, nbytes);
      (void) memcpy(value, gbp->pos, attcount);
      nbytes -= attcount;
      gbp->pos = (void *)((char *)gbp->pos + attcount);
      gbp->index += attcount;
      value = (void *)((char *)value + attcount);
      bufremain -= attcount;
    } else {
      status = hdr_fetch(gbp);
      if(status != NC_NOERR) 
        return status;
      bufremain = gbp->size;
    }
  }
 
  if (padding > 0) {
    memset(pad, 0, X_ALIGN-1);
    if (memcmp(gbp->pos, pad, padding) != 0) 
      return EINVAL;
    gbp->pos = (void *)((char *)gbp->pos + padding);
    gbp->index += padding;
  }

  return NC_NOERR;
}

static int
hdr_get_NC_attr(bufferinfo *gbp, NC_attr **attrpp) {
  NC_string *strp;
  int status;
  nc_type type; 
  MPI_Offset nelems=0;
  NC_attr *attrp;

  status = hdr_get_NC_string(gbp, &strp);
  if(status != NC_NOERR)
    return status;

  status = hdr_get_nc_type(gbp, &type);
  if(status != NC_NOERR) {
    ncmpii_free_NC_string(strp);
    return status;
  }

  status = hdr_get_size_t(gbp, &nelems); 
  if(status != NC_NOERR) {
    ncmpii_free_NC_string(strp);
    return status;
  }

  attrp = ncmpii_new_x_NC_attr(strp, type, nelems);
  if(attrp == NULL) {
    ncmpii_free_NC_string(strp);
    return status;
  }

  status = hdr_get_NC_attrV(gbp, attrp);
  if(status != NC_NOERR) {
    ncmpii_free_NC_attr(attrp); /* frees strp */ 
    return status;
  }

  *attrpp = attrp; 
  
  return NC_NOERR; 
}

static int
hdr_get_NC_attrarray(bufferinfo *gbp, NC_attrarray *ncap){
  int status;
  NCtype type = NC_UNSPECIFIED;
  NC_attr **app, **end;
  MPI_Offset tmp;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL);

  status = hdr_get_NCtype(gbp, &type);
  if(status != NC_NOERR)
    return status; 

  status = hdr_get_size_t(gbp, &tmp);
  if (status != NC_NOERR) return status;
  ncap->ndefined = tmp;

  if(ncap->ndefined == 0) {
    if (type != NC_ATTRIBUTE && type != NC_UNSPECIFIED)
      return EINVAL;
  } else {
    if(type != NC_ATTRIBUTE)
      return EINVAL;

    ncap->value = (NC_attr **) NCI_Malloc(ncap->ndefined * sizeof(NC_attr *));
    if(ncap->value == NULL)
      return NC_ENOMEM;
    ncap->nalloc = ncap->ndefined; 

    app = ncap->value;
    end = &app[ncap->ndefined];
    for( /*NADA*/; app < end; app++) {
      status = hdr_get_NC_attr(gbp, app);
      if (status != NC_NOERR) {
        ncap->ndefined = app - ncap->value;
        ncmpii_free_NC_attrarrayV(ncap);
        return status;
      }
    }
  }
  
  return NC_NOERR;
}

static int
hdr_get_NC_var(bufferinfo *gbp, NC_var **varpp) {
  NC_string *strp;
  int status;
  MPI_Offset ndims=0, dim;
  MPI_Offset tmp_dimids=0;
  NC_var *varp;
  int i;

  status = hdr_get_NC_string(gbp, &strp);
  if(status != NC_NOERR)
    return status;
   
  status = hdr_get_size_t(gbp, &ndims);
  if(status != NC_NOERR) {
     ncmpii_free_NC_string(strp); 
     return status;
  }
 
  varp = ncmpii_new_x_NC_var(strp, ndims);
  if(varp == NULL) {
    ncmpii_free_NC_string(strp);
    return NC_ENOMEM;
  }

  for (dim = 0; dim < ndims; dim++ ) {
    status = hdr_check_buffer(gbp, (gbp->version == 5 ? 8 : 4));
    if(status != NC_NOERR) {
      ncmpii_free_NC_var(varp);
      return status;
    }
    status = hdr_get_size_t(gbp, &tmp_dimids);
    varp->dimids[dim] = (int)tmp_dimids;
    if(status != NC_NOERR) {
     return status;
    }

/*  
    if (gbp->version == 5) {
    status = ncmpix_getn_long_long((const void **)(&gbp->pos), 
                              1, (MPI_Offset*)varp->dimids + dim);
    } else {
    status = ncmpix_getn_int_int((const void **)(&gbp->pos), 
                              1, (int*) varp->dimids + dim);
    }

    if(status != NC_NOERR) {
      ncmpii_free_NC_var(varp);
      return status;
    }
*/
  }

  for (i=0; i< varp->ndims; i++){
  }

  status = hdr_get_NC_attrarray(gbp, &varp->attrs);
  if(status != NC_NOERR) {
    ncmpii_free_NC_var(varp);
    return status;
  }

  status = hdr_get_nc_type(gbp, &varp->type);
  if(status != NC_NOERR) {
    ncmpii_free_NC_var(varp);
    return status;
  } 

  status = hdr_get_size_t(gbp, &varp->len);
  if(status != NC_NOERR) {
    ncmpii_free_NC_var(varp);
    return status;
  }
  
  status = hdr_check_buffer(gbp, (gbp->version == 1 ? 4 : 8));
  if(status != NC_NOERR) {
    ncmpii_free_NC_var(varp);
    return status;
  }
  status = ncmpix_get_size_t((const void **)&gbp->pos,
                        &varp->begin, (gbp->version == 1 ? 4 : 8));
  if(status != NC_NOERR) {
    ncmpii_free_NC_var(varp);
    return status;
  }

  *varpp = varp;
  return NC_NOERR;
}

static int
hdr_get_NC_vararray(bufferinfo *gbp, NC_vararray *ncap) {
  int status;
  NCtype type = NC_UNSPECIFIED;
  NC_var **vpp, **end;
  MPI_Offset tmp;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL); 

  status = hdr_get_NCtype(gbp, &type);
  if(status != NC_NOERR)
    return status;
 
  status = hdr_get_size_t(gbp, &tmp);
  if(status != NC_NOERR)
    return status;
  ncap->ndefined = tmp; /* number of defined variables allowed < 2^32 */
 
  if(ncap->ndefined == 0) {
    if (type != NC_VARIABLE && type != NC_UNSPECIFIED)
      return EINVAL;
  } else {
    if(type != NC_VARIABLE)
      return EINVAL;
    ncap->value = (NC_var **) NCI_Malloc(ncap->ndefined * sizeof(NC_var *));
    if(ncap->value == NULL)
      return NC_ENOMEM; 
    ncap->nalloc = ncap->ndefined;
    vpp = ncap->value;
    end = &vpp[ncap->ndefined];


    for( /*NADA*/; vpp < end; vpp++) {
      status = hdr_get_NC_var(gbp, vpp);
      if (status != NC_NOERR) {
        ncap->ndefined = vpp - ncap->value;
        ncmpii_free_NC_vararrayV(ncap);
        return status;
      }
    }
  }

  return NC_NOERR;
}


int
ncmpii_hdr_get_NC(NC *ncp) {
  int status;
  bufferinfo getbuf;
  schar magic[sizeof(ncmagic)];
  MPI_Offset nrecs = 0;

  assert(ncp != NULL);

  /* Initialize the get buffer */

  getbuf.nciop = ncp->nciop;
  getbuf.offset = 0;     /* read from start of the file */

  /* CDF-5's minimum header size is 4 bytes more than CDF-1 and CDF-2's */
  getbuf.size = _RNDUP( MAX(MIN_NC_XSZ+4, ncp->chunk), X_ALIGN );
  if (getbuf.size > NC_DEFAULT_CHUNKSIZE)
    getbuf.size = NC_DEFAULT_CHUNKSIZE;

  getbuf.pos = getbuf.base = (void *)NCI_Malloc(getbuf.size);
  getbuf.index = 0;

  status = hdr_fetch(&getbuf);
  
  /* Get the header from get buffer */

  memset(magic, 0, sizeof(magic));
  status = ncmpix_getn_schar_schar(
          (const void **)(&getbuf.pos), sizeof(magic), magic);
  getbuf.index += sizeof(magic);
  /* don't need to worry about CDF-1 or CDF-2 
   *     if the first bits are not 'CDF-'  */
  if(memcmp(magic, ncmagic, sizeof(ncmagic)-1) != 0) {
    NCI_Free(getbuf.base);
    return NC_ENOTNC;
  }
  /* check version number in last byte of magic */
  if (magic[sizeof(ncmagic)-1] == 0x1) {
      getbuf.version = 1;
      fSet(ncp->flags, NC_32BIT);
  } else if (magic[sizeof(ncmagic)-1] == 0x2) {
      getbuf.version = 2;
      fSet(ncp->flags, NC_64BIT_OFFSET);
      if (sizeof(MPI_Offset) != 8) {
          /* take the easy way out: if we can't support all CDF-2
           * files, return immediately */
          NCI_Free(getbuf.base);
          return NC_ESMALL;
      }
  } else if (magic[sizeof(ncmagic)-1] == 0x5) {
      getbuf.version = 5;
      fSet(ncp->flags, NC_64BIT_DATA);
      if (sizeof(MPI_Offset) != 8) {
          NCI_Free(getbuf.base);
          return NC_ESMALL;
      }
  } else {
      NCI_Free(getbuf.base);
      return NC_ENOTNC;
  }
  

  status = hdr_check_buffer(&getbuf, (getbuf.version == 1) ? 4 : 8);
  if(status != NC_NOERR) {
    NCI_Free(getbuf.base);
    return status;
  } 
  status = ncmpix_get_size_t((const void **)(&getbuf.pos), &nrecs, (getbuf.version == 5) ? 8 : 4);

  if (getbuf.version == 5) {
    getbuf.index += X_SIZEOF_LONG;
  } else {
    getbuf.index += X_SIZEOF_SIZE_T;
  }
  if(status != NC_NOERR) {
    NCI_Free(getbuf.base);
    return status;
  }
  ncp->numrecs = nrecs;

  assert((char *)getbuf.pos < (char *)getbuf.base + getbuf.size);

  status = hdr_get_NC_dimarray(&getbuf, &ncp->dims);
  if(status != NC_NOERR) {
    NCI_Free(getbuf.base);
    return status;
  }
    
 

  status = hdr_get_NC_attrarray(&getbuf, &ncp->attrs); 
  if(status != NC_NOERR) {
    NCI_Free(getbuf.base);
    return status;
  }


  status = hdr_get_NC_vararray(&getbuf, &ncp->vars);
  if(status != NC_NOERR) {
    NCI_Free(getbuf.base);
    return status; 
  }
 
  
  ncp->xsz = ncmpii_hdr_len_NC(ncp);
  status = ncmpii_NC_computeshapes(ncp);
  NCI_Free(getbuf.base);
  
  
  return status;
}

/* End Of get NC */

#define METADATA_CONSISTENCY_CHECK
/* TODO: this should be set at the configure time by a user option */

#define WARN_STR "Warning (inconsistent metadata):"

static int
ncmpii_comp_dims(NC_dimarray *nc_dim1,
                 NC_dimarray *nc_dim2)
{
    int i;
    if (nc_dim1->ndefined != nc_dim2->ndefined) {
        fprintf(stderr,"Error: number of dimensions defined is inconsistent %d != %d\n",
                nc_dim1->ndefined, nc_dim2->ndefined);
        return NC_EDIMS_NELEMS_MULTIDEFINE;
    }

    for (i=0; i<nc_dim1->ndefined; i++) {
#ifdef METADATA_CONSISTENCY_CHECK
        NC_string *name1, *name2;
        name1 = nc_dim1->value[i]->name;
        name2 = nc_dim2->value[i]->name;

        if (name1->nchars != name2->nchars ||
            strncmp(name1->cp, name2->cp, name1->nchars) != 0)
            printf("%s dimension name %s != %s\n", WARN_STR,
                   name1->cp,name2->cp);
#endif
        if (nc_dim1->value[i]->size != nc_dim2->value[i]->size) {
            /* inconsistency in dimension size is fatal */
            fprintf(stderr,"Error: dimension %s's size inconsistent %lld != %lld\n",
                    nc_dim1->value[i]->name->cp,nc_dim1->value[i]->size,
                    nc_dim2->value[i]->size);
            return NC_EDIMS_SIZE_MULTIDEFINE;
        }
    }
    return NC_NOERR;
}

static int
ncmpii_comp_attrs(NC_attrarray *nc_attr1,
                  NC_attrarray *nc_attr2)
{
    int        i, j, num;
    int       *ia, *ib;
    short int *sa, *sb;
    float     *fa, *fb;
    double    *da, *db;

    if (nc_attr1->ndefined != nc_attr2->ndefined) {
        printf("%s number of attributes (%d != %d)\n", WARN_STR,
               nc_attr1->ndefined, nc_attr2->ndefined);
        return NC_NOERR;
        /* no need to compare further */
    }

    for (i=0; i<nc_attr1->ndefined; i++) {
        NC_attr *v1 = nc_attr1->value[i];
        NC_attr *v2 = nc_attr2->value[i];

        if (v1->name->nchars != v2->name->nchars ||
            strncmp(v1->name->cp, v2->name->cp, v1->name->nchars) != 0)
            printf("%s attribute name (%s != %s)\n", WARN_STR,
                   v1->name->cp,v2->name->cp);

        if (v1->xsz != v2->xsz)
            printf("%s attribute \"%s\" size (%lld != %lld)\n", WARN_STR,
                   v1->name->cp,lld(v1->xsz),lld(v2->xsz));

        if (v1->type != v2->type)
            printf("%s attribute \"%s\" type (%d != %d)\n", WARN_STR,
                   v1->name->cp,v1->type,v2->type);

        if (v1->nelems != v2->nelems)
            printf("%s attribute \"%s\" length (%lld != %lld)\n", WARN_STR,
                   v1->name->cp,lld(v1->nelems),lld(v2->nelems));

        num = MIN(v1->nelems, v2->nelems);
        switch (v1->type) {
            case NC_CHAR:
            case NC_BYTE:
                if (strncmp(v1->xvalue, v2->xvalue, num))
                    printf("%s attribute \"%s\" BYTE/CHAR (%s != %s)\n", WARN_STR,
                           v1->name->cp,(char*)v1->xvalue,(char*)v2->xvalue);
                break;
            case NC_SHORT:
                sa = v1->xvalue;
                sb = v2->xvalue;
                for (j=0; j<num; j++) {
                    if (sa[j] != sb[j]) {
                        printf("%s attribute \"%s\" SHORT (%d != %d)\n", WARN_STR,
                               v1->name->cp,sa[j],sb[j]);
                        break;
                    }
                }
                break;
            case NC_INT:
                ia = v1->xvalue;
                ib = v2->xvalue;
                for (j=0; j<num; j++) {
                    if (ia[j] != ib[j]) {
                        printf("%s attribute \"%s\" INT (%d != %d)\n", WARN_STR,
                               v1->name->cp,ia[j],ib[j]);
                        break;
                    }
                }
                break;
            case NC_FLOAT:
                fa = v1->xvalue;
                fb = v2->xvalue;
                for (j=0; j<num; j++) {
			/* floating-point inequality here but we genuinely do
			 * expect all processors to set bit-for-bit identical
			 * headers */
                    if (fa[j] != fb[j]) {
                        printf("%s attribute \"%s\" FLOAT (%f != %f)\n", WARN_STR,
                               v1->name->cp,fa[j],fb[j]);
                        break;
                    }
                }
                break;
            case NC_DOUBLE:
                da = v1->xvalue;
                db = v2->xvalue;
                for (j=0; j<num; j++) {
			/* floating-point inequality here but we genuinely do
			 * expect all processors to set bit-for-bit identical
			 * headers */
                    if (da[j] != db[j]) {
                        printf("%s attribute \"%s\" DOUBLE (%f != %f)\n", WARN_STR,
                               v1->name->cp,da[j],db[j]);
                        break;
                    }
                }
                break;
            default: break;
        }
    }
    return NC_NOERR;
}

static int
ncmpii_comp_vars(NC_vararray *nc_var1,
                 NC_vararray *nc_var2)
{
    int i, j;
    if (nc_var1->ndefined != nc_var2->ndefined) {
        fprintf(stderr,"Error: number of defined variables is inconsistent %d != %d\n",
                nc_var1->ndefined, nc_var2->ndefined);
        return NC_EVARS_NELEMS_MULTIDEFINE;
    }

    for (i=0; i<nc_var1->ndefined; i++) {
        NC_var *v1 = nc_var1->value[i];
        NC_var *v2 = nc_var2->value[i];

#ifdef METADATA_CONSISTENCY_CHECK
        if (v1->name->nchars != v2->name->nchars ||
            strncmp(v1->name->cp, v2->name->cp, v1->name->nchars) != 0)
            printf("%s variable name %s != %s\n", WARN_STR,
                   v1->name->cp,v2->name->cp);
#endif
        if (v1->ndims != v2->ndims) {
            fprintf(stderr,"Error: variable %s's ndims is inconsistent %d != %d\n",
                    v1->name->cp, (int)v1->ndims, (int)v2->ndims);
            return NC_EVARS_NDIMS_MULTIDEFINE;
        }

        for (j=0; j<v1->ndims; j++) {
            if (v1->dimids[j] != v2->dimids[j]) {
                fprintf(stderr,"Error: variable %s's %dth dim ID is inconsistent %d != %d\n",
                        v1->name->cp, j, v1->dimids[j], v2->dimids[j]);
                return NC_EVARS_DIMIDS_MULTIDEFINE;
            }
        }

        if (v1->type != v2->type) {
            fprintf(stderr,"Error: variable %s's type is inconsistent %d != %d\n",
                    v1->name->cp, v1->type, v2->type);
            return NC_EVARS_TYPE_MULTIDEFINE;
        }

        if (v1->len != v2->len) {
            fprintf(stderr,"Error: variable %s's len is inconsistent %lld != %lld\n",
                    v1->name->cp, v1->len, v2->len);
            return NC_EVARS_LEN_MULTIDEFINE;
        }

        if (v1->begin != v2->begin) {
            fprintf(stderr,"Error: variable %s's begin is inconsistent %lld != %lld\n",
                    v1->name->cp, v1->begin, v2->begin);
            return NC_EVARS_BEGIN_MULTIDEFINE;
        }

#ifdef METADATA_CONSISTENCY_CHECK
        /* compare variable's attributes */
        ncmpii_comp_attrs(&(v1->attrs), &(v2->attrs));
#endif
    }
    return NC_NOERR;
}

int
ncmpii_hdr_check_NC(bufferinfo *getbuf, /* header from root */
                    NC         *ncp) {
    int rank, status;
    schar magic[sizeof(ncmagic)];
    MPI_Offset nrecs=0, chunksize=NC_DEFAULT_CHUNKSIZE;
    NC *temp_ncp;
    char *root_magic;

    assert(ncp != NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* check header's magic */
    memset(magic, 0, sizeof(magic));
    status = ncmpix_getn_schar_schar(
            (const void **)&getbuf->pos, sizeof(magic), magic);
    getbuf->index += sizeof(magic);
    /* don't need to worry about CDF-1 or CDF-2
     * if the first bits are not 'CDF-'  */
    if (memcmp(magic, ncmagic, sizeof(ncmagic)-1) != 0)
        return NC_ENOTNC;

    temp_ncp = ncmpii_new_NC(&chunksize);

    /* consistency of magic numbers should have already been checked during
     * ncmpi_create()
     */
#ifdef CHECK_MAGIC_CONSISTENCY
    if (magic[sizeof(ncmagic)-1] == 0x5)
        root_magic = "CDF-5";
    else if (magic[sizeof(ncmagic)-1] == 0x2)
        root_magic = "CDF-2";
    else if (magic[sizeof(ncmagic)-1] == 0x1)
        root_magic = "CDF-1";
    else
        return NC_ENOTNC;

    /* check version number in last byte of magic */
    if (ncp->flags & NC_64BIT_DATA) {
        if (magic[sizeof(ncmagic)-1] != 0x5) {
            fprintf(stderr,"Error: file format inconsistent (local=%s, root=%s)\n", "CDF-5", root_magic);
            return NC_ECMODE;
        }
        getbuf->version = 5;
        fSet(temp_ncp->flags, NC_64BIT_DATA);

        /* shouldn't this check have already been done? */
        if (sizeof(MPI_Offset) != 8)
            return NC_ESMALL;
    }
    else if (ncp->flags & NC_64BIT_OFFSET) {
        if (magic[sizeof(ncmagic)-1] != 0x2) {
            fprintf(stderr,"Error: file format inconsistent (local=%s, root=%s)\n", "CDF-2", root_magic);
            return NC_ECMODE;
        }
        getbuf->version = 2;
        fSet(temp_ncp->flags, NC_64BIT_OFFSET);

        /* shouldn't this check have already been done? */
        if (sizeof(MPI_Offset) != 8)
            /* take the easy way out: if we can't support all CDF-2
             * files, return immediately */
            return NC_ESMALL;
    }
    else {
        if (magic[sizeof(ncmagic)-1] != 0x1) {
            fprintf(stderr,"Error: file format inconsistent (local=%s, root=%s)\n", "CDF-1", root_magic);
            return NC_ECMODE;
        }
        getbuf->version = 1;
    }
#else
    if (magic[sizeof(ncmagic)-1] == 0x5) {
        getbuf->version = 5;
        fSet(temp_ncp->flags, NC_64BIT_DATA);
        if (sizeof(MPI_Offset) != 8)
            return NC_ESMALL;
    }
    else if (magic[sizeof(ncmagic)-1] == 0x2) {
        getbuf->version = 2;
        fSet(temp_ncp->flags, NC_64BIT_OFFSET);
        if (sizeof(MPI_Offset) != 8)
            return NC_ESMALL;
    }
    else if (magic[sizeof(ncmagic)-1] == 0x1)
        getbuf->version = 1;
    else
        return NC_ENOTNC;
#endif

    status = hdr_check_buffer(getbuf, (getbuf->version == 1) ? 4 : 8);
    if (status != NC_NOERR)
        return status;

    status = ncmpix_get_size_t((const void **)(&getbuf->pos), &nrecs,
                               (getbuf->version == 5) ? 8 : 4);

    if (getbuf->version == 5)
        getbuf->index += X_SIZEOF_LONG;
    else
        getbuf->index += X_SIZEOF_SIZE_T;

    if (status != NC_NOERR)
        return status;

    temp_ncp->numrecs = nrecs;

    if (temp_ncp->numrecs != ncp->numrecs) {
        fprintf(stderr,"Error: number of record variables is inconsistent %lld != %lld\n",
                temp_ncp->numrecs, ncp->numrecs);
        return NC_ENUMRECS_MULTIDEFINE;
    }

    assert((char *)getbuf->pos < (char *)getbuf->base + getbuf->size);

    status = hdr_get_NC_dimarray(getbuf, &temp_ncp->dims);
    if (status != NC_NOERR)
        return status;

    status = ncmpii_comp_dims(&temp_ncp->dims, &ncp->dims);
    if (status != NC_NOERR)
        return status;

    status = hdr_get_NC_attrarray(getbuf, &temp_ncp->attrs);
    if (status != NC_NOERR)
        return status;

#ifdef METADATA_CONSISTENCY_CHECK
    status = ncmpii_comp_attrs(&temp_ncp->attrs, &ncp->attrs);
    if (status != NC_NOERR)
        return status;
#endif

    status = hdr_get_NC_vararray(getbuf, &temp_ncp->vars);
    if (status != NC_NOERR)
        return status;

    status = ncmpii_comp_vars(&temp_ncp->vars, &ncp->vars);
    if (status != NC_NOERR)
        return status;

    temp_ncp->xsz = ncmpii_hdr_len_NC(temp_ncp);
    status = ncmpii_NC_computeshapes(temp_ncp);
  
    ncmpii_free_NC(temp_ncp);
    return status;
}

