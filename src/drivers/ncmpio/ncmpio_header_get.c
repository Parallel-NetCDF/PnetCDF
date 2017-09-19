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
#include <ncx.h>
#include "ncmpio_NC.h"

#define NC_MAGIC_LEN 4

/*----< compute_var_shape() >------------------------------------------------*/
/* Recompute the shapes of all variables: shape, xsz, and len
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

#if 0
/*----< hdr_len_NC_name() >--------------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_name(const NC_string *ncstrp,
                int              sizeof_NON_NEG)     /* NON_NEG */
{
    /* netCDF file format:
     * name       = nelems namestring
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
#endif

/*----< hdr_len_NC_dim() >---------------------------------------------------*/
inline static MPI_Offset
hdr_len_NC_dim(const NC_dim *dimp,
               int           sizeof_NON_NEG)     /* NON_NEG */
{
    /* netCDF file format:
     *  ...
     * dim        = name dim_length
     * dim_length = NON_NEG
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;

    assert(dimp != NULL);

    sz  = sizeof_NON_NEG + _RNDUP(dimp->name_len, X_ALIGN); /* name */
    sz += sizeof_NON_NEG;                                   /* dim_length */

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

    sz  = sizeof_NON_NEG + _RNDUP(attrp->name_len, X_ALIGN); /* name */
    sz += X_SIZEOF_NC_TYPE;                                  /* nc_type */
    sz += sizeof_NON_NEG;                                    /* nelems */
    sz += attrp->xsz;                                        /* [values ...] */

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
    sz  = sizeof_NON_NEG + _RNDUP(varp->name_len, X_ALIGN);   /* name */
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
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_read_at");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD)
        }
        else {
#ifdef _USE_MPI_GET_COUNT
            int get_size; /* actual read amount can be smaller */
            MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
            gbp->get_size += get_size;
#else
            gbp->get_size += gbp->size;
#endif
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

/*----< hdr_check_buffer() >-------------------------------------------------*/
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
/* NC_tag is 32-bit integer and can be the followings:
 *     ZERO (NC_UNSPECIFIED)
 *     NC_DIMENSION
 *     NC_ATTRIBUTE
 *     NC_VARIABLE
 */
inline static int
hdr_get_NC_tag(bufferinfo *gbp,
               NC_tag     *tagp)
{
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
                nc_type    *xtypep)
{
    /* nc_type is 4-byte integer, X_SIZEOF_INT */
    int status;
    uint xtype;

    status = hdr_check_buffer(gbp, X_SIZEOF_INT);
    if (status != NC_NOERR) return status;

    status = ncmpix_get_uint32((const void**)(&gbp->pos), &xtype);
    if (status != NC_NOERR) return status;

    if (xtype != NC_CHAR    &&
        xtype != NC_BYTE    &&
        xtype != NC_UBYTE   &&
        xtype != NC_SHORT   &&
        xtype != NC_USHORT  &&
        xtype != NC_INT     &&
        xtype != NC_UINT    &&
        xtype != NC_FLOAT   &&
        xtype != NC_DOUBLE  &&
        xtype != NC_INT64   &&
        xtype != NC_UINT64
       )
        DEBUG_RETURN_ERROR(NC_EBADTYPE)

    *xtypep = (nc_type) xtype;
    return NC_NOERR;
}

/*----< hdr_get_NC_name() >---------------------------------------------------*/
static int
hdr_get_NC_name(bufferinfo  *gbp,
                char       **namep)
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
    MPI_Aint pos_addr, base_addr;
    MPI_Offset nchars, padding, bufremain, strcount;

    *namep = NULL;
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

    /* Allocate a NC_string structure large enough to hold nchars characters.
     * Note nchars is strlen(namestring) without terminal character.
     */
    *namep = (char*)NCI_Malloc((size_t)nchars + 1);
    if (*namep == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    (*namep)[nchars] = '\0'; /* add terminal character */

    /* X_SIZEOF_CHAR is defined as 1 in classical CDF formats
    padding = _RNDUP(X_SIZEOF_CHAR * nchars, X_ALIGN) - X_SIZEOF_CHAR * nchars;
    */
    padding = _RNDUP(nchars, X_ALIGN) - nchars;
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    bufremain = gbp->size - (pos_addr - base_addr);
    cpos = *namep;

    /* get namestring with padding (the space in file allocated for string
     * namestring is upward aligned with 4 bytes */
    while (nchars > 0) {
        if (bufremain > 0) {
            strcount = MIN(bufremain, nchars);
            if (strcount != (size_t)strcount)
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            memcpy(cpos, gbp->pos, (size_t)strcount);
            nchars -= strcount;
            gbp->pos = (void *)((char *)gbp->pos + strcount);
            cpos += strcount;
            bufremain -= strcount;
        } else {
            err = hdr_fetch(gbp);
            if (err != NC_NOERR) {
                NCI_Free(*namep);
                *namep = NULL;
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
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header corrupted, non-zero padding found\n",__FILE__,__func__,__LINE__);
#endif
            NCI_Free(*namep);
            *namep = NULL;
            DEBUG_RETURN_ERROR(NC_ENOTNC)
        }
        gbp->pos = (void *)((char *)gbp->pos + padding);
    }

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
    char *name;
    NC_dim *dimp;

    *dimpp = NULL;

    /* get name */
    status = hdr_get_NC_name(gbp, &name);
    if (status != NC_NOERR) return status;

    /* allocate and initialize NC_dim object */
    dimp = (NC_dim*) NCI_Malloc(sizeof(NC_dim));
    if (dimp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    dimp->name     = name;
    dimp->name_len = strlen(name);

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
    if (status != NC_NOERR) { /* frees dim */
        NCI_Free(dimp->name);
        NCI_Free(dimp);
        return status;
    }

    *dimpp = dimp;
    return NC_NOERR;
}

/*----< hdr_get_NC_dimarray() >----------------------------------------------*/
/* For CDF-5 format, nelems (number of dimensions) is of type non-negative
 * INT64. However, argument ndims/dimid in all PnetCDF/NetCDF APIs are of type
 * int. Thus, we only need to use type int for internal metadata, ndefined. If
 * nelems in the input file is more than NC_MAX_DIMS, then it violates the
 * format specifications (NC_ENOTNC).
 */
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
    int i, status, ndefined=0;
    size_t alloc_size;
    NC_tag tag = NC_UNSPECIFIED;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* read NC_tag (NC_DIMENSION or ZERO) from gbp buffer */
    status = hdr_get_NC_tag(gbp, &tag);
    if (status != NC_NOERR) return status;

    /* read nelems (number of dimensions) from gbp buffer */
    if (gbp->version < 5) { /* nelems is <non-negative INT> */
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        if (status == NC_NOERR && tmp > NC_MAX_DIMS)
            /* cannot be more than max number of dimensions */
            status = NC_ENOTNC;
        else
            ndefined = (int)tmp;
    }
    else { /* nelems is <non-negative INT64> */
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        if (status == NC_NOERR && tmp > NC_MAX_DIMS)
            /* cannot be more than max number of dimensions */
            status = NC_ENOTNC;
        else
            ndefined = (int)tmp;
    }
    if (status != NC_NOERR) return status;

    /* Now ndefined is in between 0 and NC_MAX_DIMS */
    ncap->ndefined = ndefined;

    ncap->unlimited_id = -1;

    /* From the CDF file format specification, the tag is either NC_DIMENSION
     * or ABSENT (ZERO), but we follow NetCDF library to skip checking the tag
     * when ndefined is zero.
     */
    if (ndefined == 0) return NC_NOERR;

    /* Now, ndefined > 0, tag must be NC_DIMENSION */
    if (tag != NC_DIMENSION) {
#ifdef PNETCDF_DEBUG
        fprintf(stderr,"Error in file %s func %s line %d: NetCDF header corrupted, expecting tag NC_DIMENSION but got %d\n",__FILE__,__func__,__LINE__,tag);
#endif
        DEBUG_RETURN_ERROR(NC_ENOTNC)
    }

    alloc_size = _RNDUP(ncap->ndefined, NC_ARRAY_GROWBY);
    ncap->value = (NC_dim**) NCI_Malloc(alloc_size * sizeof(NC_dim*));
    if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    for (i=0; i<ndefined; i++) {
        status = hdr_get_NC_dim(gbp, ncap->value + i);
        if (status != NC_NOERR) { /* error: fail to get the next dim */
            ncap->ndefined = i; /* update to no. successful defined */
            ncmpio_free_NC_dimarray(ncap);
            return status;
        }
        if (ncap->value[i]->size == NC_UNLIMITED)
            ncap->unlimited_id = i; /* ID of unlimited dimension */
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
    int xsz;
    void *value = attrp->xvalue;
    MPI_Offset nbytes, padding, bufremain, attcount;
    MPI_Aint pos_addr, base_addr;

    ncmpii_xlen_nc_type(attrp->xtype, &xsz);
    nbytes = attrp->nelems * xsz;
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
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header corrupted, non-zero padding found\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_RETURN_ERROR(NC_ENOTNC)
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
    char *name;
    nc_type type;
    MPI_Offset nelems;
    NC_attr *attrp;

    /* get name */
    status = hdr_get_NC_name(gbp, &name);
    if (status != NC_NOERR) return status;

    /* get nc_type */
    status = hdr_get_nc_type(gbp, &type);
    if (status != NC_NOERR) return status;

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
    if (status != NC_NOERR) return status;

    /* allocate space for attribute object */
    status = ncmpio_new_NC_attr(name, type, nelems, &attrp);
    if (status != NC_NOERR) {
        NCI_Free(name);
        return status;
    }

    /* get [values ...] */
    status = hdr_get_NC_attrV(gbp, attrp);
    if (status != NC_NOERR) {
        if (attrp->xvalue != NULL) NCI_Free(attrp->xvalue);
        NCI_Free(attrp->name);
        NCI_Free(attrp);
        return status;
    }

    *attrpp = attrp;
    return NC_NOERR;
}

/*----< hdr_get_NC_attrarray() >---------------------------------------------*/
/* For CDF-5 format, nelems (number of attributes) is of type non-negative
 * INT64. However, argument nattrs in all PnetCDF/NetCDF APIs are of type int.
 * Thus, we only need to use type int for internal metadata, ndefined. If
 * nelems in the input file is more than NC_MAX_ATTRS, then it violates the
 * format specifications (NC_ENOTNC).
 */
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
    int i, status, ndefined=0;
    size_t alloc_size;
    NC_tag tag = NC_UNSPECIFIED;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* read NC_tag (NC_ATTRIBUTE or ZERO) from gbp buffer */
    status = hdr_get_NC_tag(gbp, &tag);
    if (status != NC_NOERR) return status;

    /* read nelems (number of attributes) from gbp buffer */
    if (gbp->version < 5) { /* nelems is <non-negative INT> */
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        if (status == NC_NOERR && tmp > NC_MAX_ATTRS)
            /* cannot be more than max number of attributes */
            status = NC_ENOTNC;
        else
            ndefined = (int)tmp;
    }
    else { /* nelems is <non-negative INT64> */
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        if (status == NC_NOERR && tmp > NC_MAX_ATTRS)
            /* cannot be more than max number of attributes */
            status = NC_ENOTNC;
        else
            ndefined = (int)tmp;
    }
    if (status != NC_NOERR) return status;

    /* Now ndefined is in between 0 and NC_MAX_ATTRS */
    ncap->ndefined = ndefined;

    /* From the CDF file format specification, the tag is either NC_ATTRIBUTE
     * or ABSENT (ZERO), but we follow NetCDF library to skip checking the tag
     * when ndefined is zero.
     */
    if (ndefined == 0) return NC_NOERR;

    /* Now, ndefined > 0, tag must be NC_ATTRIBUTE */
    if (tag != NC_ATTRIBUTE) {
#ifdef PNETCDF_DEBUG
        fprintf(stderr,"Error in file %s func %s line %d: NetCDF header corrupted, expecting tag NC_ATTRIBUTE but got %d\n",__FILE__,__func__,__LINE__,tag);
#endif
        DEBUG_RETURN_ERROR(NC_ENOTNC)
    }

    alloc_size = _RNDUP(ncap->ndefined, NC_ARRAY_GROWBY);
    ncap->value = (NC_attr**)NCI_Malloc(alloc_size * sizeof(NC_attr*));
    if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* get [attr ...] */
    for (i=0; i<ndefined; i++) {
        status = hdr_get_NC_attr(gbp, ncap->value + i);
        if (status != NC_NOERR) { /* Error: fail to get the next att */
            ncap->ndefined = i; /* update to no. successful defined */
            ncmpio_free_NC_attrarray(ncap);
            return status;
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
    char *name;
    MPI_Offset ndims=0, dim;
    MPI_Offset tmp_dimids=0;
    NC_var *varp;

    /* get name */
    status = hdr_get_NC_name(gbp, &name);
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
    if (status != NC_NOERR) return status;

    if (ndims != (int)ndims) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    /* allocate space for var object */
    varp = ncmpio_new_NC_var(name, (int)ndims);
    if (varp == NULL) {
        NCI_Free(name);
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
    status = hdr_get_nc_type(gbp, &varp->xtype);
    if (status != NC_NOERR) {
        ncmpio_free_NC_var(varp);
        return status;
    }
    ncmpii_xlen_nc_type(varp->xtype, &varp->xsz);

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
/* For CDF-5 format, nelems (number of variables) is of type non-negative
 * INT64. However, argument nvars/varid in all PnetCDF/NetCDF APIs are of type
 * int. Thus, we only need to use type int for internal metadata, ndefined. If
 * nelems in the input file is more than NC_MAX_VARS, then it violates the
 * format specifications (NC_ENOTNC).
 */
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
    int i, status, ndefined=0;
    size_t alloc_size;
    NC_tag tag = NC_UNSPECIFIED;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* read NC_tag (NC_VARIABLE or ZERO) from gbp buffer */
    status = hdr_get_NC_tag(gbp, &tag);
    if (status != NC_NOERR) return status;

    /* read nelems (number of variables) from gbp buffer */
    if (gbp->version < 5) { /* nelems is <non-negative INT> */
        uint tmp;
        status = hdr_get_uint32(gbp, &tmp);
        if (status == NC_NOERR && tmp > NC_MAX_VARS)
            /* cannot be more than max number of attributes */
            status = NC_ENOTNC;
        else
            ndefined = (int)tmp;
    }
    else { /* nelems is <non-negative INT64> */
        uint64 tmp;
        status = hdr_get_uint64(gbp, &tmp);
        if (status == NC_NOERR && tmp > NC_MAX_VARS)
            /* cannot be more than max number of attributes */
            status = NC_ENOTNC;
        else
            ndefined = (int)tmp;
    }
    if (status != NC_NOERR) return status;

    /* Now ndefined is in between 0 and NC_MAX_ATTRS */
    ncap->ndefined = ndefined;

    /* From the CDF file format specification, the tag is either NC_VARIABLE
     * or ABSENT (ZERO), but we follow NetCDF library to skip checking the tag
     * when ndefined is zero.
     */
    if (ndefined == 0) return NC_NOERR;

    /* Now, ndefined > 0, tag must be NC_VARIABLE */
    if (tag != NC_VARIABLE) {
#ifdef PNETCDF_DEBUG
        fprintf(stderr,"Error in file %s func %s line %d: NetCDF header corrupted, expecting tag NC_VARIABLE but got %d\n",__FILE__,__func__,__LINE__,tag);
#endif
        DEBUG_RETURN_ERROR(NC_ENOTNC)
    }

    alloc_size = _RNDUP(ncap->ndefined, NC_ARRAY_GROWBY);
    ncap->value = (NC_var**) NCI_Malloc(alloc_size * sizeof(NC_var*));
    if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* get [var ...] */
    for (i=0; i<ndefined; i++) {
        status = hdr_get_NC_var(gbp, ncap->value + i);
        if (status != NC_NOERR) { /* Error: fail to get the next var */
            ncap->ndefined = i; /* update to no. successful defined */
            ncmpio_free_NC_vararray(ncap);
            return status;
        }
        ncap->value[i]->varid = i;
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

    xlen  = NC_MAGIC_LEN;                                                    /* magic */
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
    char magic[NC_MAGIC_LEN];
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
    status = ncmpix_getn_text((const void**)(&getbuf.pos), NC_MAGIC_LEN, magic);
    if (status != NC_NOERR) return status;

    /* check if the first three bytes are 'C','D','F' */
    if (memcmp(magic, "CDF", 3) != 0) {
        /* check if is HDF5 file */
        char signature[8], *hdf5_signature="\211HDF\r\n\032\n";
        ncmpix_getn_text((const void**)(&getbuf.pos), 8, signature);
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
    if (magic[3] == 0x1) {
        getbuf.version = ncp->format = 1;
    } else if (magic[3] == 0x2) {
        getbuf.version = ncp->format = 2;
#if SIZEOF_MPI_OFFSET < 8
        /* take the easy way out: if we can't support all CDF-2
         * files, return immediately */
        NCI_Free(getbuf.base);
        DEBUG_RETURN_ERROR(NC_ESMALL)
#endif
    } else if (magic[3] == 0x5) {
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

    /* Recompute the shapes of all variables (shape, xsz, len)
     * Sets ncp->begin_var to start of first variable.
     * Sets ncp->begin_rec to start of first record variable.
     */
    status = compute_var_shape(ncp);
    if (status != NC_NOERR) goto fn_exit;

    /* Check whether variable sizes are legal for the given file format */
    status = ncmpio_NC_check_vlens(ncp);
    if (status != NC_NOERR) goto fn_exit;

    /* Check whether variable begins are in an increasing order.
     * Adding this check here is necessary for detecting corrupted metadata. */
    status = ncmpio_NC_check_voffs(ncp);
    if (status != NC_NOERR) goto fn_exit;

fn_exit:
    ncp->get_size += getbuf.get_size;
    NCI_Free(getbuf.base);

    return status;
}

