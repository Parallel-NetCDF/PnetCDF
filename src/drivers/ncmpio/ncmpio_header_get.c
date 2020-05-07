/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <string.h>  /* memcpy(), memcmp(), memset(), memmove() */
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
static MPI_Offset
hdr_len_NC_name(const NC_string *ncstrp, int sizeof_NON_NEG)
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
static MPI_Offset
hdr_len_NC_dim(const NC_dim *dimp, int sizeof_NON_NEG)
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
static MPI_Offset
hdr_len_NC_dimarray(const NC_dimarray *ncap, int sizeof_NON_NEG)
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
static MPI_Offset
hdr_len_NC_attr(const NC_attr *attrp, int sizeof_NON_NEG)
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
static MPI_Offset
hdr_len_NC_attrarray(const NC_attrarray *ncap, int sizeof_NON_NEG)
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
static MPI_Offset
hdr_len_NC_var(const NC_var *varp,
               int           sizeof_off_t,    /* OFFSET */
               int           sizeof_NON_NEG)  /* NON_NEG */
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
static MPI_Offset
hdr_len_NC_vararray(const NC_vararray *ncap,
                    int                sizeof_NON_NEG, /* NON_NEG */
                    int                sizeof_off_t)   /* OFFSET */
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
/* Fetch the next header chunk. The chunk buffer, pointed by gbp->base, is of
 * size 'gbp->size' bytes. Be careful not to overwrite leftover (yet to be
 * used) data in the buffer before fetching a new chunk.
 */
static int
hdr_fetch(bufferinfo *gbp) {
    int rank, err=NC_NOERR, mpireturn;

    assert(gbp->base != NULL);

    MPI_Comm_rank(gbp->comm, &rank);
    if (rank == 0) {
        char *readBuf;
        size_t slack, readLen;
        MPI_Status mpistatus;

        /* any leftover data in the buffer */
        slack = gbp->size - (gbp->pos - gbp->base);
        if (slack == gbp->size) slack = 0;

        /* When gbp->size == (gbp->pos - gbp->base), all data in the buffer has
         * been consumed. If not, then read additional header of size
         * (gbp->size - slack) into a contiguous buffer, pointed by gbp->base +
         * slack.
         */

        readBuf = gbp->base;
        readLen = gbp->size;
        if (slack > 0) { /* move slack to beginning of the buffer, gbp->base */
            memmove(gbp->base, gbp->pos, slack);
            readBuf += slack;
            readLen -= slack;
        }

        /* explicitly initialize mpistatus object to 0. For zero-length read,
         * MPI_Get_count may report incorrect result for some MPICH version,
         * due to the uninitialized MPI_Status object passed to MPI-IO calls.
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));

        /* fileview is already entire file visible and MPI_File_read_at does
           not change the file pointer */
        TRACE_IO(MPI_File_read_at)(gbp->collective_fh, gbp->offset, readBuf,
                                   readLen, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_read_at");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD)
        }
        else {
            /* Obtain the actual read amount. It may be smaller than readLen,
             * when the remaining file size is smaller than read chunk size.
             */
            int get_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
            gbp->get_size += get_size;

            /* If actual read amount is shorter than readLen, then we zero-out
             * the remaining buffer. This is because the MPI_Bcast below
             * broadcasts a buffer of a fixed size, gbp->size. Without zeroing
             * out, valgrind will complain about the uninitialized values.
             */
            if (get_size < readLen)
                memset(readBuf + get_size, 0, readLen - get_size);
        }
        /* only root process reads file header, keeps track of current read
         * file pointer location */
        gbp->offset += readLen;
    }

    if (gbp->safe_mode == 1) {
        TRACE_COMM(MPI_Bcast)(&err, 1, MPI_INT, 0, gbp->comm);
        if (err != NC_NOERR) return err;
    }

    /* broadcast root's read (full or partial header) to other processes */
    TRACE_COMM(MPI_Bcast)(gbp->base, gbp->size, MPI_BYTE, 0, gbp->comm);

    gbp->pos = gbp->base;

    return err;
}

/*----< hdr_get_uint32() >---------------------------------------------------*/
/* in CDF-1 format, all integers are 32-bit
 * in CDF-2 format, only variable begin (starting file offset) is 64-bit
 * in CDF-5 format, both variable's begin and size are 64-bit
 */
static int
hdr_get_uint32(bufferinfo *gbp, uint *xp)
{
    int err;

    if (gbp->pos + 4 > gbp->end) {
        err = hdr_fetch(gbp);
        if (err != NC_NOERR) return err;
    }

    err = ncmpix_get_uint32((const void **)(&gbp->pos), xp);
    return err;
}

/*----< hdr_get_uint64() >---------------------------------------------------*/
/* in CDF-1 format, all integers are 32-bit
 * in CDF-2 format, only variable begin (starting file offset) is 64-bit
 * in CDF-5 format, both variable's begin and size are 64-bit
 */
static int
hdr_get_uint64(bufferinfo *gbp, uint64 *xp)
{
    int err;

    if (gbp->pos + 8 > gbp->end) {
        err = hdr_fetch(gbp);
        if (err != NC_NOERR) return err;
    }

    err = ncmpix_get_uint64((const void **)(&gbp->pos), xp);
    return err;
}

/*----< hdr_get_NC_tag() >---------------------------------------------------*/
/* NC_tag is 32-bit integer and can be the followings:
 *     ZERO (NC_UNSPECIFIED)
 *     NC_DIMENSION
 *     NC_ATTRIBUTE
 *     NC_VARIABLE
 */
static int
hdr_get_NC_tag(bufferinfo *gbp, NC_tag *tagp)
{
    int err;
    uint tag;

    if (gbp->pos + 4 > gbp->end) {
        err = hdr_fetch(gbp);
        if (err != NC_NOERR) return err;
    }

    /* get an external unsigned 4-byte integer from the file */
    err = ncmpix_get_uint32((const void **)(&gbp->pos), &tag);
    if (err != NC_NOERR) return err;

    *tagp = (NC_tag) tag;
    return NC_NOERR;
}

/*----< hdr_get_nc_type() >--------------------------------------------------*/
static int
hdr_get_nc_type(bufferinfo *gbp, nc_type *xtypep)
{
    /* nc_type is 4-byte integer */
    int err;
    uint xtype;

    if (gbp->pos + 4 > gbp->end) {
        err = hdr_fetch(gbp);
        if (err != NC_NOERR) return err;
    }

    err = ncmpix_get_uint32((const void **)(&gbp->pos), &xtype);
    if (err != NC_NOERR) return err;

    /* check if xtype is within legal ranges of CDF-1/2/5 formats */
    if (xtype < NC_BYTE)
        DEBUG_RETURN_ERROR(NC_EBADTYPE)

    if (gbp->version < 5) {
        if (xtype > NC_DOUBLE)
            DEBUG_RETURN_ERROR(NC_EBADTYPE)
    }
    else if (xtype > NC_UINT64)
        DEBUG_RETURN_ERROR(NC_EBADTYPE)

    *xtypep = (nc_type) xtype;
    return NC_NOERR;
}

/*----< hdr_get_NC_name() >--------------------------------------------------*/
static int
hdr_get_NC_name(bufferinfo *gbp, char **namep)
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
    int err=NC_NOERR, nchars, padding, bufremain, strcount;
    char *cpos;

    *namep = NULL;

    /* get nelems (string length of name) */
    if (gbp->version < 5) {
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        if (err != NC_NOERR) return err;
        if (tmp > NC_MAX_NAME) DEBUG_RETURN_ERROR(NC_EMAXNAME)
        nchars = (int)tmp;
    }
    else {
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        if (err != NC_NOERR) return err;
        if (tmp > NC_MAX_NAME) DEBUG_RETURN_ERROR(NC_EMAXNAME)
        nchars = (int)tmp;
    }

    /* Allocate a NC_string structure large enough to hold nchars characters.
     * Note nchars is strlen(namestring) without terminal character.
     */
    *namep = (char*) NCI_Malloc((size_t)nchars + 1);
    if (*namep == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    (*namep)[nchars] = '\0'; /* add terminal character */

    /* X_SIZEOF_CHAR is defined as 1 in classical CDF formats
    padding = _RNDUP(X_SIZEOF_CHAR * nchars, X_ALIGN) - X_SIZEOF_CHAR * nchars;
    */
    padding = _RNDUP(nchars, X_ALIGN) - nchars;

    bufremain = gbp->size - (gbp->pos - gbp->base);

    cpos = *namep;

    /* get namestring with padding (the space in file allocated for string
     * namestring is upward aligned with 4 bytes */
    while (nchars > 0) {
        if (bufremain > 0) {
            strcount = MIN(bufremain, nchars);
            memcpy(cpos, gbp->pos, (size_t)strcount);
            nchars -= strcount;
            gbp->pos += strcount;
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
        if (gbp->pos + padding > gbp->end) {
            err = hdr_fetch(gbp);
            if (err != NC_NOERR) return err;
        }

        /* CDF specification: Header padding uses null (\x00) bytes.
         * However, prior to version 4.5.0, NetCDF did not implement this
         * specification entirely. In particular, it has never enforced the
         * null-byte padding for attribute values (it has for others, such as
         * names of dimension, variables, and attributes.) It also appears that
         * files created by SciPy NetCDF module or NetCDF Java module, both
         * developed independent from NetCDF-C, also fail to respect this
         * padding specification.  This becomes a problem for PnetCDF to read
         * such netCDF files, because PnetCDF enforces the header padding from
         * its very first release.  The files violating the padding
         * specification will not be readable by PnetCDF of all releases prior
         * to 1.9.0 and error code NC_EINVAL or NC_ENOTNC will be thrown when
         * opening such files.  Note if the sizes of all attribute values of
         * your files are aligned with 4-byte boundaries, then the files are
         * readable by PnetCDF.  In order to keep the files in question
         * readable by PnetCDF, checking for null-byte padding has been
         * disabled in 1.9.0. But, we keep this checking in ncvalidator, a
         * utility program that can report whether a CDF file violates the file
         * format specification, including this null-byte padding. See r3516
         * and discussion in NetCDF Github issue
         * https://github.com/Unidata/netcdf-c/issues/657.
         */
#ifdef ENABLE_NULL_BYTE_HEADER_PADDING
        char pad[X_ALIGN-1];
        memset(pad, 0, X_ALIGN-1);
        if (memcmp(gbp->pos, pad, (size_t)padding) != 0) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header corrupted, non-zero padding found\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_ASSIGN_ERROR(err, NC_ENULLPAD) /* not a fatal error */
        }
#endif
        gbp->pos += padding;
    }

    return err;
}

/*----< hdr_get_NC_dim() >---------------------------------------------------*/
static int
hdr_get_NC_dim(bufferinfo *gbp, int unlimited_id, NC_dim **dimpp)
{
    /* netCDF file format:
     *  ...
     * dim        = name  dim_length
     * dim_length = NON_NEG
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    int err, status=NC_NOERR;
    char *name;
    NC_dim *dimp;
    MPI_Offset dim_length;

    *dimpp = NULL;

    /* get name */
    err = hdr_get_NC_name(gbp, &name);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) return err;

    /* get dim_length */
    if (gbp->version < 5) {
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        dim_length = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        dim_length = (MPI_Offset)tmp;
    }
    if (err != NC_NOERR) { /* free space allocated for name */
        NCI_Free(name);
        return err;
    }

    /* check if unlimited_id already set */
    if (unlimited_id != -1 && dim_length == 0) {
        NCI_Free(name);
        return NC_EUNLIMIT;
    }

    /* allocate and initialize NC_dim object */
    dimp = (NC_dim*) NCI_Malloc(sizeof(NC_dim));
    if (dimp == NULL) {
        NCI_Free(name);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    dimp->name     = name;
    dimp->name_len = strlen(name);
    dimp->size     = dim_length;

    *dimpp = dimp;

    return status;
}

/*----< hdr_get_NC_dimarray() >----------------------------------------------*/
/* For CDF-5 format, nelems (number of dimensions) is of type non-negative
 * INT64. However, argument ndims/dimid in all PnetCDF/NetCDF APIs are of type
 * int. Thus, we only need to use type int for internal metadata, ndefined. If
 * nelems in the input file is more than NC_MAX_DIMS, then it violates the
 * format specifications (NC_EMAXDIMS).
 */
static int
hdr_get_NC_dimarray(bufferinfo *gbp, NC_dimarray *ncap)
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
    int i, err, status=NC_NOERR, ndefined=0;
    size_t alloc_size;
    NC_tag tag = NC_UNSPECIFIED;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* read NC_tag (NC_DIMENSION or ZERO) from gbp buffer */
    err = hdr_get_NC_tag(gbp, &tag);
    if (err != NC_NOERR) return err;

    /* read nelems (number of dimensions) from gbp buffer */
    if (gbp->version < 5) { /* nelems is <non-negative INT> */
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        if (err != NC_NOERR) return err;
        /* cannot be more than max number of dimensions */
        if (tmp > NC_MAX_DIMS) DEBUG_RETURN_ERROR(NC_EMAXDIMS)
        ndefined = (int)tmp;
    }
    else { /* nelems is <non-negative INT64> */
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        if (err != NC_NOERR) return err;
        /* cannot be more than max number of dimensions */
        if (tmp > NC_MAX_DIMS) DEBUG_RETURN_ERROR(NC_EMAXDIMS)
        ndefined = (int)tmp;
    }

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
    ncap->value = (NC_dim**) NCI_Calloc(alloc_size, sizeof(NC_dim*));
    if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    for (i=0; i<ndefined; i++) {
        err = hdr_get_NC_dim(gbp, ncap->unlimited_id, ncap->value + i);
        if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
        else if (err != NC_NOERR) { /* error: fail to get the next dim */
            ncmpio_free_NC_dimarray(ncap);
            return err;
        }
        if (ncap->value[i]->size == NC_UNLIMITED)
            ncap->unlimited_id = i; /* ID of unlimited dimension */
    }

    return status;
}

/*----< hdr_get_NC_attrV() >-------------------------------------------------*/
static int
hdr_get_NC_attrV(bufferinfo *gbp, NC_attr *attrp)
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
    int err=NC_NOERR, xsz, padding, bufremain;
    void *value = attrp->xvalue;
    MPI_Offset nbytes;

    ncmpii_xlen_nc_type(attrp->xtype, &xsz);
    nbytes = attrp->nelems * xsz;
    padding = attrp->xsz - nbytes;

    bufremain = gbp->size - (gbp->pos - gbp->base);
    /* gbp->size is the read chunk size, which is of type 4-byte int.
     * thus bufremain should be less than INT_MAX */

    /* get values */
    while (nbytes > 0) {
        if (bufremain > 0) {
            int attcount = MIN(nbytes, bufremain);
            memcpy(value, gbp->pos, attcount);
            nbytes -= attcount;
            gbp->pos += attcount;
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
        if (gbp->pos + padding > gbp->end) {
            err = hdr_fetch(gbp);
            if (err != NC_NOERR) return err;
        }

        /* CDF specification: Header padding uses null (\x00) bytes.
         * However, prior to version 4.5.0, NetCDF did not implement this
         * specification entirely. In particular, it has never enforced the
         * null-byte padding for attribute values (it has for others, such as
         * names of dimension, variables, and attributes.) It also appears that
         * files created by SciPy NetCDF module or NetCDF Java module, both
         * developed independent from NetCDF-C, also fail to respect this
         * padding specification.  This becomes a problem for PnetCDF to read
         * such netCDF files, because PnetCDF enforces the header padding from
         * its very first release.  The files violating the padding
         * specification will not be readable by PnetCDF of all releases prior
         * to 1.9.0 and error code NC_EINVAL or NC_ENOTNC will be thrown when
         * opening such files.  Note if the sizes of all attribute values of
         * your files are aligned with 4-byte boundaries, then the files are
         * readable by PnetCDF.  In order to keep the files in question
         * readable by PnetCDF, checking for null-byte padding has been
         * disabled in 1.9.0. But, we keep this checking in ncvalidator, a
         * utility program that can report whether a CDF file violates the file
         * format specification, including this null-byte padding. See r3516
         * and discussion in NetCDF Github issue
         * https://github.com/Unidata/netcdf-c/issues/657.
         */
#ifdef ENABLE_NULL_BYTE_HEADER_PADDING
        char pad[X_ALIGN-1];
        memset(pad, 0, X_ALIGN-1);
        if (memcmp(gbp->pos, pad, (size_t)padding) != 0) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"Error in file %s func %s line %d: NetCDF header corrupted, non-zero padding found\n",__FILE__,__func__,__LINE__);
#endif
            DEBUG_ASSIGN_ERROR(err, NC_ENULLPAD)
        }
#endif
        gbp->pos += padding;
    }
    return err;
}

/*----< hdr_get_NC_attr() >--------------------------------------------------*/
static int
hdr_get_NC_attr(bufferinfo *gbp, NC_attr **attrpp)
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     * nc_type = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * nelems  = NON_NEG       // number of elements in following sequence
     * NON_NEG = <non-negative INT> |  // CDF-1 and CDF-2
     *           <non-negative INT64>  // CDF-5
     */
    int err, status=NC_NOERR;
    char *name;
    nc_type type;
    MPI_Offset nelems;
    NC_attr *attrp;

    /* get name */
    err = hdr_get_NC_name(gbp, &name);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) return err;

    /* get nc_type */
    err = hdr_get_nc_type(gbp, &type);
    if (err != NC_NOERR) {
        NCI_Free(name);
        return err;
    }

    /* get nelems */
    if (gbp->version < 5) {
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        nelems = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        nelems = (MPI_Offset)tmp;
    }
    if (err != NC_NOERR) {
        NCI_Free(name);
        return err;
    }

    /* allocate space for attribute object, name will be assigned in attrp */
    err = ncmpio_new_NC_attr(name, type, nelems, &attrp);
    if (err != NC_NOERR) {
        NCI_Free(name);
        return err;
    }

    /* get [values ...] */
    err = hdr_get_NC_attrV(gbp, attrp);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) {
        ncmpio_free_NC_attr(attrp);
        NCI_Free(attrp);
        return err;
    }

    *attrpp = attrp;

    return status;
}

/*----< hdr_get_NC_attrarray() >---------------------------------------------*/
/* For CDF-5 format, nelems (number of attributes) is of type non-negative
 * INT64. However, argument nattrs in all PnetCDF/NetCDF APIs are of type int.
 * Thus, we only need to use type int for internal metadata, ndefined. If
 * nelems in the input file is more than NC_MAX_ATTRS, then it violates the
 * format specifications (NC_EMAXATTS).
 */
static int
hdr_get_NC_attrarray(bufferinfo *gbp, NC_attrarray *ncap)
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
    int i, err, status=NC_NOERR, ndefined=0;
    size_t alloc_size;
    NC_tag tag = NC_UNSPECIFIED;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* read NC_tag (NC_ATTRIBUTE or ZERO) from gbp buffer */
    err = hdr_get_NC_tag(gbp, &tag);
    if (err != NC_NOERR) return err;

    /* read nelems (number of attributes) from gbp buffer */
    if (gbp->version < 5) { /* nelems is <non-negative INT> */
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        if (err != NC_NOERR) return err;
        /* cannot be more than max number of attributes */
        if (tmp > NC_MAX_ATTRS) DEBUG_RETURN_ERROR(NC_EMAXATTS)
        ndefined = (int)tmp;
    }
    else { /* nelems is <non-negative INT64> */
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        if (err != NC_NOERR) return err;
        /* cannot be more than max number of attributes */
        if (tmp > NC_MAX_ATTRS) DEBUG_RETURN_ERROR(NC_EMAXATTS)
        ndefined = (int)tmp;
    }

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
    ncap->value = (NC_attr**) NCI_Calloc(alloc_size, sizeof(NC_attr*));
    if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* get [attr ...] */
    for (i=0; i<ndefined; i++) {
        err = hdr_get_NC_attr(gbp, ncap->value + i);
        if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
        else if (err != NC_NOERR) { /* Error: fail to get the next att */
            ncmpio_free_NC_attrarray(ncap);
            return err;
        }
    }

    return status;
}

/*----< hdr_get_NC_var() >---------------------------------------------------*/
static int
hdr_get_NC_var(bufferinfo  *gbp,
               NC_var     **varpp,
               int          f_ndims) /* no. dimensions defined in file */
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
    int dim, ndims, err, status=NC_NOERR;
    char *name;
    NC_var *varp;

    /* get name */
    err = hdr_get_NC_name(gbp, &name);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) return err;

    /* nelems (number of dimensions) */
    if (gbp->version < 5) {
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        if (err != NC_NOERR) {
            NCI_Free(name);
            return err;
        }
        /* cannot be more than NC_MAX_VAR_DIMS */
        if (tmp > NC_MAX_VAR_DIMS) {
            NCI_Free(name);
            DEBUG_RETURN_ERROR(NC_EMAXDIMS)
        }
        ndims = (int)tmp;
    }
    else {
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        if (err != NC_NOERR) {
            NCI_Free(name);
            return err;
        }
        /* cannot be more than NC_MAX_VAR_DIMS */
        if (tmp > NC_MAX_VAR_DIMS) {
            NCI_Free(name);
            DEBUG_RETURN_ERROR(NC_EMAXDIMS)
        }
        ndims = (int)tmp;
    }

    /* allocate space for NC_var object */
    varp = ncmpio_new_NC_var(name, ndims);
    if (varp == NULL) {
        NCI_Free(name);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    /* get [dimid ...] */
    for (dim=0; dim<ndims; dim++) {
        if (gbp->version < 5) {
            uint tmp;
            err = hdr_get_uint32(gbp, &tmp);
            if (err != NC_NOERR) break;
            /* dimid should be < f_ndims (num of dimensions defined in file) */
            if (tmp >= f_ndims) {
                DEBUG_ASSIGN_ERROR(err, NC_EBADDIM)
                goto fn_exit;
            }
            varp->dimids[dim] = (int)tmp;
        }
        else {
            uint64 tmp;
            err = hdr_get_uint64(gbp, &tmp);
            if (err != NC_NOERR) break;
            /* dimid should be < f_ndims (num of dimensions defined in file) */
            if (tmp >= f_ndims) {
                DEBUG_ASSIGN_ERROR(err, NC_EBADDIM)
                goto fn_exit;
            }
            varp->dimids[dim] = (int)tmp;
        }
    }

    /* get vatt_list */
    err = hdr_get_NC_attrarray(gbp, &varp->attrs);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) goto fn_exit;

    /* get nc_type */
    err = hdr_get_nc_type(gbp, &varp->xtype);
    if (err != NC_NOERR) goto fn_exit;

    ncmpii_xlen_nc_type(varp->xtype, &varp->xsz);

    /* get vsize */
    if (gbp->version < 5) {
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        varp->len = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        varp->len = (MPI_Offset)tmp;
    }
    if (err != NC_NOERR) goto fn_exit;

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
       header. Its value will be recalculated immediately after read from file.
     */

    /* get begin */
    if (gbp->version == 1) {
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        varp->begin = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        varp->begin = (MPI_Offset)tmp;
    }

fn_exit:
    if (err != NC_NOERR)
        ncmpio_free_NC_var(varp);
    else
        *varpp = varp;

    return (err == NC_NOERR) ? status : err;
}

/*----< hdr_get_NC_vararray() >----------------------------------------------*/
/* For CDF-5 format, nelems (number of variables) is of type non-negative
 * INT64. However, argument nvars/varid in all PnetCDF/NetCDF APIs are of type
 * int. Thus, we only need to use type int for internal metadata, ndefined. If
 * nelems in the input file is more than NC_MAX_VARS, then it violates the
 * format specifications (NC_EMAXVARS).
 */
static int
hdr_get_NC_vararray(bufferinfo  *gbp,
                    NC_vararray *ncap,
                    int          f_ndims) /* no. dimensions defined in file */
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
    int i, err, status=NC_NOERR, ndefined=0;
    size_t alloc_size;
    NC_tag tag = NC_UNSPECIFIED;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* read NC_tag (NC_VARIABLE or ZERO) from gbp buffer */
    err = hdr_get_NC_tag(gbp, &tag);
    if (err != NC_NOERR) return err;

    /* read nelems (number of variables) from gbp buffer */
    if (gbp->version < 5) { /* nelems is <non-negative INT> */
        uint tmp;
        err = hdr_get_uint32(gbp, &tmp);
        if (err != NC_NOERR) return err;
        /* cannot be more than max number of attributes */
        if (tmp > NC_MAX_VARS) DEBUG_RETURN_ERROR(NC_EMAXVARS)
        ndefined = (int)tmp;
    }
    else { /* nelems is <non-negative INT64> */
        uint64 tmp;
        err = hdr_get_uint64(gbp, &tmp);
        if (err != NC_NOERR) return err;
        /* cannot be more than max number of attributes */
        if (tmp > NC_MAX_VARS) DEBUG_RETURN_ERROR(NC_EMAXVARS)
        ndefined = (int)tmp;
    }

    /* Now ndefined is in between 0 and NC_MAX_VARS */
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
    ncap->value = (NC_var**) NCI_Calloc(alloc_size, sizeof(NC_var*));
    if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* get [var ...] */
    for (i=0; i<ndefined; i++) {
        err = hdr_get_NC_var(gbp, ncap->value + i, f_ndims);
        if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
        else if (err != NC_NOERR) { /* Error: fail to get the next var */
            ncmpio_free_NC_vararray(ncap);
            return err;
        }
        ncap->value[i]->varid = i;
    }

    return status;
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
    int i, err, status=NC_NOERR;
    bufferinfo getbuf;
    char magic[NC_MAGIC_LEN];

    assert(ncp != NULL);

    /* Initialize the get buffer that stores the header read from the file */
    getbuf.comm          = ncp->comm;
    getbuf.collective_fh = ncp->collective_fh;
    getbuf.get_size      = 0;
    getbuf.offset        = 0;   /* read from start of the file */
    getbuf.safe_mode     = ncp->safe_mode;

    /* CDF-5's minimum header size is 4 bytes more than CDF-1 and CDF-2's */
    getbuf.size = _RNDUP( MAX(MIN_NC_XSZ+4, ncp->chunk), X_ALIGN );

    getbuf.base = (char*) NCI_Malloc((size_t)getbuf.size);
    getbuf.pos  = getbuf.base;
    getbuf.end  = getbuf.base + getbuf.size;

    /* Fetch the next header chunk. The chunk is 'gbp->size' bytes big */
    err = hdr_fetch(&getbuf);
    if (err != NC_NOERR) return err;

    /* processing the header from getbuf, the get buffer */

    /* First get the file format information, magic */
    err = ncmpix_getn_text((const void **)(&getbuf.pos), NC_MAGIC_LEN, magic);
    if (err != NC_NOERR) return err;

    /* check if the first three bytes are 'C','D','F' */
    if (memcmp(magic, "CDF", 3) != 0) {
        /* check if is HDF5 file */
        char signature[8], *hdf5_signature="\211HDF\r\n\032\n";
        ncmpix_getn_text((const void **)(&getbuf.pos), 8, signature);
        if (memcmp(signature, hdf5_signature, 8) == 0) {
            DEBUG_ASSIGN_ERROR(err, NC_ENOTNC3)
            if (ncp->safe_mode)
                fprintf(stderr,"Error: file %s is HDF5 format\n",ncp->path);
        }
        else
            DEBUG_ASSIGN_ERROR(err, NC_ENOTNC)
        goto fn_exit;
    }

    /* check version number in last byte of magic */
    if (magic[3] == 0x1) {
        getbuf.version = ncp->format = 1;
    } else if (magic[3] == 0x2) {
        getbuf.version = ncp->format = 2;
    } else if (magic[3] == 0x5) {
        getbuf.version = ncp->format = 5;
    } else {
        NCI_Free(getbuf.base);
        DEBUG_RETURN_ERROR(NC_ENOTNC) /* not a netCDF file */
    }

    /* get numrecs from getbuf into ncp */
    if (getbuf.version < 5) {
        uint tmp=0;
        err = hdr_get_uint32(&getbuf, &tmp);
        if (err != NC_NOERR) goto fn_exit;
        ncp->numrecs = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp=0;
        err = hdr_get_uint64(&getbuf, &tmp);
        if (err != NC_NOERR) goto fn_exit;
        ncp->numrecs = (MPI_Offset)tmp;
    }

    assert(getbuf.pos < getbuf.end);

    /* get dim_list from getbuf into ncp */
    err = hdr_get_NC_dimarray(&getbuf, &ncp->dims);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) goto fn_exit;

    /* get gatt_list from getbuf into ncp */
    err = hdr_get_NC_attrarray(&getbuf, &ncp->attrs);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) goto fn_exit;

    /* get var_list from getbuf into ncp */
    err = hdr_get_NC_vararray(&getbuf, &ncp->vars, ncp->dims.ndefined);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) goto fn_exit;

    /* get the un-aligned size occupied by the file header */
    ncp->xsz = ncmpio_hdr_len_NC(ncp);

    /* Recompute the shapes of all variables (shape, xsz, len)
     * Sets ncp->begin_var to start of first variable.
     * Sets ncp->begin_rec to start of first record variable.
     */
    err = compute_var_shape(ncp);
    if (err != NC_NOERR) goto fn_exit;

    /* update the total number of record variables --------------------------*/
    ncp->vars.num_rec_vars = 0;
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.num_rec_vars += IS_RECVAR(ncp->vars.value[i]);

    /* Check whether variable sizes are legal for the given file format */
    err = ncmpio_NC_check_vlens(ncp);
    if (err != NC_NOERR) goto fn_exit;

    /* Check whether variable begins are in an increasing order.
     * Adding this check here is necessary for detecting corrupted metadata. */
    err = ncmpio_NC_check_voffs(ncp);
    if (err != NC_NOERR) goto fn_exit;

fn_exit:
    ncp->get_size += getbuf.get_size;
    NCI_Free(getbuf.base);

    return (err == NC_NOERR) ? status : err;
}

