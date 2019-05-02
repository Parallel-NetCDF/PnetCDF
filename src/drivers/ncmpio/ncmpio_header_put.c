/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <string.h>  /* memcpy() */
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncx.h>
#include "ncmpio_NC.h"

/*
 * "magic number" at beginning of file: 0x43444601 (big endian)
 */
static const char ncmagic1[] = {'C', 'D', 'F', 0x01};
static const char ncmagic2[] = {'C', 'D', 'F', 0x02};
static const char ncmagic5[] = {'C', 'D', 'F', 0x05};

/*----< hdr_put_NC_name() >--------------------------------------------------*/
static int
hdr_put_NC_name(bufferinfo *pbp,
                const char *name)
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
    size_t nchars = strlen(name);

    /* copy nelems */
    if (pbp->version < 5)
        err = ncmpix_put_uint32((void**)(&pbp->pos), (uint)nchars);
    else
        err = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)nchars);
    if (err != NC_NOERR) return err;

    /* copy namestring */
    return ncmpix_pad_putn_text((void **)(&pbp->pos), (MPI_Offset)nchars, name);
}

/*----< hdr_put_NC_dim() >---------------------------------------------------*/
static int
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
static int
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
    int xsz;
    MPI_Offset padding, sz;

    /* ncmpii_xlen_nc_type() returns the element size (unaligned) of
     * attrp->xtype attrp->xsz is the aligned total size of attribute values
     */
    ncmpii_xlen_nc_type(attrp->xtype, &xsz);
    sz = attrp->nelems * xsz;
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
static int
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
    status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)attrp->xtype);
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
    status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)varp->xtype);
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
        status = ncmpix_putn_text((void **)(&putbuf.pos), sizeof(ncmagic5), ncmagic5);
    }
    else if (ncp->format == 2) {
        putbuf.version = 2;
        status = ncmpix_putn_text((void **)(&putbuf.pos), sizeof(ncmagic2), ncmagic2);
    }
    else {
        putbuf.version = 1;
        status = ncmpix_putn_text((void **)(&putbuf.pos), sizeof(ncmagic1), ncmagic1);
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
        size_t bufLen = _RNDUP(ncp->xsz, X_ALIGN);
        void *buf = NCI_Malloc(bufLen); /* header's write buffer */

        /* copy header object to write buffer */
        status = ncmpio_hdr_put_NC(ncp, buf);

        if (ncp->xsz != (int)ncp->xsz) {
            NCI_Free(buf);
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }

#ifdef _USE_MPI_GET_COUNT
        /* explicitly initialize mpistatus object to 0. For zero-length read,
         * MPI_Get_count may report incorrect result for some MPICH version,
         * due to the uninitialized MPI_Status object passed to MPI-IO calls.
         * Thus we initialize it above to work around.
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));
#endif
        TRACE_IO(MPI_File_write_at)(fh, 0, buf, (int)ncp->xsz, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_write_at");
            if (status == NC_NOERR) {
                err = (err == NC_EFILE) ? NC_EWRITE : err;
                DEBUG_ASSIGN_ERROR(status, err)
            }
        }
        else {
#ifdef _USE_MPI_GET_COUNT
            int put_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
            ncp->put_size += put_size;
#else
            ncp->put_size += ncp->xsz;
#endif
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
            ncmpii_error_mpi2nc(mpireturn,"MPI_Bcast");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
    }

    if (NC_doFsync(ncp)) { /* NC_SHARE is set */
        TRACE_IO(MPI_File_sync)(fh);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_error_mpi2nc(mpireturn,"MPI_File_sync");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
        TRACE_COMM(MPI_Barrier)(ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_error_mpi2nc(mpireturn,"MPI_Barrier");
            DEBUG_RETURN_ERROR(NC_EMPI)
        }
    }

    return status;
}

