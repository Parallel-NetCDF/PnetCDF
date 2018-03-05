/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>  /* strcasecmp() */
#include <limits.h>   /* INT_MAX */
#include <assert.h>
#include <errno.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< ncmpio_set_pnetcdf_hints() >-----------------------------------------*/
/* this is where the I/O hints designated to pnetcdf are extracted */
void ncmpio_set_pnetcdf_hints(NC *ncp, MPI_Info info)
{
    char value[MPI_MAX_INFO_VAL];
    int  flag;

    if (info == MPI_INFO_NULL) return;

    /* nc_header_align_size, nc_var_align_size, and r_align * take effect when
     * a file is created or opened and later adding more header or variable
     * data */

    /* extract PnetCDF hints from user info object */
    MPI_Info_get(info, "nc_header_align_size", MPI_MAX_INFO_VAL-1, value,
                 &flag);
    if (flag) {
        errno = 0;  /* errno must set to zero before calling strtoll */
        ncp->h_align = strtoll(value,NULL,10);
        if (errno != 0) ncp->h_align = 0;
        else if (ncp->h_align < 0) ncp->h_align = 0;
    }

    MPI_Info_get(info, "nc_var_align_size", MPI_MAX_INFO_VAL-1, value, &flag);
    if (flag) {
        errno = 0;  /* errno must set to zero before calling strtoll */
        ncp->v_align = strtoll(value,NULL,10);
        if (errno != 0) ncp->v_align = 0;
        else if (ncp->v_align < 0) ncp->v_align = 0;
    }

    MPI_Info_get(info, "nc_record_align_size", MPI_MAX_INFO_VAL-1, value,
                 &flag);
    if (flag) {
        errno = 0;  /* errno must set to zero before calling strtoll */
        ncp->r_align = strtoll(value,NULL,10);
        if (errno != 0) ncp->r_align = 0;
        else if (ncp->r_align < 0) ncp->r_align = 0;
    }

    /* get header reading chunk size from info */
    MPI_Info_get(info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1, value,
                 &flag);
    if (flag) {
        errno = 0;  /* errno must set to zero before calling strtoll */
        ncp->chunk = (int) strtol(value,NULL,10);
        if (errno != 0) ncp->chunk = 0;
        else if (ncp->chunk < 0) ncp->chunk = 0;
    }

    /* hint on setting in-place byte swap (matters only for Little Endian) */
    MPI_Info_get(info, "nc_in_place_swap", MPI_MAX_INFO_VAL-1, value, &flag);
    if (flag) {
        if (strcasecmp(value, "enable") == 0) {
            fClr(ncp->flags, NC_MODE_SWAP_OFF);
            fSet(ncp->flags, NC_MODE_SWAP_ON);
        }
        else if (strcasecmp(value, "disable") == 0) {
            fClr(ncp->flags, NC_MODE_SWAP_ON);
            fSet(ncp->flags, NC_MODE_SWAP_OFF);
        }
        else if (strcasecmp(value, "auto") == 0) {
            fClr(ncp->flags, NC_MODE_SWAP_ON);
            fClr(ncp->flags, NC_MODE_SWAP_OFF);
        }
    }

#ifdef ENABLE_SUBFILING
    MPI_Info_get(info, "pnetcdf_subfiling", MPI_MAX_INFO_VAL-1, value, &flag);
    if (flag && strcasecmp(value, "enable") == 0)
        ncp->subfile_mode = 1;

    MPI_Info_get(info, "nc_num_subfiles", MPI_MAX_INFO_VAL-1, value, &flag);
    if (flag) {
        errno = 0;
        ncp->num_subfiles = strtoll(value,NULL,10);
        if (errno != 0) ncp->num_subfiles = 0;
        else if (ncp->num_subfiles < 0) ncp->num_subfiles = 0;
    }
    if (ncp->subfile_mode == 0) ncp->num_subfiles = 0;
#endif
}

/*----< ncmpio_first_offset() >-----------------------------------------------*/
/* Returns the file offset of the first variable element accessed by this
 * request. Note zero-length request should never call this subroutine.
 */
int
ncmpio_first_offset(const NC         *ncp,
                    const NC_var     *varp,
                    const MPI_Offset  start[],  /* [varp->ndims] */
                    MPI_Offset       *offset)   /* OUT: file offset */
{
    int i, ndims;

    ndims = varp->ndims; /* number of dimensions of this variable */

    if (ndims == 0) { /* scalar variable */
        *offset = varp->begin;
        return NC_NOERR;
    }

    *offset = 0;
    if (IS_RECVAR(varp)) {
        if (ndims > 1) *offset += start[ndims - 1];
        for (i=1; i<ndims-1; i++)
            *offset += start[i] * varp->dsizes[i+1];
        *offset *= varp->xsz;  /* multiply element size */
        *offset += start[0] * ncp->recsize;
    }
    else {
        if (ndims > 1) *offset += start[0] * varp->dsizes[1];
        for (i=1; i<ndims-1; i++)
            *offset += start[i] * varp->dsizes[i+1];
        *offset += start[ndims-1];
        *offset *= varp->xsz;  /* multiply element size */
    }

    *offset += varp->begin; /* beginning file offset of this variable */

    return NC_NOERR;
}

/*----< ncmpio_last_offset() >-----------------------------------------------*/
/* Returns the file offset of the last variable element accessed by this
 * request.
 * If count is NULL, this is equivalent to find the starting offset of this
 * request. Note zero-length request should never call this subroutine.
 */
int
ncmpio_last_offset(const NC         *ncp,
                   const NC_var     *varp,
                   const MPI_Offset  start[],   /* [varp->ndims] */
                   const MPI_Offset  count[],   /* [varp->ndims] */
                   const MPI_Offset  stride[],  /* [varp->ndims] */
                   MPI_Offset       *offset_ptr) /* OUT: file offset */
{
    int i, ndims;
    MPI_Offset offset, *last_indx=NULL;

    offset = varp->begin; /* beginning file offset of this variable */
    ndims  = varp->ndims; /* number of dimensions of this variable */

    if (ndims == 0) { /* scalar variable */
        *offset_ptr = varp->begin;
        return NC_NOERR;
    }

    /* when count == NULL, this is called from a var API */
    if (count != NULL) {
        last_indx = (MPI_Offset*) NCI_Malloc((size_t)ndims * SIZEOF_MPI_OFFSET);

        if (stride != NULL) {
            for (i=0; i<ndims; i++) {
                assert(count[i] > 0);
                last_indx[i] = start[i] + (count[i] - 1) * stride[i];
            }
        }
        else { /* stride == NULL */
            for (i=0; i<ndims; i++) {
                assert(count[i] > 0);
                last_indx[i] = start[i] + count[i] - 1;
            }
        }
    }
    else { /* when count == NULL stride is of no use */
        last_indx = (MPI_Offset*) start;
    }

    /* check NC_EINVALCOORDS and NC_EEDGE already done in dispatchers/var_getput.m4 */
#if 0
    /* check whether last_indx is valid */

    int firstDim = 0;
    /* check NC_EINVALCOORDS for record dimension */
    if (varp->shape[0] == NC_UNLIMITED) {
        if (ncp->format < 5 && last_indx[0] > NC_MAX_UINT) { /* CDF-1 and 2 */
            if (count != NULL) NCI_Free(last_indx);
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
        }
        /* for record variable, [0] is the NC_UNLIMITED dimension */
        if (fIsSet(reqMode, NC_REQ_RD) && last_indx[0] >= ncp->numrecs) {
            /* read cannot go beyond current numrecs */
            if (count != NULL) NCI_Free(last_indx);
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
        }
        firstDim = 1; /* done for checking the record dimension */
    }
    /* continue to check NC_EINVALCOORDS for the rest dimensions */
    for (i=firstDim; i<ndims; i++) {
        if (last_indx[i] < 0 || last_indx[i] >= varp->shape[i]) {
            if (count != NULL) NCI_Free(last_indx);
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
        }
    }
#endif

    if (varp->shape[0] == NC_UNLIMITED)
        offset += last_indx[0] * ncp->recsize;
    else
        offset += last_indx[ndims-1] * varp->xsz;

    if (ndims > 1) {
        if (IS_RECVAR(varp))
            offset += last_indx[ndims - 1] * varp->xsz;
        else
            offset += last_indx[0] * varp->dsizes[1] * varp->xsz;

        for (i=1; i<ndims-1; i++)
            offset += last_indx[i] * varp->dsizes[i+1] * varp->xsz;
    }

    if (count != NULL) NCI_Free(last_indx);

    *offset_ptr = offset;
    return NC_NOERR;
}

/*----< ncmpio_access_range() >----------------------------------------------*/
/* Returns the file offsets of access range of this request: starting file
 * offset and end offset (exclusive).
 * Note zero-length request should never call this subroutine.
 */
int
ncmpio_access_range(const NC         *ncp,
                    const NC_var     *varp,
                    const MPI_Offset  start[],   /* [varp->ndims] */
                    const MPI_Offset  count[],   /* [varp->ndims] */
                    const MPI_Offset  stride[],  /* [varp->ndims] */
                    MPI_Offset       *start_off, /* OUT: start offset */
                    MPI_Offset       *end_off)   /* OUT: end   offset */
{
    int i, ndims;
    MPI_Offset *last_indx=NULL;

    ndims = varp->ndims; /* number of dimensions of this variable */

    if (ndims == 0) { /* scalar variable */
        *start_off = varp->begin; /* beginning file offset of this variable */
        *end_off   = varp->begin + varp->xsz;
        return NC_NOERR;
    }

    assert(start != NULL);
    assert(count != NULL);

    /* find the last access index in each dimension */
    last_indx = (MPI_Offset*) NCI_Malloc((size_t)ndims * SIZEOF_MPI_OFFSET);
    if (stride != NULL) {
        for (i=0; i<ndims; i++)
            last_indx[i] = start[i] + (count[i] - 1) * stride[i];
    }
    else { /* stride == NULL */
        for (i=0; i<ndims; i++)
            last_indx[i] = start[i] + count[i] - 1;
    }

    if (IS_RECVAR(varp)) {
        *start_off = 0;
        *end_off = varp->begin + last_indx[0] * ncp->recsize;
        if (ndims > 1) {
            *start_off += start[ndims - 1];
            *end_off   += last_indx[ndims - 1] * varp->xsz;
        }
        for (i=1; i<ndims-1; i++) {
            *start_off += start[i] * varp->dsizes[i+1];
            *end_off   += last_indx[i] * varp->dsizes[i+1] * varp->xsz;
        }
        *start_off *= varp->xsz;   /* multiply element size */
        *start_off += start[0] * ncp->recsize;
        *start_off += varp->begin; /* beginning file offset of this variable */
    }
    else {
        *start_off = 0;
        *end_off = varp->begin + last_indx[ndims-1] * varp->xsz;
        if (ndims > 1) {
            *start_off += start[0] * varp->dsizes[1];
            *end_off += last_indx[0] * varp->dsizes[1] * varp->xsz;
        }
        for (i=1; i<ndims-1; i++) {
            *start_off += start[i] * varp->dsizes[i+1];
            *end_off   += last_indx[i] * varp->dsizes[i+1] * varp->xsz;
        }
        *start_off += start[ndims-1];
        *start_off *= varp->xsz;   /* multiply element size */
        *start_off += varp->begin; /* beginning file offset of this variable */
    }
    NCI_Free(last_indx);

    return NC_NOERR;
}

/*----< ncmpio_pack_xbuf() >-------------------------------------------------*/
/* Pack user buffer, buf, into xbuf, when buftype is non-contiguous or imap
 * is non-contiguous, or type-casting is needed. The immediate buffers, lbuf
 * and cbuf, may be allocated and freed within this subroutine. We try to reuse
 * the intermediate buffers as much as possible. Below describe such design.
 *
 * When called from bput APIs: (abuf means attached buffer pool)
 *     if contig && no imap && no convert
 *         buf   ==   lbuf   ==   cbuf    ==     xbuf memcpy-> abuf
 *                                               abuf
 *     if contig && no imap &&    convert
 *         buf   ==   lbuf   ==   cbuf convert-> xbuf == abuf
 *                                               abuf
 *     if contig &&    imap && no convert
 *         buf   ==   lbuf pack-> cbuf    ==     xbuf == abuf
 *                                abuf
 *     if contig &&    imap &&    convert
 *         buf   ==   lbuf pack-> cbuf convert-> xbuf == abuf
 *                                               abuf
 *  if noncontig && no imap && no convert
 *         buf pack-> lbuf   ==   cbuf    ==     xbuf == abuf
 *                    abuf
 *  if noncontig && no imap &&    convert
 *         buf pack-> lbuf   ==   cbuf convert-> xbuf == abuf
 *                                               abuf
 *  if noncontig &&    imap && no convert
 *         buf pack-> lbuf pack-> cbuf    ==     xbuf == abuf
 *                                abuf
 *  if noncontig &&    imap &&    convert
 *         buf pack-> lbuf pack-> cbuf convert-> xbuf == abuf
 *                                               abuf
 *
 * When called from put/iput APIs:
 *  if contig && no imap && no convert && xbuf_size > NC_BYTE_SWAP_BUFFER_SIZE
 *         buf   ==   lbuf   ==   cbuf    ==     xbuf
 *  if contig && no imap && no convert && xbuf_size <= NC_BYTE_SWAP_BUFFER_SIZE
 *         buf   ==   lbuf   ==   cbuf    ==     xbuf
 *                                                   malloc + memcpy
 *  if contig && no imap &&    convert
 *         buf   ==   lbuf   ==   cbuf convert-> xbuf
 *                                               malloc
 *  if contig &&    imap && no convert
 *         buf   ==   lbuf pack-> cbuf    ==     xbuf
 *                                malloc
 *  if contig &&    imap &&    convert
 *         buf   ==   lbuf pack-> cbuf convert-> xbuf
 *                                malloc         malloc
 *  if noncontig && no imap && no convert
 *         buf pack-> lbuf   ==   cbuf    ==     xbuf
 *                    malloc
 *  if noncontig && no imap &&    convert
 *         buf pack-> lbuf   ==   cbuf convert-> xbuf
 *                    malloc                     malloc
 *  if noncontig &&    imap && no convert
 *         buf pack-> lbuf pack-> cbuf    ==     xbuf
 *                    malloc      malloc
 *  if noncontig &&    imap &&    convert
 *         buf pack-> lbuf pack-> cbuf convert-> xbuf
 *                    malloc      malloc         malloc
 */
int
ncmpio_pack_xbuf(int           fmt,    /* NC_FORMAT_CDF2 NC_FORMAT_CDF5 etc. */
                 NC_var       *varp,
                 MPI_Offset    bufcount,
                 MPI_Datatype  buftype,
                 int           buftype_is_contig,
                 MPI_Offset    nelems, /* no. elements in etype in buf */
                 MPI_Datatype  etype,  /* element type in buftype */
                 MPI_Datatype  imaptype,
                 int           need_convert,
                 int           need_swap,
                 size_t        xbuf_size,
                 void         *buf,    /* user buffer */
                 void         *xbuf)   /* already allocated, in external type */
{
    int err=NC_NOERR, el_size, position;
    void *lbuf=NULL, *cbuf=NULL;
    MPI_Offset ibuf_size;

    /* check byte size of buf (internal representation) */
    MPI_Type_size(etype, &el_size);
    ibuf_size = nelems * el_size;
    if (ibuf_size > INT_MAX) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    /* Step 1: if buftype is not contiguous, i.e. a noncontiguous MPI
     * derived datatype, pack buf into a contiguous buffer, lbuf,
     */
    if (!buftype_is_contig) { /* buftype is not contiguous */
        if (imaptype == MPI_DATATYPE_NULL && !need_convert)
            /* in this case, lbuf will later become xbuf */
            lbuf = xbuf;
        else {
            /* in this case, allocate lbuf and it will be freed before
             * constructing xbuf */
            lbuf = NCI_Malloc((size_t)ibuf_size);
            if (lbuf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        }

        /* pack buf into lbuf based on buftype */
        if (bufcount > INT_MAX) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        position = 0;
        MPI_Pack(buf, (int)bufcount, buftype, lbuf, (int)ibuf_size,
                 &position, MPI_COMM_SELF);
    }
    else /* for contiguous case, we reuse buf */
        lbuf = buf;

    /* Step 2: if imap is non-contiguous, pack lbuf to cbuf */
    if (imaptype != MPI_DATATYPE_NULL) { /* true varm */
        if (!need_convert)
            /* in this case, cbuf will later become xbuf */
            cbuf = xbuf;
        else {
            /* in this case, allocate cbuf and cbuf will be freed before
             * constructing xbuf */
            cbuf = NCI_Malloc((size_t)ibuf_size);
            if (cbuf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        }

        /* pack lbuf to cbuf based on imaptype */
        position = 0;
        MPI_Pack(lbuf, 1, imaptype, cbuf, (int)ibuf_size, &position,
                 MPI_COMM_SELF);
        MPI_Type_free(&imaptype);

        /* lbuf is no longer needed */
        if (lbuf != buf) NCI_Free(lbuf);
    }
    else /* not a true varm call: reuse lbuf */
        cbuf = lbuf;

    /* Step 3: type-convert and byte-swap cbuf to xbuf, and xbuf will be
     * used in MPI write function to write to file
     */
    if (need_convert) {
        /* user buf type does not match nc var type defined in file */
        void *fillp; /* fill value in internal representation */

        /* find the fill value */
        fillp = NCI_Malloc((size_t)varp->xsz);
        ncmpio_inq_var_fill(varp, fillp);

        /* datatype conversion + byte-swap from cbuf to xbuf */
        switch(varp->xtype) {
            case NC_BYTE:
                err = ncmpii_putn_NC_BYTE(fmt,xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_UBYTE:
                err = ncmpii_putn_NC_UBYTE(xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_SHORT:
                err = ncmpii_putn_NC_SHORT(xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_USHORT:
                err = ncmpii_putn_NC_USHORT(xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_INT:
                err = ncmpii_putn_NC_INT(xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_UINT:
                err = ncmpii_putn_NC_UINT(xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_FLOAT:
                err = ncmpii_putn_NC_FLOAT(xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_DOUBLE:
                err = ncmpii_putn_NC_DOUBLE(xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_INT64:
                err = ncmpii_putn_NC_INT64(xbuf,cbuf,nelems,etype,fillp);
                break;
            case NC_UINT64:
                err = ncmpii_putn_NC_UINT64(xbuf,cbuf,nelems,etype,fillp);
                break;
            default:
                err = NC_EBADTYPE; /* this never happens */
                break;
        }
        /* The only error codes returned from the above switch block are
	 * NC_EBADTYPE or NC_ERANGE. Bad varp->xtype and itype have been sanity
	 * checked at the dispatchers, so NC_EBADTYPE is not possible. Thus,
	 * the only possible error is NC_ERANGE.  NC_ERANGE can be caused by
	 * one or more elements of buf that is out of range representable by
	 * the external data type, it is not considered a fatal error. This
	 * request must continue to finish.
         */

        NCI_Free(fillp);
        if (cbuf != buf) NCI_Free(cbuf);
    }
    else {
        if (cbuf == buf && xbuf != buf)
            memcpy(xbuf, cbuf, (size_t)xbuf_size);

        if (need_swap) /* perform array in-place byte swap on xbuf */
            ncmpii_in_swapn(xbuf, nelems, varp->xsz);
    }
    return err;
}

/*----< ncmpio_unpack_xbuf() >-----------------------------------------------*/
/* Unpack xbuf into user buffer, buf, when type-casting is needed, imap is
 * non-contiguous, or buftype is non-contiguous. The immediate buffers, cbuf
 * and lbuf, may be allocated and freed within this subroutine. We try to reuse
 * the intermediate buffers as much as possible. Below describe such design.
 *
 * When called from get/iget APIs:
 *  if no convert && imap contig    && buftype contig
 *        xbuf  ==   cbuf   ==       lbuf    ==      buf
 *  if no convert && imap contig    && buftype noncontig
 *        xbuf  ==   cbuf   ==       lbuf  unpack->  buf
 *        malloc
 *  if no convert && imap noncontig && buftype contig
 *        xbuf  ==   cbuf  unpack->  lbuf    ==      buf
 *        malloc
 *  if no convert && imap noncontig && buftype noncontig
 *        xbuf  ==   cbuf  unpack->  lbuf  unpack->  buf
 *        malloc                     malloc
 *  if    convert && imap contig    && buftype contig
 *        xbuf  convert->  cbuf   ==       lbuf    ==      buf
 *        malloc
 *  if    convert && imap contig    && buftype noncontig
 *        xbuf  convert->  cbuf   ==       lbuf  unpack->  buf
 *        malloc           malloc
 *  if    convert && imap noncontig && buftype contig
 *        xbuf  convert->  cbuf  unpack->  lbuf    ==      buf
 *        malloc           malloc
 *  if    convert && imap noncontig && buftype noncontig
 *        xbuf  convert->  cbuf  unpack->  lbuf  unpack->  buf
 *        malloc           malloc          malloc
 */
int
ncmpio_unpack_xbuf(int           fmt,   /* NC_FORMAT_CDF2 NC_FORMAT_CDF5 etc. */
                   NC_var       *varp,
                   MPI_Offset    bufcount,
                   MPI_Datatype  buftype,
                   int           buftype_is_contig,
                   MPI_Offset    nelems, /* no. elements in etype in buf */
                   MPI_Datatype  etype,  /* element type in buftype */
                   MPI_Datatype  imaptype,
                   int           need_convert,
                   int           need_swap,
                   void         *buf,  /* user buffer */
                   void         *xbuf) /* already allocated, in external type */
{
    int err=NC_NOERR, el_size, position;
    void *lbuf=NULL, *cbuf=NULL;
    MPI_Offset ibuf_size;

    /* check byte size of buf (internal representation) */
    MPI_Type_size(etype, &el_size);
    ibuf_size = nelems * el_size;
    if (ibuf_size > INT_MAX) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    /* Step 1: type-convert and byte-swap xbuf to cbuf, and xbuf contains data
     * read from file
     */
    if (need_convert) {
        /* user buf type does not match nc var type defined in file */

        if (buftype_is_contig && imaptype == MPI_DATATYPE_NULL)
            /* both imap and buftype are contiguous */
            cbuf = buf;
        else {
            cbuf = NCI_Malloc(ibuf_size);
            if (cbuf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        }

        /* datatype conversion + byte-swap from xbuf to cbuf */
        switch(varp->xtype) {
            case NC_BYTE:
                err = ncmpii_getn_NC_BYTE(fmt,xbuf,cbuf,nelems,etype);
                break;
            case NC_UBYTE:
                err = ncmpii_getn_NC_UBYTE(xbuf,cbuf,nelems,etype);
                break;
            case NC_SHORT:
                err = ncmpii_getn_NC_SHORT(xbuf,cbuf,nelems,etype);
                break;
            case NC_USHORT:
                err = ncmpii_getn_NC_USHORT(xbuf,cbuf,nelems,etype);
                break;
            case NC_INT:
                err = ncmpii_getn_NC_INT(xbuf,cbuf,nelems,etype);
                break;
            case NC_UINT:
                err = ncmpii_getn_NC_UINT(xbuf,cbuf,nelems,etype);
                break;
            case NC_FLOAT:
                err = ncmpii_getn_NC_FLOAT(xbuf,cbuf,nelems,etype);
                break;
            case NC_DOUBLE:
                err = ncmpii_getn_NC_DOUBLE(xbuf,cbuf,nelems,etype);
                break;
            case NC_INT64:
                err = ncmpii_getn_NC_INT64(xbuf,cbuf,nelems,etype);
                break;
            case NC_UINT64:
                err = ncmpii_getn_NC_UINT64(xbuf,cbuf,nelems,etype);
                break;
            default:
                err = NC_EBADTYPE; /* this never happens */
                break;
        }
        /* The only error codes returned from the above switch block are
	 * NC_EBADTYPE or NC_ERANGE. Bad varp->xtype and itype have been sanity
	 * checked at the dispatchers, so NC_EBADTYPE is not possible. Thus,
	 * the only possible error is NC_ERANGE.  NC_ERANGE can be caused by
	 * one or more elements of buf that is out of range representable by
	 * the external data type, it is not considered a fatal error. This
	 * request must continue to finish.
         */
    }
    else {
        if (need_swap) /* perform array in-place byte swap on xbuf */
            ncmpii_in_swapn(xbuf, nelems, varp->xsz);
        cbuf = xbuf;
    }

    /* Step 2: if imap is non-contiguous, unpack cbuf to lbuf */
    /* determine whether we can use cbuf as lbuf */
    if (imaptype != MPI_DATATYPE_NULL && !buftype_is_contig) {
        /* a true varm and buftype is not contiguous: we need a separate
         * buffer, lbuf, to unpack cbuf to lbuf using imaptype, and later
         * unpack lbuf to buf using buftype.
         * In this case, cbuf cannot be buf and lbuf cannot be buf.
         */
        lbuf = NCI_Malloc((size_t)ibuf_size);
        if (lbuf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    else if (imaptype == MPI_DATATYPE_NULL) /* not varm */
        lbuf = cbuf;
    else /* varm and buftype are contiguous */
        lbuf = buf;

    /* unpacked cbuf into lbuf based on imap -------------------------------*/
    if (imaptype != MPI_DATATYPE_NULL) {
        /* unpack cbuf to lbuf based on imaptype */
        position = 0;
        MPI_Unpack(cbuf, (int)ibuf_size, &position, lbuf, 1, imaptype,
                   MPI_COMM_SELF);
        MPI_Type_free(&imaptype);
        /* done with cbuf */
        if (cbuf != xbuf) NCI_Free(cbuf);
    }

    /* unpacked lbuf into buf based on buftype -----------------------------*/
    if (!buftype_is_contig) {
        if (bufcount > INT_MAX) {
            if (err == NC_NOERR)
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
        }
        else {
            position = 0;
            MPI_Unpack(lbuf, (int)ibuf_size, &position, buf, (int)bufcount,
                       buftype, MPI_COMM_SELF);
            /* done with lbuf */
            if (lbuf != buf && lbuf != xbuf) NCI_Free(lbuf);
        }
    }

    return err;
}

