/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>   /* strtoll() is first introducted in C99 */
#include <string.h>   /* strcpy() */
#include <strings.h>  /* strcasecmp() */
#include <limits.h>   /* INT_MAX */
#include <assert.h>
#include <errno.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< ncmpio_set_pnetcdf_hints() >-----------------------------------------*/
/* this is where the I/O hints designated to pnetcdf are extracted and their
 * default values are set.
 */
void ncmpio_set_pnetcdf_hints(NC *ncp,
                              MPI_Info user_info,
                              MPI_Info info_used)
{
    char value[MPI_MAX_INFO_VAL];
    int  flag;

    if (user_info == MPI_INFO_NULL) flag = 0;

    /* Note info_used cannot be MPI_INFO_NULL, as it is returned from a call to
     * MPI_File_get_info()
     */
    assert(info_used != MPI_INFO_NULL);

    /* nc_header_align_size, nc_var_align_size, and r_align take effect when
     * a file is created, or opened and later adding more metadata or variable
     * data */

    if (user_info != MPI_INFO_NULL) {
        /* aligns the size of header extent of a newly created file */
        MPI_Info_get(user_info, "nc_header_align_size", MPI_MAX_INFO_VAL-1, value,
                     &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            ncp->h_align = strtoll(value, NULL, 10);
            if (errno != 0) ncp->h_align = 0;
            else if (ncp->h_align < 0) ncp->h_align = 0;
        }
    }
    if (!flag) sprintf(value, "%d", FILE_ALIGNMENT_DEFAULT);
    MPI_Info_set(info_used, "nc_header_align_size", value);

    if (user_info != MPI_INFO_NULL) {
        /* aligns starting file offsets of individual fixed-size variables */
        MPI_Info_get(user_info, "nc_var_align_size", MPI_MAX_INFO_VAL-1, value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            ncp->fx_v_align = strtoll(value, NULL, 10);
            if (errno != 0) ncp->fx_v_align = 0;
            else if (ncp->fx_v_align < 0) ncp->fx_v_align = 0;
        }
    }
    if (!flag) sprintf(value, "%d", FILE_ALIGNMENT_DEFAULT);
    MPI_Info_set(info_used, "nc_var_align_size", value);

    if (user_info != MPI_INFO_NULL) {
        /* aligns starting file offset of the record variable section */
        MPI_Info_get(user_info, "nc_record_align_size", MPI_MAX_INFO_VAL-1, value,
                     &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            ncp->r_align = strtoll(value, NULL, 10);
            if (errno != 0) ncp->r_align = 0;
            else if (ncp->r_align < 0) ncp->r_align = 0;
        }
    }
    if (!flag) sprintf(value, "%d", FILE_ALIGNMENT_DEFAULT);
    MPI_Info_set(info_used, "nc_record_align_size", value);

    if (user_info != MPI_INFO_NULL) {
        /* header reading chunk size */
        MPI_Info_get(user_info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1, value,
                     &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            ncp->chunk = (int) strtol(value, NULL, 10);
            if (errno != 0) ncp->chunk = 0;
            else if (ncp->chunk < 0) ncp->chunk = 0;
        }
    }
    if (!flag) sprintf(value, "%d", NC_DEFAULT_CHUNKSIZE);
    MPI_Info_set(info_used, "nc_header_read_chunk_size", value);

    if (user_info != MPI_INFO_NULL) {
        /* setting in-place byte swap (matters only for Little Endian) */
        MPI_Info_get(user_info, "nc_in_place_swap", MPI_MAX_INFO_VAL-1, value, &flag);
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
    }
    if (!flag) strcpy(value, "auto");
    MPI_Info_set(info_used, "nc_in_place_swap", value);

    if (user_info != MPI_INFO_NULL) {
	/* temporal buffer size used to pack noncontiguous aggregated user
         * buffers when calling ncmpi_wait/wait_all, Default 16 MiB
         */
        MPI_Info_get(user_info, "nc_ibuf_size", MPI_MAX_INFO_VAL-1, value, &flag);
        if (flag) {
            MPI_Offset ibuf_size;
            errno = 0;  /* errno must set to zero before calling strtoll */
            ibuf_size = strtoll(value, NULL, 10);
            if (errno == 0 && ncp->ibuf_size > 0) ncp->ibuf_size = ibuf_size;
        }
    }
    if (!flag) sprintf(value, "%d", NC_DEFAULT_IBUF_SIZE);
    MPI_Info_set(info_used, "nc_ibuf_size", value);

#ifdef ENABLE_SUBFILING
    if (user_info != MPI_INFO_NULL) {
        MPI_Info_get(user_info, "pnetcdf_subfiling", MPI_MAX_INFO_VAL-1, value, &flag);
        if (flag) {
            if (strcasecmp(value, "enable") == 0)
                ncp->subfile_mode = 1;
        }
    }
    if (!flag) strcpy(value, "disable");
    MPI_Info_set(info_used, "pnetcdf_subfiling", value);

    if (user_info != MPI_INFO_NULL) {
        MPI_Info_get(user_info, "nc_num_subfiles", MPI_MAX_INFO_VAL-1, value, &flag);
        if (flag) {
            errno = 0;
            ncp->num_subfiles = strtoll(value, NULL, 10);
            if (errno != 0) ncp->num_subfiles = 0;
            else if (ncp->num_subfiles < 0) ncp->num_subfiles = 0;
        }
    }
    if (!flag) strcpy(value, "0");
    MPI_Info_set(info_used, "nc_num_subfiles", value);

    if (ncp->subfile_mode == 0) ncp->num_subfiles = 0;
#else
    MPI_Info_set(info_used, "pnetcdf_subfiling", "disable");
    MPI_Info_set(info_used, "nc_num_subfiles", "0");
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
                   const MPI_Offset  start[],    /* [varp->ndims] */
                   const MPI_Offset  count[],    /* [varp->ndims] */
                   const MPI_Offset  stride[],   /* [varp->ndims] */
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

    /* Note check NC_EINVALCOORDS and NC_EEDGE already done in
     * dispatchers/var_getput.m4 */

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
                 MPI_Offset    nelems, /* no. elements in buf */
                 MPI_Datatype  itype,  /* element type in buftype */
                 int           el_size,/* size of itype */
                 MPI_Datatype  imaptype,
                 int           need_convert,
                 int           need_swap,
                 size_t        xbuf_size,
                 void         *buf,    /* user buffer */
                 void         *xbuf)   /* already allocated, in external type */
{
    int err=NC_NOERR, position, free_lbuf=0, free_cbuf=0;
    void *lbuf=NULL, *cbuf=NULL;
    MPI_Offset ibuf_size;

    /* check byte size of buf (internal representation) */
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
            free_lbuf = 1;
        }

        if (buf != lbuf) {
            /* pack buf into lbuf based on buftype */
            if (bufcount > INT_MAX) {
                if (free_lbuf) NCI_Free(lbuf);
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            }
            position = 0;
            MPI_Pack(buf, (int)bufcount, buftype, lbuf, (int)ibuf_size,
                     &position, MPI_COMM_SELF);
        }
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
            if (cbuf == NULL) {
                if (free_lbuf) NCI_Free(lbuf);
                DEBUG_RETURN_ERROR(NC_ENOMEM)
            }
            free_cbuf = 1;
        }

        /* pack lbuf to cbuf based on imaptype */
        position = 0;
        MPI_Pack(lbuf, 1, imaptype, cbuf, (int)ibuf_size, &position,
                 MPI_COMM_SELF);
        MPI_Type_free(&imaptype);

        /* lbuf is no longer needed */
        if (free_lbuf) {
            NCI_Free(lbuf);
            free_lbuf = 0;
        }
    }
    else /* not a true varm call: reuse lbuf */
        cbuf = lbuf;

    /* Step 3: type-convert and byte-swap cbuf to xbuf, and xbuf will be
     * used in MPI write function to write to file
     */
    if (need_convert) {
        /* user buf type does not match nc var type defined in file */
        char  tmpbuf[8];
        void *fillp=tmpbuf; /* fill value in internal representation */

        /* find the fill value */
        ncmpio_inq_var_fill(varp, fillp);

        /* datatype conversion + byte-swap from cbuf to xbuf */
        switch(varp->xtype) {
            case NC_BYTE:
                err = ncmpii_putn_NC_BYTE(fmt,xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_UBYTE:
                err = ncmpii_putn_NC_UBYTE(xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_SHORT:
                err = ncmpii_putn_NC_SHORT(xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_USHORT:
                err = ncmpii_putn_NC_USHORT(xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_INT:
                err = ncmpii_putn_NC_INT(xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_UINT:
                err = ncmpii_putn_NC_UINT(xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_FLOAT:
                err = ncmpii_putn_NC_FLOAT(xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_DOUBLE:
                err = ncmpii_putn_NC_DOUBLE(xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_INT64:
                err = ncmpii_putn_NC_INT64(xbuf,cbuf,nelems,itype,fillp);
                break;
            case NC_UINT64:
                err = ncmpii_putn_NC_UINT64(xbuf,cbuf,nelems,itype,fillp);
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
        if (free_cbuf) NCI_Free(cbuf);
        if (free_lbuf) NCI_Free(lbuf);
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
                   MPI_Offset    nelems, /* no. elements in buf */
                   MPI_Datatype  itype,  /* element type in buftype */
                   MPI_Datatype  imaptype,
                   int           need_convert,
                   int           need_swap,
                   void         *buf,  /* user buffer */
                   void         *xbuf) /* already allocated, in external type */
{
    int err=NC_NOERR, el_size, position, free_lbuf=0, free_cbuf=0;
    void *lbuf=NULL, *cbuf=NULL;
    MPI_Offset ibuf_size;

    /* check byte size of buf (internal representation) */
    MPI_Type_size(itype, &el_size);
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
            free_cbuf = 1;
        }

        /* datatype conversion + byte-swap from xbuf to cbuf */
        switch(varp->xtype) {
            case NC_BYTE:
                err = ncmpii_getn_NC_BYTE(fmt,xbuf,cbuf,nelems,itype);
                break;
            case NC_UBYTE:
                err = ncmpii_getn_NC_UBYTE(xbuf,cbuf,nelems,itype);
                break;
            case NC_SHORT:
                err = ncmpii_getn_NC_SHORT(xbuf,cbuf,nelems,itype);
                break;
            case NC_USHORT:
                err = ncmpii_getn_NC_USHORT(xbuf,cbuf,nelems,itype);
                break;
            case NC_INT:
                err = ncmpii_getn_NC_INT(xbuf,cbuf,nelems,itype);
                break;
            case NC_UINT:
                err = ncmpii_getn_NC_UINT(xbuf,cbuf,nelems,itype);
                break;
            case NC_FLOAT:
                err = ncmpii_getn_NC_FLOAT(xbuf,cbuf,nelems,itype);
                break;
            case NC_DOUBLE:
                err = ncmpii_getn_NC_DOUBLE(xbuf,cbuf,nelems,itype);
                break;
            case NC_INT64:
                err = ncmpii_getn_NC_INT64(xbuf,cbuf,nelems,itype);
                break;
            case NC_UINT64:
                err = ncmpii_getn_NC_UINT64(xbuf,cbuf,nelems,itype);
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
        if (lbuf == NULL) {
            if (free_cbuf) NCI_Free(cbuf);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
        free_lbuf = 1;
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
    }

    /* unpacked lbuf into buf based on buftype -----------------------------*/
    if (!buftype_is_contig && lbuf != buf) {
        /* no need unpack when buftype is used in MPI_File_read (lbuf == buf) */
        if (bufcount > INT_MAX) {
            if (err == NC_NOERR)
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
        }
        else {
            position = 0;
            MPI_Unpack(lbuf, (int)ibuf_size, &position, buf, (int)bufcount,
                       buftype, MPI_COMM_SELF);
        }
    }
    if (free_cbuf) NCI_Free(cbuf);
    if (free_lbuf) NCI_Free(lbuf);

    return err;
}

