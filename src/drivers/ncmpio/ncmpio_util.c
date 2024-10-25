/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>   /* strtoll() is first introduced in C99 */
#include <string.h>   /* strcpy() strdup() */
#include <strings.h>  /* strcasecmp() */
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

    ncp->env_h_align = 0;
    if (user_info != MPI_INFO_NULL) {
        /* aligns the size of header extent of a newly created file */
        MPI_Info_get(user_info, "nc_header_align_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            ncp->env_h_align = strtoll(value, NULL, 10);
            if (errno != 0) ncp->env_h_align = 0;
            else if (ncp->env_h_align < 0) ncp->env_h_align = 0;
        }
    }
    if (ncp->env_h_align == 0)
        sprintf(value, "%d", FILE_ALIGNMENT_DEFAULT);
    else
        sprintf(value, "%lld", ncp->env_h_align);
    MPI_Info_set(info_used, "nc_header_align_size", value);

    ncp->env_v_align = 0;
    if (user_info != MPI_INFO_NULL) {
        /* aligns starting file offsets of individual fixed-size variables */
        MPI_Info_get(user_info, "nc_var_align_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            ncp->env_v_align = strtoll(value, NULL, 10);
            if (errno != 0) ncp->env_v_align = 0;
            else if (ncp->env_v_align < 0) ncp->env_v_align = 0;
        }
    }
    if (ncp->env_v_align == 0)
        sprintf(value, "%d", FILE_ALIGNMENT_DEFAULT);
    else
        sprintf(value, "%lld", ncp->env_v_align);
    MPI_Info_set(info_used, "nc_var_align_size", value);

    ncp->env_r_align = 0;
    if (user_info != MPI_INFO_NULL) {
        /* aligns starting file offset of the record variable section */
        MPI_Info_get(user_info, "nc_record_align_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            ncp->env_r_align = strtoll(value, NULL, 10);
            if (errno != 0) ncp->env_r_align = 0;
            else if (ncp->env_r_align < 0) ncp->env_r_align = 0;
        }
    }
    if (ncp->env_r_align == 0)
        sprintf(value, "%d", FILE_ALIGNMENT_DEFAULT);
    else
        sprintf(value, "%lld", ncp->env_r_align);
    MPI_Info_set(info_used, "nc_record_align_size", value);

    ncp->chunk = PNC_DEFAULT_CHUNKSIZE;
    if (user_info != MPI_INFO_NULL) {
        /* header reading chunk size */
        MPI_Info_get(user_info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            int chunk;
            errno = 0;  /* errno must set to zero before calling strtoll */
            chunk = atoi(value);
            if (errno != 0) ncp->chunk = 0;
            else if (ncp->chunk < 0)
                ncp->chunk = 0;
            else if (chunk > NC_MAX_INT) /* limit to NC_MAX_INT */
                ncp->chunk = NC_MAX_INT;
        }
    }
    sprintf(value, "%d", ncp->chunk);
    MPI_Info_set(info_used, "nc_header_read_chunk_size", value);

    strcpy(value, "auto");
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
    MPI_Info_set(info_used, "nc_in_place_swap", value);

    ncp->ibuf_size = PNC_DEFAULT_IBUF_SIZE;
    if (user_info != MPI_INFO_NULL) {
	/* temporal buffer size used to pack noncontiguous aggregated user
         * buffers when calling ncmpi_wait/wait_all, Default 16 MiB
         */
        MPI_Info_get(user_info, "nc_ibuf_size", MPI_MAX_INFO_VAL-1, value,
                     &flag);
        if (flag) {
            MPI_Offset ibuf_size;
            errno = 0;  /* errno must set to zero before calling strtoll */
            ibuf_size = strtoll(value, NULL, 10);
            if (errno == 0 && ncp->ibuf_size > 0) ncp->ibuf_size = ibuf_size;
        }
    }
    sprintf(value, "%lld", ncp->ibuf_size);
    MPI_Info_set(info_used, "nc_ibuf_size", value);

#ifdef ENABLE_SUBFILING
    ncp->subfile_mode = 0;
    if (user_info != MPI_INFO_NULL) {
        MPI_Info_get(user_info, "pnetcdf_subfiling", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            if (strcasecmp(value, "enable") == 0)
                ncp->subfile_mode = 1;
        }
    }
    if (ncp->subfile_mode)
        MPI_Info_set(info_used, "pnetcdf_subfiling", "enable");
    else
        MPI_Info_set(info_used, "pnetcdf_subfiling", "disable");

    ncp->num_subfiles = 0;
    if (user_info != MPI_INFO_NULL) {
        MPI_Info_get(user_info, "nc_num_subfiles", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;
            ncp->num_subfiles = atoi(value);
            if (errno != 0) ncp->num_subfiles = 0;
            else if (ncp->num_subfiles < 0) ncp->num_subfiles = 0;
        }
    }
    sprintf(value, "%d", ncp->num_subfiles);
    MPI_Info_set(info_used, "nc_num_subfiles", value);

    if (ncp->subfile_mode == 0) ncp->num_subfiles = 0;
#else
    MPI_Info_set(info_used, "pnetcdf_subfiling", "disable");
    MPI_Info_set(info_used, "nc_num_subfiles", "0");
#endif

    if (user_info != MPI_INFO_NULL) {
        /* If romio_no_indep_rw is set to true, let all processes participate
         * the read/write file header using MPI collective APIs, where only
         * rank 0 has non-zero request count.
         */
        MPI_Info_get(user_info, "romio_no_indep_rw", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            if (strcasecmp(value, "true") == 0)
                fSet((ncp)->flags, NC_HCOLL);
        }
    }

    ncp->dims.hash_size = PNC_HSIZE_DIM;
    if (user_info != MPI_INFO_NULL) {
        /* Hash table size for dimensions */
        MPI_Info_get(user_info, "nc_hash_size_dim", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling atoi */
            ncp->dims.hash_size = atoi(value);
            if (errno != 0 || ncp->dims.hash_size < 0)
                ncp->dims.hash_size = PNC_HSIZE_DIM;
        }
    }
    sprintf(value, "%d", ncp->dims.hash_size);
    MPI_Info_set(info_used, "nc_hash_size_dim", value);

    ncp->vars.hash_size = PNC_HSIZE_VAR;
    if (user_info != MPI_INFO_NULL) {
        /* Hash table size for variables */
        MPI_Info_get(user_info, "nc_hash_size_var", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling atoi */
            ncp->vars.hash_size = atoi(value);
            if (errno != 0 || ncp->vars.hash_size < 0)
                ncp->vars.hash_size = PNC_HSIZE_VAR;
        }
    }
    sprintf(value, "%d", ncp->vars.hash_size);
    MPI_Info_set(info_used, "nc_hash_size_var", value);

    ncp->attrs.hash_size = PNC_HSIZE_GATTR;
    if (user_info != MPI_INFO_NULL) {
        /* Hash table size for global attributes */
        MPI_Info_get(user_info, "nc_hash_size_gattr", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling atoi */
            ncp->attrs.hash_size = atoi(value);
            if (errno != 0 || ncp->attrs.hash_size < 0)
                ncp->attrs.hash_size = PNC_HSIZE_GATTR;
        }
    }
    sprintf(value, "%d", ncp->attrs.hash_size);
    MPI_Info_set(info_used, "nc_hash_size_gattr", value);

    ncp->hash_size_attr = PNC_HSIZE_VATTR;
    if (user_info != MPI_INFO_NULL) {
        /* Hash table size for non-global attributes */
        MPI_Info_get(user_info, "nc_hash_size_vattr", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling atoi */
            ncp->hash_size_attr = atoi(value);
            if (errno != 0 || ncp->hash_size_attr < 0)
                ncp->hash_size_attr = PNC_HSIZE_VATTR;
        }
    }
    sprintf(value, "%d", ncp->hash_size_attr);
    MPI_Info_set(info_used, "nc_hash_size_vattr", value);

    ncp->num_aggrs_per_node = 0;
    if (user_info != MPI_INFO_NULL) {
        /* Hash table size for non-global attributes */
        MPI_Info_get(user_info, "nc_num_aggrs_per_node", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling atoi */
            ncp->num_aggrs_per_node = atoi(value);
            if (errno != 0 || ncp->num_aggrs_per_node < 0)
                ncp->num_aggrs_per_node = 0;
        }
    }
    sprintf(value, "%d", ncp->num_aggrs_per_node);
    MPI_Info_set(info_used, "nc_num_aggrs_per_node", value);
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
    int err=NC_NOERR, mpireturn, free_lbuf=0, free_cbuf=0;
    void *lbuf=NULL, *cbuf=NULL;
    MPI_Offset ibuf_size;

    /* check byte size of buf (internal representation) */
    ibuf_size = nelems * el_size;

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
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count position = 0;
            mpireturn = MPI_Pack_c(buf, (MPI_Count)bufcount, buftype, lbuf,
                                   (MPI_Count)ibuf_size, &position, MPI_COMM_SELF);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_error_mpi2nc(mpireturn, "MPI_Pack_c");
#else
            int position = 0;
            if (bufcount > NC_MAX_INT || ibuf_size > NC_MAX_INT) {
                if (free_lbuf) NCI_Free(lbuf);
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            }
            mpireturn = MPI_Pack(buf, (int)bufcount, buftype, lbuf, (int)ibuf_size,
                                 &position, MPI_COMM_SELF);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_error_mpi2nc(mpireturn, "MPI_Pack");
#endif
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
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count position = 0;
        mpireturn = MPI_Pack_c(lbuf, 1, imaptype, cbuf, (MPI_Count)ibuf_size,
                               &position, MPI_COMM_SELF);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Pack_c");
#else
        int position = 0;
        if (ibuf_size > NC_MAX_INT)
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        mpireturn = MPI_Pack(lbuf, 1, imaptype, cbuf, (int)ibuf_size,
                             &position, MPI_COMM_SELF);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Pack");
#endif
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
                DEBUG_ASSIGN_ERROR(err, NC_EBADTYPE) /* this never happens */
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
    int err=NC_NOERR, mpireturn, el_size, free_lbuf=0, free_cbuf=0;
    void *lbuf=NULL, *cbuf=NULL;
    MPI_Offset ibuf_size;

    /* check byte size of buf (internal representation) */
    MPI_Type_size(itype, &el_size); /* itype is MPI primitive datatype */
    ibuf_size = nelems * el_size;

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
                DEBUG_ASSIGN_ERROR(err, NC_EBADTYPE) /* this never happens */
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
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count position = 0;
        mpireturn = MPI_Unpack_c(cbuf, (MPI_Count)ibuf_size, &position, lbuf,
                                 1, imaptype, MPI_COMM_SELF);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Unpack_c");
#else
        int position = 0;
        if (ibuf_size > NC_MAX_INT) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

        mpireturn = MPI_Unpack(cbuf, (int)ibuf_size, &position, lbuf,
                               1, imaptype, MPI_COMM_SELF);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Unpack");
#endif
        MPI_Type_free(&imaptype);
    }

    /* unpacked lbuf into buf based on buftype -----------------------------*/
    if (!buftype_is_contig && lbuf != buf) {
        /* no need unpack when buftype is used in MPI_File_read (lbuf == buf) */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count position = 0;
        mpireturn = MPI_Unpack_c(lbuf, (MPI_Count)ibuf_size, &position, buf,
                     (MPI_Count)bufcount, buftype, MPI_COMM_SELF);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Unpack_c");
#else
        if (bufcount > NC_MAX_INT) {
            if (err == NC_NOERR)
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
        }
        else {
            int position = 0;
            if (ibuf_size > NC_MAX_INT)
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            mpireturn = MPI_Unpack(lbuf, (int)ibuf_size, &position, buf,
                                   (int)bufcount, buftype, MPI_COMM_SELF);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_error_mpi2nc(mpireturn, "MPI_Unpack");
        }
#endif
    }
    if (free_cbuf) NCI_Free(cbuf);
    if (free_lbuf) NCI_Free(lbuf);

    return err;
}

/*----< ncmpio_construct_aggr_list() >---------------------------------------*/
/* Find information about MPI processes and their affinity to compute node.
 * Determine whether self process is a intra-node aggregator.
 * For an aggregator, construct rank IDs of assigned non-aggregators.
 * For a non-aggregator, calculate the rank ID of its assigned aggregator.
 */
int ncmpio_construct_aggr_list(NC *ncp)
{
    char my_procname[MPI_MAX_PROCESSOR_NAME], **all_procnames=NULL;
    int i, j, k, rank, nprocs, my_procname_len, num_nodes, root=0;
    int *node_ids=NULL, *all_procname_lens=NULL, *nprocs_per_node;
    int avg_nprocs_per_node;

    MPI_Comm_size(ncp->comm, &nprocs);
    MPI_Comm_rank(ncp->comm, &rank);

    /* Collect info about compute nodes in order to select I/O aggregators.
     * Note my_procname is null character terminated, but my_procname_len
     * does not include the null character.
     */
    MPI_Get_processor_name(my_procname, &my_procname_len);
    my_procname_len++; /* to include terminate null character */

#ifdef TEST_ALONE
printf("%d: ---------  my_procname=%s\n",rank,my_procname);
/*
if (rank < (nprocs/2)+(nprocs%2)) {sprintf(my_procname,"dummy.0"); my_procname_len=strlen(my_procname)+1;}
else {sprintf(my_procname,"dummy.1"); my_procname_len=strlen(my_procname)+1;}
*/
#endif

    if (rank == root) {
        /* root collects all procnames */
        all_procnames = (char **) NCI_Malloc(sizeof(char*) * nprocs);
        if (all_procnames == NULL)
            DEBUG_RETURN_ERROR(NC_ENOMEM)

        all_procname_lens = (int *) NCI_Malloc(sizeof(int) * nprocs);
        if (all_procname_lens == NULL) {
            NCI_Free(all_procnames);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
    }
    /* gather process name lengths from all processes first */
    MPI_Gather(&my_procname_len, 1, MPI_INT, all_procname_lens, 1, MPI_INT,
               root, ncp->comm);

    if (rank == root) {
        int *disp;
        size_t alloc_size = 0;

        for (i=0; i<nprocs; i++)
            alloc_size += all_procname_lens[i];

        all_procnames[0] = (char *) NCI_Malloc(alloc_size);
        if (all_procnames[0] == NULL) {
            NCI_Free(all_procname_lens);
            NCI_Free(all_procnames);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }

        /* Construct displacement array for the MPI_Gatherv, as each process
         * may have a different length for its process name.
         */
        disp = (int *) NCI_Malloc(sizeof(int) * nprocs);
        disp[0] = 0;
        for (i=1; i<nprocs; i++) {
            all_procnames[i] = all_procnames[i - 1] + all_procname_lens[i - 1];
            disp[i] = disp[i - 1] + all_procname_lens[i - 1];
        }

        /* gather all process names */
        MPI_Gatherv(my_procname, my_procname_len, MPI_CHAR,
                    all_procnames[0], all_procname_lens, disp, MPI_CHAR,
                    root, ncp->comm);

        NCI_Free(disp);
        NCI_Free(all_procname_lens);
    } else
        /* send process name to root */
        MPI_Gatherv(my_procname, my_procname_len, MPI_CHAR,
                    NULL, NULL, NULL, MPI_CHAR, root, ncp->comm);

    /* each MPI process's compute node ID */
    node_ids = (int *) NCI_Malloc(sizeof(int) * (nprocs + 1));

    if (rank == root) {
        /* all_procnames[] can tell us the number of nodes and number of
         * processes per node.
         */
        char **node_names;
        int last;

        /* array of pointers pointing to unique host names (compute nodes) */
        node_names = (char **) NCI_Malloc(sizeof(char*) * nprocs);

        /* number of MPI processes running on each node */
        nprocs_per_node = (int *) NCI_Malloc(sizeof(int) * nprocs);

        /* calculate nprocs_per_node[] and node_ids[] */
        last = 0;
        num_nodes = 0; /* number of unique compute nodes */
        for (i=0; i<nprocs; i++) {
            k = last;
            for (j=0; j<num_nodes; j++) {
                /* check if [i] has already appeared in [] */
                if (!strcmp(all_procnames[i], node_names[k])) { /* found */
                    node_ids[i] = k;
                    break;
                }
                k = (k == num_nodes - 1) ? 0 : k + 1;
            }
            if (j < num_nodes)  /* found, next iteration, start with node n */
                last = k;
            else {      /* not found, j == num_nodes, add a new node */
                node_names[j] = strdup(all_procnames[i]);
                nprocs_per_node[j] = 1;
                node_ids[i] = j;
                last = j;
                num_nodes++;
            }
        }
        /* num_nodes is now the number of compute nodes (unique node names) */

        avg_nprocs_per_node = nprocs_per_node[0];
        for (i=1; i<num_nodes; i++) {
            avg_nprocs_per_node += nprocs_per_node[i];
        }
        avg_nprocs_per_node /= num_nodes;

        NCI_Free(nprocs_per_node);

        /* also bcast avg_nprocs_per_node */
        node_ids[nprocs] = avg_nprocs_per_node;

        for (i=0; i<num_nodes; i++)
            free(node_names[i]); /* allocated by strdup() */
        NCI_Free(node_names);
        NCI_Free(all_procnames[0]);
        NCI_Free(all_procnames);
    }

    MPI_Bcast(node_ids, nprocs+1, MPI_INT, root, ncp->comm);

    /* TODO: what is the threshold to enable intra-node aggregation?
     * 64 process per nodes?
     */
    avg_nprocs_per_node = node_ids[nprocs];
    // ncp->aggregation = (avg_nprocs_per_node > 64) ? 1 : 0;
    ncp->num_non_aggrs = 0;

    /* my_node_id is this rank's node ID, which can be used to calculate the
     * number of processes in the same ran and their rank IDs.
     */
    int my_node_id = node_ids[rank];

    /* nprocs_my_node: the number of processes in my nodes
     * ranks_my_node[]: rank IDs of all processes in my node and it can be
     *                  used to select aggregators, e.g. based on a
     *                  ratio of 1/8, i.e. every 8 ranks on the same node is
     *                  selected to be an aggregator
     */
    int my_rank_index;
    int *ranks_my_node = (int*) NCI_Malloc(sizeof(int) * nprocs);
    for (j=0, i=0; i<nprocs; i++) {
        if (node_ids[i] == my_node_id) {
            if (i == rank)
                my_rank_index = j;
            ranks_my_node[j++] = i;
        }
    }
    int nprocs_my_node = j;
    /* Now, ranks_my_node[my_rank_index] == rank */

#ifdef TEST_ALONE
if (nprocs ==2) printf("%d: node_ids=%d %d ranks_my_node=%d %d\n",rank,node_ids[0],node_ids[1],ranks_my_node[0],ranks_my_node[1]);
else if (nprocs ==5) printf("%d: node_ids=%d %d %d %d %d ranks_my_node=%d %d %d %d %d\n",rank,
node_ids[0],node_ids[1],node_ids[2],node_ids[3],node_ids[4],
ranks_my_node[0],ranks_my_node[1],ranks_my_node[2],ranks_my_node[3],ranks_my_node[4]);
printf("%d: my_node_id=%d nprocs_my_node=%d my_rank_index=%d\n",rank,my_node_id,nprocs_my_node,my_rank_index);
#endif

    NCI_Free(node_ids);

    /* each rank needs to know:
     *   the aggregator assigned to it (rank ID)
     * each aggregator need to know:
     *   the non-aggregators assigned to it (rank IDs)
     *   the number of non-aggregators assigned to it
     */

    int naggrs_my_node = MIN(ncp->num_aggrs_per_node, nprocs_my_node);
    int num_non_aggrs = nprocs_my_node / naggrs_my_node;
    if (nprocs_my_node % naggrs_my_node) num_non_aggrs++;
    if (num_non_aggrs == 1) //<= naggrs_my_node)
        ncp->aggregation = 0;
    else {
        ncp->my_aggr = ranks_my_node[my_rank_index - my_rank_index % num_non_aggrs];
        ncp->isAggr = (ncp->my_aggr == rank) ? 1 : 0;
        if (ncp->isAggr) {
            ncp->num_non_aggrs = MIN(num_non_aggrs, nprocs_my_node - my_rank_index);
            if (ncp->num_non_aggrs == 1)
                ncp->aggregation = 0;
            else
                memcpy(ncp->nonaggr_ranks,
                       ranks_my_node + my_rank_index,
                       sizeof(int) * num_non_aggrs);
        }
    }

    NCI_Free(ranks_my_node);

#ifdef TEST_ALONE
printf("%d: ncp->aggregation=%d nprocs_my_node=%d naggrs_my_node=%d my_aggr=%d isAggr=%d num_non_aggrs=%d\n", rank, ncp->aggregation, nprocs_my_node,naggrs_my_node, ncp->my_aggr, ncp->isAggr, ncp->num_non_aggrs);
#endif

#if 0
    NUM_AGGRS_PER_NODE = 16

    nprocs_my_node = 128
    naggrs_my_node = MIN(NUM_AGGRS_PER_NODE, nprocs_my_node); 16
    num_non_aggrs = nprocs_my_node / naggrs_my_node; 128 /16 = 8

    nprocs_my_node = 122
    naggrs_my_node = MIN(NUM_AGGRS_PER_NODE, nprocs_my_node); 16
    num_non_aggrs = ceil(nprocs_my_node / naggrs_my_node); = 8
    first 15 aggrs are assigned 8 nonaggregators
    last aggr is assigned 2 nonaggregators
#endif

    /* TODO: handle boundary condition. For example, the last node may have
     *       not enough number of processes
     */

    /* For automatically determine Whether to enable intra-node write aggregation,
     * this should be done right before each collective write call.
     * 1. obtain hint cb_noddes, and striping_unit
     * 2. calculate aggregate access region
     * If in each round of two-phase I/O,
     * the number of senders to each cb_nodes is very large, then intra-node
     * aggregation should be enabled. avg_nprocs_per_node can be a factor for
     * determining whether to enable intra-node aggregation.
     *
     * When intra-node write aggregation is enabled, processes on the same node
     * will be divided into groups. The number of groups is the number of
     * aggregators on that node. The rank IDs of each group must be established.
     */

    return NC_NOERR;
}

