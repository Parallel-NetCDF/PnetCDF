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
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncx.h>
#include "ncmpio_NC.h"


#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
/*----< check_recsize_too_big() >--------------------------------------------*/
/* Because recsize (sum of single record of all record variables) will be used
 * to as the "offset stride" when constructing the file view, we must ensure
 * there is no integer overflow. Note records of all variables are interleaved,
 * for example record i of all record variables followed by record i+1 of all
 * record variables, the stride across two records of a variable can be too
 * large for a 4-byte integer to represent.
 */
static int
check_recsize_too_big(MPI_Offset recsize)
{
    int ret = NC_NOERR;

    if (recsize != (MPI_Aint)recsize) {
        fprintf(stderr, "Type overflow: unable to read/write multiple records in this dataset\non this platform. Please either access records of this record variable\none-at-a-time or run on a 64 bit platform\n");
        DEBUG_ASSIGN_ERROR(ret, NC_ESMALL)
    }
    /* the assert here might be harsh, but without it, users will get corrupt
     * data. Now, we just skip this request to avoid this assertion. */
    /* assert (recsize == (MPI_Aint)recsize); */
    return ret;
}
#endif

/*----< is_request_contiguous() >-------------------------------------------*/
/* check if a request, represented by start[] and count[], is contiguous in
 * file.
 */
static int
is_request_contiguous(int               isRecVar,
                      int               numRecVars,
                      int               ndims,
                      const MPI_Offset *shape,
                      const MPI_Offset *start,
                      const MPI_Offset *count)
{
    int i, j, most_sig_dim;

    if (ndims == 0) return 1; /* this variable is a scalar */

    for (i=0; i<ndims; i++)
         if (count[i] == 0) /* zero length request */
             return 1;

    /* most-significant dim. record dimension */
    most_sig_dim = 0;

    if (isRecVar) {
        /* if there are more than one record variable and when count[0] > 1,
         * then this request is noncontiguous.
         */
        if (numRecVars > 1) {
            if (count[0] > 1) return 0;

            /* continue to check from dimension ndims-1 up to dimension 1 */
            most_sig_dim = 1;
        }
        /* if there is only one record variable, then we need to check
         * dimension ndims-1 up to dimension 0 */
    }

    for (i=ndims-1; i>most_sig_dim; i--) {
        /* find the first count[i] that is not equal to the entire dimension */
        if (count[i] < shape[i]) {
            /* the request is contiguous only if count[i-1], count[i-2], ...
             * count[0] are all 1s. */
            for (j=i-1; j>=most_sig_dim; j--) {
                if (count[j] > 1)
                    return 0;
            }
            break;
        }
#if 0
        else { /* count[i] == shape[i] */
            /* when accessing the entire dimension, start[i] must be 0 */
            if (start[i] != 0) return 0; /* actually this should be error */
        }
#endif
    }
    return 1;
}

/*----< type_create_subarray64() >-------------------------------------------*/
/* This subroutine is to achieve the same result as MPI_Type_create_subarray()
 * but it takes arguments in type of MPI_Offset, instead of int. It also
 * checked for any possible 4-byte integer overflow.
 */
static int
type_create_subarray64(int               ndims,
                       const MPI_Offset *array_of_sizes,    /* [ndims] */
                       const MPI_Offset *array_of_subsizes, /* [ndims] */
                       const MPI_Offset *array_of_starts,   /* [ndims] */
                       int               order,
                       MPI_Datatype      oldtype,
                       MPI_Datatype     *newtype)
{
    int i, err=NC_NOERR, mpireturn;

    if (ndims == 0) DEBUG_RETURN_ERROR(NC_EDIMMETA)

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *sizes, *subsizes, *starts;

    sizes = (MPI_Count*) NCI_Malloc((size_t)ndims * 3 * sizeof(MPI_Count));
    subsizes = sizes    + ndims;
    starts   = subsizes + ndims;
    for (i=0; i<ndims; i++) {
        sizes[i]    = (MPI_Count)array_of_sizes[i];
        subsizes[i] = (MPI_Count)array_of_subsizes[i];
        starts[i]   = (MPI_Count)array_of_starts[i];
    }
    mpireturn = MPI_Type_create_subarray_c(ndims, sizes, subsizes, starts,
                                           order, oldtype, newtype);
    if (mpireturn != MPI_SUCCESS)
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_subarray_c");

    NCI_Free(sizes);
    return err;

#else
    int big_int, blklens[3] = {1, 1, 1};
    MPI_Datatype type1, type2;
    MPI_Aint extent, size, array_size, stride, disps[3];
    MPI_Aint lb;

    /* check if any of the dimensions is larger than 2^31-1 */
    big_int = 0;
    for (i=0; i<ndims; i++) {
        if (array_of_sizes[i] > NC_MAX_INT || array_of_starts[i] > NC_MAX_INT) {
            big_int = 1;
            break;
        }
    }

    if (big_int == 0) {
        int *sizes, *subsizes, *starts;
        /* none of dimensions > 2^31-1, we can safely use
         * MPI_Type_create_subarray */
        sizes = (int*) NCI_Malloc((size_t)ndims * 3 * SIZEOF_INT);
        subsizes = sizes    + ndims;
        starts   = subsizes + ndims;
        for (i=0; i<ndims; i++) {
            sizes[i]    = (int)array_of_sizes[i];
            subsizes[i] = (int)array_of_subsizes[i];
            starts[i]   = (int)array_of_starts[i];
        }
        mpireturn = MPI_Type_create_subarray(ndims, sizes, subsizes, starts,
                                             order, oldtype, newtype);
        if (mpireturn != MPI_SUCCESS)
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_subarray");

        NCI_Free(sizes);
        return err;
    }

    /* now big_int = 1 */

    /* at least one dimension is of size > 2^31-1 and we cannot use
     * MPI_Type_create_subarray() to create the newtype,
     * as its arguments array_of_sizes[] and array_of_starts[] are of
     * type int. One solution is to use a combination of
     * MPI_Type_create_hvector(), MPI_Type_create_hindexed(),
     * MPI_Type_create_resized(), and MPI_Type_struct(), as one
     * of their arguments, stride and indices[], are of type MPI_Aint
     * (a possible 8-byte integer) that can be used to store the value
     * of the dimension whose size is > 2^31-1. However, on a machine
     * where MPI_Aint is 4-byte integer, those MPI_Aint arguments will
     * cause overflow.
     */
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
    DEBUG_RETURN_ERROR(NC_EAINT_TOO_SMALL)
#endif

    MPI_Type_get_extent(oldtype, &lb, &extent);
    array_size = extent;
    for (i=0; i<ndims; i++) array_size *= array_of_sizes[i];

    if (ndims == 1) {
        /* blklens argument in MPI_Type_create_hindexed() is of type int */
        if (array_of_subsizes[0] > NC_MAX_INT) /* check int overflow */
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        blklens[1] = (int)array_of_subsizes[0];
        disps[1] = extent * array_of_starts[0];

        /* take advantage of disps argument is of type MPI_Aint */
        err = MPI_Type_create_hindexed(1, &blklens[1], &disps[1], oldtype, &type1);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_hindexed");
        MPI_Type_commit(&type1);

        /* add holes in the beginning and tail of type1 */
        err = MPI_Type_create_resized(type1, 0, array_size, newtype);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_resized");
        MPI_Type_free(&type1);

        return NC_NOERR;
    }
    /* now, ndims > 1 */

    /* first create a datatype for the least 2 significant dimensions */

    /* count and blocklength arguments in MPI_Type_create_hvector() are of
     * type int. We need to check for integer overflow */
    int count, blocklength;

    if (array_of_subsizes[ndims-2] > NC_MAX_INT ||
        array_of_subsizes[ndims-1] > NC_MAX_INT) /* check int overflow */
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    count = (int)array_of_subsizes[ndims-2];
    blocklength = (int)array_of_subsizes[ndims-1];
    stride = array_of_sizes[ndims-1] * extent;

    err = MPI_Type_create_hvector(count, blocklength, stride, oldtype, &type1);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");

    MPI_Type_commit(&type1);

    /* now iterate through the rest dimensions */
    for (i=ndims-3; i>=0; i--) {
        if (array_of_subsizes[i] > NC_MAX_INT) /* check int overflow */
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

        count = (int)array_of_subsizes[i];
        stride *= array_of_sizes[i+1];

        err = MPI_Type_create_hvector(count, 1, stride, type1, &type2);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");

        MPI_Type_commit(&type2);
        MPI_Type_free(&type1);
        type1 = type2;
    }

    /* disps[0] is the displacement to the beginning of the global array */
    disps[0] = 0;

    /* disps[1] is the first byte displacement of the subarray */
    disps[1] = array_of_starts[ndims-1] * extent;
    size = extent;
    for (i=ndims-2; i>=0; i--) {
        size *= array_of_sizes[i+1];
        disps[1] += size * array_of_starts[i];
    }

    /* disps[2] is the size of the global array */
    disps[2] = array_size;

    /* make filetype the same as calling MPI_Type_create_subarray() */
    /* adjust LB and UB without using MPI_LB or MPI_UB */
    err = MPI_Type_create_hindexed(1, blklens, &disps[1], type1, &type2);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_hindexed");
    MPI_Type_commit(&type2);
    err = MPI_Type_create_resized(type2, disps[0], disps[2], newtype);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_resized");
    MPI_Type_free(&type2);
    MPI_Type_free(&type1);

    return NC_NOERR;
#endif
}

/*----< filetype_create_vara() >--------------------------------------------*/
static int
filetype_create_vara(const NC         *ncp,
                     const NC_var     *varp,
                     const MPI_Offset *start,
                     const MPI_Offset *count,
                     MPI_Offset       *offset_ptr,         /* OUT */
                     MPI_Datatype     *filetype_ptr,       /* OUT */
                     int              *is_filetype_contig) /* OUT */
{
    int          status, err;
    MPI_Offset   offset;
    MPI_Datatype filetype, xtype;

    *offset_ptr   = varp->begin;
    *filetype_ptr = MPI_BYTE;
    if (is_filetype_contig != NULL) *is_filetype_contig = 1;

    /* when varp is a scalar, no need to create a filetype */
    if (varp->ndims == 0) return NC_NOERR;

    /* if the request is contiguous in file, no need to create a filetype */
    if (is_request_contiguous(IS_RECVAR(varp), ncp->vars.num_rec_vars,
                              varp->ndims, varp->shape, start, count)) {
        /* find the starting file offset of this request */
        status = ncmpio_first_offset(ncp, varp, start, &offset);
        *offset_ptr = offset;
        return status;
    }

    /* hereinafter fileview is noncontiguous, i.e. filetype != MPI_BYTE */
    if (is_filetype_contig != NULL) *is_filetype_contig = 0;
    offset = varp->begin;

    /* element MPI data type of variable */
    xtype = ncmpii_nc2mpitype(varp->xtype);

    /* previously, request size has been checked and it must > 0 */
    if (IS_RECVAR(varp)) {
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count blocklength;
#else
        int blocklength;
#endif
        MPI_Datatype rectype=MPI_BYTE;

#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
        /* check overflow only if MPI_Aint is smaller than MPI_Offset */
        status = check_recsize_too_big(ncp->recsize);
        if (status != NC_NOERR) return status;
#endif
        offset += start[0] * ncp->recsize;

        if (varp->ndims > 1) {
            /* when ndims > 1, we first construct a subarray type for a
             * single record, i.e. for dimensions 1 ... ndims-1 */
            status = type_create_subarray64(varp->ndims-1, varp->shape+1,
                                            count+1, start+1, MPI_ORDER_C,
                                            xtype, &rectype);
            if (status != NC_NOERR) return status;

            MPI_Type_commit(&rectype);
            blocklength = 1;
        }
        else { /* no subarray datatype is needed */
            blocklength = varp->xsz;
        }

#ifdef HAVE_MPI_LARGE_COUNT
        /* concatenate number of count[0] subarray types into filetype */
        err = MPI_Type_create_hvector_c(count[0], blocklength, ncp->recsize,
                                        rectype, &filetype);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector_c");
#else
        /* check overflow, because 1st argument of hvector is of type int */
        if (count[0] > NC_MAX_INT) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

        /* concatenate number of count[0] subarray types into filetype */
        err = MPI_Type_create_hvector((int)count[0], blocklength, ncp->recsize,
                                      rectype, &filetype);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");
#endif

        if (rectype != MPI_BYTE) MPI_Type_free(&rectype);
    }
    else { /* for non-record variable, just create a subarray datatype */
        status = type_create_subarray64(varp->ndims, varp->shape, count, start,
                                        MPI_ORDER_C, xtype, &filetype);
        if (status != NC_NOERR) return status;
    }
    MPI_Type_commit(&filetype);

    *offset_ptr   = offset;
    *filetype_ptr = filetype;

    return NC_NOERR;
}

/*----< stride_flatten() >----------------------------------------------------*/
/* flatten a stride request into a list of offset-length pairs stored in
 * blocklens[] and disps[]
 */
static int
stride_flatten(const NC_var      *varp,
               MPI_Offset         recsize,  /* record size */
               const MPI_Offset  *start,    /* [ndim] starts of subarray */
               const MPI_Offset  *count,    /* [ndim] counts of subarray */
               const MPI_Offset  *stride,   /* [ndim] strides of subarray */
               MPI_Offset        *nblocks,  /* OUT: number of blocks */
               MPI_Offset       **blocklens,/* OUT: length of each block */
               MPI_Aint         **disps)    /* OUT: displacement of each block */
{
    int i, j, k, isRecVar, ndims, el_size;
    MPI_Offset *dimlen, seg_len, nstride, array_len, off, subarray_len;

    *blocklens = NULL;
    *disps = NULL;

    /* whether record variable */
    isRecVar = IS_RECVAR(varp);

    ndims = varp->ndims;

    /* scalar variables have been handled before this subroutine is called */
    assert (ndims > 0);

    /* array element size */
    el_size = varp->xsz;

    /* dimension lengths */
    dimlen = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * ndims);

    for (i=0; i<ndims; i++) dimlen[i] = varp->shape[i];
    /* for record variable, set dimlen[0] to the record size */
    if (isRecVar) dimlen[0] = recsize;

    /* calculate the number of offset-length pairs */
    *nblocks = (stride[ndims-1] == 1) ? 1 : count[ndims-1];
    for (i=0; i<ndims-1; i++) *nblocks *= count[i];
    if (*nblocks == 0) return 1;

    *blocklens = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * (*nblocks));
    *disps = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * (*nblocks));

    /* the length of all contiguous segments are of the same size */
    seg_len  = (stride[ndims-1] == 1) ? count[ndims-1] : 1;
    seg_len *= el_size;
    nstride  = (stride[ndims-1] == 1) ? 1 : count[ndims-1];

    for (i=0; i<*nblocks; i++) (*blocklens)[i] = seg_len;

    /* set the offset-length pairs for the lowest dimension */
    (*disps)[0] = start[ndims-1] * el_size;
    for (k=1; k<nstride; k++)
        (*disps)[k] = (*disps)[k-1] + stride[ndims-1] * el_size;

    /* done with the lowest dimension */
    ndims--;

    subarray_len = nstride;
    array_len = el_size;
    /* for higher dimensions */
    while (ndims > 0) {
        /* array_len is global array size from lowest up to ndims */
        array_len *= dimlen[ndims];

        /* off is the global array offset for this dimension, ndims-1
         * For record variable, dimlen[0] is the sum of single record sizes
         */
        if (ndims == 1 && isRecVar) off = start[0] * dimlen[0];
        else off = start[ndims-1] * array_len;

        /* update all offsets from lowest up to dimension ndims-1 */
        for (j=0; j<subarray_len; j++) (*disps)[j] += off;

        /* update each successive subarray of dimension ndims-1
         * For record variable, dimlen[0] is the sum of single record sizes
         */
        if (ndims == 1 && isRecVar) off = stride[0] * dimlen[0];
        else off = stride[ndims-1] * array_len;

        for (i=1; i<count[ndims-1]; i++) {
            for (j=0; j<subarray_len; j++)
                (*disps)[k++] = (*disps)[j] + off;

            if (ndims == 1 && isRecVar) off += stride[0] * dimlen[0];
            else off += stride[ndims-1] * array_len;
        }
        ndims--;  /* move to next higher dimension */
        subarray_len *= count[ndims];
    }
    NCI_Free(dimlen);

    return 1;
}

/*----< ncmpio_filetype_create_vars() >--------------------------------------*/
int
ncmpio_filetype_create_vars(const NC         *ncp,
                            const NC_var     *varp,
                            const MPI_Offset *start,
                            const MPI_Offset *count,
                            const MPI_Offset *stride,
                            MPI_Offset       *offset_ptr,         /* OUT */
                            MPI_Datatype     *filetype_ptr,       /* OUT */
                            int              *is_filetype_contig) /* OUT */
{
    int           err=NC_NOERR, mpireturn, dim, isLargeReq;
    MPI_Aint     *disps;
    MPI_Offset    i, nblocks, nelems, *blocklens;
    MPI_Datatype  filetype=MPI_BYTE;

    if (stride == NULL)
        return filetype_create_vara(ncp, varp, start, count, offset_ptr,
                                    filetype_ptr, is_filetype_contig);

    /* check if a true vars (skip stride[] when count[] == 1) */
    for (dim=0; dim<varp->ndims; dim++)
        if (count[dim] > 1 && stride[dim] > 1)
            break;

    if (dim == varp->ndims) /* not a true vars */
        return filetype_create_vara(ncp, varp, start, count, offset_ptr,
                                    filetype_ptr, is_filetype_contig);

    /* now stride[] indicates a non-contiguous fileview */

    /* calculate request amount */
    nelems = 1;
    for (dim=0; dim<varp->ndims; dim++) nelems *= count[dim];

    /* starting file offset of the variable */
    *offset_ptr = varp->begin;

    /* when nelems == 0 or varp is a scalar, i.e. varp->ndims == 0, no need to
     * create a filetype
     */
    if (varp->ndims == 0 || nelems == 0) {
        *filetype_ptr = MPI_BYTE;
        if (is_filetype_contig != NULL) *is_filetype_contig = 1;
        return NC_NOERR;
    }

    /* hereinafter fileview is noncontiguous, i.e. filetype != MPI_BYTE */
    if (is_filetype_contig != NULL) *is_filetype_contig = 0;

    /* flatten stride access into list of offsets and lengths stored in
     * disps[] and blocklens[], respectively.
     */
    stride_flatten(varp, ncp->recsize, start, count, stride, &nblocks,
                   &blocklens, &disps);

    /* the flattened list allows one single call to hindexed constructor.
     * We cannot use MPI_Type_indexed because displacement for the record
     * dimension may not be a multiple of varp->xtype
     */
    isLargeReq = 0;
    if (nblocks > NC_MAX_INT)
        isLargeReq = 1;
    else {
        for (i=0; i<nblocks; i++) {
            if (blocklens[i] > NC_MAX_INT) {
                isLargeReq = 1;
                break;
            }
        }
    }

    if (isLargeReq) {
#ifdef HAVE_MPI_TYPE_CREATE_HINDEXED_C
        MPI_Count *blocklens_c, *disps_c;
        blocklens_c = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * nelems);
        disps_c     = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * nelems);
        for (i=0; i<nelems; i++) {
            blocklens_c[i] = blocklens[i];
            disps_c[i] = disps[i];
        }
        mpireturn = MPI_Type_create_hindexed_c(nblocks, blocklens_c, disps_c,
                                               MPI_BYTE, &filetype);
        NCI_Free(blocklens_c);
        NCI_Free(disps_c);
#else
        *filetype_ptr = MPI_DATATYPE_NULL;
        DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
        mpireturn = MPI_SUCCESS;
#endif
    }
    else {
        int *blocklens_i;
        blocklens_i = (int*) NCI_Malloc(sizeof(int) * nelems);
        for (i=0; i<nelems; i++)
            blocklens_i[i] = (int)blocklens[i];
        mpireturn = MPI_Type_create_hindexed((int)nblocks, blocklens_i, disps,
                                             MPI_BYTE, &filetype);
        NCI_Free(blocklens_i);
    }

    if (disps != NULL) NCI_Free(disps);
    if (blocklens != NULL) NCI_Free(blocklens);

    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");
        *filetype_ptr = MPI_DATATYPE_NULL;
    }
    else {
        MPI_Type_commit(&filetype);
        *filetype_ptr = filetype;
    }

    return err;
}

/*----< ncmpio_file_set_view() >---------------------------------------------*/
/* This function handles the special case for root process for setting its
 * file view: to keeps the whole file header visible to the root process. This
 * is because the root process may update the number of records or attributes
 * into the file header while in data mode. In PnetCDF design, only root
 * process can read/write the file header.
 * This function is collective if called in collective data mode
 */
int
ncmpio_file_set_view(const NC     *ncp,
                     MPI_File      fh,
                     MPI_Offset   *offset,  /* IN/OUT */
                     MPI_Datatype  filetype)
{
    int rank, err, mpireturn, status=NC_NOERR;

    if (filetype == MPI_BYTE) {
        /* filetype is a contiguous space, make the whole file visible */
        TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE,
                                    "native", MPI_INFO_NULL);
        return NC_NOERR;
    }

    MPI_Comm_rank(ncp->comm, &rank);
    if (rank == 0) {
        /* prepend the whole file header to filetype */
        MPI_Datatype root_filetype=MPI_BYTE, ftypes[2];
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count blocklens[2];
        MPI_Count disps[2];
        blocklens[0] = ncp->begin_var;
#else
        int blocklens[2];
        MPI_Aint disps[2];

        /* check if header size > 2^31 */
        if (ncp->begin_var > NC_MAX_INT) {
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW);
            goto err_out;
        }

        blocklens[0] = (int)ncp->begin_var;
#endif

        /* first block is the header extent */
            disps[0] = 0;
           ftypes[0] = MPI_BYTE;

        /* second block is filetype, the subarray request(s) to the variable */
        blocklens[1] = 1;
            disps[1] = *offset;
           ftypes[1] = filetype;

#if !defined(HAVE_MPI_LARGE_COUNT) && (SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET)
        if (*offset > NC_MAX_INT) {
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW);
            goto err_out;
        }
#endif

#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_create_struct_c(2, blocklens, disps, ftypes,
                                             &root_filetype);
#else
        mpireturn = MPI_Type_create_struct(2, blocklens, disps, ftypes,
                                           &root_filetype);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_struct");
            if (status == NC_NOERR) status = err;
        }
        MPI_Type_commit(&root_filetype);

#ifndef HAVE_MPI_LARGE_COUNT
err_out:
#endif
        TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, root_filetype, "native",
                                    MPI_INFO_NULL);
        if (root_filetype != MPI_BYTE)
            MPI_Type_free(&root_filetype);

        /* now update the explicit offset to be used in MPI-IO call later */
        *offset = ncp->begin_var;
    }
    else {
        TRACE_IO(MPI_File_set_view)(fh, *offset, MPI_BYTE, filetype, "native",
                                    MPI_INFO_NULL);
        /* the explicit offset is already set in fileview */
        *offset = 0;
    }
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_set_view");
        if (status == NC_NOERR) status = err;
    }

    return status;
}

