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
#include <limits.h>  /* INT_MAX */
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
inline static int
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
        /* if there are more than one record variable, then the record
           dimensions, count[0] must == 1. For now, we assume there
           are more than one record variable.
           TODO: we may need an API to inquire how many record variables
           are defined */
        if (numRecVars > 1) { /* more than one record variable */
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

#ifndef HAVE_MPI_TYPE_CREATE_SUBARRAY
/*----< type_create_subarray() >---------------------------------------------*/
/* this is to be used when MPI_Type_create_subarray() is not available,
 * typically for MPI-1 implementation only */
static int
type_create_subarray(int           ndims,
                     const int    *array_of_sizes,    /* [ndims] */
                     const int    *array_of_subsizes, /* [ndims] */
                     const int    *array_of_starts,   /* [ndims] */
                     int           order,
                     MPI_Datatype  oldtype,
                     MPI_Datatype *newtype)
{
    int i, err, blklens[3] = {1, 1, 1};
    MPI_Datatype type1, type2;
    MPI_Aint extent, size, array_size, stride, disps[3];

    if (ndims == 0) DEBUG_RETURN_ERROR(NC_EDIMMETA)

#ifdef HAVE_MPI_TYPE_GET_EXTENT
    MPI_Aint lb;
    MPI_Type_get_extent(oldtype, &lb, &extent);
#else
    MPI_Type_extent(oldtype, &extent);
#endif
    array_size = extent;
    for (i=0; i<ndims; i++) array_size *= array_of_sizes[i];

    if (ndims == 1) {
        /* blklens argument in MPI_Type_create_hindexed() is of type int */
        blklens[1] = array_of_subsizes[0];
        disps[1] = extent * array_of_starts[0];

#if defined (HAVE_MPI_TYPE_CREATE_HINDEXED) && defined(HAVE_MPI_TYPE_CREATE_RESIZED)
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
#else
        /* add holes in the beginning and tail of oldtype */
        MPI_Datatype types[3];
        types[0] = MPI_LB; types[1] = oldtype; types[2] = MPI_UB;
        disps[0] = 0;                          disps[2] = array_size;
        err = MPI_Type_struct(3, blklens, disps, types, newtype);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_struct");
#endif
        return NC_NOERR;
    }
    /* now, ndims > 1 */

    /* first create a datatype for the least 2 significant dimensions */

    /* blklens argument in MPI_Type_create_hvector() is of type int */
    blklens[0] = array_of_subsizes[ndims-1];
    stride = array_of_sizes[ndims-1] * extent;
#ifdef HAVE_MPI_TYPE_CREATE_HVECTOR
    err = MPI_Type_create_hvector(array_of_subsizes[ndims-2], blklens[0],
                                  stride, oldtype, &type1);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");
#else
    err = MPI_Type_hvector(array_of_subsizes[ndims-2], blklens[0],
                           stride, oldtype, &type1);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_hvector");
#endif
    MPI_Type_commit(&type1);

    /* now iterate through the rest dimensions */
    for (i=ndims-3; i>=0; i--) {
        stride *= array_of_sizes[i+1];
#ifdef HAVE_MPI_TYPE_CREATE_HVECTOR
        err = MPI_Type_create_hvector(array_of_subsizes[i], 1, stride, type1, &type2);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");
#else
        err = MPI_Type_hvector(array_of_subsizes[i], 1, stride, type1, &type2);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_hvector");
#endif
        MPI_Type_commit(&type2);
        MPI_Type_free(&type1);
        type1 = type2;
    }

    /* disps[1] is the first byte displacement of the subarray */
    disps[1] = array_of_starts[ndims-1] * extent;
    size = 1;
    for (i=ndims-2; i>=0; i--) {
        size *= array_of_sizes[i+1];
        disps[1] += size * array_of_starts[i];
    }

    /* disps[2] is the size of the global array */
    disps[2] = array_size;

    /* disps[0] is the beginning of the global array */
    disps[0] = 0;

    /* make filetype the same as calling MPI_Type_create_subarray() */
    blklens[0] = 1;
#if defined (HAVE_MPI_TYPE_CREATE_HINDEXED) && defined(HAVE_MPI_TYPE_CREATE_RESIZED)
    /* adjust LB and UB without using MPI_LB or MPI_UB */
    err = MPI_Type_create_hindexed(1, blklens, &disps[1], type1, &type2);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_hindexed");
    MPI_Type_commit(&type2);
    err = MPI_Type_create_resized(type2, disps[0], disps[2], newtype);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_resized");
    MPI_Type_free(&type2);
#else
    MPI_Datatype types[3];
    types[0] = MPI_LB;
    types[1] = type1;
    types[2] = MPI_UB;
    err = MPI_Type_struct(3, blklens, disps, types, newtype);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_struct");
#endif
    MPI_Type_free(&type1);

    return NC_NOERR;
}
#endif

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
    int i, err, tag, blklens[3] = {1, 1, 1};
    MPI_Datatype type1, type2;
    MPI_Aint extent, size, array_size, stride, disps[3];
#ifdef HAVE_MPI_TYPE_GET_EXTENT
    MPI_Aint lb;
#endif

    if (ndims == 0) DEBUG_RETURN_ERROR(NC_EDIMMETA)

    /* check if any of the dimensions is larger than 2^31-1 */
    tag = 0;
    for (i=0; i<ndims; i++) {
        if (array_of_sizes[i] > 2147483647) {
            tag = 1;
            break;
        }
    }

    if (tag == 0) {
        /* none of dimensions > 2^31-1, we can safely use
         * MPI_Type_create_subarray */
        int *sizes    = (int*) NCI_Malloc((size_t)ndims * 3 * SIZEOF_INT);
        int *subsizes = sizes    + ndims;
        int *starts   = subsizes + ndims;
        for (i=0; i<ndims; i++) {
            sizes[i]    = (int)array_of_sizes[i];
            subsizes[i] = (int)array_of_subsizes[i];
            starts[i]   = (int)array_of_starts[i];
            if (array_of_sizes[i]    != sizes[i] ||
                array_of_subsizes[i] != subsizes[i] ||
                array_of_starts[i]   != starts[i])
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
#ifdef HAVE_MPI_TYPE_CREATE_SUBARRAY
        err = MPI_Type_create_subarray(ndims, sizes, subsizes, starts,
                                       order, oldtype, newtype);
        NCI_Free(sizes);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_subarray");
#else
        err = type_create_subarray(ndims, sizes, subsizes, starts,
                                   order, oldtype, newtype);
        NCI_Free(sizes);
#endif
        return err;
    }

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

#ifdef HAVE_MPI_TYPE_GET_EXTENT
    MPI_Type_get_extent(oldtype, &lb, &extent);
#else
    MPI_Type_extent(oldtype, &extent);
#endif
    array_size = extent;
    for (i=0; i<ndims; i++) array_size *= array_of_sizes[i];

    if (ndims == 1) {
        /* blklens argument in MPI_Type_create_hindexed() is of type int */
        blklens[1] = (int)array_of_subsizes[0];
        if (array_of_subsizes[0] != blklens[1]) /* check int overflow */
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        disps[1] = extent * array_of_starts[0];

#if defined (HAVE_MPI_TYPE_CREATE_HINDEXED) && defined(HAVE_MPI_TYPE_CREATE_RESIZED)
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
#else
        /* add holes in the beginning and tail of oldtype */
        MPI_Datatype types[3];
        types[0] = MPI_LB; types[1] = oldtype; types[2] = MPI_UB;
        disps[0] = 0;                          disps[2] = array_size;
        err = MPI_Type_struct(3, blklens, disps, types, newtype);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_struct");
#endif
        return NC_NOERR;
    }
    /* now, ndims > 1 */

    /* first create a datatype for the least 2 significant dimensions */

    /* count and blocklength arguments in MPI_Type_create_hvector() are of
     * type int. We need to check for integer overflow */
    int count, blocklength;
    count = (int)array_of_subsizes[ndims-2];
    blocklength = (int)array_of_subsizes[ndims-1];
    if (array_of_subsizes[ndims-2] != count ||
        array_of_subsizes[ndims-1] != blocklength) /* check int overflow */
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    stride = array_of_sizes[ndims-1] * extent;
#ifdef HAVE_MPI_TYPE_CREATE_HVECTOR
    err = MPI_Type_create_hvector(count, blocklength, stride, oldtype, &type1);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");
#else
    err = MPI_Type_hvector(count, blocklength, stride, oldtype, &type1);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_hvector");
#endif
    MPI_Type_commit(&type1);

    /* now iterate through the rest dimensions */
    for (i=ndims-3; i>=0; i--) {
        count = (int)array_of_subsizes[i];
        if (array_of_subsizes[i] != count) /* check int overflow */
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        stride *= array_of_sizes[i+1];
#ifdef HAVE_MPI_TYPE_CREATE_HVECTOR
        err = MPI_Type_create_hvector(count, 1, stride, type1, &type2);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");
#else
        err = MPI_Type_hvector(count, 1, stride, type1, &type2);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_hvector");
#endif
        MPI_Type_commit(&type2);
        MPI_Type_free(&type1);
        type1 = type2;
    }

    /* disps[0] is the displacement to the beginning of the global array */
    disps[0] = 0;

    /* disps[1] is the first byte displacement of the subarray */
    disps[1] = array_of_starts[ndims-1] * extent;
    size = 1;
    for (i=ndims-2; i>=0; i--) {
        size *= array_of_sizes[i+1];
        disps[1] += size * array_of_starts[i];
    }

    /* disps[2] is the size of the global array */
    disps[2] = array_size;

    /* make filetype the same as calling MPI_Type_create_subarray() */
#if defined(HAVE_MPI_TYPE_CREATE_HINDEXED) && defined(HAVE_MPI_TYPE_CREATE_RESIZED)
    /* adjust LB and UB without using MPI_LB or MPI_UB */
    err = MPI_Type_create_hindexed(1, blklens, &disps[1], type1, &type2);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_hindexed");
    MPI_Type_commit(&type2);
    err = MPI_Type_create_resized(type2, disps[0], disps[2], newtype);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_resized");
    MPI_Type_free(&type2);
#else
    MPI_Datatype types[3];
    types[0] = MPI_LB;
    types[1] = type1;
    types[2] = MPI_UB;
    err = MPI_Type_struct(3, blklens, disps, types, newtype);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_struct");
#endif
    MPI_Type_free(&type1);

    return NC_NOERR;
}

/*----< filetype_create_vara() >--------------------------------------------*/
static int
filetype_create_vara(const NC         *ncp,
                     const NC_var     *varp,
                     const MPI_Offset *start,
                     const MPI_Offset *count,
                     int              *blocklen,           /* OUT */
                     MPI_Offset       *offset_ptr,         /* OUT */
                     MPI_Datatype     *filetype_ptr,       /* OUT */
                     int              *is_filetype_contig) /* OUT */
{
    int          dim, status, err;
    MPI_Offset   nbytes, offset;
    MPI_Datatype filetype, xtype;

    *offset_ptr   = varp->begin;
    *filetype_ptr = MPI_BYTE;
    if (is_filetype_contig != NULL) *is_filetype_contig = 1;

    /* calculate the request size */
    nbytes = varp->xsz;
    for (dim=0; dim<varp->ndims; dim++) nbytes *= count[dim];
    if (nbytes > INT_MAX) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    if (blocklen != NULL) *blocklen = (int)nbytes;

    /* when varp is a scalar or nbytes == 0, no need to create a filetype */
    if (varp->ndims == 0 || nbytes == 0) return NC_NOERR;

    /* if the request is contiguous in file, no need to create a filetype */
    if (is_request_contiguous(IS_RECVAR(varp), ncp->vars.num_rec_vars,
                              varp->ndims, varp->shape, start, count)) {
        /* find the starting file offset of this request */
        status = ncmpio_first_offset(ncp, varp, start, &offset);
        *offset_ptr   = offset;
        return status;
    }

    /* hereinafter fileview is noncontiguous, i.e. filetype != MPI_BYTE.
     * Since we will construct a filetype, blocklen is set to 1.
     */
    if (blocklen           != NULL) *blocklen           = 1;
    if (is_filetype_contig != NULL) *is_filetype_contig = 0;
    offset = varp->begin;

    /* element MPI data type of variable */
    xtype = ncmpii_nc2mpitype(varp->xtype);

    /* previously, request size has been checked and it must > 0 */
    if (IS_RECVAR(varp)) {
        int blocklength;
        MPI_Datatype rectype=MPI_BYTE;

#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
        /* check overflow only if MPI_Aint is smaller than MPI_Offset */
        status = check_recsize_too_big(ncp->recsize);
        if (status != NC_NOERR) return status;
#endif
        /* check overflow, because 1st argument of hvector is of type int */
        if (count[0] != (int) count[0]) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

        offset += start[0] * ncp->recsize;

        if (varp->ndims > 1) {
            /* when ndims > 1, we first construct a subarray type for a
             * single record, i.e. for dimensions 1 ... ndims-1 */
            MPI_Offset *shape64, *subcount64, *substart64;

            shape64 = (MPI_Offset*) NCI_Malloc((size_t)varp->ndims * 3 * SIZEOF_MPI_OFFSET);
            subcount64 = shape64    + varp->ndims;
            substart64 = subcount64 + varp->ndims;

            shape64[0]    = count[0];
            subcount64[0] = count[0];
            substart64[0] = 0;

            for (dim=1; dim<varp->ndims; dim++) {
                shape64[dim]    = varp->shape[dim];
                subcount64[dim] = count[dim];
                substart64[dim] = start[dim];
            }
            status = type_create_subarray64(varp->ndims-1, shape64+1,
                                 subcount64+1, substart64+1, MPI_ORDER_C,
                                 xtype, &rectype);
            NCI_Free(shape64);
            if (status != NC_NOERR) return status;

            MPI_Type_commit(&rectype);
            blocklength = 1;
        }
        else { /* no subarray datatype is needed */
            blocklength = varp->xsz;
        }

        /* concatenate number of count[0] subarray types into filetype */
#ifdef HAVE_MPI_TYPE_CREATE_HVECTOR
        err = MPI_Type_create_hvector((int)count[0], blocklength, ncp->recsize,
                                      rectype, &filetype);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");
#else
        err = MPI_Type_hvector((int)count[0], blocklength, ncp->recsize,
                               rectype, &filetype);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Type_hvector");
#endif
        if (rectype != MPI_BYTE) MPI_Type_free(&rectype);
    }
    else { /* for non-record variable, just create a subarray datatype */
        MPI_Offset *shape64, *subcount64, *substart64;
        shape64 = (MPI_Offset*) NCI_Malloc((size_t)varp->ndims * 3 * SIZEOF_MPI_OFFSET);
        subcount64 = shape64    + varp->ndims;
        substart64 = subcount64 + varp->ndims;

        for (dim=0; dim<varp->ndims; dim++) {
            shape64[dim]    = varp->shape[dim];
            subcount64[dim] = count[dim];
            substart64[dim] = start[dim];
        }
        status = type_create_subarray64(varp->ndims, shape64, subcount64,
                                        substart64, MPI_ORDER_C,
                                        xtype, &filetype);
        NCI_Free(shape64);
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
stride_flatten(int               isRecVar, /* whether record variable */
               int               ndim,     /* number of dimensions */
               int               el_size,  /* array element size */
               const MPI_Offset *dimlen,   /* [ndim] dimension lengths */
               const MPI_Offset *start,    /* [ndim] starts of subarray */
               const MPI_Offset *count,    /* [ndim] counts of subarray */
               const MPI_Offset *stride,   /* [ndim] strides of subarray */
               int              *nblocks,  /* OUT: number of blocks */
               int              *blocklens,/* OUT: length of each block */
               MPI_Aint         *disps)    /* OUT: displacement of each block */
{
    int i, j, k, seg_len;
    MPI_Offset nstride, array_len, off, subarray_len;

    /* scalar variables have been handled before this subroutine is called */
    assert (ndim > 0);

    /* calculate the number of offset-length pairs */
    *nblocks = (stride[ndim-1] == 1) ? 1 : count[ndim-1];
    for (i=0; i<ndim-1; i++) *nblocks *= count[i];
    if (*nblocks == 0) return 1;

    /* the length of all contiguous segments are of the same size */
    seg_len  = (stride[ndim-1] == 1) ? count[ndim-1] : 1;
    seg_len *= el_size;
    nstride  = (stride[ndim-1] == 1) ? 1 : count[ndim-1];

    for (i=0; i<*nblocks; i++) blocklens[i] = seg_len;

    /* set the offset-length pairs for the lowest dimension */
    disps[0] = start[ndim-1] * el_size;
    for (k=1; k<nstride; k++)
        disps[k] = disps[k-1] + stride[ndim-1] * el_size;

    /* done with the lowest dimension */
    ndim--;

    subarray_len = nstride;
    array_len = el_size;
    /* for higher dimensions */
    while (ndim > 0) {
        /* array_len is global array size from lowest up to ndim */
        array_len *= dimlen[ndim];

        /* off is the global array offset for this dimension, ndim-1
         * For record variable, dimlen[0] is the sum of single record sizes
         */
        if (ndim == 1 && isRecVar) off = start[0] * dimlen[0];
        else off = start[ndim-1] * array_len;

        /* update all offsets from lowest up to dimension ndim-1 */
        for (j=0; j<subarray_len; j++) disps[j] += off;

        /* update each successive subarray of dimension ndim-1
         * For record variable, dimlen[0] is the sum of single record sizes
         */
        if (ndim == 1 && isRecVar) off = stride[0] * dimlen[0];
        else off = stride[ndim-1] * array_len;

        for (i=1; i<count[ndim-1]; i++) {
            for (j=0; j<subarray_len; j++)
                disps[k++] = disps[j] + off;

            if (ndim == 1 && isRecVar) off += stride[0] * dimlen[0];
            else off += stride[ndim-1] * array_len;
        }
        ndim--;  /* move to next higher dimension */
        subarray_len *= count[ndim];
    }
    return 1;
}

/*----< ncmpio_filetype_create_vars() >--------------------------------------*/
int
ncmpio_filetype_create_vars(const NC         *ncp,
                            const NC_var     *varp,
                            const MPI_Offset *start,
                            const MPI_Offset *count,
                            const MPI_Offset *stride,
                            int              *blocklen,           /* OUT */
                            MPI_Offset       *offset_ptr,         /* OUT */
                            MPI_Datatype     *filetype_ptr,       /* OUT */
                            int              *is_filetype_contig) /* OUT */
{
    int           dim, err, nblocks, *blocklens;
    MPI_Aint     *disps;
    MPI_Offset    offset, nelems, *shape;
    MPI_Datatype  filetype=MPI_BYTE;

    if (stride == NULL)
        return filetype_create_vara(ncp, varp, start, count,
                                    blocklen, offset_ptr, filetype_ptr,
                                    is_filetype_contig);

    /* check if a true vars (skip stride[] when count[] == 1) */
    for (dim=0; dim<varp->ndims; dim++)
        if (count[dim] > 1 && stride[dim] > 1)
            break;

    if (dim == varp->ndims) /* not a true vars */
        return filetype_create_vara(ncp, varp, start, count,
                                    blocklen, offset_ptr, filetype_ptr,
                                    is_filetype_contig);

    /* now stride[] indicates a non-contiguous fileview */

    /* calculate request amount */
    nelems = 1;
    for (dim=0; dim<varp->ndims; dim++) nelems *= count[dim];

    /* when nelems == 0 or varp is a scalar, i.e. varp->ndims == 0, no need to
     * create a filetype
     */
    if (varp->ndims == 0 || nelems == 0) {
        *offset_ptr   = varp->begin;
        *filetype_ptr = MPI_BYTE;
        if (blocklen           != NULL) *blocklen           = 0;
        if (is_filetype_contig != NULL) *is_filetype_contig = 1;
        return NC_NOERR;
    }

    /* hereinafter fileview is noncontiguous, i.e. filetype != MPI_BYTE.
     * Since we will construct a filetype, blocklen is set to 1.
     */
    if (blocklen           != NULL) *blocklen           = 1;
    if (is_filetype_contig != NULL) *is_filetype_contig = 0;
    offset = varp->begin;

#if 1
    blocklens = (int*) NCI_Malloc((size_t)nelems * SIZEOF_INT);
    disps = (MPI_Aint*) NCI_Malloc((size_t)nelems * SIZEOF_MPI_AINT);
    shape = (MPI_Offset*) NCI_Malloc((size_t)varp->ndims * SIZEOF_MPI_OFFSET);

    for (dim=0; dim<varp->ndims; dim++) shape[dim] = varp->shape[dim];
    /* for record variable, set shape[0] to the record size */
    if (IS_RECVAR(varp)) shape[0] = ncp->recsize;

    /* flatten stride access into list of offsets and lengths stored in
     * disps[] and blocklens[], respectively.
     */
    stride_flatten(IS_RECVAR(varp), varp->ndims, varp->xsz, shape, start,
                   count, stride, &nblocks, blocklens, disps);
    NCI_Free(shape);

    /* the flattened list allows one single call to hindexed constructor.
     * We cannot use MPI_Type_indexed because displacement for the record
     * dimension may not be a multiple of varp->xtype
     */
    err = MPI_Type_create_hindexed(nblocks, blocklens, disps, MPI_BYTE,
                                   &filetype);
    NCI_Free(disps);
    NCI_Free(blocklens);
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, "MPI_Type_create_hindexed");

    MPI_Type_commit(&filetype);

    /* Below is an old implementation that uses a nested hinexed constructor
     * which can MPI error on depth beyond  DLOOP_MAX_DATATYPE_DEPTH.
     * The above approach avoids such problem by flattening stride access into
     * offset-length and using a single call to hindexed constructor.
     */
#else
    int ndims, *blockcounts;
    MPI_Aint *blockstride;
    MPI_Offset   stride_off;
    MPI_Datatype tmptype;

    ndims       = varp->ndims;
    blockcounts = (int*) NCI_Malloc((size_t)ndims * 2 * SIZEOF_INT);
    blocklens   = blockcounts + ndims;

    blockstride = (MPI_Aint*) NCI_Malloc((size_t)ndims * SIZEOF_MPI_AINT);

    tmptype = MPI_BYTE;

    blockcounts[ndims-1] = (int)count[ndims-1];
    /* check 4-byte integer overflow (blockcounts in MPI_Type_hvector
       is of type int while count[] is of type MPI_Offset */
    if (count[ndims-1] != blockcounts[ndims-1]) {
        NCI_Free(blockstride);
        NCI_Free(blockcounts);
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    }
    /* blocklens[] is unlikely a big value */
    blocklens[ndims-1] = varp->xsz;

    if (ndims == 1 && IS_RECVAR(varp)) {
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
        err = check_recsize_too_big(ncp->recsize);
        if (err != NC_NOERR) return err;
#endif
        stride_off = stride[ndims-1] * ncp->recsize;
        blockstride[ndims-1] = stride_off;
        offset += start[ndims-1] * ncp->recsize;
    } else {
        stride_off = stride[ndims-1] * varp->xsz;
        blockstride[ndims-1] = stride_off;
        offset += start[ndims-1] * varp->xsz;
    }
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
    /* check 4-byte integer overflow */
    if (stride_off != blockstride[ndims-1]) {
        if (blockstride) NCI_Free(blockstride);
        if (blockcounts) NCI_Free(blockcounts);
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    }
#endif

    for (dim=ndims-1; dim>=0; dim--) {
#ifdef HAVE_MPI_TYPE_CREATE_HVECTOR
        err = MPI_Type_create_hvector(blockcounts[dim], blocklens[dim],
                                      blockstride[dim], tmptype, &filetype);
        if (err != MPI_SUCCESS) {
            NCI_Free(blockstride);
            NCI_Free(blockcounts);
            return ncmpii_error_mpi2nc(err, "MPI_Type_create_hvector");
        }
#else
        err = MPI_Type_hvector(blockcounts[dim], blocklens[dim],
                               blockstride[dim], tmptype, &filetype);
        if (err != MPI_SUCCESS) {
            NCI_Free(blockstride);
            NCI_Free(blockcounts);
            return ncmpii_error_mpi2nc(err, "MPI_Type_hvector");
        }
#endif
        MPI_Type_commit(&filetype);
        if (tmptype != MPI_BYTE)
            MPI_Type_free(&tmptype);
        tmptype = filetype;

        if (dim - 1 >= 0) {
            blocklens[dim-1]  = 1;
            blockcounts[dim-1] = (int)count[dim - 1];
            /* check 4-byte integer overflow */
            if (count[dim-1] != blockcounts[dim-1]) {
                NCI_Free(blockstride);
                NCI_Free(blockcounts);
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            }

            if (dim-1 == 0 && IS_RECVAR(varp)) {
                stride_off = stride[dim-1] * ncp->recsize;
                blockstride[dim-1] = stride_off;
                offset += start[dim-1] * ncp->recsize;
            } else {
                stride_off = stride[dim-1] * varp->dsizes[dim] * varp->xsz;
                blockstride[dim-1] = stride_off;
                offset += start[dim-1] * varp->dsizes[dim] * varp->xsz;
            }
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
            /* check 4-byte integer overflow */
            if (stride_off != blockstride[dim-1]) {
                NCI_Free(blockstride);
                NCI_Free(blockcounts);
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            }
#endif
        }
    }
    NCI_Free(blockstride);
    NCI_Free(blockcounts);
#endif

    *offset_ptr   = offset;
    *filetype_ptr = filetype;

    return NC_NOERR;
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
        int blocklens[2];
        MPI_Aint disps[2];
        MPI_Datatype root_filetype, ftypes[2];

        /* first block is the header extent */
        blocklens[0] = (int)ncp->begin_var;
            disps[0] = 0;
           ftypes[0] = MPI_BYTE;

        /* check if header size > 2^31 */
        if (ncp->begin_var != blocklens[0])
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)

        /* second block is filetype, the subarray request(s) to the variable */
        blocklens[1] = 1;
            disps[1] = *offset;
           ftypes[1] = filetype;

#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
        if (*offset != disps[1]) {
            blocklens[1] = 0;
            DEBUG_ASSIGN_ERROR(status, NC_EAINT_TOO_SMALL)
        }
#endif

#ifdef HAVE_MPI_TYPE_CREATE_STRUCT
        MPI_Type_create_struct(2, blocklens, disps, ftypes, &root_filetype);
#else
        MPI_Type_struct(2, blocklens, disps, ftypes, &root_filetype);
#endif
        MPI_Type_commit(&root_filetype);

        TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, root_filetype,
                                    "native", MPI_INFO_NULL);
        MPI_Type_free(&root_filetype);

        /* now update the explicit offset to be used in MPI-IO call later */
        *offset = ncp->begin_var;
    }
    else {
        TRACE_IO(MPI_File_set_view)(fh, *offset, MPI_BYTE, filetype,
                                    "native", MPI_INFO_NULL);
        /* the explicit offset is already set in fileview */
        *offset = 0;
    }
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_set_view");
        if (status == NC_NOERR) status = err;
    }

    return status;
}

