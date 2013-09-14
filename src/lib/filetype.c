/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include <mpi.h>

#include "nc.h"
#include "ncx.h"
#include "macro.h"

/* Prototypes for functions used only in this file */
static int check_recsize_too_big(NC *ncp);

/*----< check_recsize_too_big() >--------------------------------------------*/
static int
check_recsize_too_big(NC *ncp)
{
    int ret = NC_NOERR;
    /* assertion: because recsize will be used to set up the file
     * view, we must ensure there is no overflow when specifying
     * how big a stride there is between items (think interleaved
     * records).
     *
     * note: 'recsize' is the sum of the record size of all record
     * variables in this dataset */
    if (ncp->recsize != (MPI_Aint)ncp->recsize) {
        fprintf(stderr, "Type overflow: unable to read/write multiple records in this dataset\non this platform. Please either access records of this record variable\none-at-a-time or run on a 64 bit platform\n");
        ret = NC_ESMALL;
    }
    /* the assert here might harsh, but without it, users will get corrupt
     * data.  */
    assert (ncp->recsize == (MPI_Aint)ncp->recsize);
    return ret;
}

/*----< NCcoordck() >--------------------------------------------------------*/
/*
 * Check whether 'coord' values (indices) are valid for the variable.
 * Note that even if the request size is zero, this check is enforced in both
 * netCDF and PnetCDF. Otherwise, many test cases under test directoy can fail. 
 */
int
NCcoordck(NC               *ncp,
          const NC_var     *varp,
          const MPI_Offset *coord,  /* i.e. start[] */
          const int         rw_flag)
{
    const MPI_Offset *ip;
    MPI_Offset *up;

    if (varp->ndims == 0)
        return NC_NOERR;        /* 'scalar' variable */

    if (IS_RECVAR(varp)) {
        if (! fIsSet(ncp->flags, NC_64BIT_DATA)) /* not CDF-5 */
            if (*coord > X_UINT_MAX)
                return NC_EINVALCOORDS; /* sanity check */

        /* for record variable, [0] is the NC_UNLIMITED dimension */
        if (rw_flag == READ_REQ) {
            if (coord[0] >= ncp->numrecs)
                return NC_EINVALCOORDS;
        }
        /* In collective data mode where numrecs is always kept consistent
         * across memory, then there is no need to update numrecs.
         * (If NC_SHARE is set, then numrecs is even sync-ed with file.)
         *
         * In independent data mode, numrecs in memory across processes
         * and file can be inconsistent. Even re-reading numrecs from file
         * cannot get the lastest value, because in independent mode,
         * numrecs in file is not updated (due to race condition).
         * For example, a subset of processes write a new record and at
         * the same time another set writes 2 new records. Even if NC_SHARE
         * is set, new values of numrecs cannot be written to the file,
         * because it can cause a race consition (atomic read-modify IO is
         * required to solve this problem and MPI-IO cannot do it). Simply
         * said, numrecs is not automatically kept consistent in
         * independent mode. Users must call ncmpi_sync_numrecs()
         * collectively to sync the value. So, here what PnetCDF can do
         * best is just to check numrecs against the local value.
         */

        /* skip checking the record dimension */
        ip = coord + 1;
        up = varp->shape + 1;
    }
    else {
        ip = coord;
        up = varp->shape;
    }

    for (; ip < coord + varp->ndims; ip++, up++) {
        if ( (*ip < 0) || (*ip >= *up) )
            return NC_EINVALCOORDS;
    }
    return NC_NOERR;
}


/*----< NCedgeck() >---------------------------------------------------------*/
/*
 * Check whether 'edges' are valid for the variable and 'start'
 */
/*ARGSUSED*/
int
NCedgeck(const NC         *ncp,
         const NC_var     *varp,
         const MPI_Offset *start,
         const MPI_Offset *edges)  /* i.e. count[] */
{
    const MPI_Offset *const end = start + varp->ndims;
    const MPI_Offset *shp = varp->shape;

    if (varp->ndims == 0)
        return NC_NOERR;  /* 'scalar' variable */

    if (IS_RECVAR(varp)) {
        if (start[0] < 0) return NC_EEDGE;
        start++;
        edges++;
        shp++;
    }

    /* out-of-boundary and negative starts is illegal */
    for (; start < end; start++, edges++, shp++) {
        if ( (*shp < 0) || (*edges > *shp) ||
             (*start < 0) || (*start + *edges > *shp))
            return NC_EEDGE;
    }

    return NC_NOERR;
}

/*----< NCstrideedgeck() >---------------------------------------------------*/
int
NCstrideedgeck(const NC         *ncp,
               const NC_var     *varp,
               const MPI_Offset *start,
               const MPI_Offset *edges,  /* i.e. count[] */
               const MPI_Offset *stride)
{
    const MPI_Offset *const end = start + varp->ndims;
    const MPI_Offset *shp = varp->shape; /* use MPI_Offset for now :( */

    if (varp->ndims == 0)
        return NC_NOERR;  /* 'scalar' variable */

    if (IS_RECVAR(varp)) {
        if ( *stride == 0 ) /*|| *stride >= X_INT64_T_MAX)*/
            /* cast needed for braindead systems with signed MPI_Offset */
            return NC_ESTRIDE;

        start++;
        edges++;
        shp++;
        stride++;
    }

    for (; start < end; start++, edges++, shp++, stride++) {
        if ( (*shp < 0) ||
            (*edges > *shp) ||
            (*edges > 0 && *start+1 + (*edges-1) * *stride > *shp) ||
            (*edges == 0 && *start > *shp) )
            return NC_EEDGE;

        if ( *stride == 0)/* || *stride >= X_INT64_T_MAX)*/
            /* cast needed for braindead systems with signed MPI_Offset */
            return NC_ESTRIDE;
    }

    return NC_NOERR;
}

/*----< ncmpii_get_offset() >------------------------------------------------*/
int
ncmpii_get_offset(NC               *ncp,
                  NC_var           *varp,
                  const MPI_Offset  starts[],   /* [varp->ndims] */
                  const MPI_Offset  counts[],   /* [varp->ndims] */
                  const MPI_Offset  strides[],  /* [varp->ndims] */
                  const int         rw_flag,
                  MPI_Offset       *offset_ptr) /* return file offset */
{
    /* returns the file offset of the last element of this request */
    MPI_Offset offset, *end_off=NULL;
    int status, i, ndims;

    offset = varp->begin; /* beginning file offset of this variable */
    ndims  = varp->ndims; /* number of dimensions of this variable */

    if (counts != NULL)
        end_off = (MPI_Offset*) NCI_Malloc(ndims * sizeof(MPI_Offset));

    if (counts != NULL && strides != NULL) {
        for (i=0; i<ndims; i++)
            end_off[i] = starts[i] + (counts[i] - 1) * strides[i];
    }
    else if (counts != NULL) { /* strides == NULL */
        for (i=0; i<ndims; i++)
            end_off[i] = starts[i] + counts[i] - 1;
    }
    else { /* when counts == NULL strides is of no use */
        end_off = (MPI_Offset*) starts;
    }

    status = NCcoordck(ncp, varp, end_off, rw_flag);  /* validate end_off[] */
    if (status != NC_NOERR) {
#ifdef CDEBUG
        printf("ncmpii_get_offset(): NCcoordck() fails\n");
#endif
        return status;
    }

    if (ndims > 0) {
        if (IS_RECVAR(varp))
            /* no need to check recsize here: if MPI_Offset is only 32 bits we
               will have had problems long before here */
            offset += end_off[0] * ncp->recsize;
        else
            offset += end_off[ndims-1] * varp->xsz;

        if (ndims > 1) {
            if (IS_RECVAR(varp))
                offset += end_off[ndims - 1] * varp->xsz;
            else
                offset += end_off[0] * varp->dsizes[1] * varp->xsz;

            for (i=1; i<ndims-1; i++)
                offset += end_off[i] * varp->dsizes[i+1] * varp->xsz;
        }
    }
    if (counts != NULL && end_off != NULL)
        NCI_Free(end_off);

    *offset_ptr = offset;
    return NC_NOERR;
}

/*----< ncmpii_is_request_contiguous() >-------------------------------------*/
int
ncmpii_is_request_contiguous(NC_var           *varp,
                             const MPI_Offset  starts[],
                             const MPI_Offset  counts[])
{
    /* determine whether the get/put request to this variable using
       starts[] and counts[] is contiguous in file */
    int i, j, most_sig_dim, ndims=varp->ndims;

    /* this variable is a scalar */
    if (ndims == 0) return 1;

    for (i=0; i<ndims; i++)
         if (counts[i] == 0) /* zero length request */
             return 1;

    most_sig_dim = 0; /* record dimension */
    if (IS_RECVAR(varp)) {
        /* if there are more than one record variabl, then the record
           dimensions, counts[0] must == 1. For now, we assume there
           are more than one record variable.
           TODO: we need the API ncmpi_inq_rec() as in netcdf 3.6.3
                 to know how many record variables are defined */
        if (counts[0] > 1) return 0;
        most_sig_dim = 1;
    }

    for (i=ndims-1; i>most_sig_dim; i--) {
        /* find the first counts[i] that is not the entire dimension */
        if (counts[i] < varp->shape[i]) {
            /* check dim from i-1, i-2, ..., most_sig_dim and
               their counts[] should all be 1 */
            for (j=i-1; j>=most_sig_dim; j--) {
                if (counts[j] > 1)
                    return 0;
            }
            break;
        }
        else { /* counts[i] == varp->shape[i] */
            /* when accessing the entire dimension, starts[i] must be 0 */
            if (starts[i] != 0) return 0;
        }
    }
    return 1;
}

/*----< ncmpii_vara_create_filetype() >--------------------------------------*/
int
ncmpii_vara_create_filetype(NC               *ncp,
                            NC_var           *varp,
                            const MPI_Offset *start,
                            const MPI_Offset *count,
                            int               rw_flag,
                            int              *blocklen,
                            MPI_Offset       *offset_ptr,
                            MPI_Datatype     *filetype_ptr)
{
    int          dim, status;
    MPI_Offset   offset, nelems=1;
    MPI_Datatype filetype;

    offset   = varp->begin;
    filetype = MPI_BYTE;

    /* New coordinate/edge check to fix NC_EINVALCOORDS bug */
    status = NCedgeck(ncp, varp, start, count);
    if (status != NC_NOERR ||
        (rw_flag == READ_REQ && IS_RECVAR(varp) && *start + *count > NC_get_numrecs(ncp)))
    {
        status = NCcoordck(ncp, varp, start, rw_flag);
        if (status != NC_NOERR)
            return status;
        else
            return NC_EEDGE;
    }

    for (dim=0; dim<varp->ndims; dim++) nelems *= count[dim];

    *blocklen = varp->xsz * nelems;
    /* Warning: blocklen might overflow */

    /* check if the request is contiguous in file
       if yes, there is no need to create a filetype */
    if (ncmpii_is_request_contiguous(varp, start, count)) {
        status = ncmpii_get_offset(ncp, varp, start, NULL, NULL, rw_flag, &offset);
        *offset_ptr   = offset;
        *filetype_ptr = filetype;
        return status;
    }

    /* filetype is defined only when varp is not a scalar and
       the number of requested elemenst > 0
       (varp->ndims == 0 meaning this is a scalar variable)
       Otherwise, keep filetype MPI_BYTE
     */
    if (varp->ndims == 0 || nelems == 0) {
        *offset_ptr   = offset;
        *filetype_ptr = filetype;

        return NC_NOERR;
    }

    /* hereinafter varp->ndims > 0 && nelems > 0 */
    int i, ndims, blklens[3], tag=0;
    int *shape=NULL, *subcount=NULL, *substart=NULL; /* all in bytes */
    MPI_Offset *shape64=NULL, *subcount64=NULL, *substart64=NULL;
    MPI_Offset size, disps[3];
    MPI_Datatype rectype, types[3], type1;

    ndims    = varp->ndims;
    shape    = (int*) NCI_Malloc(3 * ndims * sizeof(int));
    subcount = shape    + ndims;
    substart = subcount + ndims;

    /* previously, request size has been checked and it must > 0 */
    if (IS_RECVAR(varp)) {
        /* TODO: check MPI_Offset-to-int overflow */
        subcount[0] = count[0];
        substart[0] = 0;
        shape[0]    = subcount[0];

        if (ncp->recsize <= varp->len) { /* the only record variable */
            if (varp->ndims == 1) {
                shape[0] *= varp->xsz;
                subcount[0] *= varp->xsz;
            }
            else {
                for (dim = 1; dim < ndims-1; dim++) {
                    shape[dim]    = varp->shape[dim];
                    subcount[dim] = count[dim];
                    substart[dim] = start[dim];
                }
                shape[dim]    = varp->xsz * varp->shape[dim];
                subcount[dim] = varp->xsz * count[dim];
                substart[dim] = varp->xsz * start[dim];
            }
            offset += start[0] * ncp->recsize;

            MPI_Type_create_subarray(ndims, shape, subcount, substart,
                                     MPI_ORDER_C, MPI_BYTE, &filetype);
            MPI_Type_commit(&filetype);
        }
        else { /* more than one record variables */
            check_recsize_too_big(ncp);

            offset += start[0] * ncp->recsize;
            if (varp->ndims == 1) {
#if (MPI_VERSION < 2)
                MPI_Type_hvector(subcount[0], varp->xsz, ncp->recsize,
                                 MPI_BYTE, &filetype);
#else
                MPI_Type_create_hvector(subcount[0], varp->xsz, ncp->recsize,
                                        MPI_BYTE, &filetype);
#endif
                MPI_Type_commit(&filetype);
            }
            else {
                for (dim = 1; dim < ndims-1; dim++) {
                    shape[dim]    = varp->shape[dim];
                    subcount[dim] = count[dim];
                    substart[dim] = start[dim];
                }
                shape[dim]    = varp->xsz * varp->shape[dim];
                subcount[dim] = varp->xsz * count[dim];
                substart[dim] = varp->xsz * start[dim];

                MPI_Type_create_subarray(ndims-1, shape+1, subcount+1, substart+1,
                                         MPI_ORDER_C, MPI_BYTE, &rectype);
                MPI_Type_commit(&rectype);
#if (MPI_VERSION < 2)
                MPI_Type_hvector(subcount[0], 1, ncp->recsize, rectype,
                                 &filetype);
#else
                MPI_Type_create_hvector(subcount[0], 1, ncp->recsize, rectype,
                                        &filetype);
#endif
                MPI_Type_commit(&filetype);
                MPI_Type_free(&rectype);
            }
        }
    }
    else { /* non record variable */
        tag = 0;
        for (dim=0; dim< ndims-1; dim++) {
            if (varp->shape[dim] > 2147483647) { /* if shape > 2^31-1 */
                tag = 1;
                break;
            }
        }
        if ((varp->shape[dim]*varp->xsz)  > 2147483647)
            tag = 1;

        if (tag == 0) {
            for (dim = 0; dim < ndims-1; dim++ ) {
                shape[dim]    = varp->shape[dim];
                subcount[dim] = count[dim];
                substart[dim] = start[dim];
            }
            shape[dim]    = varp->xsz * varp->shape[dim];
            subcount[dim] = varp->xsz * count[dim];
            substart[dim] = varp->xsz * start[dim];

            MPI_Type_create_subarray(ndims, shape, subcount, substart,
                                     MPI_ORDER_C, MPI_BYTE, &filetype);
            MPI_Type_commit(&filetype);
        }
        else {
            shape64 = (MPI_Offset*) NCI_Malloc(3 * ndims * sizeof(MPI_Offset));
            subcount64 = shape64    + ndims;
            substart64 = subcount64 + ndims;

            if (ndims == 1) {  // for 64-bit support,  added July 23, 2008
                shape64[0]    = varp->shape[0];
                subcount64[0] = count[0];
                substart64[0] = start[0];

                offset += start[0]*varp->xsz;

                MPI_Type_contiguous(subcount64[0]*varp->xsz, MPI_BYTE, &type1);
                MPI_Type_commit(&type1);
#if (MPI_VERSION < 2)
                /* Why use the arguments differently between MPI-1 and 2 ? */
                MPI_Type_hvector(subcount64[0], varp->xsz, shape64[0]*varp->xsz,
                                 MPI_BYTE, &filetype);
#else
                MPI_Type_create_hvector(1, 1, shape64[0]*varp->xsz,
                                        type1, &filetype);
#endif
                MPI_Type_commit(&filetype);
                MPI_Type_free(&type1);
            }
            else {
                for (dim = 0; dim < ndims-1; dim++ ) {
                    shape64[dim]    = varp->shape[dim];
                    subcount64[dim] = count[dim];
                    substart64[dim] = start[dim];
                }
                shape64[dim]    = varp->xsz * varp->shape[dim];
                subcount64[dim] = varp->xsz * count[dim];
                substart64[dim] = varp->xsz * start[dim];

                MPI_Type_hvector(subcount64[dim-1],
                                 subcount64[dim],
                                 varp->xsz * varp->shape[dim],
                                 MPI_BYTE,
                                 &type1);
                MPI_Type_commit(&type1);

                size = shape[dim];
                for (i=dim-2; i>=0; i--) {
                    size *= shape[i+1];
                    MPI_Type_hvector(subcount64[i],
                                     1,
                                     size,
                                     type1,
                                     &filetype);
                    MPI_Type_commit(&filetype);

                    MPI_Type_free(&type1);
                    type1 = filetype;
                }
                disps[1] = substart64[dim];
                size = 1;
                for (i=dim-1; i>=0; i--) {
                    size *= shape64[i+1];
                    disps[1] += size*substart64[i];
                }
                disps[2] = 1;
                for (i=0; i<ndims; i++) disps[2] *= shape64[i];

                disps[0] = 0;
                blklens[0] = blklens[1] = blklens[2] = 1;
                types[0] = MPI_LB;
                types[1] = type1;
                types[2] = MPI_UB;

                MPI_Type_struct(3,
                                blklens,
                                (MPI_Aint*) disps,
                                types,
                                &filetype);

                MPI_Type_free(&type1);
            }
            NCI_Free(shape64);
        }
    }
    NCI_Free(shape);

    if (filetype != MPI_BYTE) *blocklen = 1;

    *offset_ptr   = offset;
    *filetype_ptr = filetype;

    return NC_NOERR;
}

/*----< ncmpii_vars_create_filetype() >--------------------------------------*/
int
ncmpii_vars_create_filetype(NC               *ncp,
                            NC_var           *varp,
                            const MPI_Offset  start[],
                            const MPI_Offset  count[],
                            const MPI_Offset  stride[],
                            int               rw_flag,
                            int              *blocklen,
                            MPI_Offset       *offset_ptr,
                            MPI_Datatype     *filetype_ptr)
{
    int          dim, status;
    MPI_Offset   offset, stride_off, nelems=1;
    MPI_Datatype filetype;

    if (stride == NULL)
        return ncmpii_vara_create_filetype(ncp, varp, start, count, rw_flag,
                                           blocklen, offset_ptr, filetype_ptr);
    offset   = varp->begin;
    filetype = MPI_BYTE;

    for (dim=0; dim<varp->ndims && stride[dim]==1; dim++) ;

    if (dim == varp->ndims) /* if stride[] all == 1 */
        return ncmpii_vara_create_filetype(ncp, varp, start, count, rw_flag,
                                           blocklen, offset_ptr, filetype_ptr);

    /* now stride[] indicates a non-contiguous fileview, we must define an
     * MPI derived data type for the fileview. */
    *blocklen = 1;

    /* New coordinate/edge check to fix NC_EINVALCOORDS bug */
    status = NCedgeck(ncp, varp, start, count);
    if ((status != NC_NOERR) ||
        (rw_flag == READ_REQ && IS_RECVAR(varp) && *start + *count > NC_get_numrecs(ncp)))
    {
        status = NCcoordck(ncp, varp, start, rw_flag);
        if (status != NC_NOERR)
            return status;
        else
            return NC_EEDGE;
    }

    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR)
        return status;

    if ( rw_flag == READ_REQ && IS_RECVAR(varp) &&
        ( (*count > 0 && *start+1 + (*count-1) * *stride > NC_get_numrecs(ncp)) ||
          (*count == 0 && *start > NC_get_numrecs(ncp)) ) )
        return NC_EEDGE;

    for (dim=0; dim<varp->ndims; dim++) nelems *= count[dim];

    /* filetype is defined only when varp is not a scalar and
       the number of requested elemenst > 0
       (varp->ndims == 0 meaning this is a scalar variable)
       Otherwise, keep filetype MPI_BYTE
     */
    if (varp->ndims > 0 && nelems > 0) {
        int ndims, *blockcounts, *blocklens; 
        MPI_Aint *blockstride;
        MPI_Datatype tmptype;

        ndims       = varp->ndims;
        blockcounts = (int*) NCI_Malloc(2 * ndims * sizeof(int));
        blocklens   = blockcounts + ndims;

        blockstride = (MPI_Aint*) NCI_Malloc(ndims * sizeof(MPI_Aint));

        tmptype = MPI_BYTE;

        blockcounts[ndims-1] = count[ndims-1];
        /* check 4-byte integer overflow (blockcounts in MPI_Type_hvector
           is of type int while count[] is of type MPI_Offset */
        if (count[ndims-1] != blockcounts[ndims-1])
            return NC_EINTOVERFLOW;
        /* blocklens[] is unlikely a big value */
        blocklens[ndims-1]  = varp->xsz;

        if (ndims == 1 && IS_RECVAR(varp)) {
            check_recsize_too_big(ncp);
            stride_off = stride[ndims-1] * ncp->recsize;
            blockstride[ndims-1] = stride_off;
            offset += start[ndims - 1] * ncp->recsize;
        } else {
            stride_off = stride[ndims-1] * varp->xsz;
            blockstride[ndims-1] = stride_off;
            offset += start[ndims-1] * varp->xsz;
        }
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
        /* check 4-byte integer overflow */
        if (stride_off != blockstride[ndims-1])
            return NC_EINTOVERFLOW;
#endif

        for (dim=ndims-1; dim>=0; dim--) {
#if (MPI_VERSION < 2)
            MPI_Type_hvector(blockcounts[dim], blocklens[dim], blockstride[dim],
                             tmptype, &filetype);
#else
            MPI_Type_create_hvector(blockcounts[dim], blocklens[dim],
                                    blockstride[dim], tmptype, &filetype);
#endif
            MPI_Type_commit(&filetype);
            if (tmptype != MPI_BYTE)
                MPI_Type_free(&tmptype);
            tmptype = filetype;

            if (dim - 1 >= 0) {
                blocklens[dim-1]  = 1;
                blockcounts[dim-1] = count[dim - 1];
                /* check 4-byte integer overflow */
                if (count[dim-1] != blockcounts[dim-1])
                    return NC_EINTOVERFLOW;

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
                if (stride_off != blockstride[dim-1])
                    return NC_EINTOVERFLOW;
#endif
            }
        }
        NCI_Free(blockstride);
        NCI_Free(blockcounts);
    }

    *offset_ptr   = offset;
    *filetype_ptr = filetype;

    return NC_NOERR;
}

