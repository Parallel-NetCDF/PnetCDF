/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>

#include <mpi.h>

#include <pnetcdf.h>
#include <pnc_debug.h>
#include <common.h>

/*@
  dtype_filter - Map a basic MPI datatype into one of the eight
  that we process natively.

  We unfortunately need to wrap all these types in case they aren't defined.

  Return:
  On success, one of MPI_CHAR, MPI_SIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_BYTE,
  MPI_SHORT, MPI_INT, MPI_LONG, MPI_FLOAT, and MPI_DOUBLE.
  On failure, MPI_DATATYPE_NULL.
@*/
static MPI_Datatype
dtype_filter(MPI_Datatype type)
{
    /* char types */
#ifdef HAVE_DECL_MPI_CHARACTER
    if (type == MPI_CHARACTER)
        return  MPI_CHAR;
#endif
    if (type == MPI_CHAR)
        return  MPI_CHAR;

    /* unsigned char types */
    if (type == MPI_UNSIGNED_CHAR)
        return  MPI_UNSIGNED_CHAR;

    /* signed char types */
    if (type == MPI_SIGNED_CHAR)
        return  MPI_SIGNED_CHAR;

#ifdef HAVE_DECL_MPI_INTEGER1
    if (type == MPI_INTEGER1)
        return  MPI_SIGNED_CHAR;
#endif
    if (type == MPI_BYTE)
        return  MPI_BYTE;

    /* 2-byte integer types (only supported if MPI_SHORT is 2-bytes). */
#if (SIZEOF_SHORT == 2)
  #ifdef HAVE_DECL_MPI_INTEGER2
    if (type == MPI_INTEGER2)
        return  MPI_SHORT;
  #endif
    if (type == MPI_SHORT)
        return  MPI_SHORT;

    if (type == MPI_UNSIGNED_SHORT)
        return  MPI_UNSIGNED_SHORT;
#endif

    /* 4-byte integer types */
  {
    MPI_Datatype int_4byte;
#if (SIZEOF_LONG == 4)
    MPI_Datatype uint_4byte;
#endif
#if (SIZEOF_INT == 4)
     int_4byte = MPI_INT;
#if (SIZEOF_LONG == 4)
    uint_4byte = MPI_UNSIGNED;
#endif
#elif (SIZEOF_SHORT == 4)
     int_4byte = MPI_SHORT;
#if (SIZEOF_LONG == 4)
    uint_4byte = MPI_UNSIGNED_SHORT;
#endif
#else
     int_4byte = MPI_DATATYPE_NULL; /* no 4-byte type? */
#if (SIZEOF_LONG == 4)
    uint_4byte = MPI_DATATYPE_NULL;
#endif
#endif

#ifdef HAVE_DECL_MPI_INTEGER
    if (type == MPI_INTEGER)
        return int_4byte;
#endif
#ifdef HAVE_DECL_MPI_INTEGER4
    if (type == MPI_INTEGER4)
        return int_4byte;
#endif
#if (SIZEOF_LONG == 4)
    if (type == MPI_LONG)
        return int_4byte;
    if (type == MPI_UNSIGNED_LONG)
        return uint_4byte;
#endif
#if (SIZEOF_INT == 4)
    if (type == MPI_INT)
        return  MPI_INT;
    if (type == MPI_UNSIGNED)
        return  MPI_UNSIGNED;
#elif (SIZEOF_SHORT == 4)
    if (type == MPI_SHORT)
        return  MPI_SHORT;
    if (type == MPI_UNSIGNED_SHORT)
        return  MPI_UNSIGNED_SHORT;
#endif
  }

     /* 8-byte integer types (only supported if MPI_INT or MPI_LONG
      * is 8-bytes).
      */
#if (SIZEOF_INT == 8) || (SIZEOF_LONG == 8)
  #ifdef HAVE_DECL_MPI_INTEGER8
    if (type == MPI_INTEGER8)
      #if (SIZEOF_INT == 8)
        return MPI_INT;
      #else
        return MPI_LONG;
      #endif
  #endif
  #ifdef HAVE_DECL_MPI_UNSIGNED_INTEGER8
    if (type == MPI_UNSIGNED_INTEGER8)
      #if (SIZEOF_UINT == 8)
        return MPI_UNSIGNED;
      #else
        return MPI_UNSIGNED_LONG;
      #endif
  #endif

  #if (SIZEOF_INT == 8)
    if (type == MPI_INT)
        return  MPI_INT;
    if (type == MPI_UNSIGNED)
        return  MPI_UNSIGNED;
  #endif
  #if (SIZEOF_LONG == 8)
    if (type == MPI_LONG)
      #if (SIZEOF_INT == 8)
        return MPI_INT;
      #else
        return MPI_LONG;
      #endif
    if (type == MPI_UNSIGNED_LONG)
      #if (SIZEOF_UINT == 8)
        return MPI_UNSIGNED;
      #else
        return MPI_UNSIGNED_LONG;
      #endif
  #endif
#endif

    /* 4-byte float types (we assume float is 4-bytes). */
#ifdef HAVE_DECL_MPI_REAL
    if (type == MPI_REAL)
        return  MPI_FLOAT;
#endif
#ifdef HAVE_DECL_MPI_REAL4
    if (type == MPI_REAL4)
        return  MPI_FLOAT;
#endif
    if (type == MPI_FLOAT)
        return  MPI_FLOAT;

    /* 8-byte float types (we assume double is 8-bytes). */
#ifdef HAVE_DECL_MPI_REAL8
    if (type == MPI_REAL8)
        return MPI_DOUBLE;
#endif
#ifdef HAVE_DECL_MPI_DOUBLE_PRECISION
    if (type == MPI_DOUBLE_PRECISION)
        return  MPI_DOUBLE;
#endif
    if (type == MPI_DOUBLE)
        return  MPI_DOUBLE;

    if (type == MPI_LONG_LONG_INT)
        return  MPI_LONG_LONG_INT;

    if (type == MPI_UNSIGNED_LONG_LONG)
        return  MPI_UNSIGNED_LONG_LONG;

/* MPI_LB and MPI_UB are deprecated since MPI-2.0 */
#if !defined(MPI_VERSION) || (MPI_VERSION < 2)
/* HP-MPI for example needs to handle MPI_LB and MPI_UB */
#ifdef HAVE_DECL_MPI_LB
    if (type == MPI_LB) {
#if SIZEOF_SIZE_T == 8
        return MPI_DOUBLE;
#elif SIZEOF_SIZE_T == 4
        return MPI_DOUBLE;
#endif
    }
#endif
#ifdef HAVE_DECL_MPI_UB
    if (type == MPI_UB) {
#if SIZEOF_SIZE_T == 8
        return MPI_DOUBLE;
#elif SIZEOF_SIZE_T == 4
        return MPI_DOUBLE;
#endif
    }
#endif
#endif

/* default */
    return MPI_DATATYPE_NULL;
}

/*@
  ncmpi_darray_get_totalblock - Count the total number of blocks assigned
  to the calling process by the darray datatype.

  Input:
. rank - rank of calling process
. ndims - number of dimensions for the array and process grid
. array_of_gsizes - Number of elements of type oldtype in each dimension
                    of global array (array of positive integers)
. array_of_distribs - Distribution of array in each dimension (array of state)
. array_of_dargs - Distribution argument in each dimension
                   (array of positive integers)
. array_of_psizes - Size of process grid in each dimension
                    (array of positive integers)

  Return:
. total number of blocks assigned from the distributed array
@*/
#ifdef HAVE_DECL_MPI_COMBINER_DARRAY
static int
darray_get_totalblks(int rank,
                     MPI_Offset ndims,
                     int array_of_gsizes[],
                     int array_of_distribs[],
                     int array_of_dargs[],
                     int array_of_psizes[])
{
  int total_blocks = 1;
  int subblocks, cycle, remain_cycle;
  int i;
  int pcoord;

  /* Process Grid ranking is always in C ORDER,
     so compute proc coordinates from last dim */

  for (i=(int)ndims-1; i>=0; i--) {
    if ( array_of_distribs[i] == MPI_DISTRIBUTE_NONE ) {
      total_blocks *= array_of_gsizes[i];
    } else {

      pcoord = rank % array_of_psizes[i];
      rank /= array_of_psizes[i];
      if ( array_of_dargs[i] == MPI_DISTRIBUTE_DFLT_DARG ) {
        subblocks = array_of_gsizes[i]/array_of_psizes[i];
        if ( subblocks*array_of_psizes[i]+pcoord < array_of_gsizes[i] )
          subblocks++;
      } else {
        cycle = array_of_dargs[i]*array_of_psizes[i];
        remain_cycle = array_of_gsizes[i] % cycle;

        subblocks = remain_cycle - pcoord*array_of_dargs[i];
        if (subblocks > array_of_dargs[i])
          subblocks = array_of_dargs[i];
        else if (subblocks < 0)
          subblocks = 0;

        subblocks += array_of_gsizes[i]/cycle * array_of_dargs[i];
      }

      if (subblocks == 0)
        return 0;
      total_blocks *= subblocks;

    }
  }

  return total_blocks;
}
#endif


/*----< ncmpii_dtype_decode() >----------------------------------------------*/
/*@
  ncmpii_dtype_decode - Decode the MPI derived datatype to get the bottom
  level primitive datatype information.

  Input:
. dtype - The MPI derived datatype to be decoded (can be predefined type).

  Output:
. ptype - The bottom level primitive datatype (only one allowed) in buftype
. el_size - The size of the ptype
. nelems - Number of elements/entries of such ptype in one buftype object
. iscontig_of_ptypes - Whether dtype is a contiguous number of ptype
@*/
int ncmpii_dtype_decode(MPI_Datatype  dtype,
                        MPI_Datatype *ptype,
                        int          *el_size,
                        MPI_Offset   *nelems,
                        int          *isderived,
                        int          *iscontig_of_ptypes)
{
    int i, ndims, status=NC_NOERR, tmpel_size;
    int num_ints, num_adds, num_dtypes, combiner;
    int *array_of_ints;
    MPI_Offset count, tmpnelems, total_blocks;
    MPI_Datatype tmpptype, *array_of_dtypes;
    MPI_Aint *array_of_adds;

    *isderived = 0;

    /* Note to PnetCDF developers. High-level APIs should never reach this
     * subroutine. Check whether it is a high-level of flexible API should
     * be done before entering here.
     */
    if (dtype == MPI_DATATYPE_NULL) {
        *nelems = 0;
        *ptype = dtype;
        *el_size = 0;
        *iscontig_of_ptypes = 1;
        return NC_NOERR;
    }

    MPI_Type_get_envelope(dtype, &num_ints, &num_adds, &num_dtypes, &combiner);

    if (
#ifdef HAVE_DECL_MPI_COMBINER_F90_INTEGER
        combiner == MPI_COMBINER_F90_INTEGER ||
#endif
#ifdef HAVE_DECL_MPI_COMBINER_F90_REAL
        combiner == MPI_COMBINER_F90_REAL ||
#endif
#ifdef HAVE_DECL_MPI_COMBINER_F90_COMPLEX
        combiner == MPI_COMBINER_F90_COMPLEX ||
#endif
        0) {
        fprintf(stderr, "FIXME: F90_INTEGER, F90_REAL or F90_COMPLEX are not supported.\n");
        DEBUG_RETURN_ERROR(NC_EUNSPTETYPE)
    }

    if (combiner == MPI_COMBINER_NAMED) { /* Predefined datatype */
        *nelems = 1;
        *ptype = dtype_filter(dtype);
        if (*ptype == MPI_DATATYPE_NULL)
            DEBUG_RETURN_ERROR(NC_EUNSPTETYPE)
        MPI_Type_size(dtype, el_size);
        *iscontig_of_ptypes = 1;
        return NC_NOERR;
    }

    array_of_ints = (int *) NCI_Malloc((size_t)num_ints * SIZEOF_INT);
    array_of_adds = (MPI_Aint *) NCI_Malloc((size_t)num_adds * SIZEOF_MPI_AINT);
    array_of_dtypes = (MPI_Datatype *) NCI_Malloc((size_t)num_dtypes * sizeof(MPI_Datatype));

    MPI_Type_get_contents(dtype, num_ints, num_adds, num_dtypes,
                          array_of_ints, array_of_adds, array_of_dtypes);

    switch (combiner) {
        /* single etype */
        case MPI_COMBINER_CONTIGUOUS:
        case MPI_COMBINER_HVECTOR:
        case MPI_COMBINER_VECTOR:
        case MPI_COMBINER_HINDEXED:
        case MPI_COMBINER_INDEXED:
#ifdef HAVE_DECL_MPI_COMBINER_DUP
        case MPI_COMBINER_DUP:
#endif
#ifdef HAVE_DECL_MPI_COMBINER_HVECTOR_INTEGER
        case MPI_COMBINER_HVECTOR_INTEGER:
#endif
#ifdef HAVE_DECL_MPI_COMBINER_INDEXED_BLOCK
        case MPI_COMBINER_INDEXED_BLOCK:
#endif
#ifdef HAVE_DECL_MPI_COMBINER_HINDEXED_INTEGER
        case MPI_COMBINER_HINDEXED_INTEGER:
#endif
#ifdef HAVE_DECL_MPI_COMBINER_SUBARRAY
        case MPI_COMBINER_SUBARRAY:
#endif
#ifdef HAVE_DECL_MPI_COMBINER_DARRAY
        case MPI_COMBINER_DARRAY:
#endif
#ifdef HAVE_DECL_MPI_COMBINER_RESIZED
        case MPI_COMBINER_RESIZED:
#endif
            status = ncmpii_dtype_decode(array_of_dtypes[0], ptype, el_size,
                                         nelems, isderived, iscontig_of_ptypes);
            if (*isderived) MPI_Type_free(array_of_dtypes);
            break;

        /* multiple etypes */
        case MPI_COMBINER_STRUCT:
#ifdef HAVE_DECL_MPI_COMBINER_STRUCT_INTEGER
        case MPI_COMBINER_STRUCT_INTEGER:
#endif
            count = array_of_ints[0];
            *el_size = 0;
            for (i=0; i<count && *el_size==0; i++) {
                /* need to skip null/marker types */
                status = ncmpii_dtype_decode(array_of_dtypes[i], ptype,
                                             el_size, nelems, isderived,
                                             iscontig_of_ptypes);
                if (status != NC_NOERR) return status;
                if (*isderived) MPI_Type_free(array_of_dtypes+i);
                if (*el_size > 0) *nelems *= array_of_ints[1+i];
            }
            for ( ; i<count; i++) {
                status = ncmpii_dtype_decode(array_of_dtypes[i], &tmpptype,
                                             &tmpel_size, &tmpnelems,
                                             isderived, iscontig_of_ptypes);
                if (status != NC_NOERR) return status;
                if (*isderived) MPI_Type_free(array_of_dtypes+i);
                if (tmpel_size > 0) {
                    if (tmpptype != *ptype)
                        DEBUG_RETURN_ERROR(NC_EMULTITYPES)
                    *nelems += tmpnelems*array_of_ints[1+i];
                }
            }
            *iscontig_of_ptypes = 0;
            break;
        default: break;
    }
    *isderived = 1;

    switch (combiner) {
        /* single etype */
        case MPI_COMBINER_CONTIGUOUS:
            total_blocks = array_of_ints[0];
            break;
        case MPI_COMBINER_HVECTOR:
        case MPI_COMBINER_VECTOR:
#ifdef HAVE_DECL_MPI_COMBINER_HVECTOR_INTEGER
        case MPI_COMBINER_HVECTOR_INTEGER:
#endif
#ifdef HAVE_DECL_MPI_COMBINER_INDEXED_BLOCK
        case MPI_COMBINER_INDEXED_BLOCK:
#endif
            *iscontig_of_ptypes = 0;
            total_blocks = (MPI_Offset)array_of_ints[0]*array_of_ints[1];
            break;
        case MPI_COMBINER_HINDEXED:
        case MPI_COMBINER_INDEXED:
#ifdef HAVE_DECL_MPI_COMBINER_HINDEXED_INTEGER
        case MPI_COMBINER_HINDEXED_INTEGER:
#endif
            *iscontig_of_ptypes = 0;
            for (i=0, total_blocks=0; i<array_of_ints[0]; i++)
                total_blocks += array_of_ints[1+i];
            break;
#ifdef HAVE_DECL_MPI_COMBINER_SUBARRAY
        case MPI_COMBINER_SUBARRAY:
            *iscontig_of_ptypes = 0;
            ndims = array_of_ints[0];
            for (i=0, total_blocks=1; i<ndims; i++)
                total_blocks *= array_of_ints[1+ndims+i];
            break;
#endif
#ifdef HAVE_DECL_MPI_COMBINER_DARRAY
        case MPI_COMBINER_DARRAY:
            *iscontig_of_ptypes = 0;
            ndims = array_of_ints[2];
            /* seldom reached, so put it in a separate function */
            total_blocks = darray_get_totalblks(array_of_ints[1],
                                                ndims,
                                                array_of_ints+3,
                                                array_of_ints+3+ndims,
                                                array_of_ints+3+2*ndims,
                                                array_of_ints+3+3*ndims);
            break;
#endif
#ifdef HAVE_DECL_MPI_COMBINER_RESIZED
        case MPI_COMBINER_RESIZED:
            *iscontig_of_ptypes = 0;
            total_blocks = 1;
            break;
#endif
        default: /* DUP etc. */
            total_blocks = 1;
            break;
    }
    *nelems *= total_blocks;

    NCI_Free(array_of_ints);
    NCI_Free(array_of_adds);
    NCI_Free(array_of_dtypes);

    return status;
}

/*----< ncmpii_buftype_decode() >--------------------------------------------*/
/* Obtain the following metadata about buftype:
 * etype:    element data type (MPI primitive type) in buftype
 * bufcount: If it is -1, then this is called from a high-level API and in
 *           this case buftype will be an MPI primitive data type.
 *           If bufcount is not -1, then this is called from a flexible API.
 * nelems:   number of etypes in user buffer
 * xnbytes:  number of bytes (in external data representation) to read/write
 *           from/to the file
 * esize:    byte size of etype
 * isContig: whether buftype is contiguous
 */
int
ncmpii_buftype_decode(int               ndims,
                      nc_type           xtype,
                      const MPI_Offset *count,
                      MPI_Offset        bufcount,
                      MPI_Datatype      buftype,
                      MPI_Datatype     *etype,    /* out */
                      int              *esize,    /* out */
                      MPI_Offset       *nelems,   /* out */
                      MPI_Offset       *xnbytes,  /* out */
                      int              *isContig) /* out */
{
    int i, xsz, err;
    MPI_Offset fnelems;

    err = ncmpii_xlen_nc_type(xtype, &xsz);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    /* fnelems is the total number of nc_type elements calculated from
     * count[]. count[] is the access count to the variable defined in
     * the netCDF file.
     */
    fnelems = 1;
    for (i=0; i<ndims; i++)
        fnelems *= count[i];

    if (bufcount == -1) { /* the subroutine is called from a high-level API */
        *nelems   = fnelems;
        *etype    = buftype; /* buftype is an MPI primitive data type */
        MPI_Type_size(buftype, esize);
        *xnbytes  = *nelems * xsz;
        *isContig = 1;
    }
    else if (buftype == MPI_DATATYPE_NULL) {
        /* This is called from a flexible API and buftype is set by user to
         * MPI_DATATYPE_NULL. In this case, bufcount is ignored and set by
         * this subroutine to a number match count[], and etype to match the
         * variable's external NC data type.
         */
        *nelems   = fnelems;
        *etype    = ncmpii_nc2mpitype(xtype);
        *esize    = xsz;
        *xnbytes  = *nelems * xsz;
        *isContig = 1;
    }
    else { /* This is called from a flexible API */
        int isderived;
        /* check some metadata of the MPI derived datatype */
        err = ncmpii_dtype_decode(buftype, etype, esize, nelems, &isderived,
                                  isContig);
        if (err != NC_NOERR) return err;

        /* make nelems the number of etype in the whole user buf */
        *nelems  *= bufcount;
        *xnbytes  = *nelems * xsz;

        /* check mismatch between nelems and fnelems */
        if (fnelems != *nelems) DEBUG_RETURN_ERROR(NC_EIOMISMATCH)
    }
    return NC_NOERR;
}

