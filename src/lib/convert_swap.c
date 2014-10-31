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
#include <arpa/inet.h>   /* htonl(), htons() */

#include <mpi.h>

#include "nc.h"
#include "ncx.h"

/* Prototypes for functions used only in this file */
#if 0
static void swapn(void *dst, const void *src, MPI_Offset nn, int xsize);
#endif

/*
 *  Datatype Mapping:
 *
 *  NETCDF    <--> MPI                    Description
 *   NC_BYTE       MPI_BYTE               signed 1-byte integer
 *   NC_CHAR       MPI_CHAR               char, text
 *   NC_SHORT      MPI_SHORT              signed 2-byte integer
 *   NC_INT        MPI_INT                signed 4-byte integer
 *   NC_FLOAT      MPI_FLOAT              single precision floating point
 *   NC_DOUBLE     MPI_DOUBLE             double precision floating point
 *   NC_UBYTE      MPI_UNSIGNED_CHAR      unsigned 1-byte int
 *   NC_USHORT     MPI_UNSIGNED_SHORT     unsigned 2-byte int
 *   NC_UINT       MPI_UNSIGNED           unsigned 4-byte int
 *   NC_INT64      MPI_LONG_LONG_INT      signed 8-byte int
 *   NC_UINT64     MPI_UNSIGNED_LONG_LONG unsigned 8-byte int
 *
 *  Assume: MPI_Datatype and nc_type are both enumerable types
 *          (this migh not conform wth MPI, as MPI_Datatype is intened to be
 *           an opaque data type.)
 *
 *  In OpenMPI, this assumption will fail
 */

MPI_Datatype
ncmpii_nc2mpitype(nc_type type)
{
    switch(type){
        case NC_BYTE :   return MPI_BYTE;
        case NC_CHAR :   return MPI_CHAR;
        case NC_SHORT :  return MPI_SHORT;
        case NC_INT :    return MPI_INT;
        case NC_FLOAT :  return MPI_FLOAT;
        case NC_DOUBLE : return MPI_DOUBLE;
        case NC_UBYTE :  return MPI_UNSIGNED_CHAR;
        case NC_USHORT : return MPI_UNSIGNED_SHORT;
        case NC_UINT :   return MPI_UNSIGNED;
        case NC_INT64 :  return MPI_LONG_LONG_INT;
        case NC_UINT64 : return MPI_UNSIGNED_LONG_LONG;
        default:         return MPI_DATATYPE_NULL;
    }
}

/*----< ncmpii_need_convert() >----------------------------------------------*/
inline int
ncmpii_need_convert(nc_type nctype,MPI_Datatype mpitype) {
    return !( (nctype == NC_CHAR   && mpitype == MPI_CHAR)           ||
              (nctype == NC_BYTE   && mpitype == MPI_BYTE)           ||
              (nctype == NC_SHORT  && mpitype == MPI_SHORT)          ||
              (nctype == NC_INT    && mpitype == MPI_INT)            ||
              (nctype == NC_INT    && mpitype == MPI_LONG &&
               X_SIZEOF_INT == SIZEOF_LONG)                          ||
              (nctype == NC_FLOAT  && mpitype == MPI_FLOAT)          ||
              (nctype == NC_DOUBLE && mpitype == MPI_DOUBLE)         ||
              (nctype == NC_UBYTE  && mpitype == MPI_UNSIGNED_CHAR)  ||
              (nctype == NC_USHORT && mpitype == MPI_UNSIGNED_SHORT) ||
              (nctype == NC_UINT   && mpitype == MPI_UNSIGNED)       ||
              (nctype == NC_INT64  && mpitype == MPI_LONG_LONG_INT)  ||
              (nctype == NC_UINT64 && mpitype == MPI_UNSIGNED_LONG_LONG)
            );
}

/*----< ncmpii_need_swap() >-------------------------------------------------*/
inline int 
ncmpii_need_swap(nc_type      nctype,
                 MPI_Datatype mpitype)
{
#ifdef WORDS_BIGENDIAN
    return 0;
#else
    if ((nctype == NC_CHAR   && mpitype == MPI_CHAR)           ||
        (nctype == NC_BYTE   && mpitype == MPI_BYTE)           ||
        (nctype == NC_UBYTE  && mpitype == MPI_UNSIGNED_CHAR))
        return 0;

    return 1;
#endif
}

#if 0
/*----< swapn() >------------------------------------------------------------*/
static void
swapn(void       *dst,
      const void *src,
      MPI_Offset  nn,
      int         xsize)
{
    int i;
    uchar *op = dst;
    const uchar *ip = src;
    while (nn-- != 0) {
        for (i=0; i<xsize; i++)
            op[i] = ip[xsize-1-i];
        op += xsize;
        ip += xsize;
    }
}
#endif

/* Endianness byte swap: done in-place */
#define SWAP(x,y) {tmp = (x); (x) = (y); (y) = tmp;}

/*----< ncmpii_swap() >-------------------------------------------------------*/
void
ncmpii_swapn(void       *dest_p,  /* destination array */
             const void *src_p,   /* source array */
             MPI_Offset  nelems,  /* number of elements in buf[] */
             int         esize)   /* byte size of each element */
{
    int  i;

    if (esize <= 1 || nelems <= 0) return;  /* no need */

    if (esize == 4) { /* this is the most common case */
              uint32_t *dest = (uint32_t*)       dest_p;
        const uint32_t *src  = (const uint32_t*) src_p;
        for (i=0; i<nelems; i++)
            dest[i] = htonl(src[i]);
    }
    else if (esize == 2) {
              uint16_t *dest =       (uint16_t*) dest_p;
        const uint16_t *src  = (const uint16_t*) src_p;
        for (i=0; i<nelems; i++)
            dest[i] = htons(src[i]);
    }
    else {
              uchar *op = (uchar*) dest_p;
        const uchar *ip = (uchar*) src_p;
        /* for esize is not 1, 2, or 4 */
        while (nelems-- > 0) {
            for (i=0; i<esize; i++)
                op[i] = ip[esize-1-i];
            op += esize;
            ip += esize;
        }
    }
}

/*----< ncmpii_in_swap() >---------------------------------------------------*/
void
ncmpii_in_swapn(void       *buf,
                MPI_Offset  nelems,  /* number of elements in buf[] */
                int         esize)   /* byte size of each element */
{
    int  i;
    uchar tmp, *op = (uchar*)buf;

    if (esize <= 1 || nelems <= 0) return;  /* no need */

    if (esize == 4) { /* this is the most common case */
        uint32_t *dest = (uint32_t*) buf;
        for (i=0; i<nelems; i++)
            dest[i] = htonl(dest[i]);
    }
    else if (esize == 2) {
        uint16_t *dest = (uint16_t*) buf;
        for (i=0; i<nelems; i++)
            dest[i] = htons(dest[i]);
    }
    else {
        /* for esize is not 1, 2, or 4 */
        while (nelems-- > 0) {
            for (i=0; i<esize/2; i++)
                SWAP(op[i], op[esize-1-i])
            op += esize;
        }
    }
}


/*----< ncmpii_x_putn_filetype() >-------------------------------------------*/
#define PUTN_NAME(ftype, btype) ncmpix_putn_##ftype##btype

#define X_PUTN_FILETYPE(ftype)                                                 \
int                                                                            \
ncmpii_x_putn_##ftype(void         *xp,      /* file buffer of type schar */   \
                      const void   *putbuf,  /* put buffer of type puttype */  \
                      MPI_Offset    nelems,                                    \
                      MPI_Datatype  puttype)                                   \
{                                                                              \
    if (puttype == MPI_CHAR || /* assume ECHAR has been checked before */      \
        puttype == MPI_BYTE)                                                   \
        return PUTN_NAME(ftype, _schar) (&xp, nelems, (const schar*)  putbuf); \
    else if (puttype == MPI_UNSIGNED_CHAR)                                     \
        return PUTN_NAME(ftype, _uchar) (&xp, nelems, (const uchar*)  putbuf); \
    else if (puttype == MPI_SHORT)                                             \
        return PUTN_NAME(ftype, _short) (&xp, nelems, (const short*)  putbuf); \
    else if (puttype == MPI_UNSIGNED_SHORT)                                    \
        return PUTN_NAME(ftype, _ushort)(&xp, nelems, (const ushort*) putbuf); \
    else if (puttype == MPI_INT)                                               \
        return PUTN_NAME(ftype, _int)   (&xp, nelems, (const int*)    putbuf); \
    else if (puttype == MPI_UNSIGNED)                                          \
        return PUTN_NAME(ftype, _uint)  (&xp, nelems, (const uint*)   putbuf); \
    else if (puttype == MPI_LONG)                                              \
        return PUTN_NAME(ftype, _long)  (&xp, nelems, (const long*)   putbuf); \
    else if (puttype == MPI_FLOAT)                                             \
        return PUTN_NAME(ftype, _float) (&xp, nelems, (const float*)  putbuf); \
    else if (puttype == MPI_DOUBLE)                                            \
        return PUTN_NAME(ftype, _double)(&xp, nelems, (const double*) putbuf); \
    else if (puttype == MPI_LONG_LONG_INT)                                     \
        return PUTN_NAME(ftype, _int64) (&xp, nelems, (const int64*)  putbuf); \
    else if (puttype == MPI_UNSIGNED_LONG_LONG)                                \
        return PUTN_NAME(ftype, _uint64)(&xp, nelems, (const uint64*) putbuf); \
    return NC_EBADTYPE;                                                        \
}

/*----< ncmpii_x_putn_schar() >----------------------------------------------*/
/*----< ncmpii_x_putn_uchar() >----------------------------------------------*/
/*----< ncmpii_x_putn_short() >----------------------------------------------*/
/*----< ncmpii_x_putn_ushort() >---------------------------------------------*/
/*----< ncmpii_x_putn_int() >------------------------------------------------*/
/*----< ncmpii_x_putn_uint() >-----------------------------------------------*/
/*----< ncmpii_x_putn_float() >----------------------------------------------*/
/*----< ncmpii_x_putn_double() >---------------------------------------------*/
/*----< ncmpii_x_putn_int64() >----------------------------------------------*/
/*----< ncmpii_x_putn_uint64() >---------------------------------------------*/
X_PUTN_FILETYPE(schar)
X_PUTN_FILETYPE(uchar)
X_PUTN_FILETYPE(short)
X_PUTN_FILETYPE(ushort)
X_PUTN_FILETYPE(int)
X_PUTN_FILETYPE(uint)
X_PUTN_FILETYPE(float)
X_PUTN_FILETYPE(double)
X_PUTN_FILETYPE(int64)
X_PUTN_FILETYPE(uint64)


/*----< ncmpii_x_getn_filetype() >-------------------------------------------*/
#define GETN_NAME(ftype, btype) ncmpix_getn_##ftype##btype

#define X_GETN_FILETYPE(ftype)                                                \
int                                                                           \
ncmpii_x_getn_##ftype(const void   *xp,      /* file buffer of type schar */  \
                      void         *getbuf,  /* get buffer of type gettype */ \
                      MPI_Offset    nelems,                                   \
                      MPI_Datatype  gettype)                                  \
{                                                                             \
    if (gettype == MPI_CHAR || /* assume ECHAR has been checked before */     \
        gettype == MPI_BYTE)                                                  \
        return GETN_NAME(ftype, _schar) (&xp, nelems, (schar*)  getbuf);      \
    else if (gettype == MPI_UNSIGNED_CHAR)                                    \
        return GETN_NAME(ftype, _uchar) (&xp, nelems, (uchar*)  getbuf);      \
    else if (gettype == MPI_SHORT)                                            \
        return GETN_NAME(ftype, _short) (&xp, nelems, (short*)  getbuf);      \
    else if (gettype == MPI_UNSIGNED_SHORT)                                   \
        return GETN_NAME(ftype, _ushort)(&xp, nelems, (ushort*) getbuf);      \
    else if (gettype == MPI_INT)                                              \
        return GETN_NAME(ftype, _int)   (&xp, nelems, (int*)    getbuf);      \
    else if (gettype == MPI_UNSIGNED)                                         \
        return GETN_NAME(ftype, _uint)  (&xp, nelems, (uint*)   getbuf);      \
    else if (gettype == MPI_LONG)                                             \
        return GETN_NAME(ftype, _long)  (&xp, nelems, (long*)   getbuf);      \
    else if (gettype == MPI_FLOAT)                                            \
        return GETN_NAME(ftype, _float) (&xp, nelems, (float*)  getbuf);      \
    else if (gettype == MPI_DOUBLE)                                           \
        return GETN_NAME(ftype, _double)(&xp, nelems, (double*) getbuf);      \
    else if (gettype == MPI_LONG_LONG_INT)                                    \
        return GETN_NAME(ftype, _int64) (&xp, nelems, (int64*)  getbuf);      \
    else if (gettype == MPI_UNSIGNED_LONG_LONG)                               \
        return GETN_NAME(ftype, _uint64)(&xp, nelems, (uint64*) getbuf);      \
    return NC_EBADTYPE;                                                       \
}
/*----< calling ncmpix_getn_schar_schar() >----------------------------------*/
/*----< calling ncmpix_getn_schar_uchar() >----------------------------------*/
/*----< calling ncmpix_getn_schar_short() >----------------------------------*/
/*----< calling ncmpix_getn_schar_ushort() >---------------------------------*/
/*----< calling ncmpix_getn_schar_int() >------------------------------------*/
/*----< calling ncmpix_getn_schar_uint() >-----------------------------------*/
/*----< calling ncmpix_getn_schar_float() >----------------------------------*/
/*----< calling ncmpix_getn_schar_double() >---------------------------------*/
/*----< calling ncmpix_getn_schar_int64() >----------------------------------*/
/*----< calling ncmpix_getn_schar_uint64() >---------------------------------*/

/*----< calling ncmpix_getn_uchar_schar() >----------------------------------*/
/*----< calling ncmpix_getn_uchar_uchar() >----------------------------------*/
/*----< calling ncmpix_getn_uchar_short() >----------------------------------*/
/*----< calling ncmpix_getn_uchar_ushort() >---------------------------------*/
/*----< calling ncmpix_getn_uchar_int() >------------------------------------*/
/*----< calling ncmpix_getn_uchar_uint() >-----------------------------------*/
/*----< calling ncmpix_getn_uchar_float() >----------------------------------*/
/*----< calling ncmpix_getn_uchar_double() >---------------------------------*/
/*----< calling ncmpix_getn_uchar_int64() >----------------------------------*/
/*----< calling ncmpix_getn_uchar_uint64() >---------------------------------*/

/*----< calling ncmpix_getn_short_schar() >----------------------------------*/
/*----< calling ncmpix_getn_short_uchar() >----------------------------------*/
/*----< calling ncmpix_getn_short_short() >----------------------------------*/
/*----< calling ncmpix_getn_short_ushort() >---------------------------------*/
/*----< calling ncmpix_getn_short_int() >------------------------------------*/
/*----< calling ncmpix_getn_short_uint() >-----------------------------------*/
/*----< calling ncmpix_getn_short_float() >----------------------------------*/
/*----< calling ncmpix_getn_short_double() >---------------------------------*/
/*----< calling ncmpix_getn_short_int64() >----------------------------------*/
/*----< calling ncmpix_getn_short_uint64() >---------------------------------*/

/*----< calling ncmpix_getn_int_schar() >------------------------------------*/
/*----< calling ncmpix_getn_int_uchar() >------------------------------------*/
/*----< calling ncmpix_getn_int_short() >------------------------------------*/
/*----< calling ncmpix_getn_int_ushort() >-----------------------------------*/
/*----< calling ncmpix_getn_int_int() >--------------------------------------*/
/*----< calling ncmpix_getn_int_uint() >-------------------------------------*/
/*----< calling ncmpix_getn_int_float() >------------------------------------*/
/*----< calling ncmpix_getn_int_double() >-----------------------------------*/
/*----< calling ncmpix_getn_int_int64() >------------------------------------*/
/*----< calling ncmpix_getn_int_uint64() >-----------------------------------*/

/*----< calling ncmpix_getn_float_schar() >----------------------------------*/
/*----< calling ncmpix_getn_float_uchar() >----------------------------------*/
/*----< calling ncmpix_getn_float_short() >----------------------------------*/
/*----< calling ncmpix_getn_float_ushort() >---------------------------------*/
/*----< calling ncmpix_getn_float_int() >------------------------------------*/
/*----< calling ncmpix_getn_float_uint() >-----------------------------------*/
/*----< calling ncmpix_getn_float_float() >----------------------------------*/
/*----< calling ncmpix_getn_float_double() >---------------------------------*/
/*----< calling ncmpix_getn_float_int64() >----------------------------------*/
/*----< calling ncmpix_getn_float_uint64() >---------------------------------*/

/*----< calling ncmpix_getn_double_schar() >---------------------------------*/
/*----< calling ncmpix_getn_double_uchar() >---------------------------------*/
/*----< calling ncmpix_getn_double_short() >---------------------------------*/
/*----< calling ncmpix_getn_double_ushort() >--------------------------------*/
/*----< calling ncmpix_getn_double_int() >-----------------------------------*/
/*----< calling ncmpix_getn_double_uint() >----------------------------------*/
/*----< calling ncmpix_getn_double_float() >---------------------------------*/
/*----< calling ncmpix_getn_double_double() >--------------------------------*/
/*----< calling ncmpix_getn_double_int64() >---------------------------------*/
/*----< calling ncmpix_getn_double_uint64() >--------------------------------*/

/*----< calling ncmpix_getn_ushort_schar() >---------------------------------*/
/*----< calling ncmpix_getn_ushort_uchar() >---------------------------------*/
/*----< calling ncmpix_getn_ushort_short() >---------------------------------*/
/*----< calling ncmpix_getn_ushort_ushort() >--------------------------------*/
/*----< calling ncmpix_getn_ushort_int() >-----------------------------------*/
/*----< calling ncmpix_getn_ushort_uint() >----------------------------------*/
/*----< calling ncmpix_getn_ushort_float() >---------------------------------*/
/*----< calling ncmpix_getn_ushort_double() >--------------------------------*/
/*----< calling ncmpix_getn_ushort_int64() >---------------------------------*/
/*----< calling ncmpix_getn_ushort_uint64() >--------------------------------*/

/*----< calling ncmpix_getn_uint_schar() >-----------------------------------*/
/*----< calling ncmpix_getn_uint_uchar() >-----------------------------------*/
/*----< calling ncmpix_getn_uint_short() >-----------------------------------*/
/*----< calling ncmpix_getn_uint_ushort() >----------------------------------*/
/*----< calling ncmpix_getn_uint_int() >-------------------------------------*/
/*----< calling ncmpix_getn_uint_uint() >------------------------------------*/
/*----< calling ncmpix_getn_uint_float() >-----------------------------------*/
/*----< calling ncmpix_getn_uint_double() >----------------------------------*/
/*----< calling ncmpix_getn_uint_int64() >-----------------------------------*/
/*----< calling ncmpix_getn_uint_uint64() >----------------------------------*/

/*----< calling ncmpix_getn_int64_schar() >----------------------------------*/
/*----< calling ncmpix_getn_int64_uchar() >----------------------------------*/
/*----< calling ncmpix_getn_int64_short() >----------------------------------*/
/*----< calling ncmpix_getn_int64_ushort() >---------------------------------*/
/*----< calling ncmpix_getn_int64_int() >------------------------------------*/
/*----< calling ncmpix_getn_int64_uint() >-----------------------------------*/
/*----< calling ncmpix_getn_int64_float() >----------------------------------*/
/*----< calling ncmpix_getn_int64_double() >---------------------------------*/
/*----< calling ncmpix_getn_int64_int64() >----------------------------------*/
/*----< calling ncmpix_getn_int64_uint64() >---------------------------------*/

/*----< calling ncmpix_getn_uint64_schar() >---------------------------------*/
/*----< calling ncmpix_getn_uint64_uchar() >---------------------------------*/
/*----< calling ncmpix_getn_uint64_short() >---------------------------------*/
/*----< calling ncmpix_getn_uint64_ushort() >--------------------------------*/
/*----< calling ncmpix_getn_uint64_int() >-----------------------------------*/
/*----< calling ncmpix_getn_uint64_uint() >----------------------------------*/
/*----< calling ncmpix_getn_uint64_float() >---------------------------------*/
/*----< calling ncmpix_getn_uint64_double() >--------------------------------*/
/*----< calling ncmpix_getn_uint64_int64() >---------------------------------*/
/*----< calling ncmpix_getn_uint64_uint64() >--------------------------------*/

/*----< ncmpii_x_getn_schar() >----------------------------------------------*/
/*----< ncmpii_x_getn_uchar() >----------------------------------------------*/
/*----< ncmpii_x_getn_short() >----------------------------------------------*/
/*----< ncmpii_x_getn_ushort() >---------------------------------------------*/
/*----< ncmpii_x_getn_int() >------------------------------------------------*/
/*----< ncmpii_x_getn_uint() >-----------------------------------------------*/
/*----< ncmpii_x_getn_float() >----------------------------------------------*/
/*----< ncmpii_x_getn_double() >---------------------------------------------*/
/*----< ncmpii_x_getn_int64() >----------------------------------------------*/
/*----< ncmpii_x_getn_uint64() >---------------------------------------------*/
X_GETN_FILETYPE(schar)
X_GETN_FILETYPE(uchar)
X_GETN_FILETYPE(short)
X_GETN_FILETYPE(ushort)
X_GETN_FILETYPE(int)
X_GETN_FILETYPE(uint)
X_GETN_FILETYPE(float)
X_GETN_FILETYPE(double)
X_GETN_FILETYPE(int64)
X_GETN_FILETYPE(uint64)

