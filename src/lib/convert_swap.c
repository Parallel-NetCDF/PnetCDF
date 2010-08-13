/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#include "nc.h"
#include "ncx.h"
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>
#include <arpa/inet.h>   /* htonl(), htons() */


/* Prototypes for functions used only in this file */
static void swapn(void *dst, const void *src, MPI_Offset nn, int xsize);

/*
 *  MAPPING:  MPI DATATYPE   <--->   NETCDF DATATYPE   DESCRIPTION
 *            MPI_UNSIGNED_CHAR      NC_BYTE           uchar integer
 *            MPI_BYTE               NC_BYTE           schar integer
 *            MPI_CHAR               NC_CHAR           char(text)
 *            MPI_SHORT              NC_SHORT          short
 *            MPI_INT                NC_INT            int
 *            MPI_FLOAT              NC_FLOAT          float
 *            MPI_DOUBLE             NC_DOUBLE         double
 *
 *
 *  Assume: MPI_Datatype and nc_type are both enumerable types
 */


/*----< ncmpii_echar() >-----------------------------------------------------*/
inline int
ncmpii_echar(nc_type      nctype,
             MPI_Datatype mpitype) {
    return ((nctype == NC_CHAR) == (mpitype != MPI_CHAR));
}

/*----< ncmpii_need_convert() >----------------------------------------------*/
inline int
ncmpii_need_convert(nc_type nctype,MPI_Datatype mpitype) {
    return !( (nctype == NC_CHAR && mpitype == MPI_CHAR) ||
              (nctype == NC_BYTE && mpitype == MPI_BYTE) ||
              (nctype == NC_BYTE && mpitype == MPI_UNSIGNED_CHAR) ||
              (nctype == NC_SHORT && mpitype == MPI_SHORT) ||
              (nctype == NC_INT && mpitype == MPI_INT) ||
              (nctype == NC_INT && mpitype == MPI_LONG &&
               X_SIZEOF_INT == SIZEOF_LONG) ||
              (nctype == NC_FLOAT && mpitype == MPI_FLOAT) ||
              (nctype == NC_DOUBLE && mpitype == MPI_DOUBLE) );
}

/*----< ncmpii_need_swap() >-------------------------------------------------*/
inline int 
ncmpii_need_swap(nc_type      nctype,
                 MPI_Datatype mpitype)
{
#ifdef WORDS_BIGENDIAN
    return 0;
#else
    return ( (nctype == NC_SHORT && mpitype == MPI_SHORT) ||
            (nctype == NC_INT && mpitype == MPI_INT) ||
            (nctype == NC_INT && mpitype == MPI_LONG &&
             X_SIZEOF_INT == SIZEOF_LONG) ||
            (nctype == NC_FLOAT && mpitype == MPI_FLOAT) ||
            (nctype == NC_DOUBLE && mpitype == MPI_DOUBLE) );
#endif
}

/*----< swapn() >------------------------------------------------------------*/
static void
swapn(void       *dst,
      const void *src,
      MPI_Offset  nn,
      int         xsize)
{
    int i;
    char *op = dst;
    const char *ip = src;
    while (nn-- != 0) {
        for (i=0; i<xsize; i++)
            op[i] = ip[xsize-1-i];
        op += xsize;
        ip += xsize;
    }
}


/* Endianness byte swap: done in-place */
#define SWAP(x,y) {tmp = (x); (x) = (y); (y) = tmp;}

/*----< ncmpii_in_swap() >---------------------------------------------------*/
void
ncmpii_in_swapn(void       *buf,
                MPI_Offset  nelems,  /* number of elements in buf[] */
                int         esize)   /* byte size of each element */
{
    int  i;
    char tmp, *op = buf;

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


/*----< ncmpii_x_putn_schar() >----------------------------------------------*/
int
ncmpii_x_putn_schar(void         *xbuf,
                    const void   *buf,
                    MPI_Offset    nelems,
                    MPI_Datatype  datatype)
{
    void *xp, *data;
    int status = NC_NOERR;

    xp = (void*) xbuf;
    data = (void*) buf;

    if (datatype == MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_SHORT)
        status = ncmpix_putn_schar_short(&xp, nelems, (const short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_putn_schar_int(&xp, nelems, (const int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_putn_schar_long(&xp, nelems, (const long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_putn_schar_float(&xp, nelems, (const float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_putn_schar_double(&xp, nelems, (const double*)data);

    return status;
}

/*----< ncmpii_x_putn_short() >----------------------------------------------*/
int
ncmpii_x_putn_short(void         *xbuf,
                    const void   *buf,
                    MPI_Offset    nelems,
                    MPI_Datatype  datatype)
{
    void *xp, *data;
    int status = NC_NOERR;
 
    xp = (void*) xbuf;
    data = (void*) buf;
 
    if (datatype == MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_UNSIGNED_CHAR)
        status = ncmpix_putn_short_uchar(&xp, nelems, (const unsigned char*)data);
    else if (datatype == MPI_BYTE)
        status = ncmpix_putn_short_schar(&xp, nelems, (const signed char*)data);
    else if (datatype == MPI_SHORT)
        status = ncmpix_putn_short_short(&xp, nelems, (const short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_putn_short_int(&xp, nelems, (const int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_putn_short_long(&xp, nelems, (const long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_putn_short_float(&xp, nelems, (const float*)data);
    else if (datatype == MPI_DOUBLE)
       status = ncmpix_putn_short_double(&xp, nelems, (const double*)data);

    return status;
} 

/*----< ncmpii_x_putn_int() >------------------------------------------------*/
int
ncmpii_x_putn_int(void         *xbuf,
                  const void   *buf,
                  MPI_Offset    nelems,
                  MPI_Datatype  datatype)
{
    void *xp, *data;
    int status = NC_NOERR;
 
    xp = (void*) xbuf;
    data = (void*) buf;
 
    if (datatype == MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_UNSIGNED_CHAR)
        status = ncmpix_putn_int_uchar(&xp, nelems, (const unsigned char*)data);
    else if (datatype == MPI_BYTE)
        status = ncmpix_putn_int_schar(&xp, nelems, (const signed char*)data);
    else if (datatype == MPI_SHORT)
        status = ncmpix_putn_int_short(&xp, nelems, (const short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_putn_int_int(&xp, nelems, (const int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_putn_int_long(&xp, nelems, (const long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_putn_int_float(&xp, nelems, (const float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_putn_int_double(&xp, nelems, (const double*)data);

    return status;
} 

/*----< ncmpii_x_putn_float() >----------------------------------------------*/
int
ncmpii_x_putn_float(void         *xbuf,
                    const void   *buf,
                    MPI_Offset    nelems,
                    MPI_Datatype  datatype)
{
    void *xp, *data;
    int status = NC_NOERR;
 
    xp = (void*) xbuf;
    data = (void*) buf;
 
    if (datatype == MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_UNSIGNED_CHAR)
        status = ncmpix_putn_float_uchar(&xp, nelems, (const unsigned char*)data);
    else if (datatype == MPI_BYTE)
        status = ncmpix_putn_float_schar(&xp, nelems, (const signed char*)data);
    else if (datatype == MPI_SHORT)
        status = ncmpix_putn_float_short(&xp, nelems, (const short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_putn_float_int(&xp, nelems, (const int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_putn_float_long(&xp, nelems, (const long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_putn_float_float(&xp, nelems, (const float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_putn_float_double(&xp, nelems, (const double*)data);

    return status;
} 

/*----< ncmpii_x_putn_double() >---------------------------------------------*/
int
ncmpii_x_putn_double(void         *xbuf,
                     const void   *buf,
                     MPI_Offset    nelems,
                     MPI_Datatype  datatype)
{
    void *xp, *data;
    int status = NC_NOERR;
 
    xp = (void*) xbuf; 
    data = (void*) buf;

    if (datatype ==MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_UNSIGNED_CHAR)
        status = ncmpix_putn_double_uchar(&xp, nelems, (const unsigned char*)data);
    else if (datatype == MPI_BYTE)
        status = ncmpix_putn_double_schar(&xp, nelems, (const signed char*)data);
    else if (datatype == MPI_SHORT)
        status = ncmpix_putn_double_short(&xp, nelems, (const short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_putn_double_int(&xp, nelems, (const int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_putn_double_long(&xp, nelems, (const long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_putn_double_float(&xp, nelems, (const float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_putn_double_double(&xp, nelems, (const double*)data);

    return status;
} 

/*----< ncmpii_x_getn_schar() >----------------------------------------------*/
int
ncmpii_x_getn_schar(const void   *xbuf,
                    void         *buf,
                    MPI_Offset    nelems,
                    MPI_Datatype  datatype)
{
    void *xp, *data;
    int status = NC_NOERR;

    xp = (void*) xbuf;
    data = (void*) buf;

    if (datatype == MPI_CHAR)
          status = NC_ECHAR;
    else if (datatype == MPI_SHORT)
        status = ncmpix_getn_schar_short((const void**)&xp, nelems, (short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_getn_schar_int((const void**)&xp, nelems, (int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_getn_schar_long((const void**)&xp, nelems, (long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_getn_schar_float((const void**)&xp, nelems, (float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_getn_schar_double((const void**)&xp, nelems, (double*)data);

    return status;
}

/*----< ncmpii_x_getn_short() >----------------------------------------------*/
int
ncmpii_x_getn_short(const void   *xbuf,
                    void         *buf,
                    MPI_Offset    nelems,
                    MPI_Datatype  datatype)
{
    char *xp, *data;
    int status = NC_NOERR;

    xp = (char*) xbuf;
    data = (char*) buf;
 
    if (datatype == MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_UNSIGNED_CHAR)
        status = ncmpix_getn_short_uchar((const void**)&xp, nelems, (unsigned char*)data);
    else if (datatype == MPI_BYTE)
        status = ncmpix_getn_short_schar((const void**)&xp, nelems, (signed char*)data);
    else if (datatype == MPI_SHORT)
        status = ncmpix_getn_short_short((const void**)&xp, nelems, (short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_getn_short_int((const void**)&xp, nelems, (int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_getn_short_long((const void**)&xp, nelems, (long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_getn_short_float((const void**)&xp, nelems, (float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_getn_short_double((const void**)&xp, nelems, (double*)data);

    return status;
} 

/*----< ncmpii_x_getn_int() >------------------------------------------------*/
int 
ncmpii_x_getn_int(const void   *xbuf,
                  void         *buf,
                  MPI_Offset    nelems,
                  MPI_Datatype  datatype)
{
    char *xp, *data;
    int status = NC_NOERR;
 
    xp = (char*) xbuf;
    data = (char*) buf;
 
    if (datatype == MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_UNSIGNED_CHAR)
        status = ncmpix_getn_int_uchar((const void**)&xp, nelems, (unsigned char*)data);
    else if (datatype == MPI_BYTE)
        status = ncmpix_getn_int_schar((const void**)&xp, nelems, (signed char*)data);
    else if (datatype == MPI_SHORT)
        status = ncmpix_getn_int_short((const void**)&xp, nelems, (short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_getn_int_int((const void**)&xp, nelems, (int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_getn_int_long((const void**)&xp, nelems, (long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_getn_int_float((const void**)&xp, nelems, (float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_getn_int_double((const void**)&xp, nelems, (double*)data);

    return status;
} 

/*----< ncmpii_x_getn_float() >----------------------------------------------*/
int
ncmpii_x_getn_float(const void   *xbuf,
                    void         *buf,
                    MPI_Offset    nelems,
                    MPI_Datatype  datatype)
{
    char *xp, *data;
    int  status = NC_NOERR;

    xp = (char*) xbuf;
    data = (char*) buf;
 
    if (datatype == MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_UNSIGNED_CHAR)
        status = ncmpix_getn_float_uchar((const void**)&xp, nelems, (unsigned char*)data);
    else if (datatype == MPI_BYTE)
        status = ncmpix_getn_float_schar((const void**)&xp, nelems, (signed char*)data);
    else if (datatype == MPI_SHORT)
        status = ncmpix_getn_float_short((const void**)&xp, nelems, (short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_getn_float_int((const void**)&xp, nelems, (int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_getn_float_long((const void**)&xp, nelems, (long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_getn_float_float((const void**)&xp, nelems, (float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_getn_float_double((const void**)&xp, nelems, (double*)data);

    return status;
}

/*----< ncmpii_x_getn_double() >---------------------------------------------*/
int
ncmpii_x_getn_double(const void   *xbuf,
                     void         *buf,
                     MPI_Offset   nelems,
                     MPI_Datatype datatype)
{
    char *xp, *data;
    int  status = NC_NOERR;
 
    xp = (char*) xbuf;
    data = (char*) buf;

    if (datatype == MPI_CHAR)
        status = NC_ECHAR;
    else if (datatype == MPI_UNSIGNED_CHAR)
        status = ncmpix_getn_double_uchar((const void**)&xp, nelems, (unsigned char*)data);
    else if (datatype == MPI_BYTE)
        status = ncmpix_getn_double_schar((const void**)&xp, nelems, (signed char*)data);
    else if (datatype == MPI_SHORT)
        status = ncmpix_getn_double_short((const void**)&xp, nelems, (short*)data);
    else if (datatype == MPI_INT)
        status = ncmpix_getn_double_int((const void**)&xp, nelems, (int*)data);
    else if (datatype == MPI_LONG)
        status = ncmpix_getn_double_long((const void**)&xp, nelems, (long*)data);
    else if (datatype == MPI_FLOAT)
        status = ncmpix_getn_double_float((const void**)&xp, nelems, (float*)data);
    else if (datatype == MPI_DOUBLE)
        status = ncmpix_getn_double_double((const void**)&xp, nelems, (double*)data);

    return status;
}

