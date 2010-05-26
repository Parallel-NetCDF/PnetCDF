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

#include "macro.h"


/* buffer layers:       
        
        User Level              buf     (user defined buffer of MPI_Datatype)
        MPI Datatype Level      cbuf    (contiguous buffer of ptype)
        NetCDF XDR Level        xbuf    (XDR I/O buffer)
*/

/*----< ncmpi_put_var1() >---------------------------------------------------*/
int
ncmpi_put_var1(int               ncid,
               int               varid,
               const MPI_Offset  start[],
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      datatype)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *count;

    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_INDEP_FH
    CHECK_VARID(varid, varp)
    GET_ONE_COUNT

    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,
                                (void*)buf, bufcount, datatype,
                                WRITE_REQ, INDEP_IO);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

#define PUT_VAR1_COMMON(datatype)                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset *count;                                           \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    CHECK_INDEP_FH                                               \
    CHECK_VARID(varid, varp)                                     \
    GET_ONE_COUNT                                                \
                                                                 \
    /* put_var1 is a special case of put_vars */                 \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                (void*)op, 1, datatype,          \
                                WRITE_REQ, INDEP_IO);            \
    if (varp->ndims > 0) NCI_Free(count);                        \
    return status;

/*----< ncmpi_put_var1_text() >----------------------------------------------*/
int
ncmpi_put_var1_text(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const char       *op)
{
    PUT_VAR1_COMMON(MPI_CHAR)
}

/*----< ncmpi_put_var1_schar() >---------------------------------------------*/
int
ncmpi_put_var1_schar(int                ncid,
                     int                varid,
                     const MPI_Offset   start[],
                     const signed char *op)
{
    PUT_VAR1_COMMON(MPI_BYTE)
}

/*----< ncmpi_put_var1_uchar() >---------------------------------------------*/
int
ncmpi_put_var1_uchar(int                  ncid,
                     int                  varid,
                     const MPI_Offset     start[],
                     const unsigned char *op)
{
    PUT_VAR1_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_put_var1_short() >---------------------------------------------*/
int
ncmpi_put_var1_short(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const short      *op)
{
    PUT_VAR1_COMMON(MPI_SHORT)
}

/*----< ncmpi_put_var1_int() >-----------------------------------------------*/
int
ncmpi_put_var1_int(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const int        *op)
{
    PUT_VAR1_COMMON(MPI_INT)
}

/*----< ncmpi_put_var1_long() >----------------------------------------------*/
int
ncmpi_put_var1_long(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const long       *op)
{
    PUT_VAR1_COMMON(MPI_LONG)
}

/*----< ncmpi_put_var1_float() >---------------------------------------------*/
int
ncmpi_put_var1_float(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const float      *op)
{
    PUT_VAR1_COMMON(MPI_FLOAT)
}

/*----< ncmpi_put_var1_double() >--------------------------------------------*/
int
ncmpi_put_var1_double(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const double     *op)
{
    PUT_VAR1_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_get_var1() >---------------------------------------------------*/
int
ncmpi_get_var1(int               ncid,
               int               varid,
               const MPI_Offset  start[],
               void             *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;
    MPI_Offset *count;

    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_INDEP_FH
    CHECK_VARID(varid, varp)
    GET_ONE_COUNT

    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,
                                buf, bufcount, datatype,
                                READ_REQ, INDEP_IO);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

#define GET_VAR1_COMMON(datatype)                                \
    int     status;                                              \
    NC     *ncp;                                                 \
    NC_var *varp;                                                \
    MPI_Offset *count;                                           \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    CHECK_INDEP_FH                                               \
    CHECK_VARID(varid, varp)                                     \
    GET_ONE_COUNT                                                \
                                                                 \
    /* get_var1 is a special case of get_vars */                 \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                ip, 1, datatype,                 \
                                READ_REQ, INDEP_IO);             \
    if (varp->ndims > 0) NCI_Free(count);                        \
    return status;

/*----< ncmpi_get_var1_text() >----------------------------------------------*/
int
ncmpi_get_var1_text(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    char             *ip)
{
    GET_VAR1_COMMON(MPI_CHAR)
}

/*----< ncmpi_get_var1_schar() >---------------------------------------------*/
int
ncmpi_get_var1_schar(int                ncid,
                     int                varid,
                     const MPI_Offset   start[],
                     signed char       *ip)
{
    GET_VAR1_COMMON(MPI_BYTE)
}

/*----< ncmpi_get_var1_uchar() >---------------------------------------------*/
int
ncmpi_get_var1_uchar(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     unsigned char    *ip)
{
    GET_VAR1_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_get_var1_short() >---------------------------------------------*/
int
ncmpi_get_var1_short(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     short            *ip)
{
    GET_VAR1_COMMON(MPI_SHORT)
}

/*----< ncmpi_get_var1_int() >-----------------------------------------------*/
int
ncmpi_get_var1_int(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   int              *ip)
{
    GET_VAR1_COMMON(MPI_INT)
}

/*----< ncmpi_get_var1_long() >----------------------------------------------*/
int
ncmpi_get_var1_long(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    long             *ip)
{
    GET_VAR1_COMMON(MPI_LONG)
}

/*----< ncmpi_get_var1_float() >---------------------------------------------*/
int
ncmpi_get_var1_float(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     float            *ip)
{
    GET_VAR1_COMMON(MPI_FLOAT)
}

/*----< ncmpi_get_var1_double() >--------------------------------------------*/
int
ncmpi_get_var1_double(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      double           *ip)
{
    GET_VAR1_COMMON(MPI_DOUBLE)
}

