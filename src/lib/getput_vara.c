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

/*----< ncmpi_put_vara() >---------------------------------------------------*/
int
ncmpi_put_vara(int               ncid,
               int               varid,
               const MPI_Offset  start[],
               const MPI_Offset  count[],
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_INDEP_FH
    CHECK_VARID(varid, varp)

    return ncmpii_getput_vars(ncp, varp, start, count, NULL,
                              (void*)buf, bufcount, datatype,
                              WRITE_REQ, INDEP_IO);
}

#define PUT_VARA_COMMON(datatype)                              \
    int         status;                                        \
    NC         *ncp;                                           \
    NC_var     *varp;                                          \
    MPI_Offset  nelems;                                        \
                                                               \
    CHECK_NCID                                                 \
    CHECK_WRITE_PERMISSION                                     \
    if (NC_indef(ncp)) return NC_EINDEFINE;                    \
    CHECK_INDEP_FH                                             \
    CHECK_VARID(varid, varp)                                   \
    GET_NUM_ELEMENTS                                           \
                                                               \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                              (void*)op, nelems, datatype,     \
                              WRITE_REQ, INDEP_IO);

/*----< ncmpi_put_vara_text() >----------------------------------------------*/
int
ncmpi_put_vara_text(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const char       *op)
{
    PUT_VARA_COMMON(MPI_CHAR)
}

/*----< ncmpi_put_vara_schar() >---------------------------------------------*/
int
ncmpi_put_vara_schar(int                ncid,
                     int                varid,
                     const MPI_Offset   start[],
                     const MPI_Offset   count[],
                     const signed char *op)
{
    PUT_VARA_COMMON(MPI_BYTE)
}

/*----< ncmpi_put_vara_uchar() >---------------------------------------------*/
int
ncmpi_put_vara_uchar(int                  ncid,
                     int                  varid,
                     const MPI_Offset     start[],
                     const MPI_Offset     count[],
                     const unsigned char *op)
{
    PUT_VARA_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_put_vara_short() >---------------------------------------------*/
int
ncmpi_put_vara_short(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const short      *op)
{
    PUT_VARA_COMMON(MPI_SHORT)
}

/*----< ncmpi_put_vara_int() >-----------------------------------------------*/
int
ncmpi_put_vara_int(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   const int        *op)
{
    PUT_VARA_COMMON(MPI_INT)
}

/*----< ncmpi_put_vara_long() >----------------------------------------------*/
int
ncmpi_put_vara_long(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const long       *op)
{
    PUT_VARA_COMMON(MPI_LONG)
}

/*----< ncmpi_put_vara_float() >---------------------------------------------*/
int
ncmpi_put_vara_float(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const float      *op)
{
    PUT_VARA_COMMON(MPI_FLOAT)
}

/*----< ncmpi_put_vara_double() >--------------------------------------------*/
int
ncmpi_put_vara_double(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const double     *op)
{
    PUT_VARA_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_get_vara() >---------------------------------------------------*/
int
ncmpi_get_vara(int               ncid,
               int               varid,
               const MPI_Offset  start[],
               const MPI_Offset  count[],
               void             *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_INDEP_FH
    CHECK_VARID(varid, varp)

    return ncmpii_getput_vars(ncp, varp, start, count, NULL,
                              buf, bufcount, datatype,
                              READ_REQ, INDEP_IO);
}

#define GET_VARA_COMMON(datatype)                              \
    int         status;                                        \
    NC         *ncp;                                           \
    NC_var     *varp;                                          \
    MPI_Offset  nelems;                                        \
                                                               \
    CHECK_NCID                                                 \
    if (NC_indef(ncp)) return NC_EINDEFINE;                    \
    CHECK_INDEP_FH                                             \
    CHECK_VARID(varid, varp)                                   \
    GET_NUM_ELEMENTS                                           \
                                                               \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                              ip, nelems, datatype,            \
                              READ_REQ, INDEP_IO);

/*----< ncmpi_get_vara_text() >----------------------------------------------*/
int
ncmpi_get_vara_text(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    char             *ip)
{
    GET_VARA_COMMON(MPI_CHAR)
}

/*----< ncmpi_get_vara_schar() >---------------------------------------------*/
int
ncmpi_get_vara_schar(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     signed char      *ip)
{
    GET_VARA_COMMON(MPI_BYTE)
}

/*----< ncmpi_get_vara_uchar() >---------------------------------------------*/
int
ncmpi_get_vara_uchar(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     unsigned char    *ip)
{
    GET_VARA_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_get_vara_short() >---------------------------------------------*/
int
ncmpi_get_vara_short(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     short            *ip)
{
    GET_VARA_COMMON(MPI_SHORT)
}

/*----< ncmpi_get_vara_int() >-----------------------------------------------*/
int
ncmpi_get_vara_int(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   int              *ip)
{
    GET_VARA_COMMON(MPI_INT)
}

/*----< ncmpi_get_vara_long() >----------------------------------------------*/
int
ncmpi_get_vara_long(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    long             *ip)
{
    GET_VARA_COMMON(MPI_LONG)
}

/*----< ncmpi_get_vara_float() >---------------------------------------------*/
int
ncmpi_get_vara_float(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     float            *ip)
{
    GET_VARA_COMMON(MPI_FLOAT)
}

/*----< ncmpi_get_vara_double() >--------------------------------------------*/
int
ncmpi_get_vara_double(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      double           *ip)
{
    GET_VARA_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_put_vara_all() >-----------------------------------------------*/
int
ncmpi_put_vara_all(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   const void       *buf,
                   MPI_Offset        bufcount,
                   MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_COLLECTIVE_FH
    CHECK_VARID(varid, varp)

    return ncmpii_getput_vars(ncp, varp, start, count, NULL,
                              (void*)buf, bufcount, datatype,
                              WRITE_REQ, COLL_IO);
}

#define PUT_VARA_ALL_COMMON(datatype)                           \
    int         status;                                         \
    NC         *ncp;                                            \
    NC_var     *varp;                                           \
    MPI_Offset  nelems;                                         \
                                                                \
    CHECK_NCID                                                  \
    CHECK_WRITE_PERMISSION                                      \
    if (NC_indef(ncp)) return NC_EINDEFINE;                     \
    CHECK_COLLECTIVE_FH                                         \
    CHECK_VARID(varid, varp)                                    \
    GET_NUM_ELEMENTS                                            \
                                                                \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,    \
                              (void*)op, nelems, datatype,      \
                              WRITE_REQ, COLL_IO);

/*----< ncmpi_put_vara_text_all() >------------------------------------------*/
int
ncmpi_put_vara_text_all(int               ncid,
                        int               varid,
                        const MPI_Offset  start[],
                        const MPI_Offset  count[],
                        const char       *op)
{
    PUT_VARA_ALL_COMMON(MPI_CHAR)
}

/*----< ncmpi_put_vara_schar_all() >-----------------------------------------*/
int
ncmpi_put_vara_schar_all(int                ncid,
                         int                varid,
                         const MPI_Offset   start[],
                         const MPI_Offset   count[],
                         const signed char *op)
{
    PUT_VARA_ALL_COMMON(MPI_BYTE)
}

/*----< ncmpi_put_vara_uchar_all() >-----------------------------------------*/
int
ncmpi_put_vara_uchar_all(int                  ncid,
                         int                  varid,
                         const MPI_Offset     start[],
                         const MPI_Offset     count[],
                         const unsigned char *op)
{
    PUT_VARA_ALL_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_put_vara_short_all() >-----------------------------------------*/
int
ncmpi_put_vara_short_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         const short      *op)
{
    PUT_VARA_ALL_COMMON(MPI_SHORT)
}

/*----< ncmpi_put_vara_int_all() >-------------------------------------------*/
int
ncmpi_put_vara_int_all(int               ncid,
                       int               varid,
                       const MPI_Offset  start[],
                       const MPI_Offset  count[],
                       const int        *op)
{
    PUT_VARA_ALL_COMMON(MPI_INT)
}

/*----< ncmpi_put_vara_long_all() >------------------------------------------*/
int
ncmpi_put_vara_long_all(int               ncid,
                        int               varid,
                        const MPI_Offset  start[],
                        const MPI_Offset  count[],
                        const long       *op)
{
    PUT_VARA_ALL_COMMON(MPI_LONG)
}

/*----< ncmpi_put_vara_float_all() >-----------------------------------------*/
int
ncmpi_put_vara_float_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         const float      *op)
{
    PUT_VARA_ALL_COMMON(MPI_FLOAT)
}

/*----< ncmpi_put_vara_double_all() >----------------------------------------*/
int
ncmpi_put_vara_double_all(int               ncid,
                          int               varid,
                          const MPI_Offset  start[],
                          const MPI_Offset  count[],
                          const double     *op)
{
    PUT_VARA_ALL_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_get_vara_all() >-----------------------------------------------*/
int
ncmpi_get_vara_all(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   void             *buf,
                   MPI_Offset        bufcount,
                   MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_COLLECTIVE_FH
    CHECK_VARID(varid, varp)

    return ncmpii_getput_vars(ncp, varp, start, count, NULL,
                              buf, bufcount, datatype,
                              READ_REQ, COLL_IO);
}

#define GET_VARA_ALL_COMMON(datatype)                          \
    int         status;                                        \
    NC         *ncp;                                           \
    NC_var     *varp;                                          \
    MPI_Offset  nelems;                                        \
                                                               \
    CHECK_NCID                                                 \
    if (NC_indef(ncp)) return NC_EINDEFINE;                    \
    CHECK_COLLECTIVE_FH                                        \
    CHECK_VARID(varid, varp)                                   \
    GET_NUM_ELEMENTS                                           \
                                                               \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                              ip, nelems, datatype,            \
                              READ_REQ, COLL_IO);

/*----< ncmpi_get_vara_text_all() >------------------------------------------*/
int
ncmpi_get_vara_text_all(int               ncid,
                        int               varid,
                        const MPI_Offset  start[],
                        const MPI_Offset  count[],
                        char             *ip)
{
    GET_VARA_ALL_COMMON(MPI_CHAR)
}

/*----< ncmpi_get_vara_schar_all() >-----------------------------------------*/
int
ncmpi_get_vara_schar_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         signed char      *ip)
{
    GET_VARA_ALL_COMMON(MPI_BYTE)
}

/*----< ncmpi_get_vara_uchar_all() >-----------------------------------------*/
int
ncmpi_get_vara_uchar_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         unsigned char    *ip)
{
    GET_VARA_ALL_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_get_vara_short_all() >-----------------------------------------*/
int
ncmpi_get_vara_short_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         short            *ip)
{
    GET_VARA_ALL_COMMON(MPI_SHORT)
}

/*----< ncmpi_get_vara_int_all() >-------------------------------------------*/
int
ncmpi_get_vara_int_all(int               ncid,
                       int               varid,
                       const MPI_Offset  start[],
                       const MPI_Offset  count[],
                       int              *ip)
{
    GET_VARA_ALL_COMMON(MPI_INT)
}

/*----< ncmpi_get_vara_long_all() >------------------------------------------*/
int
ncmpi_get_vara_long_all(int               ncid,
                        int               varid,
                        const MPI_Offset  start[],
                        const MPI_Offset  count[],
                        long             *ip)
{
    GET_VARA_ALL_COMMON(MPI_LONG)
}

/*----< ncmpi_get_vara_float_all() >-----------------------------------------*/
int
ncmpi_get_vara_float_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         float            *ip)
{
    GET_VARA_ALL_COMMON(MPI_FLOAT)
}

/*----< ncmpi_get_vara_double_all() >----------------------------------------*/
int
ncmpi_get_vara_double_all(int               ncid,
                          int               varid,
                          const MPI_Offset  start[],
                          const MPI_Offset  count[],
                          double           *ip)
{
    GET_VARA_ALL_COMMON(MPI_DOUBLE)
}

