/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

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


/* ftype is the variable's nc_type defined in file, eg. int64
 * btype is the I/O buffer's C data type, eg. long long
 * buftype is I/O bufer's MPI data type, eg. MPI_UNSIGNED_LONG_LONG
 * apitype is data type appeared in the API names, eg. ncmpi_get_vara_longlong
 */

#define PUT_VARA(iomode, collmode)                               \
int                                                              \
ncmpi_put_vara##iomode(int               ncid,                   \
                       int               varid,                  \
                       const MPI_Offset  start[],                \
                       const MPI_Offset  count[],                \
                       const void       *buf,                    \
                       MPI_Offset        bufcount,               \
                       MPI_Datatype      buftype)                \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
                                                                 \
    SANITY_CHECK(WRITE_REQ, collmode)                            \
                                                                 \
    /* put_vara is a special case of put_vars */                 \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,     \
                              (void*)buf, bufcount, buftype,     \
                              WRITE_REQ, collmode);              \
}

/*----< ncmpi_put_vara() >---------------------------------------------------*/
/*----< ncmpi_put_vara_all() >-----------------------------------------------*/
PUT_VARA(    , INDEP_IO)
PUT_VARA(_all, COLL_IO)

#define GET_VARA(iomode, collmode)                               \
int                                                              \
ncmpi_get_vara##iomode(int               ncid,                   \
                       int               varid,                  \
                       const MPI_Offset  start[],                \
                       const MPI_Offset  count[],                \
                       void             *buf,                    \
                       MPI_Offset        bufcount,               \
                       MPI_Datatype      buftype)                \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
                                                                 \
    SANITY_CHECK(READ_REQ, collmode)                             \
                                                                 \
    /* get_vara is a special case of get_vars */                 \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,     \
                              buf, bufcount, buftype,            \
                              READ_REQ, collmode);               \
}

/*----< ncmpi_get_vara() >---------------------------------------------------*/
/*----< ncmpi_get_vara_all() >-----------------------------------------------*/
GET_VARA(    , INDEP_IO)
GET_VARA(_all, COLL_IO)

#define PUT_VARA_TYPE(apitype, btype, mpitype, collmode)         \
int                                                              \
ncmpi_put_vara_##apitype(int               ncid,                 \
                         int               varid,                \
                         const MPI_Offset  start[],              \
                         const MPI_Offset  count[],              \
                         const btype      *op)                   \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems;                                          \
                                                                 \
    SANITY_CHECK(WRITE_REQ, collmode)                            \
                                                                 \
    GET_NUM_ELEMENTS                                             \
                                                                 \
    /* put_vara is a special case of put_vars */                 \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,     \
                              (void*)op, nelems, mpitype,        \
                              WRITE_REQ, collmode);              \
}

/*----< ncmpi_put_vara_text() >-----------------------------------------------*/
/*----< ncmpi_put_vara_schar() >----------------------------------------------*/
/*----< ncmpi_put_vara_uchar() >----------------------------------------------*/
/*----< ncmpi_put_vara_short() >----------------------------------------------*/
/*----< ncmpi_put_vara_ushort() >---------------------------------------------*/
/*----< ncmpi_put_vara_int() >------------------------------------------------*/
/*----< ncmpi_put_vara_uint() >-----------------------------------------------*/
/*----< ncmpi_put_vara_long() >-----------------------------------------------*/
/*----< ncmpi_put_vara_float() >----------------------------------------------*/
/*----< ncmpi_put_vara_double() >---------------------------------------------*/
/*----< ncmpi_put_vara_longlong() >-------------------------------------------*/
/*----< ncmpi_put_vara_ulonglong() >------------------------------------------*/
PUT_VARA_TYPE(text,      char,               MPI_CHAR,               INDEP_IO)
PUT_VARA_TYPE(schar,     schar,              MPI_BYTE,               INDEP_IO)
PUT_VARA_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR,      INDEP_IO)
PUT_VARA_TYPE(short,     short,              MPI_SHORT,              INDEP_IO)
PUT_VARA_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT,     INDEP_IO)
PUT_VARA_TYPE(int,       int,                MPI_INT,                INDEP_IO)
PUT_VARA_TYPE(uint,      uint,               MPI_UNSIGNED,           INDEP_IO)
PUT_VARA_TYPE(long,      long,               MPI_LONG,               INDEP_IO)
PUT_VARA_TYPE(float,     float,              MPI_FLOAT,              INDEP_IO)
PUT_VARA_TYPE(double,    double,             MPI_DOUBLE,             INDEP_IO)
PUT_VARA_TYPE(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
PUT_VARA_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)
// PUT_VARA_TYPE(string, char*,              MPI_CHAR,               INDEP_IO)
/* string is not yet supported */

/*----< ncmpi_put_vara_text_all() >-------------------------------------------*/
/*----< ncmpi_put_vara_schar_all() >------------------------------------------*/
/*----< ncmpi_put_vara_uchar_all() >------------------------------------------*/
/*----< ncmpi_put_vara_short_all() >------------------------------------------*/
/*----< ncmpi_put_vara_ushort_all() >-----------------------------------------*/
/*----< ncmpi_put_vara_int_all() >--------------------------------------------*/
/*----< ncmpi_put_vara_uint_all() >-------------------------------------------*/
/*----< ncmpi_put_vara_long_all() >-------------------------------------------*/
/*----< ncmpi_put_vara_float_all() >------------------------------------------*/
/*----< ncmpi_put_vara_double_all() >-----------------------------------------*/
/*----< ncmpi_put_vara_longlong_all() >---------------------------------------*/
/*----< ncmpi_put_vara_ulonglong_all() >--------------------------------------*/
PUT_VARA_TYPE(text_all,      char,               MPI_CHAR,              COLL_IO)
PUT_VARA_TYPE(schar_all,     schar,              MPI_BYTE,              COLL_IO)
PUT_VARA_TYPE(uchar_all,     uchar,              MPI_UNSIGNED_CHAR,     COLL_IO)
PUT_VARA_TYPE(short_all,     short,              MPI_SHORT,             COLL_IO)
PUT_VARA_TYPE(ushort_all,    ushort,             MPI_UNSIGNED_SHORT,    COLL_IO)
PUT_VARA_TYPE(int_all,       int,                MPI_INT,               COLL_IO)
PUT_VARA_TYPE(uint_all,      uint,               MPI_UNSIGNED,          COLL_IO)
PUT_VARA_TYPE(long_all,      long,               MPI_LONG,              COLL_IO)
PUT_VARA_TYPE(float_all,     float,              MPI_FLOAT,             COLL_IO)
PUT_VARA_TYPE(double_all,    double,             MPI_DOUBLE,            COLL_IO)
PUT_VARA_TYPE(longlong_all,  long long,          MPI_LONG_LONG_INT,     COLL_IO)
PUT_VARA_TYPE(ulonglong_all, unsigned long long, MPI_UNSIGNED_LONG_LONG,COLL_IO)
// PUT_VARA_TYPE(string_all, char*,              MPI_CHAR,              COLL_IO)
/* string is not yet supported */


#define GET_VARA_TYPE(apitype, btype, mpitype, collmode)         \
int                                                              \
ncmpi_get_vara_##apitype(int               ncid,                 \
                         int               varid,                \
                         const MPI_Offset  start[],              \
                         const MPI_Offset  count[],              \
                         btype            *ip)                   \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems;                                          \
                                                                 \
    SANITY_CHECK(READ_REQ, collmode)                             \
                                                                 \
    GET_NUM_ELEMENTS                                             \
                                                                 \
    /* get_vara is a special case of get_vars */                 \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,     \
                              ip, nelems, mpitype,               \
                              READ_REQ, collmode);               \
}

/*----< ncmpi_get_vara_text() >-----------------------------------------------*/
/*----< ncmpi_get_vara_schar() >----------------------------------------------*/
/*----< ncmpi_get_vara_uchar() >----------------------------------------------*/
/*----< ncmpi_get_vara_short() >----------------------------------------------*/
/*----< ncmpi_get_vara_ushort() >---------------------------------------------*/
/*----< ncmpi_get_vara_int() >------------------------------------------------*/
/*----< ncmpi_get_vara_uint() >-----------------------------------------------*/
/*----< ncmpi_get_vara_long() >-----------------------------------------------*/
/*----< ncmpi_get_vara_float() >----------------------------------------------*/
/*----< ncmpi_get_vara_double() >---------------------------------------------*/
/*----< ncmpi_get_vara_longlong() >-------------------------------------------*/
/*----< ncmpi_get_vara_ulonglong() >------------------------------------------*/
GET_VARA_TYPE(text,      char,               MPI_CHAR,               INDEP_IO)
GET_VARA_TYPE(schar,     schar,              MPI_BYTE,               INDEP_IO)
GET_VARA_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR,      INDEP_IO)
GET_VARA_TYPE(short,     short,              MPI_SHORT,              INDEP_IO)
GET_VARA_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT,     INDEP_IO)
GET_VARA_TYPE(int,       int,                MPI_INT,                INDEP_IO)
GET_VARA_TYPE(uint,      uint,               MPI_UNSIGNED,           INDEP_IO)
GET_VARA_TYPE(long,      long,               MPI_LONG,               INDEP_IO)
GET_VARA_TYPE(float,     float,              MPI_FLOAT,              INDEP_IO)
GET_VARA_TYPE(double,    double,             MPI_DOUBLE,             INDEP_IO)
GET_VARA_TYPE(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
GET_VARA_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)
// GET_VARA_TYPE(string, char*,              MPI_CHAR,               INDEP_IO)
/* string is not yet supported */

/*----< ncmpi_get_vara_text_all() >-------------------------------------------*/
/*----< ncmpi_get_vara_schar_all() >------------------------------------------*/
/*----< ncmpi_get_vara_uchar_all() >------------------------------------------*/
/*----< ncmpi_get_vara_short_all() >------------------------------------------*/
/*----< ncmpi_get_vara_ushort_all() >-----------------------------------------*/
/*----< ncmpi_get_vara_int_all() >--------------------------------------------*/
/*----< ncmpi_get_vara_uint_all() >-------------------------------------------*/
/*----< ncmpi_get_vara_long_all() >-------------------------------------------*/
/*----< ncmpi_get_vara_float_all() >------------------------------------------*/
/*----< ncmpi_get_vara_double_all() >-----------------------------------------*/
/*----< ncmpi_get_vara_longlong_all() >---------------------------------------*/
/*----< ncmpi_get_vara_ulonglong_all() >--------------------------------------*/
GET_VARA_TYPE(text_all,      char,               MPI_CHAR,              COLL_IO)
GET_VARA_TYPE(schar_all,     schar,              MPI_BYTE,              COLL_IO)
GET_VARA_TYPE(uchar_all,     uchar,              MPI_UNSIGNED_CHAR,     COLL_IO)
GET_VARA_TYPE(short_all,     short,              MPI_SHORT,             COLL_IO)
GET_VARA_TYPE(ushort_all,    ushort,             MPI_UNSIGNED_SHORT,    COLL_IO)
GET_VARA_TYPE(int_all,       int,                MPI_INT,               COLL_IO)
GET_VARA_TYPE(uint_all,      uint,               MPI_UNSIGNED,          COLL_IO)
GET_VARA_TYPE(long_all,      long,               MPI_LONG,              COLL_IO)
GET_VARA_TYPE(float_all,     float,              MPI_FLOAT,             COLL_IO)
GET_VARA_TYPE(double_all,    double,             MPI_DOUBLE,            COLL_IO)
GET_VARA_TYPE(longlong_all,  long long,          MPI_LONG_LONG_INT,     COLL_IO)
GET_VARA_TYPE(ulonglong_all, unsigned long long, MPI_UNSIGNED_LONG_LONG,COLL_IO)
// GET_VARA_TYPE(string_all, char*,              MPI_CHAR,              COLL_IO)
/* string is not yet supported */

