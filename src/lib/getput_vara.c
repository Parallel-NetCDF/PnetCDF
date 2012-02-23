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

#define PUT_VARA(fnmode, collmode)                               \
int                                                              \
ncmpi_put_vara##fnmode(int               ncid,                   \
                       int               varid,                  \
                       const MPI_Offset  start[],                \
                       const MPI_Offset  count[],                \
                       const void       *buf,                    \
                       MPI_Offset        bufcount,               \
                       MPI_Datatype      datatype)               \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
                                                                 \
    /* put_vara is a special case of put_vars */                 \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,     \
                              (void*)buf, bufcount, datatype,    \
                              WRITE_REQ, collmode);              \
}

/*----< ncmpi_put_vara() >---------------------------------------------------*/
/*----< ncmpi_put_vara_all() >-----------------------------------------------*/
PUT_VARA(    , INDEP_IO)
PUT_VARA(_all, COLL_IO)

#define GET_VARA(fnmode, collmode)                               \
int                                                              \
ncmpi_get_vara##fnmode(int               ncid,                   \
                       int               varid,                  \
                       const MPI_Offset  start[],                \
                       const MPI_Offset  count[],                \
                       void             *buf,                    \
                       MPI_Offset        bufcount,               \
                       MPI_Datatype      datatype)               \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
                                                                 \
    /* get_vara is a special case of get_vars */                 \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,     \
                              buf, bufcount, datatype,           \
                              READ_REQ, collmode);               \
}

/*----< ncmpi_get_vara() >---------------------------------------------------*/
/*----< ncmpi_get_vara_all() >-----------------------------------------------*/
GET_VARA(    , INDEP_IO)
GET_VARA(_all, COLL_IO)

#define PUT_VARA_TYPE(fntype, buftype, mpitype, collmode)        \
int                                                              \
ncmpi_put_vara_##fntype(int               ncid,                  \
                        int               varid,                 \
                        const MPI_Offset  start[],               \
                        const MPI_Offset  count[],               \
                        const buftype    *op)                    \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems;                                          \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_NUM_ELEMENTS                                             \
                                                                 \
    /* put_vara is a special case of put_vars */                 \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,     \
                              (void*)op, nelems, mpitype,        \
                              WRITE_REQ, collmode);              \
}

/*----< ncmpi_put_vara_text() >----------------------------------------------*/
/*----< ncmpi_put_vara_schar() >---------------------------------------------*/
/*----< ncmpi_put_vara_uchar() >---------------------------------------------*/
/*----< ncmpi_put_vara_short() >---------------------------------------------*/
/*----< ncmpi_put_vara_int() >-----------------------------------------------*/
/*----< ncmpi_put_vara_long() >----------------------------------------------*/
/*----< ncmpi_put_vara_float() >---------------------------------------------*/
/*----< ncmpi_put_vara_double() >--------------------------------------------*/

PUT_VARA_TYPE(text,   char,   MPI_CHAR,              INDEP_IO)
PUT_VARA_TYPE(schar,  schar,  MPI_BYTE,              INDEP_IO)
PUT_VARA_TYPE(uchar,  uchar,  MPI_UNSIGNED_CHAR,     INDEP_IO)
PUT_VARA_TYPE(short,  short,  MPI_SHORT,             INDEP_IO)
PUT_VARA_TYPE(int,    int,    MPI_INT,               INDEP_IO)
PUT_VARA_TYPE(long,   long,   MPI_LONG,              INDEP_IO)
PUT_VARA_TYPE(float,  float,  MPI_FLOAT,             INDEP_IO)
PUT_VARA_TYPE(double, double, MPI_DOUBLE,            INDEP_IO)

/*----< ncmpi_put_vara_text_all() >------------------------------------------*/
/*----< ncmpi_put_vara_schar_all() >-----------------------------------------*/
/*----< ncmpi_put_vara_uchar_all() >-----------------------------------------*/
/*----< ncmpi_put_vara_short_all() >-----------------------------------------*/
/*----< ncmpi_put_vara_int_all() >-------------------------------------------*/
/*----< ncmpi_put_vara_long_all() >------------------------------------------*/
/*----< ncmpi_put_vara_float_all() >-----------------------------------------*/
/*----< ncmpi_put_vara_double_all() >----------------------------------------*/

PUT_VARA_TYPE(text_all,   char,   MPI_CHAR,          COLL_IO)
PUT_VARA_TYPE(schar_all,  schar,  MPI_BYTE,          COLL_IO)
PUT_VARA_TYPE(uchar_all,  uchar,  MPI_UNSIGNED_CHAR, COLL_IO)
PUT_VARA_TYPE(short_all,  short,  MPI_SHORT,         COLL_IO)
PUT_VARA_TYPE(int_all,    int,    MPI_INT,           COLL_IO)
PUT_VARA_TYPE(long_all,   long,   MPI_LONG,          COLL_IO)
PUT_VARA_TYPE(float_all,  float,  MPI_FLOAT,         COLL_IO)
PUT_VARA_TYPE(double_all, double, MPI_DOUBLE,        COLL_IO)


#define GET_VARA_TYPE(fntype, buftype, mpitype, collmode)        \
int                                                              \
ncmpi_get_vara_##fntype(int               ncid,                  \
                        int               varid,                 \
                        const MPI_Offset  start[],               \
                        const MPI_Offset  count[],               \
                        buftype          *ip)                    \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems;                                          \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_NUM_ELEMENTS                                             \
                                                                 \
    /* get_vara is a special case of get_vars */                 \
    return ncmpii_getput_vars(ncp, varp, start, count, NULL,     \
                              ip, nelems, mpitype,               \
                              READ_REQ, collmode);               \
}

/*----< ncmpi_get_vara_text() >----------------------------------------------*/
/*----< ncmpi_get_vara_schar() >---------------------------------------------*/
/*----< ncmpi_get_vara_uchar() >---------------------------------------------*/
/*----< ncmpi_get_vara_short() >---------------------------------------------*/
/*----< ncmpi_get_vara_int() >-----------------------------------------------*/
/*----< ncmpi_get_vara_long() >----------------------------------------------*/
/*----< ncmpi_get_vara_float() >---------------------------------------------*/
/*----< ncmpi_get_vara_double() >--------------------------------------------*/

GET_VARA_TYPE(text,   char,   MPI_CHAR,              INDEP_IO)
GET_VARA_TYPE(schar,  schar,  MPI_BYTE,              INDEP_IO)
GET_VARA_TYPE(uchar,  uchar,  MPI_UNSIGNED_CHAR,     INDEP_IO)
GET_VARA_TYPE(short,  short,  MPI_SHORT,             INDEP_IO)
GET_VARA_TYPE(int,    int,    MPI_INT,               INDEP_IO)
GET_VARA_TYPE(long,   long,   MPI_LONG,              INDEP_IO)
GET_VARA_TYPE(float,  float,  MPI_FLOAT,             INDEP_IO)
GET_VARA_TYPE(double, double, MPI_DOUBLE,            INDEP_IO)

/*----< ncmpi_get_vara_text_all() >------------------------------------------*/
/*----< ncmpi_get_vara_schar_all() >-----------------------------------------*/
/*----< ncmpi_get_vara_uchar_all() >-----------------------------------------*/
/*----< ncmpi_get_vara_short_all() >-----------------------------------------*/
/*----< ncmpi_get_vara_int_all() >-------------------------------------------*/
/*----< ncmpi_get_vara_long_all() >------------------------------------------*/
/*----< ncmpi_get_vara_float_all() >-----------------------------------------*/
/*----< ncmpi_get_vara_double_all() >----------------------------------------*/

GET_VARA_TYPE(text_all,   char,   MPI_CHAR,          COLL_IO)
GET_VARA_TYPE(schar_all,  schar,  MPI_BYTE,          COLL_IO)
GET_VARA_TYPE(uchar_all,  uchar,  MPI_UNSIGNED_CHAR, COLL_IO)
GET_VARA_TYPE(short_all,  short,  MPI_SHORT,         COLL_IO)
GET_VARA_TYPE(int_all,    int,    MPI_INT,           COLL_IO)
GET_VARA_TYPE(long_all,   long,   MPI_LONG,          COLL_IO)
GET_VARA_TYPE(float_all,  float,  MPI_FLOAT,         COLL_IO)
GET_VARA_TYPE(double_all, double, MPI_DOUBLE,        COLL_IO)

