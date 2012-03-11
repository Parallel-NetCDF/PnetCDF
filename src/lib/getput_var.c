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

#define PUT_VAR(fnmode, collmode)                                \
int                                                              \
ncmpi_put_var##fnmode(int           ncid,                        \
                      int           varid,                       \
                      const void   *buf,                         \
                      MPI_Offset    bufcount,                    \
                      MPI_Datatype  datatype)                    \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset *start, *count;                                   \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_FULL_DIMENSIONS                                          \
                                                                 \
    /* put_var is a special case of put_vars */                  \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                (void*)buf, bufcount, datatype,  \
                                WRITE_REQ, collmode);            \
    if (varp->ndims > 0) NCI_Free(start);                        \
                                                                 \
    return status;                                               \
}

/*----< ncmpi_put_var() >-----------------------------------------------------*/
/*----< ncmpi_put_var_all() >-------------------------------------------------*/
PUT_VAR(    , INDEP_IO)
PUT_VAR(_all, COLL_IO)


#define GET_VAR(fnmode, collmode)                                \
int                                                              \
ncmpi_get_var##fnmode(int           ncid,                        \
                      int           varid,                       \
                      void         *buf,                         \
                      MPI_Offset    bufcount,                    \
                      MPI_Datatype  datatype)                    \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset *start, *count;                                   \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_FULL_DIMENSIONS                                          \
                                                                 \
    /* get_var is a special case of get_vars */                  \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                buf, bufcount, datatype,         \
                                READ_REQ, collmode);             \
    if (varp->ndims > 0) NCI_Free(start);                        \
                                                                 \
    return status;                                               \
}

/*----< ncmpi_get_var() >----------------------------------------------------*/
/*----< ncmpi_get_var_all() >------------------------------------------------*/
GET_VAR(    , INDEP_IO)
GET_VAR(_all, COLL_IO)


#define PUT_VAR_TYPE(fntype, buftype, mpitype, collmode)         \
int                                                              \
ncmpi_put_var_##fntype(int            ncid,                      \
                       int            varid,                     \
                       const buftype *op)                        \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems, *start, *count;                          \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_TOTAL_NUM_ELEMENTS                                       \
    GET_FULL_DIMENSIONS                                          \
                                                                 \
    /* put_var is a special case of put_vars */                  \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                (void*)op, nelems, mpitype,      \
                                WRITE_REQ, collmode);            \
    if (varp->ndims > 0) NCI_Free(start);                        \
                                                                 \
    return status;                                               \
}

/*----< ncmpi_put_var_text() >------------------------------------------------*/
/*----< ncmpi_put_var_schar() >-----------------------------------------------*/
/*----< ncmpi_put_var_uchar() >-----------------------------------------------*/
/*----< ncmpi_put_var_short() >-----------------------------------------------*/
/*----< ncmpi_put_var_ushort() >----------------------------------------------*/
/*----< ncmpi_put_var_int() >-------------------------------------------------*/
/*----< ncmpi_put_var_uint() >------------------------------------------------*/
/*----< ncmpi_put_var_long() >------------------------------------------------*/
/*----< ncmpi_put_var_float() >-----------------------------------------------*/
/*----< ncmpi_put_var_double() >----------------------------------------------*/
/*----< ncmpi_put_var_longlong() >--------------------------------------------*/
/*----< ncmpi_put_var_ulonglong() >-------------------------------------------*/
PUT_VAR_TYPE(text,      char,               MPI_CHAR,               INDEP_IO)
PUT_VAR_TYPE(schar,     schar,              MPI_BYTE,               INDEP_IO)
PUT_VAR_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR,      INDEP_IO)
PUT_VAR_TYPE(short,     short,              MPI_SHORT,              INDEP_IO)
PUT_VAR_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT,     INDEP_IO)
PUT_VAR_TYPE(int,       int,                MPI_INT,                INDEP_IO)
PUT_VAR_TYPE(uint,      uint,               MPI_UNSIGNED,           INDEP_IO)
PUT_VAR_TYPE(long,      long,               MPI_LONG,               INDEP_IO)
PUT_VAR_TYPE(float,     float,              MPI_FLOAT,              INDEP_IO)
PUT_VAR_TYPE(double,    double,             MPI_DOUBLE,             INDEP_IO)
PUT_VAR_TYPE(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
PUT_VAR_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)
// PUT_VAR_TYPE(string, char*,              MPI_CHAR,               INDEP_IO)
/* string is not yet supported */

/*----< ncmpi_put_var_text_all() >--------------------------------------------*/
/*----< ncmpi_put_var_schar_all() >-------------------------------------------*/
/*----< ncmpi_put_var_uchar_all() >-------------------------------------------*/
/*----< ncmpi_put_var_short_all() >-------------------------------------------*/
/*----< ncmpi_put_var_ushort_all() >------------------------------------------*/
/*----< ncmpi_put_var_int_all() >---------------------------------------------*/
/*----< ncmpi_put_var_uint_all() >--------------------------------------------*/
/*----< ncmpi_put_var_long_all() >--------------------------------------------*/
/*----< ncmpi_put_var_float_all() >-------------------------------------------*/
/*----< ncmpi_put_var_double_all() >------------------------------------------*/
/*----< ncmpi_put_var_longlong_all() >----------------------------------------*/
/*----< ncmpi_put_var_ulonglong_all() >---------------------------------------*/
PUT_VAR_TYPE(text_all,      char,               MPI_CHAR,               COLL_IO)
PUT_VAR_TYPE(schar_all,     schar,              MPI_BYTE,               COLL_IO)
PUT_VAR_TYPE(uchar_all,     uchar,              MPI_UNSIGNED_CHAR,      COLL_IO)
PUT_VAR_TYPE(short_all,     short,              MPI_SHORT,              COLL_IO)
PUT_VAR_TYPE(ushort_all,    ushort,             MPI_UNSIGNED_SHORT,     COLL_IO)
PUT_VAR_TYPE(int_all,       int,                MPI_INT,                COLL_IO)
PUT_VAR_TYPE(uint_all,      uint,               MPI_UNSIGNED,           COLL_IO)
PUT_VAR_TYPE(long_all,      long,               MPI_LONG,               COLL_IO)
PUT_VAR_TYPE(float_all,     float,              MPI_FLOAT,              COLL_IO)
PUT_VAR_TYPE(double_all,    double,             MPI_DOUBLE,             COLL_IO)
PUT_VAR_TYPE(longlong_all,  long long,          MPI_LONG_LONG_INT,      COLL_IO)
PUT_VAR_TYPE(ulonglong_all, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO)
// PUT_VAR_TYPE(string_all, char*,              MPI_CHAR,               COLL_IO)
/* string is not yet supported */


#define GET_VAR_TYPE(fntype, buftype, mpitype, collmode)         \
int                                                              \
ncmpi_get_var_##fntype(int      ncid,                            \
                       int      varid,                           \
                       buftype *ip)                              \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems, *start, *count;                          \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_TOTAL_NUM_ELEMENTS                                       \
    GET_FULL_DIMENSIONS                                          \
                                                                 \
    /* get_var is a special case of get_vars */                  \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                ip, nelems, mpitype,             \
                                READ_REQ, collmode);             \
    if (varp->ndims > 0) NCI_Free(start);                        \
                                                                 \
    return status;                                               \
}

/*----< ncmpi_get_var_text() >------------------------------------------------*/
/*----< ncmpi_get_var_schar() >-----------------------------------------------*/
/*----< ncmpi_get_var_uchar() >-----------------------------------------------*/
/*----< ncmpi_get_var_short() >-----------------------------------------------*/
/*----< ncmpi_get_var_ushort() >----------------------------------------------*/
/*----< ncmpi_get_var_int() >-------------------------------------------------*/
/*----< ncmpi_get_var_uint() >------------------------------------------------*/
/*----< ncmpi_get_var_long() >------------------------------------------------*/
/*----< ncmpi_get_var_float() >-----------------------------------------------*/
/*----< ncmpi_get_var_double() >----------------------------------------------*/
/*----< ncmpi_get_var_longlong() >--------------------------------------------*/
/*----< ncmpi_get_var_ulonglong() >-------------------------------------------*/
GET_VAR_TYPE(text,      char,               MPI_CHAR,               INDEP_IO)
GET_VAR_TYPE(schar,     schar,              MPI_BYTE,               INDEP_IO)
GET_VAR_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR,      INDEP_IO)
GET_VAR_TYPE(short,     short,              MPI_SHORT,              INDEP_IO)
GET_VAR_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT,     INDEP_IO)
GET_VAR_TYPE(int,       int,                MPI_INT,                INDEP_IO)
GET_VAR_TYPE(uint,      uint,               MPI_UNSIGNED,           INDEP_IO)
GET_VAR_TYPE(long,      long,               MPI_LONG,               INDEP_IO)
GET_VAR_TYPE(float,     float,              MPI_FLOAT,              INDEP_IO)
GET_VAR_TYPE(double,    double,             MPI_DOUBLE,             INDEP_IO)
GET_VAR_TYPE(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
GET_VAR_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)
// GET_VAR_TYPE(string, char*,              MPI_CHAR,               INDEP_IO)
/* string is not yet supported */

/*----< ncmpi_get_var_text_all() >--------------------------------------------*/
/*----< ncmpi_get_var_schar_all() >-------------------------------------------*/
/*----< ncmpi_get_var_uchar_all() >-------------------------------------------*/
/*----< ncmpi_get_var_short_all() >-------------------------------------------*/
/*----< ncmpi_get_var_ushort_all() >------------------------------------------*/
/*----< ncmpi_get_var_int_all() >---------------------------------------------*/
/*----< ncmpi_get_var_uint_all() >--------------------------------------------*/
/*----< ncmpi_get_var_long_all() >--------------------------------------------*/
/*----< ncmpi_get_var_float_all() >-------------------------------------------*/
/*----< ncmpi_get_var_double_all() >------------------------------------------*/
/*----< ncmpi_get_var_longlong_all() >----------------------------------------*/
/*----< ncmpi_get_var_ulonglong_all() >---------------------------------------*/
GET_VAR_TYPE(text_all,      char,               MPI_CHAR,               COLL_IO)
GET_VAR_TYPE(schar_all,     schar,              MPI_BYTE,               COLL_IO)
GET_VAR_TYPE(uchar_all,     uchar,              MPI_UNSIGNED_CHAR,      COLL_IO)
GET_VAR_TYPE(short_all,     short,              MPI_SHORT,              COLL_IO)
GET_VAR_TYPE(ushort_all,    ushort,             MPI_UNSIGNED_SHORT,     COLL_IO)
GET_VAR_TYPE(int_all,       int,                MPI_INT,                COLL_IO)
GET_VAR_TYPE(uint_all,      uint,               MPI_UNSIGNED,           COLL_IO)
GET_VAR_TYPE(long_all,      long,               MPI_LONG,               COLL_IO)
GET_VAR_TYPE(float_all,     float,              MPI_FLOAT,              COLL_IO)
GET_VAR_TYPE(double_all,    double,             MPI_DOUBLE,             COLL_IO)
GET_VAR_TYPE(longlong_all,  long long,          MPI_LONG_LONG_INT,      COLL_IO)
GET_VAR_TYPE(ulonglong_all, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO)
// GET_VAR_TYPE(string_all, char*,              MPI_CHAR,               COLL_IO)
/* string is not yet supported */

