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

#define PUT_VAR1_TYPE(fntype, buftype, mpitype, collmode)        \
int                                                              \
ncmpi_put_var1_##fntype(int               ncid,                  \
                        int               varid,                 \
                        const MPI_Offset  start[],               \
                        const buftype    *op)                    \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset *count;                                           \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_ONE_COUNT                                                \
                                                                 \
    /* put_var1 is a special case of put_vars */                 \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                (void*)op, 1, mpitype,           \
                                WRITE_REQ, collmode);            \
    if (varp->ndims > 0) NCI_Free(count);                        \
                                                                 \
    return status;                                               \
}

/*----< ncmpi_put_var1_text() >-----------------------------------------------*/
/*----< ncmpi_put_var1_schar() >----------------------------------------------*/
/*----< ncmpi_put_var1_uchar() >----------------------------------------------*/
/*----< ncmpi_put_var1_short() >----------------------------------------------*/
/*----< ncmpi_put_var1_ushort() >---------------------------------------------*/
/*----< ncmpi_put_var1_int() >------------------------------------------------*/
/*----< ncmpi_put_var1_uint() >-----------------------------------------------*/
/*----< ncmpi_put_var1_long() >-----------------------------------------------*/
/*----< ncmpi_put_var1_float() >----------------------------------------------*/
/*----< ncmpi_put_var1_double() >---------------------------------------------*/
/*----< ncmpi_put_var1_longlong() >-------------------------------------------*/
/*----< ncmpi_put_var1_ulonglong() >------------------------------------------*/
PUT_VAR1_TYPE(text,      char,               MPI_CHAR,               INDEP_IO)
PUT_VAR1_TYPE(schar,     schar,              MPI_BYTE,               INDEP_IO)
PUT_VAR1_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR,      INDEP_IO)
PUT_VAR1_TYPE(short,     short,              MPI_SHORT,              INDEP_IO)
PUT_VAR1_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT,     INDEP_IO)
PUT_VAR1_TYPE(int,       int,                MPI_INT,                INDEP_IO)
PUT_VAR1_TYPE(uint,      uint,               MPI_UNSIGNED,           INDEP_IO)
PUT_VAR1_TYPE(long,      long,               MPI_LONG,               INDEP_IO)
PUT_VAR1_TYPE(float,     float,              MPI_FLOAT,              INDEP_IO)
PUT_VAR1_TYPE(double,    double,             MPI_DOUBLE,             INDEP_IO)
PUT_VAR1_TYPE(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
PUT_VAR1_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)
// PUT_VAR1_TYPE(string, char*,              MPI_CHAR,               INDEP_IO)
/* string is not yet supported */

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

#define GET_VAR1_TYPE(fntype, buftype, mpitype, collmode)        \
int                                                              \
ncmpi_get_var1_##fntype(int               ncid,                  \
                        int               varid,                 \
                        const MPI_Offset  start[],               \
                        buftype          *ip)                    \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset *count;                                           \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_ONE_COUNT                                                \
                                                                 \
    /* get_var1 is a special case of get_vars */                 \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                ip, 1, mpitype,                  \
                                READ_REQ, collmode);             \
    if (varp->ndims > 0) NCI_Free(count);                        \
                                                                 \
    return status;                                               \
}

/*----< ncmpi_get_var1_text() >-----------------------------------------------*/
/*----< ncmpi_get_var1_schar() >----------------------------------------------*/
/*----< ncmpi_get_var1_uchar() >----------------------------------------------*/
/*----< ncmpi_get_var1_short() >----------------------------------------------*/
/*----< ncmpi_get_var1_ushort() >---------------------------------------------*/
/*----< ncmpi_get_var1_int() >------------------------------------------------*/
/*----< ncmpi_get_var1_uint() >-----------------------------------------------*/
/*----< ncmpi_get_var1_long() >-----------------------------------------------*/
/*----< ncmpi_get_var1_float() >----------------------------------------------*/
/*----< ncmpi_get_var1_double() >---------------------------------------------*/
/*----< ncmpi_get_var1_longlong() >-------------------------------------------*/
/*----< ncmpi_get_var1_ulonglong() >------------------------------------------*/
GET_VAR1_TYPE(text,      char,               MPI_CHAR,               INDEP_IO)
GET_VAR1_TYPE(schar,     schar,              MPI_BYTE,               INDEP_IO)
GET_VAR1_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR,      INDEP_IO)
GET_VAR1_TYPE(short,     short,              MPI_SHORT,              INDEP_IO)
GET_VAR1_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT,     INDEP_IO)
GET_VAR1_TYPE(int,       int,                MPI_INT,                INDEP_IO)
GET_VAR1_TYPE(uint,      uint,               MPI_UNSIGNED,           INDEP_IO)
GET_VAR1_TYPE(long,      long,               MPI_LONG,               INDEP_IO)
GET_VAR1_TYPE(float,     float,              MPI_FLOAT,              INDEP_IO)
GET_VAR1_TYPE(double,    double,             MPI_DOUBLE,             INDEP_IO)
GET_VAR1_TYPE(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
GET_VAR1_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)
// GET_VAR1_TYPE(string, char*,              MPI_CHAR,               INDEP_IO)
/* string is not yet supported */

