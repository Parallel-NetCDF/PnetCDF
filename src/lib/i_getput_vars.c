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

#include "ncmpidtype.h"
#include "macro.h"


/* buffer layers:       
        
        User Level              buf     (user defined buffer of MPI_Datatype)
        MPI Datatype Level      cbuf    (contiguous buffer of ptype)
        NetCDF XDR Level        xbuf    (XDR I/O buffer)
*/

/*----< ncmpi_iput_vars() >--------------------------------------------------*/
int
ncmpi_iput_vars(int               ncid,
                int               varid,
                const MPI_Offset  start[],
                const MPI_Offset  count[],
                const MPI_Offset  stride[],
                const void       *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      datatype,
                int              *reqid)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    *reqid = NC_REQ_NULL;
    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, NULL,
                               (void*)buf, bufcount, datatype, reqid,
                               WRITE_REQ);
}

#define IPUT_VARS_TYPE(fntype, buftype, mpitype)                         \
int                                                                      \
ncmpi_iput_vars_##fntype(int               ncid,                         \
                         int               varid,                        \
                         const MPI_Offset  start[],                      \
                         const MPI_Offset  count[],                      \
                         const MPI_Offset  stride[],                     \
                         const buftype    *op,                           \
                         int              *reqid)                        \
{                                                                        \
    int         status;                                                  \
    NC         *ncp;                                                     \
    NC_var     *varp;                                                    \
    MPI_Offset  nelems;                                                  \
                                                                         \
    *reqid = NC_REQ_NULL;                                                \
    CHECK_NCID                                                           \
    CHECK_WRITE_PERMISSION                                               \
    if (NC_indef(ncp)) return NC_EINDEFINE;                              \
    CHECK_VARID(varid, varp)                                             \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    GET_NUM_ELEMENTS                                                     \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, NULL,    \
                               (void*)op, nelems, mpitype, reqid,        \
                               WRITE_REQ);                               \
}

/*----< ncmpi_iput_vars_text() >----------------------------------------------*/
/*----< ncmpi_iput_vars_schar() >---------------------------------------------*/
/*----< ncmpi_iput_vars_uchar() >---------------------------------------------*/
/*----< ncmpi_iput_vars_short() >---------------------------------------------*/
/*----< ncmpi_iput_vars_ushort() >--------------------------------------------*/
/*----< ncmpi_iput_vars_int() >-----------------------------------------------*/
/*----< ncmpi_iput_vars_uint() >----------------------------------------------*/
/*----< ncmpi_iput_vars_long() >----------------------------------------------*/
/*----< ncmpi_iput_vars_float() >---------------------------------------------*/
/*----< ncmpi_iput_vars_double() >--------------------------------------------*/
/*----< ncmpi_iput_vars_longlong() >------------------------------------------*/
/*----< ncmpi_iput_vars_ulonglong() >-----------------------------------------*/
IPUT_VARS_TYPE(text,      char,               MPI_CHAR)
IPUT_VARS_TYPE(schar,     schar,              MPI_BYTE)
IPUT_VARS_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
IPUT_VARS_TYPE(short,     short,              MPI_SHORT)
IPUT_VARS_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
IPUT_VARS_TYPE(int,       int,                MPI_INT)
IPUT_VARS_TYPE(uint,      uint,               MPI_UNSIGNED)
IPUT_VARS_TYPE(long,      long,               MPI_LONG)
IPUT_VARS_TYPE(float,     float,              MPI_FLOAT)
IPUT_VARS_TYPE(double,    double,             MPI_DOUBLE)
IPUT_VARS_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
IPUT_VARS_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// IPUT_VARS_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

/*----< ncmpi_iget_vars() >--------------------------------------------------*/
int
ncmpi_iget_vars(int               ncid,
                int               varid,
                const MPI_Offset  start[],
                const MPI_Offset  count[],
                const MPI_Offset  stride[],
                void             *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      datatype,
                int              *reqid)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    *reqid = NC_REQ_NULL;
    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, NULL,
                               buf, bufcount, datatype, reqid, READ_REQ);
}

#define IGET_VARS_TYPE(fntype, buftype, mpitype)                         \
int                                                                      \
ncmpi_iget_vars_##fntype(int               ncid,                         \
                         int               varid,                        \
                         const MPI_Offset  start[],                      \
                         const MPI_Offset  count[],                      \
                         const MPI_Offset  stride[],                     \
                         buftype          *ip,                           \
                         int              *reqid)                        \
{                                                                        \
    int         status;                                                  \
    NC         *ncp;                                                     \
    NC_var     *varp;                                                    \
    MPI_Offset  nelems;                                                  \
                                                                         \
    *reqid = NC_REQ_NULL;                                                \
    CHECK_NCID                                                           \
    if (NC_indef(ncp)) return NC_EINDEFINE;                              \
    CHECK_VARID(varid, varp)                                             \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    GET_NUM_ELEMENTS                                                     \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, NULL,    \
                               ip, nelems, mpitype, reqid, READ_REQ);    \
}


/*----< ncmpi_iget_vars_text() >----------------------------------------------*/
/*----< ncmpi_iget_vars_schar() >---------------------------------------------*/
/*----< ncmpi_iget_vars_uchar() >---------------------------------------------*/
/*----< ncmpi_iget_vars_short() >---------------------------------------------*/
/*----< ncmpi_iget_vars_ushort() >--------------------------------------------*/
/*----< ncmpi_iget_vars_int() >-----------------------------------------------*/
/*----< ncmpi_iget_vars_uint() >----------------------------------------------*/
/*----< ncmpi_iget_vars_long() >----------------------------------------------*/
/*----< ncmpi_iget_vars_float() >---------------------------------------------*/
/*----< ncmpi_iget_vars_double() >--------------------------------------------*/
/*----< ncmpi_iget_vars_longlong() >------------------------------------------*/
/*----< ncmpi_iget_vars_ulonglong() >-----------------------------------------*/
IGET_VARS_TYPE(text,      char,               MPI_CHAR)
IGET_VARS_TYPE(schar,     schar,              MPI_BYTE)
IGET_VARS_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
IGET_VARS_TYPE(short,     short,              MPI_SHORT)
IGET_VARS_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
IGET_VARS_TYPE(int,       int,                MPI_INT)
IGET_VARS_TYPE(uint,      uint,               MPI_UNSIGNED)
IGET_VARS_TYPE(long,      long,               MPI_LONG)
IGET_VARS_TYPE(float,     float,              MPI_FLOAT)
IGET_VARS_TYPE(double,    double,             MPI_DOUBLE)
IGET_VARS_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
IGET_VARS_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// IGET_VARS_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

