/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

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


/* ftype is the variable's nc_type defined in file, eg. int64
 * btype is the I/O buffer's C data type, eg. long long
 * buftype is I/O bufer's MPI data type, eg. MPI_UNSIGNED_LONG_LONG
 * apitype is data type appeared in the API names, eg. ncmpi_get_vara_longlong
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
                MPI_Datatype      buftype,
                int              *reqid)
{
    int     status;
    NC     *ncp;
    NC_var *varp=NULL;

    *reqid = NC_REQ_NULL;
    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_COLL_IO, status)

    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, NULL,
                               (void*)buf, bufcount, buftype, reqid,
                               WRITE_REQ, 0);
}

#define IPUT_VARS_TYPE(apitype, btype, buftype)                          \
int                                                                      \
ncmpi_iput_vars_##apitype(int               ncid,                        \
                          int               varid,                       \
                          const MPI_Offset  start[],                     \
                          const MPI_Offset  count[],                     \
                          const MPI_Offset  stride[],                    \
                          const btype      *op,                          \
                          int              *reqid)                       \
{                                                                        \
    int         status;                                                  \
    NC         *ncp;                                                     \
    NC_var     *varp=NULL;                                               \
    MPI_Offset  nelems;                                                  \
                                                                         \
    *reqid = NC_REQ_NULL;                                                \
    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_COLL_IO, status)      \
                                                                         \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    GET_NUM_ELEMENTS(nelems)                                             \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, NULL,    \
                               (void*)op, nelems, buftype, reqid,        \
                               WRITE_REQ, 0);                            \
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
                MPI_Datatype      buftype,
                int              *reqid)
{
    int     status;
    NC     *ncp;
    NC_var *varp=NULL;

    *reqid = NC_REQ_NULL;
    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_COLL_IO, status)

    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;
    if (IS_RECVAR(varp) &&
        start[0] + count[0] > NC_get_numrecs(ncp)) return NC_EEDGE;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, NULL,
                               buf, bufcount, buftype, reqid, READ_REQ, 0);
}

#define IGET_VARS_TYPE(apitype, btype, buftype)                          \
int                                                                      \
ncmpi_iget_vars_##apitype(int               ncid,                        \
                          int               varid,                       \
                          const MPI_Offset  start[],                     \
                          const MPI_Offset  count[],                     \
                          const MPI_Offset  stride[],                    \
                          btype            *ip,                          \
                          int              *reqid)                       \
{                                                                        \
    int         status;                                                  \
    NC         *ncp;                                                     \
    NC_var     *varp=NULL;                                               \
    MPI_Offset  nelems;                                                  \
                                                                         \
    *reqid = NC_REQ_NULL;                                                \
    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_COLL_IO, status)       \
                                                                         \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    if (IS_RECVAR(varp) &&                                               \
        start[0] + count[0] > NC_get_numrecs(ncp)) return NC_EEDGE;      \
    GET_NUM_ELEMENTS(nelems)                                             \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, NULL,    \
                               ip, nelems, buftype, reqid, READ_REQ, 0); \
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

