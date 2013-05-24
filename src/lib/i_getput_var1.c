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

#include "ncmpidtype.h"
#include "macro.h"


/* ftype is the variable's nc_type defined in file, eg. int64
 * btype is the I/O buffer's C data type, eg. long long
 * buftype is I/O bufer's MPI data type, eg. MPI_UNSIGNED_LONG_LONG
 * apitype is data type appeared in the API names, eg. ncmpi_get_vara_longlong
 */

/*----< ncmpi_iput_var1() >--------------------------------------------------*/
int
ncmpi_iput_var1(int               ncid,
                int               varid,
                const MPI_Offset *start,
                const void       *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      buftype,
                int              *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp=NULL;
    MPI_Offset *count;

    *reqid = NC_REQ_NULL;
    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_COLL_IO, status)

    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    GET_ONE_COUNT(count)

    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,
                                 (void*)buf, bufcount, buftype, reqid,
                                 WRITE_REQ, 0);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

#define IPUT_VAR1_TYPE(apitype, btype, buftype)                         \
int                                                                     \
ncmpi_iput_var1_##apitype(int               ncid,                       \
                          int               varid,                      \
                          const MPI_Offset  start[],                    \
                          const btype      *op,                         \
                          int              *reqid)                      \
{                                                                       \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp=NULL;                                              \
    MPI_Offset *count;                                                  \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_COLL_IO, status)     \
                                                                        \
    status = NCcoordck(ncp, varp, start);                               \
    if (status != NC_NOERR) return status;                              \
    GET_ONE_COUNT(count)                                                \
                                                                        \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 (void*)op, 1, buftype, reqid,          \
                                 WRITE_REQ, 0);                         \
    if (varp->ndims > 0) NCI_Free(count);                               \
    return status;                                                      \
}

/*----< ncmpi_iput_var1_text() >----------------------------------------------*/
/*----< ncmpi_iput_var1_schar() >---------------------------------------------*/
/*----< ncmpi_iput_var1_uchar() >---------------------------------------------*/
/*----< ncmpi_iput_var1_short() >---------------------------------------------*/
/*----< ncmpi_iput_var1_ushort() >--------------------------------------------*/
/*----< ncmpi_iput_var1_int() >-----------------------------------------------*/
/*----< ncmpi_iput_var1_uint() >----------------------------------------------*/
/*----< ncmpi_iput_var1_long() >----------------------------------------------*/
/*----< ncmpi_iput_var1_float() >---------------------------------------------*/
/*----< ncmpi_iput_var1_double() >--------------------------------------------*/
/*----< ncmpi_iput_var1_longlong() >------------------------------------------*/
/*----< ncmpi_iput_var1_ulonglong() >-----------------------------------------*/
IPUT_VAR1_TYPE(text,      char,               MPI_CHAR)
IPUT_VAR1_TYPE(schar,     schar,              MPI_BYTE)
IPUT_VAR1_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
IPUT_VAR1_TYPE(short,     short,              MPI_SHORT)
IPUT_VAR1_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
IPUT_VAR1_TYPE(int,       int,                MPI_INT)
IPUT_VAR1_TYPE(uint,      uint,               MPI_UNSIGNED)
IPUT_VAR1_TYPE(long,      long,               MPI_LONG)
IPUT_VAR1_TYPE(float,     float,              MPI_FLOAT)
IPUT_VAR1_TYPE(double,    double,             MPI_DOUBLE)
IPUT_VAR1_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
IPUT_VAR1_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// IPUT_VAR1_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */


/*----< ncmpi_iget_var1() >--------------------------------------------------*/
int
ncmpi_iget_var1(int               ncid,
                int               varid,
                const MPI_Offset *start,
                void             *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      buftype,
                int              *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp=NULL;
    MPI_Offset *count;

    *reqid = NC_REQ_NULL;
    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_COLL_IO, status)

    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    if (IS_RECVAR(varp) && start[0] + 1 > NC_get_numrecs(ncp)) return NC_EEDGE;
    GET_ONE_COUNT(count)

    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL, buf,
                                 bufcount, buftype, reqid, READ_REQ, 0);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

#define IGET_VAR1_TYPE(apitype, btype, buftype)                         \
int                                                                     \
ncmpi_iget_var1_##apitype(int               ncid,                       \
                          int               varid,                      \
                          const MPI_Offset  start[],                    \
                          btype            *ip,                         \
                          int              *reqid)                      \
{                                                                       \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp=NULL;                                              \
    MPI_Offset *count;                                                  \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_COLL_IO, status)      \
                                                                        \
    status = NCcoordck(ncp, varp, start);                               \
    if (status != NC_NOERR) return status;                              \
    if (IS_RECVAR(varp) &&                                              \
        start[0] + 1 > NC_get_numrecs(ncp)) return NC_EEDGE;            \
    GET_ONE_COUNT(count)                                                \
                                                                        \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 ip, 1, buftype, reqid, READ_REQ, 0);   \
    if (varp->ndims > 0) NCI_Free(count);                               \
    return status;                                                      \
}

/*----< ncmpi_iget_var1_text() >----------------------------------------------*/
/*----< ncmpi_iget_var1_schar() >---------------------------------------------*/
/*----< ncmpi_iget_var1_uchar() >---------------------------------------------*/
/*----< ncmpi_iget_var1_short() >---------------------------------------------*/
/*----< ncmpi_iget_var1_ushort() >--------------------------------------------*/
/*----< ncmpi_iget_var1_int() >-----------------------------------------------*/
/*----< ncmpi_iget_var1_uint() >----------------------------------------------*/
/*----< ncmpi_iget_var1_long() >----------------------------------------------*/
/*----< ncmpi_iget_var1_float() >---------------------------------------------*/
/*----< ncmpi_iget_var1_double() >--------------------------------------------*/
/*----< ncmpi_iget_var1_longlong() >------------------------------------------*/
/*----< ncmpi_iget_var1_ulonglong() >-----------------------------------------*/
IGET_VAR1_TYPE(text,      char,               MPI_CHAR)
IGET_VAR1_TYPE(schar,     schar,              MPI_BYTE)
IGET_VAR1_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
IGET_VAR1_TYPE(short,     short,              MPI_SHORT)
IGET_VAR1_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
IGET_VAR1_TYPE(int,       int,                MPI_INT)
IGET_VAR1_TYPE(uint,      uint,               MPI_UNSIGNED)
IGET_VAR1_TYPE(long,      long,               MPI_LONG)
IGET_VAR1_TYPE(float,     float,              MPI_FLOAT)
IGET_VAR1_TYPE(double,    double,             MPI_DOUBLE)
IGET_VAR1_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
IGET_VAR1_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// IGET_VAR1_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

