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

/*----< ncmpi_bput_varm() >--------------------------------------------------*/
int
ncmpi_bput_varm(int               ncid,
                int               varid,
                const MPI_Offset  start[],
                const MPI_Offset  count[],
                const MPI_Offset  stride[],
                const MPI_Offset  imap[],
                const void       *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      buftype,
                int              *reqid)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    *reqid = NC_REQ_NULL;
    CHECK_NCID
    if (ncp->abuf == NULL) return NC_ENULLABUF;
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    status = NCstrideedgeck(ncp, varp, start, count, stride);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap,
                               (void*)buf, bufcount, buftype, reqid,
                               WRITE_REQ, 1);
}

#define BPUT_VARM_TYPE(apitype, btype, buftype)                          \
int                                                                      \
ncmpi_bput_varm_##apitype(int               ncid,                        \
                          int               varid,                       \
                          const MPI_Offset  start[],                     \
                          const MPI_Offset  count[],                     \
                          const MPI_Offset  stride[],                    \
                          const MPI_Offset  imap[],                      \
                          const btype      *op,                          \
                          int              *reqid)                       \
{                                                                        \
    int         status;                                                  \
    NC         *ncp;                                                     \
    NC_var     *varp;                                                    \
    MPI_Offset  nelems;                                                  \
                                                                         \
    *reqid = NC_REQ_NULL;                                                \
    CHECK_NCID                                                           \
    if (ncp->abuf == NULL) return NC_ENULLABUF;                          \
    CHECK_WRITE_PERMISSION                                               \
    if (NC_indef(ncp)) return NC_EINDEFINE;                              \
    CHECK_VARID(varid, varp)                                             \
    status = NCcoordck(ncp, varp, start);                                \
    if (status != NC_NOERR) return status;                               \
    status = NCstrideedgeck(ncp, varp, start, count, stride);            \
    if (status != NC_NOERR) return status;                               \
    GET_NUM_ELEMENTS                                                     \
                                                                         \
    return ncmpii_igetput_varm(ncp, varp, start, count, stride, imap,    \
                               (void*)op, nelems, buftype, reqid,        \
                               WRITE_REQ, 1);                            \
}

/*----< ncmpi_bput_varm_text() >----------------------------------------------*/
/*----< ncmpi_bput_varm_schar() >---------------------------------------------*/
/*----< ncmpi_bput_varm_uchar() >---------------------------------------------*/
/*----< ncmpi_bput_varm_short() >---------------------------------------------*/
/*----< ncmpi_bput_varm_ushort() >--------------------------------------------*/
/*----< ncmpi_bput_varm_int() >-----------------------------------------------*/
/*----< ncmpi_bput_varm_uint() >----------------------------------------------*/
/*----< ncmpi_bput_varm_long() >----------------------------------------------*/
/*----< ncmpi_bput_varm_float() >---------------------------------------------*/
/*----< ncmpi_bput_varm_double() >--------------------------------------------*/
/*----< ncmpi_bput_varm_longlong() >------------------------------------------*/
/*----< ncmpi_bput_varm_ulonglong() >-----------------------------------------*/
BPUT_VARM_TYPE(text,      char,               MPI_CHAR)
BPUT_VARM_TYPE(schar,     schar,              MPI_BYTE)
BPUT_VARM_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
BPUT_VARM_TYPE(short,     short,              MPI_SHORT)
BPUT_VARM_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
BPUT_VARM_TYPE(int,       int,                MPI_INT)
BPUT_VARM_TYPE(uint,      uint,               MPI_UNSIGNED)
BPUT_VARM_TYPE(long,      long,               MPI_LONG)
BPUT_VARM_TYPE(float,     float,              MPI_FLOAT)
BPUT_VARM_TYPE(double,    double,             MPI_DOUBLE)
BPUT_VARM_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
BPUT_VARM_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// BPUT_VARM_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

