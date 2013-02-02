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
 * apitype is data type appeared in the API names, eg. ncmpi_put_vara_longlong
 */

/*----< ncmpi_bput_var1() >--------------------------------------------------*/
int
ncmpi_bput_var1(int               ncid,
                int               varid,
                const MPI_Offset *start,
                const void       *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      buftype,
                int              *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *count;

    *reqid = NC_REQ_NULL;
    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    CHECK_WRITE_PERMISSION

    if (ncp->abuf == NULL) return NC_ENULLABUF;
    GET_ONE_COUNT
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;

    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,
                                 (void*)buf, bufcount, buftype, reqid,
                                 WRITE_REQ, 1);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

#define BPUT_VAR1_TYPE(apitype, btype, buftype)                         \
int                                                                     \
ncmpi_bput_var1_##apitype(int               ncid,                       \
                          int               varid,                      \
                          const MPI_Offset  start[],                    \
                          const btype      *op,                         \
                          int              *reqid)                      \
{                                                                       \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp;                                                   \
    MPI_Offset *count;                                                  \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    CHECK_NCID                                                          \
    if (NC_indef(ncp)) return NC_EINDEFINE;                             \
    CHECK_VARID(varid, varp)                                            \
    CHECK_WRITE_PERMISSION                                              \
                                                                        \
    if (ncp->abuf == NULL) return NC_ENULLABUF;                         \
    status = NCcoordck(ncp, varp, start);                               \
    if (status != NC_NOERR) return status;                              \
    GET_ONE_COUNT                                                       \
                                                                        \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 (void*)op, 1, buftype, reqid,          \
                                 WRITE_REQ, 1);                         \
    if (varp->ndims > 0) NCI_Free(count);                               \
    return status;                                                      \
}

/*----< ncmpi_bput_var1_text() >----------------------------------------------*/
/*----< ncmpi_bput_var1_schar() >---------------------------------------------*/
/*----< ncmpi_bput_var1_uchar() >---------------------------------------------*/
/*----< ncmpi_bput_var1_short() >---------------------------------------------*/
/*----< ncmpi_bput_var1_ushort() >--------------------------------------------*/
/*----< ncmpi_bput_var1_int() >-----------------------------------------------*/
/*----< ncmpi_bput_var1_uint() >----------------------------------------------*/
/*----< ncmpi_bput_var1_long() >----------------------------------------------*/
/*----< ncmpi_bput_var1_float() >---------------------------------------------*/
/*----< ncmpi_bput_var1_double() >--------------------------------------------*/
/*----< ncmpi_bput_var1_longlong() >------------------------------------------*/
/*----< ncmpi_bput_var1_ulonglong() >-----------------------------------------*/
BPUT_VAR1_TYPE(text,      char,               MPI_CHAR)
BPUT_VAR1_TYPE(schar,     schar,              MPI_BYTE)
BPUT_VAR1_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
BPUT_VAR1_TYPE(short,     short,              MPI_SHORT)
BPUT_VAR1_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
BPUT_VAR1_TYPE(int,       int,                MPI_INT)
BPUT_VAR1_TYPE(uint,      uint,               MPI_UNSIGNED)
BPUT_VAR1_TYPE(long,      long,               MPI_LONG)
BPUT_VAR1_TYPE(float,     float,              MPI_FLOAT)
BPUT_VAR1_TYPE(double,    double,             MPI_DOUBLE)
BPUT_VAR1_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
BPUT_VAR1_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// BPUT_VAR1_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */



