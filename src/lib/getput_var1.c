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

/*----< ncmpi_put_var1() >---------------------------------------------------*/
int
ncmpi_put_var1(int               ncid,
               int               varid,
               const MPI_Offset  start[],
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype)
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
                                (void*)buf, bufcount, buftype,
                                WRITE_REQ, INDEP_IO);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

/* ftype is the variable's nc_type defined in file
 * btype is the I/O buffer's data type
 */
#define PUT_VAR1_TYPE(apitype, btype, mpitype)                   \
int                                                              \
ncmpi_put_var1_##apitype(int               ncid,                 \
                         int               varid,                \
                         const MPI_Offset  start[],              \
                         const btype      *op)                   \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset *count;                                           \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    CHECK_INDEP_FH                                               \
    CHECK_VARID(varid, varp)                                     \
    GET_ONE_COUNT                                                \
                                                                 \
    /* put_var1 is a special case of put_vars */                 \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                (void*)op, 1, mpitype,           \
                                WRITE_REQ, INDEP_IO);            \
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
PUT_VAR1_TYPE(text,      char,               MPI_CHAR)
PUT_VAR1_TYPE(schar,     schar,              MPI_BYTE)
PUT_VAR1_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
PUT_VAR1_TYPE(short,     short,              MPI_SHORT)
PUT_VAR1_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
PUT_VAR1_TYPE(int,       int,                MPI_INT)
PUT_VAR1_TYPE(uint,      uint,               MPI_UNSIGNED)
PUT_VAR1_TYPE(long,      long,               MPI_LONG)
PUT_VAR1_TYPE(float,     float,              MPI_FLOAT)
PUT_VAR1_TYPE(double,    double,             MPI_DOUBLE)
PUT_VAR1_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
PUT_VAR1_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// PUT_VAR1_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

/*----< ncmpi_get_var1() >---------------------------------------------------*/
int
ncmpi_get_var1(int               ncid,
               int               varid,
               const MPI_Offset  start[],
               void             *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      buftype)
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
                                buf, bufcount, buftype,
                                READ_REQ, INDEP_IO);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

/* ftype is the variable's nc_type defined in file
 * btype is the I/O buffer's data type
 */
#define GET_VAR1_TYPE(apitype, btype, mpitype)                   \
int                                                              \
ncmpi_get_var1_##apitype(int               ncid,                 \
                         int               varid,                \
                         const MPI_Offset  start[],              \
                         btype            *ip)                   \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset *count;                                           \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    CHECK_INDEP_FH                                               \
    CHECK_VARID(varid, varp)                                     \
    GET_ONE_COUNT                                                \
                                                                 \
    /* get_var1 is a special case of get_vars */                 \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                ip, 1, mpitype,                  \
                                READ_REQ, INDEP_IO);             \
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
GET_VAR1_TYPE(text,      char,               MPI_CHAR)
GET_VAR1_TYPE(schar,     schar,              MPI_BYTE)
GET_VAR1_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
GET_VAR1_TYPE(short,     short,              MPI_SHORT)
GET_VAR1_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
GET_VAR1_TYPE(int,       int,                MPI_INT)
GET_VAR1_TYPE(uint,      uint,               MPI_UNSIGNED)
GET_VAR1_TYPE(long,      long,               MPI_LONG)
GET_VAR1_TYPE(float,     float,              MPI_FLOAT)
GET_VAR1_TYPE(double,    double,             MPI_DOUBLE)
GET_VAR1_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
GET_VAR1_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// GET_VAR1_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

