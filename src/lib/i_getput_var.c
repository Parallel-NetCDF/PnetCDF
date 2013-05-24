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

/*----< ncmpi_iput_var() >---------------------------------------------------*/
int
ncmpi_iput_var(int           ncid,
               int           varid,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int          *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp=NULL;
    MPI_Offset *start, *count;

    *reqid = NC_REQ_NULL;
    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_COLL_IO, status)

    GET_FULL_DIMENSIONS(start, count)

    /* iput_var is a special case of iput_varm */
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,
                                 (void*)buf, bufcount, buftype, reqid,
                                 WRITE_REQ, 0);
    if (varp->ndims > 0) NCI_Free(start);

    return status;
}


#define IPUT_VAR_TYPE(apitype, btype, buftype)                          \
int                                                                     \
ncmpi_iput_var_##apitype(int          ncid,                             \
                         int          varid,                            \
                         const btype *op,                               \
                         int         *reqid)                            \
{                                                                       \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp=NULL;                                              \
    MPI_Offset  nelems, *start, *count;                                 \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_COLL_IO, status)     \
                                                                        \
    GET_TOTAL_NUM_ELEMENTS(nelems)                                      \
    GET_FULL_DIMENSIONS(start, count)                                   \
                                                                        \
    /* iput_var is a special case of iput_varm */                       \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 (void*)op, nelems, buftype, reqid,     \
                                 WRITE_REQ, 0);                         \
    if (varp->ndims > 0) NCI_Free(start);                               \
                                                                        \
    return status;                                                      \
}

/*----< ncmpi_iput_var_text() >-----------------------------------------------*/
/*----< ncmpi_iput_var_schar() >----------------------------------------------*/
/*----< ncmpi_iput_var_uchar() >----------------------------------------------*/
/*----< ncmpi_iput_var_short() >----------------------------------------------*/
/*----< ncmpi_iput_var_ushort() >---------------------------------------------*/
/*----< ncmpi_iput_var_int() >------------------------------------------------*/
/*----< ncmpi_iput_var_uint() >-----------------------------------------------*/
/*----< ncmpi_iput_var_long() >-----------------------------------------------*/
/*----< ncmpi_iput_var_float() >----------------------------------------------*/
/*----< ncmpi_iput_var_double() >---------------------------------------------*/
/*----< ncmpi_iput_var_longlong() >-------------------------------------------*/
/*----< ncmpi_iput_var_ulonglong() >------------------------------------------*/
IPUT_VAR_TYPE(text,      char,               MPI_CHAR)
IPUT_VAR_TYPE(schar,     schar,              MPI_BYTE)
IPUT_VAR_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
IPUT_VAR_TYPE(short,     short,              MPI_SHORT)
IPUT_VAR_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
IPUT_VAR_TYPE(int,       int,                MPI_INT)
IPUT_VAR_TYPE(uint,      uint,               MPI_UNSIGNED)
IPUT_VAR_TYPE(long,      long,               MPI_LONG)
IPUT_VAR_TYPE(float,     float,              MPI_FLOAT)
IPUT_VAR_TYPE(double,    double,             MPI_DOUBLE)
IPUT_VAR_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
IPUT_VAR_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// IPUT_VAR_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */


/*----< ncmpi_iget_var() >---------------------------------------------------*/
int
ncmpi_iget_var(int           ncid,
               int           varid,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype,
               int          *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp=NULL;
    MPI_Offset *start, *count;

    *reqid = NC_REQ_NULL;
    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_COLL_IO, status)

    GET_FULL_DIMENSIONS(start, count)

    /* iget_var is a special case of iget_varm */
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,
                                 buf, bufcount, buftype, reqid,
                                 READ_REQ, 0);
    if (varp->ndims > 0) NCI_Free(start);

    return status;
}

#define IGET_VAR_TYPE(apitype, btype, buftype)                          \
int                                                                     \
ncmpi_iget_var_##apitype(int    ncid,                                   \
                         int    varid,                                  \
                         btype *ip,                                     \
                         int   *reqid)                                  \
{                                                                       \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp=NULL;                                              \
    MPI_Offset  nelems, *start, *count;                                 \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_COLL_IO, status)      \
                                                                        \
    GET_TOTAL_NUM_ELEMENTS(nelems)                                      \
    GET_FULL_DIMENSIONS(start, count)                                   \
                                                                        \
    /* iget_var is a special case of iget_vara */                       \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 ip, nelems, buftype, reqid,            \
                                 READ_REQ, 0);                          \
    if (varp->ndims > 0) NCI_Free(start);                               \
                                                                        \
    return status;                                                      \
}

/*----< ncmpi_iget_var_text() >-----------------------------------------------*/
/*----< ncmpi_iget_var_schar() >----------------------------------------------*/
/*----< ncmpi_iget_var_uchar() >----------------------------------------------*/
/*----< ncmpi_iget_var_short() >----------------------------------------------*/
/*----< ncmpi_iget_var_ushort() >---------------------------------------------*/
/*----< ncmpi_iget_var_int() >------------------------------------------------*/
/*----< ncmpi_iget_var_uint() >-----------------------------------------------*/
/*----< ncmpi_iget_var_long() >-----------------------------------------------*/
/*----< ncmpi_iget_var_float() >----------------------------------------------*/
/*----< ncmpi_iget_var_double() >---------------------------------------------*/
/*----< ncmpi_iget_var_longlong() >-------------------------------------------*/
/*----< ncmpi_iget_var_ulonglong() >------------------------------------------*/
IGET_VAR_TYPE(text,      char,               MPI_CHAR)
IGET_VAR_TYPE(schar,     schar,              MPI_BYTE)
IGET_VAR_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR)
IGET_VAR_TYPE(short,     short,              MPI_SHORT)
IGET_VAR_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT)
IGET_VAR_TYPE(int,       int,                MPI_INT)
IGET_VAR_TYPE(uint,      uint,               MPI_UNSIGNED)
IGET_VAR_TYPE(long,      long,               MPI_LONG)
IGET_VAR_TYPE(float,     float,              MPI_FLOAT)
IGET_VAR_TYPE(double,    double,             MPI_DOUBLE)
IGET_VAR_TYPE(longlong,  long long,          MPI_LONG_LONG_INT)
IGET_VAR_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
// IGET_VAR_TYPE(string, char*,              MPI_CHAR)
/* string is not yet supported */

