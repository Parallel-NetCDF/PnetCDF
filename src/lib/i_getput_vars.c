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

#define IPUT_VARA_COMMON(datatype)                                       \
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
                               (void*)op, nelems, datatype, reqid,       \
                               WRITE_REQ);

/*----< ncmpi_iput_vars_uchar() >--------------------------------------------*/
int
ncmpi_iput_vars_uchar(int                  ncid,
                      int                  varid,
                      const MPI_Offset     start[],
                      const MPI_Offset     count[],
                      const MPI_Offset     stride[],
                      const unsigned char *op,
                      int                 *reqid)
{
    IPUT_VARA_COMMON(MPI_UNSIGNED_CHAR);
}

/*----< ncmpi_iput_vars_schar() >--------------------------------------------*/
int
ncmpi_iput_vars_schar(int                ncid,
                      int                varid,
                      const MPI_Offset   start[],
                      const MPI_Offset   count[],
                      const MPI_Offset   stride[],
                      const signed char *op,
                      int               *reqid)
{
    IPUT_VARA_COMMON(MPI_BYTE);
}

/*----< ncmpi_iput_vars_text() >---------------------------------------------*/
int
ncmpi_iput_vars_text(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     const char       *op,
                     int              *reqid)
{
    IPUT_VARA_COMMON(MPI_CHAR);
}

/*----< ncmpi_iput_vars_short() >--------------------------------------------*/
int
ncmpi_iput_vars_short(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const MPI_Offset  stride[],
                      const short      *op,
                      int              *reqid)
{
    IPUT_VARA_COMMON(MPI_SHORT);
} 

/*----< ncmpi_iput_vars_int() >----------------------------------------------*/
int
ncmpi_iput_vars_int(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    const int        *op,
                    int              *reqid)
{
    IPUT_VARA_COMMON(MPI_INT);
}

/*----< ncmpi_iput_vars_long() >---------------------------------------------*/
int
ncmpi_iput_vars_long(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     const long       *op,
                     int              *reqid)
{
    IPUT_VARA_COMMON(MPI_LONG);
}

/*----< ncmpi_iput_vars_float() >--------------------------------------------*/
int
ncmpi_iput_vars_float(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const MPI_Offset  stride[],
                      const float      *op,
                      int              *reqid)
{
    IPUT_VARA_COMMON(MPI_FLOAT);
}

/*----< ncmpi_iput_vars_double() >-------------------------------------------*/
int
ncmpi_iput_vars_double(int               ncid,
                       int               varid,
                       const MPI_Offset  start[],
                       const MPI_Offset  count[],
                       const MPI_Offset  stride[],
                       const double     *op,
                       int              *reqid)
{
    IPUT_VARA_COMMON(MPI_DOUBLE);
}

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

#define IGET_VARA_COMMON(datatype)                                       \
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
                               ip, nelems, datatype, reqid, READ_REQ);

/*----< ncmpi_iget_vars_uchar() >--------------------------------------------*/
int
ncmpi_iget_vars_uchar(int               ncid,
                      int               varid,
                      const MPI_Offset  start[], 
                      const MPI_Offset  count[],
                      const MPI_Offset  stride[],
                      unsigned char    *ip,
                      int              *reqid)
{
    IGET_VARA_COMMON(MPI_UNSIGNED_CHAR);
}

/*----< ncmpi_iget_vars_schar() >--------------------------------------------*/
int
ncmpi_iget_vars_schar(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const MPI_Offset  stride[],
                      signed char      *ip,       
                      int              *reqid)
{
    IGET_VARA_COMMON(MPI_BYTE);
}

/*----< ncmpi_iget_vars_text() >---------------------------------------------*/
int
ncmpi_iget_vars_text(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     char             *ip,       
                     int              *reqid)
{
    IGET_VARA_COMMON(MPI_CHAR);
}

/*----< ncmpi_iget_vars_short() >--------------------------------------------*/
int
ncmpi_iget_vars_short(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const MPI_Offset  stride[],
                      short            *ip,       
                      int              *reqid)
{
    IGET_VARA_COMMON(MPI_SHORT);
} 

/*----< ncmpi_iget_vars_int() >----------------------------------------------*/
int
ncmpi_iget_vars_int(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    int              *ip,       
                    int              *reqid)
{
    IGET_VARA_COMMON(MPI_INT);
}

/*----< ncmpi_iget_vars_long() >---------------------------------------------*/
int
ncmpi_iget_vars_long(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     long             *ip,       
                     int              *reqid)
{
    IGET_VARA_COMMON(MPI_LONG);
}

/*----< ncmpi_iget_vars_float() >--------------------------------------------*/
int
ncmpi_iget_vars_float(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const MPI_Offset  stride[],
                      float            *ip,
                      int              *reqid)
{
    IGET_VARA_COMMON(MPI_FLOAT);
}

/*----< ncmpi_iget_vars_double() >-------------------------------------------*/
int
ncmpi_iget_vars_double(int               ncid,
                       int               varid,
                       const MPI_Offset  start[],
                       const MPI_Offset  count[],
                       const MPI_Offset  stride[],
                       double           *ip,
                       int              *reqid)
{
    IGET_VARA_COMMON(MPI_DOUBLE);
}

