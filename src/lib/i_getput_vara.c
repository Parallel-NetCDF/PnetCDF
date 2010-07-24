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

/*----< ncmpi_iput_vara() >--------------------------------------------------*/
int
ncmpi_iput_vara(int               ncid,
                int               varid,
                const MPI_Offset *start,
                const MPI_Offset *count,
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
    status = NCedgeck(ncp, varp, start, count);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,
                               (void*)buf, bufcount, datatype, reqid,
                               WRITE_REQ);
}

#define IPUT_VARA_COMMON(datatype)                                     \
    int         status;                                                \
    NC         *ncp;                                                   \
    NC_var     *varp;                                                  \
    MPI_Offset  nelems;                                                \
                                                                       \
    *reqid = NC_REQ_NULL;                                              \
    CHECK_NCID                                                         \
    CHECK_WRITE_PERMISSION                                             \
    if (NC_indef(ncp)) return NC_EINDEFINE;                            \
    CHECK_VARID(varid, varp)                                           \
    status = NCcoordck(ncp, varp, start);                              \
    if (status != NC_NOERR) return status;                             \
    status = NCedgeck(ncp, varp, start, count);                        \
    if (status != NC_NOERR) return status;                             \
    GET_NUM_ELEMENTS                                                   \
                                                                       \
    return ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,    \
                               (void*)op, nelems, datatype, reqid,     \
                               WRITE_REQ);

/*----< ncmpi_iput_vara_uchar() >--------------------------------------------*/
int
ncmpi_iput_vara_uchar(int                  ncid,
                      int                  varid,
                      const MPI_Offset     start[],
                      const MPI_Offset     count[],
                      const unsigned char *op,
                      int                 *reqid)
{
    IPUT_VARA_COMMON(MPI_UNSIGNED_CHAR);
}

/*----< ncmpi_iput_vara_schar() >--------------------------------------------*/
int
ncmpi_iput_vara_schar(int                ncid,
                      int                varid,
                      const MPI_Offset   start[],
                      const MPI_Offset   count[],
                      const signed char *op,
                      int               *reqid)
{
    IPUT_VARA_COMMON(MPI_BYTE);
}

/*----< ncmpi_iput_vara_text() >---------------------------------------------*/
int
ncmpi_iput_vara_text(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const char       *op,
                     int              *reqid)
{
    IPUT_VARA_COMMON(MPI_CHAR);
}

/*----< ncmpi_iput_vara_short() >--------------------------------------------*/
int
ncmpi_iput_vara_short(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const short      *op,
                      int              *reqid)
{
    IPUT_VARA_COMMON(MPI_SHORT);
} 

/*----< ncmpi_iput_vara_int() >----------------------------------------------*/
int
ncmpi_iput_vara_int(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const int        *op,
                    int              *reqid)
{
    IPUT_VARA_COMMON(MPI_INT);
}

/*----< ncmpi_iput_vara_long() >---------------------------------------------*/
int
ncmpi_iput_vara_long(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const long       *op,
                     int              *reqid)
{
    IPUT_VARA_COMMON(MPI_LONG);
}

/*----< ncmpi_iput_vara_float() >--------------------------------------------*/
int
ncmpi_iput_vara_float(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const float      *op,
                      int              *reqid)
{
    IPUT_VARA_COMMON(MPI_FLOAT);
}

/*----< ncmpi_iput_vara_double() >-------------------------------------------*/
int
ncmpi_iput_vara_double(int              ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const double     *op,
                      int              *reqid)
{
    IPUT_VARA_COMMON(MPI_DOUBLE);
}

/*----< ncmpi_iget_vara() >--------------------------------------------------*/
int
ncmpi_iget_vara(int               ncid,
                int               varid,
                const MPI_Offset *start,
                const MPI_Offset *count,
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
    status = NCedgeck(ncp, varp, start, count);
    if (status != NC_NOERR) return status;

    return ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL, buf,
                               bufcount, datatype, reqid, READ_REQ);
}

#define IGET_VARA_COMMON(datatype)                                     \
    int         status;                                                \
    NC         *ncp;                                                   \
    NC_var     *varp;                                                  \
    MPI_Offset  nelems;                                                \
                                                                       \
    *reqid = NC_REQ_NULL;                                              \
    CHECK_NCID                                                         \
    if (NC_indef(ncp)) return NC_EINDEFINE;                            \
    CHECK_VARID(varid, varp)                                           \
    status = NCcoordck(ncp, varp, start);                              \
    if (status != NC_NOERR) return status;                             \
    status = NCedgeck(ncp, varp, start, count);                        \
    if (status != NC_NOERR) return status;                             \
    GET_NUM_ELEMENTS                                                   \
                                                                       \
    return ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,    \
                               ip, nelems, datatype, reqid, READ_REQ);

/*----< ncmpi_iget_vara_uchar() >--------------------------------------------*/
int
ncmpi_iget_vara_uchar(int               ncid,
                      int               varid,
                      const MPI_Offset  start[], 
                      const MPI_Offset  count[],
                      unsigned char    *ip,
                      int              *reqid)
{
    IGET_VARA_COMMON(MPI_UNSIGNED_CHAR);
}

/*----< ncmpi_iget_vara_schar() >--------------------------------------------*/
int
ncmpi_iget_vara_schar(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      signed char      *ip,       
                      int              *reqid)
{
    IGET_VARA_COMMON(MPI_BYTE);
}

/*----< ncmpi_iget_vara_text() >---------------------------------------------*/
int
ncmpi_iget_vara_text(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     char             *ip,       
                     int              *reqid)
{
    IGET_VARA_COMMON(MPI_CHAR);
}

/*----< ncmpi_iget_vara_short() >--------------------------------------------*/
int
ncmpi_iget_vara_short(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      short            *ip,       
                      int              *reqid)
{
    IGET_VARA_COMMON(MPI_SHORT);
} 

/*----< ncmpi_iget_vara_int() >----------------------------------------------*/
int
ncmpi_iget_vara_int(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    int              *ip,       
                    int              *reqid)
{
    IGET_VARA_COMMON(MPI_INT);
}

/*----< ncmpi_iget_vara_long() >---------------------------------------------*/
int
ncmpi_iget_vara_long(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     long             *ip,       
                     int              *reqid)
{
    IGET_VARA_COMMON(MPI_LONG);
}

/*----< ncmpi_iget_vara_float() >--------------------------------------------*/
int
ncmpi_iget_vara_float(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      float            *ip,
                      int              *reqid)
{
    IGET_VARA_COMMON(MPI_FLOAT);
}

/*----< ncmpi_iget_vara_double() >-------------------------------------------*/
int
ncmpi_iget_vara_double(int               ncid,
                       int               varid,
                       const MPI_Offset  start[],
                       const MPI_Offset  count[],
                       double           *ip,
                       int              *reqid)
{
    IGET_VARA_COMMON(MPI_DOUBLE);
}

