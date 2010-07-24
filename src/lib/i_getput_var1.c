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

/*----< ncmpi_iput_var1() >--------------------------------------------------*/
int
ncmpi_iput_var1(int               ncid,
                int               varid,
                const MPI_Offset *start,
                const void       *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      datatype,
                int              *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *count;

    *reqid = NC_REQ_NULL;
    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    GET_ONE_COUNT
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;


    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,
                                 (void*)buf, bufcount, datatype, reqid,
                                 WRITE_REQ);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

#define IPUT_VAR1_COMMON(datatype)                                      \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp;                                                   \
    MPI_Offset *count;                                                  \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    CHECK_NCID                                                          \
    CHECK_WRITE_PERMISSION                                              \
    if (NC_indef(ncp)) return NC_EINDEFINE;                             \
    CHECK_VARID(varid, varp)                                            \
    status = NCcoordck(ncp, varp, start);                               \
    if (status != NC_NOERR) return status;                              \
    GET_ONE_COUNT                                                       \
                                                                        \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 (void*)op, 1, datatype, reqid,         \
                                 WRITE_REQ);                            \
    if (varp->ndims > 0) NCI_Free(count);                               \
    return status;


/*----< ncmpi_iput_var1_uchar() >--------------------------------------------*/
int
ncmpi_iput_var1_uchar(int                  ncid,
                      int                  varid,
                      const MPI_Offset     start[],
                      const unsigned char *op,
                      int                 *reqid)
{
    IPUT_VAR1_COMMON(MPI_UNSIGNED_CHAR);
}

/*----< ncmpi_iput_var1_schar() >--------------------------------------------*/
int
ncmpi_iput_var1_schar(int                ncid,
                      int                varid,
                      const MPI_Offset   start[],
                      const signed char *op,
                      int               *reqid)
{
    IPUT_VAR1_COMMON(MPI_BYTE);
}

/*----< ncmpi_iput_var1_text() >---------------------------------------------*/
int
ncmpi_iput_var1_text(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const char       *op,
                     int              *reqid)
{
    IPUT_VAR1_COMMON(MPI_CHAR);
}

/*----< ncmpi_iput_var1_short() >--------------------------------------------*/
int
ncmpi_iput_var1_short(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const short      *op,
                      int              *reqid)
{
    IPUT_VAR1_COMMON(MPI_SHORT);
} 

/*----< ncmpi_iput_var1_int() >----------------------------------------------*/
int
ncmpi_iput_var1_int(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const int        *op,
                    int              *reqid)
{
    IPUT_VAR1_COMMON(MPI_INT);
}

/*----< ncmpi_iput_var1_long() >---------------------------------------------*/
int
ncmpi_iput_var1_long(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const long       *op,
                     int              *reqid)
{
    IPUT_VAR1_COMMON(MPI_LONG);
}

/*----< ncmpi_iput_var1_float() >--------------------------------------------*/
int
ncmpi_iput_var1_float(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const float      *op,
                      int              *reqid)
{
    IPUT_VAR1_COMMON(MPI_FLOAT);
}

/*----< ncmpi_iput_var1_double() >-------------------------------------------*/
int
ncmpi_iput_var1_double(int              ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const double     *op,
                      int              *reqid)
{
    IPUT_VAR1_COMMON(MPI_DOUBLE);
}

/*----< ncmpi_iget_var1() >--------------------------------------------------*/
int
ncmpi_iget_var1(int               ncid,
                int               varid,
                const MPI_Offset *start,
                void             *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      datatype,
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
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR) return status;
    GET_ONE_COUNT

    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL, buf,
                                 bufcount, datatype, reqid, READ_REQ);
    if (varp->ndims > 0) NCI_Free(count);
    return status;
}

#define IGET_VAR1_COMMON(datatype)                                      \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp;                                                   \
    MPI_Offset *count;                                                  \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    CHECK_NCID                                                          \
    if (NC_indef(ncp)) return NC_EINDEFINE;                             \
    CHECK_VARID(varid, varp)                                            \
    status = NCcoordck(ncp, varp, start);                               \
    if (status != NC_NOERR) return status;                              \
    GET_ONE_COUNT                                                       \
                                                                        \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 ip, 1, datatype, reqid, READ_REQ);     \
    if (varp->ndims > 0) NCI_Free(count);                               \
    return status;

/*----< ncmpi_iget_var1_uchar() >--------------------------------------------*/
int
ncmpi_iget_var1_uchar(int               ncid,
                      int               varid,
                      const MPI_Offset  start[], 
                      unsigned char    *ip,
                      int              *reqid)
{
    IGET_VAR1_COMMON(MPI_UNSIGNED_CHAR);
}

/*----< ncmpi_iget_var1_schar() >--------------------------------------------*/
int
ncmpi_iget_var1_schar(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      signed char      *ip,       
                      int              *reqid)
{
    IGET_VAR1_COMMON(MPI_BYTE);
}

/*----< ncmpi_iget_var1_text() >---------------------------------------------*/
int
ncmpi_iget_var1_text(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     char             *ip,       
                     int              *reqid)
{
    IGET_VAR1_COMMON(MPI_CHAR);
}

/*----< ncmpi_iget_var1_short() >--------------------------------------------*/
int
ncmpi_iget_var1_short(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      short            *ip,       
                      int              *reqid)
{
    IGET_VAR1_COMMON(MPI_SHORT);
} 

/*----< ncmpi_iget_var1_int() >----------------------------------------------*/
int
ncmpi_iget_var1_int(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    int              *ip,       
                    int              *reqid)
{
    IGET_VAR1_COMMON(MPI_INT);
}

/*----< ncmpi_iget_var1_long() >---------------------------------------------*/
int
ncmpi_iget_var1_long(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     long             *ip,       
                     int              *reqid)
{
    IGET_VAR1_COMMON(MPI_LONG);
}

/*----< ncmpi_iget_var1_float() >--------------------------------------------*/
int
ncmpi_iget_var1_float(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      float            *ip,
                      int              *reqid)
{
    IGET_VAR1_COMMON(MPI_FLOAT);
}

/*----< ncmpi_iget_var1_double() >-------------------------------------------*/
int
ncmpi_iget_var1_double(int               ncid,
                       int               varid,
                       const MPI_Offset  start[],
                       double           *ip,
                       int              *reqid)
{
    IGET_VAR1_COMMON(MPI_DOUBLE);
}

