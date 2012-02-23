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

#define IPUT_VARA_TYPE(fntype, buftype, mpitype)                       \
int                                                                    \
ncmpi_iput_vara_##fntype(int               ncid,                       \
                         int               varid,                      \
                         const MPI_Offset  start[],                    \
                         const MPI_Offset  count[],                    \
                         const buftype    *op,                         \
                         int              *reqid)                      \
{                                                                      \
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
                               (void*)op, nelems, mpitype, reqid,      \
                               WRITE_REQ);                             \
}

/*----< ncmpi_iput_vara_text() >---------------------------------------------*/
/*----< ncmpi_iput_vara_schar() >--------------------------------------------*/
/*----< ncmpi_iput_vara_uchar() >--------------------------------------------*/
/*----< ncmpi_iput_vara_short() >--------------------------------------------*/
/*----< ncmpi_iput_vara_int() >----------------------------------------------*/
/*----< ncmpi_iput_vara_long() >---------------------------------------------*/
/*----< ncmpi_iput_vara_float() >--------------------------------------------*/
/*----< ncmpi_iput_vara_double() >-------------------------------------------*/

IPUT_VARA_TYPE(text,   char,   MPI_CHAR)
IPUT_VARA_TYPE(schar,  schar,  MPI_BYTE)
IPUT_VARA_TYPE(uchar,  uchar,  MPI_UNSIGNED_CHAR)
IPUT_VARA_TYPE(short,  short,  MPI_SHORT)
IPUT_VARA_TYPE(int,    int,    MPI_INT)
IPUT_VARA_TYPE(long,   long,   MPI_LONG)
IPUT_VARA_TYPE(float,  float,  MPI_FLOAT)
IPUT_VARA_TYPE(double, double, MPI_DOUBLE)


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

#define IGET_VARA_TYPE(fntype, buftype, mpitype)                       \
int                                                                    \
ncmpi_iget_vara_##fntype(int               ncid,                       \
                         int               varid,                      \
                         const MPI_Offset  start[],                    \
                         const MPI_Offset  count[],                    \
                         buftype          *ip,                         \
                         int              *reqid)                      \
{                                                                      \
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
                               ip, nelems, mpitype, reqid, READ_REQ);  \
}

/*----< ncmpi_iget_vara_text() >---------------------------------------------*/
/*----< ncmpi_iget_vara_schar() >--------------------------------------------*/
/*----< ncmpi_iget_vara_uchar() >--------------------------------------------*/
/*----< ncmpi_iget_vara_short() >--------------------------------------------*/
/*----< ncmpi_iget_vara_int() >----------------------------------------------*/
/*----< ncmpi_iget_vara_long() >---------------------------------------------*/
/*----< ncmpi_iget_vara_float() >--------------------------------------------*/
/*----< ncmpi_iget_vara_double() >-------------------------------------------*/

IGET_VARA_TYPE(text,   char,   MPI_CHAR)
IGET_VARA_TYPE(schar,  schar,  MPI_BYTE)
IGET_VARA_TYPE(uchar,  uchar,  MPI_UNSIGNED_CHAR)
IGET_VARA_TYPE(short,  short,  MPI_SHORT)
IGET_VARA_TYPE(int,    int,    MPI_INT)
IGET_VARA_TYPE(long,   long,   MPI_LONG)
IGET_VARA_TYPE(float,  float,  MPI_FLOAT)
IGET_VARA_TYPE(double, double, MPI_DOUBLE)

