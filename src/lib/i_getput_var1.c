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

#define IPUT_VAR1_TYPE(fntype, buftype, mpitype)                        \
int                                                                     \
ncmpi_iput_var1_##fntype(int               ncid,                        \
                         int               varid,                       \
                         const MPI_Offset  start[],                     \
                         const buftype    *op,                          \
                         int              *reqid)                       \
{                                                                       \
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
                                 (void*)op, 1, mpitype, reqid,          \
                                 WRITE_REQ);                            \
    if (varp->ndims > 0) NCI_Free(count);                               \
    return status;                                                      \
}

/*----< ncmpi_iput_var1_text() >---------------------------------------------*/
/*----< ncmpi_iput_var1_schar() >--------------------------------------------*/
/*----< ncmpi_iput_var1_uchar() >--------------------------------------------*/
/*----< ncmpi_iput_var1_short() >--------------------------------------------*/
/*----< ncmpi_iput_var1_int() >----------------------------------------------*/
/*----< ncmpi_iput_var1_long() >---------------------------------------------*/
/*----< ncmpi_iput_var1_float() >--------------------------------------------*/
/*----< ncmpi_iput_var1_double() >-------------------------------------------*/

IPUT_VAR1_TYPE(text,   char,   MPI_CHAR)
IPUT_VAR1_TYPE(schar,  schar,  MPI_BYTE)
IPUT_VAR1_TYPE(uchar,  uchar,  MPI_UNSIGNED_CHAR)
IPUT_VAR1_TYPE(short,  short,  MPI_SHORT)
IPUT_VAR1_TYPE(int,    int,    MPI_INT)
IPUT_VAR1_TYPE(long,   long,   MPI_LONG)
IPUT_VAR1_TYPE(float,  float,  MPI_FLOAT)
IPUT_VAR1_TYPE(double, double, MPI_DOUBLE)


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

#define IGET_VAR1_TYPE(fntype, buftype, mpitype)                        \
int                                                                     \
ncmpi_iget_var1_##fntype(int               ncid,                        \
                         int               varid,                       \
                         const MPI_Offset  start[],                     \
                         buftype          *ip,                          \
                         int              *reqid)                       \
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
    status = NCcoordck(ncp, varp, start);                               \
    if (status != NC_NOERR) return status;                              \
    GET_ONE_COUNT                                                       \
                                                                        \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 ip, 1, mpitype, reqid, READ_REQ);      \
    if (varp->ndims > 0) NCI_Free(count);                               \
    return status;                                                      \
}

/*----< ncmpi_iget_var1_text() >---------------------------------------------*/
/*----< ncmpi_iget_var1_uchar() >--------------------------------------------*/
/*----< ncmpi_iget_var1_short() >--------------------------------------------*/
/*----< ncmpi_iget_var1_schar() >--------------------------------------------*/
/*----< ncmpi_iget_var1_int() >----------------------------------------------*/
/*----< ncmpi_iget_var1_long() >---------------------------------------------*/
/*----< ncmpi_iget_var1_float() >--------------------------------------------*/
/*----< ncmpi_iget_var1_double() >-------------------------------------------*/

IGET_VAR1_TYPE(text,   char,   MPI_CHAR)
IGET_VAR1_TYPE(schar,  schar,  MPI_BYTE)
IGET_VAR1_TYPE(uchar,  uchar,  MPI_UNSIGNED_CHAR)
IGET_VAR1_TYPE(short,  short,  MPI_SHORT)
IGET_VAR1_TYPE(int,    int,    MPI_INT)
IGET_VAR1_TYPE(long,   long,   MPI_LONG)
IGET_VAR1_TYPE(float,  float,  MPI_FLOAT)
IGET_VAR1_TYPE(double, double, MPI_DOUBLE)

