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

/*----< ncmpi_iput_var() >---------------------------------------------------*/
int
ncmpi_iput_var(int           ncid,
               int           varid,
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  datatype,
               int          *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *start, *count;

    *reqid = NC_REQ_NULL;
    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    GET_FULL_DIMENSIONS

    /* iput_var is a special case of iput_vara */
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,
                                 (void*)buf, bufcount, datatype, reqid,
                                 WRITE_REQ);
    if (varp->ndims > 0) NCI_Free(start);

    return status;
}

#define IPUT_VAR_COMMON(datatype)                                       \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp;                                                   \
    MPI_Offset  nelems, *start, *count;                                 \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    CHECK_NCID                                                          \
    CHECK_WRITE_PERMISSION                                              \
    if (NC_indef(ncp)) return NC_EINDEFINE;                             \
    CHECK_VARID(varid, varp)                                            \
    GET_TOTAL_NUM_ELEMENTS                                              \
    GET_FULL_DIMENSIONS                                                 \
                                                                        \
    /* iput_var is a special case of iput_vara */                       \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 (void*)op, nelems, datatype, reqid,    \
                                 WRITE_REQ);                            \
    if (varp->ndims > 0) NCI_Free(start);                               \
                                                                        \
    return status;

/*----< ncmpi_iput_var_uchar() >---------------------------------------------*/
int
ncmpi_iput_var_uchar(int                  ncid,
                     int                  varid,
                     const unsigned char *op,
                     int                 *reqid)
{
    IPUT_VAR_COMMON(MPI_UNSIGNED_CHAR);
}

/*----< ncmpi_iput_var_schar() >---------------------------------------------*/
int
ncmpi_iput_var_schar(int                ncid,
                     int                varid,
                     const signed char *op,
                     int               *reqid)
{
    IPUT_VAR_COMMON(MPI_BYTE);
}

/*----< ncmpi_iput_var_text() >----------------------------------------------*/
int
ncmpi_iput_var_text(int         ncid,
                    int         varid,
                    const char *op,
                    int         *reqid)
{
    IPUT_VAR_COMMON(MPI_CHAR);
}

/*----< ncmpi_iput_var_short() >---------------------------------------------*/
int
ncmpi_iput_var_short(int          ncid,
                     int          varid,
                     const short *op,
                     int         *reqid)
{
    IPUT_VAR_COMMON(MPI_SHORT);
}

/*----< ncmpi_iput_var_int() >-----------------------------------------------*/
int
ncmpi_iput_var_int(int        ncid,
                   int        varid,
                   const int *op,
                   int       *reqid)
{
    IPUT_VAR_COMMON(MPI_INT);
} 

/*----< ncmpi_iput_var_long() >----------------------------------------------*/
int
ncmpi_iput_var_long(int         ncid,
                    int         varid,
                    const long *op,
                    int        *reqid)
{
    IPUT_VAR_COMMON(MPI_LONG);
}

/*----< ncmpi_iput_var_float() >---------------------------------------------*/
int
ncmpi_iput_var_float(int          ncid,
                     int          varid,
                     const float *op,
                     int         *reqid)
{
    IPUT_VAR_COMMON(MPI_FLOAT);
} 

/*----< ncmpi_iput_var_double() >--------------------------------------------*/
int
ncmpi_iput_var_double(int           ncid,
                      int           varid,
                      const double *op,
                      int          *reqid)
{
    IPUT_VAR_COMMON(MPI_DOUBLE);
} 

/*----< ncmpi_iget_var() >---------------------------------------------------*/
int
ncmpi_iget_var(int           ncid,
               int           varid,
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  datatype,
               int          *reqid)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *start, *count;

    *reqid = NC_REQ_NULL;
    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_VARID(varid, varp)
    GET_FULL_DIMENSIONS

    /* iget_var is a special case of iget_vara */
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,
                                 buf, bufcount, datatype, reqid,
                                 READ_REQ);
    if (varp->ndims > 0) NCI_Free(start);

    return status;
}

#define IGET_VAR_COMMON(datatype)                                       \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp;                                                   \
    MPI_Offset  nelems, *start, *count;                                 \
                                                                        \
    *reqid = NC_REQ_NULL;                                               \
    CHECK_NCID                                                          \
    if (NC_indef(ncp)) return NC_EINDEFINE;                             \
    CHECK_VARID(varid, varp)                                            \
    GET_TOTAL_NUM_ELEMENTS                                              \
    GET_FULL_DIMENSIONS                                                 \
                                                                        \
    /* iget_var is a special case of iget_vara */                       \
    status = ncmpii_igetput_varm(ncp, varp, start, count, NULL, NULL,   \
                                 ip, nelems, datatype, reqid,           \
                                 READ_REQ);                             \
    if (varp->ndims > 0) NCI_Free(start);                               \
                                                                        \
    return status;

/*----< ncmpi_iget_var_uchar() >---------------------------------------------*/
int
ncmpi_iget_var_uchar(int            ncid,
                     int            varid,
                     unsigned char *ip,
                     int           *reqid)
{
    IGET_VAR_COMMON(MPI_UNSIGNED_CHAR);
}

/*----< ncmpi_iget_var_schar() >---------------------------------------------*/
int
ncmpi_iget_var_schar(int          ncid,
                     int          varid,
                     signed char *ip,
                     int         *reqid)
{
    IGET_VAR_COMMON(MPI_BYTE);
}

/*----< ncmpi_iget_var_text() >----------------------------------------------*/
int
ncmpi_iget_var_text(int   ncid,
                    int   varid,
                    char *ip,
                    int  *reqid)
{
    IGET_VAR_COMMON(MPI_CHAR);
}

/*----< ncmpi_iget_var_short() >---------------------------------------------*/
int
ncmpi_iget_var_short(int    ncid,
                     int    varid,
                     short *ip,
                     int   *reqid)
{
    IGET_VAR_COMMON(MPI_SHORT);
}

/*----< ncmpi_iget_var_int() >-----------------------------------------------*/
int
ncmpi_iget_var_int(int  ncid,
                   int  varid,
                   int *ip,
                   int *reqid)
{
    IGET_VAR_COMMON(MPI_INT);
} 

/*----< ncmpi_iget_var_long() >----------------------------------------------*/
int
ncmpi_iget_var_long(int   ncid,
                    int   varid,
                    long *ip,
                    int  *reqid)
{
    IGET_VAR_COMMON(MPI_LONG);
}

/*----< ncmpi_iget_var_float() >---------------------------------------------*/
int
ncmpi_iget_var_float(int    ncid,
                     int    varid,
                     float *ip,
                     int   *reqid)
{
    IGET_VAR_COMMON(MPI_FLOAT);
} 

/*----< ncmpi_iget_var_double() >--------------------------------------------*/
int
ncmpi_iget_var_double(int     ncid,
                      int     varid,
                      double *ip,
                      int    *reqid)
{
    IGET_VAR_COMMON(MPI_DOUBLE);
} 

