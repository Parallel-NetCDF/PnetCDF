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

#include "macro.h"


/* buffer layers:       
        
        User Level              buf     (user defined buffer of MPI_Datatype)
        MPI Datatype Level      cbuf    (contiguous buffer of ptype)
        NetCDF XDR Level        xbuf    (XDR I/O buffer)
*/

/*----< ncmpi_put_var() >----------------------------------------------------*/
int
ncmpi_put_var(int           ncid,
              int           varid,
              const void   *buf,
              MPI_Offset    bufcount,
              MPI_Datatype  datatype)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *start, *count;

    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_INDEP_FH
    CHECK_VARID(varid, varp)
    GET_FULL_DIMENSIONS

    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,
                                (void*)buf, bufcount, datatype,
                                WRITE_REQ, INDEP_IO);
    if (varp->ndims > 0) NCI_Free(start);

    return status;
}

#define PUT_VAR_COMMON(datatype)                                 \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems, *start, *count;                          \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    CHECK_INDEP_FH                                               \
    CHECK_VARID(varid, varp)                                     \
    GET_TOTAL_NUM_ELEMENTS                                       \
    GET_FULL_DIMENSIONS                                          \
                                                                 \
    /* put_var is a special case of put_vara */                  \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                (void*)op, nelems, datatype,     \
                                WRITE_REQ, INDEP_IO);            \
    if (varp->ndims > 0) NCI_Free(start);                        \
                                                                 \
    return status;

/*----< ncmpi_put_var_text() >-----------------------------------------------*/
int
ncmpi_put_var_text(int         ncid,
                   int         varid,
                   const char *op)
{
    PUT_VAR_COMMON(MPI_CHAR)
}

/*----< ncmpi_put_var_schar() >----------------------------------------------*/
int
ncmpi_put_var_schar(int                ncid,
                    int                varid,
                    const signed char *op)
{
    PUT_VAR_COMMON(MPI_BYTE)
}

/*----< ncmpi_put_var_uchar() >----------------------------------------------*/
int
ncmpi_put_var_uchar(int                  ncid,
                    int                  varid,
                    const unsigned char *op)
{
    PUT_VAR_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_put_var_short() >----------------------------------------------*/
int
ncmpi_put_var_short(int          ncid,
                    int          varid,
                    const short *op)
{
    PUT_VAR_COMMON(MPI_SHORT)
}

/*----< ncmpi_put_var_int() >------------------------------------------------*/
int
ncmpi_put_var_int(int        ncid,
                  int        varid,
                  const int *op)
{
    PUT_VAR_COMMON(MPI_INT)
}

/*----< ncmpi_put_var_long() >-----------------------------------------------*/
int
ncmpi_put_var_long(int         ncid,
                   int         varid,
                   const long *op)
{
    PUT_VAR_COMMON(MPI_LONG)
}

/*----< ncmpi_put_var_float() >----------------------------------------------*/
int
ncmpi_put_var_float(int          ncid,
                    int          varid,
                    const float *op)
{
    PUT_VAR_COMMON(MPI_FLOAT)
}

/*----< ncmpi_put_var_double() >---------------------------------------------*/
int
ncmpi_put_var_double(int           ncid,
                     int           varid,
                     const double *op)
{
    PUT_VAR_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_get_var() >----------------------------------------------------*/
int
ncmpi_get_var(int           ncid,
              int           varid,
              void         *buf,
              MPI_Offset    bufcount,
              MPI_Datatype  datatype)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *start, *count;

    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_INDEP_FH
    CHECK_VARID(varid, varp)
    GET_FULL_DIMENSIONS

    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,
                                buf, bufcount, datatype,
                                READ_REQ, INDEP_IO);
    if (varp->ndims > 0) NCI_Free(start);

    return status;
}

#define GET_VAR_COMMON(datatype)                                 \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems, *start, *count;                          \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    CHECK_INDEP_FH                                               \
    CHECK_VARID(varid, varp)                                     \
    GET_TOTAL_NUM_ELEMENTS                                       \
    GET_FULL_DIMENSIONS                                          \
                                                                 \
    /* get_var is a special case of get_vars */                  \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                ip, nelems, datatype,            \
                                READ_REQ, INDEP_IO);             \
    if (varp->ndims > 0) NCI_Free(start);                        \
                                                                 \
    return status;


/*----< ncmpi_get_var_text() >-----------------------------------------------*/
int
ncmpi_get_var_text(int   ncid,
                   int   varid,
                   char *ip)
{
    GET_VAR_COMMON(MPI_CHAR)
}

/*----< ncmpi_get_var_schar() >----------------------------------------------*/
int
ncmpi_get_var_schar(int          ncid,
                    int          varid,
                    signed char *ip)
{
    GET_VAR_COMMON(MPI_BYTE)
}

/*----< ncmpi_get_var_uchar() >----------------------------------------------*/
int
ncmpi_get_var_uchar(int            ncid,
                    int            varid,
                    unsigned char *ip)
{
    GET_VAR_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_get_var_short() >----------------------------------------------*/
int
ncmpi_get_var_short(int    ncid,
                    int    varid,
                    short *ip)
{
    GET_VAR_COMMON(MPI_SHORT)
}

/*----< ncmpi_get_var_int() >------------------------------------------------*/
int
ncmpi_get_var_int(int  ncid,
                  int  varid,
                  int *ip)
{
    GET_VAR_COMMON(MPI_INT)
}

/*----< ncmpi_get_var_long() >-----------------------------------------------*/
int
ncmpi_get_var_long(int   ncid,
                   int   varid,
                   long *ip)
{
    GET_VAR_COMMON(MPI_LONG)
}

/*----< ncmpi_get_var_float() >----------------------------------------------*/
int
ncmpi_get_var_float(int    ncid,
                    int    varid,
                    float *ip)
{
    GET_VAR_COMMON(MPI_FLOAT)
}

/*----< ncmpi_get_var_double() >---------------------------------------------*/
int
ncmpi_get_var_double(int     ncid,
                     int     varid,
                     double *ip)
{
    GET_VAR_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_put_var_all() >------------------------------------------------*/
int
ncmpi_put_var_all(int           ncid,
                  int           varid,
                  const void   *buf,
                  MPI_Offset    bufcount,
                  MPI_Datatype  datatype)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *start, *count;

    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_COLLECTIVE_FH
    CHECK_VARID(varid, varp)
    GET_FULL_DIMENSIONS

    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,
                                (void*)buf, bufcount, datatype,
                                WRITE_REQ, COLL_IO);
    if (varp->ndims > 0) NCI_Free(start);

    return status;
}

#define PUT_VAR_ALL_COMMON(datatype)                              \
    int         status;                                           \
    NC         *ncp;                                              \
    NC_var     *varp;                                             \
    MPI_Offset  nelems, *start, *count;                           \
                                                                  \
    CHECK_NCID                                                    \
    CHECK_WRITE_PERMISSION                                        \
    if (NC_indef(ncp)) return NC_EINDEFINE;                       \
    CHECK_COLLECTIVE_FH                                           \
    CHECK_VARID(varid, varp)                                      \
    GET_TOTAL_NUM_ELEMENTS                                        \
    GET_FULL_DIMENSIONS                                           \
                                                                  \
    /* put_var is a special case of put_vars */                   \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,    \
                                (void*)op, nelems, datatype,      \
                                WRITE_REQ, COLL_IO);              \
    if (varp->ndims > 0) NCI_Free(start);                         \
                                                                  \
    return status;


/*----< ncmpi_put_var_text_all() >-------------------------------------------*/
int
ncmpi_put_var_text_all(int         ncid,
                       int         varid,
                       const char *op)
{
    PUT_VAR_ALL_COMMON(MPI_CHAR)
}

/*----< ncmpi_put_var_schar_all() >------------------------------------------*/
int
ncmpi_put_var_schar_all(int                ncid,
                        int                varid,
                        const signed char *op)
{
    PUT_VAR_ALL_COMMON(MPI_BYTE)
}

/*----< ncmpi_put_var_uchar_all() >------------------------------------------*/
int
ncmpi_put_var_uchar_all(int                  ncid,
                        int                  varid,
                        const unsigned char *op)
{
    PUT_VAR_ALL_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_put_var_short_all() >------------------------------------------*/
int
ncmpi_put_var_short_all(int          ncid,
                        int          varid,
                        const short *op)
{
    PUT_VAR_ALL_COMMON(MPI_SHORT)
}

/*----< ncmpi_put_var_int_all() >--------------------------------------------*/
int
ncmpi_put_var_int_all(int        ncid,
                      int        varid,
                      const int *op)
{
    PUT_VAR_ALL_COMMON(MPI_INT)
}

/*----< ncmpi_put_var_long_all() >-------------------------------------------*/
int
ncmpi_put_var_long_all(int         ncid,
                       int         varid,
                       const long *op)
{
    PUT_VAR_ALL_COMMON(MPI_LONG)
}

/*----< ncmpi_put_var_float_all() >------------------------------------------*/
int
ncmpi_put_var_float_all(int          ncid,
                        int          varid,
                        const float *op)
{
    PUT_VAR_ALL_COMMON(MPI_FLOAT)
}

/*----< ncmpi_put_var_double_all() >-----------------------------------------*/
int
ncmpi_put_var_double_all(int           ncid,
                         int           varid,
                         const double *op)
{
    PUT_VAR_ALL_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_get_var_all() >------------------------------------------------*/
int
ncmpi_get_var_all(int           ncid,
                  int           varid,
                  void         *buf,
                  MPI_Offset    bufcount,
                  MPI_Datatype  datatype)
{
    int         status;
    NC         *ncp;
    NC_var     *varp;
    MPI_Offset *start, *count;

    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_COLLECTIVE_FH
    CHECK_VARID(varid, varp)
    GET_FULL_DIMENSIONS

    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,
                                buf, bufcount, datatype,
                                READ_REQ, COLL_IO);
    if (varp->ndims > 0) NCI_Free(start);

    return status;
}

#define GET_VAR_ALL_COMMON(datatype)                             \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems, *start, *count;                          \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    CHECK_COLLECTIVE_FH                                          \
    CHECK_VARID(varid, varp)                                     \
    GET_TOTAL_NUM_ELEMENTS                                       \
    GET_FULL_DIMENSIONS                                          \
                                                                 \
    /* get_var is a special case of get_vars */                  \
    status = ncmpii_getput_vars(ncp, varp, start, count, NULL,   \
                                ip, nelems, datatype,            \
                                READ_REQ, COLL_IO);              \
    if (varp->ndims > 0) NCI_Free(start);                        \
                                                                 \
    return status;

/*----< ncmpi_get_var_text_all() >-------------------------------------------*/
int
ncmpi_get_var_text_all(int   ncid,
                       int   varid,
                       char *ip)
{
    GET_VAR_ALL_COMMON(MPI_CHAR)
}

/*----< ncmpi_get_var_schar_all() >------------------------------------------*/
int
ncmpi_get_var_schar_all(int          ncid,
                        int          varid,
                        signed char *ip)
{
    GET_VAR_ALL_COMMON(MPI_BYTE)
}

/*----< ncmpi_get_var_uchar_all() >------------------------------------------*/
int
ncmpi_get_var_uchar_all(int            ncid,
                        int            varid,
                        unsigned char *ip)
{
    GET_VAR_ALL_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_get_var_short_all() >------------------------------------------*/
int
ncmpi_get_var_short_all(int    ncid,
                        int    varid,
                        short *ip)
{
    GET_VAR_ALL_COMMON(MPI_SHORT)
}

/*----< ncmpi_get_var_int_all() >--------------------------------------------*/
int
ncmpi_get_var_int_all(int  ncid,
                      int  varid,
                      int *ip)
{
    GET_VAR_ALL_COMMON(MPI_INT)
}

/*----< ncmpi_get_var_long_all() >-------------------------------------------*/
int
ncmpi_get_var_long_all(int   ncid,
                       int   varid,
                       long *ip)
{
    GET_VAR_ALL_COMMON(MPI_LONG)
}

/*----< ncmpi_get_var_float_all() >------------------------------------------*/
int
ncmpi_get_var_float_all(int    ncid,
                        int    varid,
                        float *ip)
{
    GET_VAR_ALL_COMMON(MPI_FLOAT)
}

/*----< ncmpi_get_var_double_all() >-----------------------------------------*/
int
ncmpi_get_var_double_all(int     ncid,
                         int     varid,
                         double *ip)
{
    GET_VAR_ALL_COMMON(MPI_DOUBLE)
}

