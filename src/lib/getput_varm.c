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

/* varm: there maybe two layer of memory layout (remapping):
         one is specified by MPI derived datatype,
         the other is specified by imap[],
         it's encouraged to use only one option of them,
         though using both of them are supported.

   user buffer:                         |--------------------------|

   mpi derived datatype view:           |------|  |------|  |------|
                
   logic (contig) memory datastream:       |------|------|------|

   imap view:                              |--| |--|    |--| |--|

   contig I/O datastream (internal represent): |--|--|--|--|

   These two layers of memory layout will both be represented in MPI 
   derived datatype, and if double layers of memory layout is used, 
   we need to elimilate the upper one passed in MPI_Datatype parameter
   from the user, by repacking it to logic contig memory datastream view.
*/


/*----< ncmpi_put_varm() >---------------------------------------------------*/
int
ncmpi_put_varm(int               ncid,
               int               varid,
               const MPI_Offset  start[],
               const MPI_Offset  count[],
               const MPI_Offset  stride[],
               const MPI_Offset  imap[],
               const void       *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_INDEP_FH
    CHECK_VARID(varid, varp)

    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap,
                              (void*)buf, bufcount, datatype,
                              WRITE_REQ, INDEP_IO);
}

#define PUT_VARM_COMMON(datatype)                                         \
    int         status;                                                   \
    NC         *ncp;                                                      \
    NC_var     *varp;                                                     \
    MPI_Offset  nelems;                                                   \
                                                                          \
    CHECK_NCID                                                            \
    CHECK_WRITE_PERMISSION                                                \
    if (NC_indef(ncp)) return NC_EINDEFINE;                               \
    CHECK_INDEP_FH                                                        \
    CHECK_VARID(varid, varp)                                              \
    GET_NUM_ELEMENTS                                                      \
                                                                          \
    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap,      \
                              (void*)op, nelems, datatype,                \
                              WRITE_REQ, INDEP_IO);

/*----< ncmpi_put_varm_text() >----------------------------------------------*/
int
ncmpi_put_varm_text(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    const MPI_Offset  imap[],
                    const char       *op)
{
    PUT_VARM_COMMON(MPI_CHAR)
}

/*----< ncmpi_put_varm_schar() >---------------------------------------------*/
int
ncmpi_put_varm_schar(int                ncid,
                     int                varid,
                     const MPI_Offset   start[],
                     const MPI_Offset   count[],
                     const MPI_Offset   stride[],
                     const MPI_Offset   imap[],
                     const signed char *op)
{
    PUT_VARM_COMMON(MPI_BYTE)
}

/*----< ncmpi_put_varm_uchar() >---------------------------------------------*/
int
ncmpi_put_varm_uchar(int                  ncid,
                     int                  varid,
                     const MPI_Offset     start[],
                     const MPI_Offset     count[],
                     const MPI_Offset     stride[],
                     const MPI_Offset     imap[],
                     const unsigned char *op)
{
    PUT_VARM_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_put_varm_short() >---------------------------------------------*/
int
ncmpi_put_varm_short(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     const MPI_Offset  imap[],
                     const short      *op)
{
    PUT_VARM_COMMON(MPI_SHORT)
}

/*----< ncmpi_put_varm_int() >-----------------------------------------------*/
int
ncmpi_put_varm_int(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   const MPI_Offset  stride[],
                   const MPI_Offset  imap[],
                   const int        *op)
{
    PUT_VARM_COMMON(MPI_INT)
}

/*----< ncmpi_put_varm_long() >----------------------------------------------*/
int
ncmpi_put_varm_long(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    const MPI_Offset  imap[],
                    const long       *op)
{
    PUT_VARM_COMMON(MPI_LONG)
}

/*----< ncmpi_put_varm_float() >---------------------------------------------*/
int
ncmpi_put_varm_float(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     const MPI_Offset  imap[],
                     const float      *op)
{
    PUT_VARM_COMMON(MPI_FLOAT)
}

/*----< ncmpi_put_varm_double() >--------------------------------------------*/
int
ncmpi_put_varm_double(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const MPI_Offset  stride[],
                      const MPI_Offset  imap[],
                      const double     *op)
{
    PUT_VARM_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_get_varm() >---------------------------------------------------*/
int
ncmpi_get_varm(int               ncid,
               int               varid,
               const MPI_Offset  start[],
               const MPI_Offset  count[],
               const MPI_Offset  stride[],
               const MPI_Offset  imap[],
               void             *buf,
               MPI_Offset        bufcount,
               MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_INDEP_FH
    CHECK_VARID(varid, varp)

    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap, buf,
                              bufcount, datatype, READ_REQ, INDEP_IO);
}

#define GET_VARM_COMMON(datatype)                                       \
    int         status;                                                 \
    NC         *ncp;                                                    \
    NC_var     *varp;                                                   \
    MPI_Offset  nelems;                                                 \
                                                                        \
    CHECK_NCID                                                          \
    if (NC_indef(ncp)) return NC_EINDEFINE;                             \
    CHECK_INDEP_FH                                                      \
    CHECK_VARID(varid, varp)                                            \
    GET_NUM_ELEMENTS                                                    \
                                                                        \
    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap,    \
                              (void*)ip, nelems, datatype,              \
                              READ_REQ, INDEP_IO);

/*----< ncmpi_get_varm_text() >----------------------------------------------*/
int
ncmpi_get_varm_text(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    const MPI_Offset  imap[],
                    char             *ip)
{
    GET_VARM_COMMON(MPI_CHAR)
}

/*----< ncmpi_get_varm_schar() >---------------------------------------------*/
int
ncmpi_get_varm_schar(int                ncid,
                     int                varid,
                     const MPI_Offset   start[],
                     const MPI_Offset   count[],
                     const MPI_Offset   stride[],
                     const MPI_Offset   imap[],
                     signed char       *ip)
{
    GET_VARM_COMMON(MPI_BYTE)
}

/*----< ncmpi_get_varm_uchar() >---------------------------------------------*/
int
ncmpi_get_varm_uchar(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     const MPI_Offset  imap[],
                     unsigned char    *ip)
{
    GET_VARM_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_get_varm_short() >---------------------------------------------*/
int
ncmpi_get_varm_short(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     const MPI_Offset  imap[],
                     short            *ip)
{
    GET_VARM_COMMON(MPI_SHORT)
}

/*----< ncmpi_get_varm_int() >-----------------------------------------------*/
int
ncmpi_get_varm_int(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   const MPI_Offset  stride[],
                   const MPI_Offset  imap[],
                   int              *ip)
{
    GET_VARM_COMMON(MPI_INT)
}

/*----< ncmpi_get_varm_long() >----------------------------------------------*/
int
ncmpi_get_varm_long(int               ncid,
                    int               varid,
                    const MPI_Offset  start[],
                    const MPI_Offset  count[],
                    const MPI_Offset  stride[],
                    const MPI_Offset  imap[],
                    long             *ip)
{
    GET_VARM_COMMON(MPI_LONG)
}

/*----< ncmpi_get_varm_float() >---------------------------------------------*/
int
ncmpi_get_varm_float(int               ncid,
                     int               varid,
                     const MPI_Offset  start[],
                     const MPI_Offset  count[],
                     const MPI_Offset  stride[],
                     const MPI_Offset  imap[],
                     float            *ip)
{
    GET_VARM_COMMON(MPI_FLOAT)
}

/*----< ncmpi_get_varm_double() >--------------------------------------------*/
int
ncmpi_get_varm_double(int               ncid,
                      int               varid,
                      const MPI_Offset  start[],
                      const MPI_Offset  count[],
                      const MPI_Offset  stride[],
                      const MPI_Offset  imap[],
                      double           *ip)
{
    GET_VARM_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_put_varm_all() >-----------------------------------------------*/
int
ncmpi_put_varm_all(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   const MPI_Offset  stride[],
                   const MPI_Offset  imap[],
                   const void       *buf,
                   MPI_Offset        bufcount,
                   MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    CHECK_NCID
    CHECK_WRITE_PERMISSION
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_COLLECTIVE_FH
    CHECK_VARID(varid, varp)

    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap,
                              (void*)buf, bufcount, datatype,
                              WRITE_REQ, COLL_IO);
}

#define PUT_VARM_ALL_COMMON(datatype)                                  \
    int         status;                                                \
    NC         *ncp;                                                   \
    NC_var     *varp;                                                  \
    MPI_Offset  nelems;                                                \
                                                                       \
    CHECK_NCID                                                         \
    CHECK_WRITE_PERMISSION                                             \
    if (NC_indef(ncp)) return NC_EINDEFINE;                            \
    CHECK_COLLECTIVE_FH                                                \
    CHECK_VARID(varid, varp)                                           \
    GET_NUM_ELEMENTS                                                   \
    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap,   \
                              (void*)op, nelems, datatype,             \
                              WRITE_REQ, COLL_IO);

/*----< ncmpi_put_varm_text_all() >------------------------------------------*/
int
ncmpi_put_varm_text_all(int               ncid,
                        int               varid,
                        const MPI_Offset  start[],
                        const MPI_Offset  count[],
                        const MPI_Offset  stride[],
                        const MPI_Offset  imap[],
                        const char       *op)
{
    PUT_VARM_ALL_COMMON(MPI_CHAR)
}

/*----< ncmpi_put_varm_schar_all() >-----------------------------------------*/
int
ncmpi_put_varm_schar_all(int                ncid,
                         int                varid,
                         const MPI_Offset   start[],
                         const MPI_Offset   count[],
                         const MPI_Offset   stride[],
                         const MPI_Offset   imap[],
                         const signed char *op)
{
    PUT_VARM_ALL_COMMON(MPI_BYTE)
}

/*----< ncmpi_put_varm_uchar_all() >-----------------------------------------*/
int
ncmpi_put_varm_uchar_all(int                  ncid,
                         int                  varid,
                         const MPI_Offset     start[],
                         const MPI_Offset     count[],
                         const MPI_Offset     stride[],
                         const MPI_Offset     imap[],
                         const unsigned char *op)
{
    PUT_VARM_ALL_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_put_varm_short_all() >-----------------------------------------*/
int
ncmpi_put_varm_short_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         const MPI_Offset  stride[],
                         const MPI_Offset  imap[],
                         const short      *op)
{
    PUT_VARM_ALL_COMMON(MPI_SHORT)
}

/*----< ncmpi_put_varm_int_all() >-------------------------------------------*/
int
ncmpi_put_varm_int_all(int               ncid,
                       int               varid,
                       const MPI_Offset  start[],
                       const MPI_Offset  count[],
                       const MPI_Offset  stride[],
                       const MPI_Offset  imap[],
                       const int        *op)
{
    PUT_VARM_ALL_COMMON(MPI_INT)
}

/*----< ncmpi_put_varm_long_all() >------------------------------------------*/
int
ncmpi_put_varm_long_all(int               ncid,
                        int               varid,
                        const MPI_Offset  start[],
                        const MPI_Offset  count[],
                        const MPI_Offset  stride[],
                        const MPI_Offset  imap[],
                        const long       *op)
{
    PUT_VARM_ALL_COMMON(MPI_LONG)
}

/*----< ncmpi_put_varm_float_all() >-----------------------------------------*/
int
ncmpi_put_varm_float_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         const MPI_Offset  stride[],
                         const MPI_Offset  imap[],
                         const float      *op)
{
    PUT_VARM_ALL_COMMON(MPI_FLOAT)
}

/*----< ncmpi_put_varm_double_all() >----------------------------------------*/
int
ncmpi_put_varm_double_all(int               ncid,
                          int               varid,
                          const MPI_Offset  start[],
                          const MPI_Offset  count[],
                          const MPI_Offset  stride[],
                          const MPI_Offset  imap[],
                          const double     *op)
{
    PUT_VARM_ALL_COMMON(MPI_DOUBLE)
}

/*----< ncmpi_get_varm_all() >-----------------------------------------------*/
int
ncmpi_get_varm_all(int               ncid,
                   int               varid,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   const MPI_Offset  stride[],
                   const MPI_Offset  imap[],
                   void             *buf,
                   MPI_Offset        bufcount,
                   MPI_Datatype      datatype)
{
    int     status;
    NC     *ncp;
    NC_var *varp;

    CHECK_NCID
    if (NC_indef(ncp)) return NC_EINDEFINE;
    CHECK_COLLECTIVE_FH
    CHECK_VARID(varid, varp)

    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap, buf,
                              bufcount, datatype, READ_REQ, COLL_IO);
}

#define GET_VARM_ALL_COMMON(datatype)                                 \
    int         status;                                               \
    NC         *ncp;                                                  \
    NC_var     *varp;                                                 \
    MPI_Offset  nelems;                                               \
                                                                      \
    CHECK_NCID                                                        \
    if (NC_indef(ncp)) return NC_EINDEFINE;                           \
    CHECK_COLLECTIVE_FH                                               \
    CHECK_VARID(varid, varp)                                          \
    GET_NUM_ELEMENTS                                                  \
                                                                      \
    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap,  \
                              (void*)ip, nelems,datatype,             \
                              READ_REQ, COLL_IO);

/*----< ncmpi_get_varm_text_all() >------------------------------------------*/
int
ncmpi_get_varm_text_all(int               ncid,
                        int               varid,
                        const MPI_Offset  start[],
                        const MPI_Offset  count[],
                        const MPI_Offset  stride[],
                        const MPI_Offset  imap[],
                        char             *ip)
{
    GET_VARM_ALL_COMMON(MPI_CHAR)
}

/*----< ncmpi_get_varm_schar_all() >-----------------------------------------*/
int
ncmpi_get_varm_schar_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         const MPI_Offset  stride[],
                         const MPI_Offset  imap[],
                         signed char      *ip)
{
    GET_VARM_ALL_COMMON(MPI_BYTE)
}

/*----< ncmpi_get_varm_uchar_all() >-----------------------------------------*/
int
ncmpi_get_varm_uchar_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         const MPI_Offset  stride[],
                         const MPI_Offset  imap[],
                         unsigned char    *ip)
{
    GET_VARM_ALL_COMMON(MPI_UNSIGNED_CHAR)
}

/*----< ncmpi_get_varm_short_all() >-----------------------------------------*/
int
ncmpi_get_varm_short_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         const MPI_Offset  stride[],
                         const MPI_Offset  imap[],
                         short            *ip)
{
    GET_VARM_ALL_COMMON(MPI_SHORT)
}

/*----< ncmpi_get_varm_int_all() >-------------------------------------------*/
int
ncmpi_get_varm_int_all(int               ncid,
                       int               varid,
                       const MPI_Offset  start[],
                       const MPI_Offset  count[],
                       const MPI_Offset  stride[],
                       const MPI_Offset  imap[],
                       int              *ip)
{
    GET_VARM_ALL_COMMON(MPI_INT)
}

/*----< ncmpi_get_varm_long_all() >------------------------------------------*/
int
ncmpi_get_varm_long_all(int               ncid,
                        int               varid,
                        const MPI_Offset  start[],
                        const MPI_Offset  count[],
                        const MPI_Offset  stride[],
                        const MPI_Offset  imap[],
                        long             *ip)
{
    GET_VARM_ALL_COMMON(MPI_LONG)
}

/*----< ncmpi_get_varm_float_all() >-----------------------------------------*/
int
ncmpi_get_varm_float_all(int               ncid,
                         int               varid,
                         const MPI_Offset  start[],
                         const MPI_Offset  count[],
                         const MPI_Offset  stride[],
                         const MPI_Offset  imap[],
                         float            *ip)
{
    GET_VARM_ALL_COMMON(MPI_FLOAT)
}

/*----< ncmpi_get_varm_double_all() >----------------------------------------*/
int
ncmpi_get_varm_double_all(int               ncid,
                          int               varid,
                          const MPI_Offset  start[],
                          const MPI_Offset  count[],
                          const MPI_Offset  stride[],
                          const MPI_Offset  imap[],
                          double           *ip)
{
    GET_VARM_ALL_COMMON(MPI_DOUBLE)
}

/*----< ncmpii_getput_varm() >-----------------------------------------------*/
int
ncmpii_getput_varm(NC               *ncp,
                   NC_var           *varp,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   const MPI_Offset  stride[],
                   const MPI_Offset  imap[],
                   void             *buf,
                   MPI_Offset        bufcount,
                   MPI_Datatype      datatype,
                   int               rw_flag,    /* WRITE_REQ or READ_REQ */
                   int               io_method)  /* COLL_IO or INDEP_IO */
{
    void *lbuf=NULL, *cbuf=NULL;
    int err, status, warning; /* err is for API abort and status is not */
    int dim, imap_contig_blocklen, el_size, iscontig_of_ptypes;
    MPI_Offset lnelems, cnelems;
    MPI_Datatype ptype, tmptype, imaptype;

    /* "API error" will abort this API call, but not the entire program */
    err = status = warning = NC_NOERR;

    if (imap == NULL || varp->ndims == 0) {
        /* when imap == NULL, no mapping, same as vars.
           when varp->ndims == 0, reduced to scalar var, only one value
           at one fixed place */
        status = ncmpii_getput_vars(ncp, varp, start, count, stride, buf,
                                    bufcount, datatype, rw_flag, io_method);
        return status;
    }

    imap_contig_blocklen = 1;
    dim = varp->ndims;
    /* test each dim's contiguity until the 1st non-contiguous dim is reached */
    while ( --dim >= 0 && imap_contig_blocklen == imap[dim] ) {
        if (count[dim] < 0) { /* API error */
            err = NC_ENEGATIVECNT;
            goto err_check;
        }
        imap_contig_blocklen *= count[dim];
    }

    if (dim == -1) { /* imap is a contiguous layout */
        status = ncmpii_getput_vars(ncp, varp, start, count, stride, buf,
                                    bufcount, datatype, rw_flag, io_method);
        return status;
    }
    /* else imap gives non-contiguous layout, and need pack/unpack */

    CHECK_DATATYPE(datatype, ptype, el_size, lnelems, iscontig_of_ptypes)

    if (!iscontig_of_ptypes) {
        /* handling for derived datatype: pack into a contiguous buffer */
        lnelems *= bufcount;
        lbuf = NCI_Malloc(lnelems*el_size);
        status = ncmpii_data_repack((void*)buf, bufcount, datatype,
                                    lbuf, lnelems, ptype);
        if (status != NC_NOERR) { /* API error */
            NCI_Free(lbuf);
            err = ((warning != NC_NOERR) ? warning : status);
            goto err_check;
        }
    } else {
        lbuf = (void*)buf;
    }

    if (count[dim] < 0) { /* API error */
        if (!iscontig_of_ptypes && lbuf != NULL)
            NCI_Free(lbuf);
        err = ((warning != NC_NOERR) ? warning : NC_ENEGATIVECNT);
        goto err_check;
    }
    MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
                    ptype, &imaptype);
    MPI_Type_commit(&imaptype);
    cnelems = imap_contig_blocklen*count[dim];
    for (dim--; dim>=0; dim--) {
        if (count[dim] < 0) { /* API error */
            if (!iscontig_of_ptypes && lbuf != NULL)
                NCI_Free(lbuf);
            err = ((warning != NC_NOERR) ? warning : NC_ENEGATIVECNT);
            goto err_check;
        }
#if (MPI_VERSION < 2)
        MPI_Type_hvector(count[dim], 1, imap[dim]*el_size, imaptype, &tmptype);
#else
        MPI_Type_create_hvector(count[dim], 1, (MPI_Aint)imap[dim]*el_size,
                                imaptype, &tmptype);
#endif
        MPI_Type_free(&imaptype);
        MPI_Type_commit(&tmptype);
        imaptype = tmptype;
        cnelems *= count[dim];
    }

    cbuf = (void*) NCI_Malloc(cnelems*el_size);

err_check:
#define CHECK_FATAL_ERROR_COLLECTIVELY
#ifdef CHECK_FATAL_ERROR_COLLECTIVELY
    /* check API error from any proc before going into a collective call */
    if (io_method == COLL_IO) {
        int global_err;
        MPI_Allreduce(&err, &global_err, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
        if (global_err != NC_NOERR) return err;
    }
    else
#endif
        if (err != NC_NOERR) return err;

    if (rw_flag == READ_REQ) {
        status = ncmpii_getput_vars(ncp, varp, start, count, stride, cbuf,
                                    cnelems, ptype, rw_flag, io_method);
        if (status != NC_NOERR) {
            NCI_Free(cbuf);
            if (!iscontig_of_ptypes && lbuf != NULL)
                NCI_Free(lbuf);
            if (status == NC_ERANGE && warning == NC_NOERR)
                warning = status; /* to satisfy the nc_test logic */
            else
                return ((warning != NC_NOERR) ? warning : status);
        }

        /* layout cbuf to lbuf based on imap */
        status = ncmpii_data_repack(cbuf, cnelems, ptype, lbuf, 1, imaptype);
        if (status != NC_NOERR)
            return ((warning != NC_NOERR) ? warning : status);

        if (!iscontig_of_ptypes) {
            /* repack it back, like a read-modify-write operation */
            status = ncmpii_data_repack(lbuf, lnelems, ptype,
                                        (void *)buf, bufcount, datatype);
            if (status != NC_NOERR)
                return ((warning != NC_NOERR) ? warning : status);
        }
    } else { /* WRITE_REQ */
        /* layout lbuf to cbuf based on imap */
        status = ncmpii_data_repack(lbuf, 1, imaptype, cbuf, cnelems, ptype);
        if (status != NC_NOERR) {
            NCI_Free(cbuf);
            if (!iscontig_of_ptypes && lbuf != NULL)
                NCI_Free(lbuf);
            return ((warning != NC_NOERR) ? warning : status);
        }

        status = ncmpii_getput_vars(ncp, varp, start, count, stride, cbuf,
                                    cnelems, ptype, rw_flag, io_method);
    }

    if (imaptype != MPI_DATATYPE_NULL)
        MPI_Type_free(&imaptype);

    if (!iscontig_of_ptypes && lbuf != NULL)
        NCI_Free(lbuf);

    if (cbuf != NULL)
        NCI_Free(cbuf);

    return ((warning != NC_NOERR) ? warning : status);
}
