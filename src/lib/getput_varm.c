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

#define PUT_VARM(fnmode, collmode)                                   \
int                                                                  \
ncmpi_put_varm##fnmode(int               ncid,                       \
                       int               varid,                      \
                       const MPI_Offset  start[],                    \
                       const MPI_Offset  count[],                    \
                       const MPI_Offset  stride[],                   \
                       const MPI_Offset  imap[],                     \
                       const void       *buf,                        \
                       MPI_Offset        bufcount,                   \
                       MPI_Datatype      datatype)                   \
{                                                                    \
    int         status;                                              \
    NC         *ncp;                                                 \
    NC_var     *varp;                                                \
                                                                     \
    CHECK_NCID                                                       \
    CHECK_WRITE_PERMISSION                                           \
    if (NC_indef(ncp)) return NC_EINDEFINE;                          \
    if (collmode == INDEP_IO)                                        \
        CHECK_INDEP_FH                                               \
    else /* collmode == COLL_IO */                                   \
        CHECK_COLLECTIVE_FH                                          \
    CHECK_VARID(varid, varp)                                         \
                                                                     \
    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap, \
                              (void*)buf, bufcount, datatype,        \
                              WRITE_REQ, collmode);                  \
}

/*----< ncmpi_put_varm() >---------------------------------------------------*/
/*----< ncmpi_put_varm_all() >-----------------------------------------------*/
PUT_VARM(    , INDEP_IO)
PUT_VARM(_all, COLL_IO)

#define GET_VARM(fnmode, collmode)                                   \
int                                                                  \
ncmpi_get_varm##fnmode(int               ncid,                       \
                       int               varid,                      \
                       const MPI_Offset  start[],                    \
                       const MPI_Offset  count[],                    \
                       const MPI_Offset  stride[],                   \
                       const MPI_Offset  imap[],                     \
                       void             *buf,                        \
                       MPI_Offset        bufcount,                   \
                       MPI_Datatype      datatype)                   \
{                                                                    \
    int         status;                                              \
    NC         *ncp;                                                 \
    NC_var     *varp;                                                \
                                                                     \
    CHECK_NCID                                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                          \
    if (collmode == INDEP_IO)                                        \
        CHECK_INDEP_FH                                               \
    else /* collmode == COLL_IO */                                   \
        CHECK_COLLECTIVE_FH                                          \
    CHECK_VARID(varid, varp)                                         \
                                                                     \
    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap, \
                              buf, bufcount, datatype,               \
                              READ_REQ, collmode);                   \
}

/*----< ncmpi_get_varm() >---------------------------------------------------*/
/*----< ncmpi_get_varm_all() >-----------------------------------------------*/
GET_VARM(    , INDEP_IO)
GET_VARM(_all, COLL_IO)

#define PUT_VARM_TYPE(fntype, buftype, mpitype, collmode)            \
int                                                                  \
ncmpi_put_varm_##fntype(int               ncid,                      \
                        int               varid,                     \
                        const MPI_Offset  start[],                   \
                        const MPI_Offset  count[],                   \
                        const MPI_Offset  stride[],                  \
                        const MPI_Offset  imap[],                    \
                        const buftype    *op)                        \
{                                                                    \
    int         status;                                              \
    NC         *ncp;                                                 \
    NC_var     *varp;                                                \
    MPI_Offset  nelems;                                              \
                                                                     \
    CHECK_NCID                                                       \
    CHECK_WRITE_PERMISSION                                           \
    if (NC_indef(ncp)) return NC_EINDEFINE;                          \
    if (collmode == INDEP_IO)                                        \
        CHECK_INDEP_FH                                               \
    else /* collmode == COLL_IO */                                   \
        CHECK_COLLECTIVE_FH                                          \
    CHECK_VARID(varid, varp)                                         \
    GET_NUM_ELEMENTS                                                 \
                                                                     \
    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap, \
                              (void*)op, nelems, mpitype,            \
                              WRITE_REQ, collmode);                  \
}

/*----< ncmpi_put_varm_text() >----------------------------------------------*/
/*----< ncmpi_put_varm_schar() >---------------------------------------------*/
/*----< ncmpi_put_varm_uchar() >---------------------------------------------*/
/*----< ncmpi_put_varm_short() >---------------------------------------------*/
/*----< ncmpi_put_varm_int() >-----------------------------------------------*/
/*----< ncmpi_put_varm_long() >----------------------------------------------*/
/*----< ncmpi_put_varm_float() >---------------------------------------------*/
/*----< ncmpi_put_varm_double() >--------------------------------------------*/
/*----< ncmpi_get_varm() >---------------------------------------------------*/

PUT_VARM_TYPE(text,   char,   MPI_CHAR,              INDEP_IO)
PUT_VARM_TYPE(schar,  schar,  MPI_BYTE,              INDEP_IO)
PUT_VARM_TYPE(uchar,  uchar,  MPI_UNSIGNED_CHAR,     INDEP_IO)
PUT_VARM_TYPE(short,  short,  MPI_SHORT,             INDEP_IO)
PUT_VARM_TYPE(int,    int,    MPI_INT,               INDEP_IO)
PUT_VARM_TYPE(long,   long,   MPI_LONG,              INDEP_IO)
PUT_VARM_TYPE(float,  float,  MPI_FLOAT,             INDEP_IO)
PUT_VARM_TYPE(double, double, MPI_DOUBLE,            INDEP_IO)

/*----< ncmpi_put_varm_text_all() >------------------------------------------*/
/*----< ncmpi_put_varm_schar_all() >-----------------------------------------*/
/*----< ncmpi_put_varm_uchar_all() >-----------------------------------------*/
/*----< ncmpi_put_varm_short_all() >-----------------------------------------*/
/*----< ncmpi_put_varm_int_all() >-------------------------------------------*/
/*----< ncmpi_put_varm_long_all() >------------------------------------------*/
/*----< ncmpi_put_varm_float_all() >-----------------------------------------*/
/*----< ncmpi_put_varm_double_all() >----------------------------------------*/

PUT_VARM_TYPE(text_all,   char,   MPI_CHAR,          COLL_IO)
PUT_VARM_TYPE(schar_all,  schar,  MPI_BYTE,          COLL_IO)
PUT_VARM_TYPE(uchar_all,  uchar,  MPI_UNSIGNED_CHAR, COLL_IO)
PUT_VARM_TYPE(short_all,  short,  MPI_SHORT,         COLL_IO)
PUT_VARM_TYPE(int_all,    int,    MPI_INT,           COLL_IO)
PUT_VARM_TYPE(long_all,   long,   MPI_LONG,          COLL_IO)
PUT_VARM_TYPE(float_all,  float,  MPI_FLOAT,         COLL_IO)
PUT_VARM_TYPE(double_all, double, MPI_DOUBLE,        COLL_IO)


#define GET_VARM_TYPE(fntype, buftype, mpitype, collmode)            \
int                                                                  \
ncmpi_get_varm_##fntype(int               ncid,                      \
                        int               varid,                     \
                        const MPI_Offset  start[],                   \
                        const MPI_Offset  count[],                   \
                        const MPI_Offset  stride[],                  \
                        const MPI_Offset  imap[],                    \
                        buftype          *ip)                        \
{                                                                    \
    int         status;                                              \
    NC         *ncp;                                                 \
    NC_var     *varp;                                                \
    MPI_Offset  nelems;                                              \
                                                                     \
    CHECK_NCID                                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                          \
    if (collmode == INDEP_IO)                                        \
        CHECK_INDEP_FH                                               \
    else /* collmode == COLL_IO */                                   \
        CHECK_COLLECTIVE_FH                                          \
    CHECK_VARID(varid, varp)                                         \
    GET_NUM_ELEMENTS                                                 \
                                                                     \
    return ncmpii_getput_varm(ncp, varp, start, count, stride, imap, \
                              ip, nelems, mpitype,                   \
                              READ_REQ, collmode);                   \
}

/*----< ncmpi_get_varm_text() >----------------------------------------------*/
/*----< ncmpi_get_varm_schar() >---------------------------------------------*/
/*----< ncmpi_get_varm_uchar() >---------------------------------------------*/
/*----< ncmpi_get_varm_short() >---------------------------------------------*/
/*----< ncmpi_get_varm_int() >-----------------------------------------------*/
/*----< ncmpi_get_varm_long() >----------------------------------------------*/
/*----< ncmpi_get_varm_float() >---------------------------------------------*/
/*----< ncmpi_get_varm_double() >--------------------------------------------*/

GET_VARM_TYPE(text,   char,   MPI_CHAR,              INDEP_IO)
GET_VARM_TYPE(schar,  schar,  MPI_BYTE,              INDEP_IO)
GET_VARM_TYPE(uchar,  uchar,  MPI_UNSIGNED_CHAR,     INDEP_IO)
GET_VARM_TYPE(short,  short,  MPI_SHORT,             INDEP_IO)
GET_VARM_TYPE(int,    int,    MPI_INT,               INDEP_IO)
GET_VARM_TYPE(long,   long,   MPI_LONG,              INDEP_IO)
GET_VARM_TYPE(float,  float,  MPI_FLOAT,             INDEP_IO)
GET_VARM_TYPE(double, double, MPI_DOUBLE,            INDEP_IO)

/*----< ncmpi_get_varm_text_all() >------------------------------------------*/
/*----< ncmpi_get_varm_schar_all() >-----------------------------------------*/
/*----< ncmpi_get_varm_uchar_all() >-----------------------------------------*/
/*----< ncmpi_get_varm_short_all() >-----------------------------------------*/
/*----< ncmpi_get_varm_int_all() >-------------------------------------------*/
/*----< ncmpi_get_varm_long_all() >------------------------------------------*/
/*----< ncmpi_get_varm_float_all() >-----------------------------------------*/
/*----< ncmpi_get_varm_double_all() >----------------------------------------*/

GET_VARM_TYPE(text_all,   char,   MPI_CHAR,          COLL_IO)
GET_VARM_TYPE(schar_all,  schar,  MPI_BYTE,          COLL_IO)
GET_VARM_TYPE(uchar_all,  uchar,  MPI_UNSIGNED_CHAR, COLL_IO)
GET_VARM_TYPE(short_all,  short,  MPI_SHORT,         COLL_IO)
GET_VARM_TYPE(int_all,    int,    MPI_INT,           COLL_IO)
GET_VARM_TYPE(long_all,   long,   MPI_LONG,          COLL_IO)
GET_VARM_TYPE(float_all,  float,  MPI_FLOAT,         COLL_IO)
GET_VARM_TYPE(double_all, double, MPI_DOUBLE,        COLL_IO)

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
    imaptype = MPI_DATATYPE_NULL;

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
        if (rw_flag == WRITE_REQ) { /* only write needs this packing */
            status = ncmpii_data_repack((void*)buf, bufcount, datatype,
                                        lbuf, lnelems, ptype);
            if (status != NC_NOERR) { /* API error */
                err = ((warning != NC_NOERR) ? warning : status);
                goto err_check;
            }
        }
    } else {
        lbuf = (void*)buf;
    }

    if (count[dim] < 0) { /* API error */
        err = ((warning != NC_NOERR) ? warning : NC_ENEGATIVECNT);
        goto err_check;
    }
    MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
                    ptype, &imaptype);
    MPI_Type_commit(&imaptype);
    cnelems = imap_contig_blocklen*count[dim];
    for (dim--; dim>=0; dim--) {
        if (count[dim] < 0) { /* API error */
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
    /* check API error from any proc before going into a collective call.
     * optimization: to avoid MPI_Allreduce to check parameters at
     * every call, we assume caller does the right thing most of the
     * time.  If caller passed in bad parameters, we'll still conduct a
     * zero-byte operation (everyone has to participate in the
     * collective I/O call) but return error */
    if (err != NC_NOERR) {
        if (io_method == COLL_IO) {
            MPI_Offset *zeros;
            zeros = (MPI_Offset *) NCI_Calloc(varp->ndims, sizeof(MPI_Offset));
            ncmpii_getput_vars(ncp, varp, zeros, zeros, NULL, buf,
                               0, MPI_BYTE, rw_flag, io_method);
            NCI_Free(zeros);
        }
        if (lbuf != NULL && lbuf != buf)
            NCI_Free(lbuf);
        if (imaptype != MPI_DATATYPE_NULL)
            MPI_Type_free(&imaptype);
        return err;
    }

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
