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


/* ftype is the variable's nc_type defined in file, eg. int64
 * btype is the I/O buffer's C data type, eg. long long
 * buftype is I/O bufer's MPI data type, eg. MPI_UNSIGNED_LONG_LONG
 * apitype is data type appeared in the API names, eg. ncmpi_get_vara_longlong
 */

/* buffer layers:       
        
        User Level              buf     (user defined buffer of MPI_Datatype)
        MPI Datatype Level      cbuf    (contiguous buffer of ptype)
        NetCDF XDR Level        xbuf    (XDR I/O buffer)
*/

#define PUT_VARS(iomode, collmode)                               \
int                                                              \
ncmpi_put_vars##iomode(int               ncid,                   \
                       int               varid,                  \
                       const MPI_Offset  start[],                \
                       const MPI_Offset  count[],                \
                       const MPI_Offset  stride[],               \
                       const void       *buf,                    \
                       MPI_Offset        bufcount,               \
                       MPI_Datatype      buftype)                \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
                                                                 \
    return ncmpii_getput_vars(ncp, varp, start, count, stride,   \
                              (void*)buf, bufcount, buftype,     \
                              WRITE_REQ, collmode);              \
}

/*----< ncmpi_put_vars() >---------------------------------------------------*/
/*----< ncmpi_put_vars_all() >-----------------------------------------------*/
PUT_VARS(    , INDEP_IO)
PUT_VARS(_all, COLL_IO)

#define GET_VARS(iomode, collmode)                               \
int                                                              \
ncmpi_get_vars##iomode(int               ncid,                   \
                       int               varid,                  \
                       const MPI_Offset  start[],                \
                       const MPI_Offset  count[],                \
                       const MPI_Offset  stride[],               \
                       void             *buf,                    \
                       MPI_Offset        bufcount,               \
                       MPI_Datatype      buftype)                \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
                                                                 \
    return ncmpii_getput_vars(ncp, varp, start, count, stride,   \
                              buf, bufcount, buftype,            \
                              READ_REQ, collmode);               \
}

/*----< ncmpi_get_vars() >---------------------------------------------------*/
/*----< ncmpi_get_vars_all() >-----------------------------------------------*/
GET_VARS(    , INDEP_IO)
GET_VARS(_all, COLL_IO)


#define PUT_VARS_TYPE(apitype, btype, mpitype, collmode)         \
int                                                              \
ncmpi_put_vars_##apitype(int               ncid,                 \
                         int               varid,                \
                         const MPI_Offset  start[],              \
                         const MPI_Offset  count[],              \
                         const MPI_Offset  stride[],             \
                         const btype      *op)                   \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems;                                          \
                                                                 \
    CHECK_NCID                                                   \
    CHECK_WRITE_PERMISSION                                       \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_NUM_ELEMENTS                                             \
                                                                 \
    return ncmpii_getput_vars(ncp, varp, start, count, stride,   \
                              (void*)op, nelems, mpitype,        \
                              WRITE_REQ, collmode);              \
}

/*----< ncmpi_put_vars_text() >-----------------------------------------------*/
/*----< ncmpi_put_vars_schar() >----------------------------------------------*/
/*----< ncmpi_put_vars_uchar() >----------------------------------------------*/
/*----< ncmpi_put_vars_short() >----------------------------------------------*/
/*----< ncmpi_put_vars_ushort() >---------------------------------------------*/
/*----< ncmpi_put_vars_int() >------------------------------------------------*/
/*----< ncmpi_put_vars_uint() >-----------------------------------------------*/
/*----< ncmpi_put_vars_long() >-----------------------------------------------*/
/*----< ncmpi_put_vars_float() >----------------------------------------------*/
/*----< ncmpi_put_vars_double() >---------------------------------------------*/
/*----< ncmpi_put_vars_longlong() >-------------------------------------------*/
/*----< ncmpi_put_vars_ulonglong() >------------------------------------------*/
PUT_VARS_TYPE(text,      char,               MPI_CHAR,               INDEP_IO)
PUT_VARS_TYPE(schar,     schar,              MPI_BYTE,               INDEP_IO)
PUT_VARS_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR,      INDEP_IO)
PUT_VARS_TYPE(short,     short,              MPI_SHORT,              INDEP_IO)
PUT_VARS_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT,     INDEP_IO)
PUT_VARS_TYPE(int,       int,                MPI_INT,                INDEP_IO)
PUT_VARS_TYPE(uint,      uint,               MPI_UNSIGNED,           INDEP_IO)
PUT_VARS_TYPE(long,      long,               MPI_LONG,               INDEP_IO)
PUT_VARS_TYPE(float,     float,              MPI_FLOAT,              INDEP_IO)
PUT_VARS_TYPE(double,    double,             MPI_DOUBLE,             INDEP_IO)
PUT_VARS_TYPE(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
PUT_VARS_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)
// PUT_VARS_TYPE(string, char*,              MPI_CHAR,               INDEP_IO)
/* string is not yet supported */

/*----< ncmpi_put_vars_text_all() >-------------------------------------------*/
/*----< ncmpi_put_vars_schar_all() >------------------------------------------*/
/*----< ncmpi_put_vars_uchar_all() >------------------------------------------*/
/*----< ncmpi_put_vars_short_all() >------------------------------------------*/
/*----< ncmpi_put_vars_ushort_all() >-----------------------------------------*/
/*----< ncmpi_put_vars_int_all() >--------------------------------------------*/
/*----< ncmpi_put_vars_uint_all() >-------------------------------------------*/
/*----< ncmpi_put_vars_long_all() >-------------------------------------------*/
/*----< ncmpi_put_vars_float_all() >------------------------------------------*/
/*----< ncmpi_put_vars_double_all() >-----------------------------------------*/
/*----< ncmpi_put_vars_longlong_all() >---------------------------------------*/
/*----< ncmpi_put_vars_ulonglong_all() >--------------------------------------*/
PUT_VARS_TYPE(text_all,      char,               MPI_CHAR,              COLL_IO)
PUT_VARS_TYPE(schar_all,     schar,              MPI_BYTE,              COLL_IO)
PUT_VARS_TYPE(uchar_all,     uchar,              MPI_UNSIGNED_CHAR,     COLL_IO)
PUT_VARS_TYPE(short_all,     short,              MPI_SHORT,             COLL_IO)
PUT_VARS_TYPE(ushort_all,    ushort,             MPI_UNSIGNED_SHORT,    COLL_IO)
PUT_VARS_TYPE(int_all,       int,                MPI_INT,               COLL_IO)
PUT_VARS_TYPE(uint_all,      uint,               MPI_UNSIGNED,          COLL_IO)
PUT_VARS_TYPE(long_all,      long,               MPI_LONG,              COLL_IO)
PUT_VARS_TYPE(float_all,     float,              MPI_FLOAT,             COLL_IO)
PUT_VARS_TYPE(double_all,    double,             MPI_DOUBLE,            COLL_IO)
PUT_VARS_TYPE(longlong_all,  long long,          MPI_LONG_LONG_INT,     COLL_IO)
PUT_VARS_TYPE(ulonglong_all, unsigned long long, MPI_UNSIGNED_LONG_LONG,COLL_IO)
// PUT_VARS_TYPE(string_all, char*,              MPI_CHAR,              COLL_IO)
/* string is not yet supported */


#define GET_VARS_TYPE(apitype, btype, mpitype, collmode)         \
int                                                              \
ncmpi_get_vars_##apitype(int               ncid,                 \
                         int               varid,                \
                         const MPI_Offset  start[],              \
                         const MPI_Offset  count[],              \
                         const MPI_Offset  stride[],             \
                         btype            *ip)                   \
{                                                                \
    int         status;                                          \
    NC         *ncp;                                             \
    NC_var     *varp;                                            \
    MPI_Offset  nelems;                                          \
                                                                 \
    CHECK_NCID                                                   \
    if (NC_indef(ncp)) return NC_EINDEFINE;                      \
    if (collmode == INDEP_IO)                                    \
        CHECK_INDEP_FH                                           \
    else /* collmode == COLL_IO */                               \
        CHECK_COLLECTIVE_FH                                      \
    CHECK_VARID(varid, varp)                                     \
    GET_NUM_ELEMENTS                                             \
                                                                 \
    return ncmpii_getput_vars(ncp, varp, start, count, stride,   \
                              ip, nelems, mpitype,               \
                              READ_REQ, collmode);               \
}

/*----< ncmpi_get_vars_text() >-----------------------------------------------*/
/*----< ncmpi_get_vars_schar() >----------------------------------------------*/
/*----< ncmpi_get_vars_uchar() >----------------------------------------------*/
/*----< ncmpi_get_vars_short() >----------------------------------------------*/
/*----< ncmpi_get_vars_ushort() >---------------------------------------------*/
/*----< ncmpi_get_vars_int() >------------------------------------------------*/
/*----< ncmpi_get_vars_uint() >-----------------------------------------------*/
/*----< ncmpi_get_vars_long() >-----------------------------------------------*/
/*----< ncmpi_get_vars_float() >----------------------------------------------*/
/*----< ncmpi_get_vars_double() >---------------------------------------------*/
/*----< ncmpi_get_vars_longlong() >-------------------------------------------*/
/*----< ncmpi_get_vars_ulonglong() >------------------------------------------*/
GET_VARS_TYPE(text,      char,               MPI_CHAR,               INDEP_IO)
GET_VARS_TYPE(schar,     schar,              MPI_BYTE,               INDEP_IO)
GET_VARS_TYPE(uchar,     uchar,              MPI_UNSIGNED_CHAR,      INDEP_IO)
GET_VARS_TYPE(short,     short,              MPI_SHORT,              INDEP_IO)
GET_VARS_TYPE(ushort,    ushort,             MPI_UNSIGNED_SHORT,     INDEP_IO)
GET_VARS_TYPE(int,       int,                MPI_INT,                INDEP_IO)
GET_VARS_TYPE(uint,      uint,               MPI_UNSIGNED,           INDEP_IO)
GET_VARS_TYPE(long,      long,               MPI_LONG,               INDEP_IO)
GET_VARS_TYPE(float,     float,              MPI_FLOAT,              INDEP_IO)
GET_VARS_TYPE(double,    double,             MPI_DOUBLE,             INDEP_IO)
GET_VARS_TYPE(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
GET_VARS_TYPE(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)
// GET_VARS_TYPE(string, char*,              MPI_CHAR,               INDEP_IO)
/* string is not yet supported */

/*----< ncmpi_get_vars_text_all() >-------------------------------------------*/
/*----< ncmpi_get_vars_schar_all() >------------------------------------------*/
/*----< ncmpi_get_vars_uchar_all() >------------------------------------------*/
/*----< ncmpi_get_vars_short_all() >------------------------------------------*/
/*----< ncmpi_get_vars_ushort_all() >-----------------------------------------*/
/*----< ncmpi_get_vars_int_all() >--------------------------------------------*/
/*----< ncmpi_get_vars_uint_all() >-------------------------------------------*/
/*----< ncmpi_get_vars_long_all() >-------------------------------------------*/
/*----< ncmpi_get_vars_float_all() >------------------------------------------*/
/*----< ncmpi_get_vars_double_all() >-----------------------------------------*/
/*----< ncmpi_get_vars_longlong_all() >---------------------------------------*/
/*----< ncmpi_get_vars_ulonglong_all() >--------------------------------------*/
GET_VARS_TYPE(text_all,      char,               MPI_CHAR,              COLL_IO)
GET_VARS_TYPE(schar_all,     schar,              MPI_BYTE,              COLL_IO)
GET_VARS_TYPE(uchar_all,     uchar,              MPI_UNSIGNED_CHAR,     COLL_IO)
GET_VARS_TYPE(short_all,     short,              MPI_SHORT,             COLL_IO)
GET_VARS_TYPE(ushort_all,    ushort,             MPI_UNSIGNED_SHORT,    COLL_IO)
GET_VARS_TYPE(int_all,       int,                MPI_INT,               COLL_IO)
GET_VARS_TYPE(uint_all,      uint,               MPI_UNSIGNED,          COLL_IO)
GET_VARS_TYPE(long_all,      long,               MPI_LONG,              COLL_IO)
GET_VARS_TYPE(float_all,     float,              MPI_FLOAT,             COLL_IO)
GET_VARS_TYPE(double_all,    double,             MPI_DOUBLE,            COLL_IO)
GET_VARS_TYPE(longlong_all,  long long,          MPI_LONG_LONG_INT,     COLL_IO)
GET_VARS_TYPE(ulonglong_all, unsigned long long, MPI_UNSIGNED_LONG_LONG,COLL_IO)
// GET_VARS_TYPE(string_all, char*,              MPI_CHAR,              COLL_IO)
/* string is not yet supported */


/* for write case, buf needs to swapped back if swapped previously */
#define FINAL_CLEAN_UP {                                                       \
    if (is_buf_swapped) /* byte-swap back to buf's original contents */        \
        ncmpii_in_swapn(buf, fnelems, ncmpix_len_nctype(varp->type));          \
                                                                               \
    if (xbuf != NULL && xbuf != cbuf) NCI_Free(xbuf);                          \
    if (cbuf != NULL && cbuf !=  buf) NCI_Free(cbuf);                          \
}

/* The principle of buffer management is:

   for put_vars:
     1. pack buf to cbuf based on buftype
     2. type convert and byte swap cbuf to xbuf
     3. write from xbuf
     4. byte swap the buf back to its original, if it is swapped
     5. free up temp buffers, cbuf and xbuf if they were allocated separately

   for get_vars:
     1. allocate cbuf
     2. allocate xbuf
     3. read to xbuf
     4. type convert and byte swap xbuf to cbuf
     5. unpack cbuf to buf based on buftype
     6. free up temp buffers, cbuf and xbuf if they were allocated separately
*/

/*----< ncmpii_getput_vars() >-----------------------------------------------*/
int
ncmpii_getput_vars(NC               *ncp,
                   NC_var           *varp,
                   const MPI_Offset  start[],
                   const MPI_Offset  count[],
                   const MPI_Offset  stride[],
                   void             *buf,
                   MPI_Offset        bufcount,
                   MPI_Datatype      buftype,  /* data type of the bufer */
                   int               rw_flag,
                   int               io_method)
{
    void *xbuf=NULL, *cbuf=NULL;
    int el_size, buftype_is_contig, mpireturn, need_swap, is_buf_swapped=0;
    int isderived, mpi_err;
    int warning, err, status; /* err is for API abort and status is not */
    MPI_Offset fnelems, bnelems, nbytes, offset=0;
    MPI_Status mpistatus;
    MPI_Datatype ptype, filetype=MPI_BYTE;
    MPI_File fh;

    /* "API error" will abort this API call, but not the entire program */
    err = status = warning = mpi_err = NC_NOERR;

    if (varp->ndims > 0) {
        assert(start != NULL);
        assert(count != NULL);
    }

    if (io_method == COLL_IO)
        fh = ncp->nciop->collective_fh;
    else
        fh = ncp->nciop->independent_fh;

    /* find the ptype (primitive MPI data type) from buftype
     * el_size is the element size of ptype
     * bnelems is the total number of ptype elements in the I/O buffer, buf
     * fnelems is the number of nc variable elements in nc_type
     * nbytes is the amount of read/write in bytes
     */
    err = ncmpii_dtype_decode(buftype, &ptype, &el_size, &bnelems,
                              &isderived, &buftype_is_contig);
    /* bnelems now is the number of ptype in a buftype */
    if (err != NC_NOERR) goto err_check;

    err = NCMPII_ECHAR(varp->type, ptype);
    if (err != NC_NOERR) goto err_check;   

    CHECK_NELEMS(varp, fnelems, count, bnelems, bufcount, nbytes)
    /* bnelems now is the number of ptype in the whole buf */
    /* warning is set in CHECK_NELEMS() */
    need_swap = ncmpii_need_swap(varp->type, ptype);

    if (bnelems == 0) { /* if this process has nothing to read/write */
        if (io_method == INDEP_IO)
            return NCcoordck(ncp, varp, start);
#ifdef ZERO_COUNT_IGNORE_OTHER_ERRORS
        else
        /* for collective I/O, even bnelems == 0, must go on to participate
           the collective calls: MPI_File_set_view and collective read/write */
            goto err_check;
#endif
    }

    /* cbuf is the "contiguous" buffer that all I/O data are packed into a
     * contiguous memory space pointed by cbuf */
    if (!buftype_is_contig) {
        /* pack buf into cbuf, a contiguous buffer */
        cbuf = (void*) NCI_Malloc(bnelems * el_size);
        if (rw_flag == WRITE_REQ) {
            err = ncmpii_data_repack((void*)buf, bufcount, buftype,
                                     cbuf, bnelems, ptype);
            if (err != NC_NOERR) /* API error */
                goto err_check;
        }
    } else {
        cbuf = (void*) buf;
    }

    /* xbuf is the buffer whose data has been converted into the external
     * data type, ready to read/written from/to the netCDF file
     * Now, type convert and byte swap cbuf into xbuf
     */
    if ( ncmpii_need_convert(varp->type, ptype) ) {
        /* allocate new buffer for data type conversion */
        xbuf = NCI_Malloc(nbytes);

        if (rw_flag == WRITE_REQ) {
            /* automatic numeric datatype conversion + swap if necessary
               and only xbuf could be byte-swapped, not cbuf */
            DATATYPE_PUT_CONVERT(varp->type, xbuf, cbuf, bnelems, ptype)
            /* status may be set at DATATYPE_PUT_CONVERT() */
        }
    } else if (need_swap) {
        if (rw_flag == WRITE_REQ) { /* perform array in-place byte swap */
            ncmpii_in_swapn(cbuf, fnelems, ncmpix_len_nctype(varp->type));
            is_buf_swapped = (cbuf == buf) ? 1 : 0;
            /* is_buf_swapped indicates if the contents of the original user
             * buffer, buf, have been changed, i.e. byte swapped. */
        }
        xbuf = cbuf;
    } else {
        /* else, no type conversion or byte swap */
        xbuf = cbuf;
    }

    /* if record variables are too big (so big that we cannot store the
     * stride between records in an MPI_Aint, for example) then we will
     * have to process this one record at a time.  
     */

    /* check if the request is contiguous in file */
    if (stride == NULL && ncmpii_is_request_contiguous(varp, start, count)) {
        err = NCedgeck(ncp, varp, start, count);

        if (err != NC_NOERR ||
            (rw_flag == READ_REQ && IS_RECVAR(varp) &&
             start[0] + count[0] > NC_get_numrecs(ncp))) { /* API error */
            err = NCcoordck(ncp, varp, start);
            if (err != NC_NOERR) {
                goto err_check;
            }
            else {
                err = NC_EEDGE;
                goto err_check;
            }
        }

        /* this is a contiguous file access, no need to set filetype */
        err = ncmpii_get_offset(ncp, varp, start, NULL, NULL, &offset);
        /* if start[] is out of defined size, then this will return
         * NC_EINVALCOORDS error */
        if (err != NC_NOERR) /* API error */
            goto err_check;
    }
    else {
        /* this request is non-contiguous in file, create the mpi file type */
        err = ncmpii_vars_create_filetype(ncp, varp, start, count, stride,
                                          rw_flag, &offset, &filetype);
        if (err != NC_NOERR) /* API error */
            goto err_check;
    }

err_check:
    /* check API error from any proc before going into a collective call.
     * optimization: to avoid MPI_Allreduce to check parameters at
     * every call, we assume caller does the right thing most of the
     * time.  If caller passed in bad parameters, we'll still conduct a
     * zero-byte operation (everyone has to participate in the
     * collective I/O call) but return error */
    if (err != NC_NOERR) {
        /* release allocated resources */
        if (filetype != MPI_BYTE)
            MPI_Type_free(&filetype);
        filetype = MPI_BYTE;

        if (io_method == INDEP_IO) {
            FINAL_CLEAN_UP  /* swap back the data and free buffers */
            return err;
        }
        else /* io_method == COLL_IO */
	    nbytes = offset = 0;
	    /* we want all processors to return from this collective call.
             * Note that we expect the application to check for errors, but if
             * we hang because one process returned early and other processors
             * are in MPI_File_set_view, then that's not good either.  */
    }

    /* MPI_File_set_view is a collective if (io_method == COLL_IO) */
    mpireturn = MPI_File_set_view(fh, offset, MPI_BYTE, filetype,
                                  "native", MPI_INFO_NULL);
    CHECK_MPI_ERROR(mpireturn, "MPI_File_set_view", NC_EFILE)

    if (filetype != MPI_BYTE)
        MPI_Type_free(&filetype);

    if (rw_flag == WRITE_REQ) {
        if (io_method == COLL_IO) {
            mpireturn = MPI_File_write_all(fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_write_all", NC_EWRITE)
        }
        else { /* io_method == INDEP_IO */
            mpireturn = MPI_File_write(fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_write", NC_EWRITE)
        }
    }
    else {  /* rw_flag == READ_REQ */
        if (io_method == COLL_IO) {
            mpireturn = MPI_File_read_all(fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_read_all", NC_EREAD)
        }
        else { /* io_method == INDEP_IO */
            mpireturn = MPI_File_read(fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_read", NC_EREAD)
        }
    }

    /* reset the file view so the entire file is visible again */
    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    if (err != NC_NOERR) { /* no need to go further */
        FINAL_CLEAN_UP  /* swap back the data and free buffers */
        return err;
    }

    if (bnelems == 0) /* no need to go further */
        return ((warning != NC_NOERR) ? warning
                                      : ((status != NC_NOERR) ? status
                                                              : mpi_err));

    /* only bnelems > 0 needs to proceed the following */
    if (rw_flag == READ_REQ) {
        if ( ncmpii_need_convert(varp->type, ptype) ) {
            /* type conversion + swap from xbuf to cbuf*/
            DATATYPE_GET_CONVERT(varp->type, xbuf, cbuf, bnelems, ptype)
        } else if (need_swap) {
            /* perform array in-place byte swap from xbuf to cbuf */
            ncmpii_in_swapn(cbuf, fnelems, ncmpix_len_nctype(varp->type));
        }

        if (!buftype_is_contig) {
            /* unpack cbuf to buf using buftype */
            err = ncmpii_data_repack(cbuf, bnelems, ptype,
                                     (void*)buf, bufcount, buftype);
            if (err != NC_NOERR) {
                FINAL_CLEAN_UP  /* swap back the data and free buffers */
                return err;
            }
        }
    }
    else {
        if (status == NC_NOERR && IS_RECVAR(varp)) {
            /* update header's number of records in memory */
            MPI_Offset newnumrecs;
            if (stride == NULL)
                newnumrecs = start[0] + count[0];
            else
                newnumrecs = start[0] + (count[0] - 1) * stride[0] + 1;

            if (io_method == INDEP_IO) {
                ncp->numrecs = newnumrecs;
                set_NC_ndirty(ncp);
            }
            else
                ncmpii_update_numrecs(ncp, newnumrecs);
        }
    }

    FINAL_CLEAN_UP  /* swap back the data and free buffers */

    return ((warning != NC_NOERR) ? warning
                                  : ((status != NC_NOERR) ? status
                                                          : mpi_err));
}
