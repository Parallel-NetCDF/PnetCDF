/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset() */

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< ncmpio_read_write() >------------------------------------------------*/
int
ncmpio_read_write(NC           *ncp,
                  int           rw_flag,     /* NC_REQ_WR or NC_REQ_RD */
                  int           coll_indep,  /* NC_REQ_COLL or NC_REQ_INDEP */
                  MPI_Offset    offset,
                  MPI_Offset    buf_count,
                  MPI_Datatype  buf_type,
                  void         *buf,
                  int           buftype_is_contig)
{
    int status=NC_NOERR, err=NC_NOERR, mpireturn;
    MPI_Status mpistatus;
    MPI_File fh;
    MPI_Offset req_size;

#ifdef HAVE_MPI_TYPE_SIZE_C
    MPI_Count btype_size;
    /* MPI_Type_size_c is introduced in MPI 4.0 */
    mpireturn = MPI_Type_size_c(buf_type, &btype_size);
#elif defined(HAVE_MPI_TYPE_SIZE_X)
    MPI_Count btype_size;
    /* MPI_Type_size_x is introduced in MPI 3.0 */
    mpireturn = MPI_Type_size_x(buf_type, &btype_size);
#else
    int btype_size;
    mpireturn = MPI_Type_size(buf_type, &btype_size);
#endif
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size");
        /* return the first encountered error if there is any */
        err = (err == NC_EFILE) ? NC_EREAD : err;
    }
    else if (btype_size == MPI_UNDEFINED)
        DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)

    if (err != NC_NOERR) {
        if (coll_indep == NC_REQ_COLL) {
            DEBUG_ASSIGN_ERROR(status, err)
            /* write nothing, but participate the collective call */
            buf_count = 0;
        }
        else
            DEBUG_RETURN_ERROR(err)
    }

    /* request size in bytes, may be > NC_MAX_INT */
    req_size = buf_count * btype_size;

    /* explicitly initialize mpistatus object to 0. For zero-length read,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around.
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (coll_indep == NC_REQ_COLL)
        fh = ncp->collective_fh;
    else
        fh = ncp->independent_fh;

    if (rw_flag == NC_REQ_RD) {
        void         *xbuf=buf;
        int           xlen=(int)buf_count;
        MPI_Datatype  xbuf_type=buf_type;

        if (buf_count > NC_MAX_INT) {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_contiguous_c((MPI_Count)buf_count, buf_type, &xbuf_type);
            MPI_Type_commit(&xbuf_type);
            xlen = 1;
#else
            if (coll_indep == NC_REQ_COLL)
                DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
            else
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
#endif
        }
        else if (buf_count > 0 && !buftype_is_contig &&
                 req_size <= ncp->ibuf_size) {
            /* if read buffer is noncontiguous and size is < ncp->ibuf_size,
             * allocate a temporary buffer and use it to read, as some MPI,
             * e.g. Cray on KNL, can be significantly slow when read buffer is
             * noncontiguous.
             */
            if (req_size > NC_MAX_INT) {
                MPI_Type_contiguous((int)buf_count, buf_type, &xbuf_type);
                MPI_Type_commit(&xbuf_type);
                xlen = 1;
            }
            else {
                xbuf_type = MPI_BYTE;
                xlen = (int)req_size;
            }
            xbuf = NCI_Malloc((size_t)req_size);
        }

        if (coll_indep == NC_REQ_COLL) {
            TRACE_IO(MPI_File_read_at_all)(fh, offset, xbuf, xlen, xbuf_type,
                                           &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_read_at_all");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) {
                    err = (err == NC_EFILE) ? NC_EREAD : err;
                    DEBUG_ASSIGN_ERROR(status, err)
                }
            }
        } else {
            TRACE_IO(MPI_File_read_at)(fh, offset, xbuf, xlen, xbuf_type,
                                       &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_read_at");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) {
                    err = (err == NC_EFILE) ? NC_EREAD : err;
                    DEBUG_RETURN_ERROR(err)
                }
            }
        }
        if (mpireturn == MPI_SUCCESS) {
            /* update the number of bytes read since file open */
#ifdef HAVE_MPI_GET_COUNT_C
            MPI_Count get_size;
            MPI_Get_count_c(&mpistatus, MPI_BYTE, &get_size);
            ncp->get_size += get_size;
#else
            int get_size;
            mpireturn = MPI_Get_count(&mpistatus, xbuf_type, &get_size);
            if (mpireturn != MPI_SUCCESS || get_size == MPI_UNDEFINED)
                ncp->get_size += req_size;
            else {
#ifdef HAVE_MPI_TYPE_SIZE_X
                /* MPI_Type_size_x is introduced in MPI 3.0 */
                mpireturn = MPI_Type_size_x(xbuf_type, &btype_size);
#else
                mpireturn = MPI_Type_size(xbuf_type, &btype_size);
#endif
                if (mpireturn != MPI_SUCCESS || get_size == MPI_UNDEFINED)
                    ncp->get_size += req_size;
                else
                    ncp->get_size += btype_size * get_size;
            }
#endif
        }
        if (xbuf != buf) { /* unpack contiguous xbuf to noncontiguous buf */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count pos=0;
            MPI_Unpack_c(xbuf, xlen, &pos, buf, (MPI_Count)buf_count, buf_type,
                         MPI_COMM_SELF);
#else
            int pos=0;
            MPI_Unpack(xbuf, xlen, &pos, buf, (int)buf_count, buf_type,
                       MPI_COMM_SELF);
#endif
            NCI_Free(xbuf);
        }
        if (xbuf_type != buf_type && xbuf_type != MPI_BYTE)
            MPI_Type_free(&xbuf_type);
    } else { /* NC_REQ_WR */
        void         *xbuf=buf;
        int           xlen=(int)buf_count;
        MPI_Datatype  xbuf_type=buf_type;

        if (buf_count > NC_MAX_INT) {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_contiguous_c((MPI_Count)buf_count, buf_type, &xbuf_type);
            MPI_Type_commit(&xbuf_type);
            xlen = 1;
#else
            if (coll_indep == NC_REQ_COLL)
                DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
            else
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
#endif
        }
        else if (buf_count > 0 && !buftype_is_contig &&
                 req_size <= ncp->ibuf_size) {
            /* if write buffer is noncontiguous and size is < ncp->ibuf_size,
             * allocate a temporary buffer and use it to write, as some MPI,
             * e.g. Cray on KNL, can be significantly slow when write buffer is
             * noncontiguous.
             */
            if (req_size > NC_MAX_INT) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count pos=0;
                xbuf = NCI_Malloc(req_size);
                MPI_Pack_c(buf, (MPI_Count)buf_count, buf_type, xbuf,
                           (MPI_Count)req_size, &pos, MPI_COMM_SELF);
                MPI_Type_contiguous_c((MPI_Count)req_size, MPI_BYTE, &xbuf_type);
                MPI_Type_commit(&xbuf_type);
                xlen = 1;
#else
                /* skip packing write data into a temp buffer */
                xlen = (int)buf_count;
                xbuf_type = buf_type;
#endif
            }
            else {
                int pos=0;
                xlen = (int)req_size;
                xbuf = NCI_Malloc(xlen);
                MPI_Pack(buf, (int)buf_count, buf_type, xbuf, xlen, &pos,
                         MPI_COMM_SELF);
                xbuf_type = MPI_BYTE;
            }
        }

        if (coll_indep == NC_REQ_COLL) {
            TRACE_IO(MPI_File_write_at_all)(fh, offset, xbuf, xlen, xbuf_type,
                                            &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_write_at_all");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) {
                    err = (err == NC_EFILE) ? NC_EWRITE : err;
                    DEBUG_ASSIGN_ERROR(status, err)
                }
            }
        } else {
            TRACE_IO(MPI_File_write_at)(fh, offset, xbuf, xlen, xbuf_type,
                                        &mpistatus);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_write_at");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) {
                    err = (err == NC_EFILE) ? NC_EWRITE : err;
                    DEBUG_RETURN_ERROR(err)
                }
            }
        }
        if (mpireturn == MPI_SUCCESS) {
            /* update the number of bytes written since file open */
#ifdef HAVE_MPI_GET_COUNT_C
            MPI_Count put_size;
            MPI_Get_count_c(&mpistatus, MPI_BYTE, &put_size);
            ncp->put_size += put_size;
#else
            int put_size;
            mpireturn = MPI_Get_count(&mpistatus, xbuf_type, &put_size);
            if (mpireturn != MPI_SUCCESS || put_size == MPI_UNDEFINED)
                ncp->put_size += req_size;
            else {
#ifdef HAVE_MPI_TYPE_SIZE_X
                /* MPI_Type_size_x is introduced in MPI 3.0 */
                mpireturn = MPI_Type_size_x(xbuf_type, &btype_size);
#else
                mpireturn = MPI_Type_size(xbuf_type, &btype_size);
#endif
                if (mpireturn != MPI_SUCCESS || put_size == MPI_UNDEFINED)
                    ncp->put_size += req_size;
                else
                    ncp->put_size += btype_size * put_size;
            }
#endif
        }
        if (xbuf != buf) NCI_Free(xbuf);
        if (xbuf_type != buf_type && xbuf_type != MPI_BYTE)
            MPI_Type_free(&xbuf_type);
    }

    return status;
}

