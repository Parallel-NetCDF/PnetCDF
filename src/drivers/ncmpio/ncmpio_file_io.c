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
#include <limits.h> /* INT_MAX */

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
                  int           len,
                  MPI_Datatype  buf_type,
                  void         *buf,
                  int           buftype_is_contig)
{
    int status=NC_NOERR, mpireturn, err;
    MPI_Status mpistatus;
    MPI_File fh;
#if MPI_VERSION >= 3
    MPI_Count req_size;
    /* MPI_Type_size_x is introduced in MPI 3.0 */
    MPI_Type_size_x(buf_type, &req_size);
#else
    int req_size;
    MPI_Type_size(buf_type, &req_size);
#endif

    /* request size in bytes */
    req_size *= len;

#ifndef ENABLE_LARGE_SINGLE_REQ
#if MPI_VERSION >= 3
    if (req_size > INT_MAX)
        /* I/O request size > 2 GiB, ROMIO currently does not support a single
         * read/write call of amount > 2 GiB
         */
#else
    if (req_size < 0)
        /* In MPI 2.x and prior, argument "size" in MPI_Type_size is defined
         * as of type int. When int overflow occurs, the returned value in
         * "size" argument may be a negative. This means the aggregated request
         * size > 2 GiB. However, ROMIO currently does not support a single
         * request with amount > 2 GiB
         */
#endif
    {
        if (ncp->safe_mode) {
            if (rw_flag == NC_REQ_RD)
                printf("Error at %s at %d: size of read request (%lld) > INT_MAX\n",
                       __FILE__,__LINE__,(long long)req_size);
            else
                printf("Error at %s at %d: size of write request (%lld) > INT_MAX\n",
                       __FILE__,__LINE__,(long long)req_size);
        }
        if (coll_indep == NC_REQ_INDEP) DEBUG_RETURN_ERROR(NC_EMAX_REQ)
        DEBUG_ASSIGN_ERROR(status, NC_EMAX_REQ)
        buf_type = MPI_BYTE;
        len = 0; /* allow this process to participate collective call */
    }
#endif

#ifdef _USE_MPI_GET_COUNT
    /* explicitly initialize mpistatus object to 0. For zero-length read,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around.
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));
#endif

    if (coll_indep == NC_REQ_COLL)
        fh = ncp->collective_fh;
    else
        fh = ncp->independent_fh;

    if (rw_flag == NC_REQ_RD) {
        void         *xbuf=buf;
        int           xlen=len;
        MPI_Datatype  xbuf_type=buf_type;

        /* if the read buffer is noncontiguous and size is < ncp->ibuf_size,
         * allocate a temporary buffer and use it to read, as some MPI, e.g.
         * Cray on KNL, can be significantly slow when read buffer is
         * noncontiguous.
         */
        if (len > 0 && !buftype_is_contig && req_size <= ncp->ibuf_size) {
            xlen = req_size;
            xbuf = NCI_Malloc(xlen);
            xbuf_type = MPI_BYTE;
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
                    DEBUG_ASSIGN_ERROR(status, err)
                }
            }
        }
        if (mpireturn == MPI_SUCCESS) {
            /* update the number of bytes read since file open */
#ifdef _USE_MPI_GET_COUNT
            int get_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
            ncp->get_size += get_size;
#else
            ncp->get_size += req_size;
#endif
        }
        if (xbuf != buf) { /* unpack contiguous xbuf to noncontiguous buf */
            int pos=0;
            MPI_Unpack(xbuf, xlen, &pos, buf, len, buf_type, MPI_COMM_SELF);
            NCI_Free(xbuf);
        }
    } else { /* NC_REQ_WR */
        void         *xbuf=buf;
        int           xlen=len;
        MPI_Datatype  xbuf_type=buf_type;

        /* if the write buffer is noncontiguous and size is < ncp->ibuf_size,
         * allocate a temporary buffer and use it to write, as some MPI, e.g.
         * Cray on KNL, can be significantly slow when write buffer is
         * noncontiguous.
         */
        if (len > 0 && !buftype_is_contig && req_size <= ncp->ibuf_size) {
            int pos=0;
            xlen = req_size;
            xbuf = NCI_Malloc(xlen);
            MPI_Pack(buf, len, buf_type, xbuf, xlen, &pos, MPI_COMM_SELF);
            xbuf_type = MPI_BYTE;
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
                    DEBUG_ASSIGN_ERROR(status, err)
                }
            }
        }
        if (mpireturn == MPI_SUCCESS) {
            /* update the number of bytes written since file open */
#ifdef _USE_MPI_GET_COUNT
            int put_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
            ncp->put_size += put_size;
#else
            ncp->put_size += req_size;
#endif
        }
        if (xbuf != buf) NCI_Free(xbuf);
    }

    return status;
}

