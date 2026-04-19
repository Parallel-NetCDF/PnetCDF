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
#include <fcntl.h>  /* O_RDONLY, O_RDWR */

#include <assert.h>

#include <mpi.h>

#ifdef ENABLE_GIO
#include <gio.h>
#endif

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< get_count() >--------------------------------------------------------*/
/* This subroutine is independent. On success, the number of bytes read/written
 * is returned (zero indicates nothing was read/written). Like POSIX read()/
 * write(), it is not an error if this number is smaller than the number of
 * bytes requested. On error, a negative value, an NC error code, is returned.
 *
 * This subroutine is called only when using MPI-IO driver.
 */
static
MPI_Offset get_count(MPI_Status   *mpistatus,
                     MPI_Datatype  datatype)
{
    int mpireturn;

    if (datatype == MPI_DATATYPE_NULL) return 0;

#ifdef HAVE_MPI_TYPE_SIZE_C
    MPI_Count type_size;
    /* MPI_Type_size_c is introduced in MPI 4.0 */
    MPI_Type_size_c(datatype, &type_size);
#elif defined(HAVE_MPI_TYPE_SIZE_X)
    MPI_Count type_size;
    /* MPI_Type_size_x is introduced in MPI 3.0 */
    MPI_Type_size_x(datatype, &type_size);
#else
    int type_size;
    MPI_Type_size(datatype, &type_size);
#endif

#ifdef HAVE_MPI_GET_COUNT_C
    MPI_Count count;
    mpireturn = MPI_Get_count_c(mpistatus, datatype, &count);
#else
    int count;
    mpireturn = MPI_Get_count(mpistatus, datatype, &count);
#endif

    if (mpireturn != MPI_SUCCESS || count == MPI_UNDEFINED)
        /* In case of partial read/write, MPI_Get_elements() is supposed to be
         * called to obtain the number of type map elements actually
         * read/written in order to calculate the true read/write amount. Below
         * skips this step and simply returns the partial read/write amount.
         * See an example usage of MPI_Get_count() in Example 5.12 from MPI
         * standard document.
         */
        return NC_EFILE;

    return (MPI_Offset)count * type_size;
}

/*----< ncmpio_file_read_at() >----------------------------------------------*/
/* This function is independent.
 * TODO: move check count against MAX_INT and call _c API.
 *
 * This subroutine is used by vard APIs only.
 */
MPI_Offset
ncmpio_file_read_at(NC         *ncp,
                    MPI_Offset  offset,
                    void       *buf,
                    PNCIO_View  buf_view)
{
    if (ncp->driver != PNC_DRIVER_MPIIO) assert(0);

    int err=NC_NOERR;
    MPI_Offset amnt=0;

    /* For zero-sized requests, independent read can return now. */
    if (buf_view.size == 0) return 0;

    if (ncp->driver == PNC_DRIVER_MPIIO) {
        char *mpi_name;
        int mpireturn;
        MPI_File fh;
        MPI_Status mpistatus;

        /* Explicitly initialize mpistatus object to 0. For zero-length
         * read/write, MPI_Get_count() may report incorrect result for some
         * earlier MPICH versions, due to the uninitialized MPI_Status object
         * passed to MPI-IO calls. Thus we initialize it above to work around.
         * See MPICH ticket: https://trac.mpich.org/projects/mpich/ticket/2332
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

        if (fh == MPI_FILE_NULL) return 0;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count count = (buf_view.type == MPI_BYTE) ? buf_view.size : 1;

        TRACE_IO(MPI_File_read_at_c, (fh, offset, buf, count, buf_view.type,
                                      &mpistatus));
#else
        int count = (buf_view.type == MPI_BYTE) ? buf_view.size : 1;

        if (buf_view.size > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buffer size="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,buf_view.size);
#endif
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
        TRACE_IO(MPI_File_read_at, (fh, offset, buf, count, buf_view.type,
                                    &mpistatus));
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD)
        }

        /* update the number of bytes read since file open */
        if (err == NC_NOERR)
            amnt = get_count(&mpistatus, buf_view.type);
    }
#ifdef ENABLE_GIO
    else if (ncp->driver == PNC_DRIVER_GIO) {
        amnt = GIO_read_at(ncp->gio_fh, buf, ncp->file_view.count,
                                             ncp->file_view.off,
                                             ncp->file_view.len,
                                             buf_view.count,
                                             buf_view.off,
                                             buf_view.len);
    }
#endif
    else
        amnt = NC_EDRIVER;

    /* update the number of bytes read since file open */
    if (amnt >= 0) ncp->get_size += amnt;
    /* else: ignore if error, as this error is not fatal */

    return amnt;
}

/*----< ncmpio_file_read_at_all() >------------------------------------------*/
/* This function is collective.
 * This subroutine is used by vard APIs only.
 */
MPI_Offset
ncmpio_file_read_at_all(NC         *ncp,
                        MPI_Offset  offset,
                        void       *buf,
                        PNCIO_View  buf_view)
{
    if (ncp->driver != PNC_DRIVER_MPIIO) assert(0);

    int err=NC_NOERR;
    MPI_Offset amnt=0;

    if (ncp->driver == PNC_DRIVER_MPIIO) {
        char *mpi_name;
        int mpireturn;
        MPI_File fh;
        MPI_Status mpistatus;

        /* Explicitly initialize mpistatus object to 0. For zero-length
         * read/write, MPI_Get_count() may report incorrect result for some
         * earlier MPICH versions, due to the uninitialized MPI_Status object
         * passed to MPI-IO calls. Thus we initialize it above to work around.
         * See MPICH ticket: https://trac.mpich.org/projects/mpich/ticket/2332
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

        if (fh == MPI_FILE_NULL) return 0;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count count = (buf_view.type == MPI_BYTE) ? buf_view.size : 1;

        TRACE_IO(MPI_File_read_at_all_c, (fh, offset, buf, count,
                                          buf_view.type, &mpistatus));
#else
        int count = (buf_view.type == MPI_BYTE) ? buf_view.size : 1;

        if (buf_view.size > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buffer size="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,buf_view.size);
#endif
            DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
            /* participate the collective call, but read nothing */
            count = 0;
        }
        TRACE_IO(MPI_File_read_at_all, (fh, offset, buf, count,
                                        buf_view.type, &mpistatus));
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD)
        }

        /* update the number of bytes read since file open */
        if (err == NC_NOERR)
            amnt = get_count(&mpistatus, buf_view.type);
    }
#ifdef ENABLE_GIO
    else if (ncp->driver == PNC_DRIVER_GIO) {
        if (ncp->num_aggrs_per_node == 0 ||  /* INA is disabled */
            ncp->comm_attr.is_ina_aggr)      /* is an INA aggregator */
        /* When INA is disabled, all processes must participate this collective
         * read. When INA is enabled, only the INA aggregators participate.
         */
        amnt = GIO_read_at(ncp->gio_fh, buf, ncp->file_view.count,
                                             ncp->file_view.off,
                                             ncp->file_view.len,
                                             buf_view.count,
                                             buf_view.off,
                                             buf_view.len);
    }
#endif
    else
        amnt = NC_EDRIVER;

    /* update the number of bytes read since file open */
    if (amnt >= 0) ncp->get_size += amnt;
    /* else: ignore if error, as this error is not fatal */

    return amnt;
}

/*----< ncmpio_file_write_at() >---------------------------------------------*/
/* This function is independent.
 * This subroutine is used by vard APIs only.
 */
MPI_Offset
ncmpio_file_write_at(NC         *ncp,
                     MPI_Offset  offset,
                     const void *buf,
                     PNCIO_View  buf_view)
{
    if (ncp->driver != PNC_DRIVER_MPIIO) assert(0);

    int err=NC_NOERR;
    MPI_Offset amnt=0;

    /* For zero-sized requests, independent read can return now. */
    if (buf_view.size == 0) return 0;

    if (ncp->driver == PNC_DRIVER_MPIIO) {
        char *mpi_name;
        int mpireturn;
        MPI_File fh;
        MPI_Status mpistatus;

        /* Explicitly initialize mpistatus object to 0. For zero-length
         * read/write, MPI_Get_count() may report incorrect result for some
         * earlier MPICH versions, due to the uninitialized MPI_Status object
         * passed to MPI-IO calls. Thus we initialize it above to work around.
         * See MPICH ticket: https://trac.mpich.org/projects/mpich/ticket/2332
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

        if (fh == MPI_FILE_NULL) return 0;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count count = (buf_view.type == MPI_BYTE) ? buf_view.size : 1;

        TRACE_IO(MPI_File_write_at_c, (fh, offset, buf, count, buf_view.type,
                                       &mpistatus));
#else
        int count = (buf_view.type == MPI_BYTE) ? buf_view.size : 1;

        if (buf_view.size > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buffer size="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,buf_view.size);
#endif
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
        TRACE_IO(MPI_File_write_at, (fh, offset, buf, count, buf_view.type,
                                     &mpistatus));
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE)
        }

        if (err == NC_NOERR)
            amnt = get_count(&mpistatus, buf_view.type);
    }
#ifdef ENABLE_GIO
    else if (ncp->driver == PNC_DRIVER_GIO) {
        amnt = GIO_write_at(ncp->gio_fh, buf, ncp->file_view.count,
                                              ncp->file_view.off,
                                              ncp->file_view.len,
                                              buf_view.count,
                                              buf_view.off,
                                              buf_view.len);
    }
#endif
    else
        amnt = NC_EDRIVER;

    /* update the number of bytes written since file open */
    if (amnt >= 0) ncp->put_size += amnt;
    /* else: ignore if error, as this error is not fatal */

    return amnt;
}

/*----< ncmpio_file_write_at_all() >-----------------------------------------*/
/* This function is collective.
 * This subroutine is used by vard APIs only.
 */
MPI_Offset
ncmpio_file_write_at_all(NC         *ncp,
                         MPI_Offset  offset,
                         const void *buf,
                         PNCIO_View  buf_view)
{
    if (ncp->driver != PNC_DRIVER_MPIIO) assert(0);

    int err=NC_NOERR;
    MPI_Offset amnt=0;

    if (ncp->driver == PNC_DRIVER_MPIIO) {
        char *mpi_name;
        int mpireturn;
        MPI_File fh;
        MPI_Status mpistatus;

        /* Explicitly initialize mpistatus object to 0. For zero-length
         * read/write, MPI_Get_count() may report incorrect result for some
         * earlier MPICH versions, due to the uninitialized MPI_Status object
         * passed to MPI-IO calls. Thus we initialize it above to work around.
         * See MPICH ticket: https://trac.mpich.org/projects/mpich/ticket/2332
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

        if (fh == MPI_FILE_NULL) return 0;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count count = (buf_view.type == MPI_BYTE) ? buf_view.size : 1;

        TRACE_IO(MPI_File_write_at_all_c, (fh, offset, buf, count,
                                           buf_view.type, &mpistatus));
#else
        int count = (buf_view.type == MPI_BYTE) ? buf_view.size : 1;

        if (buf_view.size > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buffer size="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,buf_view.size);
#endif
            DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
            /* participate the collective call, but write nothing */
            count = 0;
        }
        TRACE_IO(MPI_File_write_at_all, (fh, offset, buf, count,
                                         buf_view.type, &mpistatus));
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE)
        }

        if (err == NC_NOERR)
            amnt = get_count(&mpistatus, buf_view.type);
    }
#ifdef ENABLE_GIO
    else if (ncp->driver == PNC_DRIVER_GIO) {
        if (ncp->num_aggrs_per_node == 0 ||  /* INA is disabled */
            ncp->comm_attr.is_ina_aggr)      /* is an INA aggregator */
        /* When INA is disabled, all processes must participate this collective
         * write. When INA is enabled, only the INA aggregators participate.
         */
        amnt = GIO_write_at_all(ncp->gio_fh, buf, ncp->file_view.count,
                                                  ncp->file_view.off,
                                                  ncp->file_view.len,
                                                  buf_view.count,
                                                  buf_view.off,
                                                  buf_view.len);
    }
#endif
    else
        amnt = NC_EDRIVER;

    /* update the number of bytes written since file open */
    if (amnt >= 0) ncp->put_size += amnt;
    /* else: ignore if error, as this error is not fatal */

    return amnt;
}

/*----< ncmpio_getput_zero_req() >-------------------------------------------*/
/* This function is called when this process has zero-length I/O request and
 * must participate all the MPI collective calls involved in the collective
 * APIs and wait_all(), which include:
 * 1. setting file view,
 * 2. collective read/write,
 * 3. re-setting file view to make the entire file visible.
 *
 * This function is collective.
 *
 * This subroutine is used by vard APIs only.
 */
int
ncmpio_getput_zero_req(NC *ncp, int reqMode)
{
    if (ncp->driver != PNC_DRIVER_MPIIO) assert(0);

    int err, status=NC_NOERR;
    MPI_Offset rlen, wlen;
    PNCIO_View buf_view;

    buf_view.size = 0;

    /* When intra-node aggregation is enabled, non-INA aggregators do not
     * access the file and participate any collective read or write. Thus
     * non-aggregators can return now.
     */
    if (ncp->num_aggrs_per_node > 0 && !ncp->comm_attr.is_ina_aggr)
        return NC_NOERR;

    /* do nothing if this came from an independent API */
    if (fIsSet(reqMode, NC_REQ_INDEP)) return NC_NOERR;

    err = ncmpio_file_set_view(ncp, MPI_BYTE, 0, NULL, NULL);
    if (status == NC_NOERR) status = err;

    if (fIsSet(reqMode, NC_REQ_RD)) {
        if (ncp->nprocs > 1)
            rlen = ncmpio_file_read_at_all(ncp, 0, NULL, buf_view);
        else
            rlen = ncmpio_file_read_at(ncp, 0, NULL, buf_view);
        if (status == NC_NOERR && rlen < 0) status = (int)rlen;
    }
    else { /* write request */
        if (ncp->nprocs > 1)
            wlen = ncmpio_file_write_at_all(ncp, 0, NULL, buf_view);
        else
            wlen = ncmpio_file_write_at(ncp, 0, NULL, buf_view);
        if (status == NC_NOERR && wlen < 0) status = (int)wlen;
    }

    if (ncp->driver == PNC_DRIVER_MPIIO)
        /* Reset file view is necessary only when using MPI-IO driver. Note
         * file view is never reused in PnetCDF.
         */
        ncmpio_file_set_view(ncp, MPI_BYTE, 0, NULL, NULL);

    return status;
}

#if 0
/* This subroutine is no longer used, as it has been split into
 * ncmpio_file_read() and ncmpio_file_write().
 */
/*----< ncmpio_read_write() >------------------------------------------------*/
int
ncmpio_read_write(NC         *ncp,
                  int         rw_flag,     /* NC_REQ_WR or NC_REQ_RD */
                  MPI_Offset  offset,
                  PNCIO_View  buf_view,
                  void       *buf)
{
    char *mpi_name;
    int i, status=NC_NOERR, err=NC_NOERR, mpireturn, coll_indep;
    int to_free_buftype=0;
    MPI_Offset rlen, wlen;

    coll_indep = NC_REQ_INDEP;
    if (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
        coll_indep = NC_REQ_COLL;

    /* for zero-sized request */
    if (buf_view.size == 0) {
        if (coll_indep == NC_REQ_INDEP)
            return NC_NOERR;

        if (rw_flag == NC_REQ_RD) {
            rlen = ncmpio_file_read_at_all(ncp, 0, NULL, buf_view);
            if (rlen < 0) status = (int)rlen;
        }
        else {
            wlen = ncmpio_file_write_at_all(ncp, 0, NULL, buf_view);
            if (wlen < 0) status = (int)wlen;
        }
        goto fn_exit;
    }

    /* buf_view.count is the number of offset-length pairs */

    /* buf_view.size is in bytes, may be > NC_MAX_INT */

    if (rw_flag == NC_REQ_RD) {
        void *xbuf=buf;
        MPI_Count orig_count;

        /* Save the original buf_view.count, as it may be modified when
         * buf_view.size <= ibuf_size, original buf_view.count is required to
         * unpack the read data to user read buffer.
         */
        orig_count = buf_view.count;

#ifndef HAVE_MPI_LARGE_COUNT
        if (buf_view.size > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buffer size="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,buf_view.size);
#endif
            if (coll_indep == NC_REQ_COLL) {
                DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
                /* write nothing, but participate the collective call */
                buf_view.size = 0;
                buf_view.count = 0;
            }
            else
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
#endif

        if (buf_view.count > 1 && buf_view.size <= ncp->ibuf_size) {
            /* The only case of read buffer being noncontiguous is when
             * nonblocking API ncmpi_wait/wait_all() is called and INA is
             * disabled. If read buffer is noncontiguous and size is <
             * ncp->ibuf_size, we allocate a temporary contiguous buffer and
             * use it to read. Later it is unpacked to user buffer. As some
             * MPI, e.g. Cray on KNL, can be significantly slow when write
             * buffer is noncontiguous.
             *
             * Note ncp->ibuf_size is never > NC_MAX_INT.
             */
            xbuf = NCI_Malloc(buf_view.size);
            /* mark buf_view is contiguous */
            buf_view.type = MPI_BYTE;
            buf_view.count = 0;
        }
        else if (buf_view.count > 1 && ncp->driver == PNC_DRIVER_MPIIO) {
            /* construct a buftype */
#ifdef HAVE_MPI_LARGE_COUNT
            /* TODO: MPI_Type_create_hindexed_c
             *       buf_view.count should be of type MPI_Count
             *       buf_view.len   should be of type MPI_Count
             *       buf_view.off   should be of type MPI_Count
             */
            mpireturn = MPI_Type_create_hindexed_c(buf_view.count,
                                                   buf_view.len,
                                                   buf_view.off,
                                                   MPI_BYTE, &buf_view.type);
            mpi_name = "MPI_Type_create_hindexed_c";
#else
            MPI_Aint *disp;
#if SIZEOF_MPI_AINT == SIZEOF_MPI_OFFSET
            disp = (MPI_Aint*) buf_view.off;
#else
            disp = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * buf_view.count);
            for (j=0; j<buf_view.count; j++)
                disp[j] = (MPI_Aint)buf_view.off[j];
#endif
            mpireturn = MPI_Type_create_hindexed(buf_view.count,
                                                 buf_view.len,
                                                 disp,
                                                 MPI_BYTE, &buf_view.type);
            mpi_name = "MPI_Type_create_hindexed";
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
            NCI_Free(disp);
#endif
#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
                buf_view.size = 0;
                buf_view.count = 0;
            }
            else {
                mpireturn = MPI_Type_commit(&buf_view.type);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                    buf_view.size = 0;
                    buf_view.count = 0;
                }
                else
                    to_free_buftype = 1;
            }
        }

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
            rlen = ncmpio_file_read_at_all(ncp, offset, xbuf, buf_view);
        else
            rlen = ncmpio_file_read_at(ncp, offset, xbuf, buf_view);
        if (status == NC_NOERR && rlen < 0) status = (int)rlen;

        if (xbuf != buf) { /* unpack contiguous xbuf to noncontiguous buf */
            char *in_ptr, *out_ptr;
            in_ptr = xbuf;

            for (i=0; i<orig_count; i++) {
                out_ptr = (char*)buf + buf_view.off[i]; // - buf_view.off[0]);
                memcpy(out_ptr, in_ptr, buf_view.len[i]);
                in_ptr += buf_view.len[i];
            }
            NCI_Free(xbuf);
        }
        if (to_free_buftype)
            MPI_Type_free(&buf_view.type);

    } else { /* NC_REQ_WR */
        void *xbuf=buf;

        if (buf_view.count > 1 && buf_view.size <= ncp->ibuf_size) {
            /* The only case of write buffer being noncontiguous is when
             * nonblocking API ncmpi_wait/wait_all() is called and INA is
             * disabled. If write buffer is noncontiguous and size is <
             * ncp->ibuf_size, pack it a temporary contiguous buffer and use it
             * to write. As some MPI, e.g. Cray on KNL, can be significantly
             * slow when write buffer is noncontiguous.
             *
             * Note ncp->ibuf_size is never > NC_MAX_INT.
             */
            char *in_ptr, *out_ptr;
            xbuf = NCI_Malloc(buf_view.size);
            out_ptr = xbuf;

            for (i=0; i<buf_view.count; i++) {
                in_ptr = (char*)buf + buf_view.off[i];
                memcpy(out_ptr, in_ptr, buf_view.len[i]);
                out_ptr += buf_view.len[i];
            }
            /* mark buf_view is contiguous */
            buf_view.type = MPI_BYTE;
            buf_view.count = 0;
        }
        else if (buf_view.count > 1 && ncp->driver == PNC_DRIVER_MPIIO) {
            /* construct a buftype */
#ifdef HAVE_MPI_LARGE_COUNT
            /* TODO: MPI_Type_create_hindexed_c
             *       buf_view.count should be of type MPI_Count
             *       buf_view.len   should be of type MPI_Count
             *       buf_view.off   should be of type MPI_Count
             */
            mpireturn = MPI_Type_create_hindexed_c(buf_view.count,
                                                   buf_view.len,
                                                   buf_view.off,
                                                   MPI_BYTE, &buf_view.type);
            mpi_name = "MPI_Type_create_hindexed_c";
#else
            MPI_Aint *disp;
#if SIZEOF_MPI_AINT == SIZEOF_MPI_OFFSET
            disp = (MPI_Aint*) buf_view.off;
#else
            disp = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * buf_view.count);
            for (j=0; j<buf_view.count; j++)
                disp[j] = (MPI_Aint)buf_view.off[j];
#endif
            /* TODO: MPI_Type_create_hindexed
             *       buf_view.count should be of type int
             *       buf_view.len   should be of type int
             *       buf_view.off   should be of type MPI_Aint
             */
            mpireturn = MPI_Type_create_hindexed(buf_view.count,
                                                 buf_view.len,
                                                 disp,
                                                 MPI_BYTE, &buf_view.type);
            mpi_name = "MPI_Type_create_hindexed";
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
            NCI_Free(disp);
#endif
#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
                buf_view.size = 0;
            }
            else {
                mpireturn = MPI_Type_commit(&buf_view.type);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                    buf_view.size = 0;
                }
                else
                    to_free_buftype = 1;
            }
        }

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
            wlen = ncmpio_file_write_at_all(ncp, offset, xbuf, buf_view);
        else
            wlen = ncmpio_file_write_at(ncp, offset, xbuf, buf_view);
        if (status == NC_NOERR && wlen < 0) status = (int)wlen;

        if (xbuf != buf) NCI_Free(xbuf);
        if (to_free_buftype)
            MPI_Type_free(&buf_view.type);
    }

fn_exit:
    if (ncp->driver == PNC_DRIVER_MPIIO)
        /* Reset file view is necessary only when using MPI-IO driver. Note
         * file view is never reused in PnetCDF.
         */
        ncmpio_file_set_view(ncp, MPI_BYTE, 0, NULL, NULL);

    return status;
}
#endif

/*----< ncmpio_file_close() >------------------------------------------------*/
/*
 * This function is collective.
 */
int
ncmpio_file_close(NC *ncp)
{
    int err=NC_NOERR;

    if (ncp->driver == PNC_DRIVER_MPIIO) {
        char *mpi_name;
        int mpireturn;

        if (ncp->independent_fh != ncp->collective_fh &&
            ncp->independent_fh != MPI_FILE_NULL) {
            TRACE_IO(MPI_File_close, (&ncp->independent_fh));
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }

        if (ncp->collective_fh != MPI_FILE_NULL) {
            TRACE_IO(MPI_File_close, (&ncp->collective_fh));
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }
    }
#ifdef ENABLE_GIO
    else if (ncp->driver == PNC_DRIVER_GIO) {
        if (ncp->gio_fh != NULL)
            /* When INA is enabled, non-INA aggregators' gio_fh may be NULL */
            err = GIO_close(ncp->gio_fh);
        ncp->gio_fh = NULL;
    }
#endif
    else
        err = NC_EDRIVER;

    return err;
}

/*----< ncmpio_file_delete() >-----------------------------------------------*/
/*
 * This function is collective.
 *
 * This subroutine is called only from ncmpi_abort. When the file is being
 * created and an error occurs, the program is still in define mode. In this
 * case, the file is deleted.
 */
int
ncmpio_file_delete(NC *ncp)
{
    int err=NC_NOERR;

    if (ncp->rank == 0) {
        char *path = ncmpii_remove_file_system_type_prefix(ncp->path);
        if (ncp->driver == PNC_DRIVER_MPIIO) {
            char *mpi_name;
            int mpireturn;
#ifdef MPICH_VERSION
            /* MPICH recognizes file system type acronym prefixed to the file name */
            TRACE_IO(MPI_File_delete, ((char *)ncp->path, ncp->mpiinfo));
#else
            TRACE_IO(MPI_File_delete, (path, ncp->mpiinfo));
#endif
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }
#ifdef ENABLE_GIO
        else if (ncp->driver == PNC_DRIVER_GIO) {
            err = GIO_delete(path);
            if (err != GIO_NOERR)
                err = ncmpii_error_gio2nc(err, "GIO_delete");
        }
#endif
        else
            err = NC_EDRIVER;
    }

    if (ncp->nprocs > 1)
        MPI_Bcast(&err, 1, MPI_INT, 0, ncp->comm);

    return err;
}

/*----< ncmpio_file_sync() >-------------------------------------------------*/
/* This function must be called collectively, no matter if it is in collective
 * or independent data mode.
 */
int
ncmpio_file_sync(NC *ncp) {
    char *mpi_name;
    int mpireturn;

#ifdef ENABLE_GIO
    if (ncp->driver == PNC_DRIVER_GIO) {
        int err=NC_NOERR;
        if (ncp->gio_fh != NULL) {
            err = GIO_sync(ncp->gio_fh);
            if (err != GIO_NOERR)
                err = ncmpii_error_gio2nc(err, "GIO_sync");
        }
        return err;
    }
#endif
    if (ncp->driver != PNC_DRIVER_MPIIO)
        return NC_EDRIVER;

    /* the remaining of this subroutine are for when using MPI-IO */

    if (ncp->independent_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync, (ncp->independent_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }
    /* when nprocs == 1, ncp->collective_fh == ncp->independent_fh */
    if (ncp->nprocs == 1) return NC_NOERR;

    /* When intra-node aggregation is enabled, non-aggregator's
     * ncp->collective_fh is always MPI_FILE_NULL. When disabled,
     * ncp->collective_fh on all ranks is never MPI_FILE_NULL as collective
     * mode is default in PnetCDF.
     */
    if (ncp->collective_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync, (ncp->collective_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }

    /* Barrier is not necessary ...
      TRACE_COMM(MPI_Barrier)(ncp->comm);
     */

    return NC_NOERR;
}

/*----< ncmpio_file_set_view() >---------------------------------------------*/
/* This subroutine is collective when using MPI-IO driver. When using GIO
 * driver, this subroutine is independent.
 *
 * This subroutine is used by vard APIs and called only when using MPI-IO driver.
 */
int
ncmpio_file_set_view(NC           *ncp,
                     MPI_Datatype  filetype,
                     MPI_Aint      npairs,
#ifdef HAVE_MPI_LARGE_COUNT
                     MPI_Count    *offsets,
                     MPI_Count    *lengths
#else
                     MPI_Offset   *offsets,
                     int          *lengths
#endif
)
{
    if (ncp->driver != PNC_DRIVER_MPIIO) assert(0);

    char *mpi_name;
    int err, mpireturn, status=NC_NOERR;
    MPI_File fh;

    if (ncp->driver == PNC_DRIVER_GIO) {
        /* GIO library does not have an API for setting a file view. The
         * flattened file view's offset-length pairs are directly passed to
         * GIO_read/write APIs, because GIO does not support reuse of a
         * file view. Using flattened offset-length pairs avoids constructing
         * and flattening a filetype.
         */
        ncp->file_view.count = npairs;
        ncp->file_view.off = offsets;
        ncp->file_view.len = lengths;

        return NC_NOERR;
    }

    /* Now, ncp->driver == PNC_DRIVER_MPIIO, i.e. using MPI-IO. */
    int to_free_filetype=0;

    /* when ncp->nprocs == 1, ncp->collective_fh == ncp->independent_fh */
    fh = (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
       ? ncp->collective_fh : ncp->independent_fh;

#ifdef PNETCDF_DEBUG
    /* When using MPI-IO driver and INA disabled, all processes hold valid
     * ncp->collective_fh (and ncp->independent_fh once entering independent
     * data mode).
     *
     * When using MPI-IO driver and INA enabled, all INA aggregators processes
     * hold valid ncp->collective_fh (and ncp->independent_fh once entering
     * independent data mode). Non-INA aggregators have only valid
     * ncp->independent_fh once entering independent data mode. Its
     * ncp->collective_fh should always be MPI_FILE_NULL.
     * */
    if (ncp->num_aggrs_per_node > 0 && !ncp->comm_attr.is_ina_aggr)
        assert(ncp->collective_fh == MPI_FILE_NULL);
#endif

    if (fh == MPI_FILE_NULL) /* This process is not an INA aggregator. */
        return NC_NOERR;

    if (npairs == 0) /* zero-sized requests */
        filetype = MPI_BYTE;
    else {
#ifdef HAVE_MPI_LARGE_COUNT
        /* construct file view */
        mpireturn = MPI_Type_create_hindexed_c(npairs, lengths, offsets,
                                               MPI_BYTE, &filetype);
#else
        assert(sizeof(*offsets) == sizeof(MPI_Aint));
        /* construct file view */
        mpireturn = MPI_Type_create_hindexed(npairs, lengths,
                                             (MPI_Aint*)offsets,
                                             MPI_BYTE, &filetype);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&filetype);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else
                to_free_filetype = 1;
        }
    }

    /* PnetCDF always builds filetype using flatten offset-length pairs, so
     * argument displacement is always 0.
     */
    TRACE_IO(MPI_File_set_view, (fh, 0, MPI_BYTE, filetype, "native",
                                 MPI_INFO_NULL));
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        if (status == NC_NOERR) status = err;
    }

    if (to_free_filetype)
        MPI_Type_free(&filetype);

    return status;
}

/*----< type_create_hindexed() >---------------------------------------------*/
static
int type_create_hindexed(PNCIO_View    view,
                         MPI_Datatype *newType)
{
    char *mpi_name;
    int status=NC_NOERR, err=NC_NOERR, mpireturn;

    *newType = MPI_DATATYPE_NULL;

    /* construct a buftype */
#ifdef HAVE_MPI_LARGE_COUNT
    /* TODO: MPI_Type_create_hindexed_c
     *       view.count should be of type MPI_Count
     *       view.len   should be of type MPI_Count
     *       view.off   should be of type MPI_Count
     */
    mpireturn = MPI_Type_create_hindexed_c(view.count,
                                           view.len,
                                           view.off,
                                           MPI_BYTE, newType);
    mpi_name = "MPI_Type_create_hindexed_c";
#else
    MPI_Aint *disp;
#if SIZEOF_MPI_AINT == SIZEOF_MPI_OFFSET
    disp = (MPI_Aint*) view.off;
#else
    disp = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * view.count);
    for (j=0; j<view.count; j++)
        disp[j] = (MPI_Aint)view.off[j];
#endif
    /* TODO: MPI_Type_create_hindexed
     *       view.count should be of type int
     *       view.len   should be of type int
     *       view.off   should be of type MPI_Aint
     */
    mpireturn = MPI_Type_create_hindexed(view.count,
                                         view.len,
                                         disp,
                                         MPI_BYTE, newType);
    mpi_name = "MPI_Type_create_hindexed";
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
    NCI_Free(disp);
#endif
#endif
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        /* return the first encountered error if there is any */
        if (status == NC_NOERR) status = err;
    }
    else {
        mpireturn = MPI_Type_commit(newType);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    return status;
}

/*----< ncmpio_file_read() >-------------------------------------------------*/
/* For zero-sized requests, file_view.count == buf_view.count == 0.
 * For non-zero sized requests, file_view.count > 0 and buf_view.count > 0 and
 * both file_view's off and len should not be NULL. The same for buf_view.
 * Accumulated amount of file_view and buf_view should be equal.
 *
 * This subroutine may be called from either collective or independent data
 * mode. Argument 'coll_indep' indicates the caller intends to perform a
 * collective or independent get request, regardless of which data mode the
 * program is currently on.
 */
MPI_Offset
ncmpio_file_read(NC         *ncp,
                 int         coll_indep, /* NC_REQ_INDEP or NC_REQ_COLL */
                 void       *buf,
                 PNCIO_View  file_view,
                 PNCIO_View  buf_view)
{
    int i, status=NC_NOERR, err=NC_NOERR;
    MPI_Offset rlen=0;
    MPI_Count off_zero=0, f_amnt, b_amnt;
    PNCIO_View orig_buf_view;

    /* If zero-sized request and independent write, this rank can return now. */
    if (buf_view.count == 0 && coll_indep == NC_REQ_INDEP)
        return NC_NOERR;

// printf("\n%s at %d: file_view count %lld off[0] %lld len[0] %lld buf_view count %lld off[0] %lld len[0] %lld\n",__func__,__LINE__,file_view.count,file_view.off[0], file_view.len[0], buf_view.count,buf_view.off[0], buf_view.len[0]);
    /* buf_view.count is the number of offset-length pairs */

    /* Calculate file_view amount in bytes, may be > NC_MAX_INT */
    for (f_amnt=0, i=0; i<file_view.count; i++)
        f_amnt += file_view.len[i];

    /* Calculate buf_view amount in bytes, may be > NC_MAX_INT */
    for (b_amnt=0, i=0; i<buf_view.count; i++)
        b_amnt += buf_view.len[i];

#ifdef PNETCDF_DEBUG
    assert(f_amnt == b_amnt);
#endif

    void *xbuf=buf;

    /* Save the original buf_view.count, as it may be modified when
     * b_amnt <= ibuf_size, original buf_view is required to unpack the read
     * data to user read buffer.
     */
    orig_buf_view = buf_view;

#ifndef HAVE_MPI_LARGE_COUNT
    if (b_amnt > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
        fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buffer size="OFFFMT"\n",
                ncp->rank, __func__,__LINE__,b_amnt);
#endif
        if (coll_indep == NC_REQ_COLL) {
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
            /* write nothing, but participate the collective call */
            buf_view.count = 0;
        }
        else
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    }
#endif

    if (buf_view.count > 1 && b_amnt <= ncp->ibuf_size) {
        /* If this read buffer is noncontiguous and amount is less than
         * ncp->ibuf_size, we allocate a temporary contiguous buffer and use it
         * to read from the file. Later it is unpacked to user buffer. As some
         * MPI, e.g. Cray on KNL, can be significantly slow when the read
         * buffer is noncontiguous. The only case of read buffer being
         * noncontiguous is when nonblocking API ncmpi_wait/wait_all() is
         * called and INA is disabled.
         *
         * Note ncp->ibuf_size is never > NC_MAX_INT.
         */
        xbuf = NCI_Malloc(b_amnt);

        /* mark buf_view is contiguous */
        buf_view.count = 1;
        buf_view.off = &off_zero;
        buf_view.len = &b_amnt;
    }

    if (ncp->driver == PNC_DRIVER_MPIIO) {
        char *mpi_name;
        int set_file_view, mpireturn;
        MPI_Offset disp;
        MPI_File fh;
        MPI_Status mpistatus;
        MPI_Datatype fileType=MPI_BYTE, bufType=MPI_DATATYPE_NULL;

        memset(&mpistatus, 0, sizeof(MPI_Status));

        /* when ncp->nprocs == 1, ncp->collective_fh == ncp->independent_fh */
        fh = (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
           ? ncp->collective_fh : ncp->independent_fh;

        /* When in collective data mode, if a process makes an independent put
         * call, skip setting the file view. This must come from subroutines
         * that would like to access file header or perform data section
         * movement at ncmpi_enddef(), due to new metadata was added. In these
         * cases, their file views are always contiguous.
         *
         * In all other cases, we must set the file view.
         */
        if (!fIsSet(ncp->flags, NC_MODE_INDEP) && coll_indep == NC_REQ_INDEP) {
            set_file_view = 0;
            assert(file_view.count == 1);
        }
        else
            set_file_view = 1;

        disp = 0;
        if (file_view.count == 1) { /* no need to create fileType */
            disp = file_view.off[0];
        }
        else if (file_view.count > 1) {
            /* Construct a file type for setting the file view. */
            err = type_create_hindexed(file_view, &fileType);
            if (err != NC_NOERR && status == NC_NOERR) status = err;
        }

        if (set_file_view) {
            TRACE_IO(MPI_File_set_view, (fh, 0, MPI_BYTE, fileType,
                                         "native", MPI_INFO_NULL));
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                if (status == NC_NOERR) status = err;
            }
        }
        if (fileType != MPI_BYTE) MPI_Type_free(&fileType);

        /* Construct a derived data type for MPI_File_read_xxx call, no matter
         * if buf_view is contiguous or not.
         */
        err = type_create_hindexed(buf_view, &bufType);
        if (err != NC_NOERR && status == NC_NOERR) status = err;

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL) {
            /* call MPI collective read */
            TRACE_IO(MPI_File_read_at_all, (fh, disp, xbuf, 1, bufType,
                                            &mpistatus));
        }
        else {
            /* call MPI independent read */
            TRACE_IO(MPI_File_read_at, (fh, disp, xbuf, 1, bufType,
                                        &mpistatus));
        }

        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE)
            if (status == NC_NOERR) status = err;
        }
        else
            /* update the number of bytes read */
            rlen = get_count(&mpistatus, MPI_BYTE);
// printf("%s at %d: rlen %lld\n",__func__,__LINE__,rlen);

        if (set_file_view) /* reset the file view to entire file visible */
            MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

        if (bufType != MPI_DATATYPE_NULL) MPI_Type_free(&bufType);
    }
#ifdef ENABLE_GIO
    else if (ncp->driver == PNC_DRIVER_GIO) {

        if (ncp->gio_fh == GIO_FILE_NULL) {
            /* If this process has not opened the file, this process must be a
             * non-INA aggregator calling from an independent get API.
             */
            assert(coll_indep == NC_REQ_INDEP);

            /* Check file open mode */
            int amode = fIsSet(ncp->flags, NC_MODE_RDONLY) ? O_RDONLY : O_RDWR;
            err = GIO_open(MPI_COMM_SELF, ncp->path, amode, ncp->mpiinfo,
                           &ncp->gio_fh);
            if (err != GIO_NOERR) {
                err = ncmpii_error_gio2nc(err, "GIO_open");
                DEBUG_FOPEN_ERROR(err);
                if (status == NC_NOERR) status = err;
            }
        }

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
            rlen = GIO_read_at_all(ncp->gio_fh, xbuf,
                                   file_view.count,
                                   file_view.off,
                                   file_view.len,
                                   buf_view.count,
                                   buf_view.off,
                                   buf_view.len);
        else
            rlen = GIO_read_at(ncp->gio_fh, xbuf,
                                   file_view.count,
                                   file_view.off,
                                   file_view.len,
                                   buf_view.count,
                                   buf_view.off,
                                   buf_view.len);

        if (status == NC_NOERR && rlen < 0) status = (int)rlen;
    }
#endif
    else
        status = NC_EDRIVER;

    if (xbuf != buf) { /* unpack contiguous xbuf to noncontiguous buf */
        char *in_ptr, *out_ptr;
        in_ptr = xbuf;

        for (i=0; i<orig_buf_view.count; i++) {
            out_ptr = (char*)buf + orig_buf_view.off[i];
            memcpy(out_ptr, in_ptr, orig_buf_view.len[i]);
            in_ptr += orig_buf_view.len[i];
        }
        NCI_Free(xbuf);
    }

    /* update the number of bytes read since file open */
    if (rlen >= 0) ncp->get_size += rlen;

    return (status == NC_NOERR) ? rlen : status;
}

/*----< ncmpio_file_write() >------------------------------------------------*/
/* For zero-sized requests, file_view.count == buf_view.count == 0.
 * For non-zero sized requests, file_view.count > 0 and buf_view.count > 0 and
 * both file_view's off and len should not be NULL. The same for buf_view.
 * Accumulated amount of file_view and buf_view should be equal.
 *
 * This subroutine may be called from either collective or independent data
 * mode. Argument 'coll_indep' indicates the caller intends to perform a
 * collective or independent put request, regardless of which data mode the
 * program is currently on.
 */
MPI_Offset
ncmpio_file_write(NC         *ncp,
                  int         coll_indep, /* NC_REQ_INDEP or NC_REQ_COLL */
                  const void *buf,
                  PNCIO_View  file_view,
                  PNCIO_View  buf_view)
{
    int i, status=NC_NOERR, err=NC_NOERR;
    MPI_Offset wlen=0;
    MPI_Count off_zero=0, f_amnt, b_amnt;

    /* If zero-sized request and independent write, this rank can return now. */
    if (buf_view.count == 0 && coll_indep == NC_REQ_INDEP)
        return NC_NOERR;

    /* buf_view.count is the number of offset-length pairs */

// if (file_view.off[0]==600 && file_view.len[0]==0) assert(0);

// printf("\n%s at %d: file_view count %lld off[0] %lld len[0] %lld buf_view count %lld off[0] %lld len[0] %lld\n",__func__,__LINE__,file_view.count,file_view.off[0], file_view.len[0], buf_view.count,buf_view.off[0], buf_view.len[0]);

    /* Calculate file_view amount in bytes, may be > NC_MAX_INT */
    for (f_amnt=0, i=0; i<file_view.count; i++)
        f_amnt += file_view.len[i];

    /* Calculate buf_view amount in bytes, may be > NC_MAX_INT */
    for (b_amnt=0, i=0; i<buf_view.count; i++)
        b_amnt += buf_view.len[i];

#ifdef PNETCDF_DEBUG
    assert(f_amnt == b_amnt);
#endif

    void *xbuf= (void*)buf;

    if (buf_view.count > 1 && b_amnt <= ncp->ibuf_size) {
        /* If this write buffer is noncontiguous and amount is less than
         * ncp->ibuf_size, we pack it into a temporary contiguous buffer and
         * use it to write to the file. As some MPI, e.g. Cray on KNL, can be
         * significantly slow when the write buffer is noncontiguous. The only
         * case of write buffer being noncontiguous is when nonblocking API
         * ncmpi_wait/wait_all() is called and INA is disabled.
         *
         * Note ncp->ibuf_size is never > NC_MAX_INT.
         */
        char *in_ptr, *out_ptr;
        xbuf = NCI_Malloc(b_amnt);
        out_ptr = xbuf;

        for (i=0; i<buf_view.count; i++) {
            in_ptr = (char*)buf + buf_view.off[i];
            memcpy(out_ptr, in_ptr, buf_view.len[i]);
            out_ptr += buf_view.len[i];
        }
        /* mark buf_view is contiguous */
        buf_view.count = 1;
        buf_view.off = &off_zero;
        buf_view.len = &b_amnt;
    }

    if (ncp->driver == PNC_DRIVER_MPIIO) {
        char *mpi_name;
        int set_file_view, mpireturn;
        MPI_Offset disp;
        MPI_File fh;
        MPI_Status mpistatus;
        MPI_Datatype fileType=MPI_BYTE, bufType=MPI_DATATYPE_NULL;

        memset(&mpistatus, 0, sizeof(MPI_Status));

        /* when ncp->nprocs == 1, ncp->collective_fh == ncp->independent_fh */
        fh = (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
           ? ncp->collective_fh : ncp->independent_fh;

        /* When in collective data mode, if a process makes an independent put
         * call, skip setting the file view. This must come from subroutines
         * that would like to access file header or perform data section
         * movement at ncmpi_enddef(), due to new metadata was added. In these
         * cases, their file views are always contiguous.
         *
         * In all other cases, we must set the file view.
         */
        if (!fIsSet(ncp->flags, NC_MODE_INDEP) && coll_indep == NC_REQ_INDEP) {
            set_file_view = 0;
            assert(file_view.count == 1);
        }
        else
            set_file_view = 1;

        disp = 0;
        if (file_view.count == 1) { /* no need to create fileType */
            disp = file_view.off[0];
        }
        else if (file_view.count > 1) {
            /* Construct a file type for setting the file view. */
            err = type_create_hindexed(file_view, &fileType);
            if (err != NC_NOERR && status == NC_NOERR) status = err;
        }

        if (set_file_view) {
            TRACE_IO(MPI_File_set_view, (fh, 0, MPI_BYTE, fileType,
                                         "native", MPI_INFO_NULL));
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                if (status == NC_NOERR) status = err;
            }
        }
        if (fileType != MPI_BYTE) MPI_Type_free(&fileType);

        /* Construct a derived data type for MPI_File_write_xxx call, no matter
         * if buf_view is contiguous or not.
         */
        err = type_create_hindexed(buf_view, &bufType);
        if (err != NC_NOERR && status == NC_NOERR) status = err;

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL) {
            /* call MPI collective write */
            TRACE_IO(MPI_File_write_at_all, (fh, disp, xbuf, 1, bufType,
                                             &mpistatus));
        }
        else {
            /* call MPI independent write */
            TRACE_IO(MPI_File_write_at, (fh, disp, xbuf, 1, bufType,
                                         &mpistatus));
        }

        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE)
            if (status == NC_NOERR) status = err;
        }
        else {
            /* update the number of bytes written */
            wlen = get_count(&mpistatus, MPI_BYTE);
// printf("%s at %d: wlen %lld\n",__func__,__LINE__,wlen);
        }

        if (set_file_view) /* reset the file view to entire file visible */
            MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                              MPI_INFO_NULL);

        if (bufType != MPI_DATATYPE_NULL) MPI_Type_free(&bufType);
    }
#ifdef ENABLE_GIO
    else if (ncp->driver == PNC_DRIVER_GIO) {

        if (ncp->gio_fh == GIO_FILE_NULL) {
            /* If this process has not opened the file, this process must be a
             * non-INA aggregator calling from an independent put API.
             */
            assert(coll_indep == NC_REQ_INDEP);

            /* Check file open mode */
            int amode = fIsSet(ncp->flags, NC_MODE_RDONLY) ? O_RDONLY : O_RDWR;
            err = GIO_open(MPI_COMM_SELF, ncp->path, amode, ncp->mpiinfo,
                           &ncp->gio_fh);
            if (err != GIO_NOERR) {
                err = ncmpii_error_gio2nc(err, "GIO_open");
                DEBUG_FOPEN_ERROR(err);
                if (status == NC_NOERR) status = err;
            }
        }

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
            wlen = GIO_write_at_all(ncp->gio_fh, xbuf,
                                    file_view.count,
                                    file_view.off,
                                    file_view.len,
                                    buf_view.count,
                                    buf_view.off,
                                    buf_view.len);
        else
            wlen = GIO_write_at(ncp->gio_fh, xbuf,
                                    file_view.count,
                                    file_view.off,
                                    file_view.len,
                                    buf_view.count,
                                    buf_view.off,
                                    buf_view.len);

        if (status == NC_NOERR && wlen < 0) status = (int)wlen;
    }
#endif
    else
        status = NC_EDRIVER;

    if (xbuf != buf) NCI_Free(xbuf);

    /* update the number of bytes written since file open */
    if (wlen >= 0) ncp->put_size += wlen;

    return (status == NC_NOERR) ? wlen : status;
}

