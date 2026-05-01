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
#include <limits.h> /* INT_MAX */

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

        if (ncp->mpio_fh_indep != ncp->mpio_fh_coll &&
            ncp->mpio_fh_indep != MPI_FILE_NULL) {
            TRACE_IO(MPI_File_close, (&ncp->mpio_fh_indep));
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }

        if (ncp->mpio_fh_coll != MPI_FILE_NULL) {
            TRACE_IO(MPI_File_close, (&ncp->mpio_fh_coll));
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }
    }
#ifdef ENABLE_GIO
    else if (ncp->driver == PNC_DRIVER_GIO) {
        if (ncp->gio_fh != NULL)
            /* When INA is enabled, non-INA aggregators' gio_fh may be NULL */
            err = GIO_close(&ncp->gio_fh);
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
        if (ncp->driver == PNC_DRIVER_MPIIO) {
            char *mpi_name;
            int mpireturn;
#ifdef MPICH_VERSION
            /* MPICH recognizes file system type acronym prefixed to the file
             * name.
             */
            TRACE_IO(MPI_File_delete, ((char *)ncp->path, ncp->info));
#else
            /* Remove the file system type prefix name if there is any, because
             * some MPI libraries do not recognize such prefix. For example,
             * when path = "lustre:/home/foo/testfile.nc", remove "lustre:" to
             * make filename pointing to "/home/foo/testfile.nc".
             */
            char *path = ncmpii_remove_file_system_type_prefix(ncp->path);
            TRACE_IO(MPI_File_delete, (path, ncp->info));
#endif
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }
#ifdef ENABLE_GIO
        else if (ncp->driver == PNC_DRIVER_GIO) {
            err = GIO_delete(ncp->path);
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

    if (ncp->mpio_fh_indep != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync, (ncp->mpio_fh_indep));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }
    /* when nprocs == 1, ncp->mpio_fh_coll == ncp->mpio_fh_indep */
    if (ncp->nprocs == 1) return NC_NOERR;

    /* When intra-node aggregation is enabled, non-aggregator's
     * ncp->mpio_fh_coll is always MPI_FILE_NULL. When disabled,
     * ncp->mpio_fh_coll on all ranks is never MPI_FILE_NULL as collective
     * mode is default in PnetCDF.
     */
    if (ncp->mpio_fh_coll != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync, (ncp->mpio_fh_coll));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }

    /* Barrier is not necessary ...
      TRACE_COMM(MPI_Barrier)(ncp->comm);
     */

    return NC_NOERR;
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
    char *xbuf;
    int i, status=NC_NOERR, err=NC_NOERR;
    MPI_Offset rlen=0;
    MPI_Count off_zero=0, f_amnt, b_amnt;
    PNCIO_View orig_buf_view;

    /* If zero-sized request and independent read, this rank can return now. */
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

#if PNETCDF_DEBUG_MODE == 1
    assert(f_amnt == b_amnt);
#endif

    /* Save the original buf_view.count, as it may be modified when
     * b_amnt <= ibuf_size, original buf_view is required to unpack the read
     * data to user read buffer.
     */
    orig_buf_view = buf_view;

#ifndef HAVE_MPI_LARGE_COUNT
    if (b_amnt > NC_MAX_INT) {
#if PNETCDF_DEBUG_MODE == 1
        fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buffer size="OFFFMT"\n",
                ncp->rank, __func__,__LINE__,b_amnt);
#endif
        if (coll_indep == NC_REQ_COLL) {
            DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
            /* read nothing, but participate the collective call */
            buf_view.count = 0;
        }
        else
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
    }
#endif

    xbuf = (char*)buf;

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
        xbuf = (char*) NCI_Malloc(b_amnt);

        /* mark buf_view is contiguous */
        buf_view.count = 1;
        buf_view.off = &off_zero;
        buf_view.len = &b_amnt;
    }

    if (ncp->driver == PNC_DRIVER_MPIIO) {
        char *xbuf_ptr, *mpi_name;
        int bufCount, set_file_view, mpireturn;
        MPI_Offset disp;
        MPI_File fh;
        MPI_Status mpistatus;
        MPI_Datatype fileType=MPI_BYTE, bufType=MPI_BYTE;

        memset(&mpistatus, 0, sizeof(MPI_Status));

        /* when ncp->nprocs == 1, ncp->mpio_fh_coll == ncp->mpio_fh_indep */
        fh = (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
           ? ncp->mpio_fh_coll : ncp->mpio_fh_indep;

        /* When in collective data mode, if a process makes an independent put
         * call, skip setting the file view. This must come from subroutines
         * that would like to access file header or perform data section
         * movement at ncmpi_enddef(), due to new metadata was added. In these
         * cases, their file views are always contiguous. However, argument
         * 'disp' to be used in MPI_File_read_at() call must be set to
         * file_view.off[0].
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
            err = ncmpio_type_create_hindexed(file_view.count, file_view.off,
                                              file_view.len, &fileType);
            if (err != NC_NOERR && status == NC_NOERR) status = err;
        }

        if (set_file_view) {
            TRACE_IO(MPI_File_set_view, (fh, disp, MPI_BYTE, fileType,
                                         "native", MPI_INFO_NULL));
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                if (status == NC_NOERR) status = err;
            }
        }
        if (fileType != MPI_BYTE) MPI_Type_free(&fileType);

        /* Construct a derived data type describing user buffer data layout to
         * be used in MPI_File_read_xxx call.
         */
        bufCount = 0;
        xbuf_ptr = xbuf;
        if (buf_view.count == 1) { /* buffer view is contiguous */
            if (buf_view.len[0] <= INT_MAX)
                bufCount = (int)buf_view.len[0];
            else if (status == NC_NOERR)
                status = NC_EINTOVERFLOW;
            xbuf_ptr += buf_view.off[0];
        }
        else if (buf_view.count > 1) {
            err = ncmpio_type_create_hindexed(buf_view.count, buf_view.off,
                                              buf_view.len, &bufType);
            if (err == NC_NOERR)
                bufCount = 1;
            else if (status == NC_NOERR)
                status = err;
        }
        /* else is for zero-sized request */

        if (set_file_view) disp = 0;

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL) {
            /* call MPI collective read */
            TRACE_IO(MPI_File_read_at_all, (fh, disp, xbuf_ptr, bufCount,
                                            bufType, &mpistatus));
        }
        else {
            /* call MPI independent read */
            TRACE_IO(MPI_File_read_at, (fh, disp, xbuf_ptr, bufCount,
                                        bufType, &mpistatus));
        }

        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD)
            if (status == NC_NOERR) status = err;
        }
        else
            /* update the number of bytes read */
            rlen = get_count(&mpistatus, MPI_BYTE);

        if (set_file_view) /* reset the file view to entire file visible */
            MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                              MPI_INFO_NULL);

        if (bufType != MPI_BYTE) MPI_Type_free(&bufType);
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
            err = GIO_open(MPI_COMM_SELF, ncp->path, amode, ncp->info,
                           &ncp->gio_fh);
            if (err != GIO_NOERR) {
                err = ncmpii_error_gio2nc(err, "GIO_open");
                DEBUG_FOPEN_ERROR(err);
                if (status == NC_NOERR) status = err;
            }
        }

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
            rlen = GIO_read_all(ncp->gio_fh, xbuf,
                                file_view.count,
                                file_view.off,
                                file_view.len,
                                buf_view.count,
                                buf_view.off,
                                buf_view.len);
        else
            rlen = GIO_read(ncp->gio_fh, xbuf,
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
    char *xbuf;
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

#if PNETCDF_DEBUG_MODE == 1
    assert(f_amnt == b_amnt);
#endif

    xbuf = (char*)buf;

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
        char *xbuf_ptr, *mpi_name;
        int bufCount, set_file_view, mpireturn;
        MPI_Offset disp;
        MPI_File fh;
        MPI_Status mpistatus;
        MPI_Datatype fileType=MPI_BYTE, bufType=MPI_BYTE;

        memset(&mpistatus, 0, sizeof(MPI_Status));

        /* when ncp->nprocs == 1, ncp->mpio_fh_coll == ncp->mpio_fh_indep */
        fh = (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
           ? ncp->mpio_fh_coll : ncp->mpio_fh_indep;

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
            err = ncmpio_type_create_hindexed(file_view.count, file_view.off,
                                              file_view.len, &fileType);
            if (err != NC_NOERR && status == NC_NOERR) status = err;
        }

        if (set_file_view) {
            TRACE_IO(MPI_File_set_view, (fh, disp, MPI_BYTE, fileType,
                                         "native", MPI_INFO_NULL));
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                if (status == NC_NOERR) status = err;
            }
        }
        if (fileType != MPI_BYTE) MPI_Type_free(&fileType);

        /* Construct a derived data type describing user buffer data layout to
         * be used in MPI_File_read_xxx call.
         */
        bufCount = 0;
        xbuf_ptr = xbuf;
        if (buf_view.count == 1) { /* buffer view is contiguous */
            if (buf_view.len[0] <= INT_MAX)
                bufCount = (int)buf_view.len[0];
            else if (status == NC_NOERR)
                status = NC_EINTOVERFLOW;
            xbuf_ptr += buf_view.off[0];
        }
        else if (buf_view.count > 1) {
            err = ncmpio_type_create_hindexed(buf_view.count, buf_view.off,
                                              buf_view.len, &bufType);
            if (err == NC_NOERR)
                bufCount = 1;
            if (status == NC_NOERR)
                status = err;
        }
        /* else is for zero-sized request */

        if (set_file_view) disp = 0;

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL) {
            /* call MPI collective write */
            TRACE_IO(MPI_File_write_at_all, (fh, disp, xbuf_ptr, bufCount,
                                             bufType, &mpistatus));
        }
        else {
            /* call MPI independent write */
            TRACE_IO(MPI_File_write_at, (fh, disp, xbuf_ptr, bufCount,
                                         bufType, &mpistatus));
        }

        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE)
            if (status == NC_NOERR) status = err;
        }
        else
            /* update the number of bytes written */
            wlen = get_count(&mpistatus, MPI_BYTE);

        if (set_file_view) /* reset the file view to entire file visible */
            MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                              MPI_INFO_NULL);

        if (bufType != MPI_BYTE) MPI_Type_free(&bufType);
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
            err = GIO_open(MPI_COMM_SELF, ncp->path, amode, ncp->info,
                           &ncp->gio_fh);
            if (err != GIO_NOERR) {
                err = ncmpii_error_gio2nc(err, "GIO_open");
                DEBUG_FOPEN_ERROR(err);
                if (status == NC_NOERR) status = err;
            }
        }

        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
            wlen = GIO_write_all(ncp->gio_fh, xbuf,
                                 file_view.count,
                                 file_view.off,
                                 file_view.len,
                                 buf_view.count,
                                 buf_view.off,
                                 buf_view.len);
        else
            wlen = GIO_write(ncp->gio_fh, xbuf,
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

