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

/*----< get_count() >--------------------------------------------------------*/
/* This subroutine is independent. On success, the number of bytes read/written
 * is returned (zero indicates nothing was read/written). Like POSIX read()/
 * write(), it is not an error if this number is smaller than the number of
 * bytes requested. On error, a negative value, an NC error code, is returned.
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
/*
 * This function is independent.
 */
/* TODO: move check count against MAX_INT and call _c API */
MPI_Offset
ncmpio_file_read_at(NC         *ncp,
                    MPI_Offset  offset,
                    void       *buf,
                    PNCIO_View  buf_view)
{
    int err=NC_NOERR, mpireturn;
    MPI_Offset amnt=0;
    MPI_Status mpistatus;

    /* explicitly initialize mpistatus object to 0. For zero-length read/write,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around. See MPICH ticket:
     * https://trac.mpich.org/projects/mpich/ticket/2332
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        MPI_File fh;

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

        if (fh == MPI_FILE_NULL) return 0;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count count = (buf_view.is_contig) ? buf_view.size : 1;

        TRACE_IO(MPI_File_read_at_c, (fh, offset, buf, count, buf_view.type,
                                      &mpistatus));
#else
        int count = (buf_view.is_contig) ? buf_view.size : 1;

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
    else if (ncp->pncio_fh != NULL)
        amnt = PNCIO_File_read_at(ncp->pncio_fh, offset, buf, buf_view);

    /* update the number of bytes read since file open */
    if (amnt >= 0) ncp->get_size += amnt;
    /* else: ignore if error, as this error is not fatal */

    return amnt;
}

/*----< ncmpio_file_read_at_all() >------------------------------------------*/
/*
 * This function is collective.
 */
MPI_Offset
ncmpio_file_read_at_all(NC         *ncp,
                        MPI_Offset  offset,
                        void       *buf,
                        PNCIO_View  buf_view)
{
    int err=NC_NOERR, mpireturn;
    MPI_Offset amnt=0;
    MPI_Status mpistatus;

    /* Explicitly initialize mpistatus object to 0. For zero-length read/write,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around. See MPICH ticket:
     * https://trac.mpich.org/projects/mpich/ticket/2332
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        MPI_File fh;

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

        if (fh == MPI_FILE_NULL) return 0;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count count = (buf_view.is_contig) ? buf_view.size : 1;

        TRACE_IO(MPI_File_read_at_all_c, (fh, offset, buf, count,
                                          buf_view.type, &mpistatus));
#else
        int count = (buf_view.is_contig) ? buf_view.size : 1;

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
    else if (ncp->pncio_fh != NULL)
        amnt = PNCIO_File_read_at_all(ncp->pncio_fh, offset, buf, buf_view);

    /* update the number of bytes read since file open */
    if (amnt >= 0) ncp->get_size += amnt;
    /* else: ignore if error, as this error is not fatal */

    return amnt;
}

/*----< ncmpio_file_write_at() >---------------------------------------------*/
/*
 * This function is independent.
 */
MPI_Offset
ncmpio_file_write_at(NC         *ncp,
                     MPI_Offset  offset,
                     const void *buf,
                     PNCIO_View  buf_view)
{
    int err=NC_NOERR, mpireturn;
    MPI_Offset amnt=0;
    MPI_Status mpistatus;

    /* Explicitly initialize mpistatus object to 0. For zero-length read/write,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around. See MPICH ticket:
     * https://trac.mpich.org/projects/mpich/ticket/2332
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        MPI_File fh;

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

        if (fh == MPI_FILE_NULL) return 0;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count count = (buf_view.is_contig) ? buf_view.size : 1;

        TRACE_IO(MPI_File_write_at_c, (fh, offset, buf, count, buf_view.type,
                                       &mpistatus));
#else
        int count = (buf_view.is_contig) ? buf_view.size : 1;

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
    else if (ncp->pncio_fh != NULL)
        amnt = PNCIO_File_write_at(ncp->pncio_fh, offset, buf, buf_view);

    /* update the number of bytes written since file open */
    if (amnt >= 0) ncp->put_size += amnt;
    /* else: ignore if error, as this error is not fatal */

    return amnt;
}

/*----< ncmpio_file_write_at_all() >-----------------------------------------*/
/*
 * This function is collective.
 */
MPI_Offset
ncmpio_file_write_at_all(NC         *ncp,
                         MPI_Offset  offset,
                         const void *buf,
                         PNCIO_View  buf_view)
{
    int err=NC_NOERR, mpireturn;
    MPI_Offset amnt=0;
    MPI_Status mpistatus;

    /* explicitly initialize mpistatus object to 0. For zero-length read/write,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around. See MPICH ticket:
     * https://trac.mpich.org/projects/mpich/ticket/2332
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        MPI_File fh;

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

        if (fh == MPI_FILE_NULL) return 0;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count count = (buf_view.is_contig) ? buf_view.size : 1;

        TRACE_IO(MPI_File_write_at_all_c, (fh, offset, buf, count,
                                           buf_view.type, &mpistatus));
#else
        int count = (buf_view.is_contig) ? buf_view.size : 1;

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
    else if (ncp->pncio_fh != NULL)
        amnt = PNCIO_File_write_at_all(ncp->pncio_fh, offset, buf, buf_view);

    /* update the number of bytes written since file open */
    if (amnt >= 0) ncp->put_size += amnt;
    /* else: ignore if error, as this error is not fatal */

    return amnt;
}

/*----< ncmpio_getput_zero_req() >-------------------------------------------*/
/* This function is called when this process has zero-length I/O request and
 * must participate all the MPI collective calls involved in the collective
 * APIs and wait_all(), which include setting fileview, collective read/write,
 * another setting fileview.
 *
 * This function is collective.
 */
int
ncmpio_getput_zero_req(NC *ncp, int reqMode)
{
    int err, status=NC_NOERR;
    MPI_Offset rlen, wlen;
    PNCIO_View buf_view;

    buf_view.size = 0;

    /* When intra-node aggregation is enabled, non-aggregators do not access
     * the file.
     */
    if (ncp->num_aggrs_per_node > 0 && ncp->rank != ncp->my_aggr)
        return NC_NOERR;

    /* do nothing if this came from an independent API */
    if (fIsSet(reqMode, NC_REQ_INDEP)) return NC_NOERR;

    err = ncmpio_file_set_view(ncp, 0, MPI_BYTE, 0, NULL, NULL);
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

    /* Reset fileview. Note fileview is never reused in PnetCDF */
    ncmpio_file_set_view(ncp, 0, MPI_BYTE, 0, NULL, NULL);

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     */

    return status;
}

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
            }
            else
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
#endif

// printf("%s at %d: buf_view count=%lld type=%s size=%lld\n",__func__,__LINE__, buf_view.count, (buf_view.type==MPI_BYTE)?"MPI_BYTE":"NOT MPI_BYTE", buf_view.size);

        if (!buf_view.is_contig && buf_view.size <= ncp->ibuf_size) {
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
            buf_view.type = MPI_BYTE;
            buf_view.is_contig = 1;
        }

        if (!buf_view.is_contig && ncp->fstype == PNCIO_FSTYPE_MPIIO) {
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
            rlen = ncmpio_file_read_at_all(ncp, offset, xbuf, buf_view);
        else
            rlen = ncmpio_file_read_at(ncp, offset, xbuf, buf_view);
        if (status == NC_NOERR && rlen < 0) status = (int)rlen;

        if (xbuf != buf) { /* unpack contiguous xbuf to noncontiguous buf */
            char *in_ptr, *out_ptr;
            in_ptr = xbuf;

#if 0
 long long *wkl, nelems; int j;
 wkl = (long long*) malloc(buf_view.size);
 nelems=buf_view.size/8;
 memcpy(wkl, xbuf, nelems*8); ncmpii_in_swapn(wkl, nelems, 8);
 printf("%s at %d: nelems=%lld xbuf=(%p) ",__func__,__LINE__, nelems, xbuf);
 for (i=0; i<nelems; i++) printf(" %lld",wkl[i]);
 printf("\n");
 free(wkl);
#endif

assert(buf != NULL);
            for (i=0; i<buf_view.count; i++) {
                out_ptr = (char*)buf + buf_view.off[i]; // - buf_view.off[0]);
                memcpy(out_ptr, in_ptr, buf_view.len[i]);
                in_ptr += buf_view.len[i];
#if 0
 wkl = (long long*) malloc(buf_view.len[i]);
 nelems=buf_view.len[i]/8;
 memcpy(wkl, out_ptr, nelems*8); ncmpii_in_swapn(wkl, nelems, 8);
 printf("%s at %d: buf_view.count=%lld i=%d nelems=%lld out_ptr=(%p) ",__func__,__LINE__, buf_view.count,i,nelems, out_ptr);
 for (j=0; j<nelems; j++) printf(" %lld",wkl[j]);
 printf("\n");
 free(wkl);
#endif
            }
            NCI_Free(xbuf);
        }
        if (to_free_buftype)
            MPI_Type_free(&buf_view.type);

    } else { /* NC_REQ_WR */
        void *xbuf=buf;

        if (!buf_view.is_contig && buf_view.size <= ncp->ibuf_size) {
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
assert(buf != NULL);
// printf("%s at %d: buf_view count=%lld size=%lld\n",__func__,__LINE__, buf_view.count,buf_view.size);

#if 0
printf("%s at %d: buf = %p\n",__func__,__LINE__, buf);
printf("%s at %d: buf_view count=%lld off=%lld %lld len=%lld %lld\n",__func__,__LINE__, buf_view.count,buf_view.off[0],buf_view.off[1],buf_view.len[0],buf_view.len[1]);
int wkl[21];
#endif
            for (i=0; i<buf_view.count; i++) {
                in_ptr = (char*)buf + buf_view.off[i];
#if 0
memcpy(wkl, in_ptr, buf_view.len[i]);
ncmpii_in_swapn(wkl, buf_view.len[i]/4, 4);
printf("%s at %d: [%lld] in_ptr=(%p) %d %d %d %d %d\n",__func__,__LINE__, buf_view.len[i]/4, in_ptr, wkl[0],wkl[1],wkl[2],wkl[3],wkl[4]);
#endif
                memcpy(out_ptr, in_ptr, buf_view.len[i]);
                out_ptr += buf_view.len[i];
            }
            /* mark the xbuf is contiguous */
            buf_view.type = MPI_BYTE;
            buf_view.is_contig = 1;
#if 0
memcpy(wkl, xbuf, 84);
ncmpii_in_swapn(wkl, 21, 4);
printf("%s at %d: size=%lld xbuf=(%p) %d %d %d %d %d\n",__func__,__LINE__, buf_view.size, xbuf, wkl[0],wkl[1],wkl[2],wkl[3],wkl[4]);
printf("%s at %d: wkl[15] = %d %d %d %d %d\n",__func__,__LINE__, wkl[15],wkl[16],wkl[17],wkl[18],wkl[19]);
#endif
        }

        if (!buf_view.is_contig && ncp->fstype == PNCIO_FSTYPE_MPIIO) {
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
    /* Reset fileview. Note fileview is never reused in PnetCDF */
    ncmpio_file_set_view(ncp, 0, MPI_BYTE, 0, NULL, NULL);

    return status;
}

/*----< ncmpio_file_close() >------------------------------------------------*/
/*
 * This function is collective.
 */
int
ncmpio_file_close(NC *ncp)
{
    int err=NC_NOERR;

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
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
    else {
        /* When intra-node aggregation is enabled, only aggregators have a
         * non-NULL ncp->pncio_fh and non-aggregators has pncio_fh == NULL.
         */
        if (ncp->pncio_fh != NULL) {
            err = PNCIO_File_close(ncp->pncio_fh);
            NCI_Free(ncp->pncio_fh);
            ncp->pncio_fh = NULL;
        }
    }

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
        if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
            char *mpi_name;
            int mpireturn;
            TRACE_IO(MPI_File_delete, ((char *)ncp->path, ncp->mpiinfo));
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }
        else
            err = PNCIO_File_delete(ncp->path);
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

    if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
        if (ncp->pncio_fh == NULL)
            return NC_NOERR;
        return PNCIO_File_sync(ncp->pncio_fh);
    }

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
/* This subroutine is collective when using MPI-IO. When using internal PNCIO
 * driver, this subroutine is independent.
 */
int
ncmpio_file_set_view(const NC     *ncp,
                     MPI_Offset    disp,    /* IN/OUT */
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
    char *mpi_name;
    int err, mpireturn, status=NC_NOERR;
    MPI_File fh;

assert(filetype == MPI_BYTE);
assert(disp == 0);

    if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
        /* Skip setting fileview for ranks whose pncio_fh is NULL */
        if (ncp->pncio_fh == NULL)
            return NC_NOERR;

        /* When PnetCDF's internal PNCIO driver is used, the request has been
         * flattened into offsets and lengths. Thus passed-in filetype is not
         * constructed. Note offsets and lengths are not relative to any MPI-IO
         * fileview. They will be reused in PNCIO driver as a flattened file
         * type struct, which avoids repeated work of constructing and
         * flattening the filetype.
         */
        return PNCIO_File_set_view(ncp->pncio_fh, disp, filetype, npairs,
                                   offsets, lengths);
    }

    /* Now, ncp->fstype == PNCIO_FSTYPE_MPIIO, i.e. using MPI-IO. */
    int to_free_filetype=0;

    /* when ncp->nprocs == 1, ncp->collective_fh == ncp->independent_fh */
    fh = (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
       ? ncp->collective_fh : ncp->independent_fh;

    if (fh == MPI_FILE_NULL) /* not INA aggregator */
        return NC_NOERR;

    if (npairs == 0) /* zero-sized requests */
        filetype = MPI_BYTE;
    else {
#ifdef HAVE_MPI_LARGE_COUNT
        /* construct fileview */
        mpireturn = MPI_Type_create_hindexed_c(npairs, lengths, offsets,
                                               MPI_BYTE, &filetype);
#else
        assert(sizeof(*offsets) == sizeof(MPI_Aint));
        /* construct fileview */
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

    TRACE_IO(MPI_File_set_view, (fh, disp, MPI_BYTE, filetype, "native",
                                 MPI_INFO_NULL));
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        if (status == NC_NOERR) status = err;
assert(0);
    }

    if (to_free_filetype)
        MPI_Type_free(&filetype);

    return status;
}

/*----< ncmpio_file_open() >-------------------------------------------------*/
int
ncmpio_file_open(NC         *ncp,
                 MPI_Comm    comm,
                 const char *path,
                 int         omode,
                 MPI_Info    info)
{
    int err=NC_NOERR;

    /* open file collectively */
    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        int mpireturn;
        MPI_File fh;

        TRACE_IO(MPI_File_open, (comm, path, omode, info, &fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);

        /* Now the file has been successfully opened */
        ncp->collective_fh  = fh;
        ncp->independent_fh = (ncp->nprocs > 1) ? MPI_FILE_NULL : fh;

        /* get the I/O hints used/modified by MPI-IO */
        TRACE_IO(MPI_File_get_info, (fh, &ncp->mpiinfo));
        if (mpireturn != MPI_SUCCESS)
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }
    else { /* ncp->fstype != PNCIO_FSTYPE_MPIIO */
        ncp->pncio_fh = (PNCIO_File*) NCI_Calloc(1,sizeof(PNCIO_File));

        err = PNCIO_File_open(comm, path, omode, info, ncp->pncio_fh);
        if (err != NC_NOERR) return err;

        /* Now the file has been successfully opened, obtain the I/O hints
         * used/modified by PNCIO driver.
         */
        err = PNCIO_File_get_info(ncp->pncio_fh, &ncp->mpiinfo);
    }

    return err;
}

