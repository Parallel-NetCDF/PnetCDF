/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/errno.h>
#include <unistd.h>   /* pwrite() */

#include <mpi.h>

#include "pncio.h"

#ifdef WKL_DEBUG
int first_ost_id;
#endif

/*----< PNCIO_WriteContig() >-------------------------------------------------*/
MPI_Offset PNCIO_WriteContig(PNCIO_File *fd,
                             const void *buf,
                             MPI_Offset  w_size,
                             MPI_Offset  offset)
{
    ssize_t err = 0;
    size_t w_count;
    MPI_Offset bytes_xfered = 0;
    char *p;

    if (w_size == 0) return NC_NOERR;

// printf("%s at %d: pwrite offset=%lld w_size=%lld\n",__func__,__LINE__,offset,w_size);
#ifdef WKL_DEBUG
int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

MPI_Offset ost_id = (offset / fd->hints->striping_unit) % fd->hints->striping_factor;
    if (first_ost_id == -1) {
        first_ost_id = ost_id;
        // printf("%2d %s file %s First pwrite offset=%lld OST %d\n",rank,__func__,fd->filename,offset,first_ost_id);
    }
    else if (ost_id != first_ost_id)
        printf("%2d Error: %s pwrite offset=%lld w_size=%lld ost_id=%lld not same 1st ost %d\n",rank,__func__,offset,w_size,ost_id,first_ost_id);

printf("%s line %d: disp=%lld offset=%lld count=%ld bufType_size=%d w_size=%lld\n",__func__,__LINE__,fd->disp,offset,count,bufType_size,w_size);

    printf("%2d %s line %d pwrite offset=%lld w_size=%lld\n",rank,__func__,__LINE__,offset,w_size);
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif
    p = (char *) buf;
    while (bytes_xfered < w_size) {
        w_count = w_size - bytes_xfered;
        err = pwrite(fd->fd_sys, p, w_count, offset + bytes_xfered);
        if (err == -1)
            goto ioerr;
        if (err == 0)
            break;
        bytes_xfered += err;
        p += err;
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->write_timing[2] += MPI_Wtime() - timing;
#endif

ioerr:
    if (err == -1)
        bytes_xfered = ncmpii_error_posix2nc("pwrite");

    return bytes_xfered;
}

/*----< file_write() >-------------------------------------------------------*/
/* This is an independent call. */
static
MPI_Offset file_write(PNCIO_File *fd,
                      MPI_Offset  offset,
                      const void *buf,
                      PNCIO_View  buf_view)
{
    MPI_Offset w_len;

    if (buf_view.size == 0) /* zero-sized request */
        return NC_NOERR;

assert(fd->filetype == MPI_BYTE);

    if (buf_view.is_contig && fd->flat_file.is_contig) {
        if (fd->flat_file.count > 0) offset += fd->flat_file.off[0];
        w_len = PNCIO_WriteContig(fd, buf, buf_view.size, offset);
    }
    else if (fd->file_system == PNCIO_LUSTRE)
        w_len = PNCIO_LUSTRE_WriteStrided(fd, buf, buf_view, offset);
    else if (fd->file_system == PNCIO_UFS)
        w_len = PNCIO_GEN_WriteStrided(fd, buf, buf_view, offset);
    else
        return NC_EFSTYPE;

    return w_len; /* when w_len < 0, it is an NetCDF error code */
}

/*----< PNCIO_File_write_at() >-----------------------------------------------*/
/* This is an independent call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
MPI_Offset PNCIO_File_write_at(PNCIO_File *fh,
                               MPI_Offset  offset,
                               const void *buf,
                               PNCIO_View  buf_view)
{
    assert(fh != NULL);

    if (buf_view.size == 0) /* zero-sized request */
        return NC_NOERR;

    if (buf_view.size < 0) return NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        return NC_EPERM;

    return file_write(fh, offset, buf, buf_view);
}

/*----< PNCIO_File_write_at_all() >-------------------------------------------*/
/* This is a collective call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
MPI_Offset PNCIO_File_write_at_all(PNCIO_File *fh,
                                   MPI_Offset  offset,
                                   const void *buf,
                                   PNCIO_View  buf_view)
{
    int err=NC_NOERR;
    MPI_Offset w_len;

    assert(fh != NULL);

    if (buf_view.size < 0) err = NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY && err == NC_NOERR)
        err = NC_EPERM;

    if (fh->file_system == PNCIO_LUSTRE)
        w_len = PNCIO_LUSTRE_WriteStridedColl(fh, buf, buf_view, offset);
    else if (fh->file_system == PNCIO_UFS)
        w_len = PNCIO_GEN_WriteStridedColl(fh, buf, buf_view, offset);
    else
        return NC_EFSTYPE;

    return (err == NC_NOERR) ? w_len : err;
}


