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
#include <unistd.h>   /* pread() */

#include <mpi.h>

#include "pncio.h"

/*----< PNCIO_ReadContig() >--------------------------------------------------*/
MPI_Offset PNCIO_ReadContig(PNCIO_File *fd,
                            void       *buf,
                            MPI_Offset  r_size,
                            MPI_Offset  offset)
{
    ssize_t err = 0;
    size_t r_count;
    MPI_Offset bytes_xfered = 0;
    char *p;

// printf("%s at %d: %s pread offset=%lld r_size=%lld\n",__func__,__LINE__,fd->filename,offset,r_size);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif
    p = (char *) buf;
    while (bytes_xfered < r_size) {
        r_count = r_size - bytes_xfered;
        err = pread(fd->fd_sys, p, r_count, offset + bytes_xfered);
        if (err == -1)
            goto ioerr;
        if (err == 0)
            break;
        bytes_xfered += err;
        p += err;
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->read_timing[2] += MPI_Wtime() - timing;
#endif

ioerr:
    if (err == -1)
        bytes_xfered = ncmpii_error_posix2nc("pread");

/*
if (offset > 0) {unsigned long long wkl[4];
    memcpy(wkl, buf, sizeof(unsigned long long) * 4);
    ncmpii_in_swapn(wkl, 4, 8);
    printf("%s at %d: %s pread offset=%lld r_size=%lld wkl=%llu %lld %lld %lld\n",__func__,__LINE__,fd->filename,offset,r_size,wkl[0],wkl[1],wkl[2],wkl[3]);
}
*/

    return bytes_xfered;
}

/*----< file_read() >--------------------------------------------------------*/
/* This is an independent call. */
static
MPI_Offset file_read(PNCIO_File *fd,
                     MPI_Offset  offset, /* relative to fileview */
                     void       *buf,
                     PNCIO_View  buf_view)
{
    MPI_Offset r_len=0;

// printf("%s at %d: offset=%lld buf_view size=%lld\n",__func__,__LINE__, offset,buf_view.size);

assert(fd->filetype == MPI_BYTE);
if (fd->flat_file.count > 0) assert(offset == 0); /* not whole file visible */

    if (buf_view.size == 0) /* zero-sized request */
        return NC_NOERR;

    if (buf_view.is_contig && fd->flat_file.is_contig) {
        if (fd->flat_file.count > 0) offset += fd->flat_file.off[0];
        r_len = PNCIO_ReadContig(fd, buf, buf_view.size, offset);
    }
    else
        r_len = PNCIO_GEN_ReadStrided(fd, buf, buf_view, offset);

    return r_len;
}

/*----< PNCIO_File_read_at() >------------------------------------------------*/
/* This is an independent call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
MPI_Offset PNCIO_File_read_at(PNCIO_File *fh,
                              MPI_Offset  offset,
                              void       *buf,
                              PNCIO_View  buf_view)
{
    assert(fh != NULL);

    if (buf_view.size == 0) return NC_NOERR;

    if (buf_view.size < 0) return NC_ENEGATIVECNT;

    /* PnetCDF has only 2 modes: read-only and read-write */
    // if (fh->access_mode & MPI_MODE_RDONLY) return NC_EPERM;

    return file_read(fh, offset, buf, buf_view);
}

/*----< PNCIO_File_read_at_all() >--------------------------------------------*/
/* This is a collective call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
MPI_Offset PNCIO_File_read_at_all(PNCIO_File *fh,
                                  MPI_Offset  offset,
                                  void       *buf,
                                  PNCIO_View  buf_view)
{
    int err=NC_NOERR;
    MPI_Offset r_len;

    assert(fh != NULL);

    if (buf_view.size < 0) err = NC_ENEGATIVECNT;

    /* PnetCDF has only 2 modes: read-only and read-write */
    // if (fh->access_mode & MPI_MODE_RDONLY && st == NC_NOERR) st = NC_EPERM;

    r_len = PNCIO_GEN_ReadStridedColl(fh, buf, buf_view, offset);

    return (err == NC_NOERR) ?  r_len :  err;
}

