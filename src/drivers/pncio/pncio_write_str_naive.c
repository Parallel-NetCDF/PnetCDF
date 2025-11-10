/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

MPI_Offset PNCIO_GEN_WriteStrided_naive(PNCIO_File      *fd,
                                        const void      *buf,
                                        PNCIO_View  buf_view,
                                        MPI_Offset       offset)
{
    int b_index;
    MPI_Count bufsize;

    /* bwr == buffer write; fwr == file write */
    MPI_Offset bwr_size, fwr_size = 0, sum, size_in_filetype, size;
    MPI_Offset abs_off_in_filetype = 0, req_len, userbuf_off;
    MPI_Offset off, req_off, disp, end_offset = 0, start_off;
    MPI_Offset w_len, total_w_len=0;

/* PnetCDF always sets fd->filetype == MPI_BYTE */
assert(fd->filetype == MPI_BYTE);

    /* Contiguous both in buftype and filetype should have been handled in a
     * call to PNCIO_WriteContig() earlier.
     */
    assert(!(buf_view.is_contig && fd->flat_file.is_contig));

    bufsize = buf_view.size;

    if (!buf_view.is_contig && fd->flat_file.is_contig) {
        /* noncontiguous in memory, contiguous in file. */

        off = fd->disp + offset;

        start_off = off;
        end_offset = off + bufsize - 1;

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        /* for each region in the buffer, grab the data and put it in place */
        for (b_index = 0; b_index < buf_view.count; b_index++) {
            userbuf_off = buf_view.off[b_index];
            req_off = off;
            req_len = buf_view.len[b_index];

            w_len = PNCIO_WriteContig(fd, (char *) buf + userbuf_off,
                                      req_len, req_off);
            if (w_len < 0) return w_len;
            total_w_len += w_len;

            /* off is (potentially) used to save the final offset later */
            off += buf_view.len[b_index];
        }

        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);
    }
    else {      /* noncontiguous in file */
        int f_index, st_index = 0;
        MPI_Offset st_fwr_size;

        /* First we're going to calculate a set of values for use in all
         * the noncontiguous in file cases:
         * start_off - starting byte position of data in file
         * end_offset - last byte offset to be accessed in the file
         * st_index - index of block in first filetype that we will be
         *            starting in (?)
         * st_fwr_size - size of the data in the first filetype block
         *               that we will write (accounts for being part-way
         *               into writing this block of the filetype
         */

        disp = fd->disp;

        size_in_filetype = offset;

        sum = 0;
        for (f_index = 0; f_index < fd->flat_file.count; f_index++) {
            sum += fd->flat_file.len[f_index];
            if (sum > size_in_filetype) {
                st_index = f_index;
                fwr_size = sum - size_in_filetype;
                abs_off_in_filetype = fd->flat_file.off[f_index] +
                    size_in_filetype - (sum - fd->flat_file.len[f_index]);
                break;
            }
        }

        /* abs. offset in bytes in the file */
        start_off = disp + abs_off_in_filetype;

        st_fwr_size = fwr_size;

        /* start_off, st_index, and st_fwr_size are
         * all calculated at this point
         */

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_off=0 and 100 bytes to be written, end_offset=99
         */
        f_index = st_index;
        fwr_size = MIN(st_fwr_size, bufsize);
        userbuf_off = fwr_size;
        end_offset = start_off + fwr_size - 1;
        while (userbuf_off < bufsize) {
            f_index++;
            fwr_size = MIN(fd->flat_file.len[f_index],
                               bufsize - userbuf_off);
            userbuf_off += fwr_size;
            end_offset = disp + fd->flat_file.off[f_index] + fwr_size - 1;
        }

        /* End of calculations.  At this point the following values have
         * been calculated and are ready for use:
         * - start_off
         * - end_offset
         * - st_index
         * - st_fwr_size
         */

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        if (buf_view.is_contig && !fd->flat_file.is_contig) {
            /* contiguous in memory, noncontiguous in file. should be the
             * most common case.
             */

            userbuf_off = 0;
            f_index = st_index;
            off = start_off;
            fwr_size = MIN(st_fwr_size, bufsize);

            /* while there is still space in the buffer, write more data */
            while (userbuf_off < bufsize) {
                if (fwr_size) {
                    /* TYPE_UB and TYPE_LB can result in
                     * fwr_size = 0. save system call in such cases */
                    req_off = off;
                    req_len = fwr_size;

                    w_len = PNCIO_WriteContig(fd, (char *) buf + userbuf_off,
                                              req_len, req_off);
                    if (w_len < 0) return w_len;
                    total_w_len += w_len;
                }
                userbuf_off += fwr_size;
                if (userbuf_off >= bufsize) break;

                if (off + fwr_size < disp + fd->flat_file.off[f_index] +
                    fd->flat_file.len[f_index]) {
                    /* important that this value be correct, as it is
                     * used to set the offset in the fd near the end of
                     * this function.
                     */
                    off += fwr_size;
                }
                /* did not reach end of contiguous block in filetype.
                 * no more I/O needed. off is incremented by fwr_size.
                 */
                else {
                    f_index++;
assert(f_index < fd->flat_file.count);
                    off = disp + fd->flat_file.off[f_index];
                    fwr_size = MIN(fd->flat_file.len[f_index],
                                       bufsize - userbuf_off);
                }
            }
        } else {
            MPI_Offset i_offset, tmp_bufsize = 0;
            /* noncontiguous in memory as well as in file */

            b_index = 0;
            i_offset = buf_view.off[0];
            f_index = st_index;
            off = start_off;
            fwr_size = st_fwr_size;
            bwr_size = buf_view.len[0];

            /* while we haven't read size * count bytes, keep going */
            while (tmp_bufsize < bufsize) {
                MPI_Offset new_bwr_size = bwr_size, new_fwr_size = fwr_size;

                size = MIN(fwr_size, bwr_size);
                /* keep max of a single read amount <= INT_MAX */
                size = MIN(size, INT_MAX);

                if (size) {
                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;

                    w_len = PNCIO_WriteContig(fd, (char *) buf + userbuf_off,
                                              req_len, req_off);
                    if (w_len < 0) return w_len;
                    total_w_len += w_len;
                }

                tmp_bufsize += size;
                if (tmp_bufsize >= bufsize) break;

                if (size == fwr_size) {
                    f_index++;
assert(f_index < fd->flat_file.count);
                    off = disp + fd->flat_file.off[f_index];
                    new_fwr_size = fd->flat_file.len[f_index];
                    if (size != bwr_size) {
                        i_offset += size;
                        new_bwr_size -= size;
                    }
                }

                if (size == bwr_size) {
                    /* reached end of contiguous block in memory */
                    b_index++;
assert(b_index < buf_view.count);
                    i_offset = buf_view.off[b_index];
                    new_bwr_size = buf_view.len[b_index];
                    if (size != fwr_size) {
                        off += size;
                        new_fwr_size -= size;
                    }
                }
                fwr_size = new_fwr_size;
                bwr_size = new_bwr_size;
            }
        }

        /* unlock the file region if we locked it */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

    }   /* end of (else noncontiguous in file) */

    return total_w_len;
}
