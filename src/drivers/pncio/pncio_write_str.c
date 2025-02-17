/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

#define BUFFERED_WRITE {                                                      \
    if (req_off >= writebuf_off + writebuf_len) {                             \
        if (writebuf_len) {                                                   \
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,             \
                                      writebuf_off);                          \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE)  \
                    PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);   \
            if (w_len < 0) goto fn_exit;                                      \
            total_w_len += w_len;                                             \
        }                                                                     \
        writebuf_off = req_off;                                               \
        writebuf_len = MIN(max_bufsize,end_offset-writebuf_off+1);            \
        if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
        r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len, writebuf_off);   \
        if (r_len < 0) goto fn_exit;                                          \
    }                                                                         \
    write_sz = (MPI_Aint)MIN(req_len, writebuf_off+writebuf_len-req_off);     \
    memcpy(writebuf+req_off-writebuf_off, (char*)buf +userbuf_off, write_sz); \
    while (write_sz != req_len) {                                             \
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);  \
        if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);           \
        if (w_len < 0) goto fn_exit;                                          \
        total_w_len += w_len;                                                 \
        req_len -= write_sz;                                                  \
        userbuf_off += write_sz;                                              \
        writebuf_off += writebuf_len;                                         \
        writebuf_len = MIN(max_bufsize,end_offset-writebuf_off+1);            \
        if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
        r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len, writebuf_off);   \
        if (r_len < 0) goto fn_exit;                                          \
        write_sz = MIN(req_len, writebuf_len);                                \
        memcpy(writebuf, (char *)buf + userbuf_off, write_sz);                \
    }                                                                         \
}


MPI_Offset PNCIO_GEN_WriteStrided(PNCIO_File *fd,
                                  const void *buf,
                                  PNCIO_View  buf_view,
                                  MPI_Offset  offset)
{

/* offset is in units of etype relative to the filetype. */

    char *writebuf = NULL;
    int i, j, k, st_index = 0;
    MPI_Aint writebuf_len, max_bufsize, write_sz, bufsize;
    MPI_Offset i_offset, sum, num, size, abs_off_in_filetype=0;
    MPI_Offset userbuf_off, off, req_off, disp, end_offset=0;
    MPI_Offset writebuf_off, start_off, new_bwr_size, new_fwr_size;
    MPI_Offset st_fwr_size, fwr_size = 0, bwr_size, req_len;
    MPI_Offset r_len, w_len, total_w_len=0;

    /* Contiguous both in buftype and filetype should have been handled in a
     * call to PNCIO_WriteContig() earlier.
     */
    assert(!(buf_view.is_contig && fd->flat_file.is_contig));

    if (fd->hints->ds_write == PNCIO_HINT_DISABLE) {
        /* If user has disabled data sieving on reads, use naive approach
         * instead.
         */
        return PNCIO_GEN_WriteStrided_naive(fd, buf, buf_view, offset);
    }

// printf("%s at %d: offset=%lld\n",__func__,__LINE__, offset);

/* PnetCDF always set these 3 conditions */
assert(fd->filetype == MPI_BYTE);
assert(fd->flat_file.size == buf_view.size);
if (fd->flat_file.count > 0) assert(offset == 0); /* not whole file visible */

    bufsize = buf_view.size;

    /* get max_bufsize from the info object. */
    max_bufsize = fd->hints->ind_wr_buffer_size;

    if (!buf_view.is_contig && fd->flat_file.is_contig) {
        /* noncontiguous in memory, contiguous in file. */

        off = fd->disp + offset;
assert(fd->disp == 0);
        if (fd->flat_file.count > 0) off += fd->flat_file.off[0];

        start_off = off;
        end_offset = off + bufsize - 1;
        writebuf_off = off;
        writebuf = (char *) NCI_Malloc(max_bufsize);
        writebuf_len = MIN(max_bufsize, end_offset - writebuf_off + 1);

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        for (i = 0; i < buf_view.count; i++) {
            userbuf_off = buf_view.off[i];
            req_off = off;
            req_len = buf_view.len[i];

            /* BUFFERED_WRITE_WITHOUT_READ does neither read-modify-write nor
             * file lock
             */
            if (req_off >= writebuf_off + writebuf_len) {
                w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,
                                          writebuf_off);
                if (w_len < 0) goto fn_exit;
                total_w_len += w_len;
                writebuf_off = req_off;
                writebuf_len = MIN(max_bufsize,end_offset-writebuf_off+1);
            }
            write_sz = MIN(req_len, writebuf_off + writebuf_len - req_off);
            memcpy(writebuf+req_off-writebuf_off, (char*)buf +userbuf_off,
                   write_sz);
            while (write_sz != req_len) {
                w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,
                                          writebuf_off);
                if (w_len < 0) goto fn_exit;
                total_w_len += w_len;
                req_len -= write_sz;
                userbuf_off += write_sz;
                writebuf_off += writebuf_len;
                writebuf_len = MIN(max_bufsize,end_offset-writebuf_off+1);
                write_sz = MIN(req_len, writebuf_len);
                memcpy(writebuf, (char *)buf + userbuf_off, write_sz);
            }

            off += buf_view.len[i];
        }

        /* write the buffer out finally */
        if (writebuf_len) {
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (w_len >= 0) total_w_len += w_len;
        }
        else
            w_len = 0;

        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        if (w_len < 0)
            goto fn_exit;
    }
    else { /* noncontiguous in file */
        MPI_Offset size_in_filetype = offset;

        disp = fd->disp;
assert(fd->disp == 0);

        sum = 0;
        for (i = 0; i < fd->flat_file.count; i++) {
            sum += fd->flat_file.len[i];
            if (sum > size_in_filetype) {
                st_index = i;
                fwr_size = sum - size_in_filetype;
                abs_off_in_filetype = fd->flat_file.off[i] +
                    size_in_filetype - (sum - fd->flat_file.len[i]);
                break;
            }
        }

        /* abs. offset in bytes in the file */
        offset = disp + abs_off_in_filetype;

        start_off = offset;
assert(offset == abs_off_in_filetype);

// printf("%s at %d: start_off=%lld abs_off_in_filetype=%lld\n",__func__,__LINE__,start_off,abs_off_in_filetype);

        /* Write request is within single flat_file contig block. This could
         * happen, for example, with subarray types that are actually fairly
         * contiguous.
         */
        if (buf_view.is_contig && bufsize <= fwr_size) {
            /* though MPI api has an integer 'count' parameter, derived
             * datatypes might describe more bytes than can fit into an integer.
             * if we've made it this far, we can pass a count of original
             * datatypes, instead of a count of bytes (which might overflow)
             * Other WriteContig calls in this path are operating on data
             * sieving buffer */
            PNCIO_WRITE_LOCK(fd, offset, SEEK_SET, bufsize);
            w_len = PNCIO_WriteContig(fd, buf, buf_view.size, offset);
            if (w_len > 0) total_w_len += w_len;
            PNCIO_UNLOCK(fd, offset, SEEK_SET, bufsize);

            goto fn_exit;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be write, end_offset=99 */

        st_fwr_size = fwr_size;
        j = st_index;
        fwr_size = MIN(fwr_size, bufsize);
        i_offset = fwr_size;
        end_offset = offset + fwr_size - 1;
        while (i_offset < bufsize) {
            j++;
            fwr_size = MIN(fd->flat_file.len[j], bufsize - i_offset);
            i_offset += fwr_size;
            end_offset = disp + fd->flat_file.off[j] + fwr_size - 1;
        }

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        writebuf_off = 0;
        writebuf_len = 0;
        writebuf = (char *) NCI_Malloc(max_bufsize);
        memset(writebuf, -1, max_bufsize);

        if (buf_view.is_contig && !fd->flat_file.is_contig) {
            /* contiguous in memory, noncontiguous in file should be the most
             * common case.
             */
            i_offset = 0;
            j = st_index;
            off = offset;
            fwr_size = MIN(st_fwr_size, bufsize);
            while (i_offset < bufsize) {
                if (fwr_size) {
                    req_off = off;
                    req_len = fwr_size;
                    userbuf_off = i_offset;
                    BUFFERED_WRITE;
                }

                i_offset += fwr_size;
                if (i_offset >= bufsize) break;

                if (off + fwr_size < disp + fd->flat_file.off[j] +
                                            fd->flat_file.len[j])
                    off += fwr_size; /* off is incremented by fwr_size. */
                else {
                    j++;
assert(j < fd->flat_file.count);
                    off = disp + fd->flat_file.off[j];
                    fwr_size = MIN(fd->flat_file.len[j],
                                       bufsize - i_offset);
                }
            }
        } else {
            /* noncontiguous in memory as well as in file */
            k = num = 0;
            i_offset = buf_view.off[0];
            j = st_index;
            off = offset;
            fwr_size = st_fwr_size;
            bwr_size = buf_view.len[0];

            while (num < bufsize) {
                size = MIN(fwr_size, bwr_size);
                if (size) {
                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;
                    BUFFERED_WRITE;
                }

                num += size;
                if (num >= bufsize) break;

                new_fwr_size = fwr_size;
                new_bwr_size = bwr_size;

                if (size == fwr_size) {
                    j++;
assert(j < fd->flat_file.count);
                    off = disp + fd->flat_file.off[j];
                    new_fwr_size = fd->flat_file.len[j];
                    if (size != bwr_size) {
                        i_offset += size;
                        new_bwr_size -= size;
                    }
                }

                if (size == bwr_size) {
                    /* reached end of contiguous block in memory */

                    k++;
assert(k < buf_view.count);
                    i_offset = buf_view.off[k];
                    new_bwr_size = buf_view.len[k];
                    if (size != fwr_size) {
                        off += size;
                        new_fwr_size -= size;
                    }
                }
                fwr_size = new_fwr_size;
                bwr_size = new_bwr_size;
            }
        }

        /* write the buffer out finally */
        if (writebuf_len) {
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);
            if (w_len < 0) goto fn_exit;
            total_w_len += w_len;
        }
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
    }

fn_exit:
    if (writebuf != NULL)
        NCI_Free(writebuf);

    return total_w_len;
}
