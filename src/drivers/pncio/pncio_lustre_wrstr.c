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
            if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)  \
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
            if (w_len < 0) {                                                  \
                NCI_Free(writebuf);                                           \
                return w_len;                                                 \
            }                                                                 \
            total_w_len += w_len;                                             \
            writebuf_off = req_off;                                           \
        }                                                                     \
        writebuf_off = req_off;                                               \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) * stripe_size     \
                           - writebuf_off);                                   \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
        r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len, writebuf_off);   \
        if (r_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return r_len;                                                     \
        }                                                                     \
    }                                                                         \
    write_sz = (MIN(req_len, writebuf_off + writebuf_len - req_off));         \
    memcpy(writebuf + req_off - writebuf_off, (char *)buf + userbuf_off,      \
           write_sz);                                                         \
    while (write_sz != req_len) {                                             \
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);  \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);           \
        if (w_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return w_len;                                                     \
        }                                                                     \
        total_w_len += w_len;                                                 \
        req_len -= write_sz;                                                  \
        userbuf_off += write_sz;                                              \
        writebuf_off += writebuf_len;                                         \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) * stripe_size     \
                           - writebuf_off);                                   \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
        r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len, writebuf_off);   \
        if (r_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return r_len;                                                     \
        }                                                                     \
        write_sz = MIN(req_len, writebuf_len);                                \
        memcpy(writebuf, (char *)buf + userbuf_off, write_sz);                \
    }                                                                         \
}

/* this macro is used when filetype is contig and buftype is not contig.
 * it does not do a read-modify-write and does not lock
 */
#define BUFFERED_WRITE_WITHOUT_READ {                                         \
    if (req_off >= writebuf_off + writebuf_len) {                             \
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);  \
        if (w_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return w_len;                                                     \
        }                                                                     \
        total_w_len += w_len;                                                 \
        writebuf_off = req_off;                                               \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) * stripe_size     \
                           - writebuf_off);                                   \
    }                                                                         \
    write_sz = MIN(req_len, writebuf_off + writebuf_len - req_off);           \
    memcpy(writebuf + req_off - writebuf_off,                                 \
           (char *)buf + userbuf_off, write_sz);                              \
    while (write_sz != req_len) {                                             \
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);  \
        if (w_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return w_len;                                                     \
        }                                                                     \
        total_w_len += w_len;                                                 \
        req_len -= write_sz;                                                  \
        userbuf_off += write_sz;                                              \
        writebuf_off += writebuf_len;                                         \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) * stripe_size     \
                           - writebuf_off);                                   \
        write_sz = MIN(req_len, writebuf_len);                                \
        memcpy(writebuf, (char *)buf + userbuf_off, write_sz);                \
    }                                                                         \
}

MPI_Offset PNCIO_LUSTRE_WriteStrided(PNCIO_File *fd,
                                     const void *buf,
                                     PNCIO_View  buf_view,
                                     MPI_Offset  offset)
{
    char *writebuf;
    int i, j, k, st_index=0, stripe_size;
    /* offset is in units of etype relative to the filetype. */
    MPI_Offset i_offset, sum, num, size, abs_off_in_filetype=0, off, disp;
    MPI_Offset userbuf_off, req_off, end_offset=0, writebuf_off, start_off;
    MPI_Offset new_bwr_size, new_fwr_size, st_fwr_size, fwr_size=0, bwr_size;
    MPI_Offset req_len, r_len, w_len, total_w_len=0;
    MPI_Count bufsize, writebuf_len, write_sz;

    /* The case of both buftype and filetype being contiguous has gone to
     * PNCIO_WriteContig().
     */

// printf("%s at %d:\n",__func__,__LINE__);

    if (fd->hints->romio_ds_write == PNCIO_HINT_DISABLE) {
        /* if user has disabled data sieving on writes, use naive
         * approach instead.
         */
        return PNCIO_GEN_WriteStrided_naive(fd, buf, buf_view, offset);
    }


/* PnetCDF always sets these 3 conditions */
assert(fd->filetype == MPI_BYTE);
assert(fd->flat_file.size == buf_view.size);
if (fd->flat_file.count > 0) assert(offset == 0); /* not whole file visible */

    bufsize = buf_view.size;

    /* get striping info */
    stripe_size = fd->hints->striping_unit;

    if (!buf_view.is_contig && fd->flat_file.is_contig) {
        /* noncontiguous in write buffer, contiguous in file. */

        off = fd->disp + offset;
        if (fd->flat_file.count > 0) off += fd->flat_file.off[0];

        start_off = off;
        end_offset = start_off + bufsize - 1;

        /* write stripe size buffer each time */
        writebuf = (char *) NCI_Malloc(MIN(bufsize, stripe_size));
        writebuf_off = 0;
        writebuf_len = 0;

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed
         */
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, bufsize);

        for (i = 0; i < buf_view.count; i++) {
            userbuf_off = buf_view.off[i];
            req_off = off;
            req_len = buf_view.len[i];
            BUFFERED_WRITE_WITHOUT_READ;
            off += buf_view.len[i];
        }

        /* write the buffer out the last round */
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);

        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, bufsize);

        NCI_Free(writebuf);

        if (w_len < 0) return w_len;
        total_w_len += w_len;

    } else { /* contiguous buffer and non-contiguous in file */
        disp = fd->disp;
/* for non-contiguous in file, PnetCDF always uses disp == 0 */
assert(disp == 0);

        /* find the starting index in fd->flat_file offset-length pairs */
        sum = 0;
        for (i = 0; i < fd->flat_file.count; i++) {
            sum += fd->flat_file.len[i];
            if (sum > offset) {
                st_index = i;
                fwr_size = sum - offset;
                abs_off_in_filetype = fd->flat_file.off[i] +
                    offset - (sum - fd->flat_file.len[i]);
                break;
            }
        }

        /* abs. offset in bytes in the file */
        offset = disp + abs_off_in_filetype;

        start_off = offset;

        /* Write request is within single flat_file contig block. This could
         * happen, for example, with subarray types that are actually fairly
         * contiguous.
         */
        if (buf_view.is_contig && bufsize <= fwr_size) {
            req_off = start_off;
            req_len = bufsize;
            end_offset = start_off + bufsize - 1;
            writebuf = (char *) NCI_Malloc(MIN(bufsize, stripe_size));
            memset(writebuf, -1, (size_t)MIN(bufsize, stripe_size));
            writebuf_off = 0;
            writebuf_len = 0;
            userbuf_off = 0;
            BUFFERED_WRITE_WITHOUT_READ;

            /* write the buffer out the last round */
            if (fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);

            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (w_len > 0) total_w_len += w_len;

            if (fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);

            NCI_Free(writebuf);

            return total_w_len;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be write, end_offset=99 */

        st_fwr_size = fwr_size;
        j = st_index;
        i_offset = fwr_size = MIN(st_fwr_size, bufsize);
        end_offset = offset + fwr_size - 1;
        while (i_offset < bufsize) {
            j++;
assert(j < fd->flat_file.count);
            off = disp + fd->flat_file.off[j];
            fwr_size = MIN(fd->flat_file.len[j], bufsize - i_offset);
            i_offset += fwr_size;
            end_offset = off + fwr_size - 1;
        }

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        writebuf_off = 0;
        writebuf_len = 0;
        writebuf = (char *) NCI_Malloc(stripe_size);
        memset(writebuf, -1, stripe_size);

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
                    off += fwr_size;
                    /* no more I/O needed. off is incremented by fwr_size. */
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
                    /* reached end of contiguous block in file */
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

        /* write the buffer out the last round */
        if (writebuf_len) {
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);
            if (w_len < 0) return w_len;
            total_w_len += w_len;
        }
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        NCI_Free(writebuf);
    }

    return buf_view.size;
}
