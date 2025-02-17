/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <limits.h>
#include <stdarg.h> /* va_start(), va_end() */

#include <pncio.h>

/* some systems do not have pread/pwrite, or requrie XOPEN_SOURCE set higher
 * than we would like.  see #1973 */
#if (HAVE_DECL_PWRITE == 0)

#include <sys/types.h>
#include <unistd.h>

ssize_t pread(int fd, void *buf, size_t count, off_t offset);
ssize_t pwrite(int fd, const void *buf, size_t count, off_t offset);

ssize_t pread(int fd, void *buf, size_t count, off_t offset)
{
    off_t lseek_ret;
    off_t old_offset;
    ssize_t read_ret;

    old_offset = lseek(fd, 0, SEEK_CUR);
    lseek_ret = lseek(fd, offset, SEEK_SET);
    if (lseek_ret == -1)
        return lseek_ret;
    read_ret = read(fd, buf, count);
    if (read_ret < 0)
        return read_ret;
    /* man page says "file offset is not changed" */
    if ((lseek_ret = lseek(fd, old_offset, SEEK_SET)) < 0)
        return lseek_ret;

    return read_ret;
}

ssize_t pwrite(int fd, const void *buf, size_t count, off_t offset)
{
    off_t lseek_ret;
    off_t old_offset;
    ssize_t write_ret;

    old_offset = lseek(fd, 0, SEEK_CUR);
    lseek_ret = lseek(fd, offset, SEEK_SET);
    if (lseek_ret == -1)
        return lseek_ret;
    write_ret = write(fd, buf, count);
    if (write_ret < 0)
        return write_ret;
    /* man page says "file offset is not changed" */
    if ((lseek_ret = lseek(fd, old_offset, SEEK_SET)) < 0)
        return lseek_ret;

    return write_ret;
}
#endif

void PNCIO_Heap_merge(PNCIO_Access * others_req, MPI_Count * count,
                      MPI_Offset * srt_off, MPI_Count * srt_len, MPI_Count * start_pos,
                      int nprocs, int nprocs_recv, MPI_Count total_elements)
{
    typedef struct {
        MPI_Offset *off_list;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Offset *len_list;
#else
        int *len_list;
#endif
        MPI_Count nelem;
    } heap_struct;

    heap_struct *a, tmp;
    int i, j, heapsize, l, r, k, smallest;

    a = (heap_struct *) NCI_Malloc((nprocs_recv + 1) * sizeof(heap_struct));

    j = 0;
    for (i = 0; i < nprocs; i++)
        if (count[i]) {
            a[j].off_list = &(others_req[i].offsets[start_pos[i]]);
            a[j].len_list = &(others_req[i].lens[start_pos[i]]);
            a[j].nelem = count[i];
            j++;
        }

    /* build a heap out of the first element from each list, with
     * the smallest element of the heap at the root */

    heapsize = nprocs_recv;
    for (i = heapsize / 2 - 1; i >= 0; i--) {
        /* Heapify(a, i, heapsize); Algorithm from Cormen et al. pg. 143
         * modified for a heap with smallest element at root. I have
         * removed the recursion so that there are no function calls.
         * Function calls are too expensive. */
        k = i;
        for (;;) {
            l = 2 * (k + 1) - 1;
            r = 2 * (k + 1);

            if ((l < heapsize) && (*(a[l].off_list) < *(a[k].off_list)))
                smallest = l;
            else
                smallest = k;

            if ((r < heapsize) && (*(a[r].off_list) < *(a[smallest].off_list)))
                smallest = r;

            if (smallest != k) {
                tmp.off_list = a[k].off_list;
                tmp.len_list = a[k].len_list;
                tmp.nelem = a[k].nelem;

                a[k].off_list = a[smallest].off_list;
                a[k].len_list = a[smallest].len_list;
                a[k].nelem = a[smallest].nelem;

                a[smallest].off_list = tmp.off_list;
                a[smallest].len_list = tmp.len_list;
                a[smallest].nelem = tmp.nelem;

                k = smallest;
            } else
                break;
        }
    }

    for (i = 0; i < total_elements; i++) {
        /* extract smallest element from heap, i.e. the root */
        srt_off[i] = *(a[0].off_list);
        srt_len[i] = *(a[0].len_list);
        (a[0].nelem)--;

        if (!a[0].nelem) {
            a[0].off_list = a[heapsize - 1].off_list;
            a[0].len_list = a[heapsize - 1].len_list;
            a[0].nelem = a[heapsize - 1].nelem;
            heapsize--;
        } else {
            (a[0].off_list)++;
            (a[0].len_list)++;
        }

        /* Heapify(a, 0, heapsize); */
        k = 0;
        for (;;) {
            l = 2 * (k + 1) - 1;
            r = 2 * (k + 1);

            if ((l < heapsize) && (*(a[l].off_list) < *(a[k].off_list)))
                smallest = l;
            else
                smallest = k;

            if ((r < heapsize) && (*(a[r].off_list) < *(a[smallest].off_list)))
                smallest = r;

            if (smallest != k) {
                tmp.off_list = a[k].off_list;
                tmp.len_list = a[k].len_list;
                tmp.nelem = a[k].nelem;

                a[k].off_list = a[smallest].off_list;
                a[k].len_list = a[smallest].len_list;
                a[k].nelem = a[smallest].nelem;

                a[smallest].off_list = tmp.off_list;
                a[smallest].len_list = tmp.len_list;
                a[smallest].nelem = tmp.nelem;

                k = smallest;
            } else
                break;
        }
    }
    NCI_Free(a);
}

