/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

/* This file contains four functions:
 *
 * PNCIO_Calc_aggregator()
 * PNCIO_Calc_file_domains()
 * PNCIO_Calc_my_req()
 * PNCIO_Free_my_req()
 * PNCIO_Calc_others_req()
 * PNCIO_Free_others_req()
 *
 * The last three of these were originally in ad_read_coll.c, but they are
 * also shared with ad_write_coll.c.  I felt that they were better kept with
 * the rest of the shared aggregation code.
 */

/* Discussion of values available from above:
 *
 * MPI_Offset st_offsets[0..nprocs-1]
 * MPI_Offset end_offsets[0..nprocs-1]
 *    These contain a list of start and end offsets for each process in
 *    the communicator.  For example, an access at loc 10, size 10 would
 *    have a start offset of 10 and end offset of 19.
 * int nprocs
 *    number of processors in the collective I/O communicator
 * MPI_Offset min_st_offset
 * MPI_Offset fd_start[0..nprocs_for_coll-1]
 *    starting location of "file domain"; region that a given process will
 *    perform aggregation for (i.e. actually do I/O)
 * MPI_Offset fd_end[0..nprocs_for_coll-1]
 *    start + size - 1 roughly, but it can be less, or 0, in the case of
 *    uneven distributions
 */

/* PNCIO_Calc_aggregator()
 *
 * The intention here is to implement a function which provides basically
 * the same functionality as in Rajeev's original version of
 * PNCIO_Calc_my_req().  He used a ceiling division approach to assign the
 * file domains, and we use the same approach here when calculating the
 * location of an offset/len in a specific file domain.  Further we assume
 * this same distribution when calculating the rank_index, which is later
 *  used to map to a specific process rank in charge of the file domain.
 *
 * A better (i.e. more general) approach would be to use the list of file
 * domains only.  This would be slower in the case where the
 * original ceiling division was used, but it would allow for arbitrary
 * distributions of regions to aggregators.  We'd need to know the
 * nprocs_for_coll in that case though, which we don't have now.
 *
 * Note a significant difference between this function and Rajeev's old code:
 * this code doesn't necessarily return a rank in the range
 * 0..nprocs_for_coll; instead you get something in 0..nprocs.  This is a
 * result of the rank mapping; any set of ranks in the communicator could be
 * used now.
 *
 * Returns an integer representing a rank in the collective I/O communicator.
 *
 * The "len" parameter is also modified to indicate the amount of data
 * actually available in this file domain.
 */
int PNCIO_Calc_aggregator(PNCIO_File *fd,
                          MPI_Offset  off,
                          MPI_Offset  min_off,
                          MPI_Offset *len,
                          MPI_Offset  fd_size,
                          MPI_Offset *fd_end)
{
    int rank_index, rank;
    MPI_Offset avail_bytes;

    /* get an index into our array of aggregators */
    rank_index = (int) ((off - min_off + fd_size) / fd_size - 1);

    if (fd->hints->striping_unit > 0) {
        /* Implementation for file domain alignment. Note fd_end[] have been
         * aligned with file system lock boundaries when it was produced by
         * PNCIO_Calc_file_domains().
         */
        rank_index = 0;
        while (off > fd_end[rank_index])
            rank_index++;
    }

    /* we index into fd_end with rank_index, and fd_end was allocated to be no
     * bigger than fd->hins->cb_nodes.   If we ever violate that, we're
     * overrunning arrays.  Obviously, we should never ever hit this abort */
    if (rank_index >= fd->hints->cb_nodes || rank_index < 0) {
        fprintf(stderr,
                "Error in PNCIO_Calc_aggregator(): rank_index(%d) >= fd->hints->cb_nodes (%d) fd_size="OFFFMT" off="OFFFMT"\n",
                rank_index, fd->hints->cb_nodes, fd_size, off);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* remember here that even in Rajeev's original code it was the case that
     * different aggregators could end up with different amounts of data to
     * aggregate.  here we use fd_end[] to make sure that we know how much
     * data this aggregator is working with.
     *
     * the +1 is to take into account the end vs. length issue.
     */
    avail_bytes = fd_end[rank_index] + 1 - off;
    if (avail_bytes < *len) {
        /* this file domain only has part of the requested contig. region */
        *len = avail_bytes;
    }

    /* map our index to a rank */
    /* NOTE: FOR NOW WE DON'T HAVE A MAPPING...JUST DO 0..NPROCS_FOR_COLL */
    rank = fd->hints->ranklist[rank_index];

    return rank;
}

void PNCIO_Calc_file_domains(MPI_Offset  *st_offsets,
                             MPI_Offset  *end_offsets,
                             int          nprocs,
                             int          nprocs_for_coll,
                             MPI_Offset  *min_st_offset_ptr,
                             MPI_Offset **fd_start_ptr,
                             MPI_Offset **fd_end_ptr,
                             MPI_Offset  *fd_size_ptr,
                             int          striping_unit)
{
/* Divide the I/O workload among "nprocs_for_coll" processes. This is
   done by (logically) dividing the file into file domains (FDs); each
   process may directly access only its own file domain. */

    MPI_Offset min_st_offset, max_end_offset, *fd_start, *fd_end, fd_size;
    int i;

/* find min of start offsets and max of end offsets of all processes */

    min_st_offset = st_offsets[0];
    max_end_offset = end_offsets[0];

    for (i = 1; i < nprocs; i++) {
        min_st_offset = MIN(min_st_offset, st_offsets[i]);
        max_end_offset = MAX(max_end_offset, end_offsets[i]);
    }

/* determine the "file domain (FD)" of each process, i.e., the portion of
   the file that will be "owned" by each process */

/* partition the total file access range equally among nprocs_for_coll
   processes */
    fd_size = ((max_end_offset - min_st_offset + 1) + nprocs_for_coll - 1) / nprocs_for_coll;
    /* ceiling division as in HPF block distribution */

    *fd_start_ptr = (MPI_Offset *) NCI_Malloc(nprocs_for_coll * 2 * sizeof(MPI_Offset));
    *fd_end_ptr = *fd_start_ptr + nprocs_for_coll;

    fd_start = *fd_start_ptr;
    fd_end = *fd_end_ptr;

    /* Wei-keng Liao: implementation for fild domain alignment to nearest file
     * lock boundary (as specified by striping_unit hint).  Could also
     * experiment with other alignment strategies here */
    if (striping_unit > 0) {
        MPI_Offset end_off;
        int rem_front, rem_back;

        /* align fd_end[0] to the nearest file lock boundary */
        fd_start[0] = min_st_offset;
        end_off = fd_start[0] + fd_size;
        rem_front = end_off % striping_unit;
        rem_back = striping_unit - rem_front;
        if (rem_front < rem_back)
            end_off -= rem_front;
        else
            end_off += rem_back;
        fd_end[0] = end_off - 1;

        /* align fd_end[i] to the nearest file lock boundary */
        for (i = 1; i < nprocs_for_coll; i++) {
            fd_start[i] = fd_end[i - 1] + 1;
            end_off = min_st_offset + fd_size * (i + 1);
            rem_front = end_off % striping_unit;
            rem_back = striping_unit - rem_front;
            if (rem_front < rem_back)
                end_off -= rem_front;
            else
                end_off += rem_back;
            fd_end[i] = end_off - 1;
        }
        fd_end[nprocs_for_coll - 1] = max_end_offset;
    } else {    /* no hints set: do things the 'old' way */
        fd_start[0] = min_st_offset;
        fd_end[0] = min_st_offset + fd_size - 1;

        for (i = 1; i < nprocs_for_coll; i++) {
            fd_start[i] = fd_end[i - 1] + 1;
            fd_end[i] = fd_start[i] + fd_size - 1;
        }
    }

/* take care of cases in which the total file access range is not
   divisible by the number of processes. In such cases, the last
   process, or the last few processes, may have unequal load (even 0).
   For example, a range of 97 divided among 16 processes.
   Note that the division is ceiling division. */

    for (i = 0; i < nprocs_for_coll; i++) {
        if (fd_start[i] > max_end_offset)
            fd_start[i] = fd_end[i] = -1;
        if (fd_end[i] > max_end_offset)
            fd_end[i] = max_end_offset;
    }

    *fd_size_ptr = fd_size;
    *min_st_offset_ptr = min_st_offset;
}


/* PNCIO_Calc_my_req() - calculate what portions of the access requests
 * of this process are located in the file domains of various processes
 * (including this one)
 */
void PNCIO_Calc_my_req(PNCIO_File    *fd,
                       MPI_Offset     min_st_offset,
                       MPI_Offset    *fd_start,
                       MPI_Offset    *fd_end,
                       MPI_Offset     fd_size,
                       int            nprocs,
                       MPI_Count     *count_my_req_procs_ptr,
                       MPI_Count    **count_my_req_per_proc_ptr,
                       PNCIO_Access **my_req_ptr,
                       MPI_Aint     **buf_idx_ptr)
/* Possibly reconsider if buf_idx's are ok as int's, or should they be aints/offsets?
   They are used as memory buffer indices so it seems like the 2G limit is in effect */
{
    MPI_Count *count_my_req_per_proc, count_my_req_procs, l;
    MPI_Aint *buf_idx;
    int proc;
    size_t memLen, alloc_sz;
    MPI_Offset fd_len, rem_len, curr_idx, off, *off_ptr;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset *len_ptr;
#else
    int *len_ptr;
#endif
    PNCIO_Access *my_req;

    *count_my_req_per_proc_ptr = NCI_Calloc(nprocs, sizeof(MPI_Count));
    count_my_req_per_proc = *count_my_req_per_proc_ptr;
/* count_my_req_per_proc[i] gives the no. of contig. requests of this
   process in process i's file domain. calloc initializes to zero.
   I'm allocating memory of size nprocs, so that I can do an
   MPI_Alltoall later on.*/

    buf_idx = (MPI_Aint *) NCI_Malloc(nprocs * sizeof(MPI_Aint));
/* buf_idx is relevant only if buftype_is_contig.
   buf_idx[i] gives the index into user_buf where data received
   from proc. i should be placed. This allows receives to be done
   without extra buffer. This can't be done if buftype is not contig. */

    /* initialize buf_idx to -1 */
    for (int i = 0; i < nprocs; i++)
        buf_idx[i] = -1;

    /* one pass just to calculate how much space to allocate for my_req */
    for (MPI_Count i = 0; i < fd->flat_file.count; i++) {
        /* short circuit offset/len processing if len == 0
         *      (zero-byte  read/write */
        if (fd->flat_file.len[i] == 0)
            continue;
        off = fd->flat_file.off[i];
        fd_len = fd->flat_file.len[i];
        /* note: we set fd_len to be the total size of the access.  then
         * PNCIO_Calc_aggregator() will modify the value to return the
         * amount that was available from the file domain that holds the
         * first part of the access.
         */
        proc = PNCIO_Calc_aggregator(fd, off, min_st_offset, &fd_len, fd_size, fd_end);
        count_my_req_per_proc[proc]++;

        /* figure out how much data is remaining in the access (i.e. wasn't
         * part of the file domain that had the starting byte); we'll take
         * care of this data (if there is any) in the while loop below.
         */
        rem_len = fd->flat_file.len[i] - fd_len;

        while (rem_len != 0) {
            off += fd_len;      /* point to first remaining byte */
            fd_len = rem_len;   /* save remaining size, pass to calc */
            proc = PNCIO_Calc_aggregator(fd, off, min_st_offset, &fd_len,
                                         fd_size, fd_end);

            count_my_req_per_proc[proc]++;
            rem_len -= fd_len;  /* reduce remaining length by amount from fd */
        }
    }

/* now allocate space for my_req, offset, and len */

    *my_req_ptr = (PNCIO_Access *) NCI_Malloc(nprocs * sizeof(PNCIO_Access));
    my_req = *my_req_ptr;

    /* combine offsets and lens into a single regions so we can make one
     * exchange instead of two later on.  Over-allocate the 'offsets' array and
     * make 'lens' point to the over-allocated part
     */
    memLen = 0;
    for (int i = 0; i < nprocs; i++)
        memLen += count_my_req_per_proc[i];

#ifdef HAVE_MPI_LARGE_COUNT
    alloc_sz = sizeof(MPI_Offset) * 2;
    my_req[0].offsets = (MPI_Offset *) NCI_Malloc(memLen * alloc_sz);
    my_req[0].lens = my_req[0].offsets + memLen;
#else
    alloc_sz = sizeof(MPI_Offset) + sizeof(int);
    my_req[0].offsets = (MPI_Offset *) NCI_Malloc(memLen * alloc_sz);
    my_req[0].lens = (int*) (my_req[0].offsets + memLen);
#endif

    off_ptr = my_req[0].offsets;
    len_ptr = my_req[0].lens;
    count_my_req_procs = 0;
    for (int i = 0; i < nprocs; i++) {
        if (count_my_req_per_proc[i]) {
            my_req[i].offsets = off_ptr;
            off_ptr += count_my_req_per_proc[i];
            my_req[i].lens = len_ptr;
            len_ptr += count_my_req_per_proc[i];
            count_my_req_procs++;
        }
        my_req[i].count = 0;    /* will be incremented where needed
                                 * later */
    }

/* now fill in my_req */
    curr_idx = 0;
    for (MPI_Count i = 0; i < fd->flat_file.count; i++) {
        /* short circuit offset/len processing if len == 0
         *      (zero-byte  read/write */
        if (fd->flat_file.len[i] == 0)
            continue;
        off = fd->flat_file.off[i];
        fd_len = fd->flat_file.len[i];
        proc = PNCIO_Calc_aggregator(fd, off, min_st_offset, &fd_len, fd_size, fd_end);

        /* for each separate contiguous access from this process */
        if (buf_idx[proc] == -1) {
            assert(curr_idx == (MPI_Aint) curr_idx);
            buf_idx[proc] = (MPI_Aint) curr_idx;
        }

        l = my_req[proc].count;
        curr_idx += fd_len;

        rem_len = fd->flat_file.len[i] - fd_len;

        /* store the proc, offset, and len information in an array
         * of structures, my_req. Each structure contains the
         * offsets and lengths located in that process's FD,
         * and the associated count.
         */
        my_req[proc].offsets[l] = off;
        my_req[proc].lens[l] = fd_len;
        my_req[proc].count++;

        while (rem_len != 0) {
            off += fd_len;
            fd_len = rem_len;
            proc = PNCIO_Calc_aggregator(fd, off, min_st_offset, &fd_len,
                                         fd_size, fd_end);

            if (buf_idx[proc] == -1) {
                assert(curr_idx == (MPI_Aint) curr_idx);
                buf_idx[proc] = (MPI_Aint) curr_idx;
            }

            l = my_req[proc].count;
            curr_idx += fd_len;
            rem_len -= fd_len;

            my_req[proc].offsets[l] = off;
            my_req[proc].lens[l] = fd_len;
            my_req[proc].count++;
        }
    }

    *count_my_req_procs_ptr = count_my_req_procs;
    *buf_idx_ptr = buf_idx;
}

void PNCIO_Free_my_req(MPI_Count    *count_my_req_per_proc,
                       PNCIO_Access *my_req,
                       MPI_Aint     *buf_idx)
{
    NCI_Free(count_my_req_per_proc);
    NCI_Free(my_req[0].offsets);
    NCI_Free(my_req);
    NCI_Free(buf_idx);
}

void PNCIO_Calc_others_req(PNCIO_File    *fd,
                           MPI_Count      count_my_req_procs,
                           MPI_Count     *count_my_req_per_proc,
                           PNCIO_Access  *my_req,
                           int            nprocs,
                           int            myrank,
                           MPI_Count     *count_others_req_procs_ptr,
                           MPI_Count    **count_others_req_per_proc_ptr,
                           PNCIO_Access **others_req_ptr)
{
/* determine what requests of other processes lie in this process's
   file domain */

/* count_others_req_procs = number of processes whose requests lie in
   this process's file domain (including this process itself)
   count_others_req_per_proc[i] indicates how many separate contiguous
   requests of proc. i lie in this process's file domain. */

    MPI_Count *count_others_req_per_proc, count_others_req_procs;
    size_t alloc_sz;
    int i, j;
    MPI_Request *requests;
    PNCIO_Access *others_req;
    size_t memLen;
    MPI_Offset *off_ptr;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset *len_ptr;
    MPI_Count *mem_ptr;
#else
    int *len_ptr;
    MPI_Aint *mem_ptr;
#endif

/* first find out how much to send/recv and from/to whom */
    count_others_req_per_proc = NCI_Malloc(nprocs * sizeof(MPI_Count));

    MPI_Alltoall(count_my_req_per_proc, 1, MPI_COUNT,
                 count_others_req_per_proc, 1, MPI_COUNT, fd->comm);

    *others_req_ptr = (PNCIO_Access *) NCI_Malloc(nprocs * sizeof(PNCIO_Access));
    others_req = *others_req_ptr;

    memLen = 0;
    for (i = 0; i < nprocs; i++)
        memLen += count_others_req_per_proc[i];

#ifdef HAVE_MPI_LARGE_COUNT
    alloc_sz = sizeof(MPI_Offset) * 2 + sizeof(MPI_Count);
    others_req[0].offsets = (MPI_Offset *) NCI_Malloc(memLen * alloc_sz);
    others_req[0].lens = others_req[0].offsets + memLen;
    others_req[0].mem_ptrs = (MPI_Count*) (others_req[0].lens + memLen);
#else
    alloc_sz = sizeof(MPI_Offset) + sizeof(int) + sizeof(MPI_Aint);
    others_req[0].offsets = (MPI_Offset *) NCI_Malloc(memLen * alloc_sz);
    others_req[0].lens = (int *) (others_req[0].offsets + memLen);
    others_req[0].mem_ptrs = (MPI_Aint*) (others_req[0].lens + memLen);
#endif
    off_ptr = others_req[0].offsets;
    len_ptr = others_req[0].lens;
    mem_ptr = others_req[0].mem_ptrs;

    count_others_req_procs = 0;
    for (i = 0; i < nprocs; i++) {
        if (count_others_req_per_proc[i]) {
            others_req[i].count = count_others_req_per_proc[i];
            others_req[i].offsets = off_ptr;
            off_ptr += count_others_req_per_proc[i];
            others_req[i].lens = len_ptr;
            len_ptr += count_others_req_per_proc[i];
            others_req[i].mem_ptrs = mem_ptr;
            mem_ptr += count_others_req_per_proc[i];
            count_others_req_procs++;
        } else
            others_req[i].count = 0;
    }
    *count_others_req_per_proc_ptr = count_others_req_per_proc;

/* now send the calculated offsets and lengths to respective processes */

    requests = (MPI_Request *)
        NCI_Malloc((count_my_req_procs + count_others_req_procs) * 2 * sizeof(MPI_Request));

    j = 0;
    for (i = 0; i < nprocs; i++) {
        if (others_req[i].count == 0)
            continue;
        if (i == myrank) {
            /* send to self uses memcpy()C, here others_req[i].count == my_req[i].count */
            memcpy(others_req[i].offsets, my_req[i].offsets,
                   my_req[i].count * sizeof(MPI_Offset));
#ifdef HAVE_MPI_LARGE_COUNT
            memcpy(others_req[i].lens, my_req[i].lens,
                   my_req[i].count * sizeof(MPI_Offset));
#else
            memcpy(others_req[i].lens, my_req[i].lens,
                   my_req[i].count * sizeof(int));
#endif
        }
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Irecv_c(others_req[i].offsets, others_req[i].count,
                        MPI_OFFSET, i, i + myrank, fd->comm, &requests[j++]);
            MPI_Irecv_c(others_req[i].lens, others_req[i].count,
                        MPI_OFFSET, i, i + myrank, fd->comm, &requests[j++]);
#else
            assert(others_req[i].count <= 2147483647); /* overflow 4-byte int */
            MPI_Irecv(others_req[i].offsets, (int)others_req[i].count,
                      MPI_OFFSET, i, i + myrank, fd->comm, &requests[j++]);
            MPI_Irecv(others_req[i].lens, (int)others_req[i].count,
                      MPI_INT, i, i + myrank, fd->comm, &requests[j++]);
#endif
        }
    }

    for (i = 0; i < nprocs; i++) {
        if (my_req[i].count && i != myrank) {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Isend_c(my_req[i].offsets, my_req[i].count,
                        MPI_OFFSET, i, i + myrank, fd->comm, &requests[j++]);
            MPI_Isend_c(my_req[i].lens, my_req[i].count,
                        MPI_OFFSET, i, i + myrank, fd->comm, &requests[j++]);
#else
            assert(my_req[i].count <= 2147483647); /* overflow 4-byte int */
            MPI_Isend(my_req[i].offsets, (int)my_req[i].count,
                      MPI_OFFSET, i, i + myrank, fd->comm, &requests[j++]);
            MPI_Isend(my_req[i].lens, (int)my_req[i].count,
                      MPI_INT, i, i + myrank, fd->comm, &requests[j++]);
#endif
        }
    }

    if (j) {
#ifdef HAVE_MPI_STATUSES_IGNORE
        MPI_Waitall(j, requests, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *) NCI_Malloc(j * sizeof(MPI_Status));
        MPI_Waitall(j, requests, statuses);
        NCI_Free(statuses);
#endif
    }

    NCI_Free(requests);

    *count_others_req_procs_ptr = count_others_req_procs;
}

void PNCIO_Free_others_req(MPI_Count    *count_others_req_per_proc,
                           PNCIO_Access *others_req)
{
    NCI_Free(count_others_req_per_proc);
    NCI_Free(others_req[0].offsets);
    NCI_Free(others_req);
}

