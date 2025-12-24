/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdbool.h> /* type bool */

#include <pncio.h>

/* prototypes of functions used for collective reads only. */
static
MPI_Offset Read_and_exch(PNCIO_File *fd, void *buf,
                          PNCIO_View buf_view, int nprocs,
                          int myrank, PNCIO_Access *others_req,
                          MPI_Offset
                          min_st_offset, MPI_Offset fd_size,
                          MPI_Offset * fd_start, MPI_Offset * fd_end,
                          MPI_Aint * buf_idx);

static void R_Exchange_data(PNCIO_File *fd, void *buf,
                            PNCIO_View buf_view,
                            MPI_Count * send_size, MPI_Count * recv_size,
                            MPI_Count * count, MPI_Count * start_pos,
                            MPI_Count * partial_send,
                            MPI_Count * recd_from_proc, int nprocs,
                            int myrank,
                            MPI_Offset min_st_offset,
                            MPI_Offset fd_size,
                            MPI_Offset * fd_start, MPI_Offset * fd_end,
                            PNCIO_Access * others_req,
                            int iter, MPI_Aint * buf_idx,
                            MPI_Aint * actual_recved_bytes);
static void Fill_user_buffer(PNCIO_File *fd, void *buf,
                             PNCIO_View buf_view, char **recv_buf,
                             MPI_Count * recv_size,
                             MPI_Count * recd_from_proc, int nprocs,
                             MPI_Offset min_st_offset,
                             MPI_Offset fd_size, MPI_Offset * fd_start,
                             MPI_Offset * fd_end);

MPI_Offset PNCIO_GEN_ReadStridedColl(PNCIO_File *fd,
                                     void       *buf,
                                     PNCIO_View  buf_view,
                                     MPI_Offset  offset)
{
/* Uses a generalized version of the extended two-phase method described
   in "An Extended Two-Phase Method for Accessing Sections of
   Out-of-Core Arrays", Rajeev Thakur and Alok Choudhary,
   Scientific Programming, (5)4:301--317, Winter 1996.
   http://www.mcs.anl.gov/home/thakur/ext2ph.ps */

    PNCIO_Access *my_req;
    /* array of nprocs structures, one for each other process in
     * whose file domain this process's request lies */

    PNCIO_Access *others_req;
    /* array of nprocs structures, one for each other process
     * whose request lies in this process's file domain. */

    int nprocs, nprocs_for_coll, myrank;
    int interleave_count = 0;
    MPI_Count *count_my_req_per_proc, count_my_req_procs;
    MPI_Count *count_others_req_per_proc, count_others_req_procs;
    MPI_Offset start_offset, end_offset, fd_size, min_st_offset;
    MPI_Offset *st_offsets = NULL, *fd_start = NULL,
        *fd_end = NULL, *end_offsets = NULL;
    MPI_Aint *buf_idx = NULL;
    MPI_Offset r_len, total_r_len=0;

// printf("%s at %d:\n",__func__,__LINE__);

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
double curT = MPI_Wtime();
#endif

    /* number of aggregators, cb_nodes, is stored in the hints */
    nprocs_for_coll = fd->hints->cb_nodes;

    /* only check for interleaving if romio_cb_read isn't disabled */
    if (fd->hints->romio_cb_read != PNCIO_HINT_DISABLE) {
        /* For this process's request, calculate the file start and end
         * offsets. Note: end_offset points to the last byte-offset that will
         * be accessed, e.g., if start_offset=0 and 100 bytes to be read,
         * end_offset=99
         */
        if (fd->flat_file.size == 0) {
            start_offset = 0;
            end_offset = -1;
        }
        else if (fd->flat_file.count > 0) {
            start_offset = offset + fd->flat_file.off[0];
            end_offset   = fd->flat_file.off[fd->flat_file.count-1]
                         + fd->flat_file.len[fd->flat_file.count-1] - 1;
        }
        else {
            start_offset = offset;
            end_offset   = offset + fd->flat_file.size - 1;
        }

        /* each process communicates its start and end offsets to other
         * processes. The result is an array each of start and end offsets
         * stored in order of process rank. */
        st_offsets = (MPI_Offset *) NCI_Malloc(nprocs * 2 * sizeof(MPI_Offset));
        end_offsets = st_offsets + nprocs;

        MPI_Allgather(&start_offset, 1, MPI_OFFSET, st_offsets, 1, MPI_OFFSET,
                      fd->comm);
        MPI_Allgather(&end_offset, 1, MPI_OFFSET, end_offsets, 1, MPI_OFFSET,
                      fd->comm);

        /* Are the accesses of different processes interleaved? Below is a
         * rudimentary check for interleaving, but should suffice for the
         * moment. */
        for (int i = 1; i < nprocs; i++)
            if ((st_offsets[i] < end_offsets[i - 1]) &&
                (st_offsets[i] <= end_offsets[i]))
                interleave_count++;
    }

    if (fd->hints->romio_cb_read == PNCIO_HINT_DISABLE
        || (!interleave_count && (fd->hints->romio_cb_read == PNCIO_HINT_AUTO))) {
        /* switch to independent read */

        if (st_offsets != NULL) NCI_Free(st_offsets);

        if (buf_view.size == 0) return 0;

/* PnetCDF always sets this condition, i.e. when fileview is non-contiguous, offset in this call is always 0. */
if (fd->flat_file.count > 0) assert(offset == 0); /* not whole file visible */

        if (buf_view.is_contig && fd->flat_file.is_contig) {
            if (fd->flat_file.count > 0) offset += fd->flat_file.off[0];
            return PNCIO_ReadContig(fd, buf, buf_view.size, offset);
        }
        else
            return PNCIO_GEN_ReadStrided(fd, buf, buf_view, offset);
    }

    /* We're going to perform aggregation of I/O.  Here we call
     * PNCIO_Calc_file_domains() to determine what processes will handle I/O
     * to what regions.  We pass nprocs_for_coll into this function; it is
     * used to determine how many processes will perform I/O, which is also
     * the number of regions into which the range of bytes must be divided.
     * These regions are called "file domains", or FDs.
     *
     * When this function returns, fd_start, fd_end, fd_size, and
     * min_st_offset will be filled in.  fd_start holds the starting byte
     * location for each file domain.  fd_end holds the ending byte location.
     * min_st_offset holds the minimum byte location that will be accessed.
     *
     * Both fd_start[] and fd_end[] are indexed by an aggregator number; this
     * needs to be mapped to an actual rank in the communicator later.
     *
     */
    PNCIO_Calc_file_domains(st_offsets, end_offsets, nprocs, nprocs_for_coll,
                            &min_st_offset, &fd_start, &fd_end, &fd_size,
                            fd->hints->striping_unit);

    /* calculate where the portions of the access requests of this process
     * are located in terms of the file domains.  this could be on the same
     * process or on other processes.  this function fills in:
     * count_my_req_procs - number of processes (including this one) for which
     *     this process has requests in their file domain
     * count_my_req_per_proc - count of requests for each process, indexed
     *     by rank of the process
     * my_req[] - array of data structures describing the requests to be
     *     performed by each process (including self).  indexed by rank.
     * buf_idx[] - array of locations into which data can be directly moved;
     *     this is only valid for contiguous buffer case
     */
    PNCIO_Calc_my_req(fd, min_st_offset, fd_end, fd_size, nprocs,
                      &count_my_req_procs, &count_my_req_per_proc, &my_req,
                      &buf_idx);

    /* perform a collective communication in order to distribute the
     * data calculated above.  fills in the following:
     * count_others_req_procs - number of processes (including this
     *     one) which have requests in this process's file domain.
     * count_others_req_per_proc[] - number of separate contiguous
     *     requests from proc i lie in this process's file domain.
     */
    PNCIO_Calc_others_req(fd, count_my_req_procs, count_my_req_per_proc,
                          my_req, nprocs, myrank, &count_others_req_procs,
                          &count_others_req_per_proc, &others_req);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->read_timing[1] += MPI_Wtime() - curT;
#endif

    /* read data in sizes of no more than collective buffer size,
     * communicate, and fill user buf.
     */
    r_len = Read_and_exch(fd, buf, buf_view, nprocs, myrank, others_req,
                          min_st_offset, fd_size, fd_start, fd_end, buf_idx);
    if (r_len > 0) total_r_len += r_len;

    /* free all memory allocated for collective I/O */
    PNCIO_Free_my_req(count_my_req_per_proc, my_req, buf_idx);
    PNCIO_Free_others_req(count_others_req_per_proc, others_req);

    NCI_Free(st_offsets);
    NCI_Free(fd_start);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->read_timing[0] += MPI_Wtime() - curT;
#endif

    return (r_len < 0) ? r_len : total_r_len;
}

static
MPI_Offset Read_and_exch(PNCIO_File *fd, void *buf,
                         PNCIO_View buf_view, int nprocs,
                         int myrank, PNCIO_Access *others_req,
                         MPI_Offset min_st_offset, MPI_Offset fd_size,
                         MPI_Offset * fd_start, MPI_Offset * fd_end,
                         MPI_Aint * buf_idx)
{
/* Read in sizes of no more than coll_bufsize, an info parameter.
   Send data to appropriate processes.
   Place recd. data in user buf.
   The idea is to reduce the amount of extra memory required for
   collective I/O. If all data were read all at once, which is much
   easier, it would require temp space more than the size of user_buf,
   which is often unacceptable. For example, to read a distributed
   array from a file, where each local array is 8Mbytes, requiring
   at least another 8Mbytes of temp space is unacceptable. */

    int i, m, ntimes, max_ntimes;
    MPI_Offset st_loc = -1, end_loc = -1, off, done, real_off;
    char *read_buf = NULL, *tmp_buf;
    MPI_Count *curr_offlen_ptr, *count, *send_size, *recv_size;
    MPI_Count *partial_send, *recd_from_proc, *start_pos;
    /* Not convinced end_loc-st_loc couldn't be > int, so make these offsets */
    MPI_Offset real_size, size, for_curr_iter, for_next_iter;
    int rank;
    MPI_Aint coll_bufsize;
    MPI_Aint actual_recved_bytes = 0;
    MPI_Offset r_len;

/* calculate the number of reads of size coll_bufsize
   to be done by each process and the max among all processes.
   That gives the no. of communication phases as well.
   coll_bufsize is obtained from the hints object. */

    coll_bufsize = fd->hints->cb_buffer_size;

    /* grab some initial values for st_loc and end_loc */
    for (i = 0; i < nprocs; i++) {
        if (others_req[i].count) {
            st_loc = others_req[i].offsets[0];
            end_loc = others_req[i].offsets[0];
            break;
        }
    }

    /* now find the real values */
    for (i = 0; i < nprocs; i++)
        for (MPI_Count j = 0; j < others_req[i].count; j++) {
            st_loc = MIN(st_loc, others_req[i].offsets[j]);
            end_loc = MAX(end_loc, (others_req[i].offsets[j]
                                        + others_req[i].lens[j] - 1));
        }

    /* calculate ntimes, the number of times this process must perform I/O
     * operations in order to complete all the requests it has received.
     * the need for multiple I/O operations comes from the restriction that
     * we only use coll_bufsize bytes of memory for internal buffering.
     */
    if ((st_loc == -1) && (end_loc == -1)) {
        /* this process does no I/O. */
        ntimes = 0;
    } else {
        /* ntimes=ceiling_div(end_loc - st_loc + 1, coll_bufsize) */
        ntimes = (int) ((end_loc - st_loc + coll_bufsize) / coll_bufsize);
    }

    MPI_Allreduce(&ntimes, &max_ntimes, 1, MPI_INT, MPI_MAX, fd->comm);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->read_counter[0] = MAX(fd->read_counter[0], max_ntimes);
#endif

    read_buf = fd->io_buf;      /* Allocated at open time */

    curr_offlen_ptr = NCI_Calloc(nprocs * 7, sizeof(*curr_offlen_ptr));
    /* its use is explained below. calloc initializes to 0. */

    count = curr_offlen_ptr + nprocs;
    /* to store count of how many off-len pairs per proc are satisfied
     * in an iteration. */

    partial_send = count + nprocs;
    /* if only a portion of the last off-len pair is sent to a process
     * in a particular iteration, the length sent is stored here.
     * calloc initializes to 0. */

    send_size = partial_send + nprocs;
    /* total size of data to be sent to each proc. in an iteration */

    recv_size = send_size + nprocs;
    /* total size of data to be recd. from each proc. in an iteration.
     * Of size nprocs so that I can use MPI_Alltoall later. */

    recd_from_proc = recv_size + nprocs;
    /* amount of data recd. so far from each proc. Used in Fill_user_buffer.
     * initialized to 0 here. */

    start_pos = recd_from_proc + nprocs;
    /* used to store the starting value of curr_offlen_ptr[i] in
     * this iteration */

    done = 0;
    off = st_loc;
    for_curr_iter = for_next_iter = 0;

    MPI_Comm_rank(fd->comm, &rank);

    for (m = 0; m < ntimes; m++) {
        /* read buf of size coll_bufsize (or less) */
        /* go through all others_req and check if any are satisfied
         * by the current read */

        /* since MPI guarantees that displacements in filetypes are in
         * monotonically nondecreasing order, I can maintain a pointer
         * (curr_offlen_ptr) to
         * current off-len pair for each process in others_req and scan
         * further only from there. There is still a problem of filetypes
         * such as:  (1, 2, 3 are not process nos. They are just numbers for
         * three chunks of data, specified by a filetype.)
         *
         * 1  -------!--
         * 2    -----!----
         * 3       --!-----
         *
         * where ! indicates where the current read_size limitation cuts
         * through the filetype.  I resolve this by reading up to !, but
         * filling the communication buffer only for 1. I copy the portion
         * left over for 2 into a tmp_buf for use in the next
         * iteration. i.e., 2 and 3 will be satisfied in the next
         * iteration. This simplifies filling in the user's buf at the
         * other end, as only one off-len pair with incomplete data
         * will be sent. I also don't need to send the individual
         * offsets and lens along with the data, as the data is being
         * sent in a particular order. */

        /* off = start offset in the file for the data actually read in
         * this iteration
         * size = size of data read corresponding to off
         * real_off = off minus whatever data was retained in memory from
         * previous iteration for cases like 2, 3 illustrated above
         * real_size = size plus the extra corresponding to real_off
         * req_off = off in file for a particular contiguous request
         * minus what was satisfied in previous iteration
         * req_size = size corresponding to req_off */

        size = MIN(coll_bufsize, end_loc - st_loc + 1 - done);
        bool flag = false;
        for (i = 0; i < nprocs; i++) {
            if (others_req[i].count) {
                for (MPI_Count j = curr_offlen_ptr[i]; j < others_req[i].count; j++) {
                    MPI_Offset req_off;
                    if (partial_send[i]) {
                        req_off = others_req[i].offsets[j] + partial_send[i];
                    } else {
                        req_off = others_req[i].offsets[j];
                    }
                    if (req_off < off + size) {
                        flag = true;
                    }
                }
            }
        }
        if (flag) {
            /* This should be only reached by I/O aggregators only */
            r_len = PNCIO_ReadContig(fd, read_buf + for_curr_iter, size, off);
            if (r_len < 0) return r_len;
            size = r_len;
        }

        real_off = off - for_curr_iter;
        real_size = size + for_curr_iter;

        for (i = 0; i < nprocs; i++)
            count[i] = send_size[i] = 0;
        for_next_iter = 0;

        for (i = 0; i < nprocs; i++) {
            if (others_req[i].count) {
                start_pos[i] = curr_offlen_ptr[i];
                MPI_Count j = 0;
                for (j = curr_offlen_ptr[i]; j < others_req[i].count; j++) {
                    MPI_Offset req_off;
#ifdef HAVE_MPI_LARGE_COUNT
                    MPI_Offset req_len;
#else
                    int req_len;
#endif
                    if (partial_send[i]) {
                        /* this request may have been partially
                         * satisfied in the previous iteration. */
                        req_off = others_req[i].offsets[j] + partial_send[i];
                        req_len = others_req[i].lens[j] - partial_send[i];
                        partial_send[i] = 0;
                        /* modify the off-len pair to reflect this change */
                        others_req[i].offsets[j] = req_off;
                        others_req[i].lens[j] = req_len;
                    } else {
                        req_off = others_req[i].offsets[j];
                        req_len = others_req[i].lens[j];
                    }
                    if (req_off < real_off + real_size) {
                        count[i]++;
                        MPI_Aint addr;
                        MPI_Get_address(read_buf + req_off - real_off, &addr);
                        others_req[i].mem_ptrs[j] = addr;
                        send_size[i] += (MIN(real_off + real_size - req_off, req_len));

                        if (real_off + real_size - req_off < req_len) {
                            partial_send[i] = (real_off + real_size - req_off);
                            if ((j + 1 < others_req[i].count) &&
                                (others_req[i].offsets[j + 1] < real_off + real_size)) {
                                /* this is the case illustrated in the
                                 * figure above. */
                                for_next_iter = MAX(for_next_iter,
                                                        real_off + real_size -
                                                        others_req[i].offsets[j + 1]);
                                /* max because it must cover requests
                                 * from different processes */
                            }
                            break;
                        }
                    } else
                        break;
                }
                curr_offlen_ptr[i] = j;
            }
        }

        for_curr_iter = for_next_iter;

        MPI_Aint recved_bytes = 0;
        R_Exchange_data(fd, buf, buf_view, send_size, recv_size, count,
                        start_pos, partial_send, recd_from_proc, nprocs,
                        myrank, min_st_offset, fd_size, fd_start, fd_end,
                        others_req, m, buf_idx, &recved_bytes);
        actual_recved_bytes += recved_bytes;


        if (for_next_iter) {
            tmp_buf = (char *) NCI_Malloc(for_next_iter);
            memcpy(tmp_buf, read_buf + real_size - for_next_iter, for_next_iter);
            NCI_Free(fd->io_buf);
            fd->io_buf = (char *) NCI_Malloc(for_next_iter + coll_bufsize);
            memcpy(fd->io_buf, tmp_buf, for_next_iter);
            read_buf = fd->io_buf;
            NCI_Free(tmp_buf);
        }

        off += size;
        done += size;
    }

    for (i = 0; i < nprocs; i++)
        count[i] = send_size[i] = 0;
    for (m = ntimes; m < max_ntimes; m++) {
        /* nothing to send, but check for recv. */
        MPI_Aint recved_bytes = 0;
        R_Exchange_data(fd, buf, buf_view, send_size, recv_size, count,
                        start_pos, partial_send, recd_from_proc, nprocs,
                        myrank, min_st_offset, fd_size, fd_start, fd_end,
                        others_req, m, buf_idx, &recved_bytes);
        actual_recved_bytes += recved_bytes;
    }

    NCI_Free(curr_offlen_ptr);

    return actual_recved_bytes;
}

static void R_Exchange_data(PNCIO_File *fd, void *buf,
                            PNCIO_View buf_view,
                            MPI_Count * send_size, MPI_Count * recv_size,
                            MPI_Count * count, MPI_Count * start_pos,
                            MPI_Count * partial_send, MPI_Count * recd_from_proc, int nprocs,
                            int myrank,
                            MPI_Offset min_st_offset, MPI_Offset fd_size,
                            MPI_Offset * fd_start, MPI_Offset * fd_end,
                            PNCIO_Access * others_req, int iter,
                            MPI_Aint * buf_idx, MPI_Aint * actual_recved_bytes)
{
    int i, nprocs_recv, nprocs_send;
    char **recv_buf = NULL;
    size_t memLen;
    MPI_Request *requests;
    MPI_Datatype send_type;
    MPI_Status *statuses;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double curT = MPI_Wtime();
#endif

/* exchange send_size info so that each process knows how much to
   receive from whom and how much memory to allocate. */

    MPI_Alltoall(send_size, 1, MPI_COUNT, recv_size, 1, MPI_COUNT, fd->comm);

    nprocs_recv = 0;
    nprocs_send = 0;
    memLen = 0;
    for (i = 0; i < nprocs; i++) {
        memLen += recv_size[i];
        if (recv_size[i])
            nprocs_recv++;
        if (send_size[i])
            nprocs_send++;
    }

    requests = (MPI_Request *)
        NCI_Malloc((nprocs_send + nprocs_recv + 1) * sizeof(MPI_Request));
/* +1 to avoid a 0-size malloc */

/* post recvs. if buf_view.is_contig, data can be directly recd. into
   user buf at location given by buf_idx. else use recv_buf. */

    MPI_Count j = 0; // think of this as a counter of non-zero sends/recs
    if (buf_view.is_contig) {
        for (i = 0; i < nprocs; i++) {
            if (recv_size[i]) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(((char *) buf) + buf_idx[i], recv_size[i],
                            MPI_BYTE, i, 0, fd->comm, requests + j);
#else
                MPI_Irecv(((char *) buf) + buf_idx[i], recv_size[i],
                            MPI_BYTE, i, 0, fd->comm, requests + j);
#endif
                j++;
                buf_idx[i] += recv_size[i];
            }
        }
    } else {
        /* allocate memory for recv_buf and post receives */
        recv_buf = (char **) NCI_Malloc(nprocs * sizeof(char *));
        recv_buf[0] = (char *) NCI_Malloc(memLen);
        for (i = 1; i < nprocs; i++)
            recv_buf[i] = recv_buf[i - 1] + recv_size[i - 1];

        j = 0;
        for (i = 0; i < nprocs; i++) {
            if (recv_size[i]) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(recv_buf[i], recv_size[i], MPI_BYTE, i,
                            0, fd->comm, requests + j);
#else
                MPI_Irecv(recv_buf[i], recv_size[i], MPI_BYTE, i,
                            0, fd->comm, requests + j);
#endif
                j++;
            }
        }
    }

/* create derived datatypes and send data */

    j = 0;
    for (i = 0; i < nprocs; i++) {
        if (send_size[i]) {
            /* take care if the last off-len pair is a partial send */
            MPI_Offset tmp = 0;
            MPI_Count k = 0;
            if (partial_send[i]) {
                k = start_pos[i] + count[i] - 1;
                tmp = others_req[i].lens[k];
                others_req[i].lens[k] = partial_send[i];
            }
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(count[i],
                                       &(others_req[i].lens[start_pos[i]]),
                                       &(others_req[i].mem_ptrs[start_pos[i]]),
                                       MPI_BYTE, &send_type);
#else
            MPI_Type_create_hindexed(count[i],
                                     &(others_req[i].lens[start_pos[i]]),
                                     &(others_req[i].mem_ptrs[start_pos[i]]),
                                     MPI_BYTE, &send_type);
#endif
            /* absolute displacement; use MPI_BOTTOM in send */
            MPI_Type_commit(&send_type);
            MPI_Isend(MPI_BOTTOM, 1, send_type, i, 0,
                      fd->comm, requests + nprocs_recv + j);
            MPI_Type_free(&send_type);
            if (partial_send[i])
                others_req[i].lens[k] = tmp;
            j++;
        }
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->read_timing[4] += MPI_Wtime() - curT;
#endif


    /* +1 to avoid a 0-size malloc */
    statuses = (MPI_Status *) NCI_Malloc((nprocs_send + nprocs_recv + 1) * sizeof(MPI_Status));

    /* wait on the receives */
    if (nprocs_recv) {
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        curT = MPI_Wtime();
#endif
        MPI_Waitall(nprocs_recv, requests, statuses);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        if (fd->is_agg) fd->read_timing[3] += MPI_Wtime() - curT;
#endif

        *actual_recved_bytes = 0;
        j = 0;
        for (i = 0; i < nprocs; i++) {
            if (recv_size[i]) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count count_recved;
                MPI_Get_count_c(&statuses[j], MPI_BYTE, &count_recved);
#else
                int count_recved;
                MPI_Get_count(&statuses[j], MPI_BYTE, &count_recved);
#endif
                *actual_recved_bytes += count_recved;
                j++;
            }
        }

        /* if noncontiguous, to the copies from the recv buffers */
        if (!buf_view.is_contig)
            Fill_user_buffer(fd, buf, buf_view, recv_buf, recv_size,
                             recd_from_proc, nprocs, min_st_offset,
                             fd_size, fd_start, fd_end);
    }

    /* wait on the sends */
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    curT = MPI_Wtime();
#endif
#ifdef HAVE_MPI_STATUSES_IGNORE
    MPI_Waitall(nprocs_send, requests + nprocs_recv, MPI_STATUSES_IGNORE);
#else
    MPI_Waitall(nprocs_send, requests + nprocs_recv, statuses + nprocs_recv);
#endif
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->read_timing[3] += MPI_Wtime() - curT;
#endif

    NCI_Free(statuses);
    NCI_Free(requests);

    if (!buf_view.is_contig) {
        NCI_Free(recv_buf[0]);
        NCI_Free(recv_buf);
    }
}

#define BUF_INCR {                                                  \
    while (buf_incr) {                                              \
        size_in_buf = MIN(buf_incr, flat_buf_sz);                   \
        user_buf_idx += size_in_buf;                                \
        flat_buf_sz -= size_in_buf;                                 \
        buf_incr -= size_in_buf;                                    \
        if (buf_incr > 0 && flat_buf_sz == 0) {                     \
            flat_buf_idx++;                                         \
            user_buf_idx = buf_view.off[flat_buf_idx];              \
            flat_buf_sz = buf_view.len[flat_buf_idx];               \
        }                                                           \
    }                                                               \
}


#define BUF_COPY {                                                  \
    while (size) {                                                  \
        size_in_buf = MIN(size, flat_buf_sz);                       \
        memcpy(((char *) buf) + user_buf_idx,                       \
               &(recv_buf[p][recv_buf_idx[p]]), size_in_buf);       \
        recv_buf_idx[p] += size_in_buf;                             \
        user_buf_idx += size_in_buf;                                \
        flat_buf_sz -= size_in_buf;                                 \
        size -= size_in_buf;                                        \
        buf_incr -= size_in_buf;                                    \
        if (size > 0 && flat_buf_sz == 0) {                         \
            flat_buf_idx++;                                         \
            user_buf_idx = buf_view.off[flat_buf_idx];              \
            flat_buf_sz = buf_view.len[flat_buf_idx];               \
        }                                                           \
    }                                                               \
    BUF_INCR                                                        \
}

static void Fill_user_buffer(PNCIO_File *fd, void *buf,
                             PNCIO_View buf_view,
                             char **recv_buf,
                             MPI_Count * recv_size,
                             MPI_Count * recd_from_proc, int nprocs,
                             MPI_Offset min_st_offset,
                             MPI_Offset fd_size, MPI_Offset * fd_start,
                             MPI_Offset * fd_end)
{

/* this function is only called if buftype is not contig */

    int p, flat_buf_idx;
    MPI_Offset flat_buf_sz, size_in_buf, buf_incr, size;
    MPI_Offset off, user_buf_idx;
    MPI_Offset len, rem_len;
    MPI_Count *curr_from_proc, *done_from_proc, *recv_buf_idx;

/*  curr_from_proc[p] = amount of data recd from proc. p that has already
                        been accounted for so far
    done_from_proc[p] = amount of data already recd from proc. p and
                        filled into user buffer in previous iterations
    user_buf_idx = current location in user buffer
    recv_buf_idx[p] = current location in recv_buf of proc. p  */
    /* combining these three related arrays into a single memory allocation
     * (the "times 3" here) can help some highly noncontiguous workloads a bit */
    curr_from_proc = NCI_Malloc(nprocs * 3 * sizeof(*curr_from_proc));
    done_from_proc = curr_from_proc + nprocs;
    recv_buf_idx = done_from_proc + nprocs;

    for (int i = 0; i < nprocs; i++) {
        recv_buf_idx[i] = curr_from_proc[i] = 0;
        done_from_proc[i] = recd_from_proc[i];
    }

    user_buf_idx = buf_view.off[0];
    flat_buf_idx = 0;
    flat_buf_sz = buf_view.len[0];

    /* flat_buf_idx = current index into flattened buftype
     * flat_buf_sz = size of current contiguous component in
     * flattened buf */

    for (MPI_Count i = 0; i < fd->flat_file.count; i++) {
        off = fd->flat_file.off[i];
        rem_len = fd->flat_file.len[i];

        /* this request may span the file domains of more than one process */
        while (rem_len != 0) {
            len = rem_len;
            /* NOTE: len value is modified by PNCIO_Calc_aggregator() to be no
             * longer than the single region that processor "p" is responsible
             * for.
             */
            p = PNCIO_Calc_aggregator(fd, off, min_st_offset, &len, fd_size, fd_end);

            if (recv_buf_idx[p] < recv_size[p]) {
                if (curr_from_proc[p] + len > done_from_proc[p]) {
                    if (done_from_proc[p] > curr_from_proc[p]) {
                        size = MIN(curr_from_proc[p] + len - done_from_proc[p],
                                   recv_size[p] - recv_buf_idx[p]);
                        buf_incr = done_from_proc[p] - curr_from_proc[p];
                        BUF_INCR
                        buf_incr = curr_from_proc[p] + len - done_from_proc[p];
                        curr_from_proc[p] = done_from_proc[p] + size;
                        BUF_COPY
                    } else {
                        size = MIN(len, recv_size[p] - recv_buf_idx[p]);
                        buf_incr = len;
                        curr_from_proc[p] += size;
                        BUF_COPY
                    }
                } else {
                    curr_from_proc[p] += len;
                    buf_incr = len;
                    BUF_INCR
                }
            } else {
                buf_incr = len;
                BUF_INCR
            }
            off += len;
            rem_len -= len;
        }
    }
    for (int i = 0; i < nprocs; i++)
        if (recv_size[i])
            recd_from_proc[i] = curr_from_proc[i];

    NCI_Free(curr_from_proc);
}
