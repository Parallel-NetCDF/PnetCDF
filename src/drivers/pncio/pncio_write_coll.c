/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "pncio.h"

/* prototypes of functions used for collective writes only. */
static MPI_Offset Exch_and_write(PNCIO_File *fd, void *buf,
                           PNCIO_View buf_view, int nprocs, int myrank,
                           PNCIO_Access *others_req,
                           MPI_Offset min_st_offset, MPI_Offset fd_size,
                           MPI_Offset * fd_start, MPI_Offset * fd_end,
                           MPI_Aint * buf_idx);

static MPI_Offset W_Exchange_data(PNCIO_File *fd, void *buf, char *write_buf,
                            PNCIO_View buf_view,
                            MPI_Count * send_size, MPI_Count * recv_size,
                            MPI_Offset off, MPI_Count size,   /* 10 */
                            MPI_Count * count, MPI_Count * start_pos,
                            MPI_Count * partial_recv, MPI_Count *
                            sent_to_proc, int nprocs,
                            int myrank,
                            MPI_Offset min_st_offset, MPI_Offset fd_size,
                            MPI_Offset * fd_start, MPI_Offset * fd_end,
                            PNCIO_Access * others_req,
                            MPI_Count *send_buf_idx, MPI_Count *curr_to_proc,
                            MPI_Count *done_to_proc, int *hole, int iter,
                            MPI_Aint * buf_idx);

static void Fill_send_buffer(PNCIO_File *fd, void *buf,
                             PNCIO_View buf_view, char **send_buf,
                             MPI_Count *send_size, MPI_Request *requests,
                             MPI_Count *sent_to_proc, int nprocs, int myrank,
                             MPI_Offset min_st_offset,
                             MPI_Offset fd_size, MPI_Offset *fd_start,
                             MPI_Offset *fd_end, MPI_Count *send_buf_idx,
                             MPI_Count *curr_to_proc, MPI_Count *done_to_proc, int iter);

MPI_Offset PNCIO_GEN_WriteStridedColl(PNCIO_File *fd,
                                      const void *buf,
                                      PNCIO_View  buf_view,
                                      MPI_Offset  offset) /* relative to fileview */
{
    /* Uses a generalized version of the extended two-phase method described in
     * "An Extended Two-Phase Method for Accessing Sections of Out-of-Core
     * Arrays", Rajeev Thakur and Alok Choudhary, Scientific Programming,
     * (5)4:301--317, Winter 1996.
     * http://www.mcs.anl.gov/home/thakur/ext2ph.ps
     */

    PNCIO_Access *my_req;
    /* array of nprocs access structures, one for each other process in
     * whose file domain this process's request lies */

    PNCIO_Access *others_req;
    /* array of nprocs access structures, one for each other process
     * whose request lies in this process's file domain. */

    int i, nprocs, nprocs_for_coll, myrank, interleave_count=0;
    MPI_Aint *buf_idx = NULL;
    MPI_Count *count_my_req_per_proc, count_my_req_procs;
    MPI_Count *count_others_req_per_proc, count_others_req_procs;
    MPI_Offset start_offset, end_offset, fd_size, min_st_offset;
    MPI_Offset *st_offsets=NULL, *fd_start=NULL;
    MPI_Offset *fd_end=NULL, *end_offsets=NULL, w_len=0;

// printf("%s at %d: offset=%lld buf_view.size=%lld\n",__func__,__LINE__, offset,buf_view.size);

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
double curT = MPI_Wtime();
#endif

    /* the number of processes that actually perform I/O, nprocs_for_coll, is
     * stored in the hints off the PNCIO_File structure
     */
    nprocs_for_coll = fd->hints->cb_nodes;

    /* only check for interleaving if cb_write isn't disabled */
    if (fd->hints->cb_write != PNCIO_HINT_DISABLE) {
        /* For this process's request, calculate the file start and end
         * offsets. Note: end_offset points to the last byte-offset that will
         * be accessed, e.g., if start_offset=0 and 100 bytes to be read,
         * end_offset=99.
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

        /* Each process communicates its start and end offsets to other
         * processes. The result is an array each of start and end offsets
         * stored in order of process rank.
         */

        st_offsets = (MPI_Offset *) NCI_Malloc(nprocs * 2 * sizeof(MPI_Offset));
        end_offsets = st_offsets + nprocs;

        MPI_Allgather(&start_offset, 1, MPI_OFFSET, st_offsets, 1, MPI_OFFSET,
                      fd->comm);
        MPI_Allgather(&end_offset, 1, MPI_OFFSET, end_offsets, 1, MPI_OFFSET,
                      fd->comm);

        /* Are the accesses of different processes interleaved? Below is a
         * rudimentary check for interleaving, but should suffice for the
         * moment.
         */
        for (i = 1; i < nprocs; i++)
            if (st_offsets[i] < end_offsets[i - 1] &&
                st_offsets[i] <= end_offsets[i])
                interleave_count++;
    }

    if (fd->hints->cb_write == PNCIO_HINT_DISABLE ||
        (!interleave_count && (fd->hints->cb_write == PNCIO_HINT_AUTO))) {

        /* use independent accesses */
        if (fd->hints->cb_write != PNCIO_HINT_DISABLE)
            NCI_Free(st_offsets);
        if (buf_view.size == 0) return 0;

        /* offset is relative to fileview */
if (fd->flat_file.count > 0) assert(offset == 0); /* not whole file visible */

        if (buf_view.is_contig && fd->flat_file.is_contig) {
            if (fd->flat_file.count > 0) offset += fd->flat_file.off[0];
            w_len = PNCIO_WriteContig(fd, buf, buf_view.size, offset);
        }
        else
            w_len = PNCIO_GEN_WriteStrided(fd, buf, buf_view, offset);

        return w_len;
    }

// printf("%s at %d:\n",__func__,__LINE__);
/* Divide the I/O workload among "nprocs_for_coll" processes. This is
   done by (logically) dividing the file into file domains (FDs); each
   process may directly access only its own file domain. */

    PNCIO_Calc_file_domains(st_offsets, end_offsets, nprocs, nprocs_for_coll,
                            &min_st_offset, &fd_start, &fd_end, &fd_size,
                            fd->hints->striping_unit);

/* calculate what portions of the access requests of this process are
   located in what file domains */

    PNCIO_Calc_my_req(fd, min_st_offset, fd_start, fd_end, fd_size, nprocs,
                      &count_my_req_procs, &count_my_req_per_proc, &my_req,
                      &buf_idx);

/* based on everyone's my_req, calculate what requests of other
   processes lie in this process's file domain.
   count_others_req_procs = number of processes whose requests lie in
   this process's file domain (including this process itself)
   count_others_req_per_proc[i] indicates how many separate contiguous
   requests of proc. i lie in this process's file domain. */

    PNCIO_Calc_others_req(fd, count_my_req_procs, count_my_req_per_proc,
                          my_req, nprocs, myrank, &count_others_req_procs,
                          &count_others_req_per_proc, &others_req);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->write_timing[1] += MPI_Wtime() - curT;
#endif

/* exchange data and write in sizes of no more than coll_bufsize. */
    /* Cast away const'ness for the below function */
    w_len = Exch_and_write(fd, (char *) buf, buf_view, nprocs, myrank,
                           others_req, min_st_offset, fd_size, fd_start,
                           fd_end, buf_idx);

    /* If this collective write is followed by an independent write,
     * it's possible to have those subsequent writes on other processes
     * race ahead and sneak in before the read-modify-write completes.
     * We carry out a collective communication at the end here so no one
     * can start independent i/o before collective I/O completes.
     *
     * need to do some gymnastics with the error codes so that if something
     * went wrong, all processes report error, but if a process has a more
     * specific error code, we can still have that process report the
     * additional information */

    /* optimization: if only one process performing i/o, we can perform
     * a less-expensive Bcast
     */
    if (fd->hints->cb_nodes == 1)
        MPI_Bcast(&w_len, 1, MPI_OFFSET, fd->hints->ranklist[0], fd->comm);
    else
        MPI_Allreduce(MPI_IN_PLACE, &w_len, 1, MPI_OFFSET, MPI_MIN, fd->comm);

    /* free all memory allocated for collective I/O */
    PNCIO_Free_my_req(count_my_req_per_proc, my_req, buf_idx);
    PNCIO_Free_others_req(count_others_req_per_proc, others_req);

    NCI_Free(st_offsets);
    NCI_Free(fd_start);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->write_timing[0] += MPI_Wtime() - curT;
#endif

    /* w_len may not be the same as buf_view.size, because data sieving may
     * write more than requested.
     */
    return buf_view.size;
}

/* If successful, it returns the amount written. Otherwise a NetCDF error code
 * (negative value) is returned.
 */
static
MPI_Offset Exch_and_write(PNCIO_File *fd, void *buf, PNCIO_View buf_view,
                          int nprocs,
                          int myrank,
                          PNCIO_Access *others_req,
                          MPI_Offset min_st_offset, MPI_Offset fd_size,
                          MPI_Offset * fd_start, MPI_Offset * fd_end,
                          MPI_Aint * buf_idx)
{
/* Send data to appropriate processes and write in sizes of no more
   than coll_bufsize.
   The idea is to reduce the amount of extra memory required for
   collective I/O. If all data were written all at once, which is much
   easier, it would require temp space more than the size of user_buf,
   which is often unacceptable. For example, to write a distributed
   array to a file, where each local array is 8Mbytes, requiring
   at least another 8Mbytes of temp space is unacceptable. */

    /* Not convinced end_loc-st_loc couldn't be > int, so make these offsets */
    MPI_Offset size=0, w_len, total_w_len=0;
    int hole, i, m, ntimes, max_ntimes;
    MPI_Offset st_loc = -1, end_loc = -1, off, done, req_off;
    char *write_buf = NULL;
    MPI_Count *curr_offlen_ptr, *send_size, *count, req_len, *recv_size;
    MPI_Count *partial_recv, *sent_to_proc, *start_pos;
    int flag;
    MPI_Count *send_buf_idx, *curr_to_proc, *done_to_proc;
    int info_flag;
    MPI_Aint coll_bufsize;
    char *value;

    /* only I/O errors are currently reported */

/* calculate the number of writes of size coll_bufsize
   to be done by each process and the max among all processes.
   That gives the no. of communication phases as well. */

    value = (char *) NCI_Malloc((MPI_MAX_INFO_VAL + 1) * sizeof(char));
    MPI_Info_get(fd->info, "cb_buffer_size", MPI_MAX_INFO_VAL, value, &info_flag);
    coll_bufsize = atoi(value);
    NCI_Free(value);


    for (i = 0; i < nprocs; i++) {
        if (others_req[i].count) {
            st_loc = others_req[i].offsets[0];
            end_loc = others_req[i].offsets[0];
            break;
        }
    }

    for (i = 0; i < nprocs; i++)
        for (MPI_Count j = 0; j < others_req[i].count; j++) {
            st_loc = MIN(st_loc, others_req[i].offsets[j]);
            end_loc = MAX(end_loc, (others_req[i].offsets[j]
                                        + others_req[i].lens[j] - 1));
        }

/* ntimes=ceiling_div(end_loc - st_loc + 1, coll_bufsize)*/

    ntimes = (int) ((end_loc - st_loc + coll_bufsize) / coll_bufsize);

    if ((st_loc == -1) && (end_loc == -1)) {
        ntimes = 0;     /* this process does no writing. */
    }

    MPI_Allreduce(&ntimes, &max_ntimes, 1, MPI_INT, MPI_MAX, fd->comm);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->write_counter[0] = MAX(fd->write_counter[0], max_ntimes);
#endif

    write_buf = fd->io_buf;

    curr_offlen_ptr = NCI_Calloc(nprocs * 10, sizeof(*curr_offlen_ptr));
    /* its use is explained below. calloc initializes to 0. */

    count = curr_offlen_ptr + nprocs;
    /* to store count of how many off-len pairs per proc are satisfied
     * in an iteration. */

    partial_recv = count + nprocs;
    /* if only a portion of the last off-len pair is recd. from a process
     * in a particular iteration, the length recd. is stored here.
     * calloc initializes to 0. */

    send_size = partial_recv + nprocs;
    /* total size of data to be sent to each proc. in an iteration.
     * Of size nprocs so that I can use MPI_Alltoall later. */

    recv_size = send_size + nprocs;
    /* total size of data to be recd. from each proc. in an iteration. */

    sent_to_proc = recv_size + nprocs;
    /* amount of data sent to each proc so far. Used in
     * Fill_send_buffer. initialized to 0 here. */

    send_buf_idx = sent_to_proc + nprocs;
    curr_to_proc = send_buf_idx + nprocs;
    done_to_proc = curr_to_proc + nprocs;
    /* Above three are used in Fill_send_buffer */

    start_pos = done_to_proc + nprocs;
    /* used to store the starting value of curr_offlen_ptr[i] in
     * this iteration */

    done = 0;
    off = st_loc;
// printf("%s at %d: off=%lld buf_view.size=%lld ntimes=%d\n",__func__,__LINE__, off,buf_view.size,ntimes);

    for (m = 0; m < ntimes; m++) {
        /* go through all others_req and check which will be satisfied
         * by the current write */

        /* Note that MPI guarantees that displacements in filetypes are in
         * monotonically nondecreasing order and that, for writes, the
         * filetypes cannot specify overlapping regions in the file. This
         * simplifies implementation a bit compared to reads. */

        /* off = start offset in the file for the data to be written in
         * this iteration
         * size = size of data written (bytes) corresponding to off
         * req_off = off in file for a particular contiguous request
         * minus what was satisfied in previous iteration
         * req_size = size corresponding to req_off */

        /* first calculate what should be communicated */

        for (i = 0; i < nprocs; i++)
            count[i] = recv_size[i] = 0;

        size = MIN(coll_bufsize, end_loc - st_loc + 1 - done);

        for (i = 0; i < nprocs; i++) {
            if (others_req[i].count) {
                start_pos[i] = curr_offlen_ptr[i];
                MPI_Count j;
                for (j = curr_offlen_ptr[i]; j < others_req[i].count; j++) {
                    if (partial_recv[i]) {
                        /* this request may have been partially
                         * satisfied in the previous iteration. */
                        req_off = others_req[i].offsets[j] + partial_recv[i];
                        req_len = others_req[i].lens[j] - partial_recv[i];
                        partial_recv[i] = 0;
                        /* modify the off-len pair to reflect this change */
                        others_req[i].offsets[j] = req_off;
                        others_req[i].lens[j] = req_len;
                    } else {
                        req_off = others_req[i].offsets[j];
                        req_len = others_req[i].lens[j];
                    }
                    if (req_off < off + size) {
                        count[i]++;
                        if (myrank != i) {
                            MPI_Aint addr;
                            MPI_Get_address(write_buf + req_off - off, &addr);
                            others_req[i].mem_ptrs[j] = addr;
                        }
                        else
                            others_req[i].mem_ptrs[j] = req_off - off;
                        recv_size[i] += MIN(off + size - req_off, req_len);

                        if (off + size - req_off < req_len) {
                            partial_recv[i] = (off + size - req_off);

                            /* --BEGIN ERROR HANDLING-- */
                            if ((j + 1 < others_req[i].count) &&
                                (others_req[i].offsets[j + 1] < off + size)) {
                                /* This error should not happen to PnetCDF, as
                                 * fileview is checked before entering this
                                 * subroutine.
                                 */
                                fprintf(stderr, "Filetype specifies overlapping write regions (which is illegal according to the MPI-2 specification\n");
                                /* allow to continue since additional
                                 * communication might have to occur
                                 */
                                return NC_EFILE;
                            }
                            /* --END ERROR HANDLING-- */
                            break;
                        }
                    } else
                        break;
                }
                curr_offlen_ptr[i] = j;
            }
        }

        w_len = W_Exchange_data(fd, buf, write_buf, buf_view, send_size,
                                recv_size, off, size, count, start_pos,
                                partial_recv, sent_to_proc, nprocs, myrank,
                                min_st_offset, fd_size, fd_start, fd_end,
                                others_req, send_buf_idx, curr_to_proc,
                                done_to_proc, &hole, m, buf_idx);

        if (w_len < 0)
            return w_len;
        else
            total_w_len += w_len;

        flag = 0;
        for (i = 0; i < nprocs; i++)
            if (count[i])
                flag = 1;

        if (flag) {
            w_len = PNCIO_WriteContig(fd, write_buf, size, off);
            if (w_len < 0)
                return w_len;
            else
                total_w_len += w_len;
        }

        off += size;
        done += size;
    }

    for (i = 0; i < nprocs; i++)
        count[i] = recv_size[i] = 0;
    for (m = ntimes; m < max_ntimes; m++) {
        /* nothing to recv, but check for send. */
        w_len = W_Exchange_data(fd, buf, write_buf, buf_view, send_size,
                                recv_size, off, size, count, start_pos,
                                partial_recv, sent_to_proc, nprocs, myrank,
                                min_st_offset, fd_size, fd_start, fd_end,
                                others_req, send_buf_idx, curr_to_proc,
                                done_to_proc, &hole, m, buf_idx);
        if (w_len < 0)
            return w_len;
        else
            total_w_len += w_len;
    }

    NCI_Free(curr_offlen_ptr);

    return total_w_len;
}


/* Sets error_code to MPI_SUCCESS if successful, or creates an error code
 * in the case of error.
 */
static
MPI_Offset W_Exchange_data(PNCIO_File *fd, void *buf, char *write_buf,
                           PNCIO_View buf_view,
                           MPI_Count *send_size, MPI_Count *recv_size,
                           MPI_Offset off, MPI_Count size,
                           MPI_Count *count, MPI_Count * start_pos,
                           MPI_Count *partial_recv,
                           MPI_Count *sent_to_proc, int nprocs,
                           int myrank,
                           MPI_Offset min_st_offset,
                           MPI_Offset fd_size,
                           MPI_Offset * fd_start, MPI_Offset * fd_end,
                           PNCIO_Access * others_req,
                           MPI_Count * send_buf_idx, MPI_Count * curr_to_proc,
                           MPI_Count * done_to_proc, int *hole, int iter,
                           MPI_Aint *buf_idx)
{
    int i, j, nprocs_recv, nprocs_send, err=NC_NOERR;
    MPI_Count *tmp_len;
    char **send_buf = NULL;
    MPI_Request *requests, *send_req;
    MPI_Datatype *recv_types, self_recv_type = MPI_DATATYPE_NULL;
    MPI_Status *statuses, status;
    MPI_Count sum, *srt_len = NULL;
    int num_rtypes, nreqs;
    MPI_Offset *srt_off = NULL;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
double curT = MPI_Wtime();
#endif

/* exchange recv_size info so that each process knows how much to
   send to whom. */

    MPI_Alltoall(recv_size, 1, MPI_COUNT, send_size, 1, MPI_COUNT, fd->comm);

    /* create derived datatypes for recv */

    nprocs_send = 0;
    nprocs_recv = 0;
    sum = 0;
    for (i = 0; i < nprocs; i++) {
        sum += count[i];
        if (recv_size[i])
            nprocs_recv++;
        if (send_size[i])
            nprocs_send++;
    }

    recv_types = (MPI_Datatype *) NCI_Malloc((nprocs_recv + 1) * sizeof(MPI_Datatype));
    /* +1 to avoid a 0-size malloc */

    tmp_len = NCI_Malloc(nprocs * sizeof(*tmp_len));
    j = 0;
    for (i = 0; i < nprocs; i++) {
        if (recv_size[i]) {
            MPI_Datatype *dtype;
            dtype = (i != myrank) ? (recv_types + j) : (&self_recv_type);

            if (partial_recv[i]) {
                /* take care if the last off-len pair is a partial recv */
                MPI_Count k = start_pos[i] + count[i] - 1;
                tmp_len[i] = others_req[i].lens[k];
                others_req[i].lens[k] = partial_recv[i];
            }
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(count[i],
                                       &(others_req[i].lens[start_pos[i]]),
                                       &(others_req[i].mem_ptrs[start_pos[i]]),
                                       MPI_BYTE, dtype);
#else
            MPI_Type_create_hindexed(count[i],
                                     &(others_req[i].lens[start_pos[i]]),
                                     &(others_req[i].mem_ptrs[start_pos[i]]),
                                     MPI_BYTE, dtype);
#endif
            /* absolute displacements; use MPI_BOTTOM in recv */
            MPI_Type_commit(dtype);
            if (i != myrank)
                j++;
        }
    }
    num_rtypes = j;     /* number of non-self receive datatypes created */

    /* To avoid a read-modify-write, check if there are holes in the
     * data to be written. For this, merge the (sorted) offset lists
     * others_req using a heap-merge. */

/* TODO: PNCIO_Heap_merge is expensive, borrow codes from ad_lustre_wrcoll.c to skip it when possible */

    /* valgrind-detcted optimization: if there is no work on this process we do
     * not need to search for holes */
    if (sum) {
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        double timing = MPI_Wtime();
#endif
        srt_off = (MPI_Offset *) NCI_Malloc(sum * sizeof(MPI_Offset));
        srt_len = NCI_Malloc(sum * sizeof(*srt_len));

        PNCIO_Heap_merge(others_req, count, srt_off, srt_len, start_pos,
                         nprocs, nprocs_recv, sum);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        if (fd->is_agg) fd->write_timing[5] += MPI_Wtime() - timing;
#endif
    }

    /* for partial recvs, restore original lengths */
    for (i = 0; i < nprocs; i++)
        if (partial_recv[i]) {
            MPI_Count k = start_pos[i] + count[i] - 1;
            others_req[i].lens[k] = tmp_len[i];
        }
    NCI_Free(tmp_len);

    /* check if there are any holes. If yes, must do read-modify-write.
     * holes can be in three places.  'middle' is what you'd expect: the
     * processes are operating on noncontigous data.  But holes can also show
     * up at the beginning or end of the file domain (see John Bent ROMIO REQ
     * #835). Missing these holes would result in us writing more data than
     * received by everyone else. */

    *hole = 0;
    if (sum) {
        if (off != srt_off[0])  /* hole at the front */
            *hole = 1;
        else {  /* coalesce the sorted offset-length pairs */
            for (i = 1; i < sum; i++) {
                if (srt_off[i] <= srt_off[0] + srt_len[0]) {
                    MPI_Count new_len = srt_off[i] + srt_len[i] - srt_off[0];
                    if (new_len > srt_len[0])
                        srt_len[0] = new_len;
                } else
                    break;
            }
            if (i < sum || size != srt_len[0])  /* hole in middle or end */
                *hole = 1;
        }

        NCI_Free(srt_off);
        NCI_Free(srt_len);
    }

    if (nprocs_recv) {
        if (*hole) {
            MPI_Offset r_len;
            r_len = PNCIO_ReadContig(fd, write_buf, size, off);
            if (r_len < 0) return r_len;
        }
    }

    if (fd->atomicity) {
        /* nreqs is the number of Isend and Irecv to be posted */
        nreqs = (send_size[myrank]) ? (nprocs_send - 1) : nprocs_send;
        requests = (MPI_Request *) NCI_Malloc((nreqs + 1) * sizeof(MPI_Request));
        send_req = requests;
    } else {
        nreqs = nprocs_send + nprocs_recv;
        if (send_size[myrank])  /* NO send to and recv from self */
            nreqs -= 2;
        requests = (MPI_Request *) NCI_Malloc((nreqs + 1) * sizeof(MPI_Request));
        /* +1 to avoid a 0-size malloc */

        /* post receives */
        j = 0;
        for (i = 0; i < nprocs; i++) {
            if (recv_size[i] == 0)
                continue;
            if (i != myrank) {
                MPI_Irecv(MPI_BOTTOM, 1, recv_types[j], i, 0,
                          fd->comm, requests + j);
                j++;
            } else if (buf_view.is_contig) {
                /* sen/recv to/from self uses MPI_Unpack() */
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count position=0;
                MPI_Unpack_c((char *) buf + buf_idx[i], recv_size[i], &position,
                             write_buf, 1, self_recv_type, MPI_COMM_SELF);
#else
                int position = 0;
                assert(recv_size[i] < INT_MAX);
                MPI_Unpack((char *) buf + buf_idx[i], (int)recv_size[i], &position,
                           write_buf, 1, self_recv_type, MPI_COMM_SELF);
#endif
                buf_idx[i] += recv_size[i];
            }
        }
        send_req = requests + j;
    }

/* post sends. if buf_view.is_contig, data can be directly sent from
   user buf at location given by buf_idx. else use send_buf. */

    if (buf_view.is_contig) {
        j = 0;
        for (i = 0; i < nprocs; i++)
            if (send_size[i] && i != myrank) {
                assert(buf_idx[i] != -1);
#if MPI_VERSION >= 4
                MPI_Isend_c((char *) buf + buf_idx[i], send_size[i],
                            MPI_BYTE, i, 0, fd->comm, send_req + j);
#else
                MPI_Isend((char *) buf + buf_idx[i], send_size[i],
                            MPI_BYTE, i, 0, fd->comm, send_req + j);
#endif
                j++;
                buf_idx[i] += send_size[i];
            }
    } else if (nprocs_send) {
        /* buftype is not contig */
        size_t msgLen = 0;
        for (i = 0; i < nprocs; i++)
            msgLen += send_size[i];
        send_buf = (char **) NCI_Malloc(nprocs * sizeof(char *));
        send_buf[0] = (char *) NCI_Malloc(msgLen * sizeof(char));
        for (i = 1; i < nprocs; i++)
            send_buf[i] = send_buf[i - 1] + send_size[i - 1];

        Fill_send_buffer(fd, buf, buf_view, send_buf, send_size, send_req,
                         sent_to_proc, nprocs, myrank, min_st_offset, fd_size,
                         fd_start, fd_end, send_buf_idx, curr_to_proc,
                         done_to_proc, iter);

        /* the send is done in Fill_send_buffer */
    }

    if (fd->atomicity) {
        /* In atomic mode, we must use blocking receives to receive data in the
         * same increasing order of MPI process rank IDs,
         */
        j = 0;
        for (i = 0; i < nprocs; i++) {
            if (recv_size[i] == 0)
                continue;
            if (i != myrank) {
                MPI_Recv(MPI_BOTTOM, 1, recv_types[j++], i, 0,
                         fd->comm, &status);
            } else {
                /* sen/recv to/from self uses MPI_Unpack() */
                char *ptr = (buf_view.is_contig) ? (char *) buf + buf_idx[i] : send_buf[i];
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count position=0;
                MPI_Unpack_c(ptr, recv_size[i], &position, write_buf, 1, self_recv_type,
                             MPI_COMM_SELF);
#else
                int position = 0;
                assert(recv_size[i] < INT_MAX);
                MPI_Unpack(ptr, (int)recv_size[i], &position, write_buf, 1, self_recv_type,
                           MPI_COMM_SELF);
#endif
                buf_idx[i] += recv_size[i];
            }
        }
    } else if (!buf_view.is_contig && recv_size[myrank]) {
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count position=0;
        MPI_Unpack_c(send_buf[myrank], recv_size[myrank], &position, write_buf, 1, self_recv_type,
                     MPI_COMM_SELF);
#else
        int position = 0;
        assert(recv_size[myrank] < INT_MAX);
        MPI_Unpack(send_buf[myrank], (int)recv_size[myrank], &position, write_buf, 1, self_recv_type,
                   MPI_COMM_SELF);
#endif
    }

    for (i = 0; i < num_rtypes; i++)
        MPI_Type_free(recv_types + i);
    NCI_Free(recv_types);

    if (self_recv_type != MPI_DATATYPE_NULL)
        MPI_Type_free(&self_recv_type);

#ifdef HAVE_MPI_STATUSES_IGNORE
    statuses = MPI_STATUSES_IGNORE;
#else
    statuses = (MPI_Status *) NCI_Malloc(nreqs * sizeof(MPI_Status));
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->write_timing[4] += MPI_Wtime() - curT;
    curT = MPI_Wtime();
#endif
    MPI_Waitall(nreqs, requests, statuses);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->write_timing[3] += MPI_Wtime() - curT;
#endif

#ifndef HAVE_MPI_STATUSES_IGNORE
    NCI_Free(statuses);
#endif
    NCI_Free(requests);
    if (!buf_view.is_contig && nprocs_send) {
        NCI_Free(send_buf[0]);
        NCI_Free(send_buf);
    }

    return err;
}

#define BUF_INCR \
{ \
    while (buf_incr) { \
        size_in_buf = MIN(buf_incr, flat_buf_sz); \
        user_buf_idx += size_in_buf; \
        flat_buf_sz -= size_in_buf; \
        buf_incr -= size_in_buf; \
        if (buf_incr > 0 && flat_buf_sz == 0) { \
            flat_buf_idx++; \
            user_buf_idx = buf_view.off[flat_buf_idx]; \
            flat_buf_sz = buf_view.len[flat_buf_idx]; \
        } \
    } \
}

#define BUF_COPY \
{ \
    while (size) { \
        size_in_buf = MIN(size, flat_buf_sz); \
        memcpy(&(send_buf[p][send_buf_idx[p]]), \
               ((char *) buf) + user_buf_idx, size_in_buf); \
        send_buf_idx[p] += size_in_buf; \
        user_buf_idx += size_in_buf; \
        flat_buf_sz -= size_in_buf; \
        size -= size_in_buf; \
        buf_incr -= size_in_buf; \
        if (size > 0 && flat_buf_sz == 0) { \
            flat_buf_idx++; \
            user_buf_idx = buf_view.off[flat_buf_idx]; \
            flat_buf_sz = buf_view.len[flat_buf_idx]; \
        } \
    } \
    BUF_INCR \
}

static
void Fill_send_buffer(PNCIO_File *fd, void *buf,
                      PNCIO_View buf_view, char **send_buf,
                      MPI_Count * send_size,
                      MPI_Request * requests, MPI_Count * sent_to_proc,
                      int nprocs, int myrank,
                      MPI_Offset min_st_offset, MPI_Offset fd_size,
                      MPI_Offset * fd_start, MPI_Offset * fd_end,
                      MPI_Count * send_buf_idx, MPI_Count * curr_to_proc,
                      MPI_Count * done_to_proc, int iter)
{
/* this function is only called if buftype is not contig */

    int p, jj;
    MPI_Offset flat_buf_idx, flat_buf_sz, size_in_buf, buf_incr, size;
    MPI_Offset off, len, rem_len, user_buf_idx;

/*  curr_to_proc[p] = amount of data sent to proc. p that has already
    been accounted for so far
    done_to_proc[p] = amount of data already sent to proc. p in
    previous iterations
    user_buf_idx = current location in user buffer
    send_buf_idx[p] = current location in send_buf of proc. p  */

    for (MPI_Count i = 0; i < nprocs; i++) {
        send_buf_idx[i] = curr_to_proc[i] = 0;
        done_to_proc[i] = sent_to_proc[i];
    }
    jj = 0;

    user_buf_idx = buf_view.off[0];
    flat_buf_idx = 0;
    flat_buf_sz = buf_view.len[0];

    /* flat_buf_idx = current index into flattened buftype
     * flat_buf_sz = size of current contiguous component in
     * flattened buf */

    for (MPI_Count i = 0; i < fd->flat_file.count; i++) {
        off = fd->flat_file.off[i];
        rem_len = fd->flat_file.len[i];

        /*this request may span the file domains of more than one process */
        while (rem_len != 0) {
            len = rem_len;
            /* NOTE: len value is modified by PNCIO_Calc_aggregator() to be no
             * longer than the single region that processor "p" is responsible
             * for.
             */
            p = PNCIO_Calc_aggregator(fd, off, min_st_offset, &len, fd_size, fd_end);

            if (send_buf_idx[p] < send_size[p]) {
                if (curr_to_proc[p] + len > done_to_proc[p]) {
                    if (done_to_proc[p] > curr_to_proc[p]) {
                        size = MIN(curr_to_proc[p] + len -
                                       done_to_proc[p], send_size[p] - send_buf_idx[p]);
                        buf_incr = done_to_proc[p] - curr_to_proc[p];
                        BUF_INCR
                        buf_incr = curr_to_proc[p] + len - done_to_proc[p];
                        /* ok to cast: bounded by cb buffer size */
                        curr_to_proc[p] = done_to_proc[p] + size;
                        BUF_COPY
                    } else {
                        size = MIN(len, send_size[p] - send_buf_idx[p]);
                        buf_incr = len;
                        curr_to_proc[p] += size;
                        BUF_COPY
                    }
                    if (send_buf_idx[p] == send_size[p] && p != myrank) {
#if MPI_VERSION >= 4
                        MPI_Isend_c(send_buf[p], send_size[p], MPI_BYTE, p,
                                    0, fd->comm, &requests[jj++]);
#else
                        MPI_Isend(send_buf[p], send_size[p], MPI_BYTE, p,
                                    0, fd->comm, &requests[jj++]);
#endif
                    }
                } else {
                    curr_to_proc[p] += len;
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
    for (int i = 0; i < nprocs; i++) {
        if (send_size[i]) {
            sent_to_proc[i] = curr_to_proc[i];
        }
    }
}
