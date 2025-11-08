/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

static int use_alltoallw;

#ifdef HAVE_MPI_LARGE_COUNT
#define MEMCPY_UNPACK(x, inbuf, start, count, outbuf) {          \
    int _k;                                                      \
    char *_ptr = (inbuf);                                        \
    MPI_Count   *mem_ptrs = others_req[x].mem_ptrs + (start);    \
    MPI_Offset *mem_lens = others_req[x].lens     + (start);    \
    for (_k=0; _k<count; _k++) {                                 \
        memcpy((outbuf) + mem_ptrs[_k], _ptr, mem_lens[_k]);     \
        _ptr += mem_lens[_k];                                    \
    }                                                            \
}
#else
#define MEMCPY_UNPACK(x, inbuf, start, count, outbuf) {          \
    int _k;                                                      \
    char *_ptr = (inbuf);                                        \
    MPI_Aint *mem_ptrs = others_req[x].mem_ptrs + (start);       \
    int      *mem_lens = others_req[x].lens     + (start);       \
    for (_k=0; _k<count; _k++) {                                 \
        memcpy((outbuf) + mem_ptrs[_k], _ptr, mem_lens[_k]);     \
        _ptr += mem_lens[_k];                                    \
    }                                                            \
}
#endif

typedef struct {
    MPI_Count    num; /* number of elements in the above off-len list */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset *off; /* list of write offsets by this rank in round m */
    MPI_Count   *len; /* list of write lengths by this rank in round m */
#else
    MPI_Offset  *off; /* list of write offsets by this rank in round m */
    int         *len; /* list of write lengths by this rank in round m */
#endif
} off_len_list;

typedef struct {
    MPI_Count   count; /* number displacement-length pairs */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count  *disp;  /* [count]: displacement */
    MPI_Count  *len;   /* [count]: size in bytes */
#else
    MPI_Aint   *disp;  /* [count]: displacement */
    int        *len;   /* [count]: size in bytes */
#endif
} disp_len_list;

/* prototypes of functions used for collective writes only. */
static MPI_Offset LUSTRE_Exch_and_write(PNCIO_File *fd,
                                        const void *buf,
                                        PNCIO_View  buf_view,
                                        PNCIO_Access *others_req,
                                        PNCIO_Access *my_req,
                                        MPI_Offset min_st_loc,
                                        MPI_Offset max_end_loc,
                                        MPI_Offset **buf_idx);

static void LUSTRE_Fill_send_buffer(PNCIO_File *fd, const void *buf,
                                    PNCIO_View *buf_view,
                                    char **send_buf,
                                    size_t send_total_size,
                                    const MPI_Count *send_size,
                                    char **self_buf,
                                    disp_len_list *send_list);

static int Exchange_data_recv(PNCIO_File           *fd,
                              const void           *buf,
                                    char           *write_buf,
                                    char          **recv_buf,
                              const PNCIO_View     *buf_view,
                              const MPI_Count      *recv_size,
                                    MPI_Offset      range_off,
                                    MPI_Count       range_size,
                              const MPI_Count      *recv_count,
                              const MPI_Count      *start_pos,
                              const PNCIO_Access   *others_req,
                              const MPI_Offset     *buf_idx,
                                    off_len_list   *srt_off_len,
                                    disp_len_list  *recv_list);

static void Exchange_data_send(      PNCIO_File     *fd,
                               const void           *buf,
                                     char           *write_buf,
                                     char          **send_buf_ptr,
                                     PNCIO_View     *buf_view,
                               const MPI_Count      *send_size,
                                     MPI_Count       self_count,
                                     MPI_Count       start_pos,
                               const PNCIO_Access   *others_req,
                               const MPI_Offset     *buf_idx,
                                     disp_len_list  *send_list);

static
int LUSTRE_Calc_aggregator(PNCIO_File *fd,
                           MPI_Offset off,
#ifdef HAVE_MPI_LARGE_COUNT
                           MPI_Offset *len
#else
                           int        *len
#endif
)
{
    MPI_Offset avail_bytes, stripe_id;

    stripe_id = off / fd->hints->striping_unit;

    avail_bytes = (stripe_id + 1) * fd->hints->striping_unit - off;
    if (avail_bytes < *len) {
        /* The request [off, off+len) has only [off, off+avail_bytes) part
         * falling into aggregator's file domain */
        *len = avail_bytes;
    }
    /* return the index to ranklist[] */
    return (stripe_id % fd->hints->cb_nodes);
}

/*----< LUSTRE_Calc_my_req() >-----------------------------------------------*/
/* calculates what portions of the read/write requests of this process fall
 * into the file domains of all I/O aggregators.
 *   IN: fd->flat_file: this rank's flattened write requests
 *       fd->flat_file.count: number of noncontiguous offset-length file requests
 *       fd->flat_file.off[fd->flat_file.count] file offsets of individual
 *       noncontiguous requests.
 *       fd->flat_file.len[fd->flat_file.count] lengths of individual
 *       noncontiguous requests.
 *   IN: buf_is_contig: whether the write buffer is contiguous or not
 *   OUT: my_req_ptr[cb_nodes] offset-length pairs of this process's requests
 *        fall into the file domain of each aggregator
 *   OUT: buf_idx_ptr[cb_nodes] index pointing to the starting location in
 *        user_buf for data to be sent to each aggregator.
 */
static
void LUSTRE_Calc_my_req(PNCIO_File    *fd,
                        int            buf_is_contig,
                        PNCIO_Access **my_req_ptr,
                        MPI_Offset   **buf_idx)
{
    int aggr, *aggr_ranks, cb_nodes;
    MPI_Count i, l;
    size_t nelems, alloc_sz;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset rem_len, avail_len, *avail_lens;
#else
    int rem_len, avail_len, *avail_lens;
#endif
    MPI_Offset curr_idx, off;
    PNCIO_Access *my_req;

    cb_nodes = fd->hints->cb_nodes;

    /* my_req[i].count gives the number of contiguous requests of this process
     * that fall in aggregator i's file domain (not process MPI rank i).
     */
    my_req = (PNCIO_Access *) NCI_Calloc(cb_nodes, sizeof(PNCIO_Access));
    *my_req_ptr = my_req;

    /* First pass is just to calculate how much space is needed to allocate
     * my_req.
     */
#ifdef HAVE_MPI_LARGE_COUNT
    alloc_sz = sizeof(int) + sizeof(MPI_Offset);
    aggr_ranks = (int*) NCI_Malloc(alloc_sz * fd->flat_file.count);
    avail_lens = (MPI_Offset*) (aggr_ranks + fd->flat_file.count);
#else
    alloc_sz = sizeof(int) * 2;
    aggr_ranks = (int*) NCI_Malloc(alloc_sz * fd->flat_file.count);
    avail_lens = aggr_ranks + fd->flat_file.count;
#endif

    /* Note that MPI standard (MPI 3.1 Chapter 13.1.1 and MPI 4.0 Chapter
     * 14.1.1) requires that the typemap displacements of etype and
     * filetype are non-negative and monotonically non-decreasing. This
     * makes fd->flat_file.off[] to be monotonically non-decreasing.
     */

/*
Alternative: especially for when fd->flat_file.count is large
1 This rank's aggregate file access region is from start_offset to end_offset.
2 start with the 1st aggregator ID and keep assign aggregator until next stripe.
  This can avoid too many calls to LUSTRE_Calc_aggregator()
*/

    /* nelems will be the number of offset-length pairs for my_req[] */
    nelems = 0;
    for (i = 0; i < fd->flat_file.count; i++) {
        /* short circuit offset/len processing if zero-byte read/write. */
        if (fd->flat_file.len[i] == 0)
            continue;

        off = fd->flat_file.off[i];
        avail_len = fd->flat_file.len[i];
        /* LUSTRE_Calc_aggregator() modifies the value of 'avail_len' to the
         * amount that is only covered by the aggr's file domain. The remaining
         * (tail) will continue to be processed to determine to whose file
         * domain it belongs. As LUSTRE_Calc_aggregator() can be expensive for
         * large value of fd->flat_file.count, we keep a copy of the returned
         * values of 'aggr' and 'avail_len' in aggr_ranks[] and avail_lens[] to
         * be used in the next for loop (not next iteration).
         *
         * Note the returned value in 'aggr' is the index to ranklist[], i.e.
         * the 'aggr'th element of array ranklist[], rather than the
         * aggregator's MPI rank ID in fd->comm.
         */
        aggr = LUSTRE_Calc_aggregator(fd, off, &avail_len);
        aggr_ranks[i] = aggr;          /* first aggregator ID of this request */
        avail_lens[i] = avail_len;     /* length covered, may be < fd->flat_file.len[i] */
        assert(aggr >= 0 && aggr <= cb_nodes);
        my_req[aggr].count++; /* increment for aggregator aggr */
        nelems++;             /* true number of noncontiguous requests
                               * in terms of file domains */

        /* rem_len is the amount of ith offset-length pair that is not covered
         * by aggregator aggr's file domain.
         */
        rem_len = fd->flat_file.len[i] - avail_len;
        assert(rem_len >= 0);

        while (rem_len > 0) {
            off += avail_len;    /* move forward to first remaining byte */
            avail_len = rem_len; /* save remaining size, pass to calc */
            aggr = LUSTRE_Calc_aggregator(fd, off, &avail_len);
            my_req[aggr].count++;
            nelems++;
            rem_len -= avail_len;/* reduce remaining length by amount from fd */
        }
    }

    /* allocate space for buf_idx.
     * buf_idx is relevant only if buftype is contiguous. buf_idx[i] gives the
     * starting index in user_buf where data will be sent to aggregator 'i'.
     * This allows sends to be done without extra buffer.
     */
    if (buf_idx != NULL && buf_is_contig) {
        buf_idx[0] = (MPI_Offset *) NCI_Malloc(nelems * sizeof(MPI_Offset));
        for (i = 1; i < cb_nodes; i++)
            buf_idx[i] = buf_idx[i - 1] + my_req[i - 1].count;
    }

    /* allocate space for my_req and its members offsets and lens */
#ifdef HAVE_MPI_LARGE_COUNT
    alloc_sz = sizeof(MPI_Offset) * 2;
    my_req[0].offsets = (MPI_Offset*) NCI_Malloc(alloc_sz * nelems);
    my_req[0].lens    = my_req[0].offsets + my_req[0].count;
    for (i=1; i<cb_nodes; i++) {
        my_req[i].offsets = my_req[i-1].offsets + my_req[i-1].count * 2;
        my_req[i].lens    = my_req[i].offsets + my_req[i].count;
        my_req[i-1].count = 0; /* reset, will increase where needed later */
    }
    my_req[cb_nodes-1].count = 0;
#else
    alloc_sz = sizeof(MPI_Offset) + sizeof(int);
    my_req[0].offsets = (MPI_Offset*) NCI_Malloc(alloc_sz * nelems);
    my_req[0].lens    = (int*) (my_req[0].offsets + my_req[0].count);

    char *ptr = (char*) my_req[0].offsets + alloc_sz * my_req[0].count;
    for (i=1; i<cb_nodes; i++) {
        my_req[i].offsets = (MPI_Offset*)ptr;
        ptr += sizeof(MPI_Offset) * my_req[i].count;
        my_req[i].lens = (int*)ptr;
        ptr += sizeof(int) * my_req[i].count;
        my_req[i].count = 0; /* reset, will be incremented where needed later */
    }
    my_req[cb_nodes-1].count = 0;
#endif

    for (i=0; i<cb_nodes; i++)
        my_req[i].count = 0; /* reset, will be incremented where needed later */

    /* now fill in my_req */
    curr_idx = 0;
    for (i = 0; i < fd->flat_file.count; i++) {
        /* short circuit offset/len processing if zero-byte read/write. */
        if (fd->flat_file.len[i] == 0)
            continue;

        off = fd->flat_file.off[i];
        aggr = aggr_ranks[i];
        assert(aggr >= 0 && aggr <= cb_nodes);
        avail_len = avail_lens[i];

        l = my_req[aggr].count;
        if (buf_idx != NULL && buf_is_contig) {
            buf_idx[aggr][l] = curr_idx;
            curr_idx += avail_len;
        }
        rem_len = fd->flat_file.len[i] - avail_len;

        /* Each my_req[i] contains the number of this process's noncontiguous
         * requests that fall into aggregator aggr's file domain.
         * my_req[aggr].offsets[] and my_req[aggr].lens store the offsets and
         * lengths of the requests.
         */
        my_req[aggr].offsets[l] = off;
        my_req[aggr].lens[l] = avail_len;
        my_req[aggr].count++;

        while (rem_len != 0) {
            off += avail_len;
            avail_len = rem_len;
            aggr = LUSTRE_Calc_aggregator(fd, off, &avail_len);
            assert(aggr >= 0 && aggr <= cb_nodes);
            l = my_req[aggr].count;
            if (buf_idx != NULL && buf_is_contig) {
                buf_idx[aggr][l] = curr_idx;
                curr_idx += avail_len;
            }
            rem_len -= avail_len;

            my_req[aggr].offsets[l] = off;
            my_req[aggr].lens[l] = avail_len;
            my_req[aggr].count++;
        }
    }
    NCI_Free(aggr_ranks);
}

/* LUSTRE_Calc_others_req() calculates what requests from each of other
 * processes fall in this aggregator's file domain.
 *   IN: my_req[cb_nodes]: offset-length pairs of this rank's requests fall
 *       into each of aggregators
 *   OUT: count_others_req_per_proc[i]: number of noncontiguous requests of
 *        rank i that falls in this aggregator's file domain.
 *   OUT: others_req_ptr[nprocs]: requests of each of other ranks fall into
 *        this aggregator's file domain.
 */
static
void LUSTRE_Calc_others_req(PNCIO_File          *fd,
                            const PNCIO_Access  *my_req,
                            PNCIO_Access       **others_req_ptr)
{
    int i, myrank, nprocs, do_alltoallv;
    MPI_Count *count_my_req_per_proc, *count_others_req_per_proc;
    PNCIO_Access *others_req;
    size_t npairs, alloc_sz, pair_sz;

    /* first find out how much to send/recv and from/to whom */

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    others_req = (PNCIO_Access *) NCI_Malloc(nprocs * sizeof(PNCIO_Access));
    *others_req_ptr = others_req;

    /* Use my_req[i].count (the number of noncontiguous requests fall in
     * aggregator i's file domain) to set count_others_req_per_proc[j] (the
     * number of noncontiguous requests from process j fall into this
     * aggregator's file domain).
     */
    count_my_req_per_proc = (MPI_Count *) NCI_Calloc(nprocs * 2, sizeof(MPI_Count));
    count_others_req_per_proc = count_my_req_per_proc + nprocs;
    for (i=0; i<fd->hints->cb_nodes; i++)
        count_my_req_per_proc[fd->hints->ranklist[i]] = my_req[i].count;

    MPI_Alltoall(count_my_req_per_proc, 1, MPI_COUNT,
                 count_others_req_per_proc, 1, MPI_COUNT, fd->comm);

    /* calculate total number of offset-length pairs to be handled by this
     * aggregator, only aggregators will have non-zero number of pairs.
     */
    npairs = 0;
    for (i=0; i<nprocs; i++) {
        npairs += count_others_req_per_proc[i];
        others_req[i].count = count_others_req_per_proc[i];
        others_req[i].curr = 0;
    }
    NCI_Free(count_my_req_per_proc);

    /* The best communication approach for aggregators to collect offset-length
     * pairs from the non-aggregators is to allocate a single contiguous memory
     * space for my_req[] to store all its pairs of offsets and lens. The same
     * for others_req[].
     */
#ifdef HAVE_MPI_LARGE_COUNT
    pair_sz = sizeof(MPI_Offset) * 2;
    alloc_sz = pair_sz + sizeof(MPI_Count);
    others_req[0].offsets  = (MPI_Offset*) NCI_Malloc(npairs * alloc_sz);
    others_req[0].lens     = others_req[0].offsets + others_req[0].count;
    others_req[0].mem_ptrs = (MPI_Count*) (others_req[0].offsets + npairs * 2);
    for (i=1; i<nprocs; i++) {
        others_req[i].offsets  = others_req[i-1].offsets + others_req[i-1].count * 2;
        others_req[i].lens     = others_req[i].offsets + others_req[i].count;
        others_req[i].mem_ptrs = others_req[i-1].mem_ptrs + others_req[i-1].count;
    }
#else
    pair_sz = sizeof(MPI_Offset) + sizeof(int);
    alloc_sz = pair_sz + sizeof(MPI_Aint);
    others_req[0].offsets  = (MPI_Offset*) NCI_Malloc(npairs * alloc_sz);
    others_req[0].lens     = (int*) (others_req[0].offsets + others_req[0].count);
    char *ptr = (char*) others_req[0].offsets + pair_sz * npairs;
    others_req[0].mem_ptrs = (MPI_Aint*)ptr;

    ptr = (char*) others_req[0].offsets + pair_sz * others_req[0].count;
    for (i=1; i<nprocs; i++) {
        others_req[i].offsets = (MPI_Offset*)ptr;
        ptr += sizeof(MPI_Offset) * others_req[i].count;
        others_req[i].lens = (int*)ptr;
        ptr += sizeof(int) * others_req[i].count;
        others_req[i].mem_ptrs = others_req[i-1].mem_ptrs + others_req[i-1].count;
    }
#endif

    /* now send the calculated offsets and lengths to respective processes */

#ifdef CONSIDER_ALLTOALLV
    /* On Perlmutter at NERSC, when the number of processes per compute node is
     * large, using MPI_Alltoallv() instead of MPI_Isend/Irecv may avoid
     * possible hanging.  When hanging occurs, the error messages are
     * 1. RXC (0x11291:0) PtlTE 397:[Fatal] OVERFLOW buffer list exhausted
     * 2. MPICH WARNING: OFI is failing to make progress on posting a receive.
     *    MPICH suspects a hang due to completion queue exhaustion. Setting
     *    environment variable FI_CXI_DEFAULT_CQ_SIZE to a higher number might
     *    circumvent this scenario. OFI retry continuing...
     *
     * Below use a threshold of 48, number of processes per compute node.
     */
    do_alltoallv = (fd->num_nodes > 0) ? (nprocs / fd->num_nodes > 48) : 0;
#else
    do_alltoallv=0;
#endif

    if (do_alltoallv) {
        MPI_Offset *r_off_buf=NULL, *s_off_buf=NULL;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count *sendCounts, *recvCounts;
        MPI_Aint *sdispls, *rdispls;
        alloc_sz   = sizeof(MPI_Count) * 2 + sizeof(MPI_Aint) * 2;
        sendCounts = (MPI_Count*) NCI_Calloc(nprocs, alloc_sz);
        recvCounts = sendCounts + nprocs;
        sdispls    = (MPI_Aint*) (recvCounts + nprocs);
        rdispls    = sdispls + nprocs;
#else
        int *sendCounts, *recvCounts, *sdispls, *rdispls;
        alloc_sz   = sizeof(int) * 4;
        sendCounts = (int*) NCI_Calloc(nprocs, alloc_sz);
        recvCounts = sendCounts + nprocs;
        sdispls    = recvCounts + nprocs;
        rdispls    = sdispls + nprocs;
#endif

        /* prepare receive side */
        r_off_buf = others_req[0].offsets;
        for (i=0; i<nprocs; i++) {
            recvCounts[i] = others_req[i].count * pair_sz;
            /* Note all others_req[*].offsets are allocated in a single malloc(). */
            rdispls[i] = (char*)others_req[i].offsets - (char*)r_off_buf;
        }

        /* prepare send side */
        s_off_buf = my_req[0].offsets;
        for (i=0; i<fd->hints->cb_nodes; i++) {
            int dest = fd->hints->ranklist[i];
            sendCounts[dest] = my_req[i].count * pair_sz;
            /* Note all my_req[*].offsets are allocated in a single malloc(). */
            sdispls[dest] = (char*)my_req[i].offsets - (char*)s_off_buf;
        }

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Alltoallv_c(s_off_buf, sendCounts, sdispls, MPI_BYTE,
                        r_off_buf, recvCounts, rdispls, MPI_BYTE, fd->comm);
#else
        MPI_Alltoallv(s_off_buf, sendCounts, sdispls, MPI_BYTE,
                      r_off_buf, recvCounts, rdispls, MPI_BYTE, fd->comm);
#endif

        NCI_Free(sendCounts);
    }
    else { /* instead of using alltoall, use MPI_Issend and MPI_Irecv */
        int nreqs;
        MPI_Request *requests = (MPI_Request *)
            NCI_Malloc((nprocs + fd->hints->cb_nodes) * sizeof(MPI_Request));

        nreqs = 0;
        for (i = 0; i < nprocs; i++) {
            if (others_req[i].count == 0) /* nothing to receive from rank i */
                continue;

            /* Note the memory address of others_req[i].lens is right after
             * others_req[i].offsets. This allows the following recv call to
             * receive both offsets and lens in a single call.
             */
            if (i == myrank) {
                /* send to self uses memcpy(), here
                 * others_req[i].count == my_req[fd->my_cb_nodes_index].count
                 */
                memcpy(others_req[i].offsets,
                       my_req[fd->my_cb_nodes_index].offsets,
                       my_req[fd->my_cb_nodes_index].count * pair_sz);
            }
            else {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(others_req[i].offsets, others_req[i].count*pair_sz,
                          MPI_BYTE, i, 0, fd->comm, &requests[nreqs++]);
#else
                MPI_Irecv(others_req[i].offsets, others_req[i].count*pair_sz,
                          MPI_BYTE, i, 0, fd->comm, &requests[nreqs++]);
#endif
            }
        }

#ifdef WKL_DEBUG
/* WRF hangs below when calling MPI_Waitall(), at running 16 nodes, 128 ranks
 * per node on Perlmutter, when these 3 env variables are set:
 *    FI_UNIVERSE_SIZE        = 2048
 *    FI_CXI_DEFAULT_CQ_SIZE  = 524288
 *    FI_CXI_RX_MATCH_MODE    = software
 *
 * Using MPI_Alltoallv seems to be able to avoid such hanging problem. (above)
 */
// MPI_Barrier(fd->comm); /* This barrier prevents the MPI_Waitall below from hanging !!! */
#endif

        for (i=0; i<fd->hints->cb_nodes; i++) {
            if (my_req[i].count == 0 || i == fd->my_cb_nodes_index)
                continue; /* nothing to send or send to self */

            /* Note the memory address of my_req[i].lens is right after
             * my_req[i].offsets. This allows the following Issend call to
             * send both offsets and lens in a single call.
             */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Issend_c(my_req[i].offsets, my_req[i].count * pair_sz, MPI_BYTE,
                       fd->hints->ranklist[i], 0, fd->comm, &requests[nreqs++]);
#else
            MPI_Issend(my_req[i].offsets, my_req[i].count * pair_sz, MPI_BYTE,
                       fd->hints->ranklist[i], 0, fd->comm, &requests[nreqs++]);
#endif
        }

        if (nreqs) {
#ifdef HAVE_MPI_STATUSES_IGNORE
            MPI_Waitall(nreqs, requests, MPI_STATUSES_IGNORE);
#else
            MPI_Status *statuses = (MPI_Status *)
                                   NCI_Malloc(nreqs * sizeof(MPI_Status));
            MPI_Waitall(nreqs, requests, statuses);
            NCI_Free(statuses);
#endif
        }
        NCI_Free(requests);
    }
}

MPI_Offset PNCIO_LUSTRE_WriteStridedColl(PNCIO_File *fd,
                                         const void *buf,
                                         PNCIO_View  buf_view,
                                         MPI_Offset  offset)
{
    /* Uses a generalized version of the extended two-phase method described in
     * "An Extended Two-Phase Method for Accessing Sections of Out-of-Core
     * Arrays", Rajeev Thakur and Alok Choudhary, Scientific Programming,
     * (5)4:301--317, Winter 1996.
     * http://www.mcs.anl.gov/home/thakur/ext2ph.ps
     */

    int i, j, nprocs, myrank;
    int do_collect = 1, do_ex_wr;
    MPI_Offset start_offset, end_offset;
    MPI_Offset min_st_loc = -1, max_end_loc = -1;
    MPI_Offset w_len=0;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

// printf("%s %d: offset=%lld\n",__func__,__LINE__,offset);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
MPI_Barrier(fd->comm);
double curT = MPI_Wtime();
#endif

    /* fd->flat_file contains a list of starting file offsets and lengths of
     * write requests made by this rank. Similarly, buf_view contains a list of
     * offset-length pairs describing the write buffer layout.  Note as PnetCDF
     * never re-uses a fileview or buffer view.
     *
     * Note that MPI standard (MPI 3.1 Chapter 13.1.1 and MPI 4.0 Chapter
     * 14.1.1) requires that the typemap displacements of etype and filetype
     * set by the user are non-negative and monotonically non-decreasing. This
     * makes fd->flat_file.off[] to be monotonically non-decreasing.
     *
     * This rank's aggregate file access region is from start_offset to
     * end_offset. Note: end_offset points to the last byte-offset to be
     * accessed. E.g., if start_offset=0 and end_offset=99, then the aggregate
     * file access region is of size 100 bytes. If this rank has no data to
     * write, end_offset == (start_offset - 1)
     */
    MPI_Offset one_off;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset one_len;
#else
    int one_len;
#endif

    if (fd->flat_file.count == 0) { /* TODO: is fd->flat_file.count == 0? */
        /* whole file is visible */
        start_offset = offset;
        end_offset = offset + buf_view.size - 1;
        if (buf_view.size > 0) { /* no-zero sized request */
            /* setting fd->flat_file is necessary for constructing my_req */
            one_off = offset;
            one_len = buf_view.size;
            fd->flat_file.off = &one_off;
            fd->flat_file.len = &one_len;
            fd->flat_file.count = 1;
        }
    }
    else {
        start_offset = offset + fd->flat_file.off[0];
        end_offset   = fd->flat_file.off[fd->flat_file.count-1]
                     + fd->flat_file.len[fd->flat_file.count-1] - 1;
    }
/*
    else if (fd->flat_file.count > 0) {
        start_offset = offset + fd->flat_file.off[0];
        end_offset   = fd->flat_file.off[fd->flat_file.count-1]
                     + fd->flat_file.len[fd->flat_file.count-1] - 1;
    }
    else {
        start_offset = offset;
        end_offset   = offset + fd->flat_file.size - 1;
    }
*/
// if (myrank==0) printf("%s %d: fd->flat_file size=%lld count=%lld offset=%lld start_offset=%lld end_offset=%lld\n",__func__,__LINE__, fd->flat_file.size, fd->flat_file.count,offset,start_offset,end_offset);

    buf_view.idx  = 0;
    buf_view.rem = buf_view.size;
    if (buf_view.count > 1)
        buf_view.rem = buf_view.len[0];

    if (fd->hints->romio_cb_write == PNCIO_HINT_DISABLE) {
        /* collective write is explicitly disabled by user */
        do_collect = 0;
    }
    else {
        /* Calculate the aggregate access region of all ranks and check if
         * write requests are interleaved among all ranks.
         */
        int is_interleaved, large_indv_req = 1;
        MPI_Offset striping_range, st_end[2], *st_end_all = NULL;

        /* Gather starting and ending file offsets of write requests from all
         * ranks into st_end_all[]. Even indices of st_end_all[] are starting
         * offsets, and odd indices are ending offsets.
         */
        st_end[0] = start_offset;
        st_end[1] = end_offset;
        st_end_all = (MPI_Offset *) NCI_Malloc(nprocs * 2 * sizeof(MPI_Offset));
        MPI_Allgather(st_end, 2, MPI_OFFSET, st_end_all, 2, MPI_OFFSET, fd->comm);

        /* The loop below does the followings.
         * 1. Calculate this rank's aggregate access region.
         * 2. Check whether or not the requests are interleaved among all ranks.
         * 3. Check whether there are LARGE individual requests. Here, "large"
         *    means a write range is > (striping_factor * striping_unit). In
         *    this case, independent write will perform faster than collective.
         */
        striping_range = fd->hints->striping_unit * fd->hints->striping_factor;
        is_interleaved = 0;
        for (i = 0; i < nprocs * 2; i += 2) {
            if (st_end_all[i] > st_end_all[i + 1]) {
                /* process rank (i/2) has no data to write */
                continue;
            }
            min_st_loc = st_end_all[i];
            max_end_loc = st_end_all[i + 1];
            if (st_end_all[i+1] - st_end_all[i] < striping_range)
                large_indv_req = 0;
            j = i; /* j is the rank of making first non-zero request */
            i += 2;
            break;
        }
        for (; i < nprocs * 2; i += 2) {
            if (st_end_all[i] > st_end_all[i + 1]) {
                /* process rank (i/2) has no data to write */
                continue;
            }
            if (st_end_all[i] < st_end_all[j+1]) {
                /* start offset of process rank (i/2) is less than the end
                 * offset of process rank (i/2-1)
                 */
                is_interleaved = 1;
            }
            min_st_loc = MIN(st_end_all[i], min_st_loc);
            max_end_loc = MAX(st_end_all[i + 1], max_end_loc);
            if (st_end_all[i+1] - st_end_all[i] < striping_range)
                large_indv_req = 0;
            j = i;
        }
        NCI_Free(st_end_all);

// if (myrank==0) printf("%s %d: do_collect=%d is_interleaved=%d buf_view size=%lld count=%lld is_contig=%d start_offset=%lld end_offset=%lld\n",__func__,__LINE__, do_collect,is_interleaved,buf_view.size,buf_view.count,buf_view.is_contig, start_offset,end_offset);
        if (fd->hints->romio_cb_write == PNCIO_HINT_ENABLE) {
            /* explicitly enabled by user */
            do_collect = 1;
        }
        else if (fd->hints->romio_cb_write == PNCIO_HINT_AUTO) {
// if (myrank==0) printf("%s %d: large_indv_req=%d cb_nodes=%d striping_factor=%d\n",__func__,__LINE__, large_indv_req,fd->hints->cb_nodes , fd->hints->striping_factor);
            /* Check if collective write is actually necessary, only when
             * romio_cb_write hint is set to PNCIO_HINT_AUTO.
             *
             * Two typical access patterns can benefit from collective write.
             *   1) access file regions of all processes are interleaved, and
             *   2) the individual request sizes are not too big, i.e. no
             *      bigger than striping_range. Large individual requests may
             *      result in a high communication cost in order to
             *      redistribute requests from non-aggregators to I/O
             *      aggregators.
             */
            if (nprocs == 1)
                do_collect = 0;
            else if (!is_interleaved && large_indv_req &&
                     fd->hints->cb_nodes <= fd->hints->striping_factor) {
                /* do independent write, if every rank's write range >
                 * striping_range and writes are not interleaved in file
                 * space
                 */
                do_collect = 0;
            }
        }
    }

    /* If collective I/O is determined not necessary, use independent I/O */
    if (!do_collect) {

        if (buf_view.size == 0) /* zero-sized request */
            return 0;

        if (fd->flat_file.is_contig && buf_view.is_contig) {
            /* both buffer and fileview are contiguous */
            if (fd->flat_file.count > 0) offset += fd->flat_file.off[0];
#ifdef WKL_DEBUG
            printf("%s %d: SWITCH to PNCIO_WriteContig !!!\n",__func__,__LINE__);
#endif

            return PNCIO_WriteContig(fd, buf, buf_view.size, offset);
        }

#ifdef WKL_DEBUG
        printf("%s %d: SWITCH to PNCIO_LUSTRE_WriteStrided !!!\n",
                   __func__,__LINE__);
#endif

        return PNCIO_LUSTRE_WriteStrided(fd, buf, buf_view, offset);
    }

    /* Now we are using collective I/O (two-phase I/O strategy) */

#ifdef ADJUST_STRIPING_UNIT
    /* adjust striping_unit when striping_factor is twice or more than the
     * number of compute nodes. Note cb_node is set to at least
     * striping_factor, if nprocs >= striping_factor. Adjustment below is to
     * let each aggregator to write to two or more consecutive OSTs, which can
     * most likely improve the performance. This will still yield an effect of
     * any one OST receiving write requests from aggregators running on only
     * one compute node.
     */
    int orig_striping_unit = fd->hints->striping_unit;

    if (fd->hints->striping_factor >= fd->num_nodes * 2) {
        fd->hints->striping_unit *= (fd->hints->striping_factor / fd->num_nodes);

        if (fd->hints->cb_buffer_size < fd->hints->striping_unit) {
            char value[MPI_MAX_INFO_VAL + 1];

            fd->hints->cb_buffer_size = fd->hints->striping_unit;
            sprintf(value, "%d", fd->hints->cb_buffer_size);
            MPI_Info_set(fd->info, "cb_buffer_size", value);
            if (fd->is_agg) {
                NCI_Free(fd->io_buf);
                fd->io_buf = (void*) NCI_Calloc(1, fd->hints->cb_buffer_size);
            }
        }
#ifdef WKL_DEBUG
        if (myrank == 0)
            printf("Warning: %s line %d: Change striping_unit from %d to %d\n",
                   __func__, __LINE__, orig_striping_unit, fd->hints->striping_unit);
#endif
    }
#endif

    /* my_req[cb_nodes] is an array of access info, one for each I/O aggregator
     * whose file domain has this rank's request.
     */
    PNCIO_Access *my_req;

    /* others_req[nprocs] is an array of access info, one for each ranks, both
     * aggregators and non-aggregators, whose write requests fall into this
     * aggregator's file domain. others_req[] matters only for aggregators.
     */
    PNCIO_Access *others_req;
    MPI_Offset **buf_idx = NULL;

    if (buf_view.is_contig)
        buf_idx = (MPI_Offset **) NCI_Malloc(fd->hints->cb_nodes *
                                                sizeof(MPI_Offset*));

    /* Calculate the portions of this rank's write requests that fall into the
     * file domains of each I/O aggregator. No inter-process communication is
     * performed in LUSTRE_Calc_my_req().
     */
    LUSTRE_Calc_my_req(fd, buf_view.is_contig, &my_req, buf_idx);

    if (fd->hints->romio_ds_write != PNCIO_HINT_DISABLE) {
        /* When data sieving is considered, below check the current file size
         * first. If the aggregate access region of this collective write is
         * beyond the current file size, then we can safely skip the read of
         * the read-modify-write of data sieving.
         */
        if (fd->is_agg) {
            /* Obtain the current file size. Note an MPI_Allgather() has been
             * called above to calculate the aggregate access region. Thus all
             * prior independent I/O should have completed by now, so it is
             * safe to call lseek() to query the file size.
             */
            MPI_Offset cur_off, fsize;

            cur_off = lseek(fd->fd_sys, 0, SEEK_CUR);
            fsize   = lseek(fd->fd_sys, 0, SEEK_END);
            /* Ignore the error, and proceed as if file size is very large. */
#ifdef PNETCDF_DEBUG
            if (fsize == -1)
                fprintf(stderr, "%s at %d: lseek SEEK_END failed on file %s (%s)\n",
                        __func__,__LINE__, fd->filename, strerror(errno));
#endif
            fd->skip_read = (fsize >=0 && min_st_loc >= fsize);

            /* restore file pointer */
            lseek(fd->fd_sys, cur_off, SEEK_SET);
        }
    }
    else
        fd->skip_read = 1;

// if (fd->is_agg && !fd->skip_read) { MPI_Offset fsize = lseek(fd->fd_sys, 0, SEEK_END); printf("%d: %s at %d: skip_read=%d min_st_loc=%lld fsize=%lld\n",myrank,__func__,__LINE__,fd->skip_read,min_st_loc,fsize); }

    /* For aggregators, calculate the portions of all other ranks' requests
     * fall into this aggregator's file domain (note only I/O aggregators are
     * assigned file domains).
     *
     * Inter-process communication is required to construct others_req[],
     * including MPI_Alltoall, MPI_Issend, MPI_Irecv, and MPI_Waitall.
     */
    LUSTRE_Calc_others_req(fd, my_req, &others_req);

    /* Two-phase I/O: first communication phase to exchange write data from all
     * ranks to the I/O aggregators, followed by the write phase where only I/O
     * aggregators write to the file.
     *
     * Unless MPI_Alltoallw() is used (when use_alltoallw is set to 1), there
     * is no collective MPI communication beyond this point, as
     * LUSTRE_Exch_and_write() calls only MPI_Issend, MPI_Irecv, and
     * MPI_Waitall. Thus it is safe for those non-aggregators making zero-sized
     * request to skip the call.
     */

    /* if this rank has data to write, then participate exchange-and-write */
    do_ex_wr = (buf_view.size == 0) ? 0 : 1;
    use_alltoallw = 0;

#ifdef USE_MPI_ALLTOALLW
    {
        /* When num_nodes < striping_factor, using MPI_Alltoallw in
         * commit_comm_phase() is faster than MPI_Issend/MPI_Irecv ... ?
         */
        char *env_str;
        if ((env_str = getenv("PNETCDF_USE_ALLTOALLW")) != NULL)
            use_alltoallw = (strcasecmp(env_str, "true") == 0) ? 1: 0;
    }
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->write_timing[1] += MPI_Wtime() - curT;
#endif

    if (do_ex_wr || fd->is_agg)
        /* This rank participates exchange and write only when it has non-zero
         * data to write or is an I/O aggregator
         */
        w_len = LUSTRE_Exch_and_write(fd, buf, buf_view, others_req, my_req,
                                      min_st_loc, max_end_loc, buf_idx);

    /* free all memory allocated */
    NCI_Free(others_req[0].offsets);
    NCI_Free(others_req);

    if (buf_idx != NULL) {
        NCI_Free(buf_idx[0]);
        NCI_Free(buf_idx);
    }
    NCI_Free(my_req[0].offsets);
    NCI_Free(my_req);

#ifdef ADJUST_STRIPING_UNIT
    /* restore the original striping_unit */
    fd->hints->striping_unit = orig_striping_unit;
#endif

    /* If this collective write is followed by an independent write, it's
     * possible to have those subsequent writes on other processes race ahead
     * and sneak in before the read-modify-write completes.  We carry out a
     * collective communication at the end here so no one can start independent
     * I/O before collective I/O completes.
     *
     * need to do some gymnastics with the error codes so that if something
     * went wrong, all processes report error, but if a process has a more
     * specific error code, we can still have that process report the
     * additional information
     */
    /* optimization: if only one process performing I/O, we can perform
     * a less-expensive Bcast. */
    if (fd->hints->cb_nodes == 1)
        MPI_Bcast(&w_len, 1, MPI_OFFSET, fd->hints->ranklist[0], fd->comm);
    else
        MPI_Allreduce(MPI_IN_PLACE, &w_len, 1, MPI_OFFSET, MPI_MIN, fd->comm);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fd->is_agg) fd->write_timing[0] += MPI_Wtime() - curT;
#endif

    /* w_len may not be the same as buf_view.size, because data sieving may
     * write more than requested.
     */
    return buf_view.size;
}

static
void comm_phase_alltoallw(PNCIO_File    *fd,
                          disp_len_list *send_list,  /* [cb_nodes] */
                          disp_len_list *recv_list)  /* [nprocs] */
{
    /* This subroutine performs the sam communication tasks as the below
     * commit_comm_phase(), but using MPI_Alltoallw() instead of MPI_Issend and
     * MPI_Irecv.
     *
     * It creates a datatype combining all displacement-length
     * pairs in each element of send_list[]. The datatype is used when calling
     * MPI_Issend to send write data to the I/O aggregators. Similarly, it
     * creates a datatype combining all displacement-length pairs in each
     * element of recv_list[] and uses it when calling MPI_Irecv or MPI_Recv
     * to receive write data from all processes.
     */
    int i, nprocs, rank;
    size_t alloc_sz;
    MPI_Datatype *sendTypes, *recvTypes;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &rank);

    /* calculate send/recv derived types metadata */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *sendCounts, *recvCounts;
    MPI_Aint *sdispls, *rdispls;
    alloc_sz = sizeof(MPI_Count) + sizeof(MPI_Aint);
    sendCounts = (MPI_Count*) NCI_Calloc(nprocs * 2, alloc_sz);
    sdispls = (MPI_Aint*) (sendCounts + (nprocs * 2));
#else
    int *sendCounts, *recvCounts, *sdispls, *rdispls;
    alloc_sz = sizeof(int) * 2;
    sendCounts = (int*) NCI_Calloc(nprocs * 2, alloc_sz);
    sdispls = (int*) (sendCounts + (nprocs * 2));
#endif
    recvCounts = sendCounts + nprocs;
    rdispls = sdispls + nprocs;

    /* allocate send/recv derived type arrays */
    sendTypes = (MPI_Datatype*)NCI_Malloc(sizeof(MPI_Datatype) * nprocs * 2);
    recvTypes = sendTypes + nprocs;

    for (i=0; i<nprocs; i++)
        sendTypes[i] = recvTypes[i] = MPI_BYTE;

    /* prepare receive side: construct recv derived data types */
    if (fd->is_agg && recv_list != NULL) {
        for (i=0; i<nprocs; i++) {
            /* check if nothing to receive or if self */
            if (recv_list[i].count == 0 || i == rank) continue;

            recvCounts[i] = 1;

            /* combine reqs using new datatype */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(recv_list[i].count, recv_list[i].len,
                                       recv_list[i].disp, MPI_BYTE,
                                       &recvTypes[i]);
#else
            MPI_Type_create_hindexed(recv_list[i].count, recv_list[i].len,
                                     recv_list[i].disp, MPI_BYTE,
                                     &recvTypes[i]);
#endif
            MPI_Type_commit(&recvTypes[i]);
        }
    }

    /* prepare send side: construct send derived data types */
    for (i=0; i<fd->hints->cb_nodes; i++) {
        /* check if nothing to send or if self */
        if (send_list[i].count == 0 || i == fd->my_cb_nodes_index) continue;

        int dest = fd->hints->ranklist[i];
        sendCounts[dest] = 1;

        /* combine reqs using new datatype */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Type_create_hindexed_c(send_list[i].count, send_list[i].len,
                                   send_list[i].disp, MPI_BYTE,
                                   &sendTypes[dest]);
#else
        MPI_Type_create_hindexed(send_list[i].count, send_list[i].len,
                                 send_list[i].disp, MPI_BYTE,
                                 &sendTypes[dest]);
#endif
        MPI_Type_commit(&sendTypes[dest]);
    }

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Alltoallw_c(MPI_BOTTOM, sendCounts, sdispls, sendTypes,
                    MPI_BOTTOM, recvCounts, rdispls, recvTypes, fd->comm);
#else
    MPI_Alltoallw(MPI_BOTTOM, sendCounts, sdispls, sendTypes,
                  MPI_BOTTOM, recvCounts, rdispls, recvTypes, fd->comm);
#endif

    for (i=0; i<nprocs; i++) {
        if (sendTypes[i] != MPI_BYTE)
            MPI_Type_free(&sendTypes[i]);
        if (recvTypes[i] != MPI_BYTE)
            MPI_Type_free(&recvTypes[i]);
    }
    NCI_Free(sendCounts);
    NCI_Free(sendTypes);

    /* clear send_list and recv_list for future reuse */
    for (i = 0; i < fd->hints->cb_nodes; i++)
        send_list[i].count = 0;

    if (recv_list != NULL)
        for (i = 0; i < nprocs; i++)
            recv_list[i].count = 0;
}

static
void commit_comm_phase(PNCIO_File    *fd,
                       disp_len_list *send_list,  /* [cb_nodes] */
                       disp_len_list *recv_list)  /* [nprocs] */
{
    /* This subroutine creates a datatype combining all displacement-length
     * pairs in each element of send_list[]. The datatype is used when calling
     * MPI_Issend to send write data to the I/O aggregators. Similarly, it
     * creates a datatype combining all displacement-length pairs in each
     * element of recv_list[] and uses it when calling MPI_Irecv or MPI_Recv
     * to receive write data from all processes.
     */
    int i, nprocs, rank, nreqs;
    MPI_Request *reqs;
    MPI_Datatype sendType, recvType;
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    int j;
    double dtype_time=MPI_Wtime();
#endif

    if (use_alltoallw)
        return comm_phase_alltoallw(fd, send_list, recv_list);

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &rank);

    nreqs = fd->hints->cb_nodes;
    nreqs += (fd->is_agg) ? nprocs : 0;
    reqs = (MPI_Request *)NCI_Malloc(sizeof(MPI_Request) * nreqs);
    nreqs = 0;

    /* receiving part */
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    /* recv buffer type profiling */
    int nrecvs=0;
    MPI_Offset max_r_amnt=0, max_r_count=0;
#endif

    if (fd->is_agg && recv_list != NULL) {
        for (i = 0; i < nprocs; i++) {
            /* check if nothing to receive or if self */
            if (recv_list[i].count == 0 || i == rank) continue;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
            MPI_Offset r_amnt=0;
            for (j=0; j<recv_list[i].count; j++)
                r_amnt += recv_list[i].len[j];
            max_r_amnt = MAX(max_r_amnt, r_amnt);
            max_r_count = MAX(max_r_count, recv_list[i].count);
            nrecvs++;
#endif

            /* combine reqs using new datatype */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(recv_list[i].count, recv_list[i].len,
                                       recv_list[i].disp, MPI_BYTE,
                                       &recvType);
#else
            MPI_Type_create_hindexed(recv_list[i].count, recv_list[i].len,
                                     recv_list[i].disp, MPI_BYTE,
                                     &recvType);
#endif
            MPI_Type_commit(&recvType);

            if (fd->atomicity) { /* Blocking Recv */
                MPI_Status status;
                MPI_Recv(MPI_BOTTOM, 1, recvType, i, 0, fd->comm, &status);
            }
            else
                MPI_Irecv(MPI_BOTTOM, 1, recvType, i, 0, fd->comm,
                          &reqs[nreqs++]);
            MPI_Type_free(&recvType);
        }
    }

    /* send reqs */
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    /* send buffer type profiling */
    int nsends=0;
    MPI_Offset max_s_amnt=0, max_s_count=0;
#endif

    for (i = 0; i < fd->hints->cb_nodes; i++) {
        /* check if nothing to send or if self */
        if (send_list[i].count == 0 || i == fd->my_cb_nodes_index) continue;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        MPI_Offset s_amnt=0;
        for (j=0; j<send_list[i].count; j++)
            s_amnt += send_list[i].len[j];
        max_s_amnt = MAX(max_s_amnt, s_amnt);
        max_s_count = MAX(max_s_count, send_list[i].count);
        nsends++;
#endif

        /* combine reqs using new datatype */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Type_create_hindexed_c(send_list[i].count, send_list[i].len,
                                   send_list[i].disp, MPI_BYTE, &sendType);
#else
        MPI_Type_create_hindexed(send_list[i].count, send_list[i].len,
                                 send_list[i].disp, MPI_BYTE, &sendType);
#endif
        MPI_Type_commit(&sendType);

        MPI_Issend(MPI_BOTTOM, 1, sendType, fd->hints->ranklist[i], 0,
                   fd->comm, &reqs[nreqs++]);
        MPI_Type_free(&sendType);
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->write_timing[4] += MPI_Wtime() - dtype_time;

/*
    fd->write_counter[2] = MAX(fd->write_counter[2], nsends);
    fd->write_counter[3] = MAX(fd->write_counter[3], nrecvs);
    fd->write_counter[4] = MAX(fd->write_counter[4], max_r_amnt);
    fd->write_counter[5] = MAX(fd->write_counter[5], max_s_amnt);
    fd->write_counter[6] = MAX(fd->write_counter[6], max_r_count);
    fd->write_counter[7] = MAX(fd->write_counter[7], max_s_count);
*/
#endif

    if (nreqs > 0) {
#ifdef HAVE_MPI_STATUSES_IGNORE
        MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *)
                               NCI_Malloc(nreqs * sizeof(MPI_Status));
        MPI_Waitall(nreqs, reqs, statuses);
        NCI_Free(statuses);
#endif
    }

    NCI_Free(reqs);

    /* clear send_list and recv_list for future reuse */
    for (i = 0; i < fd->hints->cb_nodes; i++)
        send_list[i].count = 0;

    if (recv_list != NULL)
        for (i = 0; i < nprocs; i++)
            recv_list[i].count = 0;
}

/*----< LUSTRE_Exch_and_write() >--------------------------------------------*/
/* Each process sends all its write requests to I/O aggregators based on the
 * file domain assignment to the aggregators. In this implementation, a file is
 * first divided into stripes which are assigned to the aggregators in a
 * round-robin fashion. The "exchange" of write data from non-aggregators to
 * aggregators is carried out in 'ntimes' rounds. Each round covers an
 * aggregate file region of size equal to the file stripe size times the number
 * of I/O aggregators. The file writes are carried out in every 'nbufs'
 * iterations, where 'nbufs' == cb_buffer_size / file stripe size. This approach
 * is different from ROMIO's implementation as in MPICH 4.2.3.
 *
 * Other implementations developers are referring to the paper: Wei-keng Liao,
 * and Alok Choudhary. "Dynamically Adapting File Domain Partitioning Methods
 * for Collective I/O Based on Underlying Parallel File System Locking
 * Protocols", in The Supercomputing Conference, 2008.
 */
static
MPI_Offset LUSTRE_Exch_and_write(PNCIO_File    *fd,
                                 const void    *buf,
                                 PNCIO_View     buf_view,
                                 PNCIO_Access  *others_req,
                                 PNCIO_Access  *my_req,
                                 MPI_Offset     min_st_loc,
                                 MPI_Offset     max_end_loc,
                                 MPI_Offset   **buf_idx)
{
    char **write_buf = NULL, **recv_buf = NULL, **send_buf = NULL;
    size_t alloc_sz;
    int nprocs, myrank, nbufs, ibuf, batch_idx=0, cb_nodes, striping_unit;
    MPI_Count i, j, m, ntimes;
    MPI_Count **recv_size=NULL, **recv_count=NULL;
    MPI_Count **recv_start_pos=NULL, *send_size;
    MPI_Offset end_loc, req_off, iter_end_off, *off_list, step_size;
    MPI_Offset *this_buf_idx=NULL;
    off_len_list *srt_off_len = NULL;
    disp_len_list *send_list = NULL, *recv_list = NULL;
    MPI_Offset w_len, total_w_len=0;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    cb_nodes = fd->hints->cb_nodes;
    striping_unit = fd->hints->striping_unit;

    /* The aggregate access region (across all processes) of this collective
     * write starts from min_st_loc and ends at max_end_loc. The collective
     * write is carried out in 'ntimes' rounds of two-phase I/O. Each round
     * covers an aggregate file region of size 'step_size' written only by
     * cb_nodes number of I/O aggregators. Note non-aggregators must also
     * participate all ntimes rounds to send their requests to I/O aggregators.
     *
     * step_size = the number of I/O aggregators x striping_unit
     *
     * Note the number of write phases = ntimes / nbufs, as writes (and
     * communication) are accumulated for nbufs rounds before flushed.
     */
    step_size = (MPI_Offset)cb_nodes * striping_unit;

    /* align min_st_loc downward to the nearest file stripe boundary */
    min_st_loc -= min_st_loc % (MPI_Offset) striping_unit;

    /* ntimes is the number of rounds of two-phase I/O */
    ntimes = (max_end_loc - min_st_loc + 1) / step_size;
    if ((max_end_loc - min_st_loc + 1) % step_size)
        ntimes++;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->write_counter[0] = MAX(fd->write_counter[0], ntimes);
#endif

    /* collective buffer is divided into 'nbufs' sub-buffers. Each sub-buffer
     * is of size equal to Lustre stripe size. Write data of non-aggregators
     * are sent to aggregators and stored in aggregators' sub-buffers, one for
     * each round. All nbufs sub-buffers are altogether flushed to file every
     * nbufs rounds.
     *
     * fd->hints->cb_buffer_size, collective buffer size, for Lustre must be at
     * least striping_unit. This requirement has been checked at the file
     * open/create time when fd->io_buf is allocated.
     *
     * Note cb_buffer_size and striping_unit may also be adjusted earlier in
     * PNCIO_LUSTRE_WriteStridedColl().
     */
    nbufs = fd->hints->cb_buffer_size / striping_unit;
    assert(nbufs > 0); /* must at least 1 */

    /* in case number of rounds is less than nbufs */
    nbufs = (ntimes < nbufs) ? (int)ntimes : nbufs;

    /* off_list[m] is the starting file offset of this aggregator's write
     *     region in iteration m (file domain of iteration m). This offset
     *     may not be aligned with file stripe boundaries.
     * end_loc is the ending file offset of this aggregator's file domain.
     */
    off_list = (MPI_Offset *) NCI_Malloc(ntimes * sizeof(MPI_Offset));
    end_loc = -1;
    for (m = 0; m < ntimes; m++)
        off_list[m] = max_end_loc;
    for (i = 0; i < nprocs; i++) {
// if (myrank == 0) printf("%s at %d: others_req[%d] count=%lld\n",__func__,__LINE__, i,others_req[i].count);
        for (j = 0; j < others_req[i].count; j++) {
            req_off = others_req[i].offsets[j];
            m = (int) ((req_off - min_st_loc) / step_size);
            off_list[m] = MIN(off_list[m], req_off);
            end_loc = MAX(end_loc, (others_req[i].offsets[j] + others_req[i].lens[j] - 1));
        }
    }
// if (myrank == 0) printf("%s at %d: end_loc=%lld nbufs=%d recv_list=%s\n",__func__,__LINE__, end_loc,nbufs,(recv_list==NULL)?"NULL":"NOT NULL");

    /* Allocate displacement-length pair arrays, describing the send buffer.
     * send_list[i].count: number displacement-length pairs.
     * send_list[i].len: length in bytes.
     * send_list[i].disp: displacement (send buffer address).
     */
    send_list = (disp_len_list*) NCI_Malloc(sizeof(disp_len_list) * cb_nodes);
    for (i = 0; i < cb_nodes; i++) {
        send_list[i].count = 0;
#ifdef HAVE_MPI_LARGE_COUNT
        alloc_sz = sizeof(MPI_Count) * 2;
        send_list[i].disp = (MPI_Count*) NCI_Malloc(alloc_sz * nbufs);
        send_list[i].len  = send_list[i].disp + nbufs;
#else
        alloc_sz = sizeof(MPI_Aint) + sizeof(int);
        send_list[i].disp = (MPI_Aint*) NCI_Malloc(alloc_sz * nbufs);
        send_list[i].len  = (int*) (send_list[i].disp + nbufs);
#endif
    }

    /* end_loc >= 0 indicates this process has something to write to the file.
     * Only I/O aggregators can have end_loc > 0. write_buf is the collective
     * buffer and only matter for I/O aggregators. recv_buf is the buffer used
     * only by aggregators to receive requests from non-aggregators. Its size
     * may be larger then the file stripe size, in case when writes from
     * non-aggregators overlap. In this case, it will be realloc-ed in
     * LUSTRE_W_Exchange_data(). The received data is later copied over to
     * write_buf, whose contents will be written to file.
     */
    if (end_loc >= 0 && nbufs > 0) {
        /* Allocate displacement-length pair arrays, describing the recv buffer.
         * recv_list[i].count: number displacement-length pairs.
         * recv_list[i].len: length in bytes.
         * recv_list[i].disp: displacement (recv buffer address).
         */
        assert(fd->is_agg);

        recv_list = (disp_len_list*) NCI_Malloc(sizeof(disp_len_list) * nprocs);
        for (i = 0; i < nprocs; i++) {
            recv_list[i].count = 0;
#ifdef HAVE_MPI_LARGE_COUNT
            alloc_sz = sizeof(MPI_Count) * 2;
            recv_list[i].disp = (MPI_Count*) NCI_Malloc(alloc_sz * nbufs);
            recv_list[i].len  = recv_list[i].disp + nbufs;
#else
            alloc_sz = sizeof(MPI_Aint) + sizeof(int);
            recv_list[i].disp = (MPI_Aint*) NCI_Malloc(alloc_sz * nbufs);
            recv_list[i].len  = (int*) (recv_list[i].disp + nbufs);
#endif
        }

        /* collective buffer was allocated at file open/create. For Lustre, its
         * size must be at least striping_unit, which has been checked at the
         * time fd->io_buf is allocated.
         */
        assert(fd->io_buf != NULL);

        /* divide collective buffer into nbufs sub-buffers */
        write_buf = (char **) NCI_Malloc(nbufs * sizeof(char*));
        write_buf[0] = fd->io_buf;

        /* Similarly, receive buffer consists of nbufs sub-buffers */
        recv_buf = (char **) NCI_Malloc(nbufs * sizeof(char*));
        recv_buf[0] = (char *) NCI_Malloc(striping_unit);

        /* recv_count[j][i] is the number of off-len pairs to be received from
         * each proc i in round j
         */
        recv_count    = (MPI_Count**) NCI_Malloc(3 * nbufs * sizeof(MPI_Count*));
        recv_count[0] = (MPI_Count*)  NCI_Malloc(3 * nbufs * nprocs * sizeof(MPI_Count));

        /* recv_size[j][i] is the receive size from proc i in round j */
        recv_size = recv_count + nbufs;
        recv_size[0] = recv_count[0] + nbufs * nprocs;

        /* recv_start_pos[j][i] is the starting index of offset-length arrays
         * pointed by others_req[i].curr for remote rank i in round j
         */
        recv_start_pos = recv_size + nbufs;
        recv_start_pos[0] = recv_size[0] + nbufs * nprocs;

        for (j = 1; j < nbufs; j++) {
            write_buf[j] = write_buf[j-1] + striping_unit;
            /* recv_buf[j] may be realloc in LUSTRE_W_Exchange_data() */
            recv_buf[j]       = (char *) NCI_Malloc(striping_unit);
            recv_count[j]     = recv_count[j-1]     + nprocs;
            recv_size[j]      = recv_size[j-1]      + nprocs;
            recv_start_pos[j] = recv_start_pos[j-1] + nprocs;
        }

        /* srt_off_len consists of file offset-length pairs sorted in a
         * monotonically non-decreasing order (required by MPI-IO standard)
         * which is used when writing to the file
         */
        srt_off_len = (off_len_list*) NCI_Malloc(nbufs * sizeof(off_len_list));
    }

    /* send_buf[] will be allocated in LUSTRE_W_Exchange_data(), when the use
     * buffer is not contiguous.
     */
    send_buf = (char **) NCI_Malloc(nbufs * sizeof(char*));

    /* this_buf_idx contains indices to the user write buffer for sending this
     * rank's write data to aggregators, one for each aggregator. It is used
     * only when user buffer is contiguous.
     */
    if (buf_view.is_contig)
        this_buf_idx = (MPI_Offset *) NCI_Malloc(sizeof(MPI_Offset) * cb_nodes);

    /* array of data sizes to be sent to each aggregator in a 2-phase round */
    send_size = (MPI_Count *) NCI_Calloc(cb_nodes, sizeof(MPI_Count));

    /* min_st_loc is the beginning file offsets of the aggregate access region
     *     of this collective write, and it has been downward aligned to the
     *     nearest file stripe boundary
     * iter_end_off is the ending file offset of aggregate write region of
     *     iteration m, upward aligned to the file stripe boundary.
     */
    iter_end_off = min_st_loc + step_size;

    ibuf = 0;
    for (m = 0; m < ntimes; m++) {
        MPI_Count range_size;
        MPI_Offset range_off;

        /* Note that MPI standard (MPI 3.1 Chapter 13.1.1 and MPI 4.0 Chapter
         * 14.1.1) requires that the typemap displacements of etype and
         * filetype are non-negative and monotonically non-decreasing. This
         * simplifies implementation a bit compared to reads.
         */

        /* Calculate what should be communicated.
         *
         * First, calculate the amount to be sent to each aggregator i, at this
         * round m, by going through all offset-length pairs in my_req[i].
         *
         * iter_end_off - ending file offset of aggregate write region of this
         *                round, and upward aligned to the file stripe
         *                boundary. Note the aggregate write region of this
         *                round starts from (iter_end_off-step_size) to
         *                iter_end_off, aligned with file stripe boundaries.
         * send_size[i] - total size in bytes of this process's write data
         *                fall into aggregator i's FD in this round.
         * recv_size[m][i] - size in bytes of data to be received by this
         *                aggregator from process i in round m.
         * recv_count[m][i] - number of noncontiguous offset-length pairs from
         *                process i fall into this aggregator's write region
         *                in round m.
         */
        for (i = 0; i < cb_nodes; i++) {
            /* reset communication metadata to all 0s for this round */
            send_size[i] = 0;

            if (my_req[i].count == 0) continue;
            /* my_req[i].count is the number of this rank's offset-length pairs
             * to be sent to aggregator i
             */

            if (my_req[i].curr == my_req[i].count)
                continue; /* done with aggregator i */

            if (buf_view.is_contig)
                /* buf_idx is used only when user buffer is contiguous.
                 * this_buf_idx[i] points to the starting offset of user
                 * buffer, buf, for amount of send_size[i] to be sent to
                 * aggregator i at this round.
                 */
                this_buf_idx[i] = buf_idx[i][my_req[i].curr];

            /* calculate the send amount from this rank to aggregator i */
            for (j = my_req[i].curr; j < my_req[i].count; j++) {
                if (my_req[i].offsets[j] < iter_end_off)
                    send_size[i] += my_req[i].lens[j];
                else
                    break;
            }

            /* update my_req[i].curr to point to the jth offset-length
             * pair of my_req[i], which will be used as the first pair in the
             * next round of iteration.
             */
            my_req[i].curr = j;
        }

        /* range_off is the starting file offset of this aggregator's write
         *     region at this round (may not be aligned to stripe boundary).
         * range_size is the size (in bytes) of this aggregator's write region
         *     for this round (whose size is always <= striping_unit).
         */
        range_off = off_list[m];
        range_size = MIN(striping_unit - range_off % striping_unit,
                         end_loc - range_off + 1);

        /* Calculate the amount to be received from each process i at this
         * round, by going through all offset-length pairs of others_req[i].
         */
        if (recv_count != NULL) {
            for (i=0; i<nprocs; i++) {
                /* reset communication metadata to all 0s for this round */
                recv_count[ibuf][i] = recv_size[ibuf][i] = 0;
                recv_start_pos[ibuf][i] = 0;

                if (others_req[i].count == 0) continue;

                recv_start_pos[ibuf][i] = others_req[i].curr;
                for (j = others_req[i].curr; j < others_req[i].count; j++) {
                    if (others_req[i].offsets[j] < iter_end_off) {
                        recv_count[ibuf][i]++;
                        others_req[i].mem_ptrs[j] = others_req[i].offsets[j]
                                                  - range_off;
                        recv_size[ibuf][i] += others_req[i].lens[j];
                    } else {
                        break;
                    }
                }
                /* update others_req[i].curr to point to the jth offset-length
                 * pair of others_req[i], which will be used as the first pair
                 * in the next round of iteration.
                 */
                others_req[i].curr = j;
            }
        }
        iter_end_off += step_size;

        /* exchange phase - each process sends it's write data to I/O
         * aggregators and aggregators receive from non-aggregators.
         * Communication are MPI_Issend and MPI_Irecv only. There is no
         * collective communication. Only aggregators have non-NULL write_buf
         * and recv_buf. All processes have non-NULL send_buf.
         */
        char *wbuf = (write_buf == NULL) ? NULL : write_buf[ibuf];

        /* Exchange_data_recv() and Exchange_data_send() below perform one
         * round of communication phase and there are ntimes rounds.
         */
// printf("%s at %d: end_loc=%lld nbufs=%d recv_list=%s\n",__func__,__LINE__, end_loc,nbufs,(recv_list==NULL)?"NULL":"NOT NULL");
        if (recv_list != NULL) { /* this aggregator has something to received */
            char *rbuf = (recv_buf  == NULL) ? NULL :  recv_buf[ibuf];
            int err;

            err = Exchange_data_recv(fd,
                               buf,                /* IN: user buffer */
                               wbuf,               /* OUT: write buffer */
                               &rbuf,              /* OUT: receive buffer */
                               &buf_view,
                               recv_size[ibuf],     /* IN: changed each round */
                               range_off,           /* IN: changed each round */
                               range_size,          /* IN: changed each round */
                               recv_count[ibuf],    /* IN: changed each round */
                               recv_start_pos[ibuf],/* IN: changed each round */
                               others_req,          /* IN: changed each round */
                               this_buf_idx,        /* IN: changed each round */
                               &srt_off_len[ibuf],/* OUT: write off-len pairs */
                               recv_list);        /* OUT: recv disp-len pairs */
            if (err != NC_NOERR)
                goto over;

            /* rbuf might be realloc-ed */
            if (recv_buf != NULL) recv_buf[ibuf] = rbuf;
        }

        /* sender part */
        MPI_Count self_count, self_start_pos;
        if (recv_count == NULL) {
            self_count = 0;
            self_start_pos = 0;
        }
        else {
            self_count     = recv_count[ibuf][myrank];
            self_start_pos = recv_start_pos[ibuf][myrank];
        }
        send_buf[ibuf] = NULL;

        Exchange_data_send(fd,
                           buf,             /* IN: user buffer */
                           wbuf,            /* OUT: write buffer */
                           &send_buf[ibuf], /* OUT: send buffer */
                           &buf_view,
                           send_size,       /* IN: changed each round */
                           self_count,
                           self_start_pos,
                           others_req,      /* IN: changed each round */
                           this_buf_idx,    /* IN: changed each round */
                           send_list);      /* OUT: send disp-len pairs */

        if (m % nbufs < nbufs - 1 && m < ntimes - 1) {
            /* continue to the next round */
            ibuf++;
        }
        else {
            /* commit communication and write this batch of numBufs to file */
            int numBufs = ibuf + 1;

// printf("%s at %d: m=%d nbufs=%d ntimes=%d\n",__func__,__LINE__, m,nbufs,ntimes);
            /* reset ibuf to the first element of nbufs */
            ibuf = 0;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
            double curT = MPI_Wtime();
#endif
            /* communication phase */
            commit_comm_phase(fd, send_list, recv_list);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
            if (fd->is_agg) fd->write_timing[3] += MPI_Wtime() - curT;
#endif

            /* free send_buf allocated in LUSTRE_W_Exchange_data() */
            for (j = 0; j < numBufs; j++) {
                if (send_buf[j] != NULL) {
                    NCI_Free(send_buf[j]);
                    send_buf[j] = NULL;
                }
            }
            if (!fd->is_agg) /* non-aggregators are done for this batch */
                continue;

            if (recv_list == NULL) /*  this aggregator has nothing to write */
                continue;

            /* this aggregator unpacks the data in recv_buf[] into write_buf */
            if (end_loc >= 0) {
                for (j = 0; j < numBufs; j++) {
                    char *buf_ptr = recv_buf[j];
                    for (i = 0; i < nprocs; i++) {
                        if (recv_count[j][i] > 1 && i != myrank) {
                            /* When recv_count[j][i] == 1, this case has
                             * been taken care of earlier by receiving the
                             * message directly into write_buf.
                             */
                            MEMCPY_UNPACK(i, buf_ptr, recv_start_pos[j][i],
                                          recv_count[j][i], write_buf[j]);
                            buf_ptr += recv_size[j][i];
                        }
                    }
                }
            }

            /* this aggregator writes to numBufs number of stripes */
            for (j=0; j<numBufs; j++) {

                /* if there is no data to write in round (batch_idx + j) */
                if (srt_off_len[j].num == 0)
                    continue;

                /* range_off  starting file offset of this aggregator's write
                 *            region for this round (may not be aligned to
                 *            stripe boundary)
                 * range_size size (in bytes) of this rank's write region for
                 *            this round, <= striping_unit
                 */
                range_off = off_list[batch_idx + j];
                range_size = MIN(striping_unit - range_off % striping_unit,
                                 end_loc - range_off + 1);

                /* When srt_off_len[j].num == 1, either there is no hole in the
                 * write buffer or the file domain has been read-modify-written
                 * with the received write data. When srt_off_len[j].num > 1,
                 * data sieving is not performed and holes have been found. In
                 * this case, srt_off_len[] is the list of sorted offset-length
                 * pairs describing noncontiguous writes. Now call writes for
                 * each offset-length pair. Note the offset-length pairs
                 * (represented by srt_off_len[j].off, srt_off_len[j].len, and
                 * srt_off_len[j].num) have been coalesced in
                 * LUSTRE_W_Exchange_data().
                 */
// printf("%s at %d: num=%d\n",__func__,__LINE__, srt_off_len[j].num);
                for (i = 0; i < srt_off_len[j].num; i++) {
                    /* all write requests in this round should fall into file
                     * range of [range_off, range_off+range_size). This below
                     * assertion should never fail.
                     */
                    assert(srt_off_len[j].off[i] < range_off + range_size &&
                           srt_off_len[j].off[i] >= range_off);

// printf("%s at %d: PNCIO_WriteContig num=%d [%d] off=%lld len=%lld\n",__func__,__LINE__, srt_off_len[j].num,i,srt_off_len[j].off[i],srt_off_len[j].len[i]);
                    w_len = PNCIO_WriteContig(fd,
                                     write_buf[j] + (srt_off_len[j].off[i] - range_off),
                                     srt_off_len[j].len[i],
                                     srt_off_len[j].off[i]);
                    if (w_len < 0) goto over;
                    total_w_len += w_len;
                }
                if (srt_off_len[j].num > 0) {
                    NCI_Free(srt_off_len[j].off);
                    srt_off_len[j].num = 0;
                }
            }
            batch_idx += numBufs; /* only matters for aggregators */
        }
    }

  over:
    if (srt_off_len)
        NCI_Free(srt_off_len);
    if (write_buf != NULL)
        NCI_Free(write_buf);
    if (recv_buf != NULL) {
        for (j = 0; j < nbufs; j++)
            NCI_Free(recv_buf[j]);
        NCI_Free(recv_buf);
    }
    if (recv_count != NULL) {
        NCI_Free(recv_count[0]);
        NCI_Free(recv_count);
    }
    NCI_Free(send_size);
    NCI_Free(off_list);
    if (buf_view.is_contig)
        NCI_Free(this_buf_idx);
    if (send_buf != NULL)
        NCI_Free(send_buf);
    if (send_list != NULL) {
        for (i = 0; i < cb_nodes; i++)
            NCI_Free(send_list[i].disp);
        NCI_Free(send_list);
    }
    if (recv_list != NULL) {
        for (i = 0; i < nprocs; i++)
            NCI_Free(recv_list[i].disp);
        NCI_Free(recv_list);
    }

#ifdef WKL_DEBUG
    /* check any pending messages to be received */
    MPI_Status probe_st;
    int probe_flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, fd->comm, &probe_flag, &probe_st);
    if (probe_flag) {
        printf("ERROR ++++ MPI_Iprobe rank=%4d is_agg=%d: ---- cb_nodes=%d ntimes=%lld nbufs=%d\n",myrank,fd->is_agg,cb_nodes,ntimes,nbufs);
        fflush(stdout);
    }
#endif
    return total_w_len;
}

/* This heap-merge sort also coalesces sorted offset-length pairs whenever
 * possible.
 *
 * Heapify(a, i, heapsize); Algorithm from Cormen et al. pg. 143 modified for a
 * heap with smallest element at root. The recursion has been removed so that
 * there are no function calls. Function calls are too expensive.
 */
static
void heap_merge(const PNCIO_Access *others_req,
                const MPI_Count    *count,
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count          *srt_off,
                MPI_Count          *srt_len,
#else
                MPI_Offset         *srt_off,
                int                *srt_len,
#endif
                const MPI_Count    *start_pos,
                int                 nprocs,
                int                 nprocs_recv,
                MPI_Count          *total_elements)
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
    for (i = 0; i < nprocs; i++) {
        if (count[i]) {
            a[j].off_list = others_req[i].offsets + start_pos[i];
            a[j].len_list = others_req[i].lens + start_pos[i];
            a[j].nelem = count[i];
            j++;
        }
    }

#define SWAP(x, y, tmp) { tmp = x ; x = y ; y = tmp ; }

    heapsize = nprocs_recv;

    /* Build a heap out of the first element from each list, with the smallest
     * element of the heap at the root. The first for loop is to find and move
     * the smallest a[*].off_list[0] to a[0].
     */
    for (i = heapsize / 2 - 1; i >= 0; i--) {
        k = i;
        for (;;) {
            r = 2 * (k + 1);
            l = r - 1;
            if ((l < heapsize) && (*(a[l].off_list) < *(a[k].off_list)))
                smallest = l;
            else
                smallest = k;

            if ((r < heapsize) && (*(a[r].off_list) < *(a[smallest].off_list)))
                smallest = r;

            if (smallest != k) {
                SWAP(a[k], a[smallest], tmp);
                k = smallest;
            } else
                break;
        }
    }

    /* The heap keeps the smallest element in its first element, i.e.
     * a[0].off_list[0].
     */
    j = 0;
    for (i = 0; i < *total_elements; i++) {
        /* extract smallest element from heap, i.e. the root */
        if (j == 0 || srt_off[j - 1] + srt_len[j - 1] < *(a[0].off_list)) {
            srt_off[j] = *(a[0].off_list);
            srt_len[j] = *(a[0].len_list);
            j++;
        } else {
            /* this offset-length pair can be coalesced into the previous one */
            srt_len[j - 1] = *(a[0].off_list) + *(a[0].len_list) - srt_off[j - 1];
        }
        (a[0].nelem)--;

        if (a[0].nelem) {
            (a[0].off_list)++;
            (a[0].len_list)++;
        } else {
            a[0] = a[heapsize - 1];
            heapsize--;
        }

        /* Heapify(a, 0, heapsize); */
        k = 0;
        for (;;) {
            r = 2 * (k + 1);
            l = r - 1;
            if ((l < heapsize) && (*(a[l].off_list) < *(a[k].off_list)))
                smallest = l;
            else
                smallest = k;

            if ((r < heapsize) && (*(a[r].off_list) < *(a[smallest].off_list)))
                smallest = r;

            if (smallest != k) {
                SWAP(a[k], a[smallest], tmp);
                k = smallest;
            } else
                break;
        }
    }
    NCI_Free(a);
    *total_elements = j;
}

#define CACHE_REQ(list, nelems, buf) {   \
    MPI_Aint buf_addr;                   \
    list.len[list.count] = nelems;       \
    MPI_Get_address(buf, &buf_addr);     \
    list.disp[list.count] = buf_addr;    \
    list.count++;                        \
}

static
int Exchange_data_recv(
          PNCIO_File     *fd,
    const void           *buf,         /* user buffer */
          char           *write_buf,   /* OUT: internal buffer used to write
                                        * to file */
          char          **recv_buf,    /* OUT: [nbufs] internal buffer used to
                                        * receive from other processes */
    const PNCIO_View     *buf_view,    /* IN: flattened buffer
                                        * offset-length pairs */
    const MPI_Count      *recv_size,   /* [nprocs] recv_size[i] is amount of
                                        * this aggregator recv from rank i */
          MPI_Offset      range_off,   /* starting file offset of this
                                        * aggregator's write region */
          MPI_Count       range_size,  /* amount of this aggregator's write
                                        * region */
    const MPI_Count      *recv_count,  /* [nprocs] recv_count[i] is the number
                                        * of offset-length pairs received from
                                        * rank i */
    const MPI_Count      *start_pos,   /* [nprocs] start_pos[i] starting value
                                        * of others_req[i].curr */
    const PNCIO_Access   *others_req,  /* [nprocs] others_req[i] is rank i's
                                        * write requests fall into this
                                        * aggregator's file domain */
    const MPI_Offset     *buf_idx,      /* [cb_nodes] indices to user buffer
                                        * offsets for sending this rank's
                                        * write data to aggregator i */
          off_len_list   *srt_off_len, /* OUT: list of write offset-length
                                        * pairs of this aggregator */
          disp_len_list  *recv_list)   /* OUT: displacement-length pairs of
                                        * recv buffer */
{
    char *buf_ptr, *contig_buf;
    size_t alloc_sz;
    int i, j, nprocs, myrank, nprocs_recv, hole, build_srt_off_len;
    MPI_Count sum_recv;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    /* srt_off_len contains the file offset-length pairs to be written by this
     * aggregator at this round. The file region starts from range_off with
     * size of range_size.
     */

    srt_off_len->num = 0;
    srt_off_len->off = NULL;
    sum_recv = 0;
    nprocs_recv = 0;

    /* calculate receive metadata */
    j = -1;
    for (i = 0; i < nprocs; i++) {
        srt_off_len->num += recv_count[i];
        if (j == -1 && recv_count[i] > 0) j = i;
        sum_recv += recv_size[i];
        if (recv_size[i])
            nprocs_recv++;
    }

    if (nprocs_recv == 0) return NC_NOERR;

// MPI_Count numx = srt_off_len->num; printf("nprocs_recv=%d PNCIO_DS_WR_NAGGRS_LB=%d srt_off_len->num=%lld PNCIO_DS_WR_NPAIRS_LB=%d\n",nprocs_recv,PNCIO_DS_WR_NAGGRS_LB,srt_off_len->num,PNCIO_DS_WR_NPAIRS_LB);

    /* determine whether checking holes is necessary */
    if (srt_off_len->num == 0) {
        /* this process has nothing to receive and hence no hole */
        build_srt_off_len = 0;
        hole = 0;
    } else if (srt_off_len->num == 1) {
        build_srt_off_len = 0;
        hole = 0;
#ifdef HAVE_MPI_LARGE_COUNT
        alloc_sz = sizeof(MPI_Offset) + sizeof(MPI_Count);
        srt_off_len->off = (MPI_Offset*) NCI_Malloc(alloc_sz);
        srt_off_len->len = (MPI_Count*) (srt_off_len->off + 1);
#else
        alloc_sz = sizeof(MPI_Offset) + sizeof(int);
        srt_off_len->off = (MPI_Offset*) NCI_Malloc(alloc_sz);
        srt_off_len->len = (int*) (srt_off_len->off + 1);
#endif
        srt_off_len->off[0] = others_req[j].offsets[start_pos[j]];
        srt_off_len->len[0] = others_req[j].lens[start_pos[j]];
    } else if (fd->hints->romio_ds_write == PNCIO_HINT_ENABLE) {
        /* skip building of srt_off_len and proceed to read-modify-write */
        build_srt_off_len = 0;
        /* assuming there are holes */
        hole = 1;
    } else if (fd->hints->romio_ds_write == PNCIO_HINT_AUTO) {
        if (DO_HEAP_MERGE(nprocs_recv, srt_off_len->num)) {
            /* When the number of sorted offset-length lists or the total
             * number of offset-length pairs are too large, the heap-merge sort
             * below for building srt_off_len can become very expensive. Such
             * sorting is also used to check holes to determine whether
             * read-modify-write is necessary.
             */
            build_srt_off_len = 0;
            /* assuming there are holes */
            hole = 1;
        }
        else /* heap-merge is less expensive, proceed to build srt_off_len */
            build_srt_off_len = 1;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        if (build_srt_off_len) {
            fd->write_counter[1]++;
            fd->write_counter[2] = MAX(fd->write_counter[2], srt_off_len->num);
            fd->write_counter[3] = MAX(fd->write_counter[3], nprocs_recv);
        } else {
            fd->write_counter[4]++;
            fd->write_counter[5] = MAX(fd->write_counter[5], srt_off_len->num);
            fd->write_counter[6] = MAX(fd->write_counter[6], nprocs_recv);
        }
#endif
    } else { /* if (fd->hints->romio_ds_write == PNCIO_HINT_DISABLE) */
        /* User explicitly disable data sieving to skip read-modify-write.
         * Whether or not there is a hole is not important. However,
         * srt_off_len must be constructed to merge all others_req[] into a
         * single sorted list. This step is necessary because after this
         * subroutine returns, write data from all non-aggregators will be
         * packed into the write_buf, with a possibility of overlaps, and
         * as srt_off_len stores the coalesced offset-length pairs of
         * individual non-contiguous write requests, it is used to write them
         * to the file.
         */
        build_srt_off_len = 1;
    }

    if (build_srt_off_len) {
        /* merge all the offset-length pairs from others_req[] (already sorted
         * individually) into a single list of offset-length pairs.
         */
#ifdef HAVE_MPI_LARGE_COUNT
        alloc_sz = sizeof(MPI_Offset) + sizeof(MPI_Count);
        srt_off_len->off = (MPI_Offset*) NCI_Malloc(alloc_sz * srt_off_len->num);
        srt_off_len->len = (MPI_Count*) (srt_off_len->off + srt_off_len->num);
#else
        alloc_sz = sizeof(MPI_Offset) + sizeof(int);
        srt_off_len->off = (MPI_Offset*) NCI_Malloc(alloc_sz * srt_off_len->num);
        srt_off_len->len = (int*) (srt_off_len->off + srt_off_len->num);
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        double curT = MPI_Wtime();
#endif
        heap_merge(others_req, recv_count, srt_off_len->off, srt_off_len->len,
                   start_pos, nprocs, nprocs_recv, &srt_off_len->num);

        /* Now, (srt_off_len->off and srt_off_len->len) are in an increasing
         * order of file offsets. In addition, they are coalesced.
         */
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        fd->write_timing[5] += MPI_Wtime() - curT;
#endif
        /* whether or not there are holes */
        hole = (srt_off_len->num > 1);
    }

// printf("%s at %d: romio_ds_write=%s build_srt_off_len=%d hole=%d skip_read=%d srt_off_len->num=%lld\n",__func__,__LINE__, (fd->hints->romio_ds_write == PNCIO_HINT_ENABLE)?"ENABLE": (fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)?"DISABLE":"AUTO", build_srt_off_len,hole,fd->skip_read,srt_off_len->num);
// printf("%s at %d: romio_ds_write=%s build_srt_off_len=%d hole=%d nprocs_recv=%d(PNCIO_DS_WR_NAGGRS_LB=%d) numx=%lld(PNCIO_DS_WR_NPAIRS_LB=%d)\n",__func__,__LINE__, (fd->hints->romio_ds_write == PNCIO_HINT_ENABLE)?"ENABLE": (fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)?"DISABLE":"AUTO", build_srt_off_len,hole,nprocs_recv,PNCIO_DS_WR_NAGGRS_LB,numx,PNCIO_DS_WR_NPAIRS_LB);

    /* data sieving */
    if (fd->hints->romio_ds_write != PNCIO_HINT_DISABLE && hole) {
        if (fd->skip_read)
            memset(write_buf, 0, range_size);
        else {
            MPI_Offset r_len;
            r_len = PNCIO_ReadContig(fd, write_buf, range_size, range_off);
            if (r_len < 0) return (int)r_len;
        }

        /* Once read, holes have been filled and thus the number of
         * offset-length pairs, srt_off_len->num, becomes one.
         */
        srt_off_len->num = 1;
        if (srt_off_len->off == NULL) { /* if has not been malloc-ed yet */
#ifdef HAVE_MPI_LARGE_COUNT
            alloc_sz = sizeof(MPI_Offset) + sizeof(MPI_Count);
            srt_off_len->off = (MPI_Offset*) NCI_Malloc(alloc_sz);
            srt_off_len->len = (MPI_Count*) (srt_off_len->off + 1);
#else
            alloc_sz = sizeof(MPI_Offset) + sizeof(int);
            srt_off_len->off = (MPI_Offset*) NCI_Malloc(alloc_sz);
            srt_off_len->len = (int*) (srt_off_len->off + 1);
#endif
        }
        srt_off_len->off[0] = range_off;
        srt_off_len->len[0] = range_size;
    }

    /* It is possible sum_recv (sum of message sizes to be received) is larger
     * than the size of collective buffer, write_buf, if writes from multiple
     * remote processes overlap. Receiving messages into overlapped regions of
     * the same write_buffer may cause a problem. To avoid it, we allocate a
     * temporary buffer big enough to receive all messages into disjointed
     * regions. Earlier in LUSTRE_Exch_and_write(), write_buf is already
     * allocated with twice amount of the file stripe size, with the second
     * half to be used to receive messages. If sum_recv is smaller than file
     * stripe size, we can reuse that space. But if sum_recv is bigger (an
     * overlap case, which is rare), we allocate a separate buffer of size
     * sum_recv.
     */
    sum_recv -= recv_size[myrank];
    if (sum_recv > fd->hints->striping_unit)
        *recv_buf = (char *) NCI_Realloc(*recv_buf, sum_recv);
    contig_buf = *recv_buf;

    /* cache displacement-length pairs of receive buffer */
    buf_ptr = contig_buf;
    for (i = 0; i < nprocs; i++) {
        if (recv_size[i] == 0)
            continue;
        if (i != myrank) {
            if (recv_count[i] > 1) {
                CACHE_REQ(recv_list[i], recv_size[i], buf_ptr)
                buf_ptr += recv_size[i];
            } else {
                /* recv_count[i] is the number of noncontiguous offset-length
                 * pairs describing the write requests of rank i that fall
                 * into this aggregator's file domain. When recv_count[i] is 1,
                 * there is only one such pair, meaning the receive message is
                 * to be stored contiguously. Such message can be received
                 * directly into write_buf.
                 */
                CACHE_REQ(recv_list[i], recv_size[i],
                          write_buf + others_req[i].mem_ptrs[start_pos[i]])
            }
        } else if (buf_view->is_contig && recv_count[i] > 0) {
            /* send/recv to/from self uses memcpy(). The case when buftype is
             * not contiguous will be handled later in Exchange_data_send().
             */
            char *fromBuf = (char *) buf + buf_idx[fd->my_cb_nodes_index];
            MEMCPY_UNPACK(i, fromBuf, start_pos[i], recv_count[i], write_buf);
        }
    }
    return NC_NOERR;
}

static
void Exchange_data_send(
          PNCIO_File     *fd,
    const void           *buf,          /* user buffer */
          char           *write_buf,    /* OUT: internal buffer used to write
                                         * to file, only matter when send to
                                         * self */
          char          **send_buf_ptr, /* OUT: [cb_nodes] point to internal
                                         * send buffer */
          PNCIO_View     *buf_view,     /* IN/OUT: flattened buffer
                                         * offset-length pairs */
    const MPI_Count      *send_size,    /* [cb_nodes] send_size[i] is amount of
                                         * this rank sent to aggregator i */
          MPI_Count       self_count,   /* No. offset-length pairs sent to self
                                         * rank */
          MPI_Count       start_pos,    /* others_req[myrank].curr */
    const PNCIO_Access   *others_req,   /* [nprocs] only used when send to self,
                                         * others_req[myrank] */
    const MPI_Offset     *buf_idx,      /* [cb_nodes] indices to user buffer
                                         * for sending this rank's write data
                                         * to aggregator i */
          disp_len_list  *send_list)    /* OUT: displacement-length pairs of
                                         * send buffer */
{
    int i, myrank, cb_nodes;

    *send_buf_ptr = NULL;

    MPI_Comm_rank(fd->comm, &myrank);

    cb_nodes = fd->hints->cb_nodes;
// if (myrank==0) printf("%s at %d: cb_nodes=%d\n",__func__,__LINE__, cb_nodes);
    if (buf_view->is_contig) {
        /* If buftype is contiguous, data can be directly sent from user buf
         * at location given by buf_idx.
         */
        for (i = 0; i < cb_nodes; i++) {
// if (myrank==0 && send_size[i]) printf("%s at %d: cb_nodes=%d send_size[%d]=%lld my_cb_nodes_index=%d\n",__func__,__LINE__, cb_nodes,i,send_size[i],fd->my_cb_nodes_index);
            if (send_size[i] && i != fd->my_cb_nodes_index)
                CACHE_REQ(send_list[i], send_size[i], (char*)buf + buf_idx[i]);
        }
    } else {
        char **send_buf, *self_buf;

        /* total send size of this round */
        size_t send_total_size = 0;
        for (i = 0; i < cb_nodes; i++)
            send_total_size += send_size[i];

        if (send_total_size == 0) return;

        /* The user buffer to be used to send in this round is not contiguous,
         * allocate send_buf[], a contiguous space, copy data to send_buf,
         * including ones to be sent to self, and then use send_buf to send.
         */
        send_buf = (char **) NCI_Malloc(cb_nodes * sizeof(char *));
        send_buf[0] = (char *) NCI_Malloc(send_total_size);
        for (i = 1; i < cb_nodes; i++)
            send_buf[i] = send_buf[i - 1] + send_size[i - 1];

        LUSTRE_Fill_send_buffer(fd, buf, buf_view, send_buf,
                                send_total_size, send_size, &self_buf,
                                send_list);
        /* Send buffers must not be touched before MPI_Waitall() is completed,
         * and thus send_buf will be freed in LUSTRE_Exch_and_write()
         */

        if (fd->my_cb_nodes_index >= 0 && send_size[fd->my_cb_nodes_index] > 0) {
            /* contents of user buf that must be sent to self has been copied
             * into send_buf[fd->my_cb_nodes_index]. Now unpack it into
             * write_buf.
             */
            if (self_buf == NULL) self_buf = send_buf[fd->my_cb_nodes_index];
            MEMCPY_UNPACK(myrank, self_buf, start_pos, self_count, write_buf);
        }

        *send_buf_ptr = send_buf[0];
        NCI_Free(send_buf);
    }
}

static void LUSTRE_Fill_send_buffer(PNCIO_File       *fd,
                                    const void       *buf,
                                    PNCIO_View       *buf_view, /* IN/OUT */
                                    char            **send_buf,
                                    size_t            send_total_size,
                                    const MPI_Count  *send_size,
                                    char            **self_buf,
                                    disp_len_list    *send_list)
{
    /* this function is only called if buftype is not contiguous */
    int q, first_q=-1, isUserBuf=0;
    MPI_Count send_size_rem=0, size, copy_size=0;
    char *user_buf_ptr=NULL, *send_buf_ptr=NULL, *same_buf_ptr=NULL;
    MPI_Offset off, user_buf_idx;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset len, rem_len;
#else
    int len, rem_len;
#endif

#ifdef WKL_DEBUG
int num_memcpy=0;
#endif

    *self_buf = NULL;

    /* user_buf_idx is to the index offset to buf, indicating the starting
     * location to be copied.
     *
     * buf_view stores the offset-length pairs of the flattened user buffer
     *     data type. Note this stores offset-length pairs of the data type,
     *     and write amount can be a multiple of the data type.
     * buf_view.count: the number of pairs
     * buf_view.off[i]: the ith pair's byte offset to buf. Note the
     *     flattened offsets of user buffer type may not be sorted in an
     *     increasing order, unlike fileview which is required by MPI to be
     *     sorted in a monotonically non-decreasing order.
     * buf_view.len[i]: length of the ith pair
     * buf_view.idx: index to the offset-length pair currently being
     *     processed, incremented each round.
     * buf_view.rem: amount of data in the pair that has not been copied
     *     over, changed each round.
     */
    user_buf_idx = buf_view->off[buf_view->idx]
                 + buf_view->len[buf_view->idx]
                 - buf_view->rem;
                 /* in case data left to be copied from previous round */

    /* fd->flat_file.count: the number of noncontiguous file segments this
     *     rank writes to. Each segment i is described by fd->flat_file.offs[i]
     *     and fd->flat_file.len[i].
     * fd->flat_file.idx: the index to the fd->flat_file.offs[],
     *     fd->flat_file.len[] that have been processed in the previous round.
     * The while loop below packs write data into send buffers, send_buf[],
     * based on this rank's off-len pairs in its file view,
     */
    off     = fd->flat_file.off[fd->flat_file.idx]
            + fd->flat_file.len[fd->flat_file.idx]
            - fd->flat_file.rem;
    rem_len = fd->flat_file.rem;

// int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    while (send_total_size > 0) {
        /* this off-len request may span to more than one I/O aggregator */
// if (rank == 0) printf("rank 0 %s at %d send_total_size=%zd rem_len=%lld\n",__func__,__LINE__,send_total_size,rem_len);
        while (rem_len != 0) {
            len = rem_len;
            q = LUSTRE_Calc_aggregator(fd, off, &len);
            /* NOTE: len will be modified by PNCIO_Calc_aggregator() to be no
             * more than a file stripe unit size that aggregator "q" is
             * responsible for. Note q is not the MPI rank ID, It is the array
             * index to fd->hints->ranklist[].
             *
             * Now len is the amount of data in ith off-len pair that should be
             * sent to aggregator q. Note q can also be self. In this case,
             * data is also packed into send_buf[q] or pointed to a segment of
             * buf when the data to be packed is contiguous. send_buf[q] will
             * later be copied to write buffer in MEMCPY_UNPACK, instead of
             * calling MPI_Issend to send.
             *
             * send_size[q]: data amount of this rank needs to send to
             * aggregator q in this round.
             *
             * len and send_size[q] are all always <= striping_unit
             */

// if (rank == 0) printf("rank 0 %s at %d rem_len=%lld len=%lld first_q=%d q=%d idx=%lld\n",__func__,__LINE__,rem_len,len,first_q,q,buf_view->idx);

            if (first_q != q) {
                assert(send_size_rem == 0);
                first_q = q;
                isUserBuf = 1;
                send_size_rem = send_size[q];
                copy_size = 0;
                same_buf_ptr = (char*)buf + user_buf_idx; /* no increment */
                user_buf_ptr = same_buf_ptr; /* increment after each memcpy */
                if (send_buf != NULL)
                    send_buf_ptr = send_buf[q]; /* increment after each memcpy */
            }

            /* copy len amount of data from buf to send_buf[q] */
            size = len;

            while (size) {
                MPI_Count size_in_buf = MIN(size, buf_view->rem);
                copy_size += size_in_buf;
                user_buf_idx += size_in_buf;
                send_size_rem -= size_in_buf;
                buf_view->rem -= size_in_buf;
// if (rank == 0) printf("rank 0 %s at %d size=%lld size_in_buf=%lld copy_size=%lld rem=%ld\n",__func__,__LINE__, size, size_in_buf, copy_size,buf_view->rem);
                if (buf_view->rem == 0) { /* move on to next off-len pair */
                    if (! buf_view->is_contig) {
                        /* user buffer type is not contiguous */
                        if (send_size_rem) {
                            /* after this copy send_buf[q] is still not full */
                            isUserBuf = 0;
// if (rank == 0 && (char*)buf == (char*)user_buf_ptr) printf("rank 0 copy original buf 1 size=%lld user_buf_ptr=%p\n",copy_size,user_buf_ptr);
                            memcpy(send_buf_ptr, user_buf_ptr, copy_size);
user_buf_ptr += copy_size;
                            send_buf_ptr += copy_size;
                            copy_size = 0;
                        } else if (isUserBuf == 0) {
                            /* send_buf[q] is full and not using user buf,
                             * copy the remaining delayed data */
// if (rank == 0 && (char*)buf == (char*)user_buf_ptr) printf("rank 0 copy original buf 2 size=%lld\n",copy_size);
                            memcpy(send_buf_ptr, user_buf_ptr, copy_size);
user_buf_ptr += copy_size;
                        }
#ifdef WKL_DEBUG
num_memcpy++;
#endif
                    }
                    /* update buf_view->idx, buf_view->rem,
                     * and user_buf_idx
                     */
                        buf_view->idx++;
assert(buf_view->idx <= buf_view->count);

if (buf_view->idx < buf_view->count) {
                    user_buf_idx = buf_view->off[buf_view->idx];
                    buf_view->rem = buf_view->len[buf_view->idx];
                    user_buf_ptr = (char*) buf + user_buf_idx;
}
else assert(size - size_in_buf == 0);

                }
                else if (send_size_rem == 0 && isUserBuf == 0) {
                    /* buf_view->rem > 0, send_buf[q] is full, and not using
                     * user buf to send, copy the remaining delayed data
                     */
// if (rank == 0 && (char*)buf == (char*)user_buf_ptr) printf("rank 0 copy original buf 3 size=%lld\n",copy_size);
                    memcpy(send_buf_ptr, user_buf_ptr, copy_size);
#ifdef WKL_DEBUG
num_memcpy++;
#endif
                    user_buf_ptr += copy_size;
                }
                size -= size_in_buf;
            }

            if (send_size_rem == 0) { /* data to q is fully packed */
                first_q = -1;

                if (q != fd->my_cb_nodes_index) { /* send only if not self rank */
                    if (isUserBuf)
                        CACHE_REQ(send_list[q], send_size[q], same_buf_ptr)
                    else
                        CACHE_REQ(send_list[q], send_size[q], send_buf[q])
                }
                else if (isUserBuf) {
                    /* send buffer is also (part of) user's buf. Return the
                     * buffer pointer, so the self send data can be directly
                     * unpack from user buf to write buffer.
                     */
                    *self_buf = same_buf_ptr;
                }
            }
            /* len is the amount of data copied */
            off += len;
            rem_len -= len;
            fd->flat_file.rem -= len;
            send_total_size -= len;
            if (send_total_size == 0) break;
        }

        /* done with this off-len pair, move on to the next */
        if (fd->flat_file.rem == 0) {
            fd->flat_file.idx++;
            fd->flat_file.rem = fd->flat_file.len[fd->flat_file.idx];
        }
        off = fd->flat_file.off[fd->flat_file.idx];
        rem_len = fd->flat_file.rem;
    }

#ifdef WKL_DEBUG
if (num_memcpy> 0) printf("---- fd->flat_file.count=%lld fd->flat_file.idx=%lld buf_view->count=%lld num_memcpy=%d\n",fd->flat_file.count,fd->flat_file.idx,buf_view->count,num_memcpy);
#endif
}

