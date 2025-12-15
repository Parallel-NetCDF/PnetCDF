/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * This file contains the implementation of intra-node aggregation feature,
 * which is designed to improve performance for I/O patterns that contain many
 * noncontiguous requests interleaved among processes, with a wide aggregate
 * access region on each process that involves file stripes responsible by
 * almost all the file servers. By reducing the number of processes per node
 * to participate MPI-IO operations, this feature can effectively reduce the
 * communication contention, particularly often happened to jobs that run a
 * large the number of MPI processes per compute node.
 *
 * Users can enable this feature by setting the PnetCDF I/O hint named
 * 'nc_num_aggrs_per_node' to a positive integral value, indicating the desired
 * number of processes per compute node to be selected as the intra-node I/O
 * aggregators. Processes running on the same node are divided into groups.
 * The process with the lowest rank ID is selected as the I/O aggregator of
 * that group. Non-aggregators send their requests to their aggregators, and
 * then the aggregators make I/O requests to the file, i.e. only aggregators
 * make MPI-IO calls.
 *
 * Because communication within a node can be achieved by memory copy operation
 * and thus its cost is expected to be much lower than the inter-node
 * communication, this feature can effectively reduce the communication
 * congestion or exhaustion of message queues, due to many pending asynchronous
 * messages produced in the two-phase I/O, the strategy used to implement
 * MPI collective I/O.
 *
 * The concept of intra-node request aggregation and its performance results
 * are presented in the following paper.
 * Q. Kang, S. Lee, K. Hou, R. Ross, A. Agrawal, A. Choudhary, and W. Liao.
 * Improving MPI Collective I/O for High Volume Non-Contiguous Requests With
 * Intra-Node Aggregation. IEEE Transactions on Parallel and Distributed
 * Systems, 31(11):2682-2695, November 2020.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* strcmp() strdup() */
#include <assert.h>
#include <errno.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/* swap values of x and y */
#define SWAP1(x, y, tmp) { tmp = x ; x = y; y = tmp ; }

#ifdef HAVE_MPI_LARGE_COUNT
/* swap elements of arrays x, y, and corresponding lengths and bufAddr */
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    MPI_Count aint; \
    MPI_Count cint; \
    MPI_Count d0 = (x) - offsets; \
    MPI_Count d1 = (y) - offsets; \
    if (d0 != d1) { \
        SWAP1(*(x), *(y), cint); \
        SWAP1(lengths[d0], lengths[d1], cint); \
        if (bufAddr != NULL) \
            SWAP1(bufAddr[d0], bufAddr[d1], aint); \
    } \
}
#else
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    int int4; \
    MPI_Offset aint; \
    MPI_Offset d0 = (x) - offsets; \
    MPI_Offset d1 = (y) - offsets; \
    if (d0 != d1) { \
        SWAP1(*(x), *(y), aint); \
        SWAP1(lengths[d0], lengths[d1], int4); \
        if (bufAddr != NULL) \
            SWAP1(bufAddr[d0], bufAddr[d1], aint); \
    } \
}
#endif

#define MEDIAN(a,b,c) ((*(a) < *(b)) ? \
                      ((*(b) < *(c)) ? (b) : ((*(a) < *(c)) ? (c) : (a))) : \
                      ((*(b) > *(c)) ? (b) : ((*(a) < *(c)) ? (a) : (c))))

static
size_t bin_search(
#ifdef HAVE_MPI_LARGE_COUNT
                  MPI_Count key, MPI_Count *base,
#else
                  MPI_Offset key, MPI_Offset *base,
#endif
                  size_t nmemb);

/*----< qsort_off_len_buf() >------------------------------------------------*/
/* Sort three arrays of offsets, lengths, and buffer addresses based on array
 * offsets into an increasing order. This code is based on the qsort routine
 * from Bentley & McIlroy's "Engineering a Sort Function".
 */
static void
qsort_off_len_buf(MPI_Aint    num,
#ifdef HAVE_MPI_LARGE_COUNT
                  MPI_Count  *offsets,
                  MPI_Count  *lengths,
#else
                  MPI_Offset *offsets,
                  int        *lengths,
#endif
                  MPI_Aint   *bufAddr)
{
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *pa, *pb, *pc, *pd, *pl, *pm, *pn, cmp_result, swap_cnt, i, r;
#else
    MPI_Offset *pa, *pb, *pc, *pd, *pl, *pm, *pn, cmp_result, swap_cnt, i, r;
#endif

    while (1) {
        swap_cnt = 0;
        pm = offsets + (num / 2);
        if (num > 7) {
            pl = offsets;
            pn = offsets + (num - 1);
            if (num > 40) {
                size_t d = (num / 8);
                pl = MEDIAN(pl, pl + d, pl + 2 * d);
                pm = MEDIAN(pm - d, pm, pm + d);
                pn = MEDIAN(pn - 2 * d, pn - d, pn);
            }
            pm = MEDIAN(pl, pm, pn);
        }
        SWAP(offsets, lengths, bufAddr, offsets, pm);
        pa = pb = offsets;

        pc = pd = offsets + (num - 1);
        for (;;) {
            while (pb <= pc && (cmp_result = (*pb - *offsets)) <= 0) {
                if (cmp_result == 0) {
                    swap_cnt = 1;
                    SWAP(offsets, lengths, bufAddr, pa, pb);
                    pa++;
                }
                pb++;
            }
            while (pb <= pc && (cmp_result = (*pc - *offsets)) >= 0) {
                if (cmp_result == 0) {
                    swap_cnt = 1;
                    SWAP(offsets, lengths, bufAddr, pc, pd);
                    pd--;
                }
                pc--;
            }
            if (pb > pc)
                break;
            SWAP(offsets, lengths, bufAddr, pb, pc);
            swap_cnt = 1;
            pb++;
            pc--;
        }
        if (swap_cnt == 0) {  /* Switch to insertion sort */
            for (pm = offsets; pm < offsets + num; pm++)
                for (pl = pm; pl > offsets && (*(pl-1) > *pl); pl--)
                    SWAP(offsets, lengths, bufAddr, pl, pl-1);
            return;
        }

        pn = offsets + num;
        r = MIN(pa - offsets, pb - pa);
        for (i=0; i<r; i++) SWAP(offsets, lengths, bufAddr, offsets+i, pb-r+i)

        r = MIN(pd - pc, pn - pd - 1);
        for (i=0; i<r; i++) SWAP(offsets, lengths, bufAddr, pb+i, pn-r+i)

        if ((r = pb - pa) > 1)
            qsort_off_len_buf(r, offsets, lengths, bufAddr);
        if ((r = pd - pc) > 1) {
            /* Iterate rather than recursively call self to save stack space */
            lengths = lengths + (num - r);
            if (bufAddr != NULL)
                bufAddr = bufAddr + (num - r);
            offsets = pn - r;
            num = r;
        }
        else
            break;
    }
}

/*----< heap_merge() >-------------------------------------------------------*/
/* Heapify(a, i, heapsize); Algorithm from Cormen et al. pg. 143 modified for a
 * heap with smallest element at root. The recursion has been removed so that
 * there are no function calls. Function calls are too expensive.
 *
 * Requirement: all individual offsets lists must be already sorted !!!
 */
static
void heap_merge(int              nprocs,
                const MPI_Aint  *count,    /* [nprocs] */
                MPI_Aint         nelems,
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count       *offsets,  /* [nelems] */
                MPI_Count       *blklens,  /* [nelems] */
#else
                MPI_Offset      *offsets,  /* [nelems] */
                int             *blklens,  /* [nelems] */
#endif
                MPI_Aint        *bufAddr)  /* [nelems] */
{
    typedef struct {
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count  *off_list;
        MPI_Count  *len_list;
#else
        MPI_Offset *off_list;
        int        *len_list;
#endif
        MPI_Aint  *addr_list;
        MPI_Aint  count;
    } heap_struct;

    heap_struct *a, tmp;
    int i, j, heapsize, l, r, k, smallest;
    size_t sum;

    /* This heap_merge is not in-place, taking too much memory footprint */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *srt_off = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * nelems);
    MPI_Count *srt_len = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * nelems);
#else
    MPI_Aint *srt_off = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nelems);
    int      *srt_len = (int*)      NCI_Malloc(sizeof(int)      * nelems);
#endif
    MPI_Aint *srt_addr = NULL;

    if (bufAddr != NULL)
        srt_addr = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nelems);

    a = (heap_struct *) NCI_Calloc(nprocs, sizeof(heap_struct));

    /* there are nprocs number of lists to be merged */
    j = 0;
    sum = 0;
    for (i = 0; i < nprocs; i++) {
        if (count[i]) {
            /* each of a[j].off_list is already sorted */
            a[j].off_list = offsets + sum;
            a[j].len_list = blklens + sum;
            if (bufAddr != NULL)
                a[j].addr_list = bufAddr + sum;
            sum += count[i];
            a[j].count = count[i];
            j++;
        }
    }
    nprocs = j; /* some count[i] may be zero */

#define SWAP_HEAP(x, y, tmp) { tmp = x ; x = y ; y = tmp ; }

    heapsize = nprocs;

    /* Build a heap out of the first element from each list, with the smallest
     * element of the heap at the root. The first for loop is to find and move
     * the smallest a[*].off_list[0] to a[0].
     */
    for (i = heapsize / 2 - 1; i >= 0; i--) {
        k = i;
        for (;;) {
            r = 2 * (k + 1);
            l = r - 1;
            if (l < heapsize && a[l].off_list[0] < a[k].off_list[0])
                smallest = l;
            else
                smallest = k;

            if (r < heapsize && a[r].off_list[0] < a[smallest].off_list[0])
                smallest = r;

            if (smallest != k) {
                SWAP_HEAP(a[k], a[smallest], tmp);
                k = smallest;
            } else
                break;
        }
    }

    /* The heap keeps the smallest element in its first element, i.e.
     * a[0].off_list[0].
     */
    j = 0;
    for (i = 0; i < nelems; i++) {
        /* extract smallest element from heap, i.e. the root */
        srt_off[i] = a[0].off_list[0];
        srt_len[i] = a[0].len_list[0];
        if (bufAddr != NULL)
            srt_addr[i] = a[0].addr_list[0];
        a[0].count--;

        if (!a[0].count) {
            a[0] = a[heapsize - 1];
            heapsize--;
        } else {
            a[0].off_list++;
            a[0].len_list++;
            if (bufAddr != NULL)
                a[0].addr_list++;
        }

        /* Heapify(a, 0, heapsize); */
        k = 0;
        for (;;) {
            r = 2 * (k + 1);
            l = r - 1;
            if (l < heapsize && a[l].off_list[0] < a[k].off_list[0])
                smallest = l;
            else
                smallest = k;

            if (r < heapsize && a[r].off_list[0] < a[smallest].off_list[0])
                smallest = r;

            if (smallest != k) {
                SWAP_HEAP(a[k], a[smallest], tmp);
                k = smallest;
            } else
                break;
        }
    }

#ifdef HAVE_MPI_LARGE_COUNT
    memcpy(offsets, srt_off, sizeof(MPI_Count) * nelems);
    memcpy(blklens, srt_len, sizeof(MPI_Count) * nelems);
#else
    memcpy(offsets, srt_off, sizeof(MPI_Offset) * nelems);
    memcpy(blklens, srt_len, sizeof(int)        * nelems);
#endif
    if (bufAddr != NULL)
        memcpy(bufAddr, srt_addr, sizeof(MPI_Aint) * nelems);

    NCI_Free(a);
    if (bufAddr != NULL) NCI_Free(srt_addr);
    NCI_Free(srt_len);
    NCI_Free(srt_off);
}


/*----< ncmpio_ina_init() >--------------------------------------------------*/
/* When intra-node write aggregation is enabled, this subroutine initializes
 * the metadata to be used for intra-node communication and I/O requests.
 *
 * Processes on the same node will first be divided into groups. A process with
 * the lowest rank ID in a group is selected as the aggregator. Only the
 * aggregators call the MPI-IO functions to perform I/O to the file. Thus, this
 * subroutine must be called before MPI_File_open() and should be called only
 * once at ncmpio_create() or ncmpio_open().
 *
 * The subroutine performs the following tasks.
 * 1. Make use of the affinity of each MPI process to its compute node,
 *    represented by ncp->num_nodes and ncp->node_ids[]. These two member of
 *    ncp should have been set from a call to ncmpii_construct_node_list()
 *    earlier during ncmpio_create() and ncmpio_open().
 *    + ncp->num_nodes is the number of unique compute nodes.
 *    + ncp->node_ids[ncp->nprocs] contains node IDs for all processes.
 * 2. Divide processes into groups, select aggregators, and determine whether
 *    self process is an intra-node aggregator.
 *    + ncp->my_aggr is rank ID of my aggregator.
 *    + if (ncp->my_aggr == ncp->rank) then this rank is an aggregator.
 * 3. For an aggregator, find the number of non-aggregators assigned to it and
 *    construct a list of rank IDs of non-aggregators of its group.
 *    + ncp->num_nonaggrs is the number of non-aggregators in its group.
 * 4. For a non-aggregator, find the rank ID of its assigned aggregator.
 *    + ncp->my_aggr is rank ID of my aggregator.
 *    + ncp->nonaggr_ranks[] contains the rank IDs of assigned non-aggregators.
 * 5. Create a new MPI communicator consisting of only the aggregators only.
 *    Obtain the rank ID and total process number of the new communicator.
 *    + ncp->ina_comm contains the aggregators across all nodes.
 *    + ncp->ina_nprocs is the number of processes in intra-node communicator.
 *    + ncp->ina_rank is this process's rank ID in intra-node communicator.
 */
int
ncmpio_ina_init(NC *ncp)
{
    int i, j, mpireturn, do_io, ina_nprocs, naggrs_my_node, first_rank;
    int my_rank_index, *ranks_my_node, my_node_id, nprocs_my_node;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
    int nelems = sizeof(ncp->ina_time_put) / sizeof(ncp->ina_time_put[0]);
    ncp->ina_time_init = ncp->ina_time_flatten = 0.0;
    for (i=0; i<nelems; i++) {
        ncp->ina_time_put[i] = ncp->ina_time_get[i] = 0;
        ncp->maxmem_put[i] = ncp->maxmem_get[i] = 0;
    }
    ncp->ina_npairs_put = ncp->ina_npairs_get = 0;
#endif

    /* initialize parameters of intra-node aggregation */
    ncp->my_aggr = -1;         /* rank ID of my aggregator */
    ncp->num_nonaggrs = 0;     /* number of non-aggregators assigned */
    ncp->nonaggr_ranks = NULL; /* ranks of assigned non-aggregators */

    /* Note that ill value of ncp->num_aggrs_per_node has been checked before
     * entering this subroutine. Thus ncp->num_aggrs_per_node must be > 0.
     */

    /* ncp->node_ids[] has been established in ncmpii_construct_node_list()
     * called in ncmpio_create() or ncmpio_open() before entering this
     * subroutine. my_node_id is this rank's node ID.
     */
    my_node_id = ncp->node_ids[ncp->rank];

    /* nprocs_my_node:  the number of processes in my nodes
     * ranks_my_node[]: rank IDs of all processes in my node.
     * my_rank_index:   points to ranks_my_node[] where
     *                  ranks_my_node[my_rank_index] == ncp->rank
     */
    ranks_my_node = (int*) NCI_Malloc(sizeof(int) * ncp->nprocs);
    my_rank_index = -1;
    nprocs_my_node = 0;
    for (i=0; i<ncp->nprocs; i++) {
        if (ncp->node_ids[i] == my_node_id) {
            if (i == ncp->rank)
                my_rank_index = nprocs_my_node;
            ranks_my_node[nprocs_my_node] = i;
            nprocs_my_node++;
        }
    }
    assert(my_rank_index >= 0);
    /* Now, ranks_my_node[my_rank_index] == ncp->rank */

    /* Make sure number of aggregators in my node <= nprocs_my_node. In some
     * cases, the number of processes allocated to the last few nodes can be
     * less than others.
     */
    naggrs_my_node = MIN(ncp->num_aggrs_per_node, nprocs_my_node);

    /* For each aggregation group, calculate the number of non-aggregators,
     * ncp->num_nonaggrs. Note ncp->num_nonaggrs includes self rank.
     */
    ncp->num_nonaggrs = nprocs_my_node / naggrs_my_node;
    if (nprocs_my_node % naggrs_my_node) ncp->num_nonaggrs++;

    /* Adjust the number of non-aggregators for the last group of each node,
     * to make sure it does not go beyond nprocs_my_node.
     */
    first_rank = my_rank_index - my_rank_index % ncp->num_nonaggrs;
    ncp->num_nonaggrs = MIN(ncp->num_nonaggrs, nprocs_my_node - first_rank);

    /* Assign the first rank as the intra-node aggregator of this group and
     * set the rank ID of my aggregator for each process.
     */
    ncp->my_aggr = ranks_my_node[first_rank];

    if (ncp->num_nonaggrs == 1) {
        /* When the number of processes in this group is 1, the aggregation
         * is not performed. Note num_nonaggrs includes self rank.
         *
         * Note this does not mean intra-node aggregation is disabled. The
         * indicator of whether intra-node aggregation is enabled or disabled
         * is ncp->num_aggrs_per_node, whose value should be consistent across
         * all processes. It is possible for some groups containing only one
         * process, in which the aggregation is not necessarily performed
         * within that group.
         */
        assert(ncp->my_aggr == ncp->rank);
    }
    else if (ncp->my_aggr == ncp->rank) { /* ncp->num_nonaggrs > 1 */
        /* Construct ncp->nonaggr_ranks[], the rank IDs of non-aggregators of
         * this group. Note ncp->nonaggr_ranks[], if malloc-ed, will only be
         * freed when closing the file.
         */
        ncp->nonaggr_ranks = (int*)NCI_Malloc(sizeof(int) * ncp->num_nonaggrs);

        memcpy(ncp->nonaggr_ranks, ranks_my_node + first_rank,
               sizeof(int) * ncp->num_nonaggrs);
    }
    NCI_Free(ranks_my_node);

    /* Next step is to construct a new MPI communicator consisting of all
     * intra-node aggregators. It will later be used to call MPI_File_open(),
     * so that only aggregators call MPI-IO functions to access the file.
     *
     * When using the PnetCDF's internal PNCIO driver, we can pass a list of
     * node_ids of the new communicator to the PNCIO file handler,
     * ncp->pncio_fh, so to prevent the driver from the repeated work of
     * constructing the list of node IDs, node_ids. If using MPI-IO driver,
     * then ROMIO will do this internally again anyway.
     */

    do_io = (ncp->my_aggr == ncp->rank) ? 1 : 0;

    /* construct an array containing ranks of aggregators */
    ncp->ina_node_list = (int*) NCI_Malloc(sizeof(int) * ncp->nprocs);
    TRACE_COMM(MPI_Allgather)(&do_io, 1, MPI_INT, ncp->ina_node_list, 1,
                              MPI_INT,ncp->comm);

    /* Calculate the total number of intra-node aggregators */
    for (ina_nprocs=0, i=0; i<ncp->nprocs; i++)
        if (ncp->ina_node_list[i]) ina_nprocs++;

    /* Construct ncp->node_ids[] and ncp->ina_node_list[]. Their contents
     * depend on the layout of MPI process allocation to the compute nodes.
     * The common layouts can be two kinds:
     *   + cyclic - MPI ranks are assigned to nodes round-robin-ly,
     *   + block - MPI ranks are assigned to a node and then move on to next.
     *
     * Below uses an example of nodes=3, nprocs=10, * num_aggrs_per_node=2.
     * ncp->node_ids[] should be
     *     block  process allocation: 0,0,0,0,1,1,1,2,2,2
     *     cyclic process allocation: 0,1,2,0,1,2,0,1,2,0
     * Accordingly, ncp->ina_node_list[] can be two kinds
     *     block  process allocation: 1,0,1,0,1,0,1,1,0,1
     *     cyclic process allocation: 1,1,1,0,0,0,1,1,1,0
     */

    /* ncp->node_ids[]: node IDs of processes in the new MPI communicator.
     * ncp->ina_node_list[]: the rank IDs of the new MPI communicator.
     */
    for (j=0,i=0; i<ncp->nprocs; i++) {
        if (ncp->ina_node_list[i]) {
            ncp->ina_node_list[j] = i;
            /* Modify ncp->node_ids[] to store the node IDs of the processes in
             * the new communicator. Note ncp->node_ids[] from now on is used
             * by PnetCDF's PNCIO driver only.
             */
            ncp->node_ids[j] = ncp->node_ids[i];
            j++;
        }
    }

    /* Make MPI calls to create a new communicator. */
    MPI_Group origin_group, ina_group;
    TRACE_COMM(MPI_Comm_group)(ncp->comm, &origin_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_group");
    TRACE_COMM(MPI_Group_incl)(origin_group, ina_nprocs, ncp->ina_node_list, &ina_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_incl");
    TRACE_COMM(MPI_Comm_create)(ncp->comm, ina_group, &ncp->ina_comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_create");
    TRACE_COMM(MPI_Group_free)(&ina_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_free");
    TRACE_COMM(MPI_Group_free)(&origin_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_free");

    /* Non-aggregators will have ncp->ina_comm set to MPI_COMM_NULL */
    if (ncp->ina_comm == MPI_COMM_NULL) {
        ncp->ina_nprocs = 0;
        ncp->ina_rank = -1;
    }
    else {
        MPI_Comm_size(ncp->ina_comm, &ncp->ina_nprocs);
        MPI_Comm_rank(ncp->ina_comm, &ncp->ina_rank);
    }

    /* TODO: automatically determine whether or not to enable intra-node
     * aggregation.
     *
     * The ideal case is it can be determined right before each collective
     * write call, because only at that time, the communication pattern is
     * known. If the pattern can cause contention, then enable it. Otherwise,
     * disable it.
     *
     * Such mechanism may depends on the followings.
     *   1. MPI-IO hint cb_noddes, and striping_unit
     *   2. calculate aggregate access region
     *   3. If the number of senders to each cb_nodes is very large, then
     *      intra-node aggregation should be enabled.
     *   4. Average of nprocs_per_node across all processes may be a factor for
     *      determining whether to enable intra-node aggregation. It indicates
     *      whether the high number of processes are allocated on the same
     *      node.
     */

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncp->ina_time_init = MPI_Wtime() - timing;
#endif

    return NC_NOERR;
}

/*----< flatten_subarray() >-------------------------------------------------*/
/* Flatten a subarray request, specified by start[], count[], and stride[] into
 * a list of file offset-length pairs, offsets[] and lengths[].
 */
static int
flatten_subarray(int                ndim,       /* number of dimensions */
                 int                el_size,    /* array element size */
                 MPI_Offset         var_begin,  /* starting file offset */
                 const MPI_Offset  *dimlen,     /* [ndim] dimension lengths */
                 const MPI_Offset  *start,      /* [ndim] starts of subarray */
                 const MPI_Offset  *count,      /* [ndim] counts of subarray */
                 const MPI_Offset  *stride,     /* [ndim] strides of subarray */
                 MPI_Aint          *npairs,     /* OUT: num of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
                 MPI_Count         *offsets,    /* OUT: array of offsets */
                 MPI_Count         *lengths     /* OUT: array of lengths */
#else
                 MPI_Offset        *offsets,    /* OUT: array of offsets */
                 int               *lengths     /* OUT: array of lengths */
#endif
                                     )
{
    int i, j;
    MPI_Offset length, nstride, array_len, off, subarray_len;
    size_t idx=0, idx0;

    *npairs = 0;
    if (ndim < 0) return NC_NOERR;

    if (ndim == 0) {  /* scalar record variable */
        *npairs = 1;
        offsets[0] = var_begin;
        lengths[0] = el_size;
        return NC_NOERR;
    }

    /* TODO: check if all stride[] >= 1
       Q: Is it legal if any stride[] <= 0 ? */

    /* calculate the number of offset-length pairs */
    *npairs = (stride[ndim-1] == 1) ? 1 : count[ndim-1];
    for (i=0; i<ndim-1; i++)
        *npairs *= count[i];
    if (*npairs == 0) /* not reachable, an error if count[] == 0 */
        return NC_NOERR;

    /* length of each row of the subarray are of the same size */
    length  = (stride[ndim-1] == 1) ? count[ndim-1] : 1;
    length *= el_size;
    nstride  = (stride[ndim-1] == 1) ? 1 : count[ndim-1];

    /* set the offset-length pairs for the lowest dimension */
    off = var_begin + start[ndim-1] * el_size;
    for (i=0; i<nstride; i++) {
        offsets[idx]  = off;
        lengths[idx]  = length;
        off          += stride[ndim-1] * el_size;
        idx++;
    }
    ndim--;

    subarray_len = nstride;
    array_len = 1;
    /* for higher dimensions */
    while (ndim > 0) {
        /* array_len is global array size from lowest up to ndim */
        array_len *= dimlen[ndim];

        /* off is the global array offset for this dimension, ndim-1 */
        off = start[ndim-1] * array_len * el_size;

        /* update all offsets from lowest up to dimension ndim-1 */
        idx0 = 0;
        for (j=0; j<subarray_len; j++) {
            offsets[idx0] += off;
            idx0++;
        }

        /* update each plan subarray of dimension ndim-1 */
        off = array_len * stride[ndim-1] * el_size;
        for (i=1; i<count[ndim-1]; i++) {
            idx0 = 0;
            for (j=0; j<subarray_len; j++) {
                offsets[idx] = offsets[idx0] + off;
                lengths[idx] = length;
                idx++;
                idx0++;
            }
            off += array_len * stride[ndim-1] * el_size;
        }
        ndim--;  /* move to next higher dimension */
        subarray_len *= count[ndim];
    }

    /* check if the list can be coalesced */
    for (i=0, j=1; j<*npairs; j++) {
        if (offsets[i] + lengths[i] == offsets[j])
            lengths[i] += lengths[j];
        else {
            i++;
            if (i < j) {
                offsets[i] = offsets[j];
                lengths[i] = lengths[j];
            }
        }
    }
    *npairs = i + 1;

    return NC_NOERR;
}

/*----< flatten_req() >------------------------------------------------------*/
/* Flatten one subarray request into offset-length pairs. Arrays offsets and
 * lengths are allocated in this subroutine and need to be freed by the caller.
 */
static int
flatten_req(NC                *ncp,
            NC_var            *varp,
            const MPI_Offset  *start,
            const MPI_Offset  *count,
            const MPI_Offset  *stride,
            int               *is_incr,   /* OUT: are offsets incrementing */
            MPI_Aint          *num_pairs, /* OUT: number of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count        **off_ptr,   /* OUT: array of flattened offsets */
            MPI_Count        **len_ptr    /* OUT: array of flattened lengths */
#else
            MPI_Offset       **off_ptr,   /* OUT: array of flattened offsets */
            int              **len_ptr    /* OUT: array of flattened lengths */
#endif
                                   )
{
    int i, j, err=NC_NOERR, ndims;
    MPI_Aint num, idx;
    MPI_Offset var_begin, *shape, count0, *ones=NULL;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count  prev_end_off;
    MPI_Count *offsets;
    MPI_Count *lengths;
#else
    MPI_Offset  prev_end_off;
    MPI_Offset *offsets;
    int        *lengths;
#endif

    *num_pairs = 0;    /* total number of offset-length pairs */

    /* Count the number off-len pairs, so we can malloc a contiguous memory
     * space for storing off-len pairs
     */
    if (varp->ndims == 0) { /* scalar variable */
#ifdef HAVE_MPI_LARGE_COUNT
        offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * 2);
        lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * 2);
#else
        offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * 2);
        lengths = (int*)       NCI_Malloc(sizeof(int) * 2);
#endif
        offsets[0] = varp->begin;
        lengths[0] = varp->xsz;
        *num_pairs = 1;
        *off_ptr = offsets;
        *len_ptr = lengths;
        return NC_NOERR;
    }
    else if (varp->ndims == 1 && IS_RECVAR(varp)) { /* scalar variable */
        num = count[0];
    }
    else {
        num = 1;
        if (stride != NULL && stride[varp->ndims-1] > 1)
            num = count[varp->ndims-1];  /* count of last dimension */
        for (i=0; i<varp->ndims-1; i++)
            num *= count[i];       /* all count[] except the last dimension */
    }
    *num_pairs = num;

#ifdef HAVE_MPI_LARGE_COUNT
    offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (num+1));
    lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (num+1));
#else
    offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (num+1));
    lengths = (int*)       NCI_Malloc(sizeof(int)        * (num+1));
#endif
    *off_ptr = offsets;
    *len_ptr = lengths;

    if (stride == NULL) { /* equivalent to {1, 1, ..., 1} */
        ones = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * varp->ndims);
        for (i=0; i<varp->ndims; i++) ones[i] = 1;
    }

    ndims = varp->ndims;
    var_begin = varp->begin;
    shape = varp->shape;
    if (IS_RECVAR(varp)) {
        count0 = count[0];
        var_begin += start[0] * ncp->recsize;
        ndims--;
        start++;
        count++;
        shape++;
        if (stride != NULL) stride++;
    }
    else
        count0 = 1;

    idx = 0;
    *is_incr = 1;
    prev_end_off = -1;
    for (i=0; i<count0; i++) {
        /* flatten the request into a list of offset-length pairs */
        err = flatten_subarray(ndims, varp->xsz, var_begin, shape,
                               start, count, (stride == NULL) ? ones : stride,
                               &num,           /* OUT: num of off-len pairs */
                               offsets + idx,  /* OUT: array of offsets */
                               lengths + idx); /* OUT: array of lengths */

        if (num == 0) continue;

        /* check if offsets[] are in an increasing order */
        for (j=0; j<num; j++) {
            if (prev_end_off > offsets[idx+j])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_end_off = offsets[idx+j];
        }

        idx += num;
        assert(idx <= *num_pairs);

        if (IS_RECVAR(varp))
            var_begin += ncp->recsize;
    }
    if (ones != NULL)
        NCI_Free(ones);

    /* num_pairs may be less than originally calculated, because offset-length
     * pairs are coalesced in the call to flatten_subarray().
     */
    *num_pairs = idx;

    return err;
}

/*----< flatten_reqs() >-----------------------------------------------------*/
/* Flatten multiple subarray requests into file offset-length pairs. Arrays
 * offsets and lengths are allocated here and need to be freed by the caller.
 */
static int
flatten_reqs(NC            *ncp,
             int            reqMode,   /* IN: NC_REQ_RD or NC_REQ_WR */
             int            num_reqs,  /* IN: # requests */
             const NC_req  *reqs,      /* [num_reqs] requests */
             int           *is_incr,   /* OUT: are offsets incrementing */
             MPI_Aint      *num_pairs, /* OUT: total number of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
             MPI_Count    **off_ptr,   /* OUT: array of flattened offsets */
             MPI_Count    **len_ptr    /* OUT: array of flattened lengths */
#else
             MPI_Offset   **off_ptr,   /* OUT: array of flattened offsets */
             int          **len_ptr    /* OUT: array of flattened lengths */
#endif
                                   )
{
    int i, j, status=NC_NOERR, ndims, max_ndims=0;
    MPI_Aint num, idx;
    MPI_Offset *start, *count, *shape, *stride, *ones;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count  prev_end_off;
    MPI_Count *offsets;
    MPI_Count *lengths;
#else
    MPI_Offset  prev_end_off;
    MPI_Offset *offsets;
    int        *lengths;
#endif

    *num_pairs = 0;    /* total number of offset-length pairs */

    /* Count the number off-len pairs from reqs[], so we can malloc a
     * contiguous memory space for storing off-len pairs
     */
    for (i=0; i<num_reqs; i++) {
        /* reqs[i].npairs is the number of offset-length pairs of this request,
         * calculated in ncmpio_igetput_varm() and igetput_varn()
         */
        *num_pairs += reqs[i].npairs;
        if (fIsSet(reqMode, NC_REQ_WR))
            ndims = ncp->put_lead_list[reqs[i].lead_off].varp->ndims;
        else
            ndims = ncp->get_lead_list[reqs[i].lead_off].varp->ndims;
        max_ndims = MAX(max_ndims, ndims);
    }

    /* now we can allocate a contiguous memory space for the off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
    offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (*num_pairs+1));
    lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (*num_pairs+1));
#else
    offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (*num_pairs+1));
    lengths = (int*)       NCI_Malloc(sizeof(int)        * (*num_pairs+1));
#endif
    *off_ptr = offsets;
    *len_ptr = lengths;

    ones = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * max_ndims);
    for (i=0; i<max_ndims; i++) ones[i] = 1;

    idx = 0;
    prev_end_off = -1;
    *is_incr = 1;

    /* now re-run the loop to fill in the off-len pairs */
    for (i=0; i<num_reqs; i++) {
        MPI_Offset var_begin;
        NC_lead_req *lead;
        if (fIsSet(reqMode, NC_REQ_WR))
            lead = ncp->put_lead_list + reqs[i].lead_off;
        else
            lead = ncp->get_lead_list + reqs[i].lead_off;

        if (reqs[i].npairs == 1) {
            /* When reqs[i] contains only one offset-length pair, re-use
             * reqs[i].offset_start, which has been generated earlier at a call
             * to ncmpio_intra_node_aggregation_nreqs().
             */
            offsets[idx] = reqs[i].offset_start;
            lengths[idx] = reqs[i].nelems * lead->varp->xsz;

            /* check if offsets[] are in an increasing order */
            if (prev_end_off > offsets[idx])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_end_off = offsets[idx];
            idx++;
            continue;
        }

        ndims = lead->varp->ndims;
        if (ndims > 0) {
            start  = reqs[i].start;
            count  = start + ndims;
            stride = count + ndims;
        }
        else
            start = count = stride = NULL;

        shape = lead->varp->shape;

        /* find the starting file offset for this variable */
        var_begin = lead->varp->begin;

        /* for record variable, each reqs[] is within a record */
        if (IS_RECVAR(lead->varp)) {
            ndims--;
            start++;
            count++;
            stride++;
            shape++;
            /* find the starting file offset for this record */
            var_begin += reqs[i].start[0] * ncp->recsize;
        }

        if (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) stride = NULL;

        /* flatten each request into a list of offset-length pairs and append
         * to the end of offsets and lengths
         */
        flatten_subarray(ndims, lead->varp->xsz, var_begin, shape,
                         start, count, (stride == NULL) ? ones : stride,
                         &num,           /* OUT: number of off-len pairs */
                         offsets + idx,  /* OUT: array of offsets */
                         lengths + idx); /* OUT: array of lengths */

        /* check if offsets[] are in an increasing order */
        for (j=0; j<num; j++) {
            if (prev_end_off > offsets[idx+j])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_end_off = offsets[idx+j];
        }
        idx += num;
    }
    NCI_Free(ones);

    /* num_pairs may be less than originally calculated, because offset-length
     * pairs are coalesced in the call to flatten_subarray().
     */
    *num_pairs = idx;

    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead;
        if (fIsSet(reqMode, NC_REQ_WR))
            lead = ncp->put_lead_list + reqs[i].lead_off;
        else
            lead = ncp->get_lead_list + reqs[i].lead_off;
        if (fIsSet(lead->flag, NC_REQ_TO_FREE)) {
            NCI_Free(lead->start);
            lead->start = NULL;
        }
    }

    return status;
}

/*----< flat_buf_type() >----------------------------------------------------*/
/* Scan the nonblocking requests, pointed by reqs, and build the offset-length
 * pairs of all buffers, xbuf. Note xbuf in each nonblocking request is a
 * contiguous buffer (packed from the user buffer for the write operations).
 * For record variables, if a user request is accessing more than one record,
 * the request is split into into multiple NC_req objects, one for each record.
 */
static int
flat_buf_type(const NC      *ncp,
              int            reqMode,  /* IN: NC_REQ_RD or NC_REQ_WR */
              int            num_reqs, /* IN: # requests */
              const NC_req  *reqs,     /* IN: [num_reqs] requests */
              PNCIO_View    *buf_view, /* OUT: flattened buftype */
              void         **buf)      /* OUT: pointer to I/O buffer */
/* TODO: */
#if 1
{
    int i, j, err=NC_NOERR;
    NC_lead_req *lead;
    MPI_Aint addr, addr0;
/* buffer offset should be of type MPI_Aint. length should be size_t. */

    buf_view->type = MPI_BYTE;
    buf_view->size = 0;
    buf_view->count = 0;
    buf_view->off = NULL;
    buf_view->len = NULL;
    buf_view->is_contig = 1;

    if (num_reqs == 0)
        return NC_NOERR;

    buf_view->off = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * num_reqs);
#ifdef HAVE_MPI_LARGE_COUNT
    buf_view->len = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * num_reqs);
#else
    buf_view->len = (int*)       NCI_Malloc(sizeof(int)        * num_reqs);
#endif

#if 1
    *buf = reqs[0].xbuf;

    lead = (fIsSet(reqMode, NC_REQ_WR)) ? ncp->put_lead_list
                                        : ncp->get_lead_list;

    MPI_Get_address(lead[reqs[0].lead_off].xbuf, &addr0);
// printf("%s at %d: lead xbuf=%ld nelems=%lld\n",__func__,__LINE__, addr0,lead[reqs[0].lead_off].nelems);

// assert(reqs[0].xbuf == lead[reqs[0].lead_off].xbuf);

    /* set buf_view->off[0] and buf_view->len[0] */
    MPI_Get_address(reqs[0].xbuf, &addr0); /* displacement uses MPI_BOTTOM */
    buf_view->off[0] = 0;

    /* buf_view->len[] are in bytes */
    buf_view->len[0] = reqs[0].nelems * lead[reqs[0].lead_off].varp->xsz;
#if 0
printf("%s at %d: buf_view->len[0]=%lld nelems=%lld\n",__func__,__LINE__, buf_view->len[0],reqs[0].nelems);
j=0;
printf("%s at %d: buf_view xbuf=%ld off[%d]=%lld nelems=%lld\n",__func__,__LINE__, addr0,j,buf_view->off[j],reqs[0].nelems);
#endif


/*
int *wkl, nelems; char *xbuf;
j = 0;
wkl = (int*) malloc(buf_view->len[j]);
nelems=buf_view->len[j]/4; xbuf = (char*)reqs[j].xbuf + buf_view->off[j];
memcpy(wkl, xbuf, nelems*4); ncmpii_in_swapn(wkl, nelems, 4);
printf("%s at %d: nelems=%d off=%lld buf=(%p) ",__func__,__LINE__, nelems, buf_view->off[j], xbuf);
for (i=0; i<nelems; i++) printf(" %d",wkl[i]);
printf("\n");
free(wkl);
*/

    buf_view->size = buf_view->len[0];
    for (i=0, j=1; j<num_reqs; j++) {
        /* displacement uses MPI_BOTTOM */
        MPI_Get_address(reqs[j].xbuf, &addr);
        buf_view->off[j] = addr - addr0;

// printf("%s at %d: buf_view xbuf=%ld off[%d]=%lld nelems=%lld\n",__func__,__LINE__, addr,j,buf_view->off[j],reqs[j].nelems);

// assert(reqs[j].xbuf == lead[reqs[j].lead_off].xbuf);
        /* buf_view->len[] are in bytes */
        buf_view->len[j] = reqs[j].nelems * lead[reqs[j].lead_off].varp->xsz;

/*
wkl = (int*) malloc(buf_view->len[j]);
nelems=buf_view->len[j]/4;
xbuf = (char*)reqs[j].xbuf; // + buf_view->off[j];
xbuf = (char*)(*buf) + buf_view->off[j];
memcpy(wkl, xbuf, nelems*4); ncmpii_in_swapn(wkl, nelems, 4);
printf("%s at %d: nelems=%d off=%lld buf=(%p) ",__func__,__LINE__, nelems, buf_view->off[j], xbuf);
for (i=0; i<nelems; i++) printf(" %d",wkl[i]);
printf("\n");
free(wkl);
*/

        /* accumulate buffer type size */
        buf_view->size += buf_view->len[j];

        /* coalesce the off-len pairs */
        if (buf_view->off[i] + buf_view->len[i] == buf_view->off[j])
            buf_view->len[i] += buf_view->len[j];
        else {
            i++;
            if (i < j) {
                buf_view->off[i] = buf_view->off[j];
                buf_view->len[i] = buf_view->len[j];
            }
        }
    }
    /* After coalescing, the true number of requests may be reduced */
// printf("%s at %d: buf_view->size=%lld\n",__func__,__LINE__, buf_view->size);
#else
    /* set buf_view->off[0] and buf_view->len[0] */
    MPI_Get_address(reqs[0].xbuf, &addr); /* displacement uses MPI_BOTTOM */
    buf_view->off[0] = addr;

    lead = (fIsSet(reqMode, NC_REQ_WR)) ? ncp->put_lead_list
                                        : ncp->get_lead_list;

    /* buf_view->len[] are in bytes */
    buf_view->len[0] = reqs[0].nelems * lead[reqs[0].lead_off].varp->xsz;
 ?  *buf = lead[reqs[0].lead_off].xbuf;

    buf_view->size = buf_view->len[0];
    for (i=0, j=1; j<num_reqs; j++) {
        /* displacement uses MPI_BOTTOM */
        MPI_Get_address(reqs[j].xbuf, &addr);
        buf_view->off[j] = addr;

        /* buf_view->len[] are in bytes */
        buf_view->len[j] = reqs[j].nelems * lead[reqs[j].lead_off].varp->xsz;

        /* accumulate buffer type size */
        buf_view->size += buf_view->len[j];

        /* coalesce the off-len pairs */
        if (buf_view->off[i] + buf_view->len[i] == buf_view->off[j])
            buf_view->len[i] += buf_view->len[j];
        else {
            i++;
            if (i < j) {
                buf_view->off[i] = buf_view->off[j];
                buf_view->len[i] = buf_view->len[j];
            }
        }
    }
    /* After coalescing, the true number of requests may be reduced */
#endif

    if (i + 1 < num_reqs) {
        num_reqs = i + 1; /* num_reqs is reduced */
        buf_view->off = (MPI_Offset*)NCI_Realloc(buf_view->off,
                                       sizeof(MPI_Offset) * num_reqs);
#ifdef HAVE_MPI_LARGE_COUNT
        buf_view->len = (MPI_Offset*)NCI_Realloc(buf_view->len,
                                       sizeof(MPI_Offset) * num_reqs);
#else
        buf_view->len = (int*)       NCI_Realloc(buf_view->len,
                                       sizeof(int)        * num_reqs);
#endif
    }

    buf_view->count = num_reqs;
    buf_view->is_contig = (num_reqs <= 1);

    /* construct buf_view->type if it is noncontiguous */
    if (num_reqs > 1) {
        int mpireturn;
#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_create_hindexed_c(num_reqs, buf_view->len,
                                               buf_view->off, MPI_BYTE,
                                               &buf_view->type);
#else
        MPI_Aint *disp;
#if SIZEOF_MPI_AINT == SIZEOF_MPI_OFFSET
        disp = (MPI_Aint*) buf_view->off;
#else
        disp = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * num_reqs);
        for (j=0; j<num_reqs; j++)
            disp[j] = (MPI_Aint) buf_view->off[j];
#endif

        mpireturn = MPI_Type_create_hindexed(num_reqs, buf_view->len, disp,
                                             MPI_BYTE, &buf_view->type);
#if SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET
        NCI_Free(disp);
#endif
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");

            buf_view->type = MPI_BYTE;
            NCI_Free(buf_view->off);
            NCI_Free(buf_view->len);
            buf_view->off = NULL;
            buf_view->len = NULL;
            buf_view->count = 0;
            buf_view->size = 0;
        }
        else {
            MPI_Type_commit(&buf_view->type);
        }
    }

    return err;
}
#else
{
    int i, j, err, mpireturn, status=NC_NOERR;
    NC_lead_req *lead;
    MPI_Aint addr;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *disps, *blens;
#else
    MPI_Aint *disps;
    int *blens;
#endif

    if (num_reqs == 0) {
        buf_view->type  = MPI_BYTE;
        buf_view->count = 0;
        return NC_NOERR;
    }

#ifdef HAVE_MPI_LARGE_COUNT
    disps = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num_reqs);
    blens = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num_reqs);
#else
    disps = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint)  * num_reqs);
    blens = (int*)      NCI_Malloc(sizeof(int)       * num_reqs);
#endif

    /* set disps[0] and blens[0] */
    MPI_Get_address(reqs[0].xbuf, &addr); /* displacement uses MPI_BOTTOM */
    disps[0] = addr;

    lead = (fIsSet(reqMode, NC_REQ_WR)) ? ncp->put_lead_list
                                        : ncp->get_lead_list;

    /* blens[] are in bytes */
    blens[0] = reqs[0].nelems * lead[reqs[0].lead_off].varp->xsz;
    *buf = lead[reqs[0].lead_off].xbuf;

    for (i=0, j=1; j<num_reqs; j++) {
        /* displacement uses MPI_BOTTOM */
        MPI_Get_address(reqs[j].xbuf, &addr);
        disps[j] = addr;

        /* blens[] are in bytes */
        blens[j] = reqs[j].nelems * lead[reqs[j].lead_off].varp->xsz;

        /* coalesce the disps-blens pairs */
        if (disps[i] + blens[i] == disps[j])
            blens[i] += blens[j];
        else {
            i++;
            if (i < j) {
                disps[i] = disps[j];
                blens[i] = blens[j];
            }
        }
    }

    if (i + 1 < num_reqs) {
        num_reqs = i + 1;
#ifdef HAVE_MPI_LARGE_COUNT
        disps = (MPI_Count*)NCI_Realloc(disps, sizeof(MPI_Count) * num_reqs);
        blens = (MPI_Count*)NCI_Realloc(blens, sizeof(MPI_Count) * num_reqs);
#else
        disps = (MPI_Aint*) NCI_Realloc(disps, sizeof(MPI_Aint)  * num_reqs);
        blens = (int*)      NCI_Realloc(blens, sizeof(int)       * num_reqs);
#endif
    }

    buf_view->count = num_reqs;
    buf_view->off   = disps;
    buf_view->len   = blens;

/* TODO: below datatype construction moves into ncmpio_read_write() */
    if (num_reqs == 1) {
#if 1
buf_view->count = blens[0];
#endif
        buf_view->type = MPI_BYTE;
    }
    else {
#if 1
        /* construct buffer derived datatype */
#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_create_hindexed_c(num_reqs, blens, disps,
                                               MPI_BYTE, &buf_view->type);
#else
        mpireturn = MPI_Type_create_hindexed(num_reqs, blens, disps,
                                             MPI_BYTE, &buf_view->type);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;

            buf_view->type = MPI_BYTE;
            buf_view->count = 0;
        }
        else {
            MPI_Type_commit(&buf_view->type);
buf_view->count = 1;
        }
#endif
        *buf = NULL; /* buf_view->type is constructed using MPI_BOTTOM */
    }

#if 1
    NCI_Free(blens);
    NCI_Free(disps);
#endif
    return status;
}
#endif

/*----< ina_collect_md() >---------------------------------------------------*/
/* Within each intra-node aggregation group, the aggregator collects request
 * metadata from the non-aggregators into meta, including:
 *   1. the number of offset-length pairs on each non-aggregator
 *   2. offsets array of each non-aggregator
 *   3. lengths array of each non-aggregator
 *   4. npairs is the total number of offset-length pairs of this group.
 */
static
int ina_collect_md(NC          *ncp,
                   MPI_Aint    *meta,
#ifdef HAVE_MPI_LARGE_COUNT
                   MPI_Count  **offsets, /* OUT: may be realloc-ed */
                   MPI_Count  **lengths, /* OUT: may be realloc-ed */
#else
                   MPI_Offset **offsets, /* OUT: may be realloc-ed */
                   int        **lengths, /* OUT: may be realloc-ed */
#endif
                   MPI_Aint    *npairs)  /* OUT: total no. off-len pairs */
{
    int i, err, mpireturn, status=NC_NOERR, nreqs;
    MPI_Request *req=NULL;
    MPI_Aint num_pairs=meta[0];

    /* Aggregator collects each non-aggregator's num_pairs and bufLen */
    if (ncp->my_aggr == ncp->rank) {

        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
        nreqs = 0;
        for (i=1; i<ncp->num_nonaggrs; i++)
            TRACE_COMM(MPI_Irecv)(meta + i*3, 3, MPI_AINT,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);

        if (nreqs > 0) {
#ifdef HAVE_MPI_STATUSES_IGNORE
            TRACE_COMM(MPI_Waitall)(nreqs, req, MPI_STATUSES_IGNORE);
#else
            MPI_Status *statuses = (MPI_Status *)
                                   NCI_Malloc(nreqs * sizeof(MPI_Status));
            TRACE_COMM(MPI_Waitall)(nreqs, req, statuses);
            NCI_Free(statuses);
#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
    }
    else /* non-aggregator */
        TRACE_COMM(MPI_Send)(meta, 3, MPI_AINT, ncp->my_aggr, 0, ncp->comm);

    /* Secondly, aggregators collect offset-length pairs from all its
     * non-aggregators
     */
    if (ncp->my_aggr == ncp->rank) {
        MPI_Datatype recvType;

        /* calculate the total number of offset-length pairs to receive */
        for (*npairs=0, i=0; i<ncp->num_nonaggrs; i++) *npairs += meta[i*3];

        /* offsets and lengths have been allocated for storing this rank's
         * offsets and lengths, realloc them to receive offsets and lengths
         * from non-aggregators so they can be in a contiguous buffer.
         */
#ifdef HAVE_MPI_LARGE_COUNT
        if (*npairs > num_pairs) {
            *offsets = (MPI_Count*) NCI_Realloc(*offsets, *npairs * sizeof(MPI_Count));
            *lengths = (MPI_Count*) NCI_Realloc(*lengths, *npairs * sizeof(MPI_Count));
        }
#else
        if (*npairs > num_pairs) {
            /* realloc to store all pairs in a contiguous buffer */
            *offsets = (MPI_Offset*) NCI_Realloc(*offsets, *npairs * sizeof(MPI_Offset));
            *lengths = (int*)        NCI_Realloc(*lengths, *npairs * sizeof(int));
        }
#endif

        /* To minimize number of MPI recv calls per non-aggregator, below
         * creates a derived datatype, recvType, to combine offsets and lengths
         * into one MPI_Irecv call.
         */
        nreqs = 0;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Aint aint;
        MPI_Count bklens[2];
        MPI_Count disps[2];

        MPI_Get_address(*offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Count) * meta[0]);
        MPI_Get_address(*lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(MPI_Count) * meta[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (meta[i*3] == 0) continue;
            bklens[0] = meta[i*3] * sizeof(MPI_Count);
            bklens[1] = meta[i*3] * sizeof(MPI_Count);
            mpireturn = MPI_Type_create_hindexed_c(2, bklens, disps, MPI_BYTE,
                                                   &recvType);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&recvType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
            /* post to receive offset-length pairs from non-aggregators */
            TRACE_COMM(MPI_Irecv_c)(MPI_BOTTOM, 1, recvType,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
            MPI_Type_free(&recvType);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
        }
#else
        int bklens[2];
        MPI_Aint aint, disps[2];

        MPI_Get_address(*offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Offset) * meta[0]);
        MPI_Get_address(*lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(int) * meta[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (meta[i*3] == 0) continue;
            bklens[0] = meta[i*3] * sizeof(MPI_Offset);
            bklens[1] = meta[i*3] * sizeof(int);
            mpireturn = MPI_Type_create_hindexed(2, bklens, disps, MPI_BYTE,
                                                 &recvType);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&recvType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
            /* post to receive offset-length pairs from non-aggregators */
            TRACE_COMM(MPI_Irecv)(MPI_BOTTOM, 1, recvType,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
            MPI_Type_free(&recvType);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
        }
#endif
        if (nreqs > 0) {
#ifdef HAVE_MPI_STATUSES_IGNORE
            TRACE_COMM(MPI_Waitall)(nreqs, req, MPI_STATUSES_IGNORE);
#else
            MPI_Status *statuses = (MPI_Status *)
                                   NCI_Malloc(nreqs * sizeof(MPI_Status));
            TRACE_COMM(MPI_Waitall)(nreqs, req, statuses);
            NCI_Free(statuses);
#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        NCI_Free(req);
    }
    else if (num_pairs > 0) { /* non-aggregator */
        /* To minimize number of MPI send calls to the aggregator, below
         * creates a derived datatype, sendType, to combine offsets and lengths
         * into one MPI_Send call.
         */
        MPI_Datatype sendType;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Aint aint;
        MPI_Count bklens[2];
        MPI_Count disps[2];

        bklens[0] = meta[0] * sizeof(MPI_Count);
        bklens[1] = bklens[0];
        MPI_Get_address(*offsets, &aint);
        disps[0] = aint;
        MPI_Get_address(*lengths, &aint);
        disps[1] = aint;
        mpireturn = MPI_Type_create_hindexed_c(2, bklens, disps, MPI_BYTE,
                                               &sendType);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&sendType);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        TRACE_COMM(MPI_Send_c)(MPI_BOTTOM, 1, sendType, ncp->my_aggr, 0,
                               ncp->comm);
        MPI_Type_free(&sendType);
#else
        int bklens[2];
        MPI_Aint disps[2];

        bklens[0] = meta[0] * sizeof(MPI_Aint);
        bklens[1] = meta[0] * sizeof(int);
        MPI_Get_address(*offsets, &disps[0]);
        MPI_Get_address(*lengths, &disps[1]);
        mpireturn = MPI_Type_create_hindexed(2, bklens, disps, MPI_BYTE,
                                             &sendType);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&sendType);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        TRACE_COMM(MPI_Send)(MPI_BOTTOM, 1, sendType, ncp->my_aggr, 0,
                             ncp->comm);
        MPI_Type_free(&sendType);
#endif
    }

    return status;
}

/*----< ina_put() >----------------------------------------------------------*/
/* This subroutine implements the intra-node aggregation for write operations.
 */
static
int ina_put(NC         *ncp,
            int         is_incr,   /* if offsets are incremental */
            MPI_Aint    num_pairs, /* number of offset-length pairs */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count  *offsets,
            MPI_Count  *lengths,
#else
            MPI_Offset *offsets,
            int        *lengths,
#endif
            PNCIO_View  buf_view,
            void       *buf)       /* user buffer */
{
    int i, j, err, mpireturn, status=NC_NOERR, free_buf_view_off=0;
    char *recv_buf=NULL, *wr_buf = NULL;
    MPI_Aint npairs=0, *meta=NULL, *count=NULL, *bufAddr=NULL;
    MPI_Offset wr_amnt=0;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *off_ptr, *len_ptr;
#else
    MPI_Offset *off_ptr;
    int *len_ptr;
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double endT, startT = MPI_Wtime();
    MPI_Offset mem_max;
    // ncmpi_inq_malloc_size(&mem_max);
    ncmpi_inq_malloc_max_size(&mem_max);
    ncp->maxmem_put[0] = MAX(ncp->maxmem_put[0], mem_max);
#endif

    /* buf may be noncontiguous ! */

    /* Firstly, aggregators collect metadata from non-aggregators.
     *
     * This rank tells its aggregator how much metadata to receive from this
     * rank, by sending: the number of offset-length pairs (num_pairs) and user
     * buffer size in bytes (buf_view.size). This message size to be sent by
     * this rank is 3 MPI_Offset.
     */
    if (ncp->rank == ncp->my_aggr)
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs * 3);
    else
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * 3);

    meta[0] = num_pairs;
    meta[1] = buf_view.size;
    meta[2] = is_incr;

    /* Each aggregator first collects metadata about its offset-length pairs,
     * size of write request, and whether the offsets are in an incremental
     * order. The aggregator will gather these metadata from non-aggregators
     * assigned to it.
     * For write operation, keeping the original offset-length pairs is not
     * necessary, as they will later be sorted and coalesced before calling
     * MPI-IO or PNCIO file write.
     *
     * Once ina_collect_md() returns, this aggregator's offsets and lengths may
     * grow to include the ones from non-aggregators (appended).
     */
    if (ncp->num_nonaggrs > 1) {
        err = ina_collect_md(ncp, meta, &offsets, &lengths, &npairs);
        if (err != NC_NOERR) {
            NCI_Free(meta);
            return err;
        }
    }
    else
        npairs = num_pairs;

    /* For write operation, the non-aggregators now can start sending their
     * write data to the aggregator.
     */
    if (ncp->rank != ncp->my_aggr) { /* non-aggregator */
        if (meta[0] > 0) {
            /* Non-aggregators send write data to the aggregator */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count num = (buf_view.is_contig) ? buf_view.size : 1;
            TRACE_COMM(MPI_Send_c)(buf, num, buf_view.type, ncp->my_aggr,
                                   0, ncp->comm);
#else
            int num = (buf_view.is_contig) ? buf_view.size : 1;
            TRACE_COMM(MPI_Send)(buf, num, buf_view.type, ncp->my_aggr,
                                 0, ncp->comm);
#endif
        }

        /* Must free offsets and lengths now, as they may be realloc-ed in
         * ina_collect_md()
         */
        if (offsets != NULL) NCI_Free(offsets);
        if (lengths != NULL) NCI_Free(lengths);

        /* Non-aggregators are done here, as only aggregators call MPI-IO/PNCIO
         * functions to write data to the file. Non-aggregators do not
         * participate MPI-IO calls.
         */
        NCI_Free(meta);
        return status;
    }

    /* The remaining of this subroutine is for aggregators only */

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    // ncmpi_inq_malloc_size(&mem_max);
    ncmpi_inq_malloc_max_size(&mem_max);
    ncp->maxmem_put[1] = MAX(ncp->maxmem_put[1], mem_max);
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_put[0] += endT - startT;
    startT = endT;
#endif

    off_ptr = offsets;
    len_ptr = lengths;

    /* MPI-IO has the following requirements about filetype.
     * 1. The (flattened) displacements (of a filetype) are not required to be
     *    distinct, but they cannot be negative, and they must be monotonically
     *    non-decreasing.
     * 2. If the file is opened for writing, neither the etype nor the filetype
     *    is permitted to contain overlapping regions.
     */
    if (npairs > 0) {
        /* Now this aggregator has received all offset-length pairs from its
         * non-aggregators. At first, check if a sorting is necessary.
         */
        char *ptr;
        int nreqs, indv_sorted, do_sort, overlap;
        MPI_Request *req=NULL;
        MPI_Offset recv_amnt;

        /* check if offsets of all non-aggregators are individual sorted */
        indv_sorted = 1;
        do_sort = 0;
        for (i=-1,j=0; j<ncp->num_nonaggrs; j++) {
            if (i == -1 && meta[j*3] > 0) /* find 1st whose num_pairs > 0 */
                i = j;
            if (meta[j*3+2] == 0) { /* j's offsets are not sorted */
                indv_sorted = 0;
                do_sort = 1;
                break;
            }
        }
        /* i is the first non-aggregator whose num_pairs > 0, and
         * j is the first non-aggregator whose is_incr is false
         */
// printf("%s at %d: do_sort=%d indv_sorted=%d\n",__func__,__LINE__, do_sort,indv_sorted);

        if (i >= 0 && indv_sorted == 1) {
            /* When all ranks' offsets are individually sorted, we still need
             * to check if offsets are interleaved among all non-aggregators to
             * determine whether a sort for all offset-length pairs is
             * necessary.
             */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count prev_end_off;
#else
            MPI_Offset prev_end_off;
#endif
            assert(meta[i*3+2] == 1);

            MPI_Aint sum = meta[i*3];
            prev_end_off = off_ptr[sum-1]; /* last offset of non-aggregator i */

            /* check if the offsets are interleaved */
            for (++i; i<ncp->num_nonaggrs; i++) {
                if (meta[i*3] == 0) /* zero-sized request */
                    continue;
                assert(meta[i*3+2] == 1);

                if (prev_end_off > off_ptr[sum]) {
                    /* off_ptr[sum] is the non-aggregator i' 1st offset */
                    do_sort = 1; /* offsets are not incrementing */
                    break;
                }
                /* move on to next non-aggregator */
                sum += meta[i*3];
                prev_end_off = off_ptr[sum-1];
            }
        }

        if (do_sort && indv_sorted) {
            /* Interleaved offsets are found but individual offsets are already
             * sorted. This is commonly seen from the checkerboard domain
             * partitioning pattern. In this case, heap_merge() must be called
             * to merge all individually already-sorted offsets into one single
             * sorted offset list. Note count[] is initialized and will be used
             * in heap_merge()
             */
            count = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint)*ncp->num_nonaggrs);
            for (i=0; i<ncp->num_nonaggrs; i++) count[i] = meta[i*3];
        }

        /* Construct an array of buffer addresses containing a mapping of the
         * buffer used to receive write data from non-aggregators and the
         * buffer used to write to file. bufAddr[] is calculated based on the
         * assumption that the write buffer of this aggregator is contiguous,
         * i.e. buf_view.is_contig being 1. For non-aggregators, their write
         * data will always be received into a contiguous buffer.
         */
        bufAddr = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * npairs);
        bufAddr[0] = 0;
        for (i=1; i<npairs; i++)
            bufAddr[i] = bufAddr[i-1] + len_ptr[i-1];

        if (do_sort) {
            /* Sort offsets, lengths, bufAddr altogether, based on offsets into
             * an increasing order.
             */
            if (indv_sorted) {
                /* heap-merge() runs much faster than qsort() when individual
                 * lists have already been sorted. However, it has a much
                 * bigger memory footprint.
                 */
                heap_merge(ncp->num_nonaggrs, count, npairs, off_ptr, len_ptr,
                           bufAddr);
                NCI_Free(count);
            }
            else
                /* When some individual offsets are not sorted, we cannot use
                 * heap_merge(). Note qsort() is an in-place sorting.
                 */
                qsort_off_len_buf(npairs, off_ptr, len_ptr, bufAddr);
        }
// printf("%s at %d: do_sort=%d indv_sorted=%d\n",__func__,__LINE__, do_sort,indv_sorted);

        /* Now off_ptr and len_ptr are sorted, but overlaps may exist between
         * adjacent pairs. If this is the case, they must be coalesced.
         *
         * Below loop checks if there is overlap and calculates recv_amnt and
         * wr_amnt.
         * recv_amnt is the total amount this aggregator will receive from
         *     non-aggregators, including self. recv_amnt includes overlaps.
         * wr_amnt is recv_amnt with overlap removed.
         *
         * This loop also coalesces offset-length pairs as well as the
         * corresponding buffer addresses, so they can be used to move write
         * data around in the true write buffer.
         */
        overlap = 0;
int fake_overlap=0;
        wr_amnt = recv_amnt = len_ptr[0];
        for (i=0, j=1; j<npairs; j++) {
            recv_amnt += len_ptr[j];
            if (off_ptr[i] + len_ptr[i] >= off_ptr[j] + len_ptr[j]) {
                overlap = 1;
fake_overlap=1;
                /* segment i completely covers segment j, skip j */
                continue;
            }

            MPI_Offset gap = off_ptr[i] + len_ptr[i] - off_ptr[j];
            if (gap >= 0) { /* overlap detected, merge j into i */
                /* when gap > 0,  pairs i and j overlap
                 * when gap == 0, pairs i and j are contiguous
                 */
                if (gap > 0) overlap = 1;
if (gap >= 0) fake_overlap=1;
                wr_amnt += len_ptr[j] - gap;
                if (bufAddr[i] + len_ptr[i] == bufAddr[j] + gap) {
                    /* buffers i and j are contiguous, merge j into i */
                    len_ptr[i] += len_ptr[j] - gap;
                }
                else { /* buffers are not contiguous, reduce j's len */
                    off_ptr[i+1] = off_ptr[j] + gap;
                    len_ptr[i+1] = len_ptr[j] - gap;
                    bufAddr[i+1] = bufAddr[j] + gap;
                    i++;
                }
            }
            else { /* i and j do not overlap */
                wr_amnt += len_ptr[j];
                i++;
                if (i < j) {
                    off_ptr[i] = off_ptr[j];
                    len_ptr[i] = len_ptr[j];
                    bufAddr[i] = bufAddr[j];
                }
            }
        }
/*
if (ncp->num_nonaggrs == 1 && do_sort == 1) printf("%s at %d: overlap=%d do_sort=%d after coalesce npairs changed from %ld to %d wr_amnt=%lld recv_amnt=%lld\n",__func__,__LINE__, overlap, do_sort,npairs,i+1,wr_amnt,recv_amnt);
*/

if (fake_overlap == 0) assert(npairs == i+1);

        /* Now off_ptr[], len_ptr[], bufAddr[] are coalesced and no overlap */
        npairs = i+1;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        // ncmpi_inq_malloc_size(&mem_max);
        ncmpi_inq_malloc_max_size(&mem_max);
        ncp->maxmem_put[2] = MAX(ncp->maxmem_put[2], mem_max);

        endT = MPI_Wtime();
        ncp->ina_time_put[1] += endT - startT;
        ncp->ina_npairs_put = MAX(ncp->ina_npairs_put, npairs);
        startT = endT;
#endif

        /* Allocate receive buffer. Once write data from non-aggregators have
         * received into recv_buf, it is packed into wr_buf. Then, wr_buf is
         * used to call MPI-IO/PNCIO file write. Note the wr_buf is always
         * contiguous.
         *
         * When ncp->num_nonaggrs == 1, wr_buf is set to buf which is directly
         * passed to MPI-IO/PNCIO file write.
         *
         * If file offset-length pairs have not been re-ordered, i.e. sorted
         * and overlaps removed, and this aggregator will not receive any write
         * data from its non-aggregators, then we can use user's buffer, buf,
         * to call MPI-IO/PNCIO to write to the file, without allocating an
         * additional temporary buffer.
         */
        if (!do_sort && buf_view.size == recv_amnt && !overlap)
            recv_buf = buf;
        else
            recv_buf = (char*) NCI_Malloc(recv_amnt);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        // ncmpi_inq_malloc_size(&mem_max);
        ncmpi_inq_malloc_max_size(&mem_max);
        ncp->maxmem_put[3] = MAX(ncp->maxmem_put[3], mem_max);
#endif

        if (recv_buf != buf) {
            /* Pack this aggregator's write data into front of recv_buf */
            if (buf_view.is_contig && buf_view.type == MPI_BYTE)
                memcpy(recv_buf, buf, buf_view.size);
            else {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count pos=0;
                MPI_Count num = (buf_view.is_contig) ? buf_view.size : 1;
                MPI_Pack_c(buf, num, buf_view.type, recv_buf, buf_view.size,
                           &pos, MPI_COMM_SELF);
#else
                int pos=0;
                MPI_Count num = (buf_view.is_contig) ? buf_view.size : 1;
                MPI_Pack(buf, num, buf_view.type, recv_buf, buf_view.size,
                         &pos, MPI_COMM_SELF);
#endif
            }
        }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        endT = MPI_Wtime();
        ncp->ina_time_put[2] += endT - startT;
        startT = endT;
#endif

        /* Receive write data sent from non-aggregators. Note we cannot move
         * the posting of MPI_Irecv calls to before sorting and leave
         * MPI_Waitall() to after sorting to overlap communication with the
         * sorting, because the sorting determines the receive buffer size.
         */
        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
        ptr = recv_buf + buf_view.size;
        nreqs = 0;
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (meta[i*3 + 1] == 0) continue;
#ifdef HAVE_MPI_LARGE_COUNT
            TRACE_COMM(MPI_Irecv_c)(ptr, meta[i*3 + 1], MPI_BYTE,
                           ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
#else
            TRACE_COMM(MPI_Irecv)(ptr, meta[i*3 + 1], MPI_BYTE,
                           ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
#endif
            ptr += meta[i*3 + 1];
        }

        if (nreqs > 0) {
#ifdef HAVE_MPI_STATUSES_IGNORE
            TRACE_COMM(MPI_Waitall)(nreqs, req, MPI_STATUSES_IGNORE);
#else
            MPI_Status *statuses = (MPI_Status *)
                                   NCI_Malloc(nreqs * sizeof(MPI_Status));
            TRACE_COMM(MPI_Waitall)(nreqs, req, statuses);
            NCI_Free(statuses);
#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        NCI_Free(req);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_put[3] += endT - startT;
    startT = endT;
#endif

        /* Now all write data has been collected into recv_buf. In case of any
         * overlap, we must coalesce recv_buf into wr_buf using off_ptr[],
         * len_ptr[], and bufAddr[]. For overlapped regions, requests with
         * lower j indices win the writes to the overlapped regions.
         *
         * In case the user buffer, buf, can not be used to write to the file,
         * loop below packs recv_buf, data received from non-aggregators, into
         * wr_buf, a contiguous buffer, wr_buf, which will later be used in a
         * call to MPI-IO/PNCIO file write.
         */
        if (!do_sort && wr_amnt == recv_amnt) {
            wr_buf = recv_buf;

            if (wr_buf != buf) {
                /* If write data has been packed in wr_buf, a contiguous buffer,
                 * update buf_view before passing it to the MPI-IO/PNCIO file
                 * write.
                 */
                buf_view.size = wr_amnt;
                buf_view.type = MPI_BYTE;
                buf_view.is_contig = 1;
            }
            /* else case is when user's buffer, buf, can be used to write */
        }
        else if (buf_view.is_contig && !overlap) {
            /* Note we can reuse bufAddr[] and len_ptr[] as buf_view.off and
             * buf_view.len only when buf_view.is_contig is true, because
             * bufAddr[] is constructed based on the assumption that the write
             * buffer is contiguous.
             */
            wr_buf = recv_buf;
            buf_view.size      = wr_amnt;
            buf_view.type      = MPI_BYTE;
            buf_view.is_contig = (npairs <= 1);
            buf_view.len       = len_ptr;
            buf_view.count     = npairs;
#if SIZEOF_MPI_AINT == SIZEOF_MPI_OFFSET
            buf_view.off = (MPI_Offset*)bufAddr; /* based on recv_buf */
#else
            buf_view.off = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * npairs);
            for (j=0; j<npairs; j++)
                buf_view.off[j] = (MPI_Offset)bufAddr[j];
            free_buf_view_off = 1;
#endif
        }
        else {
            /* do_sort means buffer's offsets and lengths have been moved
             * around in order to make file offset-length pairs monotonically
             * non-decreasing. We need to re-arrange the write buffer
             * accordingly by copying write data into a temporary buffer,
             * wr_buf, and write it to the file. Note when npairs and wr_amnt
             * are large, copying write data into a contiguous buffer can be
             * expensive, making INA cost high. Although it makes the two-phase
             * I/O MPI-IO and PNCIO run faster, this memory copy cost may not
             * be worthy. Besides, the memory footprint high-water mark is
             * doubled.
             */
            wr_buf = NCI_Malloc(wr_amnt);
            ptr = wr_buf;

            for (j=0; j<npairs; j++) {
                memcpy(ptr, recv_buf + bufAddr[j], len_ptr[j]);
                ptr += len_ptr[j];
            }
            /* Write data has been packed in wr_buf, a contiguous buffer,
             * update buf_view before passing it to the MPI-IO/PNCIO file
             * write.
             */
            buf_view.size = wr_amnt;
            buf_view.type = MPI_BYTE;
            buf_view.is_contig = 1;

            if (recv_buf != buf) NCI_Free(recv_buf);
        }
    } /* if (npairs > 0) */

    NCI_Free(meta);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_put[4] += endT - startT;
#endif

    /* set the fileview */
    err = ncmpio_file_set_view(ncp, 0, MPI_BYTE, npairs, off_ptr, len_ptr);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        wr_amnt = 0;
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    // ncmpi_inq_malloc_size(&mem_max);
    ncmpi_inq_malloc_max_size(&mem_max);
    ncp->maxmem_put[4] = MAX(ncp->maxmem_put[4], mem_max);
#endif

    /* carry out write request to file */
    err = ncmpio_read_write(ncp, NC_REQ_WR, 0, buf_view, wr_buf);
    if (status == NC_NOERR) status = err;

    if (free_buf_view_off) NCI_Free(buf_view.off);
    if (wr_buf != buf)  NCI_Free(wr_buf);
    if (bufAddr != NULL) NCI_Free(bufAddr);

    /* Must free offsets and lengths now, as they may be realloc-ed in
     * ina_collect_md()
     */
    if (offsets != NULL) NCI_Free(offsets);
    if (lengths != NULL) NCI_Free(lengths);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    // ncmpi_inq_malloc_size(&mem_max);
    ncmpi_inq_malloc_max_size(&mem_max);
    ncp->maxmem_put[5] = MAX(ncp->maxmem_put[5], mem_max);
#endif

    return status;
}

static
size_t bin_search(
#ifdef HAVE_MPI_LARGE_COUNT
                  MPI_Count key, MPI_Count *base,
#else
                  MPI_Offset key, MPI_Offset *base,
#endif
                  size_t nmemb)
{
    size_t low, high;

    /* only one element */
    if (nmemb == 1)
        return (base[0] <= key) ? 0 : -1;

    /* check the 1st element */
    if (base[0] <= key && key < base[1])
        return 0;

    low = 1;
    high = nmemb - 1;

    while (low <= high) {
        size_t mid = low + (high - low) / 2;
        if (base[mid] == key)
            return mid;
        if (base[mid] < key)
            low = mid + 1;
        else
            high = mid - 1;
    }
    return (low - 1);
}

/*----< ina_get() >----------------------------------------------------------*/
/* This subroutine implements the intra-node aggregation for read operations.
 */
static
int ina_get(NC         *ncp,
            int         is_incr,   /* if offsets are incremental */
            MPI_Aint    num_pairs, /* number of offset-length pairs */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count  *offsets,
            MPI_Count  *lengths,
#else
            MPI_Offset *offsets,
            int        *lengths,
#endif
            PNCIO_View  buf_view,
            void       *buf)      /* user buffer */
{
    int i, j, err, mpireturn, status=NC_NOERR, nreqs;
    int do_sort=0, indv_sorted=1, overlap=0;
    char *rd_buf = NULL;
    MPI_Aint npairs=0, max_npairs, *meta=NULL, *count=NULL;
    MPI_Offset send_amnt=0, rd_amnt=0, off_start;
    MPI_Request *req=NULL;
    PNCIO_View rd_buf_view;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *off_ptr, *len_ptr, *orig_off_ptr, *orig_len_ptr;
    MPI_Count bufLen, *orig_offsets=NULL, *orig_lengths=NULL;
    MPI_Count *blks = NULL, *disps = NULL;
#else
    MPI_Offset *orig_offsets=NULL, *orig_off_ptr, *off_ptr;
    int bufLen, *orig_lengths=NULL, *orig_len_ptr, *len_ptr, *blks = NULL;
    MPI_Aint *disps = NULL;
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double endT, startT = MPI_Wtime();
    MPI_Offset mem_max;
    // ncmpi_inq_malloc_size(&mem_max);
    ncmpi_inq_malloc_max_size(&mem_max);
    ncp->maxmem_get[0] = MAX(ncp->maxmem_get[0], mem_max);
#endif

    bufLen = buf_view.size;

    /* Firstly, aggregators collect metadata from non-aggregators.
     *
     * This rank tells its aggregator how much metadata to receive from this
     * rank, by sending
     *     1. the number of offset-length pairs (num_pairs)
     *     2. user buffer size in bytes (bufLen).
     *     3. whether this rank's offsets are sorted in increasing order.
     * This message size to be sent by this rank is 3 MPI_Offset.
     */
    if (ncp->rank == ncp->my_aggr)
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs * 3);
    else
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * 3);

    meta[0] = num_pairs;
    meta[1] = bufLen;
    meta[2] = is_incr;

    /* Each aggregator first collects metadata about its offset-length pairs,
     * size of read request, and whether the offsets are in an incremental
     * order. The aggregator will gather these metadata from non-aggregators
     * assigned to it.
     *
     * Once ina_collect_md() returns, this aggregator's offsets and lengths may
     * grow to include the ones from non-aggregators (appended).
     */
    if (ncp->num_nonaggrs > 1) {
        err = ina_collect_md(ncp, meta, &offsets, &lengths, &npairs);
        if (err != NC_NOERR) {
            NCI_Free(meta);
            return err;
        }
    }
    else
        npairs = num_pairs;

    if (ncp->rank != ncp->my_aggr) {
        if (meta[0] > 0) {
            /* For read operation, the non-aggregators now can start receiving
             * their read data from the aggregator.
             */
            MPI_Status st;
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count num = (buf_view.is_contig) ? buf_view.size : 1;
            TRACE_COMM(MPI_Recv_c)(buf, num, buf_view.type, ncp->my_aggr, 0,
                                   ncp->comm, &st);
#else
            int num = (buf_view.is_contig) ? buf_view.size : 1;
            TRACE_COMM(MPI_Recv)(buf, num, buf_view.type, ncp->my_aggr, 0,
                                 ncp->comm, &st);
#endif
        }

        /* Must free offsets and lengths now, as they may be realloc-ed in
         * ina_collect_md()
         */
        if (offsets != NULL) NCI_Free(offsets);
        if (lengths != NULL) NCI_Free(lengths);

        /* Non-aggregators are now done, as they do not participate MPI-IO or
         * PNCIO file read.
         */
        NCI_Free(meta);
        return status;
    }

    /* The remaining of this subroutine is for aggregators only. */

    /* For read operation, the original offsets and lengths must be kept
     * untouched, because the later sorting and coalescing will mess up the
     * original order of offsets and lengths, which are needed to construct a
     * datatype when an aggregator sends read data to its non-aggregators.
     */
#ifdef HAVE_MPI_LARGE_COUNT
    orig_offsets = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * npairs);
    orig_lengths = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * npairs);
    memcpy(orig_offsets, offsets, sizeof(MPI_Count) * npairs);
    memcpy(orig_lengths, lengths, sizeof(MPI_Count) * npairs);
#else
    orig_offsets = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * npairs);
    orig_lengths = (int*)        NCI_Malloc(sizeof(int) * npairs);
    memcpy(orig_offsets, offsets, sizeof(MPI_Offset) * npairs);
    memcpy(orig_lengths, lengths, sizeof(int) * npairs);
#endif
    orig_off_ptr = orig_offsets;
    orig_len_ptr = orig_lengths;
    off_ptr      = offsets;
    len_ptr      = lengths;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    // ncmpi_inq_malloc_size(&mem_max);
    ncmpi_inq_malloc_max_size(&mem_max);
    ncp->maxmem_get[1] = MAX(ncp->maxmem_get[1], mem_max);
#endif

    /* MPI-IO has the following requirements about filetype.
     * 1. The (flattened) displacements (of a filetype) are not required to be
     *    distinct, but they cannot be negative, and they must be monotonically
     *    non-decreasing.
     * 2. If the file is opened for writing, neither the etype nor the filetype
     *    is permitted to contain overlapping regions.
     */
    if (npairs > 0) {
        /* Now this aggregator has received all offset-length pairs from its
         * non-aggregators. At first, check if a sorting is necessary.
         */

        /* check if offsets of all non-aggregators are individual sorted */
        indv_sorted = 1;
        for (i=-1,j=0; j<ncp->num_nonaggrs; j++) {
            if (i == -1 && meta[j*3] > 0) /* find 1st whose num_pairs > 0 */
                i = j;
            if (meta[j*3+2] == 0) { /* j's offsets are not sorted */
                indv_sorted = 0;
                do_sort = 1;
                break;
            }
        }
        /* i is the first non-aggregator whose num_pairs > 0
         * j is the first non-aggregator whose is_incr is false
         */

        if (i >= 0 && indv_sorted == 1) {
            /* When all ranks' offsets are individually sorted, we still need
             * to check if offsets are interleaved among all non-aggregators to
             * determine whether a sort for all offset-length pairs is
             * necessary.
             */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count prev_end_off;
#else
            MPI_Offset prev_end_off;
#endif
            assert(meta[i*3+2] == 1);

            MPI_Aint sum = meta[i*3];
            prev_end_off = off_ptr[sum-1]; /* last offset of non-aggregator i */

            /* check if the offsets are interleaved */
            for (++i; i<ncp->num_nonaggrs; i++) {
                if (meta[i*3] == 0) /* zero-sized request */
                    continue;
                assert(meta[i*3+2] == 1);
                if (prev_end_off > off_ptr[sum]) {
                    /* off_ptr[sum] is the non-aggregator i' 1st offset */
                    do_sort = 1; /* offsets are not incrementing */
                    break;
                }
                /* move on to next non-aggregator */
                sum += meta[i*3];
                prev_end_off = off_ptr[sum-1];
            }
        }

        if (do_sort && indv_sorted) {
            /* Interleaved offsets are found but individual offsets are already
             * sorted. In this case, heap_merge() is called to merge all
             * offsets into one single sorted offset list. Note count[] is
             * initialized and will be used in heap_merge()
             */
            count = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint)* ncp->num_nonaggrs);
            for (i=0; i<ncp->num_nonaggrs; i++) count[i] = meta[i*3];
        }

        /* Construct an array of buffer addresses containing a mapping of the
         * buffer used to receive write data from non-aggregators and the
         * buffer used to write to file.
         */
        if (do_sort) {
            /* Sort offsets and lengths, based on offsets into an increasing
             * order.
             */
            if (indv_sorted) {
                /* heap-merge() runs much faster than qsort() when individual
                 * lists have already been sorted. However, it has a much
                 * bigger memory footprint.
                 */
                heap_merge(ncp->num_nonaggrs, count, npairs, off_ptr, len_ptr,
                           NULL);
                NCI_Free(count);
            }
            else
                /* When some individual offsets are not sorted, we cannot use
                 * heap_merge(). Note qsort() is an in-place sorting.
                 */
                qsort_off_len_buf(npairs, off_ptr, len_ptr, NULL);
        }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        // ncmpi_inq_malloc_size(&mem_max);
        ncmpi_inq_malloc_max_size(&mem_max);
        ncp->maxmem_get[2] = MAX(ncp->maxmem_get[2], mem_max);
        ncp->ina_npairs_get = MAX(ncp->ina_npairs_get, npairs);
#endif

        /* Coalesce the offset-length pairs and calculate the total read amount
         * and send amount by this aggregator.
         */
        overlap = 0;
        send_amnt = rd_amnt = len_ptr[0];
        for (i=0, j=1; j<npairs; j++) {
            MPI_Offset gap;
            send_amnt += len_ptr[j];

            gap = off_ptr[i] + len_ptr[i] - off_ptr[j];
            if (gap >= 0) { /* overlap detected, merge j into i */
                /* when gap > 0,  pairs i and j overlap
                 * when gap == 0, pairs i and j are contiguous
                 */
                MPI_Offset i_end, j_end;

                if (gap > 0) overlap = 1;

                i_end = off_ptr[i] + len_ptr[i];
                j_end = off_ptr[j] + len_ptr[j];
                if (i_end < j_end) {
                    len_ptr[i] += j_end - i_end;
                    rd_amnt += j_end - i_end;
                }
                /* else: j is entirely covered by i */
            }
            else { /* j and i are not overlapped */
                rd_amnt += len_ptr[j];
                i++;
                if (i < j) {
                    off_ptr[i] = off_ptr[j];
                    len_ptr[i] = len_ptr[j];
                }
            }
        }

        /* update npairs after coalesce */
        npairs = i+1;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        // ncmpi_inq_malloc_size(&mem_max);
        ncmpi_inq_malloc_max_size(&mem_max);
        ncp->maxmem_get[3] = MAX(ncp->maxmem_get[3], mem_max);
#endif
    } /* if (npairs > 0) */
    /* else case: This aggregation group may not have data to read, but must
     * participate the collective MPI-IO calls.
     */

    /* set the fileview */
    err = ncmpio_file_set_view(ncp, 0, MPI_BYTE, npairs, off_ptr, len_ptr);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        rd_amnt = 0;
    }

    /* Allocate read buffer and send buffer. Once data are read from file into
     * rd_buf, it is unpacked into send_buf for each non-aggregator. send_buf
     * will be directly used to send the read request data to non-aggregators.
     *
     * Note rd_amnt may not be the same as send_amnt, as there can be overlaps
     * between adjacent offset-length pairs after sorted.
     *
     * If file offset-length pairs have not been re-ordered, i.e. sorted and
     * overlaps removed, and this aggregator will not send any read data to its
     * non-aggregators, then we can use user's buffer, buf, to call
     * MPI-IO/PNCIO to read from the file, without allocating an additional
     * temporary buffer.
     */
    if (!do_sort && buf_view.size == send_amnt && !overlap) {
        rd_buf_view = buf_view;
        rd_buf = buf;
    }
    else {
        /* Read data will be stored in a contiguous read buffer. */
        rd_buf_view.size = rd_amnt;
        rd_buf_view.type = MPI_BYTE;
        rd_buf_view.is_contig = 1;
        if (rd_amnt > 0)
            rd_buf = (char*) NCI_Malloc(rd_amnt);
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    // ncmpi_inq_malloc_size(&mem_max);
    ncmpi_inq_malloc_max_size(&mem_max);
    ncp->maxmem_get[4] = MAX(ncp->maxmem_get[4], mem_max);
    endT = MPI_Wtime();
    ncp->ina_time_get[0] += endT - startT;
#endif

    err = ncmpio_read_write(ncp, NC_REQ_RD, 0, rd_buf_view, rd_buf);
    if (status == NC_NOERR) status = err;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    // ncmpi_inq_malloc_size(&mem_max);
    ncmpi_inq_malloc_max_size(&mem_max);
    ncp->maxmem_get[5] = MAX(ncp->maxmem_get[5], mem_max);
    startT = MPI_Wtime();
#endif

    /* If sorting has been performed, the orders of off_ptr[] and len_ptr[] may
     * no longer be the same as the original ones. We must use binary search to
     * find the aggregated offset-length pair containing each non-aggregator's
     * offset-length pair to construct a send buffer datatype, a view layout to
     * the read buffer, rd_buf, so the data can be directly sent from rd_buf.
     */
    if (rd_buf != buf) {
        /* First, aggregators copy the read data to their own user buffer.
         * Note off_ptr[] is sorted in an incremental order.
         *
         * When the offset-length pairs of read buffer have been sorted or
         * the read buffer size is smaller than the total get amount, we must
         * search and copy from read buffer to self's user buffer.
         */
        char *ptr=NULL, *tmp_buf=NULL;
        size_t m=0, k, scan_off=0;

        /* If this aggregator's user buftype is contiguous, the reuse its
         * read buffer. If not, allocate a temporary buffer, copy the read
         * data over, and then unpacking it to the user buffer.
         */
        if (buf_view.is_contig)
            ptr = buf;
        else if (bufLen > 0)
            ptr = tmp_buf = (char*) NCI_Malloc(bufLen);

        for (j=0; j<num_pairs; j++) {
            /* For each offset-length pair j, find the offset-length pair
             * in rd_buf containing it. Note that if the offset-length
             * pairs are not already sorted, i.e. is_incr != 1, this
             * bin_search() below can be very expensive!
             * orig_off_ptr[] and orig_len_ptr[] are the original offsets
             *     and lengths of this rank. We cannot use offsets[] and
             *     lengths[] as they contain offset-length pair of all
             *     non-aggregators and may have been sorted.
             * off_ptr[] and len_ptr[] describe the offset-length pairs of
             *     the read buffer, rd_buf.
             */
            if (!is_incr) m = 0;
            if (npairs-m == 1) assert(off_ptr[m] <= orig_off_ptr[j]);
            k = bin_search(orig_off_ptr[j], &off_ptr[m], npairs-m);
            assert(k < npairs);
            /* k returned from bin_search is relative to m */
            k += m;

            /* When is_incr is 1, the orig_off_ptr[] are in an incremental
             * order and we can continue binary search using the index from
             * previous search. When is_incr is 0, the orig_off_ptr[] are
             * NOT in an incremental order, we must do binary search on the
             * entire off_ptr[].
             */
            if (!is_incr) scan_off = 0;
            for (; m<k; m++)
                scan_off += len_ptr[m];

            /* Note orig_off_ptr[j] and orig_len_ptr[j] must entirely
             * covered by off_ptr[k] and len_ptr[k], because off_ptr[] and
             * len_ptr[] have been coalesced.
             */

            memcpy(ptr,
                   rd_buf + (scan_off + orig_off_ptr[j] - off_ptr[k]),
                   orig_len_ptr[j]);

            ptr += orig_len_ptr[j];
        }

        /* unpack read data to user read buffer, if not done already */
        if (bufLen > 0 && !buf_view.is_contig) {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count pos=0;
            MPI_Unpack_c(tmp_buf, bufLen, &pos, buf, 1, buf_view.type,
                         MPI_COMM_SELF);
#else
            int pos=0;
            MPI_Unpack(tmp_buf, bufLen, &pos, buf, 1, buf_view.type,
                       MPI_COMM_SELF);
#endif
            NCI_Free(tmp_buf);
        }
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    ncp->ina_time_get[1] += endT - startT;
    startT = endT;
#endif

    if (ncp->num_nonaggrs == 1)
        /* In this case, communication will not be necessary. */
        goto fn_exit;

    /* Aggregators start sending read data to non-aggregators. At first,
     * allocate array_of_blocklengths[] and array_of_displacements[]
     */
    for (max_npairs=0, i=1; i<ncp->num_nonaggrs; i++)
        max_npairs = MAX(meta[3*i], max_npairs);

#ifdef HAVE_MPI_LARGE_COUNT
    blks = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * max_npairs);
    disps = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * max_npairs);
#else
    blks = (int*) NCI_Malloc(sizeof(int) * max_npairs);
    disps = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * max_npairs);
#endif

    /* Now, send data to each non-aggregator */
    req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
    nreqs = 0;
    off_start = meta[0];
    for (i=1; i<ncp->num_nonaggrs; i++) {
        /* populate disps[] and blks[] */
        MPI_Aint remote_num_pairs = meta[3*i];
        MPI_Aint remote_is_incr = meta[3*i+2];

        if (remote_num_pairs == 0) continue; /* zero sized request */

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count *off = orig_off_ptr + off_start;
        MPI_Count *len = orig_len_ptr + off_start;
#else
        MPI_Offset *off = orig_off_ptr + off_start;
        int        *len = orig_len_ptr + off_start;
#endif
        size_t k, m = 0;
        size_t scan_off = 0;
        for (j=0; j<remote_num_pairs; j++) {
            MPI_Aint addr;

            /* Find the offset-length pair in rd_buf containing this pair.
             * Note that if the offset-length pairs are not already sorted,
             * i.e. remote_is_incr == 1, this bin_search() below can be very
             * expensive!
             */
            if (!remote_is_incr) m = 0;

            if (npairs-m == 1) assert(off_ptr[m] <= off[j]);
            k = bin_search(off[j], &off_ptr[m], npairs-m);
            /* k returned from bin_search is relative to m */
            k += m;
            assert(off_ptr[k] <= off[j] && off[j] < off_ptr[k] + len_ptr[k]);

            /* When is_incr is 1, the orig_off_ptr[] are in an incremental
             * order, we can continue binary search using the index from the
             * previous search.  When is_incr is 0, the orig_off_ptr[] are NOT
             * in an incremental order, we must do binary search on the entire
             * off_ptr[].
             */
            if (!remote_is_incr) scan_off = 0;
            for (; m<k; m++)
                scan_off += len_ptr[m];
            /* Note orig_off_ptr[j] and len_ptr[j] must entirely covered by
             * off_ptr[k] and len_ptr[k], because off_ptr[] and len_ptr[] have
             * been coalesced.
             */
            char *ptr = rd_buf + (scan_off + off[j] - off_ptr[k]);
            MPI_Get_address(ptr, &addr);
            disps[j] = addr;
            blks[j] = len[j];
        }
        off_start += remote_num_pairs;

        /* Construct a send buffer MPI datatype */
        MPI_Datatype sendType;
#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_create_hindexed_c(remote_num_pairs, blks, disps,
                                               MPI_BYTE, &sendType);
#else
        mpireturn = MPI_Type_create_hindexed(remote_num_pairs, blks, disps,
                                             MPI_BYTE, &sendType);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            MPI_Type_commit(&sendType);

#ifdef HAVE_MPI_LARGE_COUNT
            TRACE_COMM(MPI_Isend_c)(MPI_BOTTOM, 1, sendType,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
#else
            TRACE_COMM(MPI_Isend)(MPI_BOTTOM, 1, sendType,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
#endif
            MPI_Type_free(&sendType);
        }
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    ncp->ina_time_get[2] += endT - startT;
    startT = endT;
#endif

    if (nreqs > 0) {
#ifdef HAVE_MPI_STATUSES_IGNORE
        TRACE_COMM(MPI_Waitall)(nreqs, req, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *)
                               NCI_Malloc(nreqs * sizeof(MPI_Status));
        TRACE_COMM(MPI_Waitall)(nreqs, req, statuses);
        NCI_Free(statuses);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    NCI_Free(blks);
    NCI_Free(disps);

fn_exit:
    /* offsets[] and lengths[] are used in PNCIO read subroutines as flattened
     * filetype. They cannot be freed before the I/O is done.
     */
    if (rd_buf != NULL && rd_buf != buf) NCI_Free(rd_buf);
    if (orig_lengths != NULL) NCI_Free(orig_lengths);
    if (orig_offsets != NULL) NCI_Free(orig_offsets);
    if (req != NULL) NCI_Free(req);
    if (meta != NULL) NCI_Free(meta);

    /* Must free offsets and lengths now, as they may be realloc-ed in
     * ina_collect_md()
     */
    if (offsets != NULL) NCI_Free(offsets);
    if (lengths != NULL) NCI_Free(lengths);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    ncp->ina_time_get[3] += endT - startT;
#endif

    return status;
}

/*----< req_compare() >------------------------------------------------------*/
/* used to sort the the string file offsets of reqs[] */
static int
req_compare(const void *a, const void *b)
{
    if (((NC_req*)a)->offset_start > ((NC_req*)b)->offset_start) return (1);
    if (((NC_req*)a)->offset_start < ((NC_req*)b)->offset_start) return (-1);
    return (0);
}

/*----< ncmpio_ina_nreqs() >-------------------------------------------------*/
/* This subroutine handles PnetCDF's requests made from non-blocking APIs,
 * which contain multiple requests to one or more variable. The input arguments
 * are described below.
 *    reqMode: NC_REQ_RD for read request and NC_REQ_WR for write.
 *    num_reqs: number of elements in array req_list.
 *    req_list[]: stores pending requests from non-blocking API calls, which is
 *                used to construct file offset-length pairs and user buffer
 *                datatype.
 *    newnumrecs: number of new records
 */
int
ncmpio_ina_nreqs(NC         *ncp,
                 int         reqMode,
                 int         num_reqs,
                 NC_req     *req_list,
                 MPI_Offset  newnumrecs)
{
    int err, status=NC_NOERR, is_incr=1;
    void *buf=NULL;
    MPI_Aint num_pairs;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *offsets=NULL, *lengths=NULL;
#else
    MPI_Offset *offsets=NULL;
    int *lengths=NULL;
#endif
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif

// printf("%s at %d: rank=%d num_aggrs_per_nod =%d my_aggr=%d num_nonaggrs=%d\n",__func__,__LINE__, ncp->rank, ncp->num_aggrs_per_node, ncp->my_aggr, ncp->num_nonaggrs);

    /* populate reqs[].offset_start, starting offset of each request */
    NC_req *reqs = req_list;
    int i, descreasing=0;
    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead;
        NC_var *varp;

        lead = (reqMode == NC_REQ_RD) ? ncp->get_lead_list
                                      : ncp->put_lead_list;
        lead += reqs[i].lead_off;
        varp = lead->varp;

        if (varp->ndims == 0) { /* scalar variable */
            reqs[i].offset_start += varp->begin;
        }
        else if (reqs[i].npairs == 1) { /* only one offset-length pair */
            MPI_Offset off = varp->begin;

            if (IS_RECVAR(varp)) off += reqs[i].start[0] * ncp->recsize;

// printf("%s at %d: num_reqs=%d reqs[%d].npairs == 1 offset_start=%lld off=%lld\n", __func__,__LINE__,num_reqs,i,reqs[i].offset_start,off);
            reqs[i].offset_start += off;
        }
        else {
            /* start/count/stride have been allocated in a contiguous array */
            MPI_Offset *count, *stride, offset_end;
            count  = reqs[i].start + varp->ndims;
            stride = (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) ? NULL :
                     count + varp->ndims;

            /* calculate access range of this request */
            ncmpio_calc_start_end(ncp, varp, reqs[i].start, count, stride,
                                  &reqs[i].offset_start, &offset_end);
        }
        /* check if offset_start are in a monotonic nondecreasing order */
        if (i > 0 && reqs[i].offset_start < reqs[i-1].offset_start)
            descreasing = 1;
    }

    /* If a decreasing order is found, sort reqs[] based on reqs[].offset_start
     * into an increasing order.
     */
    if (descreasing)
        qsort(reqs, (size_t)num_reqs, sizeof(NC_req), req_compare);

// printf("%s at %d: descreasing=%d\n",__func__,__LINE__, descreasing);

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     *     is_incr:   whether offsets are incremental
     */
    if (num_reqs > 0)
        flatten_reqs(ncp, reqMode, num_reqs, reqs, &is_incr, &num_pairs,
                     &offsets, &lengths);
    else
        num_pairs = 0;

#if 0
if (0 && num_pairs==10) printf("%s at %d: num_reqs=%d num_pairs=%ld off=%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld len=%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n",__func__,__LINE__, num_reqs, num_pairs,
offsets[0],offsets[1],offsets[2],offsets[3],offsets[4],offsets[5],
offsets[6],offsets[7],offsets[8],offsets[9],
lengths[0],lengths[1],lengths[2],lengths[3],lengths[4],lengths[5],
lengths[6],lengths[7],lengths[8],lengths[9]);

else if (num_pairs==12) printf("%s at %d: num_reqs=%d num_pairs=%ld off=%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld len=%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n",__func__,__LINE__, num_reqs, num_pairs,
offsets[0],offsets[1],offsets[2],offsets[3],offsets[4],
offsets[5],offsets[6],offsets[7],offsets[8],offsets[9],
offsets[10],offsets[11],
lengths[0],lengths[1],lengths[2],lengths[3],lengths[4],
lengths[5],lengths[6],lengths[7],lengths[8],lengths[9],
lengths[10],lengths[11]);
else if (num_pairs) printf("%s at %d: num_reqs=%d num_pairs=%ld off=%lld len=%lld\n",__func__,__LINE__, num_reqs, num_pairs,offsets[0],lengths[0]);
#endif

    /* Populate buf_view, which contains metadata of the user buffers in the
     * nonblocking requests. If buf is non-contiguous, buf to NULL and
     * buf_view.type will be a derived datatype constructed using MPI_BOTTOM.
     */
    PNCIO_View buf_view;
    err = flat_buf_type(ncp, reqMode, num_reqs, reqs, &buf_view, &buf);
    if (status == NC_NOERR) status = err;
if (num_reqs > 0) assert(buf != NULL);

#if 0
if (buf_view.count > 1) printf("%s at %d: buf_view count=%lld off=%lld %lld len=%lld %lld\n",__func__,__LINE__, buf_view.count, buf_view.off[0], buf_view.off[1], buf_view.len[0],buf_view.len[1]);
else if (buf_view.count) printf("%s at %d: buf_view count=%lld off=%lld len=%lld\n",__func__,__LINE__, buf_view.count, buf_view.off[0], buf_view.len[0]);

{int *wkl;
int nelems, j,k, xsz=4;
char *xbuf, msg[1024],str[64];
printf("%s at %d: buf_view count=%lld size=%lld\n",__func__,__LINE__, buf_view.count,buf_view.size);
    wkl = (int*) malloc(buf_view.size);
    nelems=buf_view.size/xsz;
    xbuf = buf;
    memcpy(wkl, xbuf, buf_view.size); ncmpii_in_swapn(wkl, nelems, xsz);
    sprintf(msg,"%s at %d: nelems=%d buf=(%p) ",__func__,__LINE__, nelems, xbuf);
    for (k=0; k<nelems; k++) { sprintf(str," %d",wkl[k]); strcat(msg, str);}
    printf("%s\n",msg);
    free(wkl);
}
#endif

// if (reqMode == NC_REQ_RD) printf("%s at %d: buf_view count=%lld size=%lld is_contig=%d type=%s\n",__func__,__LINE__, buf_view.count, buf_view.size, buf_view.is_contig, (buf_view.type == MPI_BYTE)?"MPI_BYTE":"NOT MPI_BYTE");

    if (req_list != NULL)
        /* All metadata in req_list have been used to construct bufType and
         * bufLen. It is now safe to release the space occupied by req_list.
         */
        NCI_Free(req_list);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_flatten += MPI_Wtime() - timing;
#endif

    int saved_my_aggr, saved_num_nonaggrs;
    saved_my_aggr = ncp->my_aggr;
    saved_num_nonaggrs = ncp->num_nonaggrs;
    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP)) {
        /* Temporarily set ncp->my_aggr and ncp->num_nonaggrs to be as if
         * self rank is an INA aggregator and the INA group size is 1.
         */
        ncp->my_aggr = ncp->rank;
        ncp->num_nonaggrs = 1;
    }

// printf("%s at %d: is_incr=%d buf=%p\n",__func__,__LINE__, is_incr,buf);
    /* perform intra-node aggregation */
    if (fIsSet(reqMode, NC_REQ_WR))
        err = ina_put(ncp, is_incr, num_pairs, offsets, lengths, buf_view, buf);
    else
        err = ina_get(ncp, is_incr, num_pairs, offsets, lengths, buf_view, buf);
    if (status == NC_NOERR) status = err;

    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP)) {
        /* restore ncp->my_aggr and ncp->num_nonaggrs */
        ncp->my_aggr = saved_my_aggr;
        ncp->num_nonaggrs = saved_num_nonaggrs;
    }

#if 0
if (fIsSet(reqMode, NC_REQ_RD))
{int *wkl;
int nelems, j,k, xsz=4;
char *xbuf, msg[1024],str[64];
printf("%s at %d: buf_view count=%lld size=%lld\n",__func__,__LINE__, buf_view.count,buf_view.size);
    wkl = (int*) malloc(buf_view.size);
    nelems=buf_view.size/xsz;
    xbuf = buf;
    memcpy(wkl, xbuf, buf_view.size); ncmpii_in_swapn(wkl, nelems, xsz);
    sprintf(msg,"%s at %d: nelems=%d buf=(%p) ",__func__,__LINE__, nelems, xbuf);
    for (k=0; k<nelems; k++) { sprintf(str," %d",wkl[k]); strcat(msg, str);}
    printf("%s\n",msg);
    free(wkl);
}
#endif

    if (buf_view.type != MPI_BYTE) MPI_Type_free(&buf_view.type);
    if (buf_view.off != NULL) NCI_Free(buf_view.off);
    if (buf_view.len != NULL) NCI_Free(buf_view.len);

    /* Update the number of records if new records have been created.
     * For nonblocking APIs, there is no way for a process to know whether
     * others write to a record variable or not. Note newnumrecs has been
     * sync-ed and always >= ncp->numrecs.
     */
    if (newnumrecs > ncp->numrecs) {
        /* update new record number in file. Note newnumrecs is already
         * sync-ed among all processes and in collective mode
         * ncp->numrecs is always sync-ed in memory among processes,
         * thus no need another MPI_Allreduce to sync it. */
        err = ncmpio_write_numrecs(ncp, newnumrecs);
        if (status == NC_NOERR) status = err;
        /* retain the first error if there is any */
        if (ncp->numrecs < newnumrecs) ncp->numrecs = newnumrecs;
    }

    return status;
}

/*----< ncmpio_ina_req() >---------------------------------------------------*/
/* This subroutine handles a single request made by blocking APIs, involving
 * only one variable. Below describe the subroutine arguments.
 *    reqMode: NC_REQ_RD for read request and NC_REQ_WR for write.
 *    varp: pointer to the variable struct.
 *    start[]: starting offsets
 *    count[]: counts along each dimension
 *    stride[]: stride along each dimension
 *    buf_len: size of I/O buffer in bytes
 *    buf: pointer to the user buffer
 */
int
ncmpio_ina_req(NC               *ncp,
               int               reqMode,
               NC_var           *varp,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               MPI_Offset        buf_len,
               void             *buf)
{
    int err, status=NC_NOERR, is_incr=1;
    MPI_Aint num_pairs;
    PNCIO_View buf_view;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *offsets=NULL, *lengths=NULL;
#else
    MPI_Offset *offsets=NULL;
    int *lengths=NULL;
#endif
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif

    /* blocking API's buffer passed here is always contiguous */
    buf_view.type = MPI_BYTE;
    buf_view.is_contig = 1;
    buf_view.size = buf_len;
    buf_view.count = 0;
    buf_view.off = NULL;
    buf_view.len = NULL;

// printf("%s at %d: buf=%s\n",__func__,__LINE__, (buf==NULL)?"NULL":"NOT NULL");
    if (buf_len == 0 || buf == NULL) {
        /* This is a zero-length request. When in collective data mode, this
         * rank must still participate collective calls. When INA is enabled,
         * this rank tells its aggregator that it has no I/O data. When INA is
         * disabled, this rank must participate other collective file call.
         */
        num_pairs = 0;
        buf_view.size = 0;
        buf_view.count = 0;
    }
    else {
        /* construct file access offset-length pairs
         *     num_pairs: total number of off-len pairs
         *     offsets:   array of flattened offsets
         *     lengths:   array of flattened lengths
         *     is_incr:   whether offsets are incremental
         */
        err = flatten_req(ncp, varp, start, count, stride, &is_incr,
                          &num_pairs, &offsets, &lengths);
        if (err != NC_NOERR) { /* make this rank zero-sized request */
            is_incr = 1;
            num_pairs = 0;
            buf_len = 0;
            buf_view.size = 0;
            buf_view.count = 0;
            if (offsets != NULL) NCI_Free(offsets);
            if (lengths != NULL) NCI_Free(lengths);
            offsets = NULL;
            lengths = NULL;
        }
        status = err;
    }
// if (num_pairs > 0) printf("%s at %d: num_pairs=%ld off=%lld len=%lld\n",__func__,__LINE__, num_pairs,offsets[0],lengths[0]);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (ncp->rank == ncp->my_aggr)
        ncp->ina_time_flatten += MPI_Wtime() - timing;
#endif

    int saved_my_aggr, saved_num_nonaggrs;
    saved_my_aggr = ncp->my_aggr;
    saved_num_nonaggrs = ncp->num_nonaggrs;
    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP)) {
        /* Temporarily set ncp->my_aggr and ncp->num_nonaggrs to be as if
         * self rank is an INA aggregator and the INA group size is 1.
         */
        ncp->my_aggr = ncp->rank;
        ncp->num_nonaggrs = 1;
    }
// if (num_pairs) printf("%s at %d: num_pairs=%ld off=%lld len=%lld\n",__func__,__LINE__, num_pairs,offsets[0],lengths[0]);
// if (buf_view.count) printf("%s at %d: buf_view count=%lld off=%lld len=%lld\n",__func__,__LINE__, buf_view.count, buf_view.off[0], buf_view.len[0]);

// printf("%s at %d: buf_view count=%lld size=%lld is_contig=%d buf=%p\n",__func__,__LINE__, buf_view.count,buf_view.size,buf_view.is_contig,buf);
    /* perform intra-node aggregation */
    if (fIsSet(reqMode, NC_REQ_WR)) {
        err = ina_put(ncp, is_incr, num_pairs, offsets, lengths, buf_view, buf);
        if (status == NC_NOERR) status = err;
    }
    else {
        err = ina_get(ncp, is_incr, num_pairs, offsets, lengths, buf_view, buf);
        if (status == NC_NOERR) status = err;
    }

#if 0
if (fIsSet(reqMode, NC_REQ_RD))
{unsigned long long *wkl; int xsz=8; // int *wkl; int xsz=4;
int nelems, j,k;
char *xbuf, msg[1024],str[64];
printf("%s at %d: buf_view count=%lld size=%lld\n",__func__,__LINE__, buf_view.count,buf_view.size);
    wkl = (unsigned long long*) malloc(buf_view.size); // wkl = (int*) malloc(buf_view.size);
    nelems=buf_view.size/xsz;
    xbuf = buf;
    memcpy(wkl, xbuf, buf_view.size); ncmpii_in_swapn(wkl, nelems, xsz);
    sprintf(msg,"%s at %d: %s nelems=%d buf=(%p) ",__func__,__LINE__, ncp->path,nelems, xbuf);
    // for (k=0; k<nelems; k++) { sprintf(str," %d",wkl[k]); strcat(msg, str);}
    for (k=0; k<nelems; k++) { sprintf(str," %llu",wkl[k]); strcat(msg, str);}
    printf("%s\n",msg);
    free(wkl);
}
#endif

    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP)) {
        /* restore ncp->my_aggr and ncp->num_nonaggrs */
        ncp->my_aggr = saved_my_aggr;
        ncp->num_nonaggrs = saved_num_nonaggrs;
    }

    return status;
}

