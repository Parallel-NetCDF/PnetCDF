/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
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

#ifdef HAVE_MPI_LARGE_COUNT
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    MPI_Count aint; \
    MPI_Count cint; \
    MPI_Count d0 = (x) - offsets; \
    MPI_Count d1 = (y) - offsets; \
    if (d0 != d1) { \
        cint = *(x) ; *(x) = *(y) ; *(y) = cint ; \
        cint = lengths[d0] ; lengths[d0] = lengths[d1] ; lengths[d1] = cint ; \
        aint = bufAddr[d0] ; bufAddr[d0] = bufAddr[d1] ; bufAddr[d1] = aint ; \
    } \
}
#else
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    int int4; \
    MPI_Aint aint; \
    MPI_Aint d0 = (x) - offsets; \
    MPI_Aint d1 = (y) - offsets; \
    if (d0 != d1) { \
        aint = *(x) ; *(x) = *(y) ; *(y) = aint ; \
        int4 = lengths[d0] ; lengths[d0] = lengths[d1] ; lengths[d1] = int4 ; \
        aint = bufAddr[d0] ; bufAddr[d0] = bufAddr[d1] ; bufAddr[d1] = aint ; \
    } \
}
#endif

#define MEDIAN(a,b,c) ((*(a) < *(b)) ? \
                      ((*(b) < *(c)) ? (b) : ((*(a) < *(c)) ? (c) : (a))) : \
                      ((*(b) > *(c)) ? (b) : ((*(a) < *(c)) ? (a) : (c))))

static void
qsort_off_len_buf(MPI_Aint num,
#ifdef HAVE_MPI_LARGE_COUNT
                  MPI_Count *offsets,
                  MPI_Count *lengths,
#else
                  MPI_Aint  *offsets,
                  int       *lengths,
#endif
                  MPI_Aint  *bufAddr)
{
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *pa, *pb, *pc, *pd, *pl, *pm, *pn, cmp_result, swap_cnt;
#else
    MPI_Aint *pa, *pb, *pc, *pd, *pl, *pm, *pn, cmp_result, swap_cnt;
#endif
    MPI_Aint i, r;

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
            /* Iterate rather than recurse to save stack space */
            lengths = lengths + (num - r);
            bufAddr = bufAddr + (num - r);
            offsets = pn - r;
            num = r;
        }
        else
            break;
    }
}

/*----< ncmpio_init_intra_node_aggr() >--------------------------------------*/
/* When intra-node write aggregation is enabled, processes on the same node
 * will be divided into groups. The number of groups is the number of
 * aggregators on that node. The rank IDs of each group must be established.
 *
 * 1. Find information about MPI processes and their affinity to compute node.
 * 2. Determine whether self process is an intra-node aggregator.
 * 3. For an aggregator, find the number of non-aggregators assigned to it and
 *    construct rank IDs of assigned non-aggregators.
 * 4. For a non-aggregator, find the rank ID of its assigned aggregator.
 */
int
ncmpio_intra_node_aggr_init(NC *ncp)
{
    char my_procname[MPI_MAX_PROCESSOR_NAME], **all_procnames=NULL;
    int i, j, k, my_procname_len, num_nodes, root=0;
    int *node_ids=NULL, *all_procname_lens=NULL, *nprocs_per_node;
    int naggrs_my_node, num_nonaggrs;
    int my_rank_index, *ranks_my_node, my_node_id, nprocs_my_node;

    /* initialize parameters of local-node aggregation */
    ncp->my_aggr = -1;         /* rank ID of my aggregator */
    ncp->num_nonaggrs = 0;     /* number of non-aggregators assigned */
    ncp->nonaggr_ranks = NULL; /* ranks of assigned non-aggregators */

    if (ncp->num_aggrs_per_node == 0 || ncp->num_aggrs_per_node == ncp->nprocs)
        /* disable intra-node aggregation */
        return NC_NOERR;

    /* allocate space for storing the rank IDs of non-aggregators assigned to
     * this rank. Note ncp->nonaggr_ranks[] will be freed when closing the
     * file, if allocated.
     */
    num_nonaggrs = ncp->nprocs / ncp->num_aggrs_per_node + 1;
    ncp->nonaggr_ranks = (int*) NCI_Malloc(sizeof(int) * num_nonaggrs);

    /* Collect info about compute nodes in order to select I/O aggregators.
     * Note my_procname is null character terminated, but my_procname_len
     * does not include the null character.
     */
    MPI_Get_processor_name(my_procname, &my_procname_len);
    my_procname_len++; /* to include terminate null character */

    if (ncp->rank == root) {
        /* root collects all procnames */
        all_procnames = (char **) NCI_Malloc(sizeof(char*) * ncp->nprocs);
        if (all_procnames == NULL)
            DEBUG_RETURN_ERROR(NC_ENOMEM)

        all_procname_lens = (int *) NCI_Malloc(sizeof(int) * ncp->nprocs);
        if (all_procname_lens == NULL) {
            NCI_Free(all_procnames);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
    }
    /* gather process name lengths from all processes first */
    MPI_Gather(&my_procname_len, 1, MPI_INT, all_procname_lens, 1, MPI_INT,
               root, ncp->comm);

    if (ncp->rank == root) {
        int *disp;
        size_t alloc_size = 0;

        for (i=0; i<ncp->nprocs; i++)
            alloc_size += all_procname_lens[i];

        all_procnames[0] = (char *) NCI_Malloc(alloc_size);
        if (all_procnames[0] == NULL) {
            NCI_Free(all_procname_lens);
            NCI_Free(all_procnames);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }

        /* Construct displacement array for the MPI_Gatherv, as each process
         * may have a different length for its process name.
         */
        disp = (int *) NCI_Malloc(sizeof(int) * ncp->nprocs);
        disp[0] = 0;
        for (i=1; i<ncp->nprocs; i++) {
            all_procnames[i] = all_procnames[i - 1] + all_procname_lens[i - 1];
            disp[i] = disp[i - 1] + all_procname_lens[i - 1];
        }

        /* gather all process names */
        MPI_Gatherv(my_procname, my_procname_len, MPI_CHAR,
                    all_procnames[0], all_procname_lens, disp, MPI_CHAR,
                    root, ncp->comm);

        NCI_Free(disp);
        NCI_Free(all_procname_lens);
    } else
        /* send process name to root */
        MPI_Gatherv(my_procname, my_procname_len, MPI_CHAR,
                    NULL, NULL, NULL, MPI_CHAR, root, ncp->comm);

    /* each MPI process's compute node ID */
    node_ids = (int *) NCI_Malloc(sizeof(int) * ncp->nprocs);

    if (ncp->rank == root) {
        /* all_procnames[] can tell us the number of nodes and number of
         * processes per node.
         */
        char **node_names;
        int last;

        /* array of pointers pointing to unique host names (compute nodes) */
        node_names = (char **) NCI_Malloc(sizeof(char*) * ncp->nprocs);

        /* number of MPI processes running on each node */
        nprocs_per_node = (int *) NCI_Malloc(sizeof(int) * ncp->nprocs);

        /* calculate nprocs_per_node[] and node_ids[] */
        last = 0;
        num_nodes = 0; /* number of unique compute nodes */
        for (i=0; i<ncp->nprocs; i++) {
            k = last;
            for (j=0; j<num_nodes; j++) {
                /* check if [i] has already appeared in [] */
                if (!strcmp(all_procnames[i], node_names[k])) { /* found */
                    node_ids[i] = k;
                    break;
                }
                k = (k == num_nodes - 1) ? 0 : k + 1;
            }
            if (j < num_nodes)  /* found, next iteration, start with node n */
                last = k;
            else {      /* not found, j == num_nodes, add a new node */
                node_names[j] = strdup(all_procnames[i]);
                nprocs_per_node[j] = 1;
                node_ids[i] = j;
                last = j;
                num_nodes++;
            }
        }
        /* num_nodes is now the number of compute nodes (unique node names) */

        NCI_Free(nprocs_per_node);

        for (i=0; i<num_nodes; i++)
            free(node_names[i]); /* allocated by strdup() */
        NCI_Free(node_names);
        NCI_Free(all_procnames[0]);
        NCI_Free(all_procnames);
    }

    MPI_Bcast(node_ids, ncp->nprocs, MPI_INT, root, ncp->comm);

    /* my_node_id is this rank's node ID */
    my_node_id = node_ids[ncp->rank];

    /* nprocs_my_node: the number of processes in my nodes
     * ranks_my_node[]: rank IDs of all processes in my node.
     * my_rank_index points to ranks_my_node[] where
     * ranks_my_node[my_rank_index] == ncp->rank
     */
    ranks_my_node = (int*) NCI_Malloc(sizeof(int) * ncp->nprocs);
    my_rank_index = -1;
    nprocs_my_node = 0;
    for (i=0; i<ncp->nprocs; i++) {
        if (node_ids[i] == my_node_id) {
            if (i == ncp->rank)
                my_rank_index = nprocs_my_node;
            ranks_my_node[nprocs_my_node] = i;
            nprocs_my_node++;
        }
    }
    assert(my_rank_index >= 0);

    /* Now, ranks_my_node[my_rank_index] == ncp->rank */

    NCI_Free(node_ids);

    /* make sure number of aggregators in my node <= nprocs_my_node */
    naggrs_my_node = MIN(ncp->num_aggrs_per_node, nprocs_my_node);

    /* calculate the number of non-aggregators assigned to an aggregator.
     * Note num_nonaggrs includes self.
     */
    num_nonaggrs = nprocs_my_node / naggrs_my_node;
    if (nprocs_my_node % naggrs_my_node) num_nonaggrs++;

    if (num_nonaggrs == 1)
        /* disable aggregation if the number of non-aggregators assigned to
         * this aggregator is 1. Note num_nonaggrs includes self. It is
         * possible for aggregation enabled or disabled on different nodes and
         * even different aggregation groups on the same node.
         *
         * Use whether ncp->my_aggr < 0 to tell if aggregation is disabled or
         * enabled.
         */
        ncp->my_aggr = -1;
    else {
        /* find the rank ID of aggregator assigned to this rank */
        ncp->my_aggr = ranks_my_node[my_rank_index - my_rank_index % num_nonaggrs];

        if (ncp->my_aggr == ncp->rank) { /* this rank is an aggregator */
            /* Set the number of non-aggregators assigned to this rank. For the
             * last group, make sure it does not go beyond nprocs_my_node.
             */
            ncp->num_nonaggrs = MIN(num_nonaggrs, nprocs_my_node - my_rank_index);
            if (ncp->num_nonaggrs == 1)
                /* disable aggregation, as this aggregation group contains only
                 * self rank
                 */
                ncp->my_aggr = -1;
            else
                /* copy the rank IDs over to ncp->nonaggr_ranks[] */
                memcpy(ncp->nonaggr_ranks,
                       ranks_my_node + my_rank_index,
                       sizeof(int) * num_nonaggrs);
        }
    }
    NCI_Free(ranks_my_node);

    if (ncp->my_aggr < 0) {
        /* free ncp->nonaggr_ranks if aggregation is not enabled */
        NCI_Free(ncp->nonaggr_ranks);
        ncp->nonaggr_ranks = NULL;
    }

    /* TODO: For automatically determine Whether to enable intra-node write
     * aggregation, this should be done right before each collective write
     * call.
     *   1. obtain hint cb_noddes, and striping_unit
     *   2. calculate aggregate access region
     * In each round of two-phase I/O, when the number of senders to each
     * cb_nodes is very large, then intra-node aggregation should be enabled.
     * Average of all nprocs_per_node may be a factor for determining whether
     * to enable intra-node aggregation. It indicates whether the high number
     * of processes are allocated on the same node.
     */

    return NC_NOERR;
}

/*----< flatten_subarray() >-------------------------------------------------*/
/* flatten a subarray request into a list of offset-length pairs */
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
                 MPI_Aint          *offsets,    /* OUT: array of offsets */
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

    return NC_NOERR;
}

/*----< flatten_req() >-----------------------------------------------------*/
/* flatten one write request into offset-length pairs.
 * offsets and lengths are allocated here and need to be freed by the caller
 */
static int
flatten_req(NC                *ncp,
            NC_var            *varp,
            const MPI_Offset  *start,
            const MPI_Offset  *count,
            const MPI_Offset  *stride,
            MPI_Aint          *num_pairs, /* OUT: number of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count        **offsets,   /* OUT: array of flattened offsets */
            MPI_Count        **lengths    /* OUT: array of flattened lengths */
#else
            MPI_Aint         **offsets,   /* OUT: array of flattened offsets */
            int              **lengths    /* OUT: array of flattened lengths */
#endif
                                   )
{
    int j, err=NC_NOERR, ndims;
    MPI_Aint num, idx;
    MPI_Offset var_begin, *shape, count0, *ones=NULL;

    *num_pairs = 0;    /* total number of offset-length pairs */

    /* Count the number off-len pairs, so we can malloc a contiguous memory
     * space for storing off-len pairs
     */
    if (varp->ndims == 0) { /* scalar variable */
#ifdef HAVE_MPI_LARGE_COUNT
        *offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count));
        *lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count));
#else
        *offsets = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint));
        *lengths = (int*)     NCI_Malloc(sizeof(int));
#endif
        (*offsets)[0] = varp->begin;
        (*lengths)[0] = varp->xsz;
        *num_pairs = 1;
        return NC_NOERR;
    }
    else if (varp->ndims == 1 && IS_RECVAR(varp)) { /* scalar variable */
        num = count[0];
    }
    else {
        num = 1;
        if (stride != NULL && stride[varp->ndims-1] > 1)
            num = count[varp->ndims-1];  /* count of last dimension */
        for (j=0; j<varp->ndims-1; j++)
            num *= count[j];       /* all count[] except the last dimension */
    }
    *num_pairs = num;

#ifdef HAVE_MPI_LARGE_COUNT
    *offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num);
    *lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num);
#else
    *offsets = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * num);
    *lengths = (int*)     NCI_Malloc(sizeof(int)      * num);
#endif

    if (stride == NULL) { /* equivalent to {1, 1, ..., 1} */
        ones = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * varp->ndims);
        for (j=0; j<varp->ndims; j++) ones[j] = 1;
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
    for (j=0; j<count0; j++) {
        /* flatten the request into a list of offset-length pairs */
        err = flatten_subarray(ndims, varp->xsz, var_begin, shape,
                               start, count, (stride == NULL) ? ones : stride,
                               &num,            /* OUT: num of off-len pairs */
                               *offsets + idx,  /* OUT: array of offsets */
                               *lengths + idx); /* OUT: array of lengths */
        idx += num;
        assert(idx <= *num_pairs);

        if (IS_RECVAR(varp))
            var_begin += ncp->recsize;
    }
    if (ones != NULL)
        NCI_Free(ones);

    return err;
}

/*----< flatten_reqs() >-----------------------------------------------------*/
/* flatten all write requests into offset-length pairs.
 * offsets and lengths are allocated here and need to be freed by the caller
 */
static int
flatten_reqs(NC            *ncp,
             int            num_reqs,  /* IN: # requests */
             const NC_req  *reqs,      /* [num_reqs] requests */
             MPI_Aint      *num_pairs, /* OUT: total number of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
             MPI_Count    **offsets,   /* OUT: array of flattened offsets */
             MPI_Count    **lengths    /* OUT: array of flattened lengths */
#else
             MPI_Aint     **offsets,   /* OUT: array of flattened offsets */
             int          **lengths    /* OUT: array of flattened lengths */
#endif
                                   )
{
    int i, j, status=NC_NOERR, ndims, max_ndims=0;
    MPI_Aint num, idx;
    MPI_Offset *start, *count, *shape, *stride, *ones;

    *num_pairs = 0;    /* total number of offset-length pairs */

    /* Count the number off-len pairs from reqs[], so we can malloc a
     * contiguous memory space for storing off-len pairs
     */
    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;
        ndims = lead->varp->ndims;
        max_ndims = MAX(max_ndims, ndims);
        if (ndims > 0) {
            start  = reqs[i].start;
            count  = start + ndims;
            stride = count + ndims;
        }
        else
            start = count = stride = NULL;

        /* for record variable, each reqs[] is within a record */
        if (IS_RECVAR(lead->varp)) {
            ndims--;
            start++;
            count++;
            stride++;
        }
        if (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) stride = NULL;

        if (ndims < 0) continue;
        if (ndims == 0) {  /* 1D record variable */
            (*num_pairs)++;
            continue;
        }
        num = 1;
        if (stride != NULL && stride[ndims-1] > 1)
            num = count[ndims-1];  /* count of last dimension */
        for (j=0; j<ndims-1; j++)
            num *= count[j];  /* all count[] except the last dimension */

        (*num_pairs) += num;
    }

    /* now we can allocate a contiguous memory space for the off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
    *offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (*num_pairs));
    *lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (*num_pairs));
#else
    *offsets = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * (*num_pairs));
    *lengths = (int*)     NCI_Malloc(sizeof(int)      * (*num_pairs));
#endif
    idx = 0;

    ones = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * max_ndims);
    for (i=0; i<max_ndims; i++) ones[i] = 1;

    /* now re-run the loop to fill in the off-len pairs */
    for (i=0; i<num_reqs; i++) {
        MPI_Offset var_begin;
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;

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

        /* flatten each request into a list of offset-length pairs and
         * append to the end of offsets and lengths
         */
        flatten_subarray(ndims, lead->varp->xsz, var_begin, shape,
                         start, count, (stride == NULL) ? ones : stride,
                         &num,            /* OUT: number of off-len pairs */
                         *offsets + idx,  /* OUT: array of offsets */
                         *lengths + idx); /* OUT: array of lengths */
        idx += num;
    }
    NCI_Free(ones);

    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;
        if (fIsSet(lead->flag, NC_REQ_TO_FREE)) {
            NCI_Free(lead->start);
            lead->start = NULL;
        }
    }

    return status;
}

/*----< construct_buf_type() >-----------------------------------------------*/
/* construct an MPI derived datatype for I/O buffers from the request list, by
 * concatenate all buffers.
 */
static int
construct_buf_type(const NC     *ncp,
                   int           num_reqs,  /* IN: # requests */
                   const NC_req *reqs,      /* [num_reqs] requests */
                   MPI_Aint     *bufLen,    /* OUT: buffer size in bytes */
                   MPI_Datatype *bufType)   /* OUT: buffer datatype */
{
    int i, err, mpireturn, status=NC_NOERR;
    NC_lead_req *lead;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *blocklens = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num_reqs);
    MPI_Count *disps     = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num_reqs);
#else
    int       *blocklens = (int*)      NCI_Malloc(sizeof(int)       * num_reqs);
    MPI_Aint  *disps     = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint)  * num_reqs);
#endif

    *bufLen = 0;
    for (i=0; i<num_reqs; i++) {
        MPI_Aint addr;

        /* displacement uses MPI_BOTTOM */
        MPI_Get_address(reqs[i].xbuf, &addr);
        disps[i] = addr;

        /* blocklens[] in bytes */
        lead = ncp->put_lead_list + reqs[i].lead_off;
        blocklens[i] = reqs[i].nelems * lead->varp->xsz;

        *bufLen += blocklens[i];
    }

    /* construct buffer derived datatype */
#ifdef HAVE_MPI_LARGE_COUNT
    mpireturn = MPI_Type_create_hindexed_c(num_reqs, blocklens, disps,
                                           MPI_BYTE, bufType);
#else
    mpireturn = MPI_Type_create_hindexed(num_reqs, blocklens, disps,
                                         MPI_BYTE, bufType);
#endif
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");
        /* return the first encountered error if there is any */
        if (status == NC_NOERR) status = err;

        *bufType = MPI_DATATYPE_NULL;
    }
    else {
        MPI_Type_commit(bufType);
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count typeSize;
        MPI_Type_size_c(*bufType, &typeSize);
#else
        int typeSize;
        MPI_Type_size(*bufType, &typeSize);
#endif
        assert(typeSize == *bufLen);
    }

    NCI_Free(blocklens);
    NCI_Free(disps);

    return status;
}

/*----< intra_node_aggregation() >-------------------------------------------*/
/* This is a collective call */
static int
intra_node_aggregation(NC           *ncp,
                       MPI_Aint      num_pairs,
#ifdef HAVE_MPI_LARGE_COUNT
                       MPI_Count    *offsets,
                       MPI_Count    *lengths,
#else
                       MPI_Aint     *offsets,
                       int          *lengths,
#endif
                       MPI_Offset    bufCount,
                       MPI_Datatype  bufType,
                       void         *buf)
{
    int i, j, err, mpireturn, status=NC_NOERR, nreqs;
    char *recv_buf=NULL, *wr_buf = NULL;
    MPI_Aint npairs=0, *msg;
    MPI_Offset offset=0, buf_count;
    MPI_Datatype recvTypes, fileType=MPI_BYTE;
    MPI_File fh;
    MPI_Request *req=NULL;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufLen;
    MPI_Type_size_c(bufType, &bufLen);
#else
    int bufLen;
    MPI_Type_size(bufType, &bufLen);
#endif
    bufLen *= bufCount;

    /* First, tell aggregator how much to receive by sending:
     * (num_pairs and bufLen). The message size to be sent by this rank
     * is num_pairs * 2 * sizeof(MPI_Offset) + bufLen
     */
    if (ncp->rank == ncp->my_aggr)
        msg = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs * 2);
    else
        msg = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * 2);

    msg[0] = num_pairs;
    msg[1] = bufLen;

    /* Aggregator collects each non-aggregator's num_pairs and bufLen */
    if (ncp->rank == ncp->my_aggr) {
        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
        nreqs = 0;
        for (i=1; i<ncp->num_nonaggrs; i++)
            MPI_Irecv(msg + i*2, 2, MPI_AINT, ncp->nonaggr_ranks[i], 0,
                      ncp->comm, &req[nreqs++]);

        mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    else { /* non-aggregator */
        MPI_Send(msg, 2, MPI_AINT, ncp->my_aggr, 0, ncp->comm);
        if (num_pairs == 0)
            NCI_Free(msg);
    }

    /* Aggregator collects offset-length pairs from non-aggregators */
    if (ncp->rank == ncp->my_aggr) {
        /* calculate the total number of offset-length pairs */
        npairs = num_pairs;
        for (i=1; i<ncp->num_nonaggrs; i++) npairs += msg[i*2];

#ifdef HAVE_MPI_LARGE_COUNT
        if (npairs > num_pairs) {
            /* realloc to store all pairs in a contiguous buffer */
            offsets = (MPI_Count*) NCI_Realloc(offsets, sizeof(MPI_Count) * npairs);
            lengths = (MPI_Count*) NCI_Realloc(lengths, sizeof(MPI_Count) * npairs);
        }
#else
        if (npairs > num_pairs) {
            /* realloc to store all pairs in a contiguous buffer */
            offsets = (MPI_Aint*) NCI_Realloc(offsets, sizeof(MPI_Aint) * npairs);
            lengths = (int*) NCI_Realloc(lengths, sizeof(int) * npairs);
        }
#endif

        nreqs = 0;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Aint aint;
        MPI_Count bklens[2];
        MPI_Count disps[2];

        MPI_Get_address(offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Count) * msg[0]);
        MPI_Get_address(lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(MPI_Count) * msg[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (msg[i*2] == 0) continue;
            bklens[0] = msg[i*2] * sizeof(MPI_Count);
            bklens[1] = msg[i*2] * sizeof(MPI_Count);
            mpireturn = MPI_Type_create_hindexed_c(2, bklens, disps, MPI_BYTE,
                                                   &recvTypes);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&recvTypes);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
            /* post to receive offset-length pairs from non-aggregators */
            MPI_Irecv_c(MPI_BOTTOM, 1, recvTypes, ncp->nonaggr_ranks[i],
                        0, ncp->comm, &req[nreqs]);
            MPI_Type_free(&recvTypes);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
            nreqs++;
        }
#else
        int bklens[2];
        MPI_Aint aint, disps[2];

        MPI_Get_address(offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Aint) * msg[0]);
        MPI_Get_address(lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(int) * msg[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (msg[i*2] == 0) continue;
            bklens[0] = msg[i*2] * sizeof(MPI_Aint);
            bklens[1] = msg[i*2] * sizeof(int);
            mpireturn = MPI_Type_create_hindexed(2, bklens, disps, MPI_BYTE,
                                                 &recvTypes);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&recvTypes);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
            /* post to receive offset-length pairs from non-aggregators */
            MPI_Irecv(MPI_BOTTOM, 1, recvTypes, ncp->nonaggr_ranks[i],
                      0, ncp->comm, &req[nreqs]);
            MPI_Type_free(&recvTypes);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
            nreqs++;
        }
#endif
        mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    else if (num_pairs > 0) { /* non-aggregator */
        /* send offset-length pairs data to the aggregator */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Aint aint;
        MPI_Count bklens[2];
        MPI_Count disps[2];

        bklens[0] = msg[0] * sizeof(MPI_Count);
        bklens[1] = bklens[0];
        MPI_Get_address(offsets, &aint);
        disps[0] = aint;
        MPI_Get_address(lengths, &aint);
        disps[1] = aint;
        mpireturn = MPI_Type_create_hindexed_c(2, bklens, disps, MPI_BYTE,
                                               &recvTypes);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&recvTypes);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        MPI_Send_c(MPI_BOTTOM, 1, recvTypes, ncp->my_aggr, 0, ncp->comm);
        MPI_Type_free(&recvTypes);
#else
        int bklens[2];
        MPI_Aint disps[2];

        bklens[0] = msg[0] * sizeof(MPI_Aint);
        bklens[1] = msg[0] * sizeof(int);
        MPI_Get_address(offsets, &disps[0]);
        MPI_Get_address(lengths, &disps[1]);
        mpireturn = MPI_Type_create_hindexed(2, bklens, disps, MPI_BYTE,
                                             &recvTypes);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&recvTypes);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        MPI_Send(MPI_BOTTOM, 1, recvTypes, ncp->my_aggr, 0, ncp->comm);
        MPI_Type_free(&recvTypes);
#endif
        NCI_Free(msg);
    }

    /*
     * TODO, define a datatype to combine sends of offset-length pairs with the
     * write data into a single send call.
     */
    nreqs = 0;
    if (ncp->rank == ncp->my_aggr) {
        /* calculate the total write account */
        buf_count = bufLen;
        for (i=1; i<ncp->num_nonaggrs; i++) buf_count += msg[i*2 + 1];

        /* Allocate receive buffer, which will be sorted into an increasing
         * order based on the file offsets. Thus, after sorting pack recv_buf
         * to wr_buf to avoid creating another buffer datatype.
         */
        if (buf_count > 0) {
            recv_buf = (char*) NCI_Malloc(buf_count);
            wr_buf = (char*) NCI_Malloc(buf_count);
        }

        /* First, pack self write data into front of the recv_buf */
        if (bufLen > 0) {
            if (bufType == MPI_BYTE)
                memcpy(recv_buf, buf, bufLen);
            else {
                void *inbuf = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count position=0;
                MPI_Count incount = (buf == NULL) ? 1 : bufCount;
                MPI_Pack_c(inbuf, incount, bufType, recv_buf, bufLen, &position,
                           MPI_COMM_SELF);
#else
                int position=0;
                int incount = (buf == NULL) ? 1 : bufCount;
                MPI_Pack(inbuf, incount, bufType, recv_buf, bufLen, &position,
                         MPI_COMM_SELF);
#endif
            }
        }

        /* post requests to receive write data from non-aggregators */
        if (buf_count > 0) {
            char *ptr = recv_buf + bufLen;
            for (i=1; i<ncp->num_nonaggrs; i++) {
                if (msg[i*2 + 1] == 0) continue;
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(ptr, msg[i*2 + 1], MPI_BYTE, ncp->nonaggr_ranks[i],
                            0, ncp->comm, &req[nreqs++]);
#else
                MPI_Irecv(ptr, msg[i*2 + 1], MPI_BYTE, ncp->nonaggr_ranks[i],
                          0, ncp->comm, &req[nreqs++]);
#endif
                ptr += msg[i*2 + 1];
            }
            mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        NCI_Free(req);
        NCI_Free(msg);
    }
    else if (bufLen > 0) {
        /* send write data to the aggregator */
        void *buf_ptr = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count num = (buf == NULL) ? 1 : bufCount;
        MPI_Send_c(buf_ptr, num, bufType, ncp->my_aggr, 0, ncp->comm);
#else
        int num = (buf == NULL) ? 1 : bufCount;
        MPI_Send(buf_ptr, num, bufType, ncp->my_aggr, 0, ncp->comm);
#endif
        NCI_Free(offsets);
        NCI_Free(lengths);
    }

    /* aggregator sorts the offset-length pairs, along with the buffer */
    if (ncp->rank == ncp->my_aggr && npairs > 0) {

        /* construct array of buffer addresses */
        MPI_Aint *bufAddr = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * npairs);
        bufAddr[0] = 0;
        for (i=1; i<npairs; i++)
            bufAddr[i] = bufAddr[i-1] + lengths[i-1];

        /* sort offsets, lengths, bufAddr altogether, based on offsets into
         * an increasing order
         */
        qsort_off_len_buf(npairs, offsets, lengths, bufAddr);

        /* merge the overlapped buffer segments, skip the overlapped regions
         * for those with higher j indices (i.e. requests with lower j indices
         * win the writes to the overlapped regions)
         */
        for (i=0, j=1; j<npairs; j++) {
            if (offsets[i] + lengths[i] >= offsets[j] + lengths[j])
                /* segment i completely covers segment j, skip j */
                continue;

            MPI_Offset gap = offsets[i] + lengths[i] - offsets[j];
            if (gap >= 0) { /* segments i and j overlaps */
                if (bufAddr[i] + lengths[i] == bufAddr[j] + gap) {
                    /* buffers i and j are contiguous, merge j to i */
                    lengths[i] = MPI_Aint_add(lengths[i], lengths[j] - gap);
                }
                else { /* buffers are not contiguous, reduce j's len */
                    offsets[i+1] = offsets[j] + gap;
                    lengths[i+1] = lengths[j] - gap;
                    bufAddr[i+1] = bufAddr[j] + gap;
                    i++;
                }
            }
            else { /* i and j do not overlap */
                i++;
                if (i < j) {
                    offsets[i] = offsets[j];
                    lengths[i] = lengths[j];
                    bufAddr[i] = bufAddr[j];
                }
            }
        }
        /* update number of pairs, now all off-len pairs are not overlapped */
        npairs = i+1;

        /* pack recv_buf, data received from non-aggregators, into wr_buf, a
         * contiguous buffer, wr_buf, which will later be used in a call to
         * MPI_File_write_all()
         */
        char *ptr = wr_buf;
        buf_count = 0;
        if (npairs > 0) {
            memcpy(ptr, recv_buf + bufAddr[0], lengths[0]);
            ptr += lengths[0];
            buf_count = lengths[0];
        }
        for (i=0, j=1; j<npairs; j++) {
            memcpy(ptr, recv_buf + bufAddr[j], lengths[j]);
            ptr += lengths[j];
            /* overlap may be found, recalculate buf_count */
            buf_count += lengths[j];

            /* coalesce the offset-length pairs */
            if (offsets[i] + lengths[i] == offsets[j]) {
                /* coalesce j into i */
                lengths[i] += lengths[j];
            }
            else {
                i++;
                if (i < j) {
                    offsets[i] = offsets[j];
                    lengths[i] = lengths[j];
                }
            }
        }
        NCI_Free(bufAddr);
        if (recv_buf != NULL) NCI_Free(recv_buf);

        /* update number of pairs, now all off-len pairs are not overlapped */
        npairs = i+1;

        if (npairs == 1) {
            /* No need to create fileType if writing to a contiguous space */
            offset = offsets[0];
        }
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            /* construct fileview */
            mpireturn = MPI_Type_create_hindexed_c(npairs, lengths, offsets,
                                                   MPI_BYTE, &fileType);

#else
            /* construct fileview */
            mpireturn = MPI_Type_create_hindexed(npairs, lengths, offsets,
                                                 MPI_BYTE, &fileType);

#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&fileType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
        }
        NCI_Free(offsets);
        NCI_Free(lengths);
    }

    if (ncp->rank != ncp->my_aggr) /* non-aggregator writes nothing */
        buf_count = 0;

    /* Only aggregators writes non-zero sized of data to the file. The
     * non-aggregators participate the collective write call with zero-length
     * write requests.
     */
    fh = ncp->collective_fh;

    /* set the MPI-IO fileview, this is a collective call */
    err = ncmpio_file_set_view(ncp, fh, &offset, fileType);
    if (fileType != MPI_BYTE) MPI_Type_free(&fileType);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        buf_count = 0;
    }

    /* call MPI_File_write_at_all */
    err = ncmpio_read_write(ncp, NC_REQ_WR, NC_REQ_COLL, offset, buf_count,
                            MPI_BYTE, wr_buf, 1);
    if (status == NC_NOERR) status = err;

    if (wr_buf != NULL) NCI_Free(wr_buf);

    return status;
}

/*----< ncmpio_intra_node_aggregation_nreqs() >------------------------------*/
/* This is a collective call */
int
ncmpio_intra_node_aggregation_nreqs(NC         *ncp,
                                    int         num_reqs,
                                    NC_req     *put_list,
                                    MPI_Offset  newnumrecs)
{
    int err, status=NC_NOERR;
    MPI_Aint bufLen, num_pairs;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *offsets=NULL, *lengths=NULL;
#else
    MPI_Aint *offsets=NULL;
    int *lengths=NULL;
#endif
    MPI_Datatype bufType=MPI_BYTE;

    assert(ncp->my_aggr >= 0);

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     */
    if (num_reqs > 0)
        flatten_reqs(ncp, num_reqs, put_list, &num_pairs, &offsets, &lengths);
    else
        num_pairs = 0;

    /* construct write buffer datatype, bufType.
     * bufLen is the buffer size in bytes
     */
    if (num_reqs > 0) {
        construct_buf_type(ncp, num_reqs, put_list, &bufLen, &bufType);
        bufLen = 1;
    }
    else
        bufLen = 0;

    if (put_list != NULL)
        NCI_Free(put_list);

    err = intra_node_aggregation(ncp, num_pairs, offsets, lengths, bufLen,
                                 bufType, NULL);
    if (status == NC_NOERR) status = err;

    /* free and reset bufType */
    if (bufType != MPI_BYTE && bufType != MPI_DATATYPE_NULL)
        MPI_Type_free(&bufType);

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

/*----< ncmpio_intra_node_aggregation() >------------------------------------*/
/* This is a collective call */
int
ncmpio_intra_node_aggregation(NC               *ncp,
                              NC_var           *varp,
                              const MPI_Offset *start,
                              const MPI_Offset *count,
                              const MPI_Offset *stride,
                              MPI_Offset        bufCount,
                              MPI_Datatype      bufType,
                              void             *buf)
{
    int err, status=NC_NOERR;
    MPI_Aint num_pairs;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *offsets=NULL, *lengths=NULL;
#else
    MPI_Aint *offsets=NULL;
    int *lengths=NULL;
#endif

    if (buf == NULL) /* zero-length request */
        return intra_node_aggregation(ncp, 0, NULL, NULL, 0, MPI_BYTE, NULL);

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     */
    err = flatten_req(ncp, varp, start, count, stride, &num_pairs, &offsets,
                      &lengths);
    if (err != NC_NOERR) {
        num_pairs = 0;
        if (offsets != NULL)
            NCI_Free(offsets);
        offsets = NULL;
    }
    status = err;

    err = intra_node_aggregation(ncp, num_pairs, offsets, lengths, bufCount,
                                 bufType, buf);
    if (status == NC_NOERR) status = err;

    return status;
}

