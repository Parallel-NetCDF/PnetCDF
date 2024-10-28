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
    int i, j, k, rank, nprocs, my_procname_len, num_nodes, root=0;
    int *node_ids=NULL, *all_procname_lens=NULL, *nprocs_per_node;
    int naggrs_my_node, num_nonaggrs;
    int my_rank_index, *ranks_my_node, my_node_id, nprocs_my_node;

    /* initialize parameters of local-node aggregation */
    ncp->aggregation = 0;      /* is intra-node aggregation enabled? */
    ncp->my_aggr = 0;          /* rank ID of my aggregator */
    ncp->num_nonaggrs = 0;     /* number of non-aggregators assigned */
    ncp->nonaggr_ranks = NULL; /* ranks of assigned non-aggregators */

    if (ncp->num_aggrs_per_node == 0)
        return NC_NOERR;

    MPI_Comm_size(ncp->comm, &nprocs);
    MPI_Comm_rank(ncp->comm, &rank);

    /* allocate space for storing the rank IDs of non-aggregators assigned to
     * this rank. Note ncp->nonaggr_ranks[] will be freed when closing the
     * file, if allocated.
     */
    num_nonaggrs = nprocs / ncp->num_aggrs_per_node + 1;
    ncp->nonaggr_ranks = (int*) NCI_Malloc(sizeof(int) * num_nonaggrs);
    ncp->aggregation = 1;

    /* Collect info about compute nodes in order to select I/O aggregators.
     * Note my_procname is null character terminated, but my_procname_len
     * does not include the null character.
     */
    MPI_Get_processor_name(my_procname, &my_procname_len);
    my_procname_len++; /* to include terminate null character */

/*
printf("%d: ---------  my_procname=%s\n",rank,my_procname);
if (rank < (nprocs/2)+(nprocs%2)) {sprintf(my_procname,"dummy.0"); my_procname_len=strlen(my_procname)+1;}
else {sprintf(my_procname,"dummy.1"); my_procname_len=strlen(my_procname)+1;}
*/
#ifdef WKL_DEBUG
#endif

    if (rank == root) {
        /* root collects all procnames */
        all_procnames = (char **) NCI_Malloc(sizeof(char*) * nprocs);
        if (all_procnames == NULL)
            DEBUG_RETURN_ERROR(NC_ENOMEM)

        all_procname_lens = (int *) NCI_Malloc(sizeof(int) * nprocs);
        if (all_procname_lens == NULL) {
            NCI_Free(all_procnames);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
    }
    /* gather process name lengths from all processes first */
    MPI_Gather(&my_procname_len, 1, MPI_INT, all_procname_lens, 1, MPI_INT,
               root, ncp->comm);

    if (rank == root) {
        int *disp;
        size_t alloc_size = 0;

        for (i=0; i<nprocs; i++)
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
        disp = (int *) NCI_Malloc(sizeof(int) * nprocs);
        disp[0] = 0;
        for (i=1; i<nprocs; i++) {
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
    node_ids = (int *) NCI_Malloc(sizeof(int) * nprocs);

    if (rank == root) {
        /* all_procnames[] can tell us the number of nodes and number of
         * processes per node.
         */
        char **node_names;
        int last;

        /* array of pointers pointing to unique host names (compute nodes) */
        node_names = (char **) NCI_Malloc(sizeof(char*) * nprocs);

        /* number of MPI processes running on each node */
        nprocs_per_node = (int *) NCI_Malloc(sizeof(int) * nprocs);

        /* calculate nprocs_per_node[] and node_ids[] */
        last = 0;
        num_nodes = 0; /* number of unique compute nodes */
        for (i=0; i<nprocs; i++) {
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

    MPI_Bcast(node_ids, nprocs, MPI_INT, root, ncp->comm);

    /* my_node_id is this rank's node ID */
    my_node_id = node_ids[rank];

    /* nprocs_my_node: the number of processes in my nodes
     * ranks_my_node[]: rank IDs of all processes in my node.
     * my_rank_index points to ranks_my_node[] where
     * ranks_my_node[my_rank_index] == rank
     */
    ranks_my_node = (int*) NCI_Malloc(sizeof(int) * nprocs);
    my_rank_index = -1;
    nprocs_my_node = 0;
    for (i=0; i<nprocs; i++) {
        if (node_ids[i] == my_node_id) {
            if (i == rank)
                my_rank_index = nprocs_my_node;
            ranks_my_node[nprocs_my_node] = i;
            nprocs_my_node++;
        }
    }
    assert(my_rank_index >= 0);

    /* Now, ranks_my_node[my_rank_index] == rank */

#ifdef WKL_DEBUG
if (nprocs ==2) printf("%d: node_ids=%d %d ranks_my_node=%d %d\n",rank,node_ids[0],node_ids[1],ranks_my_node[0],ranks_my_node[1]);
else if (nprocs ==5) printf("%d: node_ids=%d %d %d %d %d ranks_my_node=%d %d %d %d %d\n",rank,
node_ids[0],node_ids[1],node_ids[2],node_ids[3],node_ids[4],
ranks_my_node[0],ranks_my_node[1],ranks_my_node[2],ranks_my_node[3],ranks_my_node[4]);
printf("%d: my_node_id=%d nprocs_my_node=%d my_rank_index=%d\n",rank,my_node_id,nprocs_my_node,my_rank_index);
#endif

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
         */
        ncp->aggregation = 0;
    else {
        /* find the rank ID of aggregator assigned to this rank */
        ncp->my_aggr = ranks_my_node[my_rank_index - my_rank_index % num_nonaggrs];

        if (ncp->my_aggr == rank) { /* this rank is an aggregator */
            /* Set the number of non-aggregators assigned to this rank. For the
             * last group, make sure it does not go beyond nprocs_my_node.
             */
            ncp->num_nonaggrs = MIN(num_nonaggrs, nprocs_my_node - my_rank_index);
            if (ncp->num_nonaggrs == 1)
                /* this aggregation group contains only self rank */
                ncp->aggregation = 0;
            else
                /* copy the rank IDs over to ncp->nonaggr_ranks[] */
                memcpy(ncp->nonaggr_ranks,
                       ranks_my_node + my_rank_index,
                       sizeof(int) * num_nonaggrs);
        }
    }
    NCI_Free(ranks_my_node);

    if (ncp->aggregation == 0) {
        /* free ncp->nonaggr_ranks if aggregation is not enabled */
        NCI_Free(ncp->nonaggr_ranks);
        ncp->nonaggr_ranks = NULL;
    }

#ifdef WKL_DEBUG
if (ncp->aggregation && ncp->rank == ncp->my_aggr) printf("%d %s:%d ncp aggregation=%d isAggr=%d my_aggr=%d num_nonaggrs=%d ranks=%d .. %d\n",
ncp->rank,__func__,__LINE__,ncp->aggregation,(ncp->my_aggr==ncp->rank),ncp->my_aggr,ncp->num_nonaggrs,ncp->nonaggr_ranks[0],ncp->nonaggr_ranks[ncp->num_nonaggrs-1]);
#endif


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


/* C struct for breaking down a request to a list of offset-length segments */
typedef struct {
    MPI_Offset off;      /* starting file offset of the request */
    MPI_Offset len;      /* requested length in bytes starting from off */
    MPI_Aint   buf_addr; /* distance of this request's I/O buffer to the first
                            request to be merged */
} off_len;

/*----< off_compare() >------------------------------------------------------*/
/* used for sorting the offsets of the off_len array */
static int
off_compare(const void *a, const void *b)
{
    if (((off_len*)a)->off > ((off_len*)b)->off) return  1;
    if (((off_len*)a)->off < ((off_len*)b)->off) return -1;
    return 0;
}

/*----< flatten_subarray() >-------------------------------------------------*/
/* flatten a subarray request into a list of offset-length pairs */
static int
flatten_subarray(int          ndim,    /* number of dimensions */
                 int          el_size, /* array element size */
                 MPI_Offset   offset,  /* starting file offset of variable */
                 MPI_Offset  *dimlen,  /* [ndim] dimension lengths */
                 MPI_Offset  *start,   /* [ndim] starts of subarray */
                 MPI_Offset  *count,   /* [ndim] counts of subarray */
                 MPI_Offset  *stride,  /* [ndim] strides of subarray */
                 MPI_Offset  *npairs,  /* OUT: number of off-len pairs */
                 MPI_Offset  *offsets, /* OUT: array of offsets */
                 MPI_Offset  *lengths) /* OUT: array of lengths */
{
    int i, j, to_free_stride=0;
    MPI_Offset length, nstride, array_len, off, subarray_len;
    size_t idx=0, idx0;

    *npairs = 0;
    if (ndim < 0) return NC_NOERR;

    if (ndim == 0) {  /* scalar record variable */
        *npairs = 1;
        offsets[0] = offset;
        lengths[0] = el_size;
        return NC_NOERR;
    }

    if (stride == NULL) { /* equivalent to {1, 1, ..., 1} */
        stride = (MPI_Offset*) NCI_Malloc((size_t)ndim * SIZEOF_MPI_OFFSET);
        for (i=0; i<ndim; i++) stride[i] = 1;
        to_free_stride = 1;
    }

    /* TODO: check if all stride[] >= 1
       Q: Is it legal if any stride[] <= 0 ? */

    /* calculate the number of offset-length pairs */
    *npairs = (stride[ndim-1] == 1) ? 1 : count[ndim-1];

#ifdef WKL_DEBUG
for (i=0; i<ndim-1; i++) if (count[i] == 0) printf("++++++++++++++++++++++++ ERROR %d count[%d]=0\n",__LINE__,i);
#endif

    for (i=0; i<ndim-1; i++)
        *npairs *= count[i];
    if (*npairs == 0) {  /* not reachable, an error if count[] == 0 */
        if (to_free_stride) NCI_Free(stride);
        return NC_NOERR;
    }

    /* length of each row of the subarray are of the same size */
    length  = (stride[ndim-1] == 1) ? count[ndim-1] : 1;
    length *= el_size;
    nstride  = (stride[ndim-1] == 1) ? 1 : count[ndim-1];

    /* set the offset-length pairs for the lowest dimension */
    off = offset + start[ndim-1] * el_size;
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
    if (to_free_stride) NCI_Free(stride);

    return NC_NOERR;
}

/*----< flatten_reqs() >-----------------------------------------------------*/
/* flatten all write requests into offset-length pairs.
 * Variable pairs is allocated here and need to be freed by the caller
 */
static int
flatten_reqs(NC            *ncp,
             int            num_reqs,  /* IN: # requests */
             const NC_req  *reqs,      /* [num_reqs] requests */
             MPI_Aint      *num_pairs, /* OUT: total number of off-len pairs */
             MPI_Offset   **offsets,   /* OUT: array of flattened offsets */
             MPI_Offset   **lengths)   /* OUT: array of flattened lengths */
{
    int i, j, status=NC_NOERR, ndims;
    MPI_Offset idx, num, *start, *count, *shape, *stride;

    *num_pairs = 0;    /* total number of offset-length pairs */

    /* Count the number off-len pairs from reqs[], so we can malloc a
     * contiguous memory space for storing off-len pairs
     */
    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;
        ndims = lead->varp->ndims;
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
    *offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (*num_pairs) * 2);
    *lengths = *offsets + (*num_pairs);
    idx = 0;

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
                         start, count, stride,
                         &num,            /* OUT: number of off-len pairs */
                         *offsets + idx,  /* OUT: array of offsets */
                         *lengths + idx); /* OUT: array of lengths */
        idx += num;
    }

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
    else
        MPI_Type_commit(bufType);

    NCI_Free(blocklens);
    NCI_Free(disps);

    return status;
}

/*----< ncmpio_intra_node_aggregation() >------------------------------------*/
/* This is a collective call */
int
ncmpio_intra_node_aggregation(NC     *ncp,
                              int     num_reqs,
                              NC_req *put_list,
                              MPI_Offset newnumrecs)
{
    int i, j, err, mpireturn, status=NC_NOERR, nreqs;
    char *recv_buf=NULL, *wr_buf = NULL;
    MPI_Aint npairs, *msg, bufLen, num_pairs;
    MPI_Offset offset=0, buf_count, *offsets=NULL, *lengths=NULL;
    MPI_Datatype fileType=MPI_BYTE, bufType=MPI_BYTE;
    MPI_File fh;
    MPI_Request *req;

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     */
    if (num_reqs > 0)
        flatten_reqs(ncp, num_reqs, put_list, &num_pairs, &offsets, &lengths);
    else
        num_pairs = 0;

// #define WKL_DEBUG
#ifdef WKL_DEBUG
for (j=1; j<num_pairs; j++) {
    if (offsets[j-1] + lengths[j-1] > offsets[j]) {
        printf("%d: DECREASING num_pairs=%ld offsets[%d]=%lld len=%lld = %lld > offsets[%d]=%lld\n",ncp->rank,num_pairs,j-1,offsets[j-1], lengths[j-1], offsets[j-1] + lengths[j-1], j, offsets[j]);
        break;
    }
}
printf("%d: %s:%d num_reqs=%d num_pairs=%ld offsets=%lld %lld %lld %lld lengths=%lld %lld %lld %lld\n",ncp->rank,__func__,__LINE__,num_reqs,num_pairs,
offsets[0], offsets[1], offsets[2], offsets[3],
lengths[0], lengths[1], lengths[2], lengths[3]);
#endif

    /* construct write buffer datatype, bufType.
     * bufLen is the buffer size in bytes
     */
    if (num_reqs > 0)
        construct_buf_type(ncp, num_reqs, put_list, &bufLen, &bufType);
    else
        bufLen = 0;

    if (put_list != NULL)
        NCI_Free(put_list);

    /* First, tell aggregator how much to receive by sending:
     * (num_pairs and bufLen). The message size to be sent by this rank
     * is num_pairs * 2 * sizeof(MPI_Offset) + bufLen
     */
    if (ncp->rank == ncp->my_aggr) {
        nreqs = ncp->num_nonaggrs - 1;
        msg = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs * 2);
    }
    else {
        nreqs = 1;
        msg = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * 2);
    }

    msg[0] = num_pairs;
    msg[1] = bufLen;
#ifdef WKL_DEBUG
printf("%s:%d rank %d num_pairs=%ld bufLen=%ld\n",__func__,__LINE__,ncp->rank,num_pairs,bufLen);
#endif

    req = (MPI_Request*) NCI_Malloc(sizeof(MPI_Request) * nreqs);
    nreqs = 0;
    if (ncp->rank == ncp->my_aggr) {
        for (i=1; i<ncp->num_nonaggrs; i++)
            MPI_Irecv(msg + i*2, 2, MPI_AINT, ncp->nonaggr_ranks[i], 0,
                      ncp->comm, &req[nreqs++]);
    }
    else
        MPI_Isend(msg, 2, MPI_AINT, ncp->my_aggr, 0, ncp->comm, &req[nreqs++]);

    MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);

    nreqs = 0;
    if (ncp->rank == ncp->my_aggr) {
#ifdef WKL_DEBUG
for (i=0; i<ncp->num_nonaggrs; i++) printf("%s:%d from %d num_pairs=%ld bufLen=%ld\n",__func__,__LINE__,i,msg[2*i],msg[2*i+1]);
#endif
        /* calculate the total number of offset-length pairs */
        npairs = num_pairs;
        for (i=1; i<ncp->num_nonaggrs; i++) npairs += msg[i*2];

        if (npairs > 0) {
            /* realloc to store all pairs in a contiguous buffer */
            offsets = (MPI_Offset*) NCI_Realloc(offsets, sizeof(MPI_Offset)
                                                * npairs * 2);
            lengths = offsets + num_pairs;
        }

        /* post requests to receive offset-length pairs from non-aggregators */
        MPI_Offset *ptr = offsets + num_pairs * 2;
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (msg[i*2] == 0) continue;
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Irecv_c(ptr, msg[i*2] * 2, MPI_OFFSET, ncp->nonaggr_ranks[i],
                        0, ncp->comm, &req[nreqs++]);
#else
            MPI_Irecv(ptr, msg[i*2] * 2, MPI_OFFSET, ncp->nonaggr_ranks[i],
                      0, ncp->comm, &req[nreqs++]);
#endif
            ptr += msg[i*2] * 2;
        }
#ifdef WKL_DEBUG
printf("%s:%d offsets=%lld %lld %lld %lld lengths=%lld %lld %lld %lld\n",__func__,__LINE__,
offsets[0], offsets[1], offsets[2], offsets[3],
lengths[0], lengths[1], lengths[2], lengths[3]);
#endif
    }
    else if (num_pairs > 0) {
        /* send offset-length pairs data to the aggregator */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Isend_c(offsets, num_pairs * 2, MPI_OFFSET, ncp->my_aggr, 0,
                    ncp->comm, &req[nreqs++]);
#else
        MPI_Isend(offsets, num_pairs * 2, MPI_OFFSET, ncp->my_aggr, 0,
                  ncp->comm, &req[nreqs++]);
#endif
    }

    MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);

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
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count position=0;
            MPI_Pack_c(MPI_BOTTOM, 1, bufType, recv_buf, bufLen, &position,
                       MPI_COMM_SELF);
#else
            int position=0;
            MPI_Pack(MPI_BOTTOM, 1, bufType, recv_buf, bufLen, &position,
                     MPI_COMM_SELF);
#endif
        }

        /* post requests to receive write data from non-aggregators */
        char *ptr = recv_buf + bufLen;
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (msg[i*2 + 1] == 0) continue;
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Irecv_c(ptr, msg[i*2 + 1], MPI_BYTE, ncp->nonaggr_ranks[i], 0,
                        ncp->comm, &req[nreqs++]);
#else
            MPI_Irecv(ptr, msg[i*2 + 1], MPI_BYTE, ncp->nonaggr_ranks[i], 0,
                      ncp->comm, &req[nreqs++]);
#endif
            ptr += msg[i*2 + 1];
        }
    }
    else if (bufLen > 0) {
        /* send write data to the aggregator */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Isend_c(MPI_BOTTOM, 1, bufType, ncp->my_aggr, 0, ncp->comm,
                    &req[nreqs++]);
#else
        MPI_Isend(MPI_BOTTOM, 1, bufType, ncp->my_aggr, 0, ncp->comm,
                  &req[nreqs++]);
#endif
        NCI_Free(offsets);
    }

    MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
    NCI_Free(req);

    /* free and reset bufType */
    if (bufType != MPI_BYTE && bufType != MPI_DATATYPE_NULL)
        MPI_Type_free(&bufType);
    bufType = MPI_BYTE;

    /* aggregator sorts the offset-length pairs, along with the buffer */
    if (ncp->rank == ncp->my_aggr && npairs > 0) {
        off_len *segs = (off_len*)NCI_Malloc(sizeof(off_len) * npairs);
        MPI_Offset *off_ptr = offsets;
        int k=0;
        for (i=0; i<ncp->num_nonaggrs; i++) {
            MPI_Offset *len_ptr = off_ptr + msg[i*2];
            for (j=0; j<msg[i*2]; j++) {
                segs[k].off = off_ptr[j];
                segs[k].len = len_ptr[j];
                if (k == 0)
                    segs[k].buf_addr = 0;
                else
                    segs[k].buf_addr = segs[k-1].buf_addr + segs[k-1].len;
                k++;
            }
            off_ptr += msg[i*2] * 2;
        }
        qsort(segs, npairs, sizeof(off_len), off_compare);

#ifdef WKL_DEBUG
for (j=1; j<npairs; j++) {
    if (segs[j-1].off >= segs[j].off) {
        printf("zzzzzzzzzzz npairs=%ld segs[%d].off=%lld >= segs[%d].off=%lld\n",npairs,j-1,segs[j-1].off, j, segs[j].off);
        break;
    }
}
for (j=1; j<npairs; j++) {
    if (segs[j-1].off + segs[j-1].len > segs[j].off) {
        printf("ooooooooverlap npairs=%ld segs[%d].off=%lld len=%lld =%lld> segs[%d].off=%lld len=%lld\n",npairs,j-1,segs[j-1].off, segs[j-1].len, segs[j-1].off + segs[j-1].len,j, segs[j].off,segs[j].len);
        break;
    }
}

printf("%d: %s:%d npairs=%ld offsets=%lld %lld %lld %lld lengths=%lld %lld %lld %lld buf_addr=%ld %ld %ld %ld\n",ncp->rank,__func__,__LINE__,npairs,
segs[0].off, segs[1].off, segs[2].off, segs[3].off,
segs[0].len, segs[1].len, segs[2].len, segs[3].len,
segs[0].buf_addr, segs[1].buf_addr, segs[2].buf_addr, segs[3].buf_addr);
printf("%d: %s:%d npairs=%ld offsets[40]=%lld %lld %lld %lld lengths=%lld %lld %lld %lld buf_addr=%ld %ld %ld %ld\n",ncp->rank,__func__,__LINE__,npairs,
segs[40].off, segs[41].off, segs[42].off, segs[43].off,
segs[40].len, segs[41].len, segs[42].len, segs[43].len,
segs[40].buf_addr, segs[41].buf_addr, segs[42].buf_addr, segs[43].buf_addr);
#endif

        NCI_Free(offsets);

        /* merge the overlapped buffer segments, skip the overlapped regions
         * for those with higher j indices (i.e. requests with lower j indices
         * win the writes to the overlapped regions)
         */
#ifdef WKL_DEBUG
MPI_Aint wkl = segs[0].len;
for (i=0, j=1; j<npairs; j++) { wkl += segs[j].len; }
#endif
        for (i=0, j=1; j<npairs; j++) {
            if (segs[i].off + segs[i].len >= segs[j].off + segs[j].len)
                /* segment i completely covers segment j, skip j */
                continue;

            MPI_Offset gap = segs[i].off + segs[i].len - segs[j].off;
            if (gap >= 0) { /* segments i and j overlaps */
                if (MPI_Aint_add(segs[i].buf_addr, segs[i].len) ==
                    MPI_Aint_add(segs[j].buf_addr, gap)) {
                    /* buffers i and j are contiguous, merge j to i */
                    segs[i].len = MPI_Aint_add(segs[i].len, segs[j].len - gap);
                }
                else { /* buffers are not contiguous, reduce j's len */
                    segs[i+1].off      = segs[j].off + gap;
                    segs[i+1].len      = segs[j].len - gap;
                    segs[i+1].buf_addr = MPI_Aint_add(segs[j].buf_addr, gap);
                    i++;
                }
            }
            else { /* i and j do not overlap */
                i++;
                if (i < j) segs[i] = segs[j];
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
        for (i=0; i<npairs; i++) {
            memcpy(ptr, recv_buf + segs[i].buf_addr, segs[i].len);
            ptr += segs[i].len;
            /* overlap may be found, recalculate buf_count */
            buf_count += segs[i].len;
        }
        if (recv_buf != NULL) NCI_Free(recv_buf);
        bufType = MPI_BYTE;

        /* coalesce the offset-length pairs */
        for (i=0, j=1; j<npairs; j++) {
            if (segs[i].off + segs[i].len == segs[j].off) {
                /* coalesce j into i */
                segs[i].len += segs[j].len;
            }
            else {
                i++;
                if (i < j) segs[i] = segs[j];
            }
        }

#ifdef WKL_DEBUG
printf("%s: %d npairs=%ld i+1=%d\n",__func__,__LINE__,npairs, i+1);
#endif
        /* update number of segments, now all off-len pairs are not overlapped */
        npairs = i+1;

        if (npairs == 1) {
            offset = segs[0].off;
#ifdef WKL_DEBUG
printf("%s: %d npairs=%ld i+1=%d offset=%lld\n",__func__,__LINE__,npairs, i+1,offset);
#endif
        }
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            /* construct fileview */
            MPI_Count *blocklens = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * npairs);
            MPI_Count *disps     = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * npairs);
            for (i=0; i<npairs; i++) {
                disps[i]     = segs[i].off;
                blocklens[i] = segs[i].len;
            }

            mpireturn = MPI_Type_create_hindexed_c(npairs, blocklens, disps, MPI_BYTE,
                                                   &fileType);

            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
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

            NCI_Free(blocklens);
            NCI_Free(disps);
#else
            /* construct fileview */
            int      *blocklens = (int*)      NCI_Malloc(SIZEOF_INT * npairs);
            MPI_Aint *disps     = (MPI_Aint*) NCI_Malloc(SIZEOF_MPI_AINT * npairs);
            for (i=0; i<npairs; i++) {
                if (segs[i].len > NC_MAX_INT) {
                    NCI_Free(disps);
                    NCI_Free(blocklens);
                    DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
                }
                disps[i]     = segs[i].off;
                blocklens[i] = segs[i].len;
            }

            mpireturn = MPI_Type_create_hindexed(npairs, blocklens, disps, MPI_BYTE,
                                                 &fileType);

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

            NCI_Free(blocklens);
            NCI_Free(disps);
#endif
        }
        NCI_Free(segs);
    }
    NCI_Free(msg);

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
                            bufType, wr_buf, ((bufType == MPI_BYTE) ? 1 : 0));
    if (bufType != MPI_BYTE) MPI_Type_free(&bufType);
    if (status == NC_NOERR) status = err;

    if (wr_buf != NULL) NCI_Free(wr_buf);

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

