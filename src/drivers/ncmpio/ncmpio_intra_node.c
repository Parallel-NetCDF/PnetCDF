/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * This file contains the implementation of intra-node aggregation feature,
 * which is designed to improve performance for I/O patterns that contain many
 * noncontiguous requests interleaved among processes, with a wide aggregate
 * access region on each process which results in an all-to-many personalized
 * communication with each MPI I/O aggregator receiving data from all MPI
 * processes. By reducing the number of processes per node participating the
 * all-to-many communication (i.e. the 'all' part), this feature effectively
 * reduces the communication contention, which often happens to jobs that run
 * on a large the number of MPI processes per compute node.
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
#include <limits.h>   /* INT_MAX */
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
    MPI_Aint nelems;

    for (nelems=0, i=0; i<nprocs; i++) nelems += count[i];

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

    /* offsets[] and lengths[] may be coalescable, but this will be delayed
     * until this subroutine callers, flatten_req() and flatten_nreqs(), where
     * both fix-sized and record variables have been flattened.
     */

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
            PNCIO_View        *file_view) /* OUT: flattened file view */
{
    int i, j, err=NC_NOERR, ndims;
    MPI_Aint num, idx;
    MPI_Offset var_begin, *shape, count0, *ones=NULL;
    MPI_Offset prev_end_off;

    file_view->count = 0;   /* number of offset-length pairs */
    file_view->off = NULL;  /* array of offsets */
    file_view->len = NULL;  /* array of lengths */

    /* Count the number offset-length pairs */

    if (varp->ndims == 0) { /* scalar variable has 1 pair offset-length */
        file_view->count  = 1;
        file_view->off    = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset));
        file_view->len    = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset));
        file_view->off[0] = varp->begin;
        file_view->len[0] = varp->xsz;
        return NC_NOERR;
    }

    if (varp->ndims == 1 && IS_RECVAR(varp)) {
        /* scalar record variable has 1 pair offset-length */
        file_view->count = count[0];
    }
    else {
        file_view->count = 1;
        if (stride != NULL && stride[varp->ndims-1] > 1)
            /* count of last dimension */
            file_view->count = count[varp->ndims-1];

        for (i=0; i<varp->ndims-1; i++)
            /* all count[] except the last dimension */
            file_view->count *= count[i];
    }

    file_view->off = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * file_view->count);
    file_view->len = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * file_view->count);

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
                               file_view->off + idx,  /* OUT: array of offsets */
                               file_view->len + idx); /* OUT: array of lengths */
        if (err != NC_NOERR) goto err_out;
        if (num == 0) continue;

        /* check if file_view->off[] are in an increasing order */
        for (j=0; j<num; j++) {
            if (prev_end_off > file_view->off[idx+j])
                *is_incr = 0;  /* file_view->off are not incrementing */
            else
                prev_end_off = file_view->off[idx+j];
        }

        idx += num;
        assert(idx <= file_view->count);

        if (IS_RECVAR(varp))
            var_begin += ncp->recsize;
    }

    file_view->count = idx;

    /* check if the offsets-lengths can be coalesced */
    for (i=0, j=1; j<file_view->count; j++) {
        if (file_view->off[i] + file_view->len[i] == file_view->off[j])
            file_view->len[i] += file_view->len[j];
        else {
            i++;
            if (i < j) {
                file_view->off[i] = file_view->off[j];
                file_view->len[i] = file_view->len[j];
            }
        }
    }
    file_view->count = i + 1;

err_out:
    if (ones != NULL) NCI_Free(ones);
    if (err != NC_NOERR) {
        NCI_Free(file_view->off);
        NCI_Free(file_view->len);
        file_view->off = NULL;
        file_view->len = NULL;
        file_view->count = 0;
    }
    return err;
}

/*----< flatten_nreqs() >----------------------------------------------------*/
/* Flatten multiple subarray requests into file offset-length pairs. Arrays
 * offsets and lengths are allocated here and need to be freed by the caller.
 */
static int
flatten_nreqs(NC            *ncp,
              int            reqMode,   /* IN: NC_REQ_RD or NC_REQ_WR */
              int            num_reqs,  /* IN: # requests */
              const NC_req  *reqs,      /* [num_reqs] requests */
              int           *is_incr,   /* OUT: are offsets incrementing */
              PNCIO_View    *file_view) /* OUT: flattened file view */
{
    int i, j, status=NC_NOERR, ndims, max_ndims=0;
    MPI_Aint num, idx;
    MPI_Offset *start, *count, *shape, *stride, *ones;
    MPI_Offset prev_end_off;

    file_view->count = 0; /* total number of offset-length pairs */
    file_view->off = NULL;
    file_view->len = NULL;
    file_view->size = 0;

    if (num_reqs == 0) return NC_NOERR;

    /* Count the number off-len pairs from reqs[], so we can malloc a
     * contiguous memory space for storing off-len pairs
     */
    for (i=0; i<num_reqs; i++) {
        /* reqs[i].npairs is the number of offset-length pairs of this request,
         * calculated in ncmpio_igetput_varm() and igetput_varn()
         */
        file_view->count += reqs[i].npairs;
        if (fIsSet(reqMode, NC_REQ_WR))
            ndims = ncp->put_lead_list[reqs[i].lead_off].varp->ndims;
        else
            ndims = ncp->get_lead_list[reqs[i].lead_off].varp->ndims;
        max_ndims = MAX(max_ndims, ndims);
    }

    /* now we can allocate a contiguous memory space for the off-len pairs */
    file_view->off = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * file_view->count);
    file_view->len = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * file_view->count);

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
            file_view->off[idx] = reqs[i].offset_start;
            file_view->len[idx] = reqs[i].nelems * lead->varp->xsz;

            /* check if file_view->off[] are in an increasing order */
            if (prev_end_off > file_view->off[idx])
                *is_incr = 0;  /* file_view->off are not incrementing */
            else
                prev_end_off = file_view->off[idx];
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
                         file_view->off + idx,  /* OUT: array of offsets */
                         file_view->len + idx); /* OUT: array of lengths */

        /* check if file_view->off[] are in an increasing order */
        for (j=0; j<num; j++) {
            if (prev_end_off > file_view->off[idx+j])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_end_off = file_view->off[idx+j];
        }
        idx += num;
    }
    NCI_Free(ones);

    file_view->count = idx;

    /* check if the offsets-lengths can be coalesced */
    for (i=0, j=1; j<file_view->count; j++) {
        if (file_view->off[i] + file_view->len[i] == file_view->off[j])
            file_view->len[i] += file_view->len[j];
        else {
            i++;
            if (i < j) {
                file_view->off[i] = file_view->off[j];
                file_view->len[i] = file_view->len[j];
            }
        }
    }
    file_view->count = i + 1;

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
              int            num_reqs, /* IN: number of requests */
              const NC_req  *reqs,     /* IN: [num_reqs] requests */
              PNCIO_View    *buf_view, /* OUT: flattened buftype */
              void         **buf)      /* OUT: pointer to I/O buffer */
{
    int i, j, err=NC_NOERR;
    NC_lead_req *lead;
    MPI_Aint addr, addr0;
    /* TODO: buffer offset should be of type MPI_Aint. length should be size_t. */

    buf_view->size = 0;
    buf_view->count = 0;
    buf_view->off = NULL;
    buf_view->len = NULL;

    if (num_reqs == 0)
        return NC_NOERR;

    buf_view->off = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * num_reqs);
    if (buf_view->off == NULL) return NC_ENOMEM;
    buf_view->len = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * num_reqs);
    if (buf_view->len == NULL) return NC_ENOMEM;

    *buf = reqs[0].xbuf;

    lead = (fIsSet(reqMode, NC_REQ_WR)) ? ncp->put_lead_list
                                        : ncp->get_lead_list;

    MPI_Get_address(lead[reqs[0].lead_off].xbuf, &addr0);

    /* set buf_view->off[0] and buf_view->len[0] */
    MPI_Get_address(reqs[0].xbuf, &addr0); /* displacement uses MPI_BOTTOM */
    buf_view->off[0] = 0;

    /* buf_view->len[] are in bytes */
    buf_view->len[0] = reqs[0].nelems * lead[reqs[0].lead_off].varp->xsz;

    buf_view->size = buf_view->len[0];
    for (i=0, j=1; j<num_reqs; j++) {
        /* displacement uses MPI_BOTTOM */
        MPI_Get_address(reqs[j].xbuf, &addr);
        buf_view->off[j] = addr - addr0;

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
    /* After coalescing, the true number of requests may be reduced, but should
     * still be > 0.
     */

    if (i + 1 < num_reqs) {
        num_reqs = i + 1; /* num_reqs is reduced */
        buf_view->off = (MPI_Offset*)NCI_Realloc(buf_view->off,
                                       sizeof(MPI_Offset) * num_reqs);
        if (buf_view->off == NULL) return NC_ENOMEM;
        buf_view->len = (MPI_Offset*)NCI_Realloc(buf_view->len,
                                       sizeof(MPI_Offset) * num_reqs);
        if (buf_view->len == NULL) return NC_ENOMEM;
    }
    assert(num_reqs > 0);

    buf_view->count = num_reqs;

    return err;
}

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
    int i, err, rank, nprocs, mpireturn, status=NC_NOERR, nreqs;
    MPI_Request *req=NULL;
    MPI_Aint num_pairs=meta[0];
    MPI_Comm comm = ncp->comm_attr.ina_intra_comm;

#if PNETCDF_DEBUG_MODE == 1
    assert(comm != MPI_COMM_NULL);
#endif

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    /* Aggregator collects each non-aggregator's num_pairs and bufLen */
    if (rank == 0) {

        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nprocs);
        nreqs = 0;
        for (i=1; i<nprocs; i++)
            TRACE_COMM(MPI_Irecv)(meta + i*3, 3, MPI_AINT, i, 0, comm,
                                  &req[nreqs++]);

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
        TRACE_COMM(MPI_Send)(meta, 3, MPI_AINT, 0, 0, comm);

    /* Secondly, aggregators collect offset-length pairs from all its
     * non-aggregators
     */
    if (rank == 0) {
        MPI_Datatype recvType;

        /* calculate the total number of offset-length pairs to receive */
        for (*npairs=0, i=0; i<nprocs; i++) *npairs += meta[i*3];

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
            *lengths = (int*) NCI_Realloc(*lengths, *npairs * sizeof(int));
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
        for (i=1; i<nprocs; i++) {
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
            TRACE_COMM(MPI_Irecv_c)(MPI_BOTTOM, 1, recvType, i, 0, comm,
                                    &req[nreqs++]);
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
        for (i=1; i<nprocs; i++) {
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
            TRACE_COMM(MPI_Irecv)(MPI_BOTTOM, 1, recvType, i, 0, comm,
                                  &req[nreqs++]);
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
        TRACE_COMM(MPI_Send_c)(MPI_BOTTOM, 1, sendType, 0, 0, comm);
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
        TRACE_COMM(MPI_Send)(MPI_BOTTOM, 1, sendType, 0, 0, comm);
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
    char *recv_buf=NULL, *wr_buf = NULL;
    int i, j, err, mpireturn, status=NC_NOERR, rank, nprocs;
    int coalesceable=0, free_buf_view_off=0;;
    MPI_Aint npairs=0, *meta=NULL, *bufAddr=NULL;
    MPI_Aint buf_npairs=0, file_npairs=0;
    MPI_Offset wr_amnt=0;
    MPI_Comm intra_comm;
    PNCIO_View wr_buf_view=buf_view;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *off_ptr, *len_ptr, *file_len=NULL;
#else
    MPI_Offset *off_ptr;
    int *len_ptr, *file_len=NULL;
#endif

MPI_Count buf_off=0, buf_len;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double endT, startT = MPI_Wtime();
    MPI_Offset mem_max;
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[0] = MAX(pnc_ina_mem_put[0], mem_max);
#endif

    intra_comm = ncp->comm_attr.ina_intra_comm;
    if (intra_comm == MPI_COMM_NULL) { /* INA is disabled */
        nprocs = 1;
        rank = 0;
    }
    else {
        MPI_Comm_size(intra_comm, &nprocs);
        MPI_Comm_rank(intra_comm, &rank);
    }

    /* Firstly, aggregators collect metadata from non-aggregators.
     *
     * This rank tells its aggregator how much metadata to receive from this
     * rank, by sending: the number of offset-length pairs (num_pairs) and user
     * buffer size in bytes (buf_view.size). This message size to be sent by
     * this rank is 3 MPI_Offset.
     */
    if (rank == 0)
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nprocs * 3);
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
     * MPI-IO or GIO file write.
     *
     * Once ina_collect_md() returns, this aggregator's offsets and lengths may
     * grow to include the ones from non-aggregators (appended).
     */
    if (nprocs > 1) {
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
    if (rank > 0) { /* non-aggregator */
        if (buf_view.count > 0) {
            MPI_Datatype sendType=MPI_BYTE;
            /* Non-aggregators send write data to the aggregator */
            if (buf_view.count > 1) {
                err = ncmpio_type_create_hindexed(buf_view.count, buf_view.off,
                                                  buf_view.len, &sendType);
                if (status == NC_NOERR) status = err;
            }
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count num = (buf_view.count == 1) ? buf_view.size : 1;
            TRACE_COMM(MPI_Send_c)(buf, num, sendType, 0, 0, intra_comm);
#else
            int num = (buf_view.count == 1) ? buf_view.size : 1;
            TRACE_COMM(MPI_Send)(buf, num, sendType, 0, 0, intra_comm);
#endif
            if (buf_view.count > 1) MPI_Type_free(&sendType);

            /* free space allocated for buf_view */
            NCI_Free(buf_view.off);
            NCI_Free(buf_view.len);
        }

        /* free space allocated for file_view */
        if (num_pairs > 0) {
            NCI_Free(offsets);
            NCI_Free(lengths);
        }

        /* Non-aggregators are done here, as only aggregators call MPI-IO/GIO
         * functions to write data to the file. Non-aggregators do not
         * participate MPI-IO calls.
         */
        NCI_Free(meta);
        return status;
    }

    /* The remaining of this subroutine is for aggregators only */

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[1] = MAX(pnc_ina_mem_put[1], mem_max);

    endT = MPI_Wtime();
    pnc_ina_put[0] += endT - startT; /* collect MD */
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
        for (i=-1,j=0; j<nprocs; j++) {
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
            for (++i; i<nprocs; i++) {
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

        /* Construct an array of buffer addresses containing a mapping of the
         * buffer used to receive write data from non-aggregators and the
         * buffer used to write to file. bufAddr[] is calculated based on the
         * assumption that the write buffer of this aggregator is contiguous,
         * i.e. buf_view.count <= 1. For non-aggregators, their write data will
         * always be received into a contiguous buffer.
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
                /* Interleaved offsets are found but individual offsets are
                 * already sorted. This is commonly seen from the checkerboard
                 * domain partitioning pattern. In this case, heap_merge() is
                 * faster to merge all offsets into one single sorted offset
                 * list. Note count[] must be initialized, so it can be used
                 * in heap_merge()
                 */
                MPI_Aint *count;
                count = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nprocs);
                for (i=0; i<nprocs; i++)
                    count[i] = meta[i*3];

                /* heap-merge() runs much faster than qsort() when individual
                 * lists have already been sorted. However, it has a much
                 * bigger memory footprint.
                 */
                heap_merge(nprocs, count, off_ptr, len_ptr, bufAddr);
                NCI_Free(count);
            }
            else
                /* When some individual offsets are not sorted, we cannot use
                 * heap_merge(). Note qsort() is an in-place sorting.
                 */
                qsort_off_len_buf(npairs, off_ptr, len_ptr, bufAddr);
        }

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
         *
         * Note off_ptr[] has been sorted into a monotonically non-decreasing
         * order. During the sorting, bufAddr[] are moved around based on their
         * corresponding off_ptr[], and thus bufAddr[] may not be in a
         * monotonically non-decreasing order.
         */
        coalesceable = 0;
        overlap = 0;
        wr_amnt = recv_amnt = len_ptr[0];
        for (i=0, j=1; j<npairs; j++) {
            recv_amnt += len_ptr[j];
            if (off_ptr[i] + len_ptr[i] >= off_ptr[j] + len_ptr[j]) {
                overlap = 1;
                /* segment i completely covers segment j, skip j */
                continue;
            }

            MPI_Offset gap = off_ptr[i] + len_ptr[i] - off_ptr[j];
            if (gap >= 0) { /* overlap detected, merge j into i */
                /* when gap > 0,  pairs i and j overlap
                 * when gap == 0, pairs i and j are contiguous
                 */
                if (gap > 0) overlap = 1;
                wr_amnt += len_ptr[j] - gap; /* subtract overlapped amount */
                if (bufAddr[i] + len_ptr[i] == bufAddr[j] + gap) {
                    /* buffers i and j are contiguous, merge j into i and
                     * subtract overlapped amount.
                     */
                    len_ptr[i] += len_ptr[j] - gap;
                }
                else { /* buffers are not contiguous, reduce j's len */
                    coalesceable = 1;
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

        /* Now off_ptr[], len_ptr[], bufAddr[] are coalesced and no overlap */
        npairs = i+1;

        /* buf_view to be used in a call to ncmpio_read_write() later should
         * use bufAddr[] and len_ptr[], as it offset-length pairs.
         */
        buf_npairs = npairs;

        /* coalesce file_view's off_ptr[] and len_ptr[] independently from
         * buf_view's
         */
        file_npairs = npairs;
        if (coalesceable) { /* file_view can be further coalesced */
            size_t cpy_amnt;
#ifdef HAVE_MPI_LARGE_COUNT
            cpy_amnt = sizeof(MPI_Count) * file_npairs;
            file_len = (MPI_Count*) NCI_Malloc(cpy_amnt);
#else
            cpy_amnt = sizeof(int) * file_npairs;
            file_len = (int*) NCI_Malloc(cpy_amnt);
#endif
            memcpy(file_len, len_ptr, cpy_amnt);

            for (i=0, j=1; j<file_npairs; j++) {
#if PNETCDF_DEBUG_MODE == 1
                /* any overlap should have been removed from the loop above */
                assert(off_ptr[i] + file_len[i] <= off_ptr[j]);
#endif
                if (off_ptr[i] + file_len[i] == off_ptr[j])
                    /* coalesce j into i */
                    file_len[i] += file_len[j];
                else { /* i and j are not coalesceable */
                    i++;
                    if (i < j) {
                        off_ptr[i] = off_ptr[j];
                        file_len[i] = file_len[j];
                    }
                }
            }
            /* update number of offset-length pairs */
            file_npairs = i+1;

        }
        else {
            /* file_view can use the same len_ptr[] as buf_view */
            file_len = len_ptr;
        }

#if PNETCDF_DEBUG_MODE == 1
        /* check if file_view's offset-lengths have been coalesced */
        for (i=1; i<file_npairs; i++) {
            assert(file_len[i-1] > 0);
            assert(off_ptr[i-1] < off_ptr[i]);
            assert(off_ptr[i-1] + file_len[i-1] < off_ptr[i]);
        }
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        ncmpi_inq_malloc_size(&mem_max);
        // ncmpi_inq_malloc_max_size(&mem_max);
        pnc_ina_mem_put[2] = MAX(pnc_ina_mem_put[2], mem_max);
        pnc_ina_npairs_put = MAX(pnc_ina_npairs_put, npairs);

        endT = MPI_Wtime();
        pnc_ina_put[1] += endT - startT; /* sorting */
        startT = endT;
#endif

        /* Allocate receive buffer. Once write data from non-aggregators have
         * received into recv_buf, it is packed into wr_buf. Then, wr_buf is
         * used to call MPI-IO/GIO file write. Note the wr_buf is always
         * contiguous.
         *
         * When nprocs == 1, wr_buf is set to buf which is directly passed to
         * MPI-IO/GIO file write.
         *
         * If file offset-length pairs have not been re-ordered, i.e. sorted
         * and overlaps removed, and this aggregator will not receive any write
         * data from its non-aggregators, then we can use user's buffer, buf,
         * to call MPI-IO/GIO to write to the file, without allocating an
         * additional temporary buffer.
         */
        if (!do_sort && buf_view.size == recv_amnt && !overlap)
            recv_buf = buf;
        else
            recv_buf = (char*) NCI_Malloc(recv_amnt);

        if (recv_buf != buf) {
            /* Copy this aggregator's write data into front of recv_buf */
            char *recv_ptr=recv_buf;
            for (j=0; j<buf_view.count; j++) {
                memcpy(recv_ptr, buf+buf_view.off[j], buf_view.len[j]);
                recv_ptr += buf_view.len[j];
            }
        }

        /* Receive write data sent from non-aggregators. Note we cannot move
         * the posting of MPI_Irecv calls to before sorting and leave
         * MPI_Waitall() to after sorting to overlap communication with the
         * sorting, because the sorting determines the receive buffer size.
         */
        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nprocs);
        ptr = recv_buf + buf_view.size;
        nreqs = 0;
        for (i=1; i<nprocs; i++) {
            if (meta[i*3 + 1] == 0) continue;

#ifdef HAVE_MPI_LARGE_COUNT
            TRACE_COMM(MPI_Irecv_c)(ptr, meta[i*3 + 1], MPI_BYTE, i, 0,
                           intra_comm, &req[nreqs++]);
#else
            TRACE_COMM(MPI_Irecv)(ptr, meta[i*3 + 1], MPI_BYTE, i, 0,
                           intra_comm, &req[nreqs++]);
#endif
            ptr += meta[i*3 + 1];
        }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        ncmpi_inq_malloc_size(&mem_max);
        // ncmpi_inq_malloc_max_size(&mem_max);
        pnc_ina_mem_put[3] = MAX(pnc_ina_mem_put[3], mem_max);

        endT = MPI_Wtime();
        pnc_ina_put[2] += endT - startT; /* post irecv */
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
        NCI_Free(req);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        endT = MPI_Wtime();
        pnc_ina_put[3] += endT - startT; /* wait */
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
         * call to MPI-IO/GIO file write.
         */
        if (!do_sort && wr_amnt == recv_amnt) {
            wr_buf = recv_buf;

            if (wr_buf != buf) {
                /* If write data has been packed in wr_buf, a contiguous buffer,
                 * update buf_view before passing it to the MPI-IO/GIO file
                 * write.
                 */
#if 1
                buf_len = wr_amnt;
                wr_buf_view.len = &buf_len;
                wr_buf_view.off = &buf_off;
                wr_buf_view.count = 1;
#else
                wr_buf_view.size = wr_amnt;
                wr_buf_view.count = 0;
#endif
            }
            else /* user's buffer, buf, can be used to write */
                wr_buf_view = buf_view;
        }
        else if (buf_view.count <= 1 && !overlap) {
            /* Note we can reuse bufAddr[] and len_ptr[] as buf_view.off and
             * buf_view.len only when buf_view is contiguous, because bufAddr[]
             * is constructed based on the assumption that the write buffer is
             * contiguous.
             */
            wr_buf = recv_buf;
            wr_buf_view.size      = wr_amnt;
            wr_buf_view.len       = len_ptr;
            wr_buf_view.count     = buf_npairs;
#if SIZEOF_MPI_AINT == SIZEOF_MPI_OFFSET
            wr_buf_view.off = (MPI_Offset*)bufAddr; /* based on recv_buf */
#else
            wr_buf_view.off = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * buf_npairs);
            for (j=0; j<buf_npairs; j++)
                wr_buf_view.off[j] = (MPI_Offset)bufAddr[j];
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
             * I/O MPI-IO and GIO run faster, this memory copy cost may not
             * be worthy. Besides, the memory footprint high-water mark is
             * doubled.
             */
            wr_buf = NCI_Malloc(wr_amnt);
            ptr = wr_buf;

            for (j=0; j<buf_npairs; j++) {
                memcpy(ptr, recv_buf + bufAddr[j], len_ptr[j]);
                ptr += len_ptr[j];
            }
            /* Write data has been packed in wr_buf, a contiguous buffer,
             * update buf_view before passing it to the MPI-IO/GIO file
             * write.
             */
#if 1
            buf_len = wr_amnt;
            wr_buf_view.len = &buf_len;
            wr_buf_view.off = &buf_off;
            wr_buf_view.count = 1;
#else
            wr_buf_view.size = wr_amnt;
            wr_buf_view.count = 0;
#endif

            if (recv_buf != buf) NCI_Free(recv_buf);
        }
    } /* if (npairs > 0) */

    NCI_Free(meta);

#if 1
PNCIO_View file_view;
file_view.count = file_npairs;
file_view.off = off_ptr;
file_view.len = file_len;
// printf("\n%s at %d: file_view count %lld off %lld len %lld\n",__func__,__LINE__,file_npairs,off_ptr[0],file_len[0]);

#else
    /* set the fileview */
    err = ncmpio_file_set_view(ncp, MPI_BYTE, file_npairs, off_ptr, file_len);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        wr_amnt = 0;
    }
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[4] = MAX(pnc_ina_mem_put[4], mem_max);

    endT = MPI_Wtime();
    pnc_ina_put[4] += endT - startT; /* setview */
    startT = endT;
#endif

    /* carry out write request to file */
    int coll_indep = (fIsSet(ncp->flags, NC_MODE_INDEP)) ? NC_REQ_INDEP : NC_REQ_COLL;

    MPI_Offset wlen = ncmpio_file_write(ncp, coll_indep, wr_buf, file_view, wr_buf_view);
    if (wlen < 0) {
        if (status == NC_NOERR) status = (int)wlen;
        wr_amnt = 0;
    }

    if (free_buf_view_off) NCI_Free(wr_buf_view.off);
    if (wr_buf != buf)  NCI_Free(wr_buf);
    if (bufAddr != NULL) NCI_Free(bufAddr);

    /* free space allocated for file_view and buf_view */
    if (npairs > 0) {
        NCI_Free(offsets);
        NCI_Free(lengths);
    }
    if (coalesceable) NCI_Free(file_len);

    if (buf_view.count > 0) {
        NCI_Free(buf_view.off);
        NCI_Free(buf_view.len);
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[5] = MAX(pnc_ina_mem_put[5], mem_max);
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[5] = MAX(pnc_ina_mem_put[5], mem_max);

    endT = MPI_Wtime();
    pnc_ina_put[5] += endT - startT; /* write */
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
    int i, j, err, mpireturn, status=NC_NOERR, nprocs, rank, nreqs;
    int do_sort=0, indv_sorted=1, overlap=0;
    char *rd_buf = NULL;
    MPI_Aint npairs=0, max_npairs, *meta=NULL;
    MPI_Offset send_amnt=0, rd_amnt=0, off_start;
    MPI_Request *req=NULL;
    MPI_Comm intra_comm;
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
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[0] = MAX(pnc_ina_mem_get[0], mem_max);
#endif

    intra_comm = ncp->comm_attr.ina_intra_comm;
    if (intra_comm == MPI_COMM_NULL) { /* INA is disabled */
        nprocs = 1;
        rank = 0;
    }
    else {
        MPI_Comm_size(intra_comm, &nprocs);
        MPI_Comm_rank(intra_comm, &rank);
    }

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
    if (rank == 0)
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nprocs * 3);
    else
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * 3);

    meta[0] = num_pairs;
    meta[1] = bufLen;
    meta[2] = is_incr;

    /* Each aggregator first collects metadata about its own offset-length
     * pairs, size of read request, and checks whether the offsets are in an
     * incremental order. The aggregator will gather the same metadata from the
     * non-aggregators assigned to it.
     *
     * Once ina_collect_md() returns, this aggregator's offsets[] and lengths[]
     * may grow to include the ones from non-aggregators (appended).
     */
    if (nprocs > 1) {
        err = ina_collect_md(ncp, meta, &offsets, &lengths, &npairs);
        if (err != NC_NOERR) {
            NCI_Free(meta);
            return err;
        }
    }
    else
        npairs = num_pairs;

    if (rank > 0) { /* not an INA aggregator */
        /* For read operation, the non-aggregators now can start receiving
         * their read data from the aggregator.
         */
        if (buf_view.count > 0) {
            MPI_Status st;
            MPI_Datatype recvType=MPI_BYTE;

            /* Non-aggregators send write data to the aggregator */
            if (buf_view.count > 1) {
                err = ncmpio_type_create_hindexed(buf_view.count, buf_view.off,
                                                  buf_view.len, &recvType);
                if (status == NC_NOERR) status = err;
            }
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count num = (buf_view.count <= 1) ? buf_view.size : 1;
            TRACE_COMM(MPI_Recv_c)(buf, num, recvType, 0, 0, intra_comm, &st);
#else
            int num = (buf_view.count <= 1) ? buf_view.size : 1;
            TRACE_COMM(MPI_Recv)(buf, num, recvType, 0, 0, intra_comm, &st);
#endif
            if (buf_view.count > 1) MPI_Type_free(&recvType);

            /* free space allocated for buf_view */
            NCI_Free(buf_view.off);
            NCI_Free(buf_view.len);
        }

        /* free space allocated for file_view */
        if (num_pairs > 0) {
            NCI_Free(offsets);
            NCI_Free(lengths);
        }

        /* Non-INA aggregators are now done, as they do not participate MPI-IO
         * or GIO file read (neither collective nor independent).
         */
        NCI_Free(meta);
        return status;
    }

    /* The remaining of this subroutine is for INA aggregators only. */

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[1] = MAX(pnc_ina_mem_get[1], mem_max);

    endT = MPI_Wtime();
    pnc_ina_get[0] += endT - startT; /* collect MD */
    startT = endT;
#endif

    /* For read operations, the original offsets[] and lengths[] must be kept
     * untouched, because the later sorting and coalescing will mess up the
     * original order of offsets[] and lengths[], which are necessary for
     * constructing a datatype used by the INA aggregator to send data read
     * from the file to its non-INA aggregators.
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
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[2] = MAX(pnc_ina_mem_get[2], mem_max);
#endif

    /* MPI-IO has the following requirements about filetype.
     * 1. The (flattened) displacements (of a filetype) are not required to be
     *    distinct, but they cannot be negative, and they must be monotonically
     *    non-decreasing.
     * 2. If the file is opened for writing, neither the etype nor the filetype
     *    is permitted to contain overlapping regions.
     */
    if (npairs > 0) {
        /* Now this INA aggregator has received all offset-length pairs from
         * its non-aggregators. At first, it checks if a sorting is necessary.
         */

        /* Check if offsets of all non-aggregators are individually sorted */
        indv_sorted = 1;
        for (i=-1,j=0; j<nprocs; j++) {
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
            /* Even when the offsets of all non-INA aggregators and this
             * aggregator offsets are individually sorted, we still need to
             * check if offsets are interleaved. If interleaved, we must sort
             * all offset-length pairs.
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
            for (++i; i<nprocs; i++) {
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

        if (do_sort) {
            /* Sort the offsets into an increasing order. Note during sorting,
             * the length and buffer must also be moved together with their
             * corresponding offset.
             *
             * At first, construct an array of buffer addresses containing a
             * mapping of the buffer used to send read data to the non-INA
             * aggregators and the buffer used to read from the file.
             */
            if (indv_sorted) {
                /* Interleaved offsets are found but individual offsets are
                 * already sorted. This is commonly seen from the checkerboard
                 * domain partitioning pattern. In this case, heap_merge() is
                 * faster to merge all offsets into one single sorted offset
                 * list. Note count[] must be initialized, so it can be used
                 * in heap_merge()
                 */
                MPI_Aint *count;
                count = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nprocs);
                for (i=0; i<nprocs; i++)
                    count[i] = meta[i*3];

                /* heap-merge() runs much faster than qsort() when individual
                 * lists have already been sorted. However, it has a much
                 * bigger memory footprint.
                 */
                heap_merge(nprocs, count, off_ptr, len_ptr, NULL);
                NCI_Free(count);
            }
            else
                /* When some individual offsets are not sorted, we cannot use
                 * heap_merge(). Note qsort() is an in-place sorting.
                 */
                qsort_off_len_buf(npairs, off_ptr, len_ptr, NULL);
        }

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
        ncmpi_inq_malloc_size(&mem_max);
        // ncmpi_inq_malloc_max_size(&mem_max);
        pnc_ina_mem_get[2] = MAX(pnc_ina_mem_get[2], mem_max);
        pnc_ina_npairs_get = MAX(pnc_ina_npairs_get, npairs);

        endT = MPI_Wtime();
        pnc_ina_get[1] += endT - startT; /* sorting */
        startT = endT;
#endif
    } /* if (npairs > 0) */
    /* else case: This INA aggregation group has zero data to read, but this
     * aggregator must participate the collective I/O calls.
     */

#if 1
PNCIO_View file_view;
file_view.count = npairs;
file_view.off = off_ptr;
file_view.len = len_ptr;
#else
    /* set the fileview */
    err = ncmpio_file_set_view(ncp, MPI_BYTE, npairs, off_ptr, len_ptr);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        rd_amnt = 0;
    }
#endif

    /* Allocate read buffer and send buffer. Once data are read from file into
     * rd_buf, it is unpacked into send_buf for each non-aggregator. send_buf
     * will be directly used to send the read request data to non-aggregators.
     *
     * Note rd_amnt may not be the same as send_amnt, as there can be overlaps
     * between adjacent offset-length pairs after sort.
     *
     * If file offset-length pairs have not been re-ordered, i.e. sorted and
     * overlaps removed, and this aggregator will not send any read data to its
     * non-aggregators, then we can use user's buffer, buf, to call MPI-IO/GIO
     * to read from the file, without allocating an additional temporary
     * buffer.
     */
    if (!do_sort && buf_view.size == send_amnt && !overlap) {
        rd_buf_view = buf_view;
        rd_buf = buf;
    }
    else {
        /* Read data will be stored in a contiguous read buffer. */
        rd_buf_view.size = rd_amnt;
        rd_buf_view.count = 0;
        if (rd_amnt > 0)
            rd_buf = (char*) NCI_Malloc(rd_amnt);
    }

#if 1
    int coll_indep = (fIsSet(ncp->flags, NC_MODE_INDEP)) ? NC_REQ_INDEP : NC_REQ_COLL;

MPI_Count buf_off=0, buf_len=rd_buf_view.size;
if (rd_buf_view.count == 0 && rd_buf_view.size > 0) {
    rd_buf_view.count = 1;
    rd_buf_view.off = &buf_off;
    rd_buf_view.len = &buf_len;
}
    MPI_Offset rlen = ncmpio_file_read(ncp, coll_indep, rd_buf, file_view, rd_buf_view);
    if (rlen < 0) {
        if (status == NC_NOERR) status = (int)rlen;
        rd_amnt = 0;
    }
#else
    err = ncmpio_read_write(ncp, NC_REQ_RD, 0, rd_buf_view, rd_buf);
    if (status == NC_NOERR) status = err;
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[3] = MAX(pnc_ina_mem_get[3], mem_max);

    endT = MPI_Wtime();
    pnc_ina_get[2] += endT - startT;
    startT = endT;
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
        if (buf_view.count <= 1)
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
        if (bufLen > 0 && buf_view.count > 1) {
            char *buf_ptr=tmp_buf;
            for (j=0; j<buf_view.count; j++) {
                memcpy((char*)buf + buf_view.off[j], buf_ptr, buf_view.len[j]);
                buf_ptr += buf_view.len[j];
            }
            NCI_Free(tmp_buf);
        }
    }

    if (nprocs == 1)
        /* In this case, communication will not be necessary. */
        goto fn_exit;

    /* Aggregators start sending read data to non-aggregators. At first,
     * allocate array_of_blocklengths[] and array_of_displacements[]
     */
    for (max_npairs=0, i=1; i<nprocs; i++)
        max_npairs = MAX(meta[3*i], max_npairs);

#ifdef HAVE_MPI_LARGE_COUNT
    blks = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * max_npairs);
    disps = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * max_npairs);
#else
    blks = (int*) NCI_Malloc(sizeof(int) * max_npairs);
    disps = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * max_npairs);
#endif

    /* Now, send data to each non-aggregator */
    req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nprocs);
    nreqs = 0;
    off_start = meta[0];
    for (i=1; i<nprocs; i++) {
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
            TRACE_COMM(MPI_Isend_c)(MPI_BOTTOM, 1, sendType, i, 0, intra_comm,
                                    &req[nreqs++]);
#else
            TRACE_COMM(MPI_Isend)(MPI_BOTTOM, 1, sendType, i, 0, intra_comm,
                                  &req[nreqs++]);
#endif
            MPI_Type_free(&sendType);
        }
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[4] = MAX(pnc_ina_mem_get[4], mem_max);

    endT = MPI_Wtime();
    pnc_ina_get[3] += endT - startT;
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
    /* offsets[] and lengths[] are used in GIO read subroutines as flattened
     * filetype. They cannot be freed before the I/O is done.
     */
    if (rd_buf != NULL && rd_buf != buf) NCI_Free(rd_buf);
    if (orig_lengths != NULL) NCI_Free(orig_lengths);
    if (orig_offsets != NULL) NCI_Free(orig_offsets);
    if (req != NULL) NCI_Free(req);
    if (meta != NULL) NCI_Free(meta);

    /* free space allocated for file_view and buf_view */
    if (npairs > 0) {
        NCI_Free(offsets);
        NCI_Free(lengths);
    }

    if (buf_view.count > 0) {
        NCI_Free(buf_view.off);
        NCI_Free(buf_view.len);
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[5] = MAX(pnc_ina_mem_get[5], mem_max);

    endT = MPI_Wtime();
    pnc_ina_get[4] += endT - startT;
#endif

    return status;
}

/*----< req_compare() >------------------------------------------------------*/
/* This subroutine is used to sort the string file offsets of reqs[] */
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
 *
 * Note even when INA is disabled, this subroutine is still called by all
 * nonblocking put/get APIs.
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
    PNCIO_View file_view, buf_view;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif

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
        /* check if offset_start are in a monotonic non-decreasing order */
        if (i > 0 && reqs[i].offset_start < reqs[i-1].offset_start)
            descreasing = 1;
    }

    /* If a decreasing order is found, sort reqs[] based on reqs[].offset_start
     * into an increasing order.
     */
    if (descreasing)
        qsort(reqs, (size_t)num_reqs, sizeof(NC_req), req_compare);

    /* construct file_view, the file access offset-length pairs
     *   file_view.count: total number of offset_length pairs
     *   file_view.off[]: array of file offsets
     *   file_view.len[]: array of lengths
     *   is_incr: whether file_view.off[] are incremental
     */
    err = flatten_nreqs(ncp, reqMode, num_reqs, reqs, &is_incr, &file_view);
    if (status == NC_NOERR) status = err;

    /* Note offsets lengths may contain overlaps between consecutive pairs when
     * the user's requests contain overlaps. They, if exist, will be resolved
     * later in ina_put() and ina_get().
     */

    /* Populate buf_view, which contains metadata describing the user buffer
     * from the nonblocking requests.
     */
    err = flat_buf_type(ncp, reqMode, num_reqs, reqs, &buf_view, &buf);
    if (status == NC_NOERR) status = err;

#if PNETCDF_DEBUG_MODE == 1
    if (num_reqs > 0) assert(buf_view.count > 0);
#endif

    if (req_list != NULL)
        /* All metadata in req_list have been used to construct bufType and
         * bufLen. It is now safe to release the space occupied by req_list.
         */
        NCI_Free(req_list);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    pnc_ina_flatten += MPI_Wtime() - timing;
#endif

    MPI_Comm saved_ina_intra_comm;
    saved_ina_intra_comm = ncp->comm_attr.ina_intra_comm;
    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP)) {
        /* Temporarily set ncp->comm_attr.ina_intra_comm to be as if self rank
         * is an INA aggregator and the INA group size is 1.
         */
        ncp->comm_attr.ina_intra_comm = MPI_COMM_SELF;
    }

    /* perform intra-node aggregation */
    if (fIsSet(reqMode, NC_REQ_WR))
        err = ina_put(ncp, is_incr, file_view.count, file_view.off, file_view.len, buf_view, buf);
    else
        err = ina_get(ncp, is_incr, file_view.count, file_view.off, file_view.len, buf_view, buf);
    if (status == NC_NOERR) status = err;

    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP)) {
        /* restore ncp->comm_attr.ina_intra_comm */
        ncp->comm_attr.ina_intra_comm = saved_ina_intra_comm;
    }

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
 *
 * Note even when INA is disabled, this subroutine is still called by all
 * blocking put/get APIs.
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
    PNCIO_View file_view, buf_view;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif

#warning TODO: fix view count

    if (buf_len == 0 || buf == NULL) {
        /* This is a zero-length request. When in collective data mode, this
         * rank must still participate collective calls. When INA is enabled,
         * this rank tells its aggregator that it has no I/O data. When INA is
         * disabled, this rank must participate other collective file call.
         */
        file_view.count = 0;
        buf_view.count = 0;
        buf_view.size = 0;
    }
    else {
        /* construct file_view, the file access offset-length pairs
         *   file_view.count: total number of offset_length pairs
         *   file_view.off[]: array of file offsets
         *   file_view.len[]: array of lengths
         *   is_incr: whether file_view.off[] are incremental
         */
        err = flatten_req(ncp, varp, start, count, stride, &is_incr,
                          &file_view);
        if (err == NC_NOERR) {
            /* buffer passed to blocking APIs is always contiguous */
            buf_view.count = 1;
            buf_view.off = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset));
            buf_view.len = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset));
            buf_view.off[0] = 0;
            buf_view.len[0] = buf_len;
            buf_view.size = buf_len;
        }
        else { /* make this rank zero-sized request */
            file_view.count = 0;
            buf_view.count = 0;
            buf_view.size = 0;
        }
        status = err;
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    pnc_ina_flatten += MPI_Wtime() - timing;
#endif

    MPI_Comm saved_ina_intra_comm;
    saved_ina_intra_comm = ncp->comm_attr.ina_intra_comm;
    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP)) {
        /* Temporarily set ncp->comm_attr.ina_intra_comm to be as if self rank
         * is an INA aggregator and the INA group size is 1.
         */
        ncp->comm_attr.ina_intra_comm = MPI_COMM_SELF;
    }

    /* perform intra-node aggregation */
    if (fIsSet(reqMode, NC_REQ_WR)) {
        err = ina_put(ncp, is_incr, file_view.count, file_view.off, file_view.len, buf_view, buf);
        if (status == NC_NOERR) status = err;
    }
    else {
        err = ina_get(ncp, is_incr, file_view.count, file_view.off, file_view.len, buf_view, buf);
        if (status == NC_NOERR) status = err;
    }

    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP)) {
        /* restore ncp->comm_attr.ina_intra_comm */
        ncp->comm_attr.ina_intra_comm = saved_ina_intra_comm;
    }

    return status;
}

