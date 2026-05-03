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

/* swap elements of arrays x, y, and corresponding lengths and bufAddr */
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    MPI_Offset tmp; \
    MPI_Offset d0 = (x) - offsets; \
    MPI_Offset d1 = (y) - offsets; \
    if (d0 != d1) { \
        SWAP1(*(x), *(y), tmp); \
        SWAP1(lengths[d0], lengths[d1], tmp); \
        if (bufAddr != NULL) \
            SWAP1(bufAddr[d0], bufAddr[d1], tmp); \
    } \
}

#define MEDIAN(a,b,c) ((*(a) < *(b)) ? \
                      ((*(b) < *(c)) ? (b) : ((*(a) < *(c)) ? (c) : (a))) : \
                      ((*(b) > *(c)) ? (b) : ((*(a) < *(c)) ? (a) : (c))))


/*----< qsort_off_len_buf() >------------------------------------------------*/
/* Sort three arrays of offsets, lengths, and buffer addresses based on array
 * offsets into an increasing order. This code is based on the qsort routine
 * from Bentley & McIlroy's "Engineering a Sort Function".
 */
static void
qsort_off_len_buf(MPI_Offset  num,
                  MPI_Offset *offsets,
                  MPI_Offset *lengths,
                  MPI_Offset *bufAddr)
{
    MPI_Offset *pa, *pb, *pc, *pd, *pl, *pm, *pn, cmp_result, swap_cnt, i, r;

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
        for (i=0; i<r; i++)
            SWAP(offsets, lengths, bufAddr, offsets+i, pb-r+i)

        r = MIN(pd - pc, pn - pd - 1);
        for (i=0; i<r; i++)
            SWAP(offsets, lengths, bufAddr, pb+i, pn-r+i)

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
 *
 * count[] contains the number of offset-blklens-bufAddr triplets per process
 * offsets[] contains offsets of nprocs processes appending one after another
 * blklens[] contains lengths of nprocs processes
 * bufAddr[] contains buffer addresses of nprocs processes
 * nelems is the total number of offset-blklens-bufAddr triplets
 */
static
void heap_merge(int              nprocs,
                const MPI_Offset *count,    /* [nprocs] */
                MPI_Offset      *offsets,  /* IN/OUT: [nelems] */
                MPI_Offset      *blklens,  /* IN/OUT: [nelems] */
                MPI_Offset      *bufAddr)  /* IN/OUT: [nelems] */
{
    typedef struct {
        MPI_Offset *off_list;
        MPI_Offset *len_list;
        MPI_Offset *addr_list;
        MPI_Offset   count;
    } heap_struct;

    heap_struct *a, tmp;
    int i, j, heapsize, l, r, k, smallest;
    size_t sum;
    MPI_Offset nelems, *srt_addr=NULL, *srt_off, *srt_len;

    for (nelems=0, i=0; i<nprocs; i++) nelems += count[i];

    /* This heap_merge is not in-place, taking too much memory footprint */
    srt_off = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * nelems);
    srt_len = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * nelems);

    if (bufAddr != NULL)
        srt_addr = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * nelems);

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

    memcpy(offsets, srt_off, sizeof(MPI_Offset) * nelems);
    memcpy(blklens, srt_len, sizeof(MPI_Offset) * nelems);
    if (bufAddr != NULL)
        memcpy(bufAddr, srt_addr, sizeof(MPI_Offset) * nelems);

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
                 MPI_Offset         *npairs,     /* OUT: num of off-len pairs */
                 MPI_Offset        *offsets,    /* OUT: array of offsets */
                 MPI_Offset        *lengths)    /* OUT: array of lengths */
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
flatten_req(const NC          *ncp,
            const NC_var      *varp,
            const MPI_Offset  *start,
            const MPI_Offset  *count,
            const MPI_Offset  *stride,
            int               *is_incr,   /* OUT: are offsets incrementing */
            PNCIO_View        *file_view) /* OUT: flattened file view */
{
    int i, j, err=NC_NOERR, ndims;
    MPI_Offset num, idx, var_begin, prev_end_off;
    MPI_Offset *shape, count0, *ones=NULL;

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

    file_view->off = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) *
                                              file_view->count);
    file_view->len = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) *
                                              file_view->count);

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
                               &num, /* OUT: num of off-len pairs */
                               file_view->off+idx,  /* OUT: array of offsets */
                               file_view->len+idx); /* OUT: array of lengths */
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
flatten_nreqs(const NC      *ncp,
              int            reqMode,   /* IN: NC_REQ_RD or NC_REQ_WR */
              int            num_reqs,  /* IN: # requests */
              const NC_req  *reqs,      /* [num_reqs] requests */
              int           *is_incr,   /* OUT: are offsets incrementing */
              PNCIO_View    *file_view) /* OUT: flattened file view */
{
    int i, j, err=NC_NOERR, ndims, max_ndims=0;
    MPI_Offset num, idx, prev_end_off;
    MPI_Offset *start, *count, *shape, *stride, *ones=NULL;

    file_view->count = 0; /* total number of offset-length pairs */
    file_view->off = NULL;
    file_view->len = NULL;

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
        err = flatten_subarray(ndims, lead->varp->xsz, var_begin, shape,
                               start, count, (stride == NULL) ? ones : stride,
                               &num, /* OUT: number of off-len pairs */
                               file_view->off+idx,  /* OUT: array of offsets */
                               file_view->len+idx); /* OUT: array of lengths */
        if (err != NC_NOERR)
            goto err_out;

        /* check if file_view->off[] are in an increasing order */
        for (j=0; j<num; j++) {
            if (prev_end_off > file_view->off[idx+j])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_end_off = file_view->off[idx+j];
        }
        idx += num;
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
 * metadata from the non-INA aggregators into meta[], including:
 *   1. the number of offset-length pairs of each INA group member
 *   2. offsets array of each INA group member
 *   3. lengths array of each INA group member
 *   4. INA aggregator's file_view will be updated:
 *      INA aggregator's file_view->count will be the total number of
 *      offset-length pairs of this INA group.
 *      INA aggregator's file_view->off will contain all the vile_view.off of
 *      this INA group (appending).
 *      INA aggregator's file_view->len will contain all the vile_view.len of
 *      this INA group (appending).
 */
static
int ina_collect_md(const NC   *ncp,
                   MPI_Offset *meta,
                   PNCIO_View *file_view)
{
    int i, err, rank, nprocs, mpireturn, status=NC_NOERR, nreqs;
    MPI_Request *req=NULL;
    MPI_Offset num_pairs=meta[0];
    MPI_Comm comm = ncp->comm_attr.ina_intra_comm;

#if PNETCDF_DEBUG_MODE == 1
    assert(comm != MPI_COMM_NULL);
#endif

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    /* Aggregator collects metadata of all INA group members */
    if (rank == 0) {
        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nprocs);
        nreqs = 0;
        for (i=1; i<nprocs; i++)
            TRACE_COMM(MPI_Irecv)(meta + i*3, 3, MPI_OFFSET, i, 0, comm,
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
        TRACE_COMM(MPI_Send)(meta, 3, MPI_OFFSET, 0, 0, comm);

    /* Secondly, aggregators collect offset-length pairs from all its INA group
     * members.
     */
    if (rank == 0) {
        MPI_Aint aint;
        MPI_Offset bklens[2], disps[2];
        MPI_Datatype recvType;

        /* calculate the total number of offset-length pairs to receive */
        for (i=1; i<nprocs; i++) file_view->count += meta[i*3];

        /* Realloc this INA aggregator's file_view->off and file_view->len to
         * receive the non-INA aggregators' file_view so they can be in a
         * contiguous buffer.
         */
        if (file_view->count > num_pairs) {
            file_view->off = (MPI_Offset*) NCI_Realloc(file_view->off,
                             sizeof(MPI_Offset) * file_view->count);
            file_view->len = (MPI_Offset*) NCI_Realloc(file_view->len,
                             sizeof(MPI_Offset) * file_view->count);
        }

        /* To minimize number of MPI recv calls per non-aggregator, below
         * creates a derived datatype, recvType, to combine file_view->off and
         * file_view->len into one MPI_Irecv call.
         */
        nreqs = 0;
        MPI_Get_address(file_view->off, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Offset) * meta[0]);
        MPI_Get_address(file_view->len, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(MPI_Offset) * meta[0]);
        for (i=1; i<nprocs; i++) {
            if (meta[i*3] == 0) continue;
            bklens[0] = meta[i*3] * sizeof(MPI_Offset);
            bklens[1] = meta[i*3] * sizeof(MPI_Offset);
            err = ncmpio_type_create_hindexed(2, disps, bklens, &recvType);
            if (status == NC_NOERR) status = err;

            /* post to receive offset-length pairs from non-aggregators */
            TRACE_COMM(MPI_Irecv)(MPI_BOTTOM, 1, recvType, i, 0, comm,
                                  &req[nreqs++]);
            MPI_Type_free(&recvType);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
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
    }
    else if (num_pairs > 0) { /* non-aggregator */
        /* To minimize number of MPI send calls to the aggregator, below
         * creates a derived datatype, sendType, to combine file_view->off and
         * file_view->len into one MPI_Send call.
         */
        MPI_Aint aint;
        MPI_Offset bklens[2], disps[2];
        MPI_Datatype sendType;

        bklens[0] = meta[0] * sizeof(MPI_Offset);
        bklens[1] = bklens[0];
        MPI_Get_address(file_view->off, &aint);
        disps[0] = aint;
        MPI_Get_address(file_view->len, &aint);
        disps[1] = aint;
        err = ncmpio_type_create_hindexed(2, disps, bklens, &sendType);
        if (status == NC_NOERR) status = err;

        TRACE_COMM(MPI_Send)(MPI_BOTTOM, 1, sendType, 0, 0, comm);
        MPI_Type_free(&sendType);
    }

    return status;
}

/*----< ina_put() >----------------------------------------------------------*/
/* This subroutine implements the intra-node aggregation for write operations.
 * It also handles the case when INA is disabled. Note heap space allocated in
 * file_view and buf_view will be freed by end of this subroutine.
 */
static
int ina_put(NC         *ncp,
            int         is_incr,   /* if file_view.off[] are incremental */
            PNCIO_View  file_view,
            PNCIO_View  buf_view,
            void       *buf)       /* user write buffer */
{
    char *recv_buf=NULL, *wr_buf = NULL, *mpi_name;
    int i, j, err, mpireturn, status=NC_NOERR, rank, nprocs, coalesceable=0;
    MPI_Offset saved_file_view_count, *meta=NULL;
    MPI_Offset wr_amnt=0, *bufAddr=NULL, *saved_file_view_len;
    MPI_Comm intra_comm;

#if PNETCDF_PROFILING == 1
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

    /* Each aggregator's first step is to collect metadata from all INA group
     * members about their request's file offset-length pairs, write amount,
     * and whether the offsets are in an incremental order. This step is
     * necessary for the next step which is to collect all members' file
     * offset-length pairs.
     *
     * Once ina_collect_md() returns, this INA aggregator's file_view.count
     * increases, file_view.off and file_view.len have grown to include the
     * ones from all the INA group members (appending one after another).
     *
     * For write operation, keeping the original offset-length pairs is not
     * necessary, as they will later be sorted and coalesced before calling
     * MPI-IO or GIO file write.
     */
    if (rank == 0)
        meta = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * nprocs * 3);
    else
        meta = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * 3);

    meta[0] = file_view.count;
    meta[1] = buf_view.size;
    meta[2] = is_incr;

    /* The INA aggregator collects request metadata from all its INA group
     * members. Note only the INA aggregator's file_view will get updated
     * after returned from ina_collect_md().
     */
    if (nprocs > 1) {
        err = ina_collect_md(ncp, meta, &file_view);
        if (err != NC_NOERR) {
            NCI_Free(meta);
            return err;
        }
    }

    /* For write operation, the non-aggregator now can start sending their
     * write data to its aggregator.
     */
    if (rank > 0) { /* non-aggregator */
        if (buf_view.count > 0) {
            MPI_Datatype sendType=MPI_BYTE;

            /* Non-INA aggregator sends write data to its INA aggregator */
            /* If this rank has non-zero sized request, it constructs an MPI
             * derived datatype and call MPI_send to send write data from its
             * INA aggregator.
             */
            if (buf_view.count > 1) {
                err = ncmpio_type_create_hindexed(buf_view.count, buf_view.off,
                                                  buf_view.len, &sendType);
                if (status == NC_NOERR)
                    status = err;

                /* When err != NC_NOERR, sendType is set to MPI_DATATYPE_NULL
                 * which will trigger an error when calling MPI_Send().
                 */
                TRACE_COMM(MPI_Send)(buf, 1, sendType, 0, 0, intra_comm);
                mpi_name = "MPI_Send";
            }
            else { /* buf_view.count == 1 */
#ifdef HAVE_MPI_LARGE_COUNT
                TRACE_COMM(MPI_Send_c)(buf, buf_view.size, MPI_BYTE, 0, 0,
                                       intra_comm);
                mpi_name = "MPI_Send_c";
#else
                int num;
                if (buf_view.size <= INT_MAX)
                    num = (int)buf_view.size;
                else {
                    num = 1;
                    err = ncmpio_type_contiguous(buf_view.size, &sendType);
                    if (status == NC_NOERR)
                        status = err;
                }

                /* When err != NC_NOERR, sendType is set to MPI_DATATYPE_NULL
                 * which will trigger an error when calling MPI_Send().
                 */
                TRACE_COMM(MPI_Send)(buf, num, sendType, 0, 0, intra_comm);
                mpi_name = "MPI_Send";
#endif
            }

            if (mpireturn != MPI_SUCCESS && status == NC_NOERR)
                status = ncmpii_error_mpi2nc(mpireturn, mpi_name);

            if (sendType != MPI_BYTE && sendType != MPI_DATATYPE_NULL)
                MPI_Type_free(&sendType);

            /* free space allocated for buf_view */
            NCI_Free(buf_view.off);
            NCI_Free(buf_view.len);
        }

        /* free space allocated for file_view */
        if (file_view.count > 0) {
            NCI_Free(file_view.off);
            NCI_Free(file_view.len);
        }

        /* Non-INA aggregators are done here, as only INA aggregators make a
         * call to ncmpio_file_write() to write data to the file. Non-INA
         * aggregators do not participate file I/O calls, except when
         * performing independent writes, in which case intra_comm is
         * temporarily set to MPI_COMM_SELF.
         */
        NCI_Free(meta);
        return status;
    }

    /* The remaining of this subroutine is for aggregators only -------------*/

#if PNETCDF_PROFILING == 1
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[1] = MAX(pnc_ina_mem_put[1], mem_max);

    endT = MPI_Wtime();
    pnc_ina_put[0] += endT - startT; /* collect MD */
    startT = endT;
#endif

    /* Once an INA aggregator collected metadata from all its INA group
     * members, it must construct new file view and buffer view whose
     * offset-length pairs are coalesced, with overlaps removed, and abide by
     * the MPI standard requirement on offsets being in a monotonically
     * non-decreasing order.
     *
     * MPI-IO has the following requirements about filetype.
     * 1. The (flattened) displacements (of a filetype) are not required to be
     *    distinct, but they cannot be negative, and they must be monotonically
     *    non-decreasing.
     * 2. If the file is opened for writing, neither the etype nor the filetype
     *    is permitted to contain overlapping regions.
     */

    if (file_view.count == 0) goto do_write;

    /* Now this aggregator has received all offset-length pairs from its
     * non-aggregators. If this INA group makes a non-zero sized request, the
     * first step is to check if a sorting of file offsets is necessary.
     */
    char *ptr;
    int nreqs, indv_sorted, do_sort, overlap;
    MPI_Request *req=NULL;
    MPI_Offset recv_amnt;

    /* Check whether or not all INA group members' file_view.off[] are
     * individually sorted.
     */
    indv_sorted = 1;
    do_sort = 0;
    for (i=-1,j=0; j<nprocs; j++) {
        if (i == -1 && meta[j*3] > 0)
            /* i is the 1st member whose file_view.count > 0 */
            i = j;
        if (meta[j*3+2] == 0) {
            /* member j's file_view.off are not sorted */
            indv_sorted = 0;
            do_sort = 1;
            break;
        }
    }
    /* i is the first INA group member whose file_view.count > 0, and
     * j is the first INA group member whose is_incr is false
     */

    if (i >= 0 && indv_sorted == 1) {
        /* Even when file_view.off[] of all INA group members are individually
         * sorted, we still need to check if offsets are interleaved. If
         * interleaved, we must sort all offset-length pairs.
         */
        MPI_Offset prev_end_off, sum;

        assert(meta[i*3+2] == 1);

        sum = meta[i*3];

        /* prev_end_off is the last offset of INA group member i */
        prev_end_off = file_view.off[sum-1];

        /* check if the file_view.off are interleaved */
        for (++i; i<nprocs; i++) {
            if (meta[i*3] == 0) /* zero-sized request */
                continue;

            assert(meta[i*3+2] == 1);

            /* file_view.off[sum] is member i' 1st offset */
            if (prev_end_off > file_view.off[sum]) {
                do_sort = 1; /* indicate file_view.off are not incrementing */
                break;
            }
            /* move on to the next member */
            sum += meta[i*3];
            prev_end_off = file_view.off[sum-1];
        }
    }

    /* Construct an array of buffer addresses containing a mapping of the
     * buffer used to receive write data from non-aggregators and the buffer
     * used to write to file. bufAddr[] is calculated based on the assumption
     * that the write buffer of this aggregator is contiguous, i.e.
     * buf_view.count <= 1. For non-aggregators, their write data will always
     * be received into a contiguous buffer.
     */
    bufAddr = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * file_view.count);
    bufAddr[0] = 0;
    for (i=1; i<file_view.count; i++)
        bufAddr[i] = bufAddr[i-1] + file_view.len[i-1];

    if (do_sort) {
        /* Sort file_view.off, file_view.len, bufAddr altogether, based on
         * file_view.off into an increasing order. Note during sorting, the
         * length and buffer must also be moved together with their
         * corresponding offset.
         */
        if (indv_sorted) {
            /* Interleaved offsets are found but individual offsets are already
             * sorted. This is commonly seen from the checkerboard domain
             * partitioning pattern. In this case, heap_merge() is faster to
             * merge all offsets into one single sorted offset list.  Note
             * count[] must be initialized, so it can be used in heap_merge()
             */
            MPI_Offset *count;
            count = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * nprocs);
            for (i=0; i<nprocs; i++)
                count[i] = meta[i*3];

            /* heap-merge() runs much faster than qsort() when individual lists
             * have already been sorted. However, it has a much bigger memory
             * footprint.
             */
            heap_merge(nprocs, count, file_view.off, file_view.len, bufAddr);
            NCI_Free(count);
        }
        else
            /* When some individual file_view.off are not sorted, we cannot use
             * heap_merge(). Note qsort() is an in-place sorting.
             */
            qsort_off_len_buf(file_view.count, file_view.off, file_view.len,
                              bufAddr);
    }

    /* Now file_view.off and file_view.len are sorted, but overlaps may exist
     * between adjacent pairs. If this is the case, they must be coalesced.
     *
     * The loop below checks if there is an overlap and calculates recv_amnt
     * and wr_amnt.
     * recv_amnt is the total amount this aggregator will receive from its INA
     *      group members, including self. recv_amnt includes overlaps.
     * wr_amnt is recv_amnt with overlap removed.
     *
     * This loop also coalesces offset-length pairs as well as the
     * corresponding buffer addresses, so they can be used to move write data
     * around in the final write buffer.
     *
     * Note file_view.off[] has been sorted into a monotonically non-decreasing
     * order. During the sorting, bufAddr[] are moved around based on their
     * corresponding file_view.off[], and thus bufAddr[] may not be in a
     * monotonically non-decreasing order.
     */
    coalesceable = 0;
    overlap = 0;
    wr_amnt = recv_amnt = file_view.len[0];
    for (i=0, j=1; j<file_view.count; j++) {
        recv_amnt += file_view.len[j];
        if (file_view.off[i] + file_view.len[i] >=
            file_view.off[j] + file_view.len[j]) {
            /* segment i completely covers segment j, skip j */
            overlap = 1;
            continue;
        }

        MPI_Offset gap = file_view.off[i] + file_view.len[i] - file_view.off[j];
        if (gap >= 0) { /* overlap detected, merge j into i */
            /* when gap > 0,  pairs i and j overlap
             * when gap == 0, pairs i and j are contiguous
             */
            if (gap > 0) overlap = 1;
            wr_amnt += file_view.len[j] - gap; /* subtract overlapped amount */
            if (bufAddr[i] + file_view.len[i] == bufAddr[j] + gap) {
                /* buffers i and j are contiguous, merge j into i and
                 * subtract overlapped amount.
                 */
                file_view.len[i] += file_view.len[j] - gap;
            }
            else { /* buffers are not contiguous, reduce j's len */
                coalesceable = 1;
                file_view.off[i+1] = file_view.off[j] + gap;
                file_view.len[i+1] = file_view.len[j] - gap;
                bufAddr[i+1] = bufAddr[j] + gap;
                i++;
            }
        }
        else { /* i and j do not overlap */
            wr_amnt += file_view.len[j];
            i++;
            if (i < j) {
                file_view.off[i] = file_view.off[j];
                file_view.len[i] = file_view.len[j];
                bufAddr[i] = bufAddr[j];
            }
        }
    }

    /* Now file_view.off[], file_view.len[], bufAddr[] are coalesced and no
     * overlap. Update file_view.count.
     */
    file_view.count = i+1;

    /* If file_view can be further coalesced, a new set of offsets and lengths
     * must be allocated for file_view. These new offsets and lengths cannot be
     * used for buf_view, because the buffer addresses may not be coalesceable
     * even the corresponding file_view can. Thus the old offsets must be kept
     * to construct buf_view.
     *
     * Note file_view.len can be updated in place, because it will not be used
     * by buf_view.
     */

    /* buf_view to be used in a call to ncmpio_file_write() later will use
     * bufAddr[] and file_view.len[], as it offset-length pairs. Because
     * file_view.count and file_view.len may change below if coalesceable is
     * true, we now save them for later use.
     */
    saved_file_view_len = file_view.len;
    saved_file_view_count = file_view.count;

    if (coalesceable) { /* file_view can be further coalesced */
        size_t cpy_amnt;
        MPI_Offset *file_len;

        cpy_amnt = sizeof(MPI_Offset) * file_view.count;
        file_len = (MPI_Offset*) NCI_Malloc(cpy_amnt);
        memcpy(file_len, file_view.len, cpy_amnt);
        file_view.len = file_len;

        for (i=0, j=1; j<file_view.count; j++) {
#if PNETCDF_DEBUG_MODE == 1
            /* any overlap should have been removed from the loop above */
            assert(file_view.off[i] + file_view.len[i] <= file_view.off[j]);
#endif
            if (file_view.off[i] + file_view.len[i] == file_view.off[j])
                /* coalesce j into i */
                file_view.len[i] += file_view.len[j];
            else { /* i and j are not coalesceable */
                i++;
                if (i < j) {
                    file_view.off[i] = file_view.off[j];
                    file_view.len[i] = file_view.len[j];
                }
            }
        }

        /* update number of offset-length pairs */
        file_view.count = i+1;
    }

#if PNETCDF_DEBUG_MODE == 1
    /* check if file_view's offset-lengths have been coalesced */
    for (i=1; i<file_view.count; i++) {
        assert(file_view.len[i-1] > 0);
        assert(file_view.off[i-1] < file_view.off[i]);
        assert(file_view.off[i-1] + file_view.len[i-1] < file_view.off[i]);
    }
#endif

#if PNETCDF_PROFILING == 1
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[2] = MAX(pnc_ina_mem_put[2], mem_max);
    pnc_ina_npairs_put = MAX(pnc_ina_npairs_put, file_view.count);

    endT = MPI_Wtime();
    pnc_ina_put[1] += endT - startT; /* sorting */
    startT = endT;
#endif

    /* Allocate receive buffer. Once write data from non-aggregators have
     * received into recv_buf, it is packed into wr_buf. Then, wr_buf is
     * used to call MPI-IO/GIO file write. Note the wr_buf is always
     * contiguous.
     *
     * When nprocs == 1, wr_buf is set to buf and together with buf_view are
     * directly passed to ncmpio_file_write() call.
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
        mpi_name = "MPI_Irecv_c";
#else
        int num;
        MPI_Datatype recvType=MPI_BYTE;

        if (meta[i*3 + 1] <= INT_MAX)
            num = (int)meta[i*3 + 1];
        else {
            num = 1;
            err = ncmpio_type_contiguous(meta[i*3 + 1], &recvType);
            if (status == NC_NOERR)
                status = err;
        }

        /* When err != NC_NOERR, recvType is set to MPI_DATATYPE_NULL
         * which will trigger an error when calling MPI_Recv().
         */
        TRACE_COMM(MPI_Irecv)(ptr, num, recvType, i, 0, intra_comm,
                              &req[nreqs++]);
        mpi_name = "MPI_Irecv";

        if (recvType != MPI_BYTE && recvType != MPI_DATATYPE_NULL)
            MPI_Type_free(&recvType);
#endif
        if (mpireturn != MPI_SUCCESS && status == NC_NOERR)
            status = ncmpii_error_mpi2nc(mpireturn, mpi_name);

        ptr += meta[i*3 + 1];
    }

#if PNETCDF_PROFILING == 1
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

#if PNETCDF_PROFILING == 1
    endT = MPI_Wtime();
    pnc_ina_put[3] += endT - startT; /* wait */
    startT = endT;
#endif

    /* Now all write data has been collected into recv_buf. In case of any
     * overlap, we must coalesce recv_buf into wr_buf using file_view.off[],
     * file_view.len[], and bufAddr[]. For overlapped regions, requests with
     * lower j indices win the writes to the overlapped regions.
     *
     * In case the user buffer, buf, can not be used to write to the file, loop
     * below packs recv_buf, data received from non-aggregators, into wr_buf, a
     * contiguous buffer, wr_buf, which will later be used in a call to
     * ncmpio_file_write().
     */
    if (!do_sort && wr_amnt == recv_amnt) {
        wr_buf = recv_buf;

        if (wr_buf != buf) {
            /* Since write data has been packed in wr_buf, a contiguous buffer,
             * update buf_view before passing it to ncmpio_file_write().
             */
            if (buf_view.count > 0) {
                NCI_Free(buf_view.len);
                NCI_Free(buf_view.off);
            }
            buf_view.off = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset));
            buf_view.len = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset));
            buf_view.off[0] = 0;
            buf_view.len[0] = wr_amnt;
            buf_view.count = 1;
            buf_view.size = wr_amnt;
        }
#if PNETCDF_DEBUG_MODE == 1
        else {
            /* User's buffer, buf, can be used to write to the file. In this
             * case, it also means all non-INA aggregators in this group have
             * zero-sized request.
             */
            for (i=1; i<nprocs; i++)
                assert(meta[3*i+1] == 0);
        }
#endif
    }
    else if (buf_view.count <= 1 && !overlap) {
        /* Note we can reuse bufAddr[] and file_view.len[] (before it is
         * coalesced) as buf_view.off and buf_view.len only when buf_view is
         * contiguous, because bufAddr[] is constructed based on the assumption
         * that the write buffer is contiguous.
         */
        wr_buf = recv_buf;

        /* Update buf_view before passing it to ncmpio_file_write(). */
        if (buf_view.count > 0) {
            NCI_Free(buf_view.len);
            NCI_Free(buf_view.off);
        }
        buf_view.off   = bufAddr; /* based on recv_buf */
        buf_view.len   = saved_file_view_len;
        buf_view.count = saved_file_view_count;
        buf_view.size  = wr_amnt;
    }
    else {
        /* do_sort == 1 means buffer's offsets and lengths have been moved
         * around to make file_view.off[] monotonically non-decreasing. In this
         * case, we need to re-arrange the write buffer accordingly by copying
         * write data into a temporary buffer, wr_buf, and write it to the
         * file.
         *
         * Note when npairs and wr_amnt are large, copying write data into a
         * contiguous buffer can be expensive, making INA cost high. Although
         * it makes the two-phase I/O MPI-IO and GIO run faster, this memory
         * copy cost may not be worthy. Besides, the memory footprint
         * high-water mark is doubled.
         */
        wr_buf = NCI_Malloc(wr_amnt);
        ptr = wr_buf;

        /* Copy write data into wr_buf, a contiguous buffer. */
        for (j=0; j<saved_file_view_count; j++) {
            memcpy(ptr, recv_buf + bufAddr[j], saved_file_view_len[j]);
            ptr += saved_file_view_len[j];
        }

        /* saved_file_view_len can now be freed, if it is != file_view.len */
        if (saved_file_view_len != file_view.len) NCI_Free(saved_file_view_len);

        /* Update buf_view before passing it to ncmpio_file_write(). */
        if (buf_view.count > 0) {
            NCI_Free(buf_view.len);
            NCI_Free(buf_view.off);
        }
        buf_view.len = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset));
        buf_view.off = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset));
        buf_view.off[0] = 0;
        buf_view.len[0] = wr_amnt;
        buf_view.count = 1;
        buf_view.size = wr_amnt;

        if (recv_buf != buf) NCI_Free(recv_buf);
    }

do_write:
    NCI_Free(meta);

#if PNETCDF_PROFILING == 1
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[4] = MAX(pnc_ina_mem_put[4], mem_max);

    endT = MPI_Wtime();
    pnc_ina_put[4] += endT - startT; /* setview */
    startT = endT;
#endif

    /* carry out write request to file */
    int coll_indep = (fIsSet(ncp->flags, NC_MODE_INDEP)) ? NC_REQ_INDEP
                                                         : NC_REQ_COLL;

    MPI_Offset wlen;

    wlen = ncmpio_file_write(ncp, coll_indep, wr_buf, file_view, buf_view);
    if (wlen < 0) {
        if (status == NC_NOERR) status = (int)wlen;
        wr_amnt = 0;
    }

    if (wr_buf != buf) NCI_Free(wr_buf);

    /* Free bufAddr if it is not used by buf_view.off */
    if (bufAddr != NULL && bufAddr != buf_view.off) NCI_Free(bufAddr);

    /* free space allocated for file_view and buf_view */
    if (file_view.count > 0) {
        NCI_Free(file_view.off);
        /* file_view.len and buf_view.len may share the same address */
        if (file_view.len != buf_view.len) NCI_Free(file_view.len);
    }
    if (buf_view.count > 0) {
        NCI_Free(buf_view.off);
        NCI_Free(buf_view.len);
    }

#if PNETCDF_PROFILING == 1
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_put[5] = MAX(pnc_ina_mem_put[5], mem_max);

    endT = MPI_Wtime();
    pnc_ina_put[5] += endT - startT; /* write */
#endif

    return status;
}

/*----< bin_search() >-------------------------------------------------------*/
static
size_t bin_search(MPI_Offset        key,
                  MPI_Offset        nmemb,
                  const MPI_Offset *base)
{
    MPI_Offset low, high;

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
 * It also handles the case when INA is disabled. Note heap space allocated in
 * file_view and buf_view will be freed by end of this subroutine.
 */
static
int ina_get(NC         *ncp,
            int         is_incr,   /* if file_view.off[] are incremental */
            PNCIO_View  file_view,
            PNCIO_View  buf_view,
            void       *buf)       /* user read buffer */
{
    char *rd_buf = NULL;
    int i, j, err, mpireturn, status=NC_NOERR, nprocs, rank, nreqs;
    int do_sort=0, indv_sorted=1, overlap=0;
    MPI_Offset *meta=NULL, *blks=NULL, *disps=NULL;
    MPI_Offset max_npairs, send_amnt=0, rd_amnt=0, off_start;
    MPI_Request *req=NULL;
    MPI_Comm intra_comm;
    PNCIO_View rd_buf_view, orig_fview;

#if PNETCDF_PROFILING == 1
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

    /* Each aggregator's first step is to collect metadata from all INA group
     * members about their request's file offset-length pairs, write amount,
     * and whether the offsets are in an incremental order. This step is
     * necessary for the next step which is to collect all members' file
     * offset-length pairs.
     *
     * Once ina_collect_md() returns, this INA aggregator's file_view.count
     * increases, file_view.off and file_view.len have grown to include the
     * ones from all the INA group members (appending one after another).
     *
     * For read operation, the original file offset-length pairs must be kept,
     * as they are required to unpack file read data into messages for sending
     * them to the INA group members. The buf_view's offset-length pairs to be
     * passed to ncmpio_file read() will be modified to be sorted into an
     * incremental order and coalesced.
     */
    if (rank == 0)
        meta = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * nprocs * 3);
    else
        meta = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * 3);

    meta[0] = file_view.count;
    meta[1] = buf_view.size;
    meta[2] = is_incr;

    /* This INA aggregator's file_view.count must be preserved, so it can be
     * used to unpack read data into its user buffer. Below the call to
     * ina_collect_md() will increase/expand this INA aggregator's
     * file_view.count, file_view.off[] and file_view.len[].
     */
    orig_fview.count = file_view.count;

    /* The INA aggregator collects request metadata from all its INA group
     * members. Note only the INA aggregator's file_view will get updated
     * after returned from ina_collect_md().
     */
    if (nprocs > 1) {
        err = ina_collect_md(ncp, meta, &file_view);
        if (err != NC_NOERR) {
            NCI_Free(meta);
            return err;
        }
    }

    /* For read operation, the non-INA aggregator members now can start
     * receiving their read data from its aggregator.
     */
    if (rank > 0) { /* not an INA aggregator */
        if (buf_view.count > 0) {
            char *mpi_name;
            MPI_Status st;
            MPI_Datatype recvType=MPI_BYTE;

            /* If this rank has non-zero sized request, it constructs an MPI
             * derived datatype and call MPI_recv to receive read data from its
             * INA aggregator.
             */
            if (buf_view.count > 1) {
                err = ncmpio_type_create_hindexed(buf_view.count, buf_view.off,
                                                  buf_view.len, &recvType);
                if (status == NC_NOERR)
                    status = err;

                /* When err != NC_NOERR, recvType is set to MPI_DATATYPE_NULL
                 * which will trigger an error when calling MPI_Send().
                 */
                TRACE_COMM(MPI_Recv)(buf, 1, recvType, 0, 0, intra_comm, &st);
                mpi_name = "MPI_Recv";
            }
            else { /* buf_view.count == 1 */
#ifdef HAVE_MPI_LARGE_COUNT
                TRACE_COMM(MPI_Recv_c)(buf, buf_view.size, MPI_BYTE, 0, 0,
                                       intra_comm, &st);
                mpi_name = "MPI_Recv_c";
#else
                int num;
                if (buf_view.size <= INT_MAX)
                    num = (int)buf_view.size;
                else {
                    num = 1;
                    err = ncmpio_type_contiguous(buf_view.size, &recvType);
                    if (status == NC_NOERR)
                        status = err;
                }

                /* When err != NC_NOERR, recvType is set to MPI_DATATYPE_NULL
                 * which will trigger an error when calling MPI_Recv().
                 */
                TRACE_COMM(MPI_Recv)(buf, num, recvType, 0, 0, intra_comm, &st);
                mpi_name = "MPI_Recv";
#endif
            }

            if (mpireturn != MPI_SUCCESS && status == NC_NOERR)
                status = ncmpii_error_mpi2nc(mpireturn, mpi_name);

            if (recvType != MPI_BYTE && recvType != MPI_DATATYPE_NULL)
                MPI_Type_free(&recvType);

            /* free space allocated for buf_view */
            NCI_Free(buf_view.off);
            NCI_Free(buf_view.len);
        }

        /* free space allocated for file_view */
        if (file_view.count > 0) {
            NCI_Free(file_view.off);
            NCI_Free(file_view.len);
        }

        /* Non-INA aggregators are now done, as they do not participate MPI-IO
         * or GIO file read (neither collective nor independent).
         */
        NCI_Free(meta);
        return status;
    }

    /* The remaining of this subroutine is for INA aggregators only. */

#if PNETCDF_PROFILING == 1
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[1] = MAX(pnc_ina_mem_get[1], mem_max);

    endT = MPI_Wtime();
    pnc_ina_get[0] += endT - startT; /* collect MD */
    startT = endT;
#endif

    if (file_view.count == 0) {
        /* This INA aggregation group has zero data to read, but this
         * aggregator must participate the collective I/O calls.
         */
        goto do_read;
    }

    /* For read operations, the INA aggregator's current file_view.off[] and
     * file_view.len[] collected from all its INA group members must be kept
     * untouched, because the later sorting and coalescing will be performed on
     * file_view, messing up the original ones are needed to construct MPI
     * datatype for the INA aggregator to send the data read from file to its
     * INA group members.
     */
    size_t alloc_amnt = sizeof(MPI_Offset) * file_view.count;
    orig_fview.off = (MPI_Offset*) NCI_Malloc(alloc_amnt);
    orig_fview.len = (MPI_Offset*) NCI_Malloc(alloc_amnt);
    memcpy(orig_fview.off, file_view.off, alloc_amnt);
    memcpy(orig_fview.len, file_view.len, alloc_amnt);

#if PNETCDF_PROFILING == 1
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[2] = MAX(pnc_ina_mem_get[2], mem_max);
#endif

    /* Once an INA aggregator collected metadata from all its INA group
     * members, it must construct new file view and buffer view whose
     * offset-length pairs are coalesced, with overlaps removed, and abide by
     * the MPI standard requirement on offsets being in a monotonically
     * non-decreasing order.
     *
     * MPI-IO has the following requirements about filetype.
     * 1. The (flattened) displacements (of a filetype) are not required to be
     *    distinct, but they cannot be negative, and they must be monotonically
     *    non-decreasing.
     * 2. If the file is opened for writing, neither the etype nor the filetype
     *    is permitted to contain overlapping regions.
     */

    /* Now this INA aggregator has received all offset-length pairs from
     * its non-aggregators. At first, it checks if a sorting is necessary.
     *
     * First check whether or not all INA group members' file_view.off[] are
     * individually sorted.
     */
    indv_sorted = 1;
    for (i=-1,j=0; j<nprocs; j++) {
        if (i == -1 && meta[j*3] > 0)
            /* i is the 1st member whose file_view.count > 0 */
            i = j;
        if (meta[j*3+2] == 0) {
            /* member j's file_view.off are not sorted */
            indv_sorted = 0;
            do_sort = 1;
            break;
        }
    }
    /* i is the first INA group member whose file_view.count > 0, and
     * j is the first INA group member whose is_incr is false
     */

    if (i >= 0 && indv_sorted == 1) {
        /* Even when file_view.off[] of all INA group members are individually
         * sorted, we still need to check if offsets are interleaved. If
         * interleaved, we must sort all offset-length pairs.
         */
        MPI_Offset prev_end_off, sum;

        assert(meta[i*3+2] == 1);

        sum = meta[i*3];

        /* prev_end_off is the last offset of INA group member i */
        prev_end_off = file_view.off[sum-1];

        /* check if the file_view.off are interleaved */
        for (++i; i<nprocs; i++) {
            if (meta[i*3] == 0) /* zero-sized request */
                continue;

            assert(meta[i*3+2] == 1);

            if (prev_end_off > file_view.off[sum]) {
                /* file_view.off[sum] is the non-aggregator i' 1st offset */
                do_sort = 1; /* file_view.off are not incrementing */
                break;
            }
            /* move on to the next member */
            sum += meta[i*3];
            prev_end_off = file_view.off[sum-1];
        }
    }

    if (do_sort) {
        /* Sort file_view.off[] into an increasing order. Note during sorting,
         * file_view.len[] is also moved together with their corresponding
         * file_view.off[].
         */
        if (indv_sorted) {
            /* Interleaved offsets are found in the aggregated file_view.off[],
             * but individual file_view.off[] are already sorted. This is
             * commonly seen from the checkerboard domain partitioning pattern.
             * In this case, heap_merge() is faster to merge all file_view.off
             * into one single sorted offset list. Note count[] must be
             * initialized, so it can be used in heap_merge()
             */
            MPI_Offset *count;
            count = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * nprocs);
            for (i=0; i<nprocs; i++)
                count[i] = meta[i*3];

            /* heap-merge() runs much faster than qsort() when individual
             * lists have already been sorted. However, it has a much
             * bigger memory footprint.
             */
            heap_merge(nprocs, count, file_view.off, file_view.len, NULL);
            NCI_Free(count);
        }
        else
            /* When some individual file_view.off are not already sorted, we
             * cannot use heap_merge(). Note qsort() is an in-place sorting.
             */
            qsort_off_len_buf(file_view.count, file_view.off, file_view.len,
                              NULL);
    }

    /* Coalesce the file_view.off[] and file_view.len[] and calculate the total
     * read amount and the total send amounts to all INA group members by this
     * aggregator.
     */
    overlap = 0;
    send_amnt = rd_amnt = file_view.len[0];
    for (i=0, j=1; j<file_view.count; j++) {
        MPI_Offset gap;
        send_amnt += file_view.len[j];

        gap = file_view.off[i] + file_view.len[i] - file_view.off[j];
        if (gap >= 0) { /* overlap detected, merge j into i */
            /* when gap > 0,  pairs i and j overlap
             * when gap == 0, pairs i and j are contiguous
             */
            MPI_Offset i_end, j_end;

            if (gap > 0) overlap = 1;

            i_end = file_view.off[i] + file_view.len[i];
            j_end = file_view.off[j] + file_view.len[j];
            if (i_end < j_end) {
                file_view.len[i] += j_end - i_end;
                rd_amnt += j_end - i_end;
            }
            /* else: j is entirely covered by i */
        }
        else { /* j and i are not overlapped */
            rd_amnt += file_view.len[j];
            i++;
            if (i < j) {
                file_view.off[i] = file_view.off[j];
                file_view.len[i] = file_view.len[j];
            }
        }
    }

    /* update file_view.count after coalesce */
    file_view.count = i+1;

#if PNETCDF_PROFILING == 1
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[2] = MAX(pnc_ina_mem_get[2], mem_max);
    pnc_ina_npairs_get = MAX(pnc_ina_npairs_get, file_view.count);

    endT = MPI_Wtime();
    pnc_ina_get[1] += endT - startT; /* sorting */
    startT = endT;

#endif

do_read:

    /* Allocate read buffer and send buffer. Once data are read from file into
     * rd_buf. MPI derived datatypes will be constructed for sending the read
     * data to each INA group member. rd_buf will be directly used to send so
     * no additional memory allocation is necessary.
     *
     * Note rd_amnt may not be the same as send_amnt, as there can be overlaps
     * between adjacent file offset-length pairs after sorting, with overlaps
     * removed.
     *
     * If this INA aggregator does not have to send any read data to its INA
     * group members, i.e. all its non-INA aggregator members have zero-sized
     * requests, and the file offset-length pairs have not been re-ordered,
     * i.e. sorted and overlaps removed, then we can use user's buffer, buf, to
     * call ncmpio_file_read() to read from the file, without allocating an
     * additional temporary buffer.
     */
    if (!do_sort && buf_view.size == send_amnt && !overlap) {
        /* Only this INA aggregator has non-zero sized data to read and
         * its buf_view does not required to be sorted and overlaps removed.
         */
        rd_buf_view = buf_view;
        rd_buf = buf;
    }
    else {
        /* Allocate a read buffer to store data read from the file. */
        if (rd_amnt > 0)
            rd_buf = (char*) NCI_Malloc(rd_amnt);

        rd_buf_view.count = 1;
        rd_buf_view.off = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset));
        rd_buf_view.len = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset));
        rd_buf_view.off[0] = 0;
        rd_buf_view.len[0] = rd_amnt;
        rd_buf_view.size = rd_amnt;
    }

    int coll_indep = (fIsSet(ncp->flags, NC_MODE_INDEP)) ? NC_REQ_INDEP
                                                         : NC_REQ_COLL;

    MPI_Offset rlen;

    rlen = ncmpio_file_read(ncp, coll_indep, rd_buf, file_view, rd_buf_view);
    if (rlen < 0) {
        if (status == NC_NOERR) status = (int)rlen;
        rd_amnt = 0;
    }

    if (file_view.count == 0) {
        /* This INA aggregation group has zero data to read. */
        if (meta != NULL) NCI_Free(meta);
        return status;
    }

#if PNETCDF_PROFILING == 1
    ncmpi_inq_malloc_size(&mem_max);
    // ncmpi_inq_malloc_max_size(&mem_max);
    pnc_ina_mem_get[3] = MAX(pnc_ina_mem_get[3], mem_max);

    endT = MPI_Wtime();
    pnc_ina_get[2] += endT - startT;
    startT = endT;
#endif

    /* If sorting has been performed, the orders of file_view.off[] and
     * file_view.len[] may no longer be the same as the original ones. We must
     * use binary search to find the offset-length pair from file_view that
     * contains each INA group member's offset-length pair to construct a send
     * buffer datatype, a view layout to the read buffer, rd_buf, so the data
     * can be directly sent from rd_buf.
     */
    if (rd_buf != buf) {
        /* First, the INA aggregator takes case of its own read, by coping the
         * read data to its user buffer. Note file_view.off[] has been sorted
         * in an incremental order.
         *
         * When the offset-length pairs of read buffer have been sorted or the
         * read buffer size is smaller than the total read amount, we must
         * search and copy from read buffer to self's user buffer.
         */
        char *ptr=NULL, *tmp_buf=NULL;
        size_t m=0, k, scan_off=0;

        /* If this INA aggregator's user buffer is contiguous, then we have
         * used buf as the file read buffer. If not, allocate a temporary
         * buffer, copy the read data over, and then unpacking it to the user
         * buffer.
         */
        if (buf_view.count <= 1)
            ptr = buf;
        else if (buf_view.size > 0)
            ptr = tmp_buf = (char*) NCI_Malloc(buf_view.size);

        for (j=0; j<orig_fview.count; j++) {
            /* This loop unpacks read data into this INA aggregator's user
             * buffer only. Note the first 'orig_fview.count' elements of
             * orig_fview.off[] and orig_fview.len[] are this INA aggregator's
             * own request.
             *
             * Note file_view have been modified and now describes the layout
             * of rd_buf.
             *
             * For each orig_fview's offset-length pair j, find in file_view
             * the offset-length pair in rd_buf covering it. Note that if
             * orig_fview's offset-length pairs are not already sorted, i.e.
             * is_incr != 1, this bin_search() below can be very expensive!
             *
             * orig_fview.off[] and orig_fview.len[] are the original aggregated
             *     file_view.off[] and file_view.len[] of this INA aggregator.
             *     And, the first orig_fview.count of elements of offset-length
             *     pairs are this aggregator's own requests. orig_fview is no
             *     longer the same as file_view.off[] and file_view.len[].
             * file_view.off[] and file_view.len[] describe the offset-length
             *     pairs of the read buffer, rd_buf. It is a modified
             *     orig_fview with a sorting and coalescing applied.
             */
            if (!is_incr) m = 0;

            if (file_view.count-m == 1)
                assert(file_view.off[m] <= orig_fview.off[j]);

            k = bin_search(orig_fview.off[j], file_view.count - m,
                           &file_view.off[m]);

            assert(k < file_view.count);

            /* k returned from bin_search is relative to m */
            k += m;

            /* When is_incr is 1, the orig_fview.off[] are in an incremental
             * order and we can continue binary search using the index from
             * previous search. When is_incr is 0, the orig_fview.off[] are
             * NOT in an incremental order, we must do binary search on the
             * entire file_view.off[].
             */
            if (!is_incr) scan_off = 0;
            for (; m<k; m++)
                scan_off += file_view.len[m];

            /* Note orig_fview.off[j] and orig_fview.len[j] must entirely
             * covered by file_view.off[k] and file_view.len[k], because
             * file_view.off[] and file_view.len[] have been coalesced.
             */
            memcpy(ptr,
                   rd_buf + (scan_off + orig_fview.off[j] - file_view.off[k]),
                   orig_fview.len[j]);

            ptr += orig_fview.len[j];
        }

        /* unpack read data to user read buffer, if not done already */
        if (buf_view.size > 0 && buf_view.count > 1) {
            char *buf_ptr=tmp_buf;
            for (j=0; j<buf_view.count; j++) {
                memcpy((char*)buf + buf_view.off[j], buf_ptr, buf_view.len[j]);
                buf_ptr += buf_view.len[j];
            }
            NCI_Free(tmp_buf);
        }
    } /* Done with coping read data to this INA aggregator's own  buffer. */

    if (nprocs == 1)
        /* In this case, communication will not be necessary. */
        goto fn_exit;

    /* This INA aggregator starts preparing sending read data to its INA group
     * members. At first, construct an MPI derived datatype to be used in
     * MPI_Isend(), by allocating and setting array_of_blocklengths[] and
     * array_of_displacements[]
     */
    for (max_npairs=0, i=1; i<nprocs; i++)
        max_npairs = MAX(meta[3*i], max_npairs);

    blks  = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * max_npairs);
    disps = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * max_npairs);

    /* Now, send data to each INA group member */
    req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * nprocs);
    nreqs = 0;
    off_start = meta[0];
    for (i=1; i<nprocs; i++) {
        /* populate disps[] and blks[] */
        size_t k, m, scan_off;
        MPI_Offset remote_num_pairs, remote_is_incr, *off, *len;

        remote_num_pairs = meta[3*i];
        remote_is_incr = meta[3*i+2];

        if (remote_num_pairs == 0) continue; /* zero sized request */

        off = orig_fview.off + off_start;
        len = orig_fview.len + off_start;

        m = 0;
        scan_off = 0;
        for (j=0; j<remote_num_pairs; j++) {
            MPI_Aint addr;

            /* Find the offset-length pair in rd_buf containing this pair j.
             * Note that if the INA group member i's file_view's offset-length
             * pairs are not already sorted, i.e. remote_is_incr == 1, this
             * bin_search() below can be very expensive!
             */
            if (!remote_is_incr) m = 0;

            if (file_view.count-m == 1)
                assert(file_view.off[m] <= off[j]);

            k = bin_search(off[j], file_view.count-m, &file_view.off[m]);
            /* k returned from bin_search is relative to m */
            k += m;

            assert(file_view.off[k] <= off[j] &&
                   off[j] < file_view.off[k] + file_view.len[k]);

            /* When is_incr is 1, the orig_fview.off[] are in an incremental
             * order, we can continue binary search using the index from the
             * previous search. When is_incr is 0, the orig_fview.off[] are NOT
             * in an incremental order, we must do binary search on the entire
             * file_view.off[].
             */
            if (!remote_is_incr) scan_off = 0;

            for (; m<k; m++)
                scan_off += file_view.len[m];

            /* Note orig_fview.off[j] and file_view.len[j] must entirely
             * covered by file_view.off[k] and file_view.len[k], because
             * file_view.off[] and file_view.len[] have been coalesced.
             */
            char *ptr = rd_buf + (scan_off + off[j] - file_view.off[k]);
            MPI_Get_address(ptr, &addr);
            disps[j] = addr;
            blks[j] = len[j];
        }

        /* Construct a send buffer MPI datatype */
        MPI_Datatype sendType;
        err = ncmpio_type_create_hindexed(remote_num_pairs, disps, blks,
                                          &sendType);
        if (status == NC_NOERR) status = err;

        TRACE_COMM(MPI_Isend)(MPI_BOTTOM, 1, sendType, i, 0, intra_comm,
                              &req[nreqs++]);
        if (sendType != MPI_BYTE && sendType != MPI_DATATYPE_NULL)
            MPI_Type_free(&sendType);

        /* Move on to the next INA group member i */
        off_start += remote_num_pairs;
    }

#if PNETCDF_PROFILING == 1
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
    if (rd_buf != NULL && rd_buf != buf) {
        NCI_Free(rd_buf_view.len);
        NCI_Free(rd_buf_view.off);
        NCI_Free(rd_buf);
    }
    if (orig_fview.off != NULL) NCI_Free(orig_fview.off);
    if (orig_fview.len != NULL) NCI_Free(orig_fview.len);
    if (req != NULL) NCI_Free(req);
    if (meta != NULL) NCI_Free(meta);

    /* free space allocated for file_view and buf_view */
    if (file_view.count > 0) {
        NCI_Free(file_view.off);
        NCI_Free(file_view.len);
    }

    if (buf_view.count > 0) {
        NCI_Free(buf_view.off);
        NCI_Free(buf_view.len);
    }

#if PNETCDF_PROFILING == 1
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

#if PNETCDF_PROFILING == 1
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

    if (req_list != NULL) {
        /* All metadata in req_list have been used to set file_view and
         * buf_view. It is now safe to release the space occupied by req_list.
         */
        NCI_Free(req_list);
    }

#if PNETCDF_PROFILING == 1
    pnc_ina_flatten += MPI_Wtime() - timing;
#endif

    /* When a non-INA aggregator performs independent I/O, we need to
     * temporarily set ncp->comm_attr.ina_intra_comm to be MPI_COMM_SELF, as if
     * self rank is an INA aggregator and the INA group size is of 1.
     */
    MPI_Comm saved_ina_intra_comm;
    saved_ina_intra_comm = ncp->comm_attr.ina_intra_comm;
    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP))
        ncp->comm_attr.ina_intra_comm = MPI_COMM_SELF;

    /* Perform intra-node aggregation and file I/O. Note that the space
     * allocated in file_view and buf_view will be freed at the end of
     * ina_put() and ina_get().
     */
    if (fIsSet(reqMode, NC_REQ_WR))
        err = ina_put(ncp, is_incr, file_view, buf_view, buf);
    else
        err = ina_get(ncp, is_incr, file_view, buf_view, buf);
    if (status == NC_NOERR) status = err;

    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP))
        /* restore ncp->comm_attr.ina_intra_comm */
        ncp->comm_attr.ina_intra_comm = saved_ina_intra_comm;

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
               const NC_var     *varp,
               const MPI_Offset *start,
               const MPI_Offset *count,
               const MPI_Offset *stride,
               MPI_Offset        buf_len,
               void             *buf)
{
    int err, status=NC_NOERR, is_incr=1;
    PNCIO_View file_view, buf_view;

#if PNETCDF_PROFILING == 1
    double timing = MPI_Wtime();
#endif

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

#if PNETCDF_PROFILING == 1
    pnc_ina_flatten += MPI_Wtime() - timing;
#endif

    /* When a non-INA aggregator performs independent I/O, we need to
     * temporarily set ncp->comm_attr.ina_intra_comm to be MPI_COMM_SELF, as if
     * self rank is an INA aggregator and the INA group size is of 1.
     */
    MPI_Comm saved_ina_intra_comm;
    saved_ina_intra_comm = ncp->comm_attr.ina_intra_comm;
    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP))
        ncp->comm_attr.ina_intra_comm = MPI_COMM_SELF;

    /* Perform intra-node aggregation and file I/O. Note that the space
     * allocated in file_view and buf_view will be freed at the end of
     * ina_put() and ina_get().
     */
    if (fIsSet(reqMode, NC_REQ_WR))
        err = ina_put(ncp, is_incr, file_view, buf_view, buf);
    else
        err = ina_get(ncp, is_incr, file_view, buf_view, buf);

    if (status == NC_NOERR) status = err;

    if (ncp->num_aggrs_per_node == 0 || fIsSet(ncp->flags, NC_MODE_INDEP))
        /* restore ncp->comm_attr.ina_intra_comm */
        ncp->comm_attr.ina_intra_comm = saved_ina_intra_comm;

    return status;
}

