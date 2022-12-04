/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <common.h>
#include <math.h>
#include <mpi.h>
#include <ncchkio_driver.h>
#include <pnc_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ncchkio_internal.h"

int ncchkioi_get_var_cb_chunk (NC_chk *ncchkp,
							   NC_chk_var *varp,
							   const MPI_Offset *start,
							   const MPI_Offset *count,
							   const MPI_Offset *stride,
							   void *buf) {
	int err = NC_NOERR;
	int i, j, k;
	int cid;  // Chunk iterator

	MPI_Offset *ostart = NULL, *osize;
	int *tsize = NULL, *tssize, *tstart, *tsizep, *tssizep, *tstartp;  // Size for sub-array type
	MPI_Offset *citr;												   // Chunk iterator

	int *rcnt_local = NULL, *rcnt_all = NULL;  // Number of processes that writes to each chunk

	int overlapsize;	// Size of overlaping region of request and chunk
	char *cbuf = NULL;	// Intermediate continuous buffer

	int packoff;		 // Pack offset
	MPI_Datatype ptype;	 // Pack datatype

	int nread;		   // # chunks to read form file
	int *rids = NULL;  // Id of chunks to read from file

	int nsend, nrecv;							// Number of send and receive
	MPI_Request *sreqs = NULL, *rreqs = NULL;	// Send and recv req
	MPI_Status *sstats = NULL, *rstats = NULL;	// Send and recv status
	char **sbufs = NULL, **rbufs = NULL;		// Send and recv buffer
	int *rsizes = NULL;							// recv size of each message
	MPI_Message rmsg;							// Receive message

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_INIT)

	// Allocate buffering for write count
	rcnt_local = (int *)NCI_Malloc (sizeof (int) * varp->nchunk * 2);
	rcnt_all   = rcnt_local + varp->nchunk;

	// Allocate buffering for overlaping index
	tsize  = (int *)NCI_Malloc (sizeof (int) * varp->ndim * 3);
	tssize = tsize + varp->ndim;
	tstart = tssize + varp->ndim;
	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * 3);
	osize  = ostart + varp->ndim;

	// Chunk iterator
	citr = osize + varp->ndim;

	// We need to calculate the size of message of each chunk
	// This is just for allocating send buffer
	// We do so by iterating through all request and all chunks they cover
	// If we are not the owner of a chunk, we need to send message
	memset (rcnt_local, 0, sizeof (int) * varp->nchunk);
	nsend = 0;

	// Iterate through chunks
	ncchkioi_chunk_itr_init (varp, start, count, citr, &cid);
	do {
		rcnt_local[cid] = 1;

		if (varp->chunk_owner[cid] != ncchkp->rank) {
			// Count number of mnessage we need to send
			nsend++;
		}
	} while (ncchkioi_chunk_itr_next (varp, start, count, citr, &cid));

	NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_GET_CB_INIT)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SYNC)

	// Sync number of messages of each chunk
	CHK_ERR_ALLREDUCE (rcnt_local, rcnt_all, varp->nchunk, MPI_INT, MPI_SUM, ncchkp->comm);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SYNC)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_IO_INIT)

	// We need to prepare chunk in the chunk cache
	// For chunks not yet allocated, we need to read them form file collectively
	// We collect chunk id of those chunks
	// Calculate number of recv request
	// This is for all the chunks
	rids  = (int *)NCI_Malloc (sizeof (int) * varp->nmychunk);
	nread = 0;
	nrecv = 0;
	for (i = 0; i < varp->nmychunk; i++) {
		cid = varp->mychunks[i];
		// We don't need message for our own data
		nrecv += rcnt_all[cid] - rcnt_local[cid];
		// Count number of chunks we need to prepare
		// We read only chunks that is required

		if (rcnt_all[cid] || rcnt_local[cid]) {
			if (varp->chunk_cache[cid] == NULL) {
				// err = ncchkioi_cache_alloc(ncchkp, varp->chunksize, varp->chunk_cache + cid);
				if (varp->chunk_index[cid].len > 0) { rids[nread++] = cid; }
			} else {
				// ncchkioi_cache_visit(ncchkp, varp->chunk_cache[cid]);
			}
		}
	}

	NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_GET_CB)  // I/O time count separately

#ifdef PNETCDF_PROFILING
	MPI_Barrier (ncchkp->comm);
#endif
	// Decompress chunks into chunk cache
	err = ncchkioi_load_var (ncchkp, varp, nread, rids);
	CHK_ERR
	// Increase batch number to indicate allocated chunk buffer can be freed for future allocation
	(ncchkp->cache_serial)++;

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB)

	// Allocate buffer for send and recv
	// We need to accept nrecv requests and receive nsend of replies
	rreqs = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * (nrecv + nsend));
	CHK_PTR (rreqs)
	rstats = (MPI_Status *)NCI_Malloc (sizeof (MPI_Status) * (nrecv + nsend));
	CHK_PTR (rstats)
	rbufs = (char **)NCI_Malloc (sizeof (char *) * (nrecv + nsend));
	CHK_PTR (rbufs)
	rsizes = (int *)NCI_Malloc (sizeof (int) * (nrecv + nsend));
	CHK_PTR (rsizes)
	// We need to send nsend requests and reply nrecv of requests
	sbufs = (char **)NCI_Malloc (sizeof (char *) * (nrecv + nsend));
	CHK_PTR (sbufs)
	sreqs = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * (nrecv + nsend));
	CHK_PTR (sreqs)
	sstats = (MPI_Status *)NCI_Malloc (sizeof (MPI_Status) * (nrecv + nsend));
	CHK_PTR (sstats)

	// Post send
	k = 0;
	// Initialize chunk iterator
	ncchkioi_chunk_itr_init_ex (varp, start, count, citr, &cid, ostart, osize);
	// Iterate through chunks
	do {
		// We got something to send if we are not owner
		if (varp->chunk_owner[cid] != ncchkp->rank) {
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_PACK_REQ)

			// Calculate chunk overlap
			overlapsize = varp->esize;
			for (j = 0; j < varp->ndim; j++) { overlapsize *= osize[j]; }

			// Allocate buffer
			sbufs[k] = (char *)NCI_Malloc (sizeof (int) * varp->ndim * 2);	// For request
			CHK_PTR (sbufs[k])
			rbufs[k + nrecv] =
				(char *)NCI_Malloc (overlapsize);  // For reply, first nrecv are for request
			CHK_PTR (rbufs[k + nrecv])

			// Metadata
			tstartp = (int *)sbufs[k];
			packoff = varp->ndim * sizeof (int);
			tsizep	= (int *)(sbufs[k] + packoff);
			packoff += varp->ndim * sizeof (int);
			for (j = 0; j < varp->ndim; j++) {
				tstartp[j] = (int)(ostart[j] - citr[j]);
				tsizep[j]  = (int)osize[j];
			}

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_PACK_REQ)
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REQ)

			// Send request
			CHK_ERR_ISEND (sbufs[k], packoff, MPI_BYTE, varp->chunk_owner[cid], cid, ncchkp->comm,
						   sreqs + k);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REQ)
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REP)

			// Post recv reply
			CHK_ERR_IRECV (rbufs[k + nrecv], overlapsize, MPI_BYTE, varp->chunk_owner[cid],
						   cid + 1024, ncchkp->comm, rreqs + nrecv + k);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REP)

			k++;
		}
	} while (ncchkioi_chunk_itr_next_ex (varp, start, count, citr, &cid, ostart, osize));

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REQ)

	// Post recv
	k = 0;
	for (i = 0; i < varp->nmychunk; i++) {
		cid = varp->mychunks[i];
		// We are the owner of the chunk
		// Receive data from other process
		for (j = 0; j < rcnt_all[cid] - rcnt_local[cid]; j++) {
			// Get message size, including metadata
			CHK_ERR_MPROBE (MPI_ANY_SOURCE, cid, ncchkp->comm, &rmsg, rstats);
			CHK_ERR_GET_COUNT (rstats, MPI_BYTE, rsizes + k);

			// Allocate buffer
			rbufs[k] = (char *)NCI_Malloc (rsizes[k]);
			CHK_PTR (rbufs[k])

			// Post irecv
			CHK_ERR_IMRECV (rbufs[k], rsizes[k], MPI_BYTE, &rmsg, rreqs + k);
			k++;
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REQ)

	// Allocate intermediate buffer
	cbuf = (char *)NCI_Malloc (varp->chunksize);
	CHK_PTR (cbuf)

	// For each chunk we own, we need to receive incoming data
	k = 0;
	for (i = 0; i < varp->nmychunk; i++) {
		cid = varp->mychunks[i];

		NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SELF)

		// Handle our own data first if we have any
		if (rcnt_local[cid] > 0) {
			// Convert chunk id to iterator
			get_chunk_itr (varp, cid, citr);

			// Calculate overlapping region
			get_chunk_overlap (varp, citr, start, count, ostart, osize);

			// Pack type from chunk buffer to (contiguous) intermediate buffer
			for (j = 0; j < varp->ndim; j++) {
				tstart[j] = (int)(ostart[j] - citr[j]);
				tsize[j]  = varp->chunkdim[j];
				tssize[j] = (int)osize[j];
			}
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
										  varp->etype, &ptype);
			CHK_ERR_TYPE_COMMIT (&ptype);

			// Pack data into intermediate buffer
			packoff = 0;
			CHK_ERR_PACK (varp->chunk_cache[cid]->buf, 1, ptype, cbuf, varp->chunksize, &packoff,
						  ncchkp->comm);
			overlapsize = packoff;
			MPI_Type_free (&ptype);

			// Pack type from (contiguous) intermediate buffer to user buffer
			for (j = 0; j < varp->ndim; j++) {
				tstart[j] = (int)(ostart[j] - start[j]);
				tsize[j]  = (int)count[j];
			}
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
										  varp->etype, &ptype);
			CHK_ERR_TYPE_COMMIT (&ptype);

			// Pack data into user buffer
			packoff = 0;
			CHK_ERR_UNPACK (cbuf, overlapsize, &packoff, buf, 1, ptype, ncchkp->comm);
			MPI_Type_free (&ptype);
		}

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SELF)
		NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REQ)

		// Wait for all send requests related to this chunk
		// We remove the impact of -1 mark in rcnt_local[cid]
		// printf("Rank: %d, CHK_ERR_WAITALL_recv(%d, %d)\n", ncchkp->rank, rcnt_all[cid] -
		// rcnt_local[cid], k); fflush(stdout);
		CHK_ERR_WAITALL (rcnt_all[cid] - rcnt_local[cid], rreqs + k, rstats + k);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REQ)

		// Now, it is time to process data from other processes
		for (j = 0; j < varp->ndim; j++) { tsize[j] = varp->chunkdim[j]; }

		// Process data received
		// printf("nrecv = %d, rcnt_all = %d, rcnt_local = %d\n", nrecv, rcnt_all[cid],
		// rcnt_local[cid]); fflush(stdout);
		for (j = k; j < k + rcnt_all[cid] - rcnt_local[cid]; j++) {
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_UNPACK_REQ)

			// Metadata
			tstartp = (int *)rbufs[j];
			packoff = varp->ndim * sizeof (int);
			tssizep = (int *)(rbufs[j] + packoff);
			packoff += varp->ndim * sizeof (int);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_UNPACK_REQ)
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_PACK_REP)

			// Pack type
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssizep, tstartp, MPI_ORDER_C,
										  varp->etype, &ptype);
			CHK_ERR_TYPE_COMMIT (&ptype);

			// Allocate buffer
			MPI_Type_size (ptype, &overlapsize);
			sbufs[j + nsend] = (char *)NCI_Malloc (overlapsize);  // For reply

			// Data
			packoff = 0;
			CHK_ERR_PACK (varp->chunk_cache[cid]->buf, 1, ptype, sbufs[j + nsend], overlapsize,
						  &packoff, ncchkp->comm);
			MPI_Type_free (&ptype);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_PACK_REP)
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REP)

			// Send reply
			CHK_ERR_ISEND (sbufs[j + nsend], packoff, MPI_BYTE, rstats[j].MPI_SOURCE, cid + 1024,
						   ncchkp->comm, sreqs + j + nsend);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REP)
		}
		k += rcnt_all[cid] - rcnt_local[cid];

		// princbuf(ncchkp->rank, varp->chunk_cache[cid], varp->chunksize);
	}

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REQ)

	// Wait for all send request
	// printf("Rank: %d, CHK_ERR_WAITALL_send(%d, %d)\n", ncchkp->rank, nsend, 0); fflush(stdout);
	CHK_ERR_WAITALL (nsend, sreqs, sstats);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REQ)

	// Receive replies from the owners and update the user buffer
	k = 0;
	// Initialize chunk iterator
	ncchkioi_chunk_itr_init_ex (varp, start, count, citr, &cid, ostart, osize);
	// Iterate through chunks
	do {
		// We got something to recv if we are not owner
		if (varp->chunk_owner[cid] != ncchkp->rank) {
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_UNPACK_REP)

			// Pack type from recv buffer to user buffer
			for (j = 0; j < varp->ndim; j++) {
				tstart[j] = (int)(ostart[j] - start[j]);
				tsize[j]  = (int)count[j];
				tssize[j] = (int)osize[j];
			}
			// printf("Rank: %d, ostart=[%lld, %lld], osize=[%lld, %lld]\n", ncchkp->rank,
			// ostart[0], ostart[1], osize[0], osize[1]); fflush(stdout); printf("Rank: %d,
			// CHK_ERR_TYPE_CREATE_SUBARRAY4([%d, %d], [%d, %d], [%d, %d]\n", ncchkp->rank,
			// tsize[0], tsize[1], tssize[0], tssize[1], tstart[0], tstart[1]); fflush(stdout);
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
										  varp->etype, &ptype);
			// printf("Rank: %d, commit\n", ncchkp->rank); fflush(stdout);
			CHK_ERR_TYPE_COMMIT (&ptype);
			MPI_Type_size (ptype, &overlapsize);

			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REP)

			// printf("Rank: %d, wait recv, nrecv = %d, k = %d, nsend = %d\n", ncchkp->rank, nrecv,
			// k, nsend); fflush(stdout);
			// Wait for reply
			// printf("Rank: %d, MPI_Wait_recv(%d)\n", ncchkp->rank, nrecv + k); fflush(stdout);
			MPI_Wait (rreqs + nrecv + k, rstats + nrecv + k);

			NC_CHK_TIMER_STOPEX (NC_CHK_TIMER_GET_CB_RECV_REP, NC_CHK_TIMER_GET_CB_UNPACK_REP)

			// Pack data
			// printf("Rank: %d, pack\n", ncchkp->rank); fflush(stdout);
			packoff = 0;
			CHK_ERR_UNPACK (rbufs[nrecv + k], overlapsize, &packoff, buf, 1, ptype, ncchkp->comm);
			MPI_Type_free (&ptype);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_UNPACK_REP)

			k++;
		}
	} while (ncchkioi_chunk_itr_next_ex (varp, start, count, citr, &cid, ostart, osize));

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REP)

	// printf("Rank: %d, wait_final\n", ncchkp->rank); fflush(stdout);
	// Wait for all send replies
	// printf("Rank: %d, CHK_ERR_WAITALL_send(%d, %d)\n", ncchkp->rank, nrecv, nsend);
	// fflush(stdout);
	CHK_ERR_WAITALL (nrecv, sreqs + nsend, sstats + nsend);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REP)

	// printf("Rank: %d, exiting\n", ncchkp->rank); fflush(stdout);

err_out:;

	// Free buffers
	NCI_Free (rcnt_local);

	NCI_Free (rids);

	NCI_Free (tsize);

	NCI_Free (ostart);

	for (i = 0; i < nsend + nrecv; i++) {
		NCI_Free (sbufs[i]);
		NCI_Free (rbufs[i]);
	}
	NCI_Free (sreqs);
	NCI_Free (sstats);
	NCI_Free (sbufs);
	NCI_Free (rreqs);
	NCI_Free (rstats);
	NCI_Free (rbufs);
	NCI_Free (rsizes);

	if (cbuf != NULL) { NCI_Free (cbuf); }

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REP)

	return err;
}

int ncchkioi_get_var_cb_proc (NC_chk *ncchkp,
							  NC_chk_var *varp,
							  const MPI_Offset *start,
							  const MPI_Offset *count,
							  const MPI_Offset *stride,
							  void *buf) {
	int err = NC_NOERR;
	int i, j, k;
	int cid, cown;	// Chunk iterator

	MPI_Offset *ostart= NULL , *osize;
	int *tsize = NULL, *tssize, *tstart, *tssizep, *tstartp;  // Size for sub-array type
	MPI_Offset *citr;										  // Chunk iterator

	int *rcnt_local = NULL, *rcnt_all;	// Number of processes that writes to each proc

	int rrange_local[2], rrange_all[2];	 // Number of processes that writes to each chunk

	int overlapsize;	  // Size of overlaping region of request and chunk
	int max_tbuf = 0;	  // Size of intermediate buffer
	char *tbuf	 = NULL;  // Intermediate buffer

	int packoff;		 // Pack offset
	MPI_Datatype ptype;	 // Pack datatype

	int nread;		   // # chunks to read form file
	int *rids = NULL;  // Id of chunks to read from file

	int nsend, nrecv;									  // Number of send and receive
	MPI_Request *sreq = NULL, *rreq, *sreq_re, *rreq_re;  // Send and recv req
	MPI_Status *sstat = NULL, rstat, *sstat_re;			  // Send and recv status
	char **sbuf = NULL, **rbuf, **sbufp, **rbufp, **sbuf_re, **rbuf_re;	 // Send and recv buffer
	int *rsize, *ssize = NULL, *rsize_re, *ssize_re;  // recv size of each message
	int *sdst = NULL;								  // recv size of each message
	int *smap = NULL;
	MPI_Message rmsg;  // Receive message

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_INIT)

	// Allocate buffering for write count
	rcnt_local = (int *)NCI_Malloc (sizeof (int) * (ncchkp->np * 2 + varp->nchunk * 1));
    CHK_PTR(rcnt_local)
	rcnt_all   = rcnt_local + ncchkp->np;
	smap	   = rcnt_all + ncchkp->np;

	// Allocate buffering for overlaping index
	tsize  = (int *)NCI_Malloc (sizeof (int) * varp->ndim * 3);
    CHK_PTR(tsize)
	tssize = tsize + varp->ndim;
	tstart = tssize + varp->ndim;
	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * 3);
    CHK_PTR(ostart)
	osize  = ostart + varp->ndim;

	// Chunk iterator
	citr = osize + varp->ndim;

	// We need to calculate the size of message of each chunk
	// This is just for allocating send buffer
	// We do so by iterating through all request and all chunks they cover
	// If we are not the owner of a chunk, we need to send message
	memset (rcnt_local, 0, sizeof (int) * (ncchkp->np + varp->nchunk));
	nsend = 0;

	// Count total number of messages and build a map of accessed chunk to list of comm
	// datastructure
	rrange_local[0] = varp->nchunk;
	rrange_local[1] = 0;
	ncchkioi_chunk_itr_init (varp, start, count, citr, &cid);  // Initialize chunk iterator
	do {
		// Chunk owner
		cown = varp->chunk_owner[cid];

		// Mapping to skip list of send requests
		if (rcnt_local[cown] == 0 && cown != ncchkp->rank) { smap[cown] = nsend++; }
		rcnt_local[cown] = 1;  // Need to send message if not owner

		// Record lowest and highest chunk accessed
		if (rrange_local[0] > cid) { rrange_local[0] = cid; }
		if (rrange_local[1] < cid) { rrange_local[1] = cid; }
	} while (ncchkioi_chunk_itr_next (varp, start, count, citr, &cid));

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_INIT)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SYNC)

	// Sync number of messages of each chunk and access range
	CHK_ERR_ALLREDUCE (rcnt_local, rcnt_all, ncchkp->np, MPI_INT, MPI_SUM, ncchkp->comm);
	nrecv = rcnt_all[ncchkp->rank] -
			rcnt_local[ncchkp->rank];  // We don't need to receive request form self

	rrange_local[1] *= -1;
	CHK_ERR_ALLREDUCE (rrange_local, rrange_all, 2, MPI_INT, MPI_MIN, ncchkp->comm);
	rrange_all[1] *= -1;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SYNC)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_PACK_REQ)

	// Allocate data structure for messaging
	sbuf = (char **)NCI_Malloc (sizeof (char *) * (2 * nsend + nrecv));
	CHK_PTR (sbuf)
	ssize = (int *)NCI_Malloc (sizeof (int) * (nsend * 2 + nrecv * 1));
	CHK_PTR (ssize)
	sdst = ssize + (nsend + nrecv);
	sreq = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * (nsend + nrecv));
	CHK_PTR (sreq)
	sstat = (MPI_Status *)NCI_Malloc (sizeof (MPI_Status) * (nsend + nrecv));
	CHK_PTR (sstat)

	rbuf = (char **)NCI_Malloc (sizeof (char *) * (nsend + nrecv * 2));
	CHK_PTR (rbuf)
	rsize = (int *)NCI_Malloc (sizeof (int) * (nsend + nrecv));
	CHK_PTR (rsize)
	rreq = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * (nsend + nrecv));
	CHK_PTR (rreq)

	sbuf_re	 = sbuf + nsend;
	sbufp	 = sbuf_re + nrecv;
	ssize_re = ssize + nsend;
	sreq_re	 = sreq + nsend;
	sstat_re = sstat + nsend;

	rbuf_re	 = rbuf + nrecv;
	rbufp	 = rbuf_re + nsend;
	rsize_re = rsize + nrecv;
	rreq_re	 = rreq + nrecv;

	// Count size of each request
	memset (ssize, 0, sizeof (int) * nsend);
	memset (rsize_re, 0, sizeof (int) * nsend);
	ncchkioi_chunk_itr_init_ex (varp, start, count, citr, &cid, ostart,
								osize);	 // Initialize chunk iterator
	do {
		// Chunk owner
		cown = varp->chunk_owner[cid];
		if (cown != ncchkp->rank) {
			j		= smap[cown];
			sdst[j] = cown;	 // Record a reverse map by the way

			// Count overlap
			overlapsize = varp->esize;
			for (i = 0; i < varp->ndim; i++) { overlapsize *= osize[i]; }
			ssize[j] += sizeof (int) * (varp->ndim * 2 + 1);
			rsize_re[j] += overlapsize;
		}
	} while (ncchkioi_chunk_itr_next_ex (varp, start, count, citr, &cid, ostart, osize));

	// Allocate buffer for send
	for (i = 0; i < nsend; i++) {
		ssize[i] += sizeof (int);
		sbuf[i] = sbufp[i] = (char *)NCI_Malloc (ssize[i]);
		CHK_PTR (sbuf[i])
		*((int *)sbufp[i]) = rsize_re[i];
		sbufp[i] += sizeof (int);
		rbuf_re[i] = (char *)NCI_Malloc (rsize_re[i]);
		CHK_PTR (rbuf_re[i])
	}

	// Pack requests
	ncchkioi_chunk_itr_init_ex (varp, start, count, citr, &cid, ostart,
								osize);	 // Initialize chunk iterator
	do {
		// Chunk owner
		cown = varp->chunk_owner[cid];
		if (cown != ncchkp->rank) {
			j = smap[cown];

			// Metadata
			*((int *)sbufp[j]) = cid;
			sbufp[j] += sizeof (int);
			tstartp = (int *)sbufp[j];
			sbufp[j] += sizeof (int) * varp->ndim;
			tssizep = (int *)sbufp[j];
			sbufp[j] += sizeof (int) * varp->ndim;
			for (i = 0; i < varp->ndim; i++) {
				tstartp[i] = (int)(ostart[i] - citr[i]);
				tssizep[i] = (int)osize[i];
			}
		}
	} while (ncchkioi_chunk_itr_next_ex (varp, start, count, citr, &cid, ostart, osize));

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_PACK_REQ)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REQ)

	// Post send
	for (i = 0; i < nsend; i++) {
		CHK_ERR_ISEND (sbuf[i], ssize[i], MPI_BYTE, sdst[i], 0, ncchkp->comm, sreq + i);
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REQ)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REP)

	// Post receive
	for (i = 0; i < nsend; i++) {
		CHK_ERR_IRECV (rbuf_re[i], rsize_re[i], MPI_BYTE, sdst[i], 1, ncchkp->comm, rreq_re + i);
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REP)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REQ)

	// Post recv
	for (i = 0; i < nrecv; i++) {
		// Get message size, including metadata
		CHK_ERR_MPROBE (MPI_ANY_SOURCE, 0, ncchkp->comm, &rmsg, &rstat);
		CHK_ERR_GET_COUNT (&rstat, MPI_BYTE, rsize + i);

		// Allocate buffer
		rbuf[i] = rbufp[i] = (char *)NCI_Malloc (rsize[i]);
		CHK_PTR (rbuf[i])

		// Post irecv
		CHK_ERR_IMRECV (rbuf[i], rsize[i], MPI_BYTE, &rmsg, rreq + i);
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REQ)

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_IO_INIT)

	// We need to prepare chunk in the chunk cache
	// For chunks not yet allocated, we need to read them from file collectively
	// We collect chunk id of those chunks
	// Calculate number of recv request
	// This is for all the chunks
	for (j = 0; j < varp->nmychunk && varp->mychunks[j] < rrange_all[0]; j++)
		;
	for (k = j; k < varp->nmychunk && varp->mychunks[k] <= rrange_all[1]; k++)
		;
	rids  = (int *)NCI_Malloc (sizeof (int) * (k - j));
	nread = 0;
	for (i = j; i < k; i++) {
		cid = varp->mychunks[i];
						printf("checking chunk %d, size is %d\n",cid, varp->chunk_index[cid].len);
		if (varp->chunk_cache[cid] == NULL) {
			// err = ncchkioi_cache_alloc(ncchkp, varp->chunksize, varp->chunk_cache + cid);
			// varp->chunk_cache[cid] = (char*)NCI_Malloc(varp->chunksize);
			if (varp->chunk_index[cid].len > 0) { rids[nread++] = cid; printf("chunk %d need read\n",cid);}
		} else {
			// ncchkioi_cache_visit(ncchkp, varp->chunk_cache[cid]);
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_IO_INIT)

	NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_GET_CB)  // I/O time count separately

#ifdef PNETCDF_PROFILING
	MPI_Barrier (ncchkp->comm);
#endif
	// Decompress chunks into chunk cache
	err = ncchkioi_load_var (ncchkp, varp, nread, rids);
	CHK_ERR
	// Increase batch number to indicate allocated chunk buffer can be freed for future allocation
	(ncchkp->cache_serial)++;
	
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB)

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SELF)

	tbuf = (char *)NCI_Malloc (varp->chunksize);

	// Handle our own data
	ncchkioi_chunk_itr_init_ex (varp, start, count, citr, &cid, ostart,
								osize);	 // Initialize chunk iterator
	do {
		if (varp->chunk_owner[cid] == ncchkp->rank) {
			// Pack type from chunk cache to (contiguous) intermediate buffer
			for (j = 0; j < varp->ndim; j++) {
				tstart[j] = (int)(ostart[j] - citr[j]);
				tsize[j]  = varp->chunkdim[j];
				tssize[j] = (int)osize[j];
			}
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
										  varp->etype, &ptype);
			CHK_ERR_TYPE_COMMIT (&ptype);

			// Pack data into intermediate buffer
			packoff = 0;
			CHK_ERR_PACK (varp->chunk_cache[cid]->buf, 1, ptype, tbuf, varp->chunksize, &packoff,
						  ncchkp->comm);
			MPI_Type_free (&ptype);
			overlapsize = packoff;

			// Pack type from (contiguous) intermediate buffer to chunk buffer
			for (j = 0; j < varp->ndim; j++) {
				tstart[j] = (int)(ostart[j] - start[j]);
				tsize[j]  = (int)count[j];
			}
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
										  varp->etype, &ptype);
			CHK_ERR_TYPE_COMMIT (&ptype);

			// Unpack data into chunk buffer
			packoff = 0;
			CHK_ERR_UNPACK (tbuf, overlapsize, &packoff, buf, 1, ptype, ncchkp->comm);
			MPI_Type_free (&ptype);
		}
	} while (ncchkioi_chunk_itr_next_ex (varp, start, count, citr, &cid, ostart, osize));

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SELF)

	// Handle incoming requests
	for (i = 0; i < varp->ndim; i++) { tsize[i] = varp->chunkdim[i]; }
	for (i = 0; i < nrecv; i++) {
		NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REQ)

		// Will wait any provide any benefit?
		MPI_Waitany (nrecv, rreq, &j, &rstat);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REQ)

		packoff		= 0;
		ssize_re[j] = *((int *)rbufp[j]);
		rbufp[j] += sizeof (int);
		sbuf_re[j] = (char *)NCI_Malloc (ssize_re[j]);
		CHK_PTR (sbuf_re[j])
		while (rbufp[j] < rbuf[j] + rsize[j]) {
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_UNPACK_REQ)

			// Metadata
			cid = *((int *)rbufp[j]);
			rbufp[j] += sizeof (int);
			tstartp = (int *)rbufp[j];
			rbufp[j] += sizeof (int) * varp->ndim;
			tssizep = (int *)rbufp[j];
			rbufp[j] += sizeof (int) * varp->ndim;

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_UNPACK_REQ)
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_PACK_REP)

			// Pack type
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssizep, tstartp, MPI_ORDER_C,
										  varp->etype, &ptype);
			CHK_ERR_TYPE_COMMIT (&ptype);

			// Data
			CHK_ERR_PACK (varp->chunk_cache[cid]->buf, 1, ptype, sbuf_re[j], ssize_re[j], &packoff,
						  ncchkp->comm);
			MPI_Type_free (&ptype);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_PACK_REP)
		}

		NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REQ)

		// Send Response
		CHK_ERR_ISEND (sbuf_re[j], packoff, MPI_BYTE, rstat.MPI_SOURCE, 1, ncchkp->comm,
					   sreq_re + j);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REQ)
	}

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REQ)

	// Wait for all request
	CHK_ERR_WAITALL (nsend, sreq, sstat);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REQ)

	// Handle reply
	for (i = 0; i < varp->ndim; i++) { tsize[i] = count[i]; }
	for (i = 0; i < nsend; i++) {
		NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REP)

		// Will wait any provide any benefit?
		MPI_Waitany (nsend, rreq_re, &j, &rstat);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REP)
		NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_UNPACK_REP)

		sbufp[j] = sbuf[j] + sizeof (int);	// Skip reply size
		packoff	 = 0;
		while (packoff < rsize_re[j]) {
			// Retrieve metadata from the request we sent
			cid = *((int *)sbufp[j]);
			sbufp[j] += sizeof (int);
			tstartp = (int *)sbufp[j];
			sbufp[j] += sizeof (int) * varp->ndim;
			tssizep = (int *)sbufp[j];
			sbufp[j] += sizeof (int) * varp->ndim;

			// Bring back the request
			get_chunk_itr (varp, cid, citr);
			for (k = 0; k < varp->ndim; k++) { tstartp[k] += (int)(citr[k] - start[k]); }

			// Pack type
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssizep, tstartp, MPI_ORDER_C,
										  varp->etype, &ptype);
			CHK_ERR_TYPE_COMMIT (&ptype);

			// Pack data
			CHK_ERR_UNPACK (rbuf_re[j], rsize_re[j], &packoff, buf, 1, ptype, ncchkp->comm);
			MPI_Type_free (&ptype);
		}

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_UNPACK_REP)
	}

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REP)

	// Wait for all Response
	CHK_ERR_WAITALL (nrecv, sreq_re, sstat_re);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REP)

err_out:;

	// Free buffers
	NCI_Free (rcnt_local);

	NCI_Free (rids);

	NCI_Free (tsize);

	NCI_Free (ostart);

	NCI_Free (sreq);
	NCI_Free (sstat);
	NCI_Free (ssize);
	for (i = 0; i < nsend + nrecv; i++) {
		NCI_Free (sbuf[i]);
		NCI_Free (rbuf[i]);
	}
	NCI_Free (sbuf);

	NCI_Free (rreq);
	NCI_Free (rbuf);
	NCI_Free (rsize);

	if (tbuf != NULL) { NCI_Free (tbuf); }

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB)

	return err;
}
