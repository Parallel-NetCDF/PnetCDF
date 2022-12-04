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

int printbuf (int rank, char *buf, int len) {
	int i;
	printf ("Rank %d: ", rank);
	for (i = 0; i < len; i++) { printf ("%02x ", buf[i]); }
	printf ("\n");
}

int ncchkioi_put_varn_cb_chunk (NC_chk *ncchkp,
								NC_chk_var *varp,
								int nreq,
								MPI_Offset *const *starts,
								MPI_Offset *const *counts,
								MPI_Offset *const *strides,
								void **bufs) {
	int err = NC_NOERR;
	int i, j, k;
	int cid, req;  // Chunk and request iterator

	int *tsize, *tssize, *tstart, *tsizep, *tstartp;  // Size for sub-array type
	MPI_Offset *ostart, *osize;
	MPI_Offset *citr;

	int *wcnt_local, *wcnt_all;	 // Number of processes that writes to each chunk

	int nread;	// Chunks to read for background
	int *rids;

	int overlapsize;	// Size of overlaping region of request and chunk
	int max_tbuf;		// Size of intermediate buffer
	char *tbuf = NULL;	// Intermediate buffer

	int packoff;		 // Pack offset
	MPI_Datatype ptype;	 // Pack datatype

	int nsend, nrecv;			  // Number of send and receive
	MPI_Request *sreqs, *rreqs;	  // Send and recv req
	MPI_Status *sstats, *rstats;  // Send and recv status
	char **sbufs, **rbufs;		  // Send and recv buffer
	int *rsizes;				  // recv size of each message
	MPI_Message rmsg;			  // Receive message

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_INIT)

	// Allocate buffering for write count
	wcnt_local = (int *)NCI_Malloc (sizeof (int) * varp->nchunk * 2);
	wcnt_all   = wcnt_local + varp->nchunk;

	// Allocate buffering for overlaping index
	tstart = (int *)NCI_Malloc (sizeof (int) * varp->ndim * 3);
	tsize  = tstart + varp->ndim;
	tssize = tsize + varp->ndim;
	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * 3);
	osize  = ostart + varp->ndim;

	// Chunk iterator
	citr = osize + varp->ndim;

	// We need to calculate the size of message of each chunk
	// This is just for allocating send buffer
	// We do so by iterating through all request and all chunks they cover
	// If we are not the owner of a chunk, we need to send message
	memset (wcnt_local, 0, sizeof (int) * varp->nchunk);
	nsend	 = 0;
	max_tbuf = 0;
	for (req = 0; req < nreq; req++) {
		// Initialize chunk iterator
		ncchkioi_chunk_itr_init_ex (varp, starts[req], counts[req], citr, &cid, ostart,
									osize);	 // Initialize chunk iterator

		// Iterate through chunks
		do {
			// Calculate overlapping
			overlapsize = varp->esize;
			for (j = 0; j < varp->ndim; j++) { overlapsize *= osize[j]; }

			if (varp->chunk_owner[cid] != ncchkp->rank) {
				// Count number of mnessage we need to send
				if (wcnt_local[cid] == 0) { nsend++; }
				wcnt_local[cid] += overlapsize + sizeof (int) * 2 * varp->ndim;
			} else {
				// We mark covered chunk of our own to prevent unnecessary calculation of overlap
				// -1 is purely a mark, we need to add 1 back to global message count
				wcnt_local[cid] = -1;

				// Record max overlapsize so we know how large the intermediate buffer is needed
				// later
				if (max_tbuf < overlapsize) { max_tbuf = overlapsize; }
			}

		} while (
			ncchkioi_chunk_itr_next_ex (varp, starts[req], counts[req], citr, &cid, ostart, osize));
	}

	// Allocate buffer for sending
	sbufs  = (char **)NCI_Malloc (sizeof (char *) * nsend);
	sreqs  = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * nsend);
	sstats = (MPI_Status *)NCI_Malloc (sizeof (MPI_Status) * nsend);
	j	   = 0;
	// Allocate buffer for data
	for (cid = 0; cid < varp->nchunk; cid++) {
		// Count number of mnessage we need to send
		if (wcnt_local[cid] > 0) {
			// Add space for number of reqs
			sbufs[j++] = (char *)NCI_Malloc (wcnt_local[cid]);
			// We don't need message size anymore, wcnt_local is used to track number of message
			// from now on
			wcnt_local[cid] = 1;
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_INIT)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SYNC)

	// Sync number of messages of each chunk
	CHK_ERR_ALLREDUCE (wcnt_local, wcnt_all, varp->nchunk, MPI_INT, MPI_SUM, ncchkp->comm);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SYNC)

	// Calculate number of recv request
	// This is for all the chunks
	nrecv = 0;
	for (i = 0; i < varp->nmychunk; i++) {
		cid = varp->mychunks[i];
		// We don't need message for our own data
		nrecv += wcnt_all[cid] - wcnt_local[cid];
	}
	rreqs  = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * nrecv);
	rstats = (MPI_Status *)NCI_Malloc (sizeof (MPI_Status) * nrecv);
	rbufs  = (char **)NCI_Malloc (sizeof (char *) * nrecv);
	rsizes = (int *)NCI_Malloc (sizeof (int) * nrecv);

	// Post send and recv
	nrecv = 0;
	nsend = 0;
	for (cid = 0; cid < varp->nchunk; cid++) {
		if (varp->chunk_owner[cid] == ncchkp->rank) {
			NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_RECV_REQ)

			// We are the owner of the chunk
			// Receive data from other process
			for (i = 0; i < wcnt_all[cid] - wcnt_local[cid]; i++) {
				// Get message size, including metadata
				CHK_ERR_MPROBE (MPI_ANY_SOURCE, cid, ncchkp->comm, &rmsg, rstats);
				CHK_ERR_GET_COUNT (rstats, MPI_BYTE, rsizes + nrecv);

				// Allocate buffer
				rbufs[nrecv] = (char *)NCI_Malloc (rsizes[nrecv]);

				// Post irecv
				CHK_ERR_IMRECV (rbufs[nrecv], rsizes[nrecv], MPI_BYTE, &rmsg, rreqs + nrecv);
				nrecv++;
			}

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_RECV_REQ)
		} else {
			// If we any of our request overlap with this chunk, we need to send data
			// We send only 1 message for 1 chunk
			if (wcnt_local[cid] > 0) {
				NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_PACK_REQ)

				packoff = 0;
				// Get chunk iterator
				get_chunk_itr (varp, cid, citr);
				for (req = 0; req < nreq; req++) {
					// Calculate chunk overlap
					overlapsize =
						get_chunk_overlap (varp, citr, starts[req], counts[req], ostart, osize);

					// If current request have any overlap with the chunk, we pack the data and
					// metadata
					if (overlapsize > 0) {
						// Metadata
						tstartp = (int *)(sbufs[nsend] + packoff);
						packoff += varp->ndim * sizeof (int);
						tsizep = (int *)(sbufs[nsend] + packoff);
						packoff += varp->ndim * sizeof (int);
						for (j = 0; j < varp->ndim; j++) {
							tstartp[j] = (int)(ostart[j] - citr[j]);
							tsizep[j]  = (int)osize[j];
						}

						// Pack type
						for (j = 0; j < varp->ndim; j++) {
							tstart[j] = (int)(ostart[j] - starts[req][j]);
							tsize[j]  = (int)counts[req][j];
							tssize[j] = (int)osize[j];
						}
						CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart,
													  MPI_ORDER_C, varp->etype, &ptype);
						CHK_ERR_TYPE_COMMIT (&ptype);

						// Data
						CHK_ERR_PACK (bufs[req], 1, ptype, sbufs[nsend], packoff + overlapsize,
									  &packoff, ncchkp->comm);
						MPI_Type_free (&ptype);
					}
				}

				NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_PACK_REQ)
				NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SEND_REQ)

				// Send the request
				CHK_ERR_ISEND (sbufs[nsend], packoff, MPI_BYTE, varp->chunk_owner[cid], cid,
							   ncchkp->comm, sreqs + nsend);
				nsend++;

				NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SEND_REQ)
			}
		}
	}

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SEND_REQ)

	// Wait for all send
	CHK_ERR_WAITALL (nsend, sreqs, sstats);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SEND_REQ)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_INIT)

	// Preparing chunk cache
	nread = 0;
	for (i = 0; i < varp->nmychunk; i++) {
		cid = varp->mychunks[i];
		if (wcnt_all[cid] && varp->chunk_cache[cid] == NULL) {
			if (varp->chunk_index[cid].len > 0) { nread++; }
		}
	}
	rids  = (int *)NCI_Malloc (sizeof (int) * nread);
	nread = 0;
	for (i = 0; i < varp->nmychunk; i++) {
		cid = varp->mychunks[i];
		if (wcnt_all[cid] || wcnt_local[cid]) {
			if (varp->chunk_cache[cid] == NULL) {
				err = ncchkioi_cache_alloc (ncchkp, varp->chunksize, varp->chunk_cache + cid);
				CHK_ERR
				// varp->chunk_cache[cid] = (NC_chk_cache*)NCI_Malloc(varp->chunksize);
				if (varp->chunk_index[cid].len > 0) { rids[nread++] = cid; }
			} else {
				ncchkioi_cache_visit (ncchkp, varp->chunk_cache[cid]);
			}
		}
	}
	// Increase batch number to indicate allocated chunk buffer can be freed for future allocation
	(ncchkp->cache_serial)++;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_INIT)

	// Read background
	ncchkioi_load_var_bg (ncchkp, varp, nread, rids);

	// Allocate intermediate buffer
	if (max_tbuf > 0) { tbuf = (char *)NCI_Malloc (max_tbuf); }

	// For each chunk we own, we need to receive incoming data
	nrecv = 0;
	for (i = 0; i < varp->nmychunk; i++) {
		cid = varp->mychunks[i];

		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SELF)

		// Handle our own data first if we have any
		if (wcnt_local[cid] < 0) {
			for (req = 0; req < nreq; req++) {
				// Convert chunk id to iterator
				get_chunk_itr (varp, cid, citr);

				// Calculate overlapping region
				overlapsize =
					get_chunk_overlap (varp, citr, starts[req], counts[req], ostart, osize);

				// If anything overlaps
				if (overlapsize > 0) {
					// Pack type from user buffer to (contiguous) intermediate buffer
					for (j = 0; j < varp->ndim; j++) {
						tstart[j] = (int)(ostart[j] - starts[req][j]);
						tsize[j]  = (int)counts[req][j];
						tssize[j] = (int)osize[j];
					}

					CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
												  varp->etype, &ptype);
					CHK_ERR_TYPE_COMMIT (&ptype);

					// Pack data into intermediate buffer
					packoff = 0;
					CHK_ERR_PACK (bufs[req], 1, ptype, tbuf, overlapsize, &packoff, ncchkp->comm);

					MPI_Type_free (&ptype);

					// Pack type from (contiguous) intermediate buffer to chunk buffer
					for (j = 0; j < varp->ndim; j++) {
						tstart[j] = (int)(ostart[j] - citr[j]);
						tsize[j]  = varp->chunkdim[j];
					}
					CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
												  varp->etype, &ptype);
					CHK_ERR_TYPE_COMMIT (&ptype);

					// Unpack data into chunk buffer
					packoff = 0;
					CHK_ERR_UNPACK (tbuf, overlapsize, &packoff, varp->chunk_cache[cid]->buf, 1,
									ptype, ncchkp->comm);

					MPI_Type_free (&ptype);

					// Mark chunk as dirty
					varp->dirty[cid] = 1;
				}
			}
		}

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SELF)
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_RECV_REQ)

		// Now, it is time to process data from other processes

		// Wait for all send requests related to this chunk
		// We remove the impact of -1 mark in wcnt_local[cid]
		CHK_ERR_WAITALL (wcnt_all[cid] - wcnt_local[cid], rreqs + nrecv, rstats + nrecv);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_RECV_REQ)
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_UNPACK_REQ)

		// Process data received
		for (j = nrecv; j < nrecv + wcnt_all[cid] - wcnt_local[cid]; j++) {
			packoff = 0;
			while (packoff < rsizes[j]) {
				// Metadata
				tstartp = (int *)(rbufs[j] + packoff);
				packoff += varp->ndim * sizeof (int);
				tsizep = (int *)(rbufs[j] + packoff);
				packoff += varp->ndim * sizeof (int);

				// Packtype
				CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, varp->chunkdim, tsizep, tstartp,
											  MPI_ORDER_C, varp->etype, &ptype);
				CHK_ERR_TYPE_COMMIT (&ptype);

				// Data
				CHK_ERR_UNPACK (rbufs[j], rsizes[j], &packoff, varp->chunk_cache[cid]->buf, 1,
								ptype, ncchkp->comm);
				MPI_Type_free (&ptype);

				// Mark chunk as dirty
				varp->dirty[cid] = 1;
			}
		}
		nrecv += wcnt_all[cid] - wcnt_local[cid];

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_UNPACK_REQ)
	}

	// Free buffers
	NCI_Free (wcnt_local);

	NCI_Free (tstart);

	NCI_Free (ostart);

	NCI_Free (sreqs);
	NCI_Free (sstats);
	for (i = 0; i < nsend; i++) { NCI_Free (sbufs[i]); }
	NCI_Free (sbufs);

	NCI_Free (rreqs);
	NCI_Free (rstats);
	for (i = 0; i < nrecv; i++) { NCI_Free (rbufs[i]); }
	NCI_Free (rbufs);
	NCI_Free (rsizes);

	if (tbuf != NULL) { NCI_Free (tbuf); }

	if (rids != NULL) { NCI_Free (rids); }

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB)

err_out:;
	return err;
}

int ncchkioi_put_varn_cb_proc (NC_chk *ncchkp,
							   NC_chk_var *varp,
							   int nreq,
							   MPI_Offset *const *starts,
							   MPI_Offset *const *counts,
							   void **bufs) {
	int err;
	int i, j, k;
	int cid, cown;	// Chunk iterator and owner
	int req;

	MPI_Offset *ostart, *osize;
	int *tsize, *tssize, *tstart, *tssizep, *tstartp;  // Size for sub-array type
	MPI_Offset *citr;  // Bounding box for chunks overlapping my own write region

	int *wcnt_local, *wcnt_all;			 // Number of processes that writes to each chunk
	int wrange_local[2], wrange_all[2];	 // Number of processes that writes to each chunk

	int nread;	// Chunks to read for background
	int *rids;

	int overlapsize;	  // Size of overlaping region of request and chunk
	int max_tbuf = 0;	  // Size of intermediate buffer
	char *tbuf	 = NULL;  // Intermediate buffer

	int packoff;		 // Pack offset
	MPI_Datatype ptype;	 // Pack datatype

	int nsend, nrecv;						// Number of send and receive
	MPI_Request *sreq, *rreq;				// Send and recv req
	MPI_Status *sstat, rstat;				// Send and recv status
	char **sbuf, **sbufp, **rbuf, **rbufp;	// Send and recv buffer
	int *rsize, *ssize;						// recv size of each message
	int *sdst;								// recv size of each message
	int *smap;
	MPI_Message rmsg;  // Receive message

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_INIT)

	// Allocate buffering for write count
	wcnt_local = (int *)NCI_Malloc (sizeof (int) * ncchkp->np * 3);
	wcnt_all   = wcnt_local + ncchkp->np;
	smap	   = wcnt_all + ncchkp->np;

	// Allocate buffering for overlaping index
	tstart = (int *)NCI_Malloc (sizeof (int) * varp->ndim * 3);
	tssize = tstart + varp->ndim;
	tsize  = tssize + varp->ndim;
	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * 3);
	osize  = ostart + varp->ndim;

	// Chunk iterator
	citr = osize + varp->ndim;

	// We need to calculate the size of message of each processes
	// This is just for allocating send buffer
	// We do so by iterating through all request and all chunks they cover
	// If we are not the owner of a chunk, we need to send message
	memset (wcnt_local, 0, sizeof (int) * ncchkp->np);
	nsend = 0;

	// Count total number of messages and build a map of accessed chunk to list of comm
	// datastructure
	wrange_local[0] = varp->nchunk;
	wrange_local[1] = 0;
	for (req = 0; req < nreq; req++) {
		ncchkioi_chunk_itr_init (varp, starts[req], counts[req], citr,
								 &cid);	 // Initialize chunk iterator
		do {
			// Chunk owner
			cown = varp->chunk_owner[cid];

			// Mapping to skip list of send requests
			if (wcnt_local[cown] == 0 && cown != ncchkp->rank) { smap[cown] = nsend++; }
			wcnt_local[cown] = 1;  // Need to send message if not owner

			// Record lowest and highest chunk accessed
			if (wrange_local[0] > cid) { wrange_local[0] = cid; }
			if (wrange_local[1] < cid) { wrange_local[1] = cid; }
		} while (ncchkioi_chunk_itr_next (varp, starts[req], counts[req], citr, &cid));
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_INIT)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SYNC)

	// Sync number of messages of each chunk
	CHK_ERR_ALLREDUCE (wcnt_local, wcnt_all, ncchkp->np, MPI_INT, MPI_SUM, ncchkp->comm);
	wrange_local[1] *= -1;
	CHK_ERR_ALLREDUCE (wrange_local, wrange_all, 2, MPI_INT, MPI_MIN, ncchkp->comm);
	nrecv = wcnt_all[ncchkp->rank] -
			wcnt_local[ncchkp->rank];  // We don't need to receive request form self
	wrange_all[1] *= -1;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SYNC)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_PACK_REQ)

	// Allocate data structure for messaging
	sbuf  = (char **)NCI_Malloc (sizeof (char *) * nsend * 2);
	sbufp = sbuf + nsend;
	ssize = (int *)NCI_Malloc (sizeof (int) * nsend * 2);
	sdst  = ssize + nsend;
	sreq  = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * nsend);
	sstat = (MPI_Status *)NCI_Malloc (sizeof (MPI_Status) * nsend);

	rbuf  = (char **)NCI_Malloc (sizeof (char *) * nrecv * 2);
	rbufp = rbuf + nrecv;
	rsize = (int *)NCI_Malloc (sizeof (int) * nrecv);
	rreq  = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * nrecv);

	// Count size of each request
	memset (ssize, 0, sizeof (int) * nsend);
	for (req = 0; req < nreq; req++) {
		ncchkioi_chunk_itr_init_ex (varp, starts[req], counts[req], citr, &cid, ostart,
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
				ssize[j] += overlapsize + sizeof (int) * (varp->ndim * 2 + 1);
			}
		} while (
			ncchkioi_chunk_itr_next_ex (varp, starts[req], counts[req], citr, &cid, ostart, osize));
	}
	// Allocate buffer for send
	for (i = 0; i < nsend; i++) { sbuf[i] = sbufp[i] = (char *)NCI_Malloc (ssize[i]); }

	// Pack requests
	for (req = 0; req < nreq; req++) {
		ncchkioi_chunk_itr_init_ex (varp, starts[req], counts[req], citr, &cid, ostart,
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
				sbufp[j] += varp->ndim * sizeof (int);
				tssizep = (int *)sbufp[j];
				sbufp[j] += varp->ndim * sizeof (int);
				for (i = 0; i < varp->ndim; i++) {
					tstartp[i] = (int)(ostart[i] - citr[i]);
					tssizep[i] = (int)osize[i];
				}

				// Pack type
				for (i = 0; i < varp->ndim; i++) {
					tstart[i] = (int)(ostart[i] - starts[req][i]);
					tsize[i]  = (int)counts[req][i];
				}
				CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssizep, tstart, MPI_ORDER_C,
											  varp->etype, &ptype);
				CHK_ERR_TYPE_COMMIT (&ptype);

				// Data
				packoff = 0;
				CHK_ERR_PACK (bufs[req], 1, ptype, sbufp[j], ssize[j], &packoff, ncchkp->comm);
				sbufp[j] += packoff;
				MPI_Type_free (&ptype);
			}
		} while (
			ncchkioi_chunk_itr_next_ex (varp, starts[req], counts[req], citr, &cid, ostart, osize));
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_PACK_REQ)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SEND_REQ)

	// Post send
	for (i = 0; i < nsend; i++) {
		CHK_ERR_ISEND (sbuf[i], ssize[i], MPI_BYTE, sdst[i], 0, ncchkp->comm, sreq + i);
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SEND_REQ)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_RECV_REQ)

	// Post recv
	for (i = 0; i < nrecv; i++) {
		// Get message size, including metadata
		CHK_ERR_MPROBE (MPI_ANY_SOURCE, 0, ncchkp->comm, &rmsg, &rstat);
		CHK_ERR_GET_COUNT (&rstat, MPI_BYTE, rsize + i);

		// Allocate buffer
		rbuf[i] = rbufp[i] = (char *)NCI_Malloc (rsize[i]);

		// Post irecv
		CHK_ERR_IMRECV (rbuf[i], rsize[i], MPI_BYTE, &rmsg, rreq + i);
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_RECV_REQ)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_INIT)

	// Preparing chunk cache
	for (j = 0; j < varp->nmychunk && varp->mychunks[j] < wrange_all[0]; j++)
		;
	for (k = j; k < varp->nmychunk && varp->mychunks[k] <= wrange_all[1]; k++)
		;
	rids  = (int *)NCI_Malloc (sizeof (int) * (k - j));
	nread = 0;
	for (i = j; i < k; i++) {
		cid = varp->mychunks[i];
		if (varp->chunk_cache[cid] == NULL) {
			err = ncchkioi_cache_alloc (ncchkp, varp->chunksize, varp->chunk_cache + cid);
			CHK_ERR
			// varp->chunk_cache[cid] = (char*)NCI_Malloc(varp->chunksize);
			if (varp->chunk_index[cid].len > 0) { rids[nread++] = cid; }
		} else {
			ncchkioi_cache_visit (ncchkp, varp->chunk_cache[cid]);
		}
	}
	// Increase batch number to indicate allocated chunk buffer can be freed for future allocation
	(ncchkp->cache_serial)++;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_INIT)

	// Read background
	ncchkioi_load_var_bg (ncchkp, varp, nread, rids);

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SELF)

	tbuf = (char *)NCI_Malloc (varp->chunksize);

	// Handle our own data
	for (req = 0; req < nreq; req++) {
		ncchkioi_chunk_itr_init_ex (varp, starts[req], counts[req], citr, &cid, ostart,
									osize);	 // Initialize chunk iterator
		do {
			if (varp->chunk_owner[cid] == ncchkp->rank) {
				// Pack type from user buffer to (contiguous) intermediate buffer
				for (j = 0; j < varp->ndim; j++) {
					tstart[j] = (int)(ostart[j] - starts[req][j]);
					tsize[j]  = (int)counts[req][j];
					tssize[j] = (int)osize[j];
				}
				CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
											  varp->etype, &ptype);
				CHK_ERR_TYPE_COMMIT (&ptype);

				// Pack data into intermediate buffer
				packoff = 0;
				CHK_ERR_PACK (bufs[req], 1, ptype, tbuf, varp->chunksize, &packoff, ncchkp->comm);
				MPI_Type_free (&ptype);
				overlapsize = packoff;

				// Pack type from (contiguous) intermediate buffer to chunk buffer
				for (j = 0; j < varp->ndim; j++) {
					tstart[j] = (int)(ostart[j] - citr[j]);
					tsize[j]  = varp->chunkdim[j];
				}
				CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart, MPI_ORDER_C,
											  varp->etype, &ptype);
				CHK_ERR_TYPE_COMMIT (&ptype);

				// Unpack data into chunk buffer
				packoff = 0;
				CHK_ERR_UNPACK (tbuf, overlapsize, &packoff, varp->chunk_cache[cid]->buf, 1, ptype,
								ncchkp->comm);
				MPI_Type_free (&ptype);

				// Mark chunk as dirty
				varp->dirty[cid] = 1;
			}
		} while (
			ncchkioi_chunk_itr_next_ex (varp, starts[req], counts[req], citr, &cid, ostart, osize));
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SELF)

	// Handle incoming requests
	for (i = 0; i < varp->ndim; i++) { tsize[i] = varp->chunkdim[i]; }
	for (i = 0; i < nrecv; i++) {
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_RECV_REQ)

		// Will wait any provide any benefit?
		MPI_Waitany (nrecv, rreq, &j, &rstat);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_RECV_REQ)
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_UNPACK_REQ)

		while (rbufp[j] < rbuf[j] + rsize[j]) {
			// Metadata
			cid = *(int *)(rbufp[j]);
			rbufp[j] += sizeof (int);
			tstartp = (int *)rbufp[j];
			rbufp[j] += varp->ndim * sizeof (int);
			tssizep = (int *)rbufp[j];
			rbufp[j] += varp->ndim * sizeof (int);

			// Pack type
			CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssizep, tstartp, MPI_ORDER_C,
										  varp->etype, &ptype);
			CHK_ERR_TYPE_COMMIT (&ptype);

			// Data
			packoff = 0;
			CHK_ERR_UNPACK (rbufp[j], rsize[j], &packoff, varp->chunk_cache[cid]->buf, 1, ptype,
							ncchkp->comm);
			rbufp[j] += packoff;
			MPI_Type_free (&ptype);

			// Mark chunk as dirty
			varp->dirty[cid] = 1;
		}
		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_UNPACK_REQ)
	}

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SEND_REQ)

	CHK_ERR_WAITALL (nsend, sreq, sstat);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SEND_REQ)

	// Free buffers
	NCI_Free (wcnt_local);

	NCI_Free (tstart);

	NCI_Free (ostart);

	NCI_Free (sreq);
	NCI_Free (sstat);
	NCI_Free (ssize);
	for (i = 0; i < nsend; i++) { NCI_Free (sbuf[i]); }
	NCI_Free (sbuf);

	NCI_Free (rreq);
	for (i = 0; i < nrecv; i++) { NCI_Free (rbuf[i]); }
	NCI_Free (rbuf);
	NCI_Free (rsize);

	NCI_Free (tbuf);

	NCI_Free (rids);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB)

err_out:;
	return err;
}

int ncchkioi_put_varn (NC_chk *ncchkp,
					   NC_chk_var *varp,
					   int nreq,
					   MPI_Offset *const *starts,
					   MPI_Offset *const *counts,
					   const void *buf) {
	int err = NC_NOERR;
	int i, j;
	MPI_Offset rsize;
	char *bptr = (char *)buf;
	char **bufs;

	if (varp->isrec) {
		for (i = 0; i < nreq; i++) {
			if (ncchkp->recsize < starts[i][0] + counts[i][0]) {
				ncchkp->recsize = starts[i][0] + counts[i][0];
			}
		}
		CHK_ERR_ALLREDUCE (MPI_IN_PLACE, &(ncchkp->recsize), 1, MPI_LONG_LONG, MPI_MAX,
						   ncchkp->comm);  // Sync number of recs
		if (varp->dimsize[0] < ncchkp->recsize) {
			NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_PUT)
			NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_RESIZE)

			ncchkioi_var_resize (ncchkp, varp);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_RESIZE)
			NC_CHK_TIMER_START (NC_CHK_TIMER_PUT)
		}
	}

	// Calculate buffer offset of each request
	bufs = (char **)NCI_Malloc (sizeof (char *) * nreq);
	for (i = 0; i < nreq; i++) {
		bufs[i] = bptr;
		rsize	= varp->esize;
		for (j = 0; j < varp->ndim; j++) { rsize *= counts[i][j]; }
		bptr += rsize;
	}

	// Collective buffer
	switch (ncchkp->comm_unit) {
		case NC_CHK_COMM_CHUNK:
			ncchkioi_put_varn_cb_chunk (ncchkp, varp, nreq, starts, counts, NULL, (void **)bufs);
			break;
		case NC_CHK_COMM_PROC:
			ncchkioi_put_varn_cb_proc (ncchkp, varp, nreq, starts, counts, (void **)bufs);
			break;
	}

	// Write the compressed variable
	ncchkioi_save_var (ncchkp, varp);

err_out:;
	NCI_Free (bufs);

	return err;
}