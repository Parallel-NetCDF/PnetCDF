/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
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

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int ncchkioi_iget_cb_chunk (NC_chk *ncchkp, int nreq, int *reqids, int *stats) {
	int err;
	int i, j;
	int nvar;
	int vid;  // Iterators for variable id
	int *varids;
	int *nreqs;	 // Number of reqids in each variable
	int *nums;	 // Number of reqs in each varn
	int **vreqids;
	int num, maxnum = 0;
	MPI_Offset **starts, **counts, **strides;
	MPI_Offset rsize;
	char **bufs;
	NC_chk_req *req;

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_INIT)

	// Count total number of request in per variable for packed varn request
	nums  = (int *)NCI_Malloc (sizeof (int) * ncchkp->vars.cnt * 2);
	nreqs = nums + ncchkp->vars.cnt;
	memset (nums, 0, sizeof (int) * ncchkp->vars.cnt);
	memset (nreqs, 0, sizeof (int) * ncchkp->vars.cnt);
	for (i = 0; i < nreq; i++) {
		req = ncchkp->getlist.reqs + reqids[i];
		nreqs[req->varid]++;
		nums[req->varid] += req->nreq;
	}

	/* Allocate a skip list of reqids for each vriable
	 * At the same time, we find out the number of starts and counts we need to allocate
	 */
	vreqids	   = (int **)NCI_Malloc (sizeof (int *) * ncchkp->vars.cnt);
	vreqids[0] = (int *)NCI_Malloc (sizeof (int) * nreq);
	maxnum	   = 0;
	i		   = 0;
	nvar	   = 0;
	for (vid = 0; vid < ncchkp->vars.cnt; vid++) {
		if (nreqs[vid] > 0) {
			// Assign buffer to reqid skip list
			vreqids[vid] = vreqids[0] + i;
			i += nreqs[vid];

			// maximum number of starts and counts we need across all variables
			if (maxnum < nums[vid]) { maxnum = nums[vid]; }

			// Number of variable that has request to write
			nvar++;
		}
	}

	varids = (int *)NCI_Malloc (sizeof (int) * nvar);

	// Fill up the skip list
	memset (nreqs, 0, sizeof (int) * ncchkp->vars.cnt);
	for (i = 0; i < nreq; i++) {
		req										 = ncchkp->getlist.reqs + reqids[i];
		vreqids[req->varid][nreqs[req->varid]++] = reqids[i];
	}

	// Allocate parameters
	starts = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *) * maxnum * 2);
	counts = starts + maxnum;
	bufs   = (char **)NCI_Malloc (sizeof (char *) * maxnum);

	/* Pack requests variable by variable
	 */
	nvar = 0;
	for (vid = 0; vid < ncchkp->vars.cnt; vid++) {
		if (nreqs[vid] > 0) {
			// Fill varid in the skip list
			varids[nvar++] = vid;

			// Collect parameters
			num = 0;
			for (j = 0; j < nreqs[vid]; j++) {
				req = ncchkp->getlist.reqs + vreqids[vid][j];

				if (req->nreq > 1) {
					for (i = 0; i < req->nreq; i++) {
						starts[num] = req->starts[i];
						counts[num] = req->counts[i];
						bufs[num++] = req->xbufs[i];
					}
				} else {
					starts[num] = req->start;
					counts[num] = req->count;
					bufs[num++] = req->xbuf;
				}
			}

			NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_GET_CB)
			NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_GET_CB_INIT)

			// Perform collective buffering
			ncchkioi_get_varn_cb_chunk (ncchkp, ncchkp->vars.data + vid, num, starts, counts, NULL,
										(void **)bufs);

			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB)
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_INIT)
		}
	}

	// Free buffers
	NCI_Free (nums);

	NCI_Free (vreqids[0]);
	NCI_Free (vreqids);

	NCI_Free (varids);

	NCI_Free (starts);
	NCI_Free (bufs);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB)
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_INIT)

	return NC_NOERR;
}

int ncchkioi_iget_cb_proc (NC_chk *ncchkp, int nreq, int *reqids, int *stats) {
	int err;
	int i, j, k;
	int cid, cown;	// Chunk iterator
	int vid;
	int r, **reqs;

	MPI_Offset *ostart, *osize;
	int *tsize, *tssize, *tstart, *tssizep, *tstartp;  // Size for sub-array type
	MPI_Offset *citr;  // Bounding box for chunks overlapping my own write region

	int *rcnt_local, *rcnt_all;	 // Number of processes that writes to each proc

	int nread;
	int *rlo_local, *rhi_local;
	int *rlo_all, *rhi_all;
	int *rids;

	int overlapsize;	// Size of overlaping region of request and chunk
	char *tbuf = NULL;	// Intermediate buffer

	int packoff;		 // Pack offset
	MPI_Datatype ptype;	 // Pack datatype
	MPI_Offset poff;	 // Offset of buffer to pack to/ from
	int plen;

	int nsend, nrecv;											  // Number of send and receive
	MPI_Request *sreq, *rreq, *sreq_re, *rreq_re;				  // Send and recv req
	MPI_Status *sstat, rstat, *sstat_re;						  // Send and recv status
	char **sbuf, **sbufp, **rbuf, **rbufp, **sbuf_re, **rbuf_re;  // Send and recv buffer
	int *rsize, *ssize, *rsize_re, *ssize_re;					  // recv size of each message
	int *roff, *roff_re;										  // recv size of each message
	int *sdst;													  // recv size of each message
	int *smap;
	MPI_Message rmsg;  // Receive message
	NC_chk_var *varp;
	NC_chk_req *req;

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_INIT)

	// Allocate buffering for write count
	rcnt_local = (int *)NCI_Malloc (sizeof (int) * ncchkp->np * 3);
	rcnt_all   = rcnt_local + ncchkp->np;
	smap	   = rcnt_all + ncchkp->np;

	// Intermediate buffer for our own data
	tbuf = (char *)NCI_Malloc (ncchkp->max_chunk_size);

	// Allocate buffering for overlaping index
	tsize  = (int *)NCI_Malloc (sizeof (int) * ncchkp->max_ndim * 3);
	tssize = tsize + ncchkp->max_ndim;
	tstart = tssize + ncchkp->max_ndim;
	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * ncchkp->max_ndim * 3);
	osize  = ostart + ncchkp->max_ndim;

	// Chunk iterator
	citr = osize + ncchkp->max_ndim;

	// Access range
	rlo_local = (int *)NCI_Malloc (sizeof (int) * ncchkp->vars.cnt * 5);
	rhi_local = rlo_local + ncchkp->vars.cnt;
	rlo_all	  = rhi_local + ncchkp->vars.cnt;
	rhi_all	  = rlo_all + ncchkp->vars.cnt;
	rids	  = rhi_all + ncchkp->vars.cnt;

	for (i = 0; i < ncchkp->vars.cnt; i++) {
		rlo_local[i] = 2147483647;
		rhi_local[i] = -1;
	}

	// We need to calculate the size of message of each chunk
	// This is just for allocating send buffer
	// We do so by iterating through all request and all chunks they cover
	// If we are not the owner of a chunk, we need to send message
	memset (rcnt_local, 0, sizeof (int) * ncchkp->np);
	nsend = 0;

	// count total number of messages and build a map of accessed chunk to list of comm
	// datastructure
	for (i = 0; i < nreq; i++) {
		req	 = ncchkp->getlist.reqs + reqids[i];
		varp = ncchkp->vars.data + req->varid;
		for (r = 0; r < req->nreq; r++) {
			ncchkioi_chunk_itr_init (varp, req->starts[r], req->counts[r], citr,
									 &cid);	 // Initialize chunk iterator
			do {
				// Chunk owner
				cown = varp->chunk_owner[cid];

				// Mapping to skip list of send requests
				if (rcnt_local[cown] == 0 && cown != ncchkp->rank) { smap[cown] = nsend++; }
				rcnt_local[cown] = 1;  // Need to send message if not owner

				if (rlo_local[req->varid] > cid) { rlo_local[req->varid] = cid; }
				if (rhi_local[req->varid] < cid) { rhi_local[req->varid] = cid; }
			} while (ncchkioi_chunk_itr_next (varp, req->starts[r], req->counts[r], citr, &cid));
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_INIT)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SYNC)

	// Sync number of messages of each chunk
	CHK_ERR_ALLREDUCE (rcnt_local, rcnt_all, ncchkp->np, MPI_INT, MPI_SUM, ncchkp->comm);
	nrecv = rcnt_all[ncchkp->rank] -
			rcnt_local[ncchkp->rank];  // We don't need to receive request form self

#ifdef PNETCDF_PROFILING
	ncchkp->nsend += nrecv + nsend;
	ncchkp->nrecv += nrecv + nsend;
#endif

	for (i = 0; i < ncchkp->vars.cnt; i++) { rhi_local[i] *= -1; }
	CHK_ERR_ALLREDUCE (rlo_local, rlo_all, ncchkp->vars.cnt * 2, MPI_INT, MPI_MIN, ncchkp->comm);
	for (i = 0; i < ncchkp->vars.cnt; i++) { rhi_all[i] *= -1; }

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SYNC)
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_PACK_REQ)

	// Allocate data structure for messaging
	sbuf  = (char **)NCI_Malloc (sizeof (char *) * (nsend * 2 + nrecv));
	sbufp = sbuf + (nsend + nrecv);
	ssize = (int *)NCI_Malloc (sizeof (int) * (nsend * 2 + nrecv * 1));
	sdst  = ssize + (nsend + nrecv);
	sreq  = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * (nsend + nrecv));
	sstat = (MPI_Status *)NCI_Malloc (sizeof (MPI_Status) * (nsend + nrecv));
	reqs  = (int **)NCI_Malloc (sizeof (int *) * nsend);

	rbuf  = (char **)NCI_Malloc (sizeof (char *) * (nsend + nrecv * 2));
	rbufp = rbuf + (nsend + nrecv);
	rsize = (int *)NCI_Malloc (sizeof (int) * (nsend + nrecv));
	rreq  = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * (nsend + nrecv));

	sbuf_re	 = sbuf + nsend;
	ssize_re = ssize + nsend;
	sreq_re	 = sreq + nsend;
	sstat_re = sstat + nsend;

	rbuf_re	 = rbuf + nrecv;
	rsize_re = rsize + nrecv;
	rreq_re	 = rreq + nrecv;

	// req->counts[r] size of each request
	memset (ssize, 0, sizeof (int) * nsend);
	memset (rsize_re, 0, sizeof (int) * nsend);
	memset (rcnt_local, 0, sizeof (int) * nsend);
	for (i = 0; i < nreq; i++) {
		req	 = ncchkp->getlist.reqs + reqids[i];
		varp = ncchkp->vars.data + req->varid;
		for (r = 0; r < req->nreq; r++) {
			ncchkioi_chunk_itr_init_ex (varp, req->starts[r], req->counts[r], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				// Chunk owner
				cown = varp->chunk_owner[cid];
				if (cown != ncchkp->rank) {
					j		= smap[cown];
					sdst[j] = cown;	 // Record a reverse map by the way

					// Count overlap
					overlapsize = varp->esize;
					for (k = 0; k < varp->ndim; k++) { overlapsize *= osize[k]; }
					ssize[j] += sizeof (int) * (varp->ndim * 2 + 2);
					rsize_re[j] += overlapsize;
					rcnt_local[j]++;
				}
			} while (ncchkioi_chunk_itr_next_ex (varp, req->starts[r], req->counts[r], citr, &cid,
												 ostart, osize));
		}
	}

	// Allocate buffer for send
	for (i = 0; i < nsend; i++) {
		ssize[i] += sizeof (int);
#ifdef PNETCDF_DEBUG
		assert (ssize[i] >= 0);
#endif
		sbuf[i] = sbufp[i] = (char *)NCI_Malloc (ssize[i]);
		*((int *)sbufp[i]) = rsize_re[i];
		sbufp[i] += sizeof (int);
		rbuf_re[i] = (char *)NCI_Malloc (rsize_re[i]);
		reqs[i]	   = (int *)NCI_Malloc (sizeof (int) * rcnt_local[i] * 2);
#ifdef PNETCDF_PROFILING
		ncchkp->sendsize += ssize[i];
		ncchkp->recvsize += rsize_re[i];
#endif
	}

	// Pack requests
	memset (rcnt_local, 0, sizeof (int) * nsend);
	for (i = 0; i < nreq; i++) {
		req	 = ncchkp->getlist.reqs + reqids[i];
		varp = ncchkp->vars.data + req->varid;
		for (r = 0; r < req->nreq; r++) {
			ncchkioi_chunk_itr_init_ex (varp, req->starts[r], req->counts[r], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				// Chunk owner
				cown = varp->chunk_owner[cid];
				if (cown != ncchkp->rank) {
					j = smap[cown];

					// Pack metadata
					*((int *)sbufp[j]) = varp->varid;
					sbufp[j] += sizeof (int);
					*((int *)sbufp[j]) = cid;
					sbufp[j] += sizeof (int);
					tstartp = (int *)sbufp[j];
					sbufp[j] += sizeof (int) * varp->ndim;
					tssizep = (int *)sbufp[j];
					sbufp[j] += sizeof (int) * varp->ndim;
					for (k = 0; k < varp->ndim; k++) {
						tstartp[k] = (int)(ostart[k] - citr[k]);
						tssizep[k] = (int)osize[k];
					}

					// Record source of the request
					reqs[j][rcnt_local[j]++] = i;
					reqs[j][rcnt_local[j]++] = r;

#ifdef PNETCDF_PROFILING
					ncchkp->nremote++;
#endif
				}
			} while (ncchkioi_chunk_itr_next_ex (varp, req->starts[r], req->counts[r], citr, &cid,
												 ostart, osize));
		}
	}

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
#ifdef PNETCDF_PROFILING
		ncchkp->recvsize += rsize[i];
#endif
		// Post irecv
		CHK_ERR_IMRECV (rbuf[i], rsize[i], MPI_BYTE, &rmsg, rreq + i);
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REQ)

	nread = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		if (rhi_all[i] >= rlo_all[i]) {
			varp		= ncchkp->vars.data + i;
			rids[nread] = i;
			for (j = 0; j < varp->nmychunk && varp->mychunks[j] < rlo_all[i]; j++)
				;
			for (k = j; k < varp->nmychunk && varp->mychunks[k] <= rhi_all[i]; k++)
				;
			rlo_all[nread]	 = j;
			rhi_all[nread++] = k;
		}
	}

	NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_GET_CB)
#ifdef PNETCDF_PROFILING
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_BARR)
	MPI_Barrier (ncchkp->comm);
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_BARR)
#endif
	err = ncchkioi_load_nvar (ncchkp, nread, rids, rlo_all, rhi_all);
	CHK_ERR 
	// Increase batch number to indicate allocated chunk buffer can be freed for future allocation
	(ncchkp->cache_serial)++;
	
	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB)

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SELF)

	// Handle our own data
	for (i = 0; i < nreq; i++) {
		req	 = ncchkp->getlist.reqs + reqids[i];
		varp = ncchkp->vars.data + req->varid;
		for (r = 0; r < req->nreq; r++) {
			ncchkioi_chunk_itr_init_ex (varp, req->starts[r], req->counts[r], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				if (varp->chunk_owner[cid] == ncchkp->rank) {
					// Pack type from chunk cache to (contiguous) intermediate buffer
					for (j = 0; j < varp->ndim; j++) {
						tstart[j] = (int)(ostart[j] - citr[j]);
						tsize[j]  = varp->chunkdim[j];
						tssize[j] = (int)osize[j];
					}
					err =
						ncchkioi_subarray_off_len (varp->ndim, tsize, tssize, tstart, &poff, &plen);
					if (err == 0) {
						plen *= varp->esize;
						poff *= varp->esize;
						memcpy (tbuf, varp->chunk_cache[cid]->buf + poff, plen);
						overlapsize = plen;
					} else {
						CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart,
													  MPI_ORDER_C, varp->etype, &ptype);
						CHK_ERR_TYPE_COMMIT (&ptype);

						// Pack data into intermediate buffer
						packoff = 0;
						CHK_ERR_PACK (varp->chunk_cache[cid]->buf, 1, ptype, tbuf, varp->chunksize,
									  &packoff, ncchkp->comm);
						MPI_Type_free (&ptype);
						overlapsize = packoff;
					}

					// Pack type from (contiguous) intermediate buffer to chunk buffer
					for (j = 0; j < varp->ndim; j++) {
						tstart[j] = (int)(ostart[j] - req->starts[r][j]);
						tsize[j]  = (int)req->counts[r][j];
					}
					err =
						ncchkioi_subarray_off_len (varp->ndim, tsize, tssize, tstart, &poff, &plen);
					if (err == 0) {
						plen *= varp->esize;
						poff *= varp->esize;
						memcpy (req->xbufs[r] + poff, tbuf, plen);
					} else {
						CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart,
													  MPI_ORDER_C, varp->etype, &ptype);
						CHK_ERR_TYPE_COMMIT (&ptype);

						// Unpack data into chunk buffer
						packoff = 0;
						CHK_ERR_UNPACK (tbuf, overlapsize, &packoff, req->xbufs[r], 1, ptype,
										ncchkp->comm);
						MPI_Type_free (&ptype);
					}
#ifdef PNETCDF_PROFILING
					ncchkp->nlocal++;
#endif
				}
			} while (ncchkioi_chunk_itr_next_ex (varp, req->starts[r], req->counts[r], citr, &cid,
												 ostart, osize));
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SELF)

	// Handle incoming requests
	for (i = 0; i < nrecv; i++) {
		NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_RECV_REQ)

		// Will wait any provide any benefit?
		MPI_Waitany (nrecv, rreq, &j, &rstat);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_RECV_REQ)

		packoff		= 0;
		ssize_re[j] = *((int *)rbufp[j]);
		rbufp[j] += sizeof (int);
#ifdef PNETCDF_DEBUG
		assert (ssize_re[j] >= 0);
#endif
		sbuf_re[j] = (char *)NCI_Malloc (ssize_re[j]);
		while (rbufp[j] < rbuf[j] + rsize[j]) {
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_UNPACK_REQ)

			// Retrieve metadata
			vid = *((int *)rbufp[j]);
			rbufp[j] += sizeof (int);
			cid = *((int *)rbufp[j]);
			rbufp[j] += sizeof (int);
			varp	= ncchkp->vars.data + vid;
			tstartp = (int *)rbufp[j];
			rbufp[j] += sizeof (int) * varp->ndim;
			tssizep = (int *)rbufp[j];
			rbufp[j] += sizeof (int) * varp->ndim;

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_UNPACK_REQ)
			NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_PACK_REP)

			err = ncchkioi_subarray_off_len (varp->ndim, varp->chunkdim, tssizep, tstartp, &poff,
											 &plen);
			if (err == 0) {
				plen *= varp->esize;
				poff *= varp->esize;
				memcpy (sbuf_re[j] + packoff, varp->chunk_cache[cid]->buf + poff, plen);
				packoff += plen;
			} else {
				// Pack type
				CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, varp->chunkdim, tssizep, tstartp,
											  MPI_ORDER_C, varp->etype, &ptype);
				CHK_ERR_TYPE_COMMIT (&ptype);

				// Pack data
				CHK_ERR_PACK (varp->chunk_cache[cid]->buf, 1, ptype, sbuf_re[j], ssize_re[j],
							  &packoff, ncchkp->comm);
				MPI_Type_free (&ptype);
			}

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_PACK_REP)
#ifdef PNETCDF_PROFILING
			ncchkp->nreq++;
#endif
		}

		NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REQ)

		// Send Response
		CHK_ERR_ISEND (sbuf_re[j], packoff, MPI_BYTE, rstat.MPI_SOURCE, 1, ncchkp->comm,
					   sreq_re + j);
#ifdef PNETCDF_PROFILING
		ncchkp->sendsize += packoff;
#endif
		NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_GET_CB_SEND_REQ)
	}

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REQ)

	// Wait for all request
	CHK_ERR_WAITALL (nsend, sreq, sstat);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REQ)

	// Handle reply
	memset (rcnt_local, 0, sizeof (int) * nsend);
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
			vid = *((int *)sbufp[j]);
			sbufp[j] += sizeof (int);
			cid = *((int *)sbufp[j]);
			sbufp[j] += sizeof (int);
			varp	= ncchkp->vars.data + vid;
			tstartp = (int *)sbufp[j];
			sbufp[j] += sizeof (int) * varp->ndim;
			tssizep = (int *)sbufp[j];
			sbufp[j] += sizeof (int) * varp->ndim;

			k	= reqs[j][rcnt_local[j]++];
			r	= reqs[j][rcnt_local[j]++];
			req = ncchkp->getlist.reqs + reqids[k];
			get_chunk_itr (varp, cid, citr);
			for (k = 0; k < varp->ndim; k++) {
				tstartp[k] += (int)(citr[k] - req->starts[r][k]);
				tsize[k] = req->counts[r][k];
			}

			err = ncchkioi_subarray_off_len (varp->ndim, tsize, tssizep, tstartp, &poff, &plen);
			if (err == 0) {
				plen *= varp->esize;
				poff *= varp->esize;
				memcpy (req->xbufs[r] + poff, rbuf_re[j] + packoff, plen);
				packoff += plen;
			} else {
				// Pack type
				CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssizep, tstartp, MPI_ORDER_C,
											  varp->etype, &ptype);
				CHK_ERR_TYPE_COMMIT (&ptype);

				// Pack data
				CHK_ERR_UNPACK (rbuf_re[j], rsize_re[j], &packoff, req->xbufs[r], 1, ptype,
								ncchkp->comm);
				MPI_Type_free (&ptype);
			}
		}

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_UNPACK_REP)
	}

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CB_SEND_REP)

	// Wait for all Response
	CHK_ERR_WAITALL (nrecv, sreq_re, sstat_re);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB_SEND_REP)

	// Free buffers
	NCI_Free (rcnt_local);

	NCI_Free (tsize);

	NCI_Free (ostart);

	NCI_Free (tbuf);

	NCI_Free (sreq);
	NCI_Free (sstat);
	NCI_Free (ssize);
	for (i = 0; i < nsend; i++) { NCI_Free (reqs[i]); }
	for (i = 0; i < nsend + nrecv; i++) {
		NCI_Free (sbuf[i]);
		NCI_Free (rbuf[i]);
	}
	NCI_Free (sbuf);
	NCI_Free (reqs);

	NCI_Free (rreq);
	NCI_Free (rbuf);
	NCI_Free (rsize);

	NCI_Free (rlo_local);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CB)

err_out:;
	return err;
}
