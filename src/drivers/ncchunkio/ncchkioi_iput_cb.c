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

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int ncchkioi_iput_cb_chunk (NC_chk *ncchkp, int nreq, int *reqids, int *stats) {
	int err;
	int i, j;
	int vid;	 // Iterators for variable id
	int *nreqs;	 // Number of reqids in each variable
	int *nums;	 // Number of reqs in each varn
	int **vreqids;
	int num, maxnum = 0;
	MPI_Offset **starts, **counts, **strides;
	MPI_Offset rsize;
	char **bufs;
	NC_chk_req *req;

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_INIT)

	// Count total number of request in per variable for packed varn request
	nums  = (int *)NCI_Malloc (sizeof (int) * ncchkp->vars.cnt * 2);
	nreqs = nums + ncchkp->vars.cnt;
	memset (nums, 0, sizeof (int) * ncchkp->vars.cnt);
	memset (nreqs, 0, sizeof (int) * ncchkp->vars.cnt);
	for (i = 0; i < nreq; i++) {
		req = ncchkp->putlist.reqs + reqids[i];
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
	for (vid = 0; vid < ncchkp->vars.cnt; vid++) {
		if (nreqs[vid] > 0) {
			// Assign buffer to reqid skip list
			vreqids[vid] = vreqids[0] + i;
			i += nreqs[vid];

			// maximum number of starts and counts we need across all variables
			if (maxnum < nums[vid]) { maxnum = nums[vid]; }
		}
	}

	// Fill up the skip list
	memset (nreqs, 0, sizeof (int) * ncchkp->vars.cnt);
	for (i = 0; i < nreq; i++) {
		req										 = ncchkp->putlist.reqs + reqids[i];
		vreqids[req->varid][nreqs[req->varid]++] = reqids[i];
	}

	// Allocate parameters
	starts = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *) * maxnum * 2);
	counts = starts + maxnum;
	bufs   = (char **)NCI_Malloc (sizeof (char *) * maxnum);

	/* Pack requests variable by variable
	 */
	for (vid = 0; vid < ncchkp->vars.cnt; vid++) {
		if (nreqs[vid] > 0) {
			// Collect parameters
			num = 0;
			for (j = 0; j < nreqs[vid]; j++) {
				req = ncchkp->putlist.reqs + vreqids[vid][j];

				for (i = 0; i < req->nreq; i++) {
					starts[num] = req->starts[i];
					counts[num] = req->counts[i];
					bufs[num++] = req->xbufs[i];
				}
			}

			NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_PUT_CB)
			NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_PUT_CB_INIT)

			// Perform collective buffering
			ncchkioi_put_varn_cb_chunk (ncchkp, ncchkp->vars.data + vid, num, starts, counts, NULL,
										(void **)bufs);

			NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB)
			NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_INIT)
		}
	}

	// Free buffers
	NCI_Free (nums);

	NCI_Free (vreqids[0]);
	NCI_Free (vreqids);

	NCI_Free (starts);
	NCI_Free (bufs);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB)
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_INIT)

	return NC_NOERR;
}

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int ncchkioi_iput_cb_proc (NC_chk *ncchkp, int nreq, int *reqids, int *stats) {
	int err;
	int i, j, k;
	int cid, cown;	// Chunk iterator and owner
	int vid;
	int r;
	MPI_Offset *ostart = NULL, *osize;
	int *tsize, *tssize, *tstart = NULL, *tssizep, *tstartp;  // Size for sub-array type
	MPI_Offset *citr;  // Bounding box for chunks overlapping my own write region

	int *wcnt_local = NULL, *wcnt_all;	// Number of processes that writes to each chunk

	int nread;
	int *rlo_local, *rhi_local;
	int *rlo_all, *rhi_all;
	int *rids;

	int overlapsize;	// Size of overlaping region of request and chunk
	char *tbuf = NULL;	// Intermediate buffer

	int packoff;		 // Pack offset
	MPI_Offset pboff;	 // Offset of buffer to pack to/ from
	MPI_Datatype ptype;	 // Pack datatype
	int plen;

	int nsend, nrecv;						 // Number of send and receive
	MPI_Request *sreq = NULL, *rreq = NULL;	 // Send and recv req
	MPI_Status *sstat = NULL, rstat;		 // Send and recv status
	char **sbuf = NULL, **rbuf = NULL;		 // Send and recv buffer
	char **sbufp, **rbufp;					 // Send and recv buffer pointer
	int *rsize = NULL, *ssize = NULL;		 // Send and recv size of each message
	MPI_Offset totalsize;
	int *sdst;	// recv size of each message
	int *smap;
	MPI_Message rmsg;  // Receive message
	NC_chk_var *varp;
	NC_chk_req *req;

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_INIT)

	// Allocate buffering for write count
	wcnt_local = (int *)NCI_Malloc (sizeof (int) * ncchkp->np * 3);
	CHK_PTR (wcnt_local)
	wcnt_all = wcnt_local + ncchkp->np;
	smap	 = wcnt_all + ncchkp->np;

	// Intermediate buffer for our own data
	tbuf = (char *)NCI_Malloc (ncchkp->max_chunk_size);
	CHK_PTR (tbuf)

	// Allocate buffering for overlaping index
	tstart = (int *)NCI_Malloc (sizeof (int) * ncchkp->max_ndim * 3);
	CHK_PTR (tstart)
	tssize = tstart + ncchkp->max_ndim;
	tsize  = tssize + ncchkp->max_ndim;
	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * ncchkp->max_ndim * 3);
	CHK_PTR (ostart)
	osize = ostart + ncchkp->max_ndim;

	// Chunk iterator
	citr = osize + ncchkp->max_ndim;

	// Access range
	rlo_local = (int *)NCI_Malloc (sizeof (int) * ncchkp->vars.cnt * 5);
	CHK_PTR (rlo_local)
	rhi_local = rlo_local + ncchkp->vars.cnt;
	rlo_all	  = rhi_local + ncchkp->vars.cnt;
	rhi_all	  = rlo_all + ncchkp->vars.cnt;
	rids	  = rhi_all + ncchkp->vars.cnt;

	for (i = 0; i < ncchkp->vars.cnt; i++) {
		rlo_local[i] = 2147483647;
		rhi_local[i] = -1;
	}

	// We need to calculate the size of message of each processes
	// This is just for allocating send buffer
	// We do so by iterating through all request and all chunks they cover
	// If we are not the owner of a chunk, we need to send message
	memset (wcnt_local, 0, sizeof (int) * ncchkp->np);
	nsend = 0;

	// Count total number of messages and build a map of accessed chunk to list of comm
	// datastructure
	for (i = 0; i < nreq; i++) {
		req	 = ncchkp->putlist.reqs + reqids[i];
		varp = ncchkp->vars.data + req->varid;
		for (r = 0; r < req->nreq; r++) {
			ncchkioi_chunk_itr_init (varp, req->starts[r], req->counts[r], citr,
									 &cid);	 // Initialize chunk iterator
			do {
				// Chunk owner
				cown = varp->chunk_owner[cid];

				// Mapping to skip list of send requests
				if (wcnt_local[cown] == 0 && cown != ncchkp->rank) { smap[cown] = nsend++; }
				wcnt_local[cown] = 1;  // Need to send message if not owner

				if (rlo_local[req->varid] > cid) { rlo_local[req->varid] = cid; }
				if (rhi_local[req->varid] < cid) { rhi_local[req->varid] = cid; }
			} while (ncchkioi_chunk_itr_next (varp, req->starts[r], req->counts[r], citr, &cid));
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_INIT)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SYNC)

	// Sync number of messages of each chunk
	CHK_ERR_ALLREDUCE (wcnt_local, wcnt_all, ncchkp->np, MPI_INT, MPI_SUM, ncchkp->comm);
	nrecv = wcnt_all[ncchkp->rank] -
			wcnt_local[ncchkp->rank];  // We don't need to receive request form self

	for (i = 0; i < ncchkp->vars.cnt; i++) { rhi_local[i] *= -1; }
	CHK_ERR_ALLREDUCE (rlo_local, rlo_all, ncchkp->vars.cnt * 2, MPI_INT, MPI_MIN, ncchkp->comm);
	for (i = 0; i < ncchkp->vars.cnt; i++) { rhi_all[i] *= -1; }

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SYNC)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_PACK_REQ)

	// Allocate data structure for messaging
	sbuf = (char **)NCI_Malloc (sizeof (char *) * nsend * 2);
	CHK_PTR (sbuf)
	sbufp = sbuf + nsend;
	ssize = (int *)NCI_Malloc (sizeof (int) * nsend * 2);
	CHK_PTR (ssize)
	sdst = ssize + nsend;
	sreq = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * nsend);
	CHK_PTR (sreq)
	sstat = (MPI_Status *)NCI_Malloc (sizeof (MPI_Status) * nsend);
	CHK_PTR (sstat)

	rbuf = (char **)NCI_Malloc (sizeof (char *) * nrecv * 2);
	CHK_PTR (rbuf)
	rbufp = rbuf + nrecv;
	rsize = (int *)NCI_Malloc (sizeof (int) * nrecv);
	CHK_PTR (rsize)
	rreq = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * nrecv);
	CHK_PTR (rreq)

	// Count size of each request
	memset (ssize, 0, sizeof (int) * nsend);
	for (k = 0; k < nreq; k++) {
		req	 = ncchkp->putlist.reqs + reqids[k];
		varp = ncchkp->vars.data + req->varid;

		for (r = 0; r < req->nreq; r++) {
			ncchkioi_chunk_itr_init_ex (varp, req->starts[r], req->counts[r], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				// Chunk index and owner
				cown = varp->chunk_owner[cid];
				if (cown != ncchkp->rank) {
					j		= smap[cown];
					sdst[j] = cown;	 // Record a reverse map by the way

					// Count overlap
					overlapsize = varp->esize;
					for (i = 0; i < varp->ndim; i++) { overlapsize *= osize[i]; }
					ssize[j] += overlapsize + sizeof (int) * (varp->ndim * 2 + 2);
#ifdef PNETCDF_DEBUG
					if (ssize[j] < 0) { RET_ERR (NC_EAINT_TOO_SMALL) }
#endif
				}
			} while (ncchkioi_chunk_itr_next_ex (varp, req->starts[r], req->counts[r], citr, &cid,
												 ostart, osize));
		}
	}

	// Allocate buffer for send
	totalsize = 0;
	for (i = 0; i < nsend; i++) {
#ifdef PNETCDF_DEBUG
		assert (ssize[i] >= 0);
#endif
		totalsize += ssize[i];
	}
	if (nsend > 0) {
		sbuf[0] = sbufp[0] = (char *)NCI_Malloc (totalsize);
		CHK_PTR (sbuf[0])
		for (i = 1; i < nsend; i++) { sbuf[i] = sbufp[i] = sbuf[i - 1] + ssize[i - 1]; }
	}
#ifdef PNETCDF_PROFILING
	ncchkp->nsend += nsend;
	ncchkp->sendsize += totalsize;
#endif

	// Pack requests
	for (k = 0; k < nreq; k++) {
		req	 = ncchkp->putlist.reqs + reqids[k];
		varp = ncchkp->vars.data + req->varid;

		for (r = 0; r < req->nreq; r++) {
			ncchkioi_chunk_itr_init_ex (varp, req->starts[r], req->counts[r], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				// Chunk index and owner
				cown = varp->chunk_owner[cid];

				if (cown != ncchkp->rank) {
					j = smap[cown];

					// Pack metadata
					*((int *)(sbufp[j])) = req->varid;
					sbufp[j] += sizeof (int);
					*((int *)(sbufp[j])) = cid;
					sbufp[j] += sizeof (int);
					tstartp = (int *)sbufp[j];
					sbufp[j] += varp->ndim * sizeof (int);
					tssizep = (int *)sbufp[j];
					sbufp[j] += varp->ndim * sizeof (int);

					for (i = 0; i < varp->ndim; i++) {
						tstartp[i] = (int)(ostart[i] - citr[i]);
						tssizep[i] = (int)osize[i];
					}

					// Pack type from user buffer to send buffer
					for (i = 0; i < varp->ndim; i++) {
						tsize[i]  = (int)req->counts[r][i];
						tstart[i] = (int)(ostart[i] - req->starts[r][i]);
					}
					err = ncchkioi_subarray_off_len (varp->ndim, tsize, tssizep, tstart, &pboff,
													 &plen);
					if (err == 0) {
						plen *= varp->esize;
						pboff *= varp->esize;
#ifdef PNETCDF_DEBUG
						if (sbufp[j] - sbuf[j] + plen > ssize[j]) { RET_ERR (NC_EINTERNAL) }
#endif
						memcpy (sbufp[j], req->xbufs[r] + pboff, plen);
						sbufp[j] += plen;
					} else {
						CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssizep, tstart,
													  MPI_ORDER_C, varp->etype, &ptype);
						CHK_ERR_TYPE_COMMIT (&ptype);

						// Pack data
						packoff = 0;
						CHK_ERR_PACK (req->xbufs[r], 1, ptype, sbufp[j], ssize[j], &packoff,
									  ncchkp->comm);
						sbufp[j] += packoff;
						MPI_Type_free (&ptype);
					}

#ifdef PNETCDF_PROFILING
					ncchkp->nremote++;
#endif
				}
			} while (ncchkioi_chunk_itr_next_ex (varp, req->starts[r], req->counts[r], citr, &cid,
												 ostart, osize));
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_PACK_REQ)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SEND_REQ)

	// Post send
	for (i = 0; i < nsend; i++) {
		CHK_ERR_ISEND (sbuf[i], ssize[i], MPI_BYTE, sdst[i], 0, ncchkp->comm, sreq + i);
		if (ssize[i] == 24 && sdst[i] == 983) {
			printf ("rank %d: wrong ssize to 983\n", ncchkp->rank);
			abort ();
		}
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
		CHK_PTR (rbuf[i])

#ifdef PNETCDF_PROFILING
		ncchkp->recvsize += rsize[i];
#endif

		// Post irecv
		CHK_ERR_IMRECV (rbuf[i], rsize[i], MPI_BYTE, &rmsg, rreq + i);
	}
#ifdef PNETCDF_PROFILING
	ncchkp->nrecv += nrecv;
#endif

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_RECV_REQ)

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

#ifdef PNETCDF_PROFILING
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_BARR)
	MPI_Barrier (ncchkp->comm);
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_BARR)
#endif

	NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_PUT_CB)
	err = ncchkioi_load_nvar_bg (ncchkp, nread, rids, rlo_all, rhi_all);
	CHK_ERR
		// Increase batch number to indicate allocated chunk buffer can be freed for future
		// allocation
		(ncchkp->cache_serial)
	++;

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SELF)

	// Handle our own data
	for (k = 0; k < nreq; k++) {
		req	 = ncchkp->putlist.reqs + reqids[k];
		varp = ncchkp->vars.data + req->varid;

		for (r = 0; r < req->nreq; r++) {
			ncchkioi_chunk_itr_init_ex (varp, req->starts[r], req->counts[r], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				// Chunk index and owner
				if (varp->chunk_owner[cid] == ncchkp->rank) {
					// Pack type from user buffer to (contiguous) intermediate buffer
					for (j = 0; j < varp->ndim; j++) {
						tstart[j] = (int)(ostart[j] - req->starts[r][j]);
						tsize[j]  = (int)req->counts[r][j];
						tssize[j] = (int)osize[j];
					}
					err = ncchkioi_subarray_off_len (varp->ndim, tsize, tssize, tstart, &pboff,
													 &plen);
					if (err == 0) {
						plen *= varp->esize;
						pboff *= varp->esize;
						memcpy (tbuf, req->xbufs[r] + pboff, plen);
						overlapsize = plen;
					} else {
						CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart,
													  MPI_ORDER_C, varp->etype, &ptype);
						CHK_ERR_TYPE_COMMIT (&ptype);

						// Pack data into intermediate buffer
						packoff = 0;
						CHK_ERR_PACK (req->xbufs[r], 1, ptype, tbuf, varp->chunksize, &packoff,
									  ncchkp->comm);
						MPI_Type_free (&ptype);
						overlapsize = packoff;
					}

					// Pack type from (contiguous) intermediate buffer to chunk buffer
					for (j = 0; j < varp->ndim; j++) {
						tstart[j] = (int)(ostart[j] - citr[j]);
						tsize[j]  = varp->chunkdim[j];
					}
					err = ncchkioi_subarray_off_len (varp->ndim, tsize, tssize, tstart, &pboff,
													 &plen);
					if (err == 0) {
						plen *= varp->esize;
						pboff *= varp->esize;
						memcpy (varp->chunk_cache[cid]->buf + pboff, tbuf, plen);
					} else {
						CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, tsize, tssize, tstart,
													  MPI_ORDER_C, varp->etype, &ptype);
						CHK_ERR_TYPE_COMMIT (&ptype);

						// Unpack data into chunk buffer
						packoff = 0;
						CHK_ERR_UNPACK (tbuf, overlapsize, &packoff, varp->chunk_cache[cid]->buf, 1,
										ptype, ncchkp->comm);
						MPI_Type_free (&ptype);
					}

					// Mark chunk as dirty
					varp->dirty[cid] = 1;
#ifdef PNETCDF_PROFILING
					ncchkp->nlocal++;
#endif
				}
			} while (ncchkioi_chunk_itr_next_ex (varp, req->starts[r], req->counts[r], citr, &cid,
												 ostart, osize));
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SELF)

	// Handle incoming requests
	for (i = 0; i < nrecv; i++) {
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_RECV_REQ)

		// Will wait any provide any benefit?
		MPI_Waitany (nrecv, rreq, &j, &rstat);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_RECV_REQ)
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_UNPACK_REQ)

		packoff = 0;
		while (rbufp[j] - rbuf[j] < rsize[j]) {
			// Retrieve metadata
			vid = *((int *)(rbufp[j]));
			rbufp[j] += sizeof (int);
			cid = *((int *)(rbufp[j]));
			rbufp[j] += sizeof (int);
			varp	= ncchkp->vars.data + vid;
			tstartp = (int *)rbufp[j];
			rbufp[j] += varp->ndim * sizeof (int);
			tssizep = (int *)rbufp[j];
			rbufp[j] += varp->ndim * sizeof (int);

			err = ncchkioi_subarray_off_len (varp->ndim, varp->chunkdim, tssizep, tstartp, &pboff,
											 &plen);
			if (err == 0) {
				plen *= varp->esize;
				pboff *= varp->esize;
#ifdef PNETCDF_DEBUG
				if (rbufp[j] - rbuf[j] + plen > rsize[j]) { RET_ERR (NC_EINTERNAL) }
#endif
				memcpy (varp->chunk_cache[cid]->buf + pboff, rbufp[j], plen);
				rbufp[j] += plen;
			} else {
				// Pack type
				CHK_ERR_TYPE_CREATE_SUBARRAY (varp->ndim, varp->chunkdim, tssizep, tstartp,
											  MPI_ORDER_C, varp->etype, &ptype);
				CHK_ERR_TYPE_COMMIT (&ptype);

				// Pack data
				packoff = 0;
				CHK_ERR_UNPACK (rbufp[j], rsize[j], &packoff, varp->chunk_cache[cid]->buf, 1, ptype,
								ncchkp->comm);
				rbufp[j] += packoff;
				MPI_Type_free (&ptype);
			}

			// Mark chunk as dirty
			varp->dirty[cid] = 1;

#ifdef PNETCDF_PROFILING
			ncchkp->nreq++;
#endif
		}
		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_UNPACK_REQ)
	}

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_CB_SEND_REQ)

	CHK_ERR_WAITALL (nsend, sreq, sstat);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB_SEND_REQ)

err_out:;
	// Free buffers
	NCI_Free (wcnt_local);

	NCI_Free (tstart);

	NCI_Free (ostart);

	NCI_Free (sreq);
	NCI_Free (sstat);
	NCI_Free (ssize);
	if (nsend > 0) { NCI_Free (sbuf[0]); }
	NCI_Free (sbuf);

	NCI_Free (rreq);
	for (i = 0; i < nrecv; i++) { NCI_Free (rbuf[i]); }
	NCI_Free (rbuf);
	NCI_Free (rsize);

	NCI_Free (tbuf);

	NCI_Free (rlo_local);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_CB)

	return err;
}
