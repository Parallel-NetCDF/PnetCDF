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

#include "../ncmpio/ncmpio_NC.h"
#include "ncchkio_internal.h"

int ncchkioi_var_init_core (
	NC_chk *ncchkp, NC_chk_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
	int err = NC_NOERR;
	int ret;
	int i, j;
	int valid;
	MPI_Offset len;
	NC_chk_var *var;

	if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
		if (varp->chunkdim == NULL) {  // This is a new uninitialized variable
			// Init value
			varp->mychunks = NULL;	// To be added later

			// Update dimsize on rec dim
			if (ncchkp->recdim >= 0) {
				if (varp->dimsize[0] < ncchkp->recsize) { varp->dimsize[0] = ncchkp->recsize; }
			}

			// Determine its block size
			varp->chunkdim = (int *)NCI_Malloc (sizeof (int) * varp->ndim);
			varp->nchunks  = (int *)NCI_Malloc (sizeof (int) * varp->ndim);

			// First check attribute
			valid = 1;
			ret	  = ncchkp->driver->inq_att (ncchkp->ncp, varp->varid, "_chunkdim", NULL, &len);
			if (ret == NC_NOERR && len == varp->ndim) {
				ret = ncchkp->driver->get_att (ncchkp->ncp, varp->varid, "_chunkdim",
											   varp->chunkdim, MPI_INT);
				if (ret != NC_NOERR) { valid = 0; }
				// chunkdim must be at leasst 1
				for (j = 0; j < varp->ndim; j++) {
					if (varp->chunkdim[j] <= 0) {
						valid = 0;
						printf ("Warning: chunk size invalid, use default");
						break;
					}
				}
			} else {
				valid = 0;
			}

			// Now, try global default
			if ((!valid) && ncchkp->chunkdim) {
				valid = 1;
				for (i = 0; i < varp->ndim; i++) {
					if (ncchkp->chunkdim[varp->dimids[i]] > 0) {
						varp->chunkdim[i] = ncchkp->chunkdim[varp->dimids[i]];
					} else {
						valid = 0;
						break;
					}
				}
			}

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_META)

			// Still no clue, try to infer form I/O pattern (expensive)
			// If there is no I/O records, the default is just set to entire variable (only 1 chunk)
			if (!valid) {
				// Infering not supported
				err = ncchkioi_calc_chunk_size (ncchkp, varp, nreq, starts, counts);
				CHK_ERR
			}

			NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_META)

			// Calculate total # chunks, # chunks along each dim, chunksize
			varp->nchunkrec = 1;
			varp->chunksize = NC_Type_size (varp->xtype);
			for (i = 0; i < varp->ndim; i++) {	// chunkdim must be at leasst 1
				if (varp->dimsize[i] % varp->chunkdim[i] == 0) {
					varp->nchunks[i] = (int)(varp->dimsize[i] / (MPI_Offset)varp->chunkdim[i]);
				} else {
					varp->nchunks[i] = (int)(varp->dimsize[i] / (MPI_Offset)varp->chunkdim[i] + 1);
				}
				if (i > 0) { varp->nchunkrec *= varp->nchunks[i]; }
				varp->chunksize *= varp->chunkdim[i];
			}
			if (varp->isrec) {
				varp->nrec		= varp->nchunks[0];
				varp->nrecalloc = ncchkp->default_recnalloc;
				while (varp->nrecalloc < varp->nchunks[0]) {
					varp->nrecalloc *= NC_CHK_REC_MULTIPLIER;
				}
			} else {
				varp->nrec		= 1;
				varp->nrecalloc = 1;
				varp->nchunkrec *= varp->nchunks[0];
			}
			varp->nchunk	  = varp->nchunkrec * varp->nrec;
			varp->nchunkalloc = varp->nrecalloc * varp->nchunkrec;

			// Calculate number of chunks below each dimension
			varp->cidsteps				   = (int *)NCI_Malloc (sizeof (int) * varp->ndim);
			varp->cidsteps[varp->ndim - 1] = 1;
			for (i = varp->ndim - 2; i >= 0; i--) {
				varp->cidsteps[i] = varp->cidsteps[i + 1] * varp->nchunks[i + 1];
			}

			// Determine block ownership
			varp->dirty		  = (int *)NCI_Malloc (sizeof (int) * varp->nchunkalloc);
			varp->chunk_cache = (NC_chk_cache **)NCI_Malloc (sizeof (char *) * varp->nchunkalloc);
			memset (varp->chunk_cache, 0, sizeof (char *) * varp->nchunkalloc);
			memset (varp->dirty, 0, sizeof (int) * varp->nchunkalloc);

			// Block ownership to be decisded later
			varp->chunk_owner = (int *)NCI_Malloc (sizeof (int) * varp->nchunkalloc);

			// Determine block offset
			varp->chunk_index = (NC_chk_chunk_index_entry *)NCI_Malloc (
				sizeof (NC_chk_chunk_index_entry) * (varp->nchunkalloc + 1));

			// Try if there are offset recorded in attributes, it can happen after opening a file
			if (varp->isnew) {
				varp->metaoff = -1;
				;
				memset (varp->chunk_index, 0,
						sizeof (NC_chk_chunk_index_entry) * (varp->nchunk + 1));
			}

			/* Select compression driver based on attribute */
			ret = ncchkp->driver->inq_att (ncchkp->ncp, varp->varid, "_filter", NULL, &len);
			if (ret == NC_NOERR && len == 1) {
				ret = ncchkp->driver->get_att (ncchkp->ncp, varp->varid, "_filter",
											   &(varp->filter), MPI_INT);
				if (ret != NC_NOERR) { return err; }
			} else {
				varp->filter = ncchkp->default_filter;
			}
			switch (varp->filter) {
				case NC_CHK_FILTER_NONE:
					varp->filter_driver = NULL;
					break;
				case NC_CHK_FILTER_DUMMY:
					varp->filter_driver = ncchk_dummy_inq_driver ();
					break;
#ifdef ENABLE_ZLIB
				case NC_CHK_FILTER_ZLIB:
					varp->filter_driver = ncchk_zlib_inq_driver ();
					break;
#endif
#ifdef ENABLE_SZ
				case NC_CHK_FILTER_SZ:
					varp->filter_driver = ncchk_sz_inq_driver ();
					break;
#endif
				default:
					if (ncchkp->rank == 0) {
						printf ("Warning: Unknown filter driver id %d, use NC_CHK_FILTER_DUMMY\n",
								varp->filter);
					}
					varp->filter_driver = ncchk_dummy_inq_driver ();
					break;
					break;
			}

			// Update max ndim and chunksize
			if (ncchkp->max_ndim < varp->ndim) { ncchkp->max_ndim = varp->ndim; }
			if (ncchkp->max_chunk_size < varp->chunksize) {
				ncchkp->max_chunk_size = varp->chunksize;
			}

			if (ncchkp->cache_limit_hint == -1) {
				ncchkp->cache_limit += (size_t) (varp->nmychunkrec) * (size_t) (varp->chunksize);
			}
		}
	}

err_out:;
	return err;
}

int ncchkioi_var_init (
	NC_chk *ncchkp, NC_chk_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
	int err;

	err = ncchkioi_var_init_core (ncchkp, varp, nreq, starts, counts);
	CHK_ERR

	if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
		err = ncchkioi_calc_chunk_owner (ncchkp, varp, nreq, starts, counts);
		CHK_ERR
	}

err_out:;
	return err;
}

void ncchkioi_var_free (NC_chk_var *varp) {
	int i;

	if (varp->chunkdim != NULL) {
		NCI_Free (varp->dimsize);
		NCI_Free (varp->chunkdim);
		NCI_Free (varp->dimids);
		NCI_Free (varp->nchunks);
		NCI_Free (varp->cidsteps);
		NCI_Free (varp->chunk_index);
		NCI_Free (varp->chunk_owner);
		NCI_Free (varp->dirty);
		// for(i = 0; i < varp->nmychunk; i++){
		//    if (varp->chunk_cache[varp->mychunks[i]] != NULL){
		//        NCI_Free(varp->chunk_cache[varp->mychunks[i]]);
		//    }
		//}
		NCI_Free (varp->chunk_cache);
		NCI_Free (varp->mychunks);
	}
}

int ncchkioi_init_nvar_core_gather (NC_chk *ncchkp,
									int nvar,
									NC_chk_var **varps,
									int *rcnt,
									int *roff,
									MPI_Offset **starts,
									MPI_Offset **counts) {
	int err = NC_NOERR;
	int i, j;
	NC_chk_var *varp;
	ncchkioi_chunk_overlap_t *ocnt[2], *ocnt_all[2];
	size_t ocnt_size[2];
	MPI_Status stat;
	MPI_Request req;

	// Iinit vars
	ocnt_size[0] = ocnt_size[1] = 0;
	ocnt[0] = ocnt[1] = NULL;
	for (i = 0; i < nvar; i++) {
		varp = varps[i];
		j	 = i & 1;

		err = ncchkioi_var_init_core (ncchkp, varp, rcnt[i], starts + roff[i], counts + roff[i]);
		CHK_ERR

		if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
			if (varp->nchunkrec > ocnt_size[j]) {
				ocnt_size[j] = varp->nchunkrec;
				NCI_Free (ocnt[j]);
				ocnt[j] = (ncchkioi_chunk_overlap_t *)NCI_Malloc (
					sizeof (ncchkioi_chunk_overlap_t) * varp->nchunkrec * 2);
				ocnt_all[j] = ocnt[j] + varp->nchunkrec;
			}

			err = ncchkioi_calc_chunk_overlap (ncchkp, varp, rcnt[i], starts, counts, ocnt[j]);
			CHK_ERR
		}

		if ((i > 0) && (req != MPI_REQUEST_NULL)) {	 // Wait comm for prev var
			err = MPI_Wait (&req, &stat);
			ncchkioi_assign_chunk_owner (ncchkp, varps[i - 1], ocnt_all[(i - 1) & 1]);
			ncchkioi_write_chunk_ocnt (ncchkp, varps[i - 1], ocnt[(i - 1) & 1],
									   sizeof (ncchkioi_chunk_overlap_t));
		}

		if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
			err = ncchkioi_sync_ocnt_reduce (ncchkp, varp->nchunkrec, ocnt[j], ocnt_all[j], &req);
			CHK_ERR
		} else {
			req = MPI_REQUEST_NULL;
		}
	}
	// Last var
	if (req != MPI_REQUEST_NULL) {
		err = MPI_Wait (&req, &stat);
		ncchkioi_assign_chunk_owner (ncchkp, varp, ocnt_all[(i - 1) & 1]);
		ncchkioi_write_chunk_ocnt (ncchkp, varp, ocnt[(i - 1) & 1],
								   sizeof (ncchkioi_chunk_overlap_t));
	}

err_out:;
	NCI_Free (ocnt[0]);
	NCI_Free (ocnt[1]);
	return err;
}

int ncchkioi_init_nvar_core_reduce (NC_chk *ncchkp,
									int nvar,
									NC_chk_var **varps,
									int *rcnt,
									int *roff,
									MPI_Offset **starts,
									MPI_Offset **counts) {
	int err = NC_NOERR;
	int i, j;
	NC_chk_var *varp;
	ncchkioi_chunk_overlap_t *ocnt[2], *ocnt_all[2];
	size_t ocnt_size[2];
	MPI_Status stat;
	MPI_Request req;

	// Iinit vars
	ocnt_size[0] = ocnt_size[1] = 0;
	ocnt[0] = ocnt[1] = NULL;
	for (i = 0; i < nvar; i++) {
		varp = varps[i];
		j	 = i & 1;

		err = ncchkioi_var_init_core (ncchkp, varp, rcnt[i], starts + roff[i], counts + roff[i]);
		CHK_ERR

		NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_COWN)

		if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
			if (varp->nchunkrec > ocnt_size[j]) {
				ocnt_size[j] = varp->nchunkrec;
				NCI_Free (ocnt[j]);
				ocnt[j] = (ncchkioi_chunk_overlap_t *)NCI_Malloc (
					sizeof (ncchkioi_chunk_overlap_t) * varp->nchunkrec * 2);
				ocnt_all[j] = ocnt[j] + varp->nchunkrec;
			}

			err = ncchkioi_calc_chunk_overlap (ncchkp, varp, rcnt[i], starts + roff[i],
											   counts + roff[i], ocnt[j]);
			CHK_ERR
		}

		if ((i > 0) && (req != MPI_REQUEST_NULL)) {	 // Wait comm for prev var
			err = MPI_Wait (&req, &stat);
			ncchkioi_assign_chunk_owner (ncchkp, varps[i - 1], ocnt_all[(i - 1) & 1]);
			ncchkioi_write_chunk_ocnt (ncchkp, varps[i - 1], ocnt[(i - 1) & 1],
									   sizeof (ncchkioi_chunk_overlap_t));
		}

		if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
			ncchkioi_sync_ocnt_reduce (ncchkp, varp->nchunkrec, ocnt[j], ocnt_all[j], &req);
		} else {
			req = MPI_REQUEST_NULL;
		}

		NC_CHK_TIMER_STOPEX (NC_CHK_TIMER_VAR_INIT_COWN, NC_CHK_TIMER_VAR_INIT_META)
	}
	// Last var
	NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_COWN)
	if (req != MPI_REQUEST_NULL) {
		err = MPI_Wait (&req, &stat);
		ncchkioi_assign_chunk_owner (ncchkp, varp, ocnt_all[(i - 1) & 1]);
		ncchkioi_write_chunk_ocnt (ncchkp, varp, ocnt[(i - 1) & 1],
								   sizeof (ncchkioi_chunk_overlap_t));
	}
	NC_CHK_TIMER_STOPEX (NC_CHK_TIMER_VAR_INIT_COWN, NC_CHK_TIMER_VAR_INIT_META)

err_out:;
	NCI_Free (ocnt[0]);
	NCI_Free (ocnt[1]);
	return err;
}

int ncchkioi_init_nvar (NC_chk *ncchkp, int nput, int *putreqs, int nget, int *getreqs) {
	int err = NC_NOERR, ret;
	int i, j;
	int nflag;
	unsigned int *flag, *flag_all;
	int nvar;
	int *vmap;
	NC_chk_var *varp;
	NC_chk_var **varps;
	int *rcnt, *roff;
	MPI_Offset **starts, **counts;
	NC_chk_req *req;
	int nread;
	int *lens;
	MPI_Aint *fdisps, *mdisps;
	MPI_Datatype ftype, mtype;
	MPI_Status status;

	NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_META)

	CHK_ERR_ALLREDUCE (MPI_IN_PLACE, &(ncchkp->recsize), 1, MPI_LONG_LONG, MPI_MAX,
					   ncchkp->comm);  // Sync number of recs

	// Flag of touched vars
	nflag = ncchkp->vars.cnt / 32 + 1;
	flag  = (unsigned int *)NCI_Malloc (sizeof (int) * nflag * 2);
	CHK_PTR (flag)
	flag_all = flag + nflag;
	memset (flag, 0, sizeof (int) * nflag);
	for (i = 0; i < nput; i++) {
		req = ncchkp->putlist.reqs + putreqs[i];
		flag[req->varid >> 5] |= 1u << (req->varid % 32);
	}
	for (i = 0; i < nget; i++) {
		req = ncchkp->getlist.reqs + getreqs[i];
		flag[req->varid >> 5] |= 1u << (req->varid % 32);
	}

	// Sync flag
	CHK_ERR_ALLREDUCE (flag, flag_all, nflag, MPI_UNSIGNED, MPI_BOR, ncchkp->comm);

	// Build a skip list of touched vars
	nvar = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		if (flag_all[i >> 5] & (1u << (i % 32))) {
			if ((ncchkp->vars.data + i)->chunkdim == NULL) {  // If not yet inited
				nvar++;
			} else {
				flag_all[i >> 5] ^= (1u << (i % 32));
				if ((ncchkp->vars.data + i)->dimsize[0] < ncchkp->recsize) {
					ncchkioi_var_resize (ncchkp, ncchkp->vars.data + i);
				}
			}
		}
	}
	varps = (NC_chk_var **)NCI_Malloc (sizeof (NC_chk_var *) * nvar);
	CHK_PTR (varps)
	vmap = (int *)NCI_Malloc (sizeof (int) * ncchkp->vars.cnt);
	CHK_PTR (vmap)
	nvar = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		if (flag_all[i >> 5] & (1u << (i % 32))) {
			varps[nvar] = ncchkp->vars.data + i;
			vmap[i]		= nvar++;
		}
	}

	// Count reqs for each var
	roff = (int *)NCI_Malloc (sizeof (int) * (nvar + 1));
	CHK_PTR (roff)
	rcnt = (int *)NCI_Malloc (sizeof (int) * nvar);
	CHK_PTR (rcnt)
	memset (rcnt, 0, sizeof (int) * nvar);
	for (i = 0; i < nput; i++) {
		req = ncchkp->putlist.reqs + putreqs[i];
		j	= req->varid;
		if (flag_all[j >> 5] & (1u << (j % 32))) { rcnt[vmap[j]] += req->nreq; }
	}
	for (i = 0; i < nget; i++) {
		req = ncchkp->getlist.reqs + getreqs[i];
		j	= req->varid;
		if (flag_all[j >> 5] & (1u << (j % 32))) { rcnt[vmap[j]] += req->nreq; }
	}
	roff[0] = 0;
	for (i = 0; i < nvar; i++) { roff[i + 1] = roff[i] + rcnt[i]; }

	// Gather starts and counts
	starts = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *) * roff[nvar] * 2);
	CHK_PTR (starts)
	counts = starts + roff[nvar];
	memset (rcnt, 0, sizeof (int) * nvar);
	for (i = 0; i < nput; i++) {
		req = ncchkp->putlist.reqs + putreqs[i];
		j	= req->varid;
		if (flag_all[j >> 5] & (1u << (j % 32))) {
			j = vmap[req->varid];
			if (req->nreq > 1) {
				memcpy (starts + roff[j] + rcnt[j], req->starts, sizeof (MPI_Offset *) * req->nreq);
				memcpy (counts + roff[j] + rcnt[j], req->counts, sizeof (MPI_Offset *) * req->nreq);
				rcnt[j] += req->nreq;
			} else {
				starts[roff[j] + rcnt[j]]	  = req->start;
				counts[roff[j] + (rcnt[j]++)] = req->count;
			}
		}
	}
	for (i = 0; i < nget; i++) {
		req = ncchkp->getlist.reqs + getreqs[i];
		j	= req->varid;
		if (flag_all[j >> 5] & (1u << (j % 32))) {
			j = vmap[req->varid];
			if (req->nreq > 1) {
				memcpy (starts + roff[j] + rcnt[j], req->starts, sizeof (MPI_Offset *) * req->nreq);
				memcpy (counts + roff[j] + rcnt[j], req->counts, sizeof (MPI_Offset *) * req->nreq);
				rcnt[j] += req->nreq;
			} else {
				starts[roff[j] + rcnt[j]]	  = req->start;
				counts[roff[j] + (rcnt[j]++)] = req->count;
			}
		}
	}

	// Buffer for index table type
	lens = NCI_Malloc (sizeof (int) * nvar);
	CHK_PTR (lens)
	fdisps = NCI_Malloc (sizeof (MPI_Aint) * nvar * 2);
	CHK_PTR (fdisps)
	mdisps = fdisps + nvar;
	nread  = 0;

	// Iinit vars
	ncchkp->cown_size = 0;	// Reset owner penalty
	err = ncchkioi_init_nvar_core_reduce (ncchkp, nvar, varps, rcnt, roff, starts, counts);
	CHK_ERR

	// Read the index table for existing variables
	// MPI Type to load the index table for existing variables
	for (i = 0; i < nvar; i++) {
		varp = varps[i];
		if (!(varp->isnew)) {
			ret = ncchkp->driver->get_att (ncchkp->ncp, varp->varid, "_metaoffset",
										   &(varp->metaoff), MPI_LONG_LONG);
			if (ret == NC_NOERR) {
				lens[nread]		= sizeof (NC_chk_chunk_index_entry) * (varp->nchunk);
				fdisps[nread]	= varp->metaoff;
				mdisps[nread++] = (MPI_Aint) (varp->chunk_index);
			} else {
				varp->metaoff = -1;
				memset (varp->chunk_index, 0,
						sizeof (NC_chk_chunk_index_entry) * (varp->nchunk + 1));
			}
		}
	}
	if (nread) {
		ncchkioi_sort_file_offset (nread, fdisps, mdisps, lens);

		MPI_Type_create_hindexed (nread, lens, fdisps, MPI_BYTE, &ftype);
		CHK_ERR_TYPE_COMMIT (&ftype);

		MPI_Type_create_hindexed (nread, lens, mdisps, MPI_BYTE, &mtype);
		CHK_ERR_TYPE_COMMIT (&mtype);

		// Set file view
		CHK_ERR_SET_VIEW (((NC *)(ncchkp->ncp))->collective_fh, ((NC *)(ncchkp->ncp))->begin_var,
						  MPI_BYTE, ftype, "native", MPI_INFO_NULL);

		// Read data
		CHK_ERR_READ_AT_ALL (((NC *)(ncchkp->ncp))->collective_fh, 0, MPI_BOTTOM, 1, mtype,
							 &status);

		// Restore file view
		CHK_ERR_SET_VIEW (((NC *)(ncchkp->ncp))->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native",
						  MPI_INFO_NULL);

#ifdef WORDS_BIGENDIAN	// Switch back to little endian
		for (i = 0; i < nvar; i++) {
			ncchkioi_idx_in_swapn (varps[i]->chunk_index, varps[i]->nchunk + 1);
		}
#endif

		MPI_Type_free (&ftype);
		MPI_Type_free (&mtype);
	}

	NCI_Free (lens);
	NCI_Free (fdisps);

	NCI_Free (flag);
	NCI_Free (varps);
	NCI_Free (vmap);
	NCI_Free (roff);
	NCI_Free (rcnt);
	NCI_Free (starts);

err_out:;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_META)
	return err;
}
