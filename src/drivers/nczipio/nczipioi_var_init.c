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
#include <nczipio_driver.h>
#include <pnc_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../ncmpio/ncmpio_NC.h"
#include "nczipio_internal.h"

int nczipioi_var_init_core (
	NC_zip *nczipp, NC_zip_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
	int err = NC_NOERR;
	int ret;
	int i, j;
	int valid;
	MPI_Offset len;
	NC_zip_var *var;

	if (varp->varkind == NC_ZIP_VAR_COMPRESSED) {
		if (varp->chunkdim == NULL) {  // This is a new uninitialized variable
			// Update dimsize on rec dim
			if (nczipp->recdim >= 0) {
				if (varp->dimsize[0] < nczipp->recsize) { varp->dimsize[0] = nczipp->recsize; }
			}

			// Determine its block size
			varp->chunkdim = (int *)NCI_Malloc (sizeof (int) * varp->ndim);
			varp->nchunks  = (int *)NCI_Malloc (sizeof (int) * varp->ndim);

			// First check attribute
			valid = 1;
			ret	  = nczipp->driver->inq_att (nczipp->ncp, varp->varid, "_chunkdim", NULL, &len);
			if (ret == NC_NOERR && len == varp->ndim) {
				ret = nczipp->driver->get_att (nczipp->ncp, varp->varid, "_chunkdim",
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
			if ((!valid) && nczipp->chunkdim) {
				valid = 1;
				for (i = 0; i < varp->ndim; i++) {
					if (nczipp->chunkdim[varp->dimids[i]] > 0) {
						varp->chunkdim[i] = nczipp->chunkdim[varp->dimids[i]];
					} else {
						valid = 0;
						break;
					}
				}
			}

			NC_ZIP_TIMER_STOP (NC_ZIP_TIMER_VAR_INIT_META)

			// Still no clue, try to infer form I/O pattern (expensive)
			// If there is no I/O records, the default is just set to entire variable (only 1 chunk)
			if (!valid) {
				// Infering not supported
				err = nczipioi_calc_chunk_size (nczipp, varp, nreq, starts, counts);
				CHK_ERR
			}

			NC_ZIP_TIMER_START (NC_ZIP_TIMER_VAR_INIT_META)

			// Calculate total # chunks, # chunks along each dim, chunksize
			varp->nchunkrec = 1;
			varp->chunksize = NC_Type_size (varp->xtype);
			for (i = 0; i < varp->ndim; i++) {	// chunkdim must be at leasst 1
				if (varp->dimsize[i] % varp->chunkdim[i] == 0) {
					varp->nchunks[i] = (int)varp->dimsize[i] / varp->chunkdim[i];
				} else {
					varp->nchunks[i] = (int)varp->dimsize[i] / varp->chunkdim[i] + 1;
				}
				if (i > 0) { varp->nchunkrec *= varp->nchunks[i]; }
				varp->chunksize *= varp->chunkdim[i];
			}
			if (varp->isrec) {
				varp->nrec		= varp->nchunks[0];
				varp->nrecalloc = nczipp->default_recnalloc;
				while (varp->nrecalloc < varp->nchunks[0]) {
					varp->nrecalloc *= NC_ZIP_REC_MULTIPLIER;
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
			varp->chunk_cache = (NC_zip_cache **)NCI_Malloc (sizeof (char *) * varp->nchunkalloc);
			memset (varp->chunk_cache, 0, sizeof (char *) * varp->nchunkalloc);
			memset (varp->dirty, 0, sizeof (int) * varp->nchunkalloc);

			// Block ownership to be decisded later
			varp->chunk_owner = (int *)NCI_Malloc (sizeof (int) * varp->nchunkalloc);

			// Determine block offset
			varp->chunk_index = (NC_zip_chunk_index_entry *)NCI_Malloc (
				sizeof (NC_zip_chunk_index_entry) * (varp->nchunkalloc + 1));

			// Try if there are offset recorded in attributes, it can happen after opening a file
			if (varp->isnew) {
				varp->metaoff = -1;
				;
				memset (varp->chunk_index, 0,
						sizeof (NC_zip_chunk_index_entry) * (varp->nchunk + 1));
			}

			/* Select compression driver based on attribute */
			ret = nczipp->driver->inq_att (nczipp->ncp, varp->varid, "_zipdriver", NULL, &len);
			if (ret == NC_NOERR && len == 1) {
				ret = nczipp->driver->get_att (nczipp->ncp, varp->varid, "_zipdriver",
											   &(varp->zipdriver), MPI_INT);
				if (ret != NC_NOERR) { return err; }
			} else {
				varp->zipdriver = nczipp->default_zipdriver;
			}
			switch (varp->zipdriver) {
				case NC_ZIP_DRIVER_NONE:
					varp->zip = NULL;
					break;
				case NC_ZIP_DRIVER_DUMMY:
					varp->zip = nczip_dummy_inq_driver ();
					break;
#ifdef ENABLE_ZLIB
				case NC_ZIP_DRIVER_ZLIB:
					varp->zip = nczip_zlib_inq_driver ();
					break;
#endif
#ifdef ENABLE_SZ
				case NC_ZIP_DRIVER_SZ:
					varp->zip = nczip_sz_inq_driver ();
					break;
#endif
				default:
					if (nczipp->rank == 0) {
						printf ("Warning: Unknown zip driver id %d, use NC_ZIP_DRIVER_DUMMY\n",
								varp->zipdriver);
					}
					varp->zip = nczip_dummy_inq_driver ();
					break;
					break;
			}

			// Update max ndim and chunksize
			if (nczipp->max_ndim < varp->ndim) { nczipp->max_ndim = varp->ndim; }
			if (nczipp->max_chunk_size < varp->chunksize) {
				nczipp->max_chunk_size = varp->chunksize;
			}

			if (nczipp->cache_limit_hint == -1) {
				nczipp->cache_limit += (size_t) (varp->nmychunkrec) * (size_t) (varp->chunksize);
			}
		}
	}

err_out:;
	return err;
}

int nczipioi_var_init (
	NC_zip *nczipp, NC_zip_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
	int err;

	err = nczipioi_var_init_core (nczipp, varp, nreq, starts, counts);
	CHK_ERR
	err = nczipioi_calc_chunk_owner (nczipp, varp, nreq, starts, counts);
	CHK_ERR

err_out:;
	return err;
}

void nczipioi_var_free (NC_zip_var *varp) {
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

int nczipioi_init_nvar_core_gather (NC_zip *nczipp,
									int nvar,
									NC_zip_var **varps,
									int *rcnt,
									int *roff,
									MPI_Offset **starts,
									MPI_Offset **counts) {
	int err = NC_NOERR;
	int i, j;
	NC_zip_var *varp;
	nczipioi_chunk_overlap_t *ocnt[2], *ocnt_all[2];
	size_t ocnt_size[2];
	MPI_Status stat;
	MPI_Request req;

	// Iinit vars
	ocnt_size[0] = ocnt_size[1] = 0;
	ocnt[0] = ocnt[1] = NULL;
	for (i = 0; i < nvar; i++) {
		varp = varps[i];
		j	 = i & 1;

		err = nczipioi_var_init_core (nczipp, varp, rcnt[i], starts + roff[i], counts + roff[i]);
		CHK_ERR

		if (varp->nchunkrec > ocnt_size[j]) {
			ocnt_size[j] = varp->nchunkrec;
			NCI_Free (ocnt[j]);
			ocnt[j] = (nczipioi_chunk_overlap_t *)NCI_Malloc (sizeof (nczipioi_chunk_overlap_t) *
															  varp->nchunkrec * 2);
			ocnt_all[j] = ocnt[j] + varp->nchunkrec;
		}

		err = nczipioi_calc_chunk_overlap (nczipp, varp, rcnt[i], starts, counts, ocnt[j]);
		CHK_ERR

		if (i > 0) {  // Wait comm for prev var
			err = MPI_Wait (&req, &stat);
			nczipioi_assign_chunk_owner (nczipp, varps[i - 1], ocnt_all[(i - 1) & 1]);
			nczipioi_write_chunk_ocnt (nczipp, varps[i - 1], ocnt[(i - 1) & 1],
									   sizeof (nczipioi_chunk_overlap_t));
		}

		nczipioi_sync_ocnt_reduce (nczipp, varp->nchunkrec, ocnt[j], ocnt_all[j], &req);
	}
	// Last var
	err = MPI_Wait (&req, &stat);
	nczipioi_assign_chunk_owner (nczipp, varp, ocnt_all[(i - 1) & 1]);
	nczipioi_write_chunk_ocnt (nczipp, varp, ocnt[(i - 1) & 1], sizeof (nczipioi_chunk_overlap_t));

err_out:;
	NCI_Free (ocnt[0]);
	NCI_Free (ocnt[1]);
	return err;
}

int nczipioi_init_nvar_core_reduce (NC_zip *nczipp,
									int nvar,
									NC_zip_var **varps,
									int *rcnt,
									int *roff,
									MPI_Offset **starts,
									MPI_Offset **counts) {
	int err = NC_NOERR;
	int i, j;
	NC_zip_var *varp;
	nczipioi_chunk_overlap_t *ocnt[2], *ocnt_all[2];
	size_t ocnt_size[2];
	MPI_Status stat;
	MPI_Request req;

	// Iinit vars
	ocnt_size[0] = ocnt_size[1] = 0;
	ocnt[0] = ocnt[1] = NULL;
	for (i = 0; i < nvar; i++) {
		varp = varps[i];
		j	 = i & 1;

		err = nczipioi_var_init_core (nczipp, varp, rcnt[i], starts + roff[i], counts + roff[i]);
		CHK_ERR

		if (varp->nchunkrec > ocnt_size[j]) {
			ocnt_size[j] = varp->nchunkrec;
			NCI_Free (ocnt[j]);
			ocnt[j] = (nczipioi_chunk_overlap_t *)NCI_Malloc (sizeof (nczipioi_chunk_overlap_t) *
															  varp->nchunkrec * 2);
			ocnt_all[j] = ocnt[j] + varp->nchunkrec;
		}
		
		NC_ZIP_TIMER_START (NC_ZIP_TIMER_VAR_INIT_COWN)

		err = nczipioi_calc_chunk_overlap (nczipp, varp, rcnt[i], starts, counts, ocnt[j]);
		CHK_ERR

		if (i > 0) {  // Wait comm for prev var
			err = MPI_Wait (&req, &stat);
			nczipioi_assign_chunk_owner (nczipp, varps[i - 1], ocnt_all[(i - 1) & 1]);
			nczipioi_write_chunk_ocnt (nczipp, varps[i - 1], ocnt[(i - 1) & 1],
									   sizeof (nczipioi_chunk_overlap_t));
		}

		nczipioi_sync_ocnt_reduce (nczipp, varp->nchunkrec, ocnt[j], ocnt_all[j], &req);

		NC_ZIP_TIMER_STOPEX (NC_ZIP_TIMER_VAR_INIT_COWN, NC_ZIP_TIMER_VAR_INIT_META)
	}
	// Last var
	NC_ZIP_TIMER_START (NC_ZIP_TIMER_VAR_INIT_COWN)
	err = MPI_Wait (&req, &stat);
	nczipioi_assign_chunk_owner (nczipp, varp, ocnt_all[(i - 1) & 1]);
	nczipioi_write_chunk_ocnt (nczipp, varp, ocnt[(i - 1) & 1], sizeof (nczipioi_chunk_overlap_t));
	NC_ZIP_TIMER_STOPEX (NC_ZIP_TIMER_VAR_INIT_COWN, NC_ZIP_TIMER_VAR_INIT_META)

err_out:;
	NCI_Free (ocnt[0]);
	NCI_Free (ocnt[1]);
	return err;
}

int nczipioi_init_nvar (NC_zip *nczipp, int nput, int *putreqs, int nget, int *getreqs) {
	int err = NC_NOERR;
	int i, j;
	int nflag;
	unsigned int *flag, *flag_all;
	int nvar;
	int *vmap;
	NC_zip_var *varp;
	NC_zip_var **varps;
	int *rcnt, *roff;
	MPI_Offset **starts, **counts;
	NC_zip_req *req;
	int nread;
	int *lens;
	MPI_Aint *fdisps, *mdisps;
	MPI_Datatype ftype, mtype;
	MPI_Status status;

	NC_ZIP_TIMER_START (NC_ZIP_TIMER_VAR_INIT_META)

	CHK_ERR_ALLREDUCE (MPI_IN_PLACE, &(nczipp->recsize), 1, MPI_LONG_LONG, MPI_MAX,
					   nczipp->comm);  // Sync number of recs

	// Flag of touched vars
	nflag = nczipp->vars.cnt / 32 + 1;
	flag  = (unsigned int *)NCI_Malloc (sizeof (int) * nflag * 2);
	CHK_PTR (flag)
	flag_all = flag + nflag;
	memset (flag, 0, sizeof (int) * nflag);
	for (i = 0; i < nput; i++) {
		req = nczipp->putlist.reqs + putreqs[i];
		flag[req->varid >> 5] |= 1u << (req->varid % 32);
	}
	for (i = 0; i < nget; i++) {
		req = nczipp->getlist.reqs + getreqs[i];
		flag[req->varid >> 5] |= 1u << (req->varid % 32);
	}

	// Sync flag
	CHK_ERR_ALLREDUCE (flag, flag_all, nflag, MPI_UNSIGNED, MPI_BOR, nczipp->comm);

	// Build a skip list of touched vars
	nvar = 0;
	for (i = 0; i < nczipp->vars.cnt; i++) {
		if (flag_all[i >> 5] & (1u << (i % 32))) {
			if ((nczipp->vars.data + i)->chunkdim == NULL) {  // If not yet inited
				nvar++;
			} else {
				flag_all[i >> 5] ^= (1u << (i % 32));
				if ((nczipp->vars.data + i)->dimsize[0] < nczipp->recsize) {
					nczipioi_var_resize (nczipp, nczipp->vars.data + i);
				}
			}
		}
	}
	varps = (NC_zip_var **)NCI_Malloc (sizeof (NC_zip_var *) * nvar);
	CHK_PTR (varps)
	vmap = (int *)NCI_Malloc (sizeof (int) * nczipp->vars.cnt);
	CHK_PTR (vmap)
	nvar = 0;
	for (i = 0; i < nczipp->vars.cnt; i++) {
		if (flag_all[i >> 5] & (1u << (i % 32))) {
			varps[nvar] = nczipp->vars.data + i;
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
		req = nczipp->putlist.reqs + putreqs[i];
		j	= req->varid;
		if (flag_all[j >> 5] & (1u << (j % 32))) { rcnt[vmap[j]] += req->nreq; }
	}
	for (i = 0; i < nget; i++) {
		req = nczipp->getlist.reqs + getreqs[i];
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
		req = nczipp->putlist.reqs + putreqs[i];
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
		req = nczipp->getlist.reqs + getreqs[i];
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
	nczipp->cown_size=0;	// Reset owner penalty
	err = nczipioi_init_nvar_core_reduce (nczipp, nvar, varps, rcnt, roff, starts, counts);
	CHK_ERR

	// Read the index table for existing variables
	// MPI Type to load the index table for existing variables
	for (i = 0; i < nvar; i++) {
		varp = varps[i];
		if (!(varp->isnew)) {
			err = nczipp->driver->get_att (nczipp->ncp, varp->varid, "_metaoffset",
										   &(varp->metaoff), MPI_LONG_LONG);
			if (err == NC_NOERR) {
				lens[nread]		= sizeof (NC_zip_chunk_index_entry) * (varp->nchunk);
				fdisps[nread]	= varp->metaoff;
				mdisps[nread++] = (MPI_Aint) (varp->chunk_index);
			} else {
				varp->metaoff = -1;
				memset (varp->chunk_index, 0,
						sizeof (NC_zip_chunk_index_entry) * (varp->nchunk + 1));
			}
		}
	}
	if (nread) {
		nczipioi_sort_file_offset (nread, fdisps, mdisps, lens);

		MPI_Type_create_hindexed (nread, lens, fdisps, MPI_BYTE, &ftype);
		CHK_ERR_TYPE_COMMIT (&ftype);

		MPI_Type_create_hindexed (nread, lens, mdisps, MPI_BYTE, &mtype);
		CHK_ERR_TYPE_COMMIT (&mtype);

		// Set file view
		CHK_ERR_SET_VIEW (((NC *)(nczipp->ncp))->collective_fh, ((NC *)(nczipp->ncp))->begin_var,
						  MPI_BYTE, ftype, "native", MPI_INFO_NULL);

		// Read data
		CHK_ERR_READ_AT_ALL (((NC *)(nczipp->ncp))->collective_fh, 0, MPI_BOTTOM, 1, mtype,
							 &status);

		// Restore file view
		CHK_ERR_SET_VIEW (((NC *)(nczipp->ncp))->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native",
						  MPI_INFO_NULL);

#ifdef WORDS_BIGENDIAN	// Switch back to little endian
		for (i = 0; i < nvar; i++) {
			nczipioi_idx_in_swapn (varps[i]->chunk_index, varps[i]->nchunk + 1);
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

	NC_ZIP_TIMER_STOP (NC_ZIP_TIMER_VAR_INIT_META)
	return err;
}
