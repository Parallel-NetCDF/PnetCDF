/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

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

void ncchkioi_write_chunk_ocnt (NC_chk *ncchkp, NC_chk_var *varp, void *ocnt, size_t ocnt_size) {
#ifdef PNETCDF_PROFILING
	{
		int i, j;
		char *pprefix = getenv ("PNETCDF_OWNER_PREFIX");

		if (pprefix != NULL) {
			if (ncchkp->rank == 0) {
				void *ocnt_in;
				int *cown;
				MPI_Status stat;
				FILE *pfile;
				char fname[1024], ppath[1024];

				ocnt_in = NCI_Malloc (ocnt_size * varp->nchunkrec);
				cown	= NCI_Malloc (sizeof (int) * varp->nchunkrec);

				strcpy (fname, ncchkp->path);
				for (i = strlen (fname); i > 0; i--) {
					if (fname[i] == '.') {
						fname[i] = '\0';
					} else if (fname[i] == '\\' || fname[i] == '/') {
						i++;
						break;
					}
				}
				sprintf (ppath, "%s%s_owner.csv", pprefix, fname + i);
				pfile = fopen (ppath, "a");

				fprintf (pfile, "Var:, %d\n", varp->varid);
				fprintf (pfile, "Rank\\Chunk, ");
				for (j = 0; j < varp->nchunkrec; j++) { fprintf (pfile, "%d, ", j); }
				fprintf (pfile, "\nOwner, ");
				for (j = 0; j < varp->nchunk; j++) {
					fprintf (pfile, "%d, ", varp->chunk_owner[j]);
				}
				fprintf (pfile, "\n0, ");
				if (ocnt_size == sizeof (MPI_Offset)) {
					for (j = 0; j < varp->nchunkrec; j++) {
						fprintf (pfile, "%lld, ", ((MPI_Offset *)ocnt)[j]);
					}
				} else {
					for (j = 0; j < varp->nchunkrec; j++) {
						fprintf (pfile, "%lld, ", ((ncchkioi_chunk_overlap_t *)ocnt)[j].osize);
					}
				}
				fprintf (pfile, "\n");
				for (i = 1; i < ncchkp->np; i++) {
					if (ocnt_size == sizeof (MPI_Offset)) {
						MPI_Recv (ocnt_in, varp->nchunkrec, MPI_LONG_LONG, i, 0, ncchkp->comm,
								  &stat);
						fprintf (pfile, "%d, ", i);
						for (j = 0; j < varp->nchunkrec; j++) {
							fprintf (pfile, "%lld, ", ((MPI_Offset *)ocnt_in)[j]);
						}
					} else {
						MPI_Recv (ocnt_in, varp->nchunkrec, ncchkp->overlaptype, i, 0, ncchkp->comm,
								  &stat);
						fprintf (pfile, "%d, ", i);
						for (j = 0; j < varp->nchunkrec; j++) {
							fprintf (pfile, "%lld, ",
									 ((ncchkioi_chunk_overlap_t *)ocnt_in)[j].osize);
						}
					}
					fprintf (pfile, "\n");

					MPI_Recv (cown, varp->nchunkrec, MPI_INT, i, 0, ncchkp->comm, &stat);
					for (j = 0; j < varp->nchunkrec; j++) {
						if (cown[j] != varp->chunk_owner[j]) {
							printf ("Warning: cown[%d][%d] on rank %d = %d, != %d\n", varp->varid, j,
									i, cown[j], varp->chunk_owner[j]);
						}
					}
				}

				fclose (pfile);
				NCI_Free (ocnt_in);
				NCI_Free (cown);
			} else {
				if (ocnt_size == sizeof (MPI_Offset)) {
					MPI_Send (ocnt, varp->nchunkrec, MPI_LONG_LONG, 0, 0, ncchkp->comm);
				} else {
					MPI_Send (ocnt, varp->nchunkrec, ncchkp->overlaptype, 0, 0, ncchkp->comm);
				}
				MPI_Send (varp->chunk_owner, varp->nchunkrec, MPI_INT, 0, 0, ncchkp->comm);
			}
		}
	}
#endif
}

void max_osize_rank_op (void *inp, void *inoutp, int *len, MPI_Datatype *dptr) {
	int i;
	ncchkioi_chunk_overlap_t *in	= (ncchkioi_chunk_overlap_t *)inp;
	ncchkioi_chunk_overlap_t *inout = (ncchkioi_chunk_overlap_t *)inoutp;

	for (i = 0; i < *len; i++) {
		if (in->osize > inout->osize) {
			inout->osize = in->osize;
			inout->rank	 = in->rank;
		} else if ((in->osize == inout->osize) && (in->rank < inout->rank)) {
			inout->osize = in->osize;
			inout->rank	 = in->rank;
		}
		in++;
		inout++;
	}
}

int ncchkioi_calc_chunk_owner (
	NC_chk *ncchkp, NC_chk_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
	return ncchkioi_calc_chunk_owner_reduce (ncchkp, varp, nreq, starts, counts);
}

static inline void ncchkioi_rec_chunk_overlap (MPI_Offset *ostart,
											   MPI_Offset *osize,
											   MPI_Offset *citr,
											   NC_chk_var *varp,
											   MPI_Offset *ocnt,
											   NC_chk_req *reqp) {
	int i;
	int req;
	int cid;  // Chunk iterator
	MPI_Offset overlapsize;

	for (req = 0; req < reqp->nreq; req++) {
		ncchkioi_chunk_itr_init_ex (varp, reqp->starts[req], reqp->counts[req], citr, &cid, ostart,
									osize);	 // Initialize chunk iterator
		do {
			if (cid < varp->nchunkrec) {  // Count only first record
				// Count overlap
				overlapsize = 1;
				for (i = 0; i < varp->ndim; i++) { overlapsize *= osize[i]; }
				ocnt[cid] += (double)overlapsize;
				if (ocnt[cid] > varp->chunksize) { ocnt[cid] = (double)varp->chunksize; }
			}
		} while (ncchkioi_chunk_itr_next_ex (varp, reqp->starts[req], reqp->counts[req], citr, &cid,
											 ostart, osize));
	}
}

int ncchkioi_calc_chunk_overlap (NC_chk *ncchkp,
								 NC_chk_var *varp,
								 int nreq,
								 MPI_Offset **starts,
								 MPI_Offset **counts,
								 ncchkioi_chunk_overlap_t *ocnt) {
	int err = NC_NOERR;
	int i, j, k;
	int cid;  // Chunk iterator
	int req;
	MPI_Offset overlapsize;
	MPI_Offset *ostart, *osize;
	MPI_Offset *citr;  // Bounding box for chunks overlapping my own write region

	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * 3);
	CHK_PTR (ostart)
	osize = ostart + varp->ndim;
	citr  = osize + varp->ndim;

	memset (ocnt, 0, sizeof (ncchkioi_chunk_overlap_t) * varp->nchunkrec);

	// Count overlapsize of each request
	if (varp->isrec) {
		for (req = 0; req < nreq; req++) {
			ncchkioi_chunk_itr_init_ex (varp, starts[req], counts[req], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				if (cid < varp->nchunkrec) {  // Count only first record
					// Count overlap
					overlapsize = 1;
					for (i = 0; i < varp->ndim; i++) { overlapsize *= osize[i]; }
					ocnt[cid].osize += (double)overlapsize;
					if (ocnt[cid].osize > varp->chunksize) {
						ocnt[cid].osize = (double)varp->chunksize;
					}
				}
			} while (ncchkioi_chunk_itr_next_ex (varp, starts[req], counts[req], citr, &cid, ostart,
												 osize));
		}
	} else {
		for (req = 0; req < nreq; req++) {
			ncchkioi_chunk_itr_init_ex (varp, starts[req], counts[req], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				// Count overlap
				overlapsize = 1;
				for (i = 0; i < varp->ndim; i++) { overlapsize *= osize[i]; }
				ocnt[cid].osize += overlapsize;
				if (ocnt[cid].osize > varp->chunksize) { ocnt[cid].osize = varp->chunksize; }
			} while (ncchkioi_chunk_itr_next_ex (varp, starts[req], counts[req], citr, &cid, ostart,
												 osize));
		}
	}

	// First 16 bit used as noise
	for (i = 0; i < varp->nchunkrec; i++) {
		ocnt[i].rank = ncchkp->rank;
		ocnt[i].osize *= varp->esize;
		ocnt[i].osize <<= 16;
	}

	// Noise to break tie
	j = (ncchkp->rank - ncchkp->assigned_chunks) % ncchkp->np;
	if (j < 0) j += ncchkp->np;
	if (j > varp->nchunkrec) { j = varp->nchunkrec; }
	k = ncchkp->np - 1;	 // noise from 0 ~ np-1
	for (i = j; i < varp->nchunkrec; i++) {
		ocnt[i].osize += k;
		k--;
		if (k < 0) { k += ncchkp->np; }
	}
	for (i = 0; i < j; i++) {
		ocnt[i].osize += k;
		k--;
		if (k < 0) { k += ncchkp->np; }
	}
	ncchkp->assigned_chunks += varp->nchunk;

err_out:;
	NCI_Free (ostart);
	return err;
}

void ncchkioi_assign_chunk_owner (NC_chk *ncchkp,
								  NC_chk_var *varp,
								  ncchkioi_chunk_overlap_t *ocnt) {
	int i, j;
	for (i = 0; i < varp->nchunkrec; i++) { varp->chunk_owner[i] = ocnt[i].rank; }
	if (varp->isrec) {
		for (i = varp->nchunkrec; i < varp->nchunk; i += varp->nchunkrec) {
			memcpy (varp->chunk_owner + i, varp->chunk_owner, sizeof (int) * varp->nchunkrec);
		}
	}

	// Build skip list of my own chunks
	if (varp->nchunk > 0) {
		varp->nmychunkrec = 0;
		for (j = 0; j < varp->nchunkrec; j++) {
			if (varp->chunk_owner[j] == ncchkp->rank) { varp->nmychunkrec++; }
		}
		varp->nmychunk = varp->nmychunkrec * varp->nrec;
		varp->mychunks = (int *)NCI_Realloc (varp->mychunks, sizeof (int) * varp->nmychunkrec * varp->nrecalloc);
		varp->nmychunk = 0;
		for (j = 0; j < varp->nchunk; j++) {
			if (varp->chunk_owner[j] == ncchkp->rank) {
				varp->mychunks[varp->nmychunk++] = j;
				if (varp->isnew) {	// Only apply to new var, old var will be read when it is
									// needed
					// varp->chunk_cache[j] = (void*)NCI_Malloc(varp->chunksize);  // Allocate
					// buffer for blocks we own
					// memset(varp->chunk_cache[j], 0 , varp->chunksize);
				}
			}
		}
	} else {
		varp->nmychunk = varp->nmychunkrec = 0;
		varp->mychunks					   = NULL;
	}

	// Update global chunk count
	ncchkp->nmychunks += (MPI_Offset) (varp->nmychunk);
	ncchkp->cown_size +=
		(MPI_Offset) ((double)((MPI_Offset) (varp->nmychunk) * (MPI_Offset) (varp->chunksize)) *
					  ncchkp->cown_ratio);
}

int ncchkioi_sync_ocnt_reduce (NC_chk *ncchkp,
							   int nchunk,
							   ncchkioi_chunk_overlap_t *ocnt,
							   ncchkioi_chunk_overlap_t *ocnt_all,
							   MPI_Request *req) {
	int err = NC_NOERR;
	int i;

	// Construct MPI type for overlap if not already constructed
	if (ncchkp->overlaptype == MPI_DATATYPE_NULL) {
		err = MPI_Type_contiguous (sizeof (ncchkioi_chunk_overlap_t), MPI_BYTE,
								   &(ncchkp->overlaptype));
		CHK_MPIERR
		err = MPI_Type_commit (&(ncchkp->overlaptype));
		CHK_MPIERR
	}

	if (ncchkp->max_cown_op == MPI_OP_NULL) {
		err = MPI_Op_create (max_osize_rank_op, 1, &(ncchkp->max_cown_op));
		CHK_MPIERR
	}

	// Apply owner penalty
	for (i = 0; i < nchunk; i++) {
		ocnt[i].osize -= ncchkp->cown_size << 16;  // Penality for load ballance, set at 1/16
	}

	if (req) {
		CHK_ERR_IALLREDUCE (ocnt, ocnt_all, nchunk, ncchkp->overlaptype, ncchkp->max_cown_op,
							ncchkp->comm, req);
	} else {
		CHK_ERR_ALLREDUCE (ocnt, ocnt_all, nchunk, ncchkp->overlaptype, ncchkp->max_cown_op,
						   ncchkp->comm);
	}

err_out:;
	return err;
}

int ncchkioi_sync_ocnt_gather (NC_chk *ncchkp,
							   int nchunk,
							   ncchkioi_chunk_overlap_t *ocnt,
							   MPI_Offset **ocnt_all,
							   MPI_Request *req) {
	int err = NC_NOERR;

	// Construct MPI type for overlap if not already constructed
	if (ncchkp->overlaptype == MPI_DATATYPE_NULL) {
		MPI_Datatype tmptype;

		err = MPI_Type_contiguous (sizeof (MPI_Offset), MPI_BYTE, &tmptype);
		CHK_MPIERR
		err = MPI_Type_commit (&tmptype);
		CHK_MPIERR
		err = MPI_Type_create_resized (tmptype, 0, sizeof (ncchkioi_chunk_overlap_t),
									   &(ncchkp->overlaptype));
		CHK_MPIERR
		err = MPI_Type_commit (&(ncchkp->overlaptype));
		CHK_MPIERR
		err = MPI_Type_free (&tmptype);
	}

	if (req) {
		err = MPI_Igather (ocnt, nchunk, ncchkp->overlaptype, ocnt_all[0], nchunk, MPI_LONG_LONG, 0,
						   ncchkp->comm, req);
	} else {
		err = MPI_Gather (ocnt, nchunk, ncchkp->overlaptype, ocnt_all[0], nchunk, MPI_LONG_LONG, 0,
						  ncchkp->comm);
	}
	CHK_MPIERR

err_out:;
	return err;
}

int ncchkioi_sync_ocnt_gather_bcast (NC_chk *ncchkp,
									 NC_chk_var *varp,
									 MPI_Offset **ocnt_in,
									 ncchkioi_chunk_overlap_t *ocnt_all,
									 MPI_Request *req) {
	int err = NC_NOERR;
	int i, j, k;
	MPI_Offset *cown_size;

	if (ncchkp->rank == 0) {
		cown_size = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * ncchkp->np);
		memset (cown_size, 0, sizeof (MPI_Offset) * ncchkp->np);
		for (i = 0; i < varp->nchunkrec; i++) {
			ocnt_all[i].rank  = 0;
			ocnt_all[i].osize = ocnt_in[0][i];
			k				  = 0;
			for (j = 1; j < ncchkp->np; j++) {
				if (ocnt_in[j][i] - cown_size[j] > ocnt_in[k][i] - cown_size[k]) { k = j; }
			}
			cown_size[k] +=
				(MPI_Offset) ((double)(varp->chunksize) * ncchkp->cown_ratio) * varp->nrec;
			ocnt_all[i].rank  = k;
			ocnt_all[i].osize = ocnt_in[i][k];
		}
	}

	if (req) {
		err = MPI_Ibcast (ocnt_all, varp->nchunkrec, ncchkp->overlaptype, 0, ncchkp->comm, req);
	} else {
		err = MPI_Bcast (ocnt_all, varp->nchunkrec, ncchkp->overlaptype, 0, ncchkp->comm);
	}
	CHK_MPIERR

err_out:;
	return err;
}

int ncchkioi_calc_chunk_owner_reduce (
	NC_chk *ncchkp, NC_chk_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
	int err = NC_NOERR;
	int i, j, k;
	int cid;  // Chunk iterator
	int req;
	double noise, noise_step;
	MPI_Offset overlapsize;
	ncchkioi_chunk_overlap_t *ocnt, *ocnt_all;

	NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_COWN)

	ocnt = (ncchkioi_chunk_overlap_t *)NCI_Malloc (sizeof (ncchkioi_chunk_overlap_t) *
												   varp->nchunkrec * 2);
	CHK_PTR (ocnt)
	ocnt_all = ocnt + varp->nchunkrec;

	err = ncchkioi_calc_chunk_overlap (ncchkp, varp, nreq, starts, counts, ocnt);
	CHK_ERR

	if (ncchkp->exact_cown) {
		// err = ncchkioi_sync_ocnt_gather (ncchkp, varp->nchunkrec, ocnt, ocnt_all, NULL);
		// CHK_ERR
		RET_ERR (NC_ENOTSUPPORT)
	} else {
		err = ncchkioi_sync_ocnt_reduce (ncchkp, varp->nchunkrec, ocnt, ocnt_all, NULL);
		CHK_ERR
	}

	ncchkioi_assign_chunk_owner (ncchkp, varp, ocnt_all);

	ncchkioi_write_chunk_ocnt (ncchkp, varp, ocnt, sizeof (ncchkioi_chunk_overlap_t));

	NCI_Free (ocnt);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_COWN)

err_out:;
	return err;
}

static inline int ncchkioi_reduce_max_csize_n (
	NC_chk *ncchkp, int nvar, NC_chk_var **varps, MPI_Offset **ocnts, int **cowns) {
	int err = NC_NOERR;
	int i, j, k, v;
	int nchunk;
	MPI_Offset **ocnts_all[2];
	MPI_Offset *cown_size;
	MPI_Offset *ocnt, **ocnt_all;
	int *cown;
	NC_chk_var *varp;
	MPI_Request req;
	MPI_Request *bcast_reqs;
	MPI_Status stat;

	bcast_reqs = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * nvar);
	CHK_PTR (bcast_reqs)

	if (ncchkp->rank == 0) {
		// Size owned by each process
		cown_size = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * ncchkp->np);
		CHK_PTR (cown_size)
		memset (cown_size, 0, sizeof (MPI_Offset) * ncchkp->np);

		// Max #chunks across vars
		nchunk = 0;
		for (v = 0; v < nvar; v++) {
			varp = varps[v];
			if (varp->nchunkrec > nchunk) { nchunk = varp->nchunkrec; }
		}
		// Allocate 2 set of ocnts_all
		ocnts_all[0] = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *) * ncchkp->np * 2);
		CHK_PTR (ocnts_all[0])
		ocnts_all[1] = ocnts_all[0] + ncchkp->np;

		ocnts_all[0][0] = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * nchunk * ncchkp->np * 2);
		CHK_PTR (ocnts_all[0][0])
		ocnts_all[1][0] = ocnts_all[0][0] + nchunk * ncchkp->np;
		for (i = 1; i < ncchkp->np; i++) {
			ocnts_all[0][i] = ocnts_all[0][i - 1] + nchunk;
			ocnts_all[1][i] = ocnts_all[1][i - 1] + nchunk;
		}

		err = MPI_Igather (ocnt, varp->nchunkrec, MPI_LONG_LONG, ocnt_all, varp->nchunkrec,
						   MPI_LONG_LONG, 0, ncchkp->comm, &req);
		CHK_ERR

		for (v = 0; v < nvar; v++) {
			cown	 = cowns[v];
			varp	 = varps[v];
			ocnt	 = ocnts[v];
			ocnt_all = ocnts_all[v & 1];

			// Wait for comm
			err = MPI_Wait (&req, &stat);
			CHK_ERR

			// Post comm for next var
			if (v < nvar - 1) {
				err = MPI_Igather (ocnts[v + 1], varps[v + 1]->nchunkrec, MPI_LONG_LONG,
								   ocnts_all[(v + 1) & 1], varps[v + 1]->nchunkrec, MPI_LONG_LONG,
								   0, ncchkp->comm, &req);
				CHK_ERR
			}

			// Compute max rank for this var
			memset (cown, 0, sizeof (int) * varp->nchunkrec);
			for (i = 0; i < varp->nchunk; i++) {
				k = 0;
				for (j = 1; j < ncchkp->np; j++) {
					if (ocnt_all[j][i] - cown_size[j] > ocnt_all[k][i] - cown_size[k]) { k = j; }
				}
				cown_size[k] +=
					(MPI_Offset) ((double)(varp->chunksize) * ncchkp->cown_ratio) * varp->nrec;
				cown[i] = k;
			}

			// Bcast result
			err = MPI_Ibcast (cown, varp->nchunkrec, MPI_INT, 0, ncchkp->comm, bcast_reqs + v);
			CHK_ERR
		}
	} else {
		for (v = 0; v < nvar; v++) {
			// Send to rank 0
			err = MPI_Gather (ocnts[v], varps[v]->nchunkrec, MPI_LONG_LONG, NULL,
							  varps[v]->nchunkrec, MPI_LONG_LONG, 0, ncchkp->comm);
			CHK_ERR
			// Recv result
			err = MPI_Ibcast (cowns[v], varp->nchunkrec, MPI_INT, 0, ncchkp->comm, bcast_reqs + v);
			CHK_ERR
		}
	}

	err = MPI_Waitall (nvar, bcast_reqs, MPI_STATUS_IGNORE);
	CHK_ERR

	if (ncchkp->rank == 0) {
		NCI_Free (cown_size);
		NCI_Free (ocnts_all[0][0]);
		NCI_Free (ocnts_all[0]);
	}
	NCI_Free (bcast_reqs);

err_out:;
	return err;
}

int ncchkioi_calc_chunk_owner_gather (
	NC_chk *ncchkp, int nvar, NC_chk_var **varps, int nput, int *putreqs, int nget, int *getreqs) {
	int err = NC_NOERR;
	int i, j, k;
	int cid;  // Chunk iterator
	int req;
	int max_ndim;
	int nchunks;
	MPI_Offset overlapsize;
	MPI_Offset *ostart, *osize;
	MPI_Offset *citr;  // Bounding box for chunks overlapping my own write region
	MPI_Offset **ocnts;
	int **cowns;
	NC_chk_var *varp;
	int *idmap;
	NC_chk_req *reqp;

	NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_COWN)

	// Allocate buffer for overlappinp structure
	// Box of single overlap
	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * ncchkp->max_ndim * 3);
	osize  = ostart + ncchkp->max_ndim;
	citr   = osize + ncchkp->max_ndim;
	// Calculate total number of chunks to assign
	idmap	= (int *)NCI_Malloc (sizeof (int) * ncchkp->vars.cnt);
	nchunks = 0;
	for (i = 0; i < nvar; i++) {
		idmap[varps[i]->varid] = i;
		nchunks += varps[i]->nchunkrec;
	}
	// Overlap count struct
	ocnts	 = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *) * nvar);
	ocnts[0] = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * nchunks);
	cowns	 = (int **)NCI_Malloc (sizeof (int *) * nvar);
	cowns[0] = varps[0]->chunk_owner;
	for (i = 1; i < nvar; i++) {
		ocnts[i] = ocnts[i - 1] + varps[i - 1]->nchunkrec;
		cowns[i] = varps[i]->chunk_owner;
	}

	// Count overlapsize for each request
	memset (ocnts[0], 0, sizeof (ncchkioi_chunk_overlap_t) * nchunks);
	for (i = 0; i < nput; i++) {
		reqp = ncchkp->putlist.reqs + putreqs[i];
		ncchkioi_rec_chunk_overlap (ostart, osize, citr, ncchkp->vars.data + reqp->varid,
									ocnts[idmap[reqp->varid]], reqp);
	}
	for (i = 0; i < nget; i++) {
		reqp = ncchkp->getlist.reqs + getreqs[i];
		ncchkioi_rec_chunk_overlap (ostart, osize, citr, ncchkp->vars.data + reqp->varid,
									ocnts[idmap[reqp->varid]], reqp);
	}

	// Calculate the max rank
	ncchkioi_reduce_max_csize_n (ncchkp, nvar, varps, ocnts, cowns);

	// Copy owner to other records
	for (i = 0; i < nvar; i++) {
		varp = varps[i];
		if (varp->isrec) {
			for (j = varp->nchunkrec; j < varp->nchunk; j += varp->nchunkrec) {
				memcpy (varp->chunk_owner + j, varp->chunk_owner, sizeof (int) * varp->nchunkrec);
			}
		}
		ncchkioi_write_chunk_ocnt (ncchkp, varp, ocnts[i], sizeof (MPI_Offset));
	}

	NCI_Free (ostart);
	NCI_Free (ocnts[0]);
	NCI_Free (ocnts);
	NCI_Free (cowns);
	NCI_Free (idmap);

err_out:;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_COWN)

	return err;
}