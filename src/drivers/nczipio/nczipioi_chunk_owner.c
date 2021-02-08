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
#include <nczipio_driver.h>
#include <pnc_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nczipio_internal.h"

typedef struct nczipioi_chunk_overlap_t {
	MPI_Offset osize;
	int rank;
} nczipioi_chunk_overlap_t;

void max_osize_rank_op (nczipioi_chunk_overlap_t *in, nczipioi_chunk_overlap_t *inout, int *len, MPI_Datatype *dptr) {
	int i;

	for (i = 0; i < *len; i++) {
		if (in->osize > inout->osize) {
			inout->osize = in->osize;
			inout->rank	 = in->rank;
		}
		in++;
		inout++;
	}
}

int nczipioi_calc_chunk_owner (
	NC_zip *nczipp, NC_zip_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
	int err;
	int i, j, k;
	int cid;  // Chunk iterator
	int req;
	double noise, noise_step;
	MPI_Offset overlapsize;
	MPI_Offset *ostart, *osize;
	MPI_Offset *citr;  // Bounding box for chunks overlapping my own write region
	nczipioi_chunk_overlap_t *ocnt, *ocnt_all;
	MPI_Datatype overlaptype;

	NC_ZIP_TIMER_START (NC_ZIP_TIMER_INIT_COWN)

    MPI_Type_contiguous(sizeof(nczipioi_chunk_overlap_t),MPI_BYTE,&overlaptype);
    MPI_Type_commit(&overlaptype);

	ostart = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * 3);
	osize  = ostart + varp->ndim;
	citr   = osize + varp->ndim;

	ocnt	 = (nczipioi_chunk_overlap_t *)NCI_Malloc (sizeof (nczipioi_chunk_overlap_t) * varp->nchunkrec * 2);
	ocnt_all = ocnt + varp->nchunkrec;
	memset (ocnt, 0, sizeof (nczipioi_chunk_overlap_t) * varp->nchunkrec);

	// Count overlapsize of each request
	if (varp->isrec) {
		for (req = 0; req < nreq; req++) {
			nczipioi_chunk_itr_init_ex (varp, starts[req], counts[req], citr, &cid, ostart,
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
			} while (nczipioi_chunk_itr_next_ex (varp, starts[req], counts[req], citr, &cid, ostart,
												 osize));
		}
	} else {
		for (req = 0; req < nreq; req++) {
			nczipioi_chunk_itr_init_ex (varp, starts[req], counts[req], citr, &cid, ostart,
										osize);	 // Initialize chunk iterator
			do {
				// Count overlap
				overlapsize = 1;
				for (i = 0; i < varp->ndim; i++) { overlapsize *= osize[i]; }
				ocnt[cid].osize += (double)overlapsize;
				if (ocnt[cid].osize > varp->chunksize) {
					ocnt[cid].osize = (double)varp->chunksize;
				}
			} while (nczipioi_chunk_itr_next_ex (varp, starts[req], counts[req], citr, &cid, ostart,
												 osize));
		}
	}
	for (i = 0; i < varp->nchunkrec; i++) {
		ocnt[i].rank = nczipp->rank;
		ocnt[i].osize -=
			(double)(nczipp->nmychunks >> 4);  // Penality for load ballance, set at 1/16
		ocnt[i].osize <<= 16;
	}
	// Noise to break tie
	j = (nczipp->rank - nczipp->assigned_chunks) % nczipp->np;
	if (j < 0) j += nczipp->np;
	k		   = nczipp->np - 1;  // noise from 0 ~ np-1
	for (i = j; i < varp->nchunkrec; i++) {
		ocnt[i].osize += k;
		k--;
		if (k < 0) { k += nczipp->np; }
	}
	for (i = 0; i < j; i++) {
		ocnt[i].osize += k;
		k--;
		if (k < 0) { k += nczipp->np; }
	}
	nczipp->assigned_chunks += varp->nchunk;

	CHK_ERR_ALLREDUCE (ocnt, ocnt_all, varp->nchunkrec, overlaptype, MPI_MAXLOC, nczipp->comm);
	for (i = 0; i < varp->nchunkrec; i++) { varp->chunk_owner[i] = ocnt_all[i].rank; }
	if (varp->isrec) {
		for (i = varp->nchunkrec; i < varp->nchunk; i += varp->nchunkrec) {
			memcpy (varp->chunk_owner + i, varp->chunk_owner, sizeof (int) * varp->nchunkrec);
		}
	}

	MPI_Type_free(&overlaptype);

#ifdef PNETCDF_DEBUG
#ifdef PNETCDF_PROFILING
	{
		char *pprefix = getenv ("PNETCDF_OWNER_PREFIX");
		if (pprefix != NULL) {
			if (nczipp->rank == 0) {
				MPI_Status stat;
				FILE *pfile;
				char fname[1024], ppath[1024];

				strcpy (fname, nczipp->path);
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
				for (j = 0; j < varp->nchunk; j++) { fprintf (pfile, "%d, ", j); }
				fprintf (pfile, "\nOwner, ");
				for (j = 0; j < varp->nchunk; j++) {
					fprintf (pfile, "%d, ", varp->chunk_owner[j]);
				}
				fprintf (pfile, "\n0, ");
				for (j = 0; j < varp->nchunk; j++) { fprintf (pfile, "%lf, ", ocnt[j].osize); }
				fprintf (pfile, "\n");
				for (i = 1; i < nczipp->np; i++) {
					MPI_Recv (ocnt_all, varp->nchunk, MPI_2INT, i, 0, nczipp->comm, &stat);
					fprintf (pfile, "%d, ", i);
					for (j = 0; j < varp->nchunk; j++) {
						fprintf (pfile, "%lf, ", ocnt_all[j].osize);
					}
					fprintf (pfile, "\n");
				}

				fclose (pfile);
			} else {
				MPI_Send (ocnt, varp->nchunk, MPI_2INT, 0, 0, nczipp->comm);
			}
		}
	}
#endif
#endif

	NCI_Free (ostart);
	NCI_Free (ocnt);

	NC_ZIP_TIMER_STOP (NC_ZIP_TIMER_INIT_COWN)

	return NC_NOERR;
}