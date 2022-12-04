/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <common.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <ncchkio_driver.h>
#include <pnc_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ncchkio_internal.h"

MPI_Offset gcd (MPI_Offset a, MPI_Offset b) {
	if (b) {
		while ((a %= b) && (b %= a))
			;
	}
	return a + b;
}

void gcd_reduce (long long *in, long long *inout, int *len, MPI_Datatype *dptr) {
	int i;

	for (i = 0; i < *len; i++) {
		if (*inout)
			while (((*in) %= (*inout)) && ((*inout) %= (*in)))
				;
		(*inout) = (*inout) + (*in);
		in++;
		inout++;
	}
}

int smaller (const void *a, const void *b) { return (*(MPI_Offset *)b - *(MPI_Offset *)a); }

int ncchkioi_calc_chunk_size (
	NC_chk *ncchkp, NC_chk_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts) {
	int err = NC_NOERR;
	int r, i, j;
	int primes[] = {2,	3,	5,	7,	11, 13, 17, 19, 23, 29, 31, 37, 41,
					43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
	MPI_Offset *chunkdim;
	MPI_Offset **candidates;
	MPI_Offset chunksize;
	MPI_Offset ub, lb;
	MPI_Op gcd_op;

	NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_CSIZE)

	// Upper and lower bound of reasonable chunk size
	ub = (MPI_Offset)INT_MAX;  // Max chunk size supported
	lb = 1;
	for (i = 0; i < varp->ndim; i++) { lb *= varp->dimsize[i]; }
	lb /= (MPI_Offset)INT_MAX;	// Max # chunks supported
	if (lb < varp->ndim * 3) {	// Metadata should not exceed data
		lb = varp->ndim * 3;
	}
	if (lb < 1024) {  // At least 1 KiB for efficiency
		lb = 1024;
	}

	/* Infer chunk size by reqs
	 * Assume the application is doing blocked division
	 * If we set chunk dim to gcd of all access boundary, no communication required
	 * If the pattern is completely randomized, the result will likely be 1
	 */
	chunkdim = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim);
	if (nreq > 0) {
		candidates	  = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *) * varp->ndim);
		candidates[0] = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * nreq);
		for (i = 1; i < varp->ndim; i++) { candidates[i] = candidates[i - 1] + nreq; }
		for (r = 0; r < nreq; r++) {
			for (i = 0; i < varp->ndim; i++) {
				candidates[i][r] = gcd (starts[r][i], counts[r][i]);
			}
		}
		for (i = 0; i < varp->ndim; i++) {
			qsort (candidates[i], nreq, sizeof (MPI_Offset), smaller);
			chunkdim[i] = candidates[i][0];
			for (r = 1; r < nreq / 2; r++) {  // Take the top 50% to drop out fragment writes
				chunkdim[i] = gcd (chunkdim[i], candidates[i][r]);
			}
		}
	} else {
		for (i = 0; i < varp->ndim; i++) {
			chunkdim[i] = 0;  // We have no clue, listen to other processes
		}
	}

	// Global gcd
	MPI_Op_create ((MPI_User_function *)gcd_reduce, 1, &gcd_op);
	CHK_ERR_ALLREDUCE (MPI_IN_PLACE, chunkdim, varp->ndim, MPI_LONG_LONG, gcd_op, ncchkp->comm);
	MPI_Op_free (&gcd_op);

	// If we have no clue accross processes, set chunk to max
	for (i = 0; i < varp->ndim; i++) {
		if (chunkdim[i] == 0) { chunkdim[i] = varp->dimsize[i]; }
	}

	// At least 1 for rec dim
	if (varp->isrec) {
		if (chunkdim[0] == 0) { chunkdim[0] = 1; }
	}

	// Check if chunk size is resonable (not too large or too small)
	chunksize = 1;
	for (i = 0; i < varp->ndim; i++) { chunksize *= chunkdim[i]; }

	// we only support chunk size up to INT_MAX
	if (chunksize > ub) {
		// Can we find perffect split using small prime numbers?
		j = 0;
		while ((j < 25) && (chunksize > ub)) {
			r = 1;
			for (i = 0; i < varp->ndim; i++) {	// Spliting chunks along dims
				if (chunkdim[i] % primes[j] == 0) {
					chunkdim[i] /= primes[j];
					chunksize /= primes[j];
					r = 0;
				}
			}
			if (r) {  // No fit, try next prime
				j++;
			}
		}
		if (j >= 25) {	// If not, we still need to split even we need to introduce communication
						// overhead
			for (i = 0; chunksize > ub; i++) {	// Merging chunks
				chunkdim[i % varp->ndim] /= 2;
				chunksize /= 2;
			}
		}
	} else if (chunksize < lb) {  // Data smaller than metadata
		int tmp;
		int *heap;
		int hsize;

		// Build heap of smallest chunk dim
		heap = (int *)NCI_Malloc (sizeof (int) * varp->ndim);
		for (i = 0; i < varp->ndim; i++) {
			heap[i] = i;
			j		= i;
			r		= (j - 1) / 2;
			while (j > 0 && chunkdim[heap[j]] < chunkdim[heap[r]]) {
				tmp - heap[j];
				heap[j] = heap[r];
				heap[r] = tmp;
				j		= r;
				r		= (j - 1) / 2;
			}
		}

		hsize = varp->ndim;
		while (chunksize < lb && hsize > 0) {
			j = heap[0];
			if (chunkdim[j] * 2 <= varp->dimsize[j]) {	// Merge chunk along smallest dim
				chunkdim[j] *= 2;
				chunksize *= 2;
			} else {  // Already reach var dim, remove from consideration
				heap[0] = heap[--hsize];
			}
			// Heapify
			r = 0;
			i = r * 2 + 1;
			j = r * 2 + 2;
			while (i < hsize) {
				if ((j >= hsize) || (chunkdim[heap[i]] < chunkdim[heap[j]])) {
					if (chunkdim[heap[i]] < chunkdim[heap[r]]) {
						tmp		= heap[r];
						heap[r] = heap[i];
						heap[i] = tmp;
						r		= i;
					} else {
						break;
					}
				} else {
					if (chunkdim[heap[j]] < chunkdim[heap[r]]) {
						tmp		= heap[r];
						heap[r] = heap[j];
						heap[j] = tmp;
						r		= j;
					} else {
						break;
					}
				}
				i = r * 2 + 1;
				j = r * 2 + 2;
			}
		}
		NCI_Free (heap);

		// Still not enough after doing everything, just set to entire var
		if (chunksize < lb) {
			memcpy (chunkdim, varp->dimsize, sizeof (MPI_Offset) * varp->ndim);

			// At least 1 for rec dim
			if (varp->isrec) {
				if (chunkdim[0] == 0) { chunkdim[0] = 1; }
			}
		}
	}

	for (i = 0; i < varp->ndim; i++) { varp->chunkdim[i] = (int)chunkdim[i]; }

err_out:;

	NCI_Free (chunkdim);
	if (nreq > 0) {
		NCI_Free (candidates[0]);
		NCI_Free (candidates);
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_CSIZE)

	return err;
}
