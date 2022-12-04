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

static int ncchkioi_cache_evict (NC_chk *ncchkp) {
	int err=NC_NOERR;
	NC_chk_cache *target;

	target = ncchkp->cache_head;

	if (target == NULL || target->serial >= ncchkp->cache_serial) {
		printf ("Rank %d: Cache limit exceeded\n", ncchkp->rank);
		RET_ERR(NC_ENOMEM)
	}

	// Remove from list
	ncchkp->cache_head = target->next;
	if (ncchkp->cache_tail == target) { ncchkp->cache_tail = NULL; }

	ncchkp->cache_used -= target->bsize;  // Return budget
	ncchkp->cache_head = target->next;

	*(target->ref) = NULL;	// Mark as evicted
	NCI_Free (target->buf);
	NCI_Free (target);

err_out:;
	return err;
}

int ncchkioi_cache_alloc (NC_chk *ncchkp, MPI_Offset size, NC_chk_cache **ref) {
	int err = NC_NOERR;
	NC_chk_cache *target;

	// Evict cached data if no space
	if (ncchkp->cache_limit > 0) {
		while (ncchkp->cache_used + size > ncchkp->cache_limit) {
			err = ncchkioi_cache_evict (ncchkp);
			CHK_ERR
		}
	}
	ncchkp->cache_used += size;

	// Prepare cache entry
	target = (NC_chk_cache *)NCI_Malloc (sizeof (NC_chk_cache));
	if (target == NULL) { DEBUG_RETURN_ERROR (NC_ENOMEM) }
	target->bsize  = size;
	target->next   = NULL;
	target->prev   = ncchkp->cache_tail;
	target->ref	   = ref;
	target->serial = ncchkp->cache_serial;
	target->buf	   = NCI_Malloc (size);
#ifdef PNETCDF_DEBUG
	memset (target->buf, 0, size);
#endif

	// Insert to list tail
	if (ncchkp->cache_tail != NULL) {
		ncchkp->cache_tail->next = target;
	} else {
		ncchkp->cache_head = target;
	}
	ncchkp->cache_tail = target;

	// Assign reference
	*ref = target;

err_out:;
	return err;
}

void ncchkioi_cache_visit (NC_chk *ncchkp, NC_chk_cache *target) {
	if (target != ncchkp->cache_tail) {
		// Remove from list
		if (target->prev != NULL) { target->prev->next = target->next; }
		if (target->next != NULL) { target->next->prev = target->prev; }

		// Insert to list tail
		target->next			 = NULL;
		target->prev			 = ncchkp->cache_tail;
		ncchkp->cache_tail->next = target;
		ncchkp->cache_tail		 = target;
	}
}

void ncchkioi_cache_free (NC_chk *ncchkp) {
	NC_chk_cache *pre, *cur;

	cur = ncchkp->cache_head;
	while (cur != NULL) {
		pre = cur;
		cur = cur->next;
		NCI_Free (pre->buf);
		NCI_Free (pre);
	}
}
