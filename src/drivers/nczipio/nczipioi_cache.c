/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>
#include "nczipio_internal.h"

static int nczipioi_cache_evict(NC_zip *nczipp){
    NC_zip_cache *target;

    target = nczipp->cache_head;

    if (target == NULL || target->serial >= nczipp->cache_serial){
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    // Remove from list
    nczipp->cache_head = target->next;
    if (nczipp->cache_tail == target){
        nczipp->cache_tail = NULL;
    }

    nczipp->cache_used -= target->bsize;    // Return budget
    nczipp->cache_head = target->next;

    *(target->ref) = NULL;  // Mark as evicted
    NCI_Free(target->buf);
    NCI_Free(target);
}

int nczipioi_cache_alloc(NC_zip *nczipp, MPI_Offset size, NC_zip_cache **ref){
    int err;
    NC_zip_cache *target;

    // Evict cached data if no space
    if (nczipp->cache_limit > 0){
        while(nczipp->cache_used + size > nczipp->cache_limit){
            err = nczipioi_cache_evict(nczipp);
        }
    }
    nczipp->cache_used += size;

    // Prepare cache entry
    target = (NC_zip_cache*)NCI_Malloc(sizeof(NC_zip_cache));
    if (target == NULL){
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    target->bsize = size;
    target->next = NULL;
    target->prev = nczipp->cache_tail;
    target->ref = ref;
    target->serial = nczipp->cache_serial;
    target->buf = NCI_Malloc(size);
#ifdef PNETCDF_DEBUG
    memset(target->buf, 0, size);
#endif

    // Insert to list tail
    if (nczipp->cache_tail != NULL){
        nczipp->cache_tail->next = target;
    }
    else{
        nczipp->cache_head = target;
    }
    nczipp->cache_tail = target;

    // Assign reference
    *ref = target;
}

void nczipioi_cache_visit(NC_zip *nczipp, NC_zip_cache *target){
    if (target != nczipp->cache_tail){
        // Remove from list
        if (target->prev != NULL){
            target->prev->next = target->next;
        }
        if (target->next != NULL){
            target->next->prev = target->prev;
        }

        // Insert to list tail
        target->next = NULL;
        target->prev = nczipp->cache_tail;
        nczipp->cache_tail->next = target;
        nczipp->cache_tail = target;
    }
}   

void nczipioi_cache_free(NC_zip *nczipp){
    NC_zip_cache *pre, *cur;

    cur = nczipp->cache_head;
    while(cur != NULL){
        pre = cur;
        cur = cur->next;
        NCI_Free(pre->buf);
        NCI_Free(pre);
    }
}
