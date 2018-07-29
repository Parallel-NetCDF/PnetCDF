/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncbbio_driver.h>

#define LOG_BUFFER_SIZE 1024 /* Size of initial metadata buffer */
#define LOG_ARRAY_SIZE 32 /* Size of initial metadata offset list */
#define SIZE_MULTIPLIER 20    /* When metadata buffer is full, we'll NCI_Reallocate it to META_BUFFER_MULTIPLIER times the original size*/

/*
 * Initialize a variable sized buffer
 * IN   bp: buffer structure to be initialized
 */
int ncbbio_log_buffer_init(NC_bb_buffer *bp){
    bp->buffer = NCI_Malloc(LOG_BUFFER_SIZE);
    if (bp->buffer == NULL){
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }
    bp->nalloc = LOG_BUFFER_SIZE;
    bp->nused = 0;
    return NC_NOERR;
}

/*
 * Free the variable sized buffer
 * IN   bp: buffer structure to be freed
 */
void ncbbio_log_buffer_free(NC_bb_buffer * bp){
    NCI_Free(bp->buffer);
}

/*
 * Allocate space in the variable sized buffer
 * This function works as NCI_Malloc in the metadata buffer
 * IN    bp:    buffer structure
 * IN    size:    size required in the buffer
 */
char* ncbbio_log_buffer_alloc(NC_bb_buffer *bp, size_t size) {
    char* ret;

    /* Expand buffer if needed
     * bp->nused is the size currently in use
     * bp->nalloc is the size of internal buffer
     * If the remaining size is less than the required size, we reallocate the buffer
     */
    if (bp->nalloc < bp->nused + size) {
        /*
         * We don't know how large the required size is, loop until we have enough space
         * Must make sure realloc successed before increasing bp->nalloc
         */
        size_t newsize = bp->nalloc;
        while (newsize < bp->nused + size) {
            /* (new size) = (old size) * (META_BUFFER_MULTIPLIER) */
            newsize *= SIZE_MULTIPLIER;
        }
        /* ret is used to temporaryly hold the allocated buffer so we don't lose ncbbp->metadata.buffer if allocation fails */
        ret = (char*)NCI_Realloc(bp->buffer, newsize);
        /* If not enough memory */
        if (ret == NULL) {
            return ret;
        }
        /* Point to the new buffer and update nalloc */
        bp->buffer = ret;
        bp->nalloc = newsize;
    }

    /* Increase used buffer size and return the allocated space */
    ret = (char*)(((char*)bp->buffer) + bp->nused);
    bp->nused += size;

    return ret;
}

/*
 * Initialize log entry array
 * IN   ep: array to be initialized
 */
int ncbbio_log_sizearray_init(NC_bb_sizevector *sp){
    sp->values = (size_t*)NCI_Malloc(LOG_ARRAY_SIZE * sizeof(size_t));
    if (sp->values == NULL){
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }
    sp->nalloc = LOG_ARRAY_SIZE;
    sp->nused = 0;
    return NC_NOERR;
}

/*
 * Free the log entry array
 * IN   ep: array to be freed
 */
void ncbbio_log_sizearray_free(NC_bb_sizevector *sp){
    NCI_Free(sp->values);
}

/*
 * Append entry to array
 * IN    ep:    array structure
 * IN    ent:    entry to be added
 */
int ncbbio_log_sizearray_append(NC_bb_sizevector *sp, size_t size) {
    size_t *ret;

    /* Expand array if needed
     * sp->nused is the size currently in use
     * sp->nalloc is the size of internal buffer
     * If the remaining size is less than the required size, we reallocate the buffer
     */
    if (sp->nalloc < sp->nused + 1) {
        /*
         * Must make sure realloc successed before increasing sp->nalloc
         * (new size) = (old size) * (META_BUFFER_MULTIPLIER)
         */
        size_t newsize = sp->nalloc * SIZE_MULTIPLIER;
        /* ret is used to temporaryly hold the allocated buffer so we don't lose ncbbp->metadata.buffer if allocation fails */
        ret = (size_t*)NCI_Realloc(sp->values, newsize * sizeof(size_t));
        /* If not enough memory */
        if (ret == NULL) {
            return NC_ENOMEM;
        }
        /* Point to the new buffer and update nalloc */
        sp->values = ret;
        sp->nalloc = newsize;
    }

    /* Add entry to tail */
    sp->values[sp->nused++] = size;

    return NC_NOERR;
}

/*
 * Initialize vector
 * IN   vp: vector to be initialized
 */
int ncbbio_log_intvector_init(NC_bb_intvector *vp){
    vp->values = (int*)NCI_Malloc(LOG_ARRAY_SIZE * SIZEOF_INT);
    if (vp->values == NULL){
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }
    vp->nalloc = LOG_ARRAY_SIZE;
    vp->nused = 0;
    return NC_NOERR;
}

/*
 * Free the vector
 * IN   vp: vector to be freed
 */
void ncbbio_log_intvector_free(NC_bb_intvector *vp){
    NCI_Free(vp->values);
}

/*
 * Append entry to vector
 * IN    vp:    vector structure
 * IN    val:   value to be added
 */
int ncbbio_log_intvector_append(NC_bb_intvector *vp, int size) {
    int *ret;

    /* Expand array if needed
     * vp->nused is the size currently in use
     * vp->nalloc is the size of internal buffer
     * If the remaining size is less than the required size, we reallocate the buffer
     */
    if (vp->nalloc < vp->nused + 1) {
        /*
         * Must make sure realloc successed before increasing vp->nalloc
         * (new size) = (old size) * (META_BUFFER_MULTIPLIER)
         */
        size_t newsize = vp->nalloc * SIZE_MULTIPLIER;
        /* ret is used to temporaryly hold the allocated buffer so we don't lose ncbbp->metadata.buffer if allocation fails */
        ret = (int*)NCI_Realloc(vp->values, newsize * SIZEOF_INT);
        /* If not enough memory */
        if (ret == NULL) {
            return NC_ENOMEM;
        }
        /* Point to the new buffer and update nalloc */
        vp->values = ret;
        vp->nalloc = newsize;
    }

    /* Add entry to tail */
    vp->values[vp->nused++] = size;

    return NC_NOERR;
}

int ncbbio_metaidx_init(NC_bb *ncbbp) {
    NC_bb_metadataidx *ip = &(ncbbp->metaidx);

    ip->nalloc = LOG_ARRAY_SIZE;
    ip->nused = 0;
    ip->entries = (NC_bb_metadataptr*)NCI_Malloc(sizeof(NC_bb_metadataptr) * ip->nalloc);

    return NC_NOERR;
}

int ncbbio_metaidx_add(NC_bb *ncbbp, NC_bb_metadataentry *ptr) {
    NC_bb_metadataidx *ip = &(ncbbp->metaidx);
    NC_bb_metadataptr *tmp;

    if (ip->nused == ip->nalloc) {
        ip->nalloc *= SIZE_MULTIPLIER;
        tmp = (NC_bb_metadataptr*)NCI_Realloc(ip->entries, sizeof(NC_bb_metadataptr) * ip->nalloc);
        ip->entries = tmp;
    }

    ip->entries[ip->nused].ptr = ptr;
    ip->entries[ip->nused].valid = 1;
    ip->entries[ip->nused++].reqid = -1;

    return NC_NOERR;
}

int ncbbio_metaidx_free(NC_bb *ncbbp) {
    NC_bb_metadataidx *ip = &(ncbbp->metaidx);

    NCI_Free(ip->entries);

    return NC_NOERR;
}

