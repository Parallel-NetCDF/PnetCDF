/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* These are routines for allocating and deallocating memory.
   They should be called as NCI_Malloc(size) and
   NCI_Free(ptr). In macro.h, they are macro-replaced to 
   NCI_Malloc_fn(size,__LINE__,__FILE__) and 
   NCI_Free_fn(ptr,__LINE__,__FILE__).
 */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#ifdef NC_TRACK_MALLOC
#include <search.h> /* tsearch() and tdelete() */
#endif

#include <mpi.h>

/* global variables for malloc tracking */
static void       *ncmpii_mem_root;
static MPI_Offset  ncmpii_mem_alloc;
static MPI_Offset  ncmpii_max_mem_alloc;

/*----< ncmpii_init_malloc_tracking() >---------------------------------------*/
void ncmpii_init_malloc_tracking(void)
{
    ncmpii_mem_alloc     = 0;
    ncmpii_max_mem_alloc = 0;
    ncmpii_mem_root      = NULL;
}

/*----< ncmpii_inq_malloc_size() >--------------------------------------------*/
/* get the current aggregate size allocated by malloc */
int ncmpii_inq_malloc_size(MPI_Offset *size)
{
#ifdef NC_TRACK_MALLOC
    *size = ncmpii_mem_alloc;
    return 1;
#else
    return 0;
#endif
}

/*----< ncmpii_inq_max_malloc_size() >----------------------------------------*/
/* get the max watermark ever researched by malloc */
int ncmpii_inq_max_malloc_size(MPI_Offset *size)
{
#ifdef NC_TRACK_MALLOC
    *size = ncmpii_max_mem_alloc;
    return 1;
#else
    return 0;
#endif
}

#ifdef NC_TRACK_MALLOC
#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

typedef struct {
    void       *self;
    void       *buf;
    MPI_Offset  size;
} ncmpii_mem_entry;

static
int ncmpii_cmp(const void *a, const void *b) {
    ncmpii_mem_entry *fa = (ncmpii_mem_entry*)a;
    ncmpii_mem_entry *fb = (ncmpii_mem_entry*)b;

    if (fa->buf > fb->buf) return  1;
    if (fa->buf < fb->buf) return -1;
    return 0;
}

static
void walker(const void *node, const VISIT which, const int depth) {
    ncmpii_mem_entry *f;
    f = *(ncmpii_mem_entry **)node;
    printf("walker: malloc residue buf=%p size=%lld\n", f->buf, f->size);
}
#endif

/* check if there is any malloc residue */
int ncmpii_inq_malloc_walk(void)
{
#ifdef NC_TRACK_MALLOC
    /* check if malloc tree is empty */
    if (ncmpii_mem_root != NULL)
        twalk(ncmpii_mem_root, walker);
    return 1;
#else
    return 0;
#endif
}

#ifdef NC_TRACK_MALLOC
/*----< ncmpii_add_mem_entry() >----------------------------------------------*/
/* add a new malloc entry to the table */
static
void ncmpii_add_mem_entry(void   *buf,
                          MPI_Offset  size)
{
    /* use C tsearch utility */
    ncmpii_mem_entry *node = (ncmpii_mem_entry*) malloc(sizeof(ncmpii_mem_entry));
    node->self = node;
    node->buf  = buf;
    node->size = size;
    /* search and add a new item */
    void *ret = tsearch(node, &ncmpii_mem_root, ncmpii_cmp);
    if (ret == NULL) {
        printf("Error: tsearch()\n");
        return;
    }
    ncmpii_mem_alloc += size;
    ncmpii_max_mem_alloc = MAX(ncmpii_max_mem_alloc, ncmpii_mem_alloc);
}

/*----< ncmpii_del_mem_entry() >----------------------------------------------*/
/* delete a malloc entry from the table */
static
void ncmpii_del_mem_entry(void *buf)
{
    /* use C tsearch utility */
    if (ncmpii_mem_root != NULL) {
        ncmpii_mem_entry node;
        node.buf  = buf;
        ncmpii_mem_entry **ret = tfind(&node, &ncmpii_mem_root, ncmpii_cmp);
        if (ret == NULL) {
            printf("Error: tdelete() buf=%p\n", buf);
            return;
        }
        ncmpii_mem_alloc -= (*ret)->size;
        void *tmp = (*ret)->self;
        ret = tdelete(&node, &ncmpii_mem_root, ncmpii_cmp);
        if (ret == NULL) {
            printf("Error: tdelete() buf=%p\n", buf);
            return;
        }
        free(tmp);
    }
}
#endif

void *NCI_Malloc_fn(size_t size, int lineno, const char *fname);
void *NCI_Calloc_fn(size_t nelem, size_t elsize, int lineno, const char *fname);
void *NCI_Realloc_fn(void *ptr, size_t size, int lineno, const char *fname);
void  NCI_Free_fn(void *ptr, int lineno, const char *fname);


void *NCI_Malloc_fn(size_t size, int lineno, const char *fname)
{
#ifdef NC_DEBUG
    if (size == 0)
        fprintf(stderr, "Attempt to malloc zero-size in file %s, line %d\n", fname, lineno);
#endif
    if (size == 0) return NULL;
    void *new = malloc(size);
    if (!new) {
	fprintf(stderr, "malloc(%zd) failed in file %s, line %d\n", size, fname, lineno);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
#ifdef NC_TRACK_MALLOC
    ncmpii_add_mem_entry(new, size);
#endif
    return new;
}


void *NCI_Calloc_fn(size_t nelem, size_t elsize, int lineno, const char *fname)
{
#ifdef NC_DEBUG
    if (nelem == 0 || elsize == 0)
        fprintf(stderr, "Attempt to calloc zero-size in file %s, line %d\n", fname, lineno);
#endif
    if (nelem == 0 || elsize == 0) return NULL;
    void *new =calloc(nelem, elsize);
    if (!new) {
	fprintf(stderr, "calloc(%zd, %zd) failed in file %s, line %d\n", nelem, elsize, fname, lineno);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
#ifdef NC_TRACK_MALLOC
    ncmpii_add_mem_entry(new, nelem * elsize);
#endif
    return new;
}


void *NCI_Realloc_fn(void *ptr, size_t size, int lineno, const char *fname)
{
#ifdef NC_DEBUG
    if (size == 0)
        fprintf(stderr, "Attempt to realloc zero-size in file %s, line %d\n", fname, lineno);
#endif
#ifdef NC_TRACK_MALLOC
    if (ptr != NULL) ncmpii_del_mem_entry(ptr);
#endif
    if (size == 0) return NULL;
    void *new = (void *) realloc(ptr, size);
    if (!new) {
	fprintf(stderr, "realloc failed in file %s, line %d\n", fname, lineno);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
#ifdef NC_TRACK_MALLOC
    ncmpii_add_mem_entry(new, size);
#endif
    return new;
}


void NCI_Free_fn(void *ptr, int lineno, const char *fname)
{
#ifdef NC_DEBUG
    if (ptr == NULL)
	fprintf(stderr, "Attempt to free null pointer in file %s, line %d\n", fname, lineno);
#endif
    if (ptr == NULL) return;
#ifdef NC_TRACK_MALLOC
    ncmpii_del_mem_entry(ptr);
#endif
    free(ptr);
}


