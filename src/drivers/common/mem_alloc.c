/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* These are routines for allocating and deallocating heap memory dynamically.
   They should be called as
       NCI_Malloc(size)
       NCI_Calloc(nelems, esize)
       NCI_Realloc(ptr, size)
       NCI_Free(ptr)

   In macro.h, they are macro-replaced to
       NCI_Malloc_fn(size, __LINE__, __FILE__) and
       NCI_Calloc_fn(nelems, esize, __LINE__, __FILE__) and
       NCI_Realloc_fn(ptr, size, __LINE__, __func__, __FILE__)
       NCI_Free_fn(ptr,__LINE__,__FILE__).
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h> /* strcpy(), strlen() */
#include <errno.h>  /* EINVAL */
#include <assert.h>

#include <pnetcdf.h>
#include <pnc_debug.h>
#include <common.h>

/* PNC_MALLOC_TRACE is set at the configure time when --enable-debug is used */
#ifdef PNC_MALLOC_TRACE

#ifdef HAVE_SEARCH_H
#include <search.h> /* tfind(), tsearch() and tdelete() */
#endif

/* static variables for malloc tracing (initialized to 0s) */
static void   *ncmpii_mem_root;
static size_t  ncmpii_mem_alloc;
static size_t  ncmpii_max_mem_alloc;

#ifdef ENABLE_THREAD_SAFE
/* updating the binary tree used in tfind()/tsearch()/tdelete() is not
 * thread-safe, protect these subroutines with a mutex */
#include<pthread.h>
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif

#if 0
/*----< ncmpii_init_malloc_tracing() >----------------------------------------*/
void ncmpii_init_malloc_tracing(void)
{
    ncmpii_mem_alloc     = 0;
    ncmpii_max_mem_alloc = 0;
    ncmpii_mem_root      = NULL;
}
#endif

typedef struct {
    void   *self;
    void   *buf;
    size_t  size;
    int     lineno;
    char   *func;
    char   *filename;
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
    if (which == preorder || which == leaf)
        fprintf(stdout, "Warning: malloc yet to be freed (buf=%p size=%zd filename=%s func=%s line=%d)\n", f->buf, f->size, f->filename, f->func, f->lineno);
}

/*----< ncmpii_add_mem_entry() >----------------------------------------------*/
/* add a new malloc entry to the table */
static
void ncmpii_add_mem_entry(void       *buf,
                          size_t      size,
                          const int   lineno,
                          const char *func,
                          const char *filename)
{
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif

    /* use C tsearch utility */
    ncmpii_mem_entry *node = (ncmpii_mem_entry*) malloc(sizeof(ncmpii_mem_entry));
    node->self     = node;
    node->buf      = buf;
    node->size     = size;
    node->lineno   = lineno;
    node->func     = (char*)malloc(strlen(func)+1);
    node->filename = (char*)malloc(strlen(filename)+1);
    strcpy(node->func, func);
    node->func[strlen(func)] = '\0';
    strcpy(node->filename, filename);
    node->filename[strlen(filename)] = '\0';

    /* search and add a new item */
    void *ret = tsearch(node, &ncmpii_mem_root, ncmpii_cmp);
    if (ret == NULL) {
        fprintf(stderr, "Error at line %d file %s: tsearch()\n",
                __LINE__,__FILE__);
    }
    else {
        ncmpii_mem_alloc += size;
        ncmpii_max_mem_alloc = MAX(ncmpii_max_mem_alloc, ncmpii_mem_alloc);
    }
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif
}

/*----< ncmpii_del_mem_entry() >---------------------------------------------*/
/* delete a malloc entry from the table */
static
void ncmpii_del_mem_entry(void *buf)
{
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif
    /* use C tsearch utility */
    if (ncmpii_mem_root != NULL) {
        ncmpii_mem_entry node;
        node.buf  = buf;
        void *ret = tfind(&node, &ncmpii_mem_root, ncmpii_cmp);
        ncmpii_mem_entry **found = (ncmpii_mem_entry**) ret;
        if (ret == NULL) {
            fprintf(stderr, "Error at line %d file %s: tfind() buf=%p\n",
                    __LINE__,__FILE__,buf);
            goto fn_exit;
        }
        /* free space for func and filename */
        free((*found)->func);
        free((*found)->filename);

        /* subtract the space amount to be freed */
        ncmpii_mem_alloc -= (*found)->size;
        void *tmp = (*found)->self;
        ret = tdelete(&node, &ncmpii_mem_root, ncmpii_cmp);
        if (ret == NULL) {
            fprintf(stderr, "Error at line %d file %s: tdelete() buf=%p\n",
                    __LINE__,__FILE__,buf);
            goto fn_exit;
        }
        free(tmp);
    }
    else
        fprintf(stderr, "Error at line %d file %s: ncmpii_mem_root is NULL\n",
                __LINE__,__FILE__);
fn_exit:
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif
    return;
}
#endif

/*----< NCI_Malloc_fn() >-----------------------------------------------------*/
/* This subroutine is esstentially the same as calling malloc().
 * According to malloc man page, If size is 0, then malloc() returns either
 * NULL, or a unique pointer value that can later be successfully passed to
 * free(). Thus, there is no need to check whether size is zero.
 */
void *NCI_Malloc_fn(size_t      size,
                    const int   lineno,
                    const char *func,
                    const char *filename)
{
    void *buf = malloc(size);
#ifdef PNETCDF_DEBUG
    if (size > 0 && buf == NULL)
        fprintf(stderr, "malloc(%zd) failed in file %s func %s line %d\n",
                size, filename, func, lineno);
#endif
#ifdef PNC_MALLOC_TRACE
    ncmpii_add_mem_entry(buf, size, lineno, func, filename);
#endif
    return buf;
}


/*----< NCI_Calloc_fn() >-----------------------------------------------------*/
/* This subroutine is esstentially the same as calling calloc().
 * According to calloc man page, If nelem is 0, then calloc() returns either
 * NULL, or a unique pointer value that can later be successfully passed to
 * free(). Thus, there is no need to check whether nelem is zero.
 */
void *NCI_Calloc_fn(size_t      nelem,
                    size_t      elsize,
                    const int   lineno,
                    const char *func,
                    const char *filename)
{
    void *buf = calloc(nelem, elsize);
#ifdef PNETCDF_DEBUG
    if (nelem > 0 && buf == NULL)
        fprintf(stderr, "calloc(%zd, %zd) failed in file %s func %s line %d\n",
                nelem, elsize, filename, func, lineno);
#endif
#ifdef PNC_MALLOC_TRACE
    ncmpii_add_mem_entry(buf, nelem * elsize, lineno, func, filename);
#endif
    return buf;
}


/*----< NCI_Realloc_fn() >----------------------------------------------------*/
/* This subroutine is esstentially the same as calling realloc().
 * According to realloc man page, if ptr is NULL, then the call is equivalent
 * to malloc(size), for all values of size; if size is equal to zero, and ptr
 * is not NULL, then the call is equivalent to free(ptr). Unless ptr is NULL,
 * it must have been returned by an earlier call to malloc(), calloc() or
 * realloc(). If size was equal to 0, either NULL or a pointer suitable to be
 * passed to free() is returned.
 */
void *NCI_Realloc_fn(void       *ptr,
                     size_t      size,
                     const int   lineno,
                     const char *func,
                     const char *filename)
{
    if (ptr == NULL) return NCI_Malloc_fn(size, lineno, func, filename);

    if (size == 0) {
        NCI_Free_fn(ptr, lineno, func, filename);
        return NCI_Malloc_fn(0, lineno, func, filename);
    }

#ifdef PNC_MALLOC_TRACE
    ncmpii_del_mem_entry(ptr);
#endif
    void *buf = (void *) realloc(ptr, size);
#ifdef PNETCDF_DEBUG
    if (buf == NULL)
	fprintf(stderr, "realloc failed in file %s func %s line %d\n",
                filename, func, lineno);
#endif
#ifdef PNC_MALLOC_TRACE
    ncmpii_add_mem_entry(buf, size, lineno, func, filename);
#endif
    return buf;
}


/*----< NCI_Free_fn() >-------------------------------------------------------*/
/* This subroutine is esstentially the same as calling free().
 * According to free man page, free() frees the memory space pointed to by ptr,
 * which must have been returned by a previous call to malloc(), calloc() or
 * realloc(). Otherwise, or if free(ptr) has already been called before,
 * undefined  behavior occurs. If ptr is NULL, no operation is performed.
 */
void NCI_Free_fn(void       *ptr,
                 const int   lineno,
                 const char *func,
                 const char *filename)
{
    if (ptr == NULL) return;
#ifdef PNC_MALLOC_TRACE
    ncmpii_del_mem_entry(ptr);
#endif
    free(ptr);
}


/*----< ncmpi_inq_malloc_size() >--------------------------------------------*/
/* This is an independent subroutine.
 * report the current aggregate size allocated by malloc, yet to be freed */
int ncmpi_inq_malloc_size(MPI_Offset *size)
{
#ifdef PNC_MALLOC_TRACE
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif
    *size = (MPI_Offset)ncmpii_mem_alloc;
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif
    return NC_NOERR;
#else
    DEBUG_RETURN_ERROR(NC_ENOTENABLED)
#endif
}

/*----< ncmpi_inq_malloc_max_size() >----------------------------------------*/
/* This is an independent subroutine.
 * get the max watermark ever researched by malloc (aggregated amount) */
int ncmpi_inq_malloc_max_size(MPI_Offset *size)
{
#ifdef PNC_MALLOC_TRACE
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif
    *size = (MPI_Offset)ncmpii_max_mem_alloc;
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif
    return NC_NOERR;
#else
    DEBUG_RETURN_ERROR(NC_ENOTENABLED)
#endif
}

/*----< ncmpi_inq_malloc_list() >--------------------------------------------*/
/* This is an independent subroutine.
 * walk the malloc tree and print yet-to-be-freed malloc residues */
int ncmpi_inq_malloc_list(void)
{
#ifdef PNC_MALLOC_TRACE
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_lock(&lock);
#endif
    /* check if malloc tree is empty */
    if (ncmpii_mem_root != NULL)
        twalk(ncmpii_mem_root, walker);
#ifdef ENABLE_THREAD_SAFE
    pthread_mutex_unlock(&lock);
#endif
    return NC_NOERR;
#else
    DEBUG_RETURN_ERROR(NC_ENOTENABLED)
#endif
}

