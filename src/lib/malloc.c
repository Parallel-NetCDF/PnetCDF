/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* These are routines for allocating and deallocating memory.
   They should be called as NCI_Malloc(size) and
   NCI_Free(ptr). In macro.h, they are macro-replaced to 
   NCI_Malloc_fn(size,__LINE__,__FILE__) and 
   NCI_Free_fn(ptr,__LINE__,__FILE__).
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <mpi.h>

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
	fprintf(stderr, "Out of memory in file %s, line %d\n", fname, lineno);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
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
	fprintf(stderr, "Out of memory in file %s, line %d\n", fname, lineno);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return new;
}


void *NCI_Realloc_fn(void *ptr, size_t size, int lineno, const char *fname)
{
#ifdef NC_DEBUG
    if (size == 0)
        fprintf(stderr, "Attempt to realloc zero-size in file %s, line %d\n", fname, lineno);
#endif
    if (size == 0) return NULL;
    void *new = (void *) realloc(ptr, size);
    if (!new) {
	fprintf(stderr, "realloc failed in file %s, line %d\n", fname, lineno);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return new;
}


void NCI_Free_fn(void *ptr, int lineno, const char *fname)
{
#ifdef NC_DEBUG
    if (ptr == NULL)
	fprintf(stderr, "Attempt to free null pointer in file %s, line %d\n", fname, lineno);
#endif
    if (ptr == NULL) return;
    free(ptr);
}


