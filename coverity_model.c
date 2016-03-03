void* NCI_Malloc_fn(size_t size, const int lineno, const char *func, const char *filename) {
     __coverity_alloc__(size);
}

void* NCI_Calloc_fn(size_t nelem, size_t elsize, const int lineno, const char *func, const char *filename) {
     __coverity_alloc__(nelem*elsize);
}

void* NCI_Realloc_fn(void *ptr, size_t size, const int lineno, const char *func, const char *filename) {
     __coverity_alloc__(size);
}

void NCI_Free_fn(void *ptr, const int lineno, const char *func, const char *filename) {
     __coverity_free__(ptr);
}

int MPI_Abort(int comm, int errorcode) {
    __coverity_panic__();
}
