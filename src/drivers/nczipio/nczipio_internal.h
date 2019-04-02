#ifndef _nczipio_INTERNAL_H
#define _nczipio_INTERNAL_H

#include "nczipio_driver.h"

typedef struct NC_zip_vector{
    int esize;
    int size;
    int nalloc;
    char *data;
}NC_zip_vector;

// File
extern int nczipioi_init(NC_zip*);
extern int nczipioi_parse_var_info(NC_zip*);
extern int nczipioi_var_list_init(NC_zip_var_list*);
extern int nczipioi_var_list_free(NC_zip_var_list*);
extern int nczipioi_var_list_add(NC_zip_var_list*, NC_zip_var); 

// Util
extern int nczipioi_extract_hint(NC_zip*, MPI_Info);
extern int nczipioi_export_hint(NC_zip *nczipp, MPI_Info info);
extern MPI_Offset NC_Type_size(nc_type);

// Var
extern int nczipioi_var_init(NC_zip*, NC_zip_var*, int);
extern int nczipioi_var_init(NC_zip*, NC_zip_var*, int);
extern int nczipioi_load_var(NC_zip*, NC_zip_var*, int, int*);
extern int nczipioi_save_var(NC_zip*, NC_zip_var*);
extern int nczipioi_save_nvar(NC_zip*, int, int*);

// Chunks
extern int get_chunk_idx(NC_zip_var*, int*);
extern int get_chunk_cord(NC_zip_var*, int, int*);
extern int get_chunk_overlap(NC_zip_var*, int*, const MPI_Offset*, const MPI_Offset*, MPI_Offset*, MPI_Offset*);
extern int get_chunk_overlap_str(NC_zip_var*, int*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, int*, int*);
extern int nczipioi_chunk_itr_init(NC_zip_var*, MPI_Offset*, MPI_Offset*, int*, int*, int*);
extern int nczipioi_chunk_itr_init_str(NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*, int*, int*, int*);
extern int nczipioi_chunk_itr_next(NC_zip_var*, MPI_Offset*, MPI_Offset*, int*, int*, int*);
extern int nczipioi_chunk_itr_next_str(NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*, int*, int*, int*);
extern int nczipioi_chunk_itr_init_cord(NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*);
extern int nczipioi_chunk_itr_next_cord(NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*);
extern int get_chunk_overlap_cord(NC_zip_var*, MPI_Offset*, const MPI_Offset*, const MPI_Offset*, MPI_Offset*, MPI_Offset*);
extern int get_chunk_idx_cord(NC_zip_var*, MPI_Offset*);
extern int get_chunk_itr(NC_zip_var*, int, MPI_Offset*);

// Get
//extern int nczipioi_get_var_old(NC_zip*, NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*, void*);
extern int nczipioi_get_var_cb_chunk(NC_zip*, NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*, void*);
extern int nczipioi_get_var_cb_proc(NC_zip*, NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*, void*);
extern int nczipioi_get_varn(NC_zip *, NC_zip_var*, int , MPI_Offset* const *, MPI_Offset* const *, const void*);
extern int nczipioi_get_varn_cb_chunk(NC_zip*, NC_zip_var*, int nreq, MPI_Offset* const*, MPI_Offset* const*, MPI_Offset* const*, void**);
extern int nczipioi_get_varn_cb_proc(NC_zip*, NC_zip_var*, int nreq, MPI_Offset* const*, MPI_Offset* const*, void**);
extern int nczipioi_iget_var(NC_zip*, int, MPI_Offset*, MPI_Offset*, MPI_Offset*, const MPI_Offset*, void*, MPI_Offset, MPI_Datatype, int*);
extern int nczipioi_iget_varn(NC_zip*, int, int, MPI_Offset**, MPI_Offset**, void*, MPI_Offset, MPI_Datatype, int*);

// Put
//extern int nczipioi_put_var_old(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);
extern int nczipioi_put_var(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);
extern int nczipioi_put_var_cb_chunk(NC_zip*, NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*, void*);
extern int nczipioi_put_var_cb_proc(NC_zip*, NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*, void*);
extern int nczipioi_put_varn(NC_zip*, NC_zip_var*, int, MPI_Offset* const*, MPI_Offset* const *, const void*);
extern int nczipioi_put_varn_cb_chunk(NC_zip*, NC_zip_var*, int, MPI_Offset* const*, MPI_Offset* const*, MPI_Offset* const*, void**);
extern int nczipioi_put_varn_cb_proc(NC_zip*, NC_zip_var*, int, MPI_Offset* const*, MPI_Offset* const*, void**);
extern int nczipioi_iput_var(NC_zip*, int, MPI_Offset*, MPI_Offset*, MPI_Offset*, const void*, const void*, int*);
extern int nczipioi_iput_varn(NC_zip*, int, int, MPI_Offset**, MPI_Offset**, const void*, const void*, int*);
extern int nczipioi_iput_cb_chunk(NC_zip*, int, int*, int*);
extern int nczipioi_iput_cb_proc(NC_zip*, int, int*, int*);

// Nonblocking
extern int nczipioi_req_list_init(NC_zip_req_list*);
extern int nczipioi_req_list_free(NC_zip_req_list*);
extern int nczipioi_req_list_add(NC_zip_req_list*, int*);
extern int nczipioi_req_list_remove(NC_zip_req_list*, int);
extern int nczipioi_wait_put_reqs(NC_zip*, int, int*, int*);
extern int nczipioi_wait_get_reqs(NC_zip*, int, int*, int*);
extern int nczipioi_wait(NC_zip*, int, int*, int*, int);

// Vector
extern int nczipioi_vector_init(NC_zip_vector*, int);
extern int nczipioi_vector_init_ex(NC_zip_vector*, int, int);
extern void nczipioi_vector_free(NC_zip_vector*);
extern int nczipioi_vector_append(NC_zip_vector*, void*);
#endif