#ifndef _nczipio_INTERNAL_H
#define _nczipio_INTERNAL_H

#include "nczipio_driver.h"

#define NC_ZIP_DRIVER_NONE	0
#define NC_ZIP_DRIVER_DUMMY 1
#define NC_ZIP_DRIVER_ZLIB	2
#define NC_ZIP_DRIVER_SZ	3

#define NC_ZIP_DEFAULT_REC_ALLOC 1024
#define NC_ZIP_REC_MULTIPLIER	 2

#define RET_ERR(E)               \
	{                            \
		err = E;                 \
		DEBUG_TRACE_ERROR (err); \
		goto err_out;            \
	}
#define CHK_ERR                  \
	if (err != MPI_SUCCESS) {    \
		DEBUG_TRACE_ERROR (err); \
		goto err_out;            \
	}

#define CHK_PTR(P)               \
	if (!P) {                    \
		err = NC_ENOMEM;         \
		DEBUG_TRACE_ERROR (err); \
		goto err_out;            \
	}

#define CHK_ERR_WAIT(V0, V1)                         \
	err = MPI_Wait (V0, V1);                         \
	if (err != MPI_SUCCESS) {                        \
		err = ncmpii_error_mpi2nc (err, "MPI_Wait"); \
		DEBUG_RETURN_ERROR (err)                     \
	}

#define CHK_ERR_ALLREDUCE(V0, V1, V2, V3, V4, V5)         \
	err = MPI_Allreduce (V0, V1, V2, V3, V4, V5);         \
	if (err != MPI_SUCCESS) {                             \
		err = ncmpii_error_mpi2nc (err, "MPI_Allreduce"); \
		DEBUG_RETURN_ERROR (err)                          \
	}

#define CHK_ERR_IALLREDUCE(V0, V1, V2, V3, V4, V5, V6)     \
	err = MPI_Iallreduce (V0, V1, V2, V3, V4, V5, V6);     \
	if (err != MPI_SUCCESS) {                              \
		err = ncmpii_error_mpi2nc (err, "MPI_Iallreduce"); \
		DEBUG_RETURN_ERROR (err)                           \
	}

#define CHK_ERR_REDUCE(V0, V1, V2, V3, V4, V5, V6)        \
	err = MPI_Reduce (V0, V1, V2, V3, V4, V5, V6);        \
	if (err != MPI_SUCCESS) {                             \
		err = ncmpii_error_mpi2nc (err, "MPI_Allreduce"); \
		DEBUG_RETURN_ERROR (err)                          \
	}

#define CHK_ERR_GATHER(V0, V1, V2, V3, V4, V5, V6, V7) \
	err = MPI_Gather (V0, V1, V2, V3, V4, V5, V6, V7); \
	if (err != MPI_SUCCESS) {                          \
		err = ncmpii_error_mpi2nc (err, "MPI_Gather"); \
		DEBUG_RETURN_ERROR (err)                       \
	}

#define CHK_ERR_PACK(V0, V1, V2, V3, V4, V5, V6)     \
	err = MPI_Pack (V0, V1, V2, V3, V4, V5, V6);     \
	if (err != MPI_SUCCESS) {                        \
		err = ncmpii_error_mpi2nc (err, "MPI_Pack"); \
		DEBUG_RETURN_ERROR (err)                     \
	}
#ifdef PNETCDF_DEBUG
#define CHK_ERR_UNPACK(V0, V1, V2, V3, V4, V5, V6)          \
	{                                                       \
		int esize;                                          \
		MPI_Type_size (V5, &esize);                         \
		if (V1 - *((int *)(V2)) < V4 * esize) { abort (); } \
		err = MPI_Unpack (V0, V1, V2, V3, V4, V5, V6);      \
		if (err != MPI_SUCCESS) {                           \
			err = ncmpii_error_mpi2nc (err, "MPI_Unpack");  \
			DEBUG_RETURN_ERROR (err)                        \
		}                                                   \
	}
#else
#define CHK_ERR_UNPACK(V0, V1, V2, V3, V4, V5, V6)     \
	err = MPI_Unpack (V0, V1, V2, V3, V4, V5, V6);     \
	if (err != MPI_SUCCESS) {                          \
		err = ncmpii_error_mpi2nc (err, "MPI_Unpack"); \
		DEBUG_RETURN_ERROR (err)                       \
	}
#endif
#define CHK_ERR_TYPE_COMMIT(V0)                             \
	err = MPI_Type_commit (V0);                             \
	if (err != MPI_SUCCESS) {                               \
		err = ncmpii_error_mpi2nc (err, "MPI_Type_commit"); \
		DEBUG_RETURN_ERROR (err)                            \
	}
#ifdef PNETCDF_DEBUG
#define CHK_ERR_TYPE_CREATE_SUBARRAY(V0, V1, V2, V3, V4, V5, V6)                               \
	{                                                                                          \
		int d;                                                                                 \
		for (d = 0; d < V0; d++) {                                                             \
			if (V1[d] < V2[d] + V3[d]) {                                                       \
				printf (                                                                       \
					"Error: Subarray outside array at dim %d. size = %d, ssize = %d, start = " \
					"%d\n",                                                                    \
					d, V1[d], V2[d], V3[d]);                                                   \
				abort ();                                                                      \
			}                                                                                  \
			if (V2[d] <= 0) {                                                                  \
				printf ("Error: Subarray size <= 0 at dim %d. ssize = %d\n", d, V2[d]);        \
				abort ();                                                                      \
			}                                                                                  \
		}                                                                                      \
		err = MPI_Type_create_subarray (V0, V1, V2, V3, V4, V5, V6);                           \
		if (err != MPI_SUCCESS) {                                                              \
			err = ncmpii_error_mpi2nc (err, "MPI_Type_create_subarray");                       \
			DEBUG_RETURN_ERROR (err)                                                           \
		}                                                                                      \
	}
#else
#define CHK_ERR_TYPE_CREATE_SUBARRAY(V0, V1, V2, V3, V4, V5, V6)     \
	err = MPI_Type_create_subarray (V0, V1, V2, V3, V4, V5, V6);     \
	if (err != MPI_SUCCESS) {                                        \
		err = ncmpii_error_mpi2nc (err, "MPI_Type_create_subarray"); \
		DEBUG_RETURN_ERROR (err)                                     \
	}
#endif
#define CHK_ERR_WAITALL(V0, V1, V2)                     \
	err = MPI_Waitall (V0, V1, V2);                     \
	if (err != MPI_SUCCESS) {                           \
		err = ncmpii_error_mpi2nc (err, "MPI_Waitall"); \
		DEBUG_RETURN_ERROR (err)                        \
	}
#define CHK_ERR_MPROBE(V0, V1, V2, V3, V4)             \
	err = MPI_Mprobe (V0, V1, V2, V3, V4);             \
	if (err != MPI_SUCCESS) {                          \
		err = ncmpii_error_mpi2nc (err, "MPI_Mprobe"); \
		DEBUG_RETURN_ERROR (err)                       \
	}
#define CHK_ERR_GET_COUNT(V0, V1, V2)                     \
	err = MPI_Get_count (V0, V1, V2);                     \
	if (err != MPI_SUCCESS) {                             \
		err = ncmpii_error_mpi2nc (err, "MPI_Get_count"); \
		DEBUG_RETURN_ERROR (err)                          \
	}
#define CHK_ERR_IMRECV(V0, V1, V2, V3, V4)             \
	err = MPI_Imrecv (V0, V1, V2, V3, V4);             \
	if (err != MPI_SUCCESS) {                          \
		err = ncmpii_error_mpi2nc (err, "MPI_Imrecv"); \
		DEBUG_RETURN_ERROR (err)                       \
	}
#define CHK_ERR_ISEND(V0, V1, V2, V3, V4, V5, V6)     \
	err = MPI_Isend (V0, V1, V2, V3, V4, V5, V6);     \
	if (err != MPI_SUCCESS) {                         \
		err = ncmpii_error_mpi2nc (err, "MPI_Isend"); \
		DEBUG_RETURN_ERROR (err)                      \
	}
#define CHK_ERR_IRECV(V0, V1, V2, V3, V4, V5, V6)     \
	err = MPI_Irecv (V0, V1, V2, V3, V4, V5, V6);     \
	if (err != MPI_SUCCESS) {                         \
		err = ncmpii_error_mpi2nc (err, "MPI_Irecv"); \
		DEBUG_RETURN_ERROR (err)                      \
	}
#define CHK_ERR_SET_VIEW(V0, V1, V2, V3, V4, V5)              \
	err = MPI_File_set_view (V0, V1, V2, V3, V4, V5);         \
	if (err != MPI_SUCCESS) {                                 \
		err = ncmpii_error_mpi2nc (err, "MPI_File_set_view"); \
		DEBUG_RETURN_ERROR (err)                              \
	}
#define CHK_ERR_READ_AT_ALL(V0, V1, V2, V3, V4, V5)              \
	err = MPI_File_read_at_all (V0, V1, V2, V3, V4, V5);         \
	if (err != MPI_SUCCESS) {                                    \
		err = ncmpii_error_mpi2nc (err, "MPI_File_read_at_all"); \
		DEBUG_RETURN_ERROR (err)                                 \
	}
#define CHK_ERR_WRITE_AT_ALL(V0, V1, V2, V3, V4, V5)              \
	err = MPI_File_write_at_all (V0, V1, V2, V3, V4, V5);         \
	if (err != MPI_SUCCESS) {                                     \
		err = ncmpii_error_mpi2nc (err, "MPI_File_write_at_all"); \
		DEBUG_RETURN_ERROR (err)                                  \
	}
#define CHK_ALLOC(V0) \
	if (V0 == NULL) { DEBUG_RETURN_ERROR (NC_ENOMEM) }

typedef struct NC_zip_vector {
	int esize;
	int size;
	int nalloc;
	char *data;
} NC_zip_vector;

// File
extern int nczipioi_init (NC_zip *, int);
extern int nczipioi_parse_var_info (NC_zip *);
extern int nczipioi_var_list_init (NC_zip_var_list *);
extern int nczipioi_var_list_free (NC_zip_var_list *);
extern int nczipioi_var_list_add (NC_zip_var_list *);

// Util
extern int nczipioi_extract_hint (NC_zip *, MPI_Info);
extern int nczipioi_export_hint (NC_zip *, MPI_Info);
extern MPI_Offset NC_Type_size (nc_type);
extern int nczipioi_print_profile (NC_zip *);
extern void nczipioi_sort_file_offset (int, MPI_Aint *, MPI_Aint *, int *);
extern int nczipioi_update_statistics (NC_zip *);
extern int nczipioi_get_default_chunk_dim (NC_zip *);
extern int nczipioi_subarray_off_len (int, int *, int *, int *, int *, int *);
extern void nczipioi_idx_in_swapn (NC_zip_chunk_index_entry *, MPI_Offset);

// Misc
extern int nczipioi_calc_chunk_owner (NC_zip *, NC_zip_var *, int, MPI_Offset **, MPI_Offset **);
extern int nczipioi_calc_chunk_size (NC_zip *, NC_zip_var *, int, MPI_Offset **, MPI_Offset **);
extern int nczipioiconvert (void *, void *, MPI_Datatype, MPI_Datatype, int);

// Var
extern int nczipioi_var_init (NC_zip *, NC_zip_var *, int, MPI_Offset **, MPI_Offset **);
extern int nczipioi_load_var (NC_zip *, NC_zip_var *, int, int *);
extern int nczipioi_load_var_bg (NC_zip *, NC_zip_var *, int, int *);
extern int nczipioi_load_nvar (NC_zip *, int, int *, int *, int *);
extern int nczipioi_load_nvar_bg (NC_zip *, int, int *, int *, int *);
extern int nczipioi_save_var (NC_zip *, NC_zip_var *);
extern int nczipioi_save_nvar (NC_zip *, int, int *);
extern void nczipioi_var_free (NC_zip_var *);
extern int nczipioi_var_resize (NC_zip *, NC_zip_var *);
extern int nczipioi_init_nvar (NC_zip *, int, int *, int, int *);
extern int nczipioi_resize_nvar (NC_zip *, int, int *, int, int *);

// Cache
extern int nczipioi_cache_alloc (NC_zip *, MPI_Offset, NC_zip_cache **);
extern void nczipioi_cache_visit (NC_zip *, NC_zip_cache *);
extern void nczipioi_cache_free (NC_zip *);

// Chunks
extern int nczipioi_chunk_itr_init (
	NC_zip_var *, const MPI_Offset *, const MPI_Offset *, MPI_Offset *, int *);
extern int nczipioi_chunk_itr_next (
	NC_zip_var *, const MPI_Offset *, const MPI_Offset *, MPI_Offset *, int *);
extern MPI_Offset get_chunk_overlap (
	NC_zip_var *, MPI_Offset *, const MPI_Offset *, const MPI_Offset *, MPI_Offset *, MPI_Offset *);
extern int get_chunk_id (NC_zip_var *, MPI_Offset *);
extern int get_chunk_itr (NC_zip_var *, int, MPI_Offset *);
extern int nczipioi_chunk_itr_init_ex (NC_zip_var *,
									   const MPI_Offset *,
									   const MPI_Offset *,
									   MPI_Offset *,
									   int *,
									   MPI_Offset *,
									   MPI_Offset *);
extern int nczipioi_chunk_itr_next_ex (NC_zip_var *,
									   const MPI_Offset *,
									   const MPI_Offset *,
									   MPI_Offset *,
									   int *,
									   MPI_Offset *,
									   MPI_Offset *);

// Get
// extern int nczipioi_get_var_old(NC_zip*, NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*,
// void*);
extern int nczipioi_get_var_cb_chunk (
	NC_zip *, NC_zip_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int nczipioi_get_var_cb_proc (
	NC_zip *, NC_zip_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int nczipioi_get_varn (
	NC_zip *, NC_zip_var *, int, MPI_Offset *const *, MPI_Offset *const *, const void *);
extern int nczipioi_get_varn_cb_chunk (NC_zip *,
									   NC_zip_var *,
									   int,
									   MPI_Offset *const *,
									   MPI_Offset *const *,
									   MPI_Offset *const *,
									   void **);
extern int nczipioi_get_varn_cb_proc (
	NC_zip *, NC_zip_var *, int, MPI_Offset *const *, MPI_Offset *const *, void **);
extern int nczipioi_iget_var (NC_zip *,
							  int,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  void *,
							  MPI_Offset,
							  MPI_Datatype,
							  int *);
extern int nczipioi_iget_varn (NC_zip *,
							   int,
							   int,
							   MPI_Offset *const *,
							   MPI_Offset *const *,
							   void *,
							   MPI_Offset,
							   MPI_Datatype,
							   int *);
extern int nczipioi_iget_cb_chunk (NC_zip *, int, int *, int *);
extern int nczipioi_iget_cb_proc (NC_zip *, int, int *, int *);

// Put
// extern int nczipioi_put_var_old(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const
// MPI_Offset*, void*);
extern int nczipioi_put_var (
	NC_zip *, NC_zip_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int nczipioi_put_var_cb_chunk (
	NC_zip *, NC_zip_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int nczipioi_put_var_cb_proc (
	NC_zip *, NC_zip_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int nczipioi_put_varn (
	NC_zip *, NC_zip_var *, int, MPI_Offset *const *, MPI_Offset *const *, const void *);
extern int nczipioi_put_varn_cb_chunk (NC_zip *,
									   NC_zip_var *,
									   int,
									   MPI_Offset *const *,
									   MPI_Offset *const *,
									   MPI_Offset *const *,
									   void **);
extern int nczipioi_put_varn_cb_proc (
	NC_zip *, NC_zip_var *, int, MPI_Offset *const *, MPI_Offset *const *, void **);
extern int nczipioi_iput_var (NC_zip *,
							  int,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  const void *,
							  const void *,
							  int *);
extern int nczipioi_iput_varn (NC_zip *,
							   int,
							   int,
							   MPI_Offset *const *,
							   MPI_Offset *const *,
							   const void *,
							   const void *,
							   int *);
extern int nczipioi_iput_cb_chunk (NC_zip *, int, int *, int *);
extern int nczipioi_iput_cb_proc (NC_zip *, int, int *, int *);

// Nonblocking
extern int nczipioi_req_list_init (NC_zip_req_list *);
extern int nczipioi_req_list_free (NC_zip_req_list *);
extern int nczipioi_req_list_add (NC_zip_req_list *, int *);
extern int nczipioi_req_list_remove (NC_zip_req_list *, int);
extern int nczipioi_wait_put_reqs (NC_zip *, int, int *, int *);
extern int nczipioi_wait_get_reqs (NC_zip *, int, int *, int *);
extern int nczipioi_wait (NC_zip *, int, int *, int *, int);

// Vector
extern int nczipioi_vector_init (NC_zip_vector *, int);
extern int nczipioi_vector_init_ex (NC_zip_vector *, int, int);
extern void nczipioi_vector_free (NC_zip_vector *);
extern int nczipioi_vector_append (NC_zip_vector *, void *);
#endif