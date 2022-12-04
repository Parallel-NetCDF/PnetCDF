#ifndef _ncchkio_INTERNAL_H
#define _ncchkio_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "ncchkio_driver.h"
#ifdef PNETCDF_DEBUG
#include <assert.h>
#endif

#define NC_CHK_FILTER_NONE	0
#define NC_CHK_FILTER_DUMMY 1
#define NC_CHK_FILTER_ZLIB	2
#define NC_CHK_FILTER_SZ	3

#define NC_CHK_DEFAULT_REC_ALLOC 1024
#define NC_CHK_REC_MULTIPLIER	 2

#ifdef PNETCDF_DEBUG
#define DEBUG_ABORT                                             \
	{                                                           \
		char *_env_str = getenv ("PNETCDF_ABORT_ON_ERR");       \
		if (_env_str != NULL && *_env_str != '0') { abort (); } \
	}
#else
#define DEBUG_ABORT
#endif

#define RET_ERR(E)               \
	{                            \
		err = E;                 \
		DEBUG_TRACE_ERROR (err); \
		DEBUG_ABORT              \
		goto err_out;            \
	}
#define CHK_ERR            \
	if (err != NC_NOERR) { \
		DEBUG_ABORT        \
		goto err_out;      \
	}

#define CHK_MPIERR                              \
	if (err != MPI_SUCCESS) {                   \
		err = ncmpii_error_mpi2nc (err, "MPI"); \
		DEBUG_TRACE_ERROR (err);                \
		DEBUG_ABORT                             \
		goto err_out;                           \
	}

#define CHK_PTR(P)               \
	if (!P) {                    \
		err = NC_ENOMEM;         \
		DEBUG_TRACE_ERROR (err); \
		DEBUG_ABORT              \
		goto err_out;            \
	}

#define CHK_ERR_WAIT(V0, V1) \
	err = MPI_Wait (V0, V1); \
	CHK_MPIERR

#define CHK_ERR_ALLREDUCE(V0, V1, V2, V3, V4, V5) \
	err = MPI_Allreduce (V0, V1, V2, V3, V4, V5); \
	CHK_MPIERR

#define CHK_ERR_IALLREDUCE(V0, V1, V2, V3, V4, V5, V6) \
	err = MPI_Iallreduce (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR

#define CHK_ERR_REDUCE(V0, V1, V2, V3, V4, V5, V6) \
	err = MPI_Reduce (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR

#define CHK_ERR_GATHER(V0, V1, V2, V3, V4, V5, V6, V7) \
	err = MPI_Gather (V0, V1, V2, V3, V4, V5, V6, V7); \
	CHK_MPIERR

#ifdef PNETCDF_DEBUG
#define CHK_ERR_PACK(V0, V1, V2, V3, V4, V5, V6)     \
	{                                                \
		assert ((V0) != NULL);                       \
		assert ((V3) != NULL);                       \
		err = MPI_Pack (V0, V1, V2, V3, V4, V5, V6); \
		CHK_MPIERR                                   \
	}
#else
#define CHK_ERR_PACK(V0, V1, V2, V3, V4, V5, V6) \
	err = MPI_Pack (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR
#endif

#ifdef PNETCDF_DEBUG
#define CHK_ERR_UNPACK(V0, V1, V2, V3, V4, V5, V6)          \
	{                                                       \
		int esize;                                          \
		MPI_Type_size (V5, &esize);                         \
		if (V1 - *((int *)(V2)) < V4 * esize) { abort (); } \
		err = MPI_Unpack (V0, V1, V2, V3, V4, V5, V6);      \
		CHK_MPIERR                                          \
	}
#else
#define CHK_ERR_UNPACK(V0, V1, V2, V3, V4, V5, V6) \
	err = MPI_Unpack (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR
#endif

#define CHK_ERR_TYPE_COMMIT(V0) \
	err = MPI_Type_commit (V0); \
	CHK_MPIERR

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
		CHK_MPIERR                                                                             \
	}
#else
#define CHK_ERR_TYPE_CREATE_SUBARRAY(V0, V1, V2, V3, V4, V5, V6) \
	err = MPI_Type_create_subarray (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR
#endif

#define CHK_ERR_WAITALL(V0, V1, V2) \
	err = MPI_Waitall (V0, V1, V2); \
	CHK_MPIERR
#define CHK_ERR_MPROBE(V0, V1, V2, V3, V4) \
	err = MPI_Mprobe (V0, V1, V2, V3, V4); \
	CHK_MPIERR

#define CHK_ERR_GET_COUNT(V0, V1, V2) \
	err = MPI_Get_count (V0, V1, V2); \
	CHK_MPIERR

#define CHK_ERR_IMRECV(V0, V1, V2, V3, V4) \
	err = MPI_Imrecv (V0, V1, V2, V3, V4); \
	CHK_MPIERR

#ifdef PNETCDF_DEBUG
#define CHK_ERR_ISEND(V0, V1, V2, V3, V4, V5, V6) \
	assert (V1 >= 0);                             \
	err = MPI_Isend (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR
#else
#define CHK_ERR_ISEND(V0, V1, V2, V3, V4, V5, V6) \
	err = MPI_Isend (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR
#endif

#ifdef PNETCDF_DEBUG
#define CHK_ERR_IRECV(V0, V1, V2, V3, V4, V5, V6) \
	assert (V1 >= 0);                             \
	err = MPI_Irecv (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR
#else
#define CHK_ERR_IRECV(V0, V1, V2, V3, V4, V5, V6) \
	err = MPI_Irecv (V0, V1, V2, V3, V4, V5, V6); \
	CHK_MPIERR
#endif

#define CHK_ERR_SET_VIEW(V0, V1, V2, V3, V4, V5)      \
	err = MPI_File_set_view (V0, V1, V2, V3, V4, V5); \
	CHK_MPIERR

#define CHK_ERR_READ_AT_ALL(V0, V1, V2, V3, V4, V5)      \
	err = MPI_File_read_at_all (V0, V1, V2, V3, V4, V5); \
	CHK_MPIERR

#define CHK_ERR_WRITE_AT_ALL(V0, V1, V2, V3, V4, V5)      \
	err = MPI_File_write_at_all (V0, V1, V2, V3, V4, V5); \
	CHK_MPIERR

#define CHK_ALLOC(V0) \
	if (V0 == NULL) { DEBUG_RETURN_ERROR (NC_ENOMEM) }

typedef struct NC_chk_vector {
	int esize;
	int size;
	int nalloc;
	char *data;
} NC_chk_vector;

// File
extern int ncchkioi_init (NC_chk *, int);
extern int ncchkioi_parse_var_info (NC_chk *);
extern int ncchkioi_var_list_init (NC_chk_var_list *);
extern int ncchkioi_var_list_free (NC_chk_var_list *);
extern int ncchkioi_var_list_add (NC_chk_var_list *);

// Util
extern int ncchkioi_extract_hint (NC_chk *, MPI_Info);
extern int ncchkioi_export_hint (NC_chk *, MPI_Info);
extern MPI_Offset NC_Type_size (nc_type);
extern void ncchkioi_sort_file_offset (int, MPI_Aint *, MPI_Aint *, int *);
extern int ncchkioi_update_statistics (NC_chk *);
extern int ncchkioi_get_default_chunk_dim (NC_chk *);
extern int ncchkioi_subarray_off_len (int, int *, int *, int *, MPI_Offset *, int *);
extern void ncchkioi_idx_in_swapn (NC_chk_chunk_index_entry *, MPI_Offset);
#ifdef PNETCDF_PROFILING
extern int ncchkioi_print_profile (NC_chk *);
extern void ncchkioi_profile_add_time (NC_chk *ncchkp, int id, double t);
#endif

// Misc
typedef struct ncchkioi_chunk_overlap_t {
	MPI_Offset osize;
	int rank;
} ncchkioi_chunk_overlap_t;
extern int ncchkioi_init_nvar_core_reduce (NC_chk *ncchkp,
										   int nvar,
										   NC_chk_var **varps,
										   int *rcnt,
										   int *roff,
										   MPI_Offset **starts,
										   MPI_Offset **counts);
extern int ncchkioi_calc_chunk_overlap (NC_chk *ncchkp,
										NC_chk_var *varp,
										int nreq,
										MPI_Offset **starts,
										MPI_Offset **counts,
										ncchkioi_chunk_overlap_t *ocnt);
extern void ncchkioi_assign_chunk_owner (NC_chk *ncchkp,
										 NC_chk_var *varp,
										 ncchkioi_chunk_overlap_t *ocnt);
extern int ncchkioi_sync_ocnt_reduce (NC_chk *ncchkp,
									  int nchunk,
									  ncchkioi_chunk_overlap_t *ocnt,
									  ncchkioi_chunk_overlap_t *ocnt_all,
									  MPI_Request *req);
extern void ncchkioi_write_chunk_ocnt (NC_chk *ncchkp,
									   NC_chk_var *varp,
									   void *ocnt,
									   size_t ocnt_size);
extern int ncchkioi_calc_chunk_owner (NC_chk *, NC_chk_var *, int, MPI_Offset **, MPI_Offset **);
extern int ncchkioi_calc_chunk_owner_reduce (
	NC_chk *ncchkp, NC_chk_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts);
extern int ncchkioi_calc_chunk_size (NC_chk *, NC_chk_var *, int, MPI_Offset **, MPI_Offset **);
extern int ncchkioiconvert (void *, void *, MPI_Datatype, MPI_Datatype, int);

// Var
extern int ncchkioi_var_init (NC_chk *, NC_chk_var *, int, MPI_Offset **, MPI_Offset **);
extern int ncchkioi_load_var (NC_chk *, NC_chk_var *, int, int *);
extern int ncchkioi_load_var_bg (NC_chk *, NC_chk_var *, int, int *);
extern int ncchkioi_load_nvar (NC_chk *, int, int *, int *, int *);
extern int ncchkioi_load_nvar_bg (NC_chk *, int, int *, int *, int *);
extern int ncchkioi_save_var (NC_chk *, NC_chk_var *);
extern int ncchkioi_save_nvar (NC_chk *, int, int *);
extern void ncchkioi_var_free (NC_chk_var *);
extern int ncchkioi_var_resize (NC_chk *, NC_chk_var *);
extern int ncchkioi_init_nvar (NC_chk *, int, int *, int, int *);
extern int ncchkioi_resize_nvar (NC_chk *, int, int *, int, int *);

// Cache
extern int ncchkioi_cache_alloc (NC_chk *, MPI_Offset, NC_chk_cache **);
extern void ncchkioi_cache_visit (NC_chk *, NC_chk_cache *);
extern void ncchkioi_cache_free (NC_chk *);

// Chunks
extern int ncchkioi_chunk_itr_init (
	NC_chk_var *, const MPI_Offset *, const MPI_Offset *, MPI_Offset *, int *);
extern int ncchkioi_chunk_itr_next (
	NC_chk_var *, const MPI_Offset *, const MPI_Offset *, MPI_Offset *, int *);
extern MPI_Offset get_chunk_overlap (
	NC_chk_var *, MPI_Offset *, const MPI_Offset *, const MPI_Offset *, MPI_Offset *, MPI_Offset *);
extern int get_chunk_id (NC_chk_var *, MPI_Offset *);
extern int get_chunk_itr (NC_chk_var *, int, MPI_Offset *);
extern int ncchkioi_chunk_itr_init_ex (NC_chk_var *,
									   const MPI_Offset *,
									   const MPI_Offset *,
									   MPI_Offset *,
									   int *,
									   MPI_Offset *,
									   MPI_Offset *);
extern int ncchkioi_chunk_itr_next_ex (NC_chk_var *,
									   const MPI_Offset *,
									   const MPI_Offset *,
									   MPI_Offset *,
									   int *,
									   MPI_Offset *,
									   MPI_Offset *);

// Get
// extern int ncchkioi_get_var_old(NC_chk*, NC_chk_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*,
// void*);
extern int ncchkioi_get_var_cb_chunk (
	NC_chk *, NC_chk_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int ncchkioi_get_var_cb_proc (
	NC_chk *, NC_chk_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int ncchkioi_get_varn (
	NC_chk *, NC_chk_var *, int, MPI_Offset *const *, MPI_Offset *const *, const void *);
extern int ncchkioi_get_varn_cb_chunk (NC_chk *,
									   NC_chk_var *,
									   int,
									   MPI_Offset *const *,
									   MPI_Offset *const *,
									   MPI_Offset *const *,
									   void **);
extern int ncchkioi_get_varn_cb_proc (
	NC_chk *, NC_chk_var *, int, MPI_Offset *const *, MPI_Offset *const *, void **);
extern int ncchkioi_iget_var (NC_chk *,
							  int,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  void *,
							  MPI_Offset,
							  MPI_Datatype,
							  int *);
extern int ncchkioi_iget_varn (NC_chk *,
							   int,
							   int,
							   MPI_Offset *const *,
							   MPI_Offset *const *,
							   void *,
							   MPI_Offset,
							   MPI_Datatype,
							   int *);
extern int ncchkioi_iget_cb_chunk (NC_chk *, int, int *, int *);
extern int ncchkioi_iget_cb_proc (NC_chk *, int, int *, int *);

// Put
// extern int ncchkioi_put_var_old(NC_chk*, NC_chk_var*, const MPI_Offset*, const MPI_Offset*, const
// MPI_Offset*, void*);
extern int ncchkioi_put_var (
	NC_chk *, NC_chk_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int ncchkioi_put_var_cb_chunk (
	NC_chk *, NC_chk_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int ncchkioi_put_var_cb_proc (
	NC_chk *, NC_chk_var *, const MPI_Offset *, const MPI_Offset *, const MPI_Offset *, void *);
extern int ncchkioi_put_varn (
	NC_chk *, NC_chk_var *, int, MPI_Offset *const *, MPI_Offset *const *, const void *);
extern int ncchkioi_put_varn_cb_chunk (NC_chk *,
									   NC_chk_var *,
									   int,
									   MPI_Offset *const *,
									   MPI_Offset *const *,
									   MPI_Offset *const *,
									   void **);
extern int ncchkioi_put_varn_cb_proc (
	NC_chk *, NC_chk_var *, int, MPI_Offset *const *, MPI_Offset *const *, void **);
extern int ncchkioi_iput_var (NC_chk *,
							  int,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  const MPI_Offset *,
							  const void *,
							  const void *,
							  int *);
extern int ncchkioi_iput_varn (NC_chk *,
							   int,
							   int,
							   MPI_Offset *const *,
							   MPI_Offset *const *,
							   const void *,
							   const void *,
							   int *);
extern int ncchkioi_iput_cb_chunk (NC_chk *, int, int *, int *);
extern int ncchkioi_iput_cb_proc (NC_chk *, int, int *, int *);

// Nonblocking
extern int ncchkioi_req_list_init (NC_chk_req_list *);
extern int ncchkioi_req_list_free (NC_chk_req_list *);
extern int ncchkioi_req_list_add (NC_chk_req_list *, int *);
extern int ncchkioi_req_list_remove (NC_chk_req_list *, int);
extern int ncchkioi_wait_put_reqs (NC_chk *, int, int *, int *);
extern int ncchkioi_wait_get_reqs (NC_chk *, int, int *, int *);
extern int ncchkioi_wait (NC_chk *, int, int *, int *, int);

// Vector
extern int ncchkioi_vector_init (NC_chk_vector *, int);
extern int ncchkioi_vector_init_ex (NC_chk_vector *, int, int);
extern void ncchkioi_vector_free (NC_chk_vector *);
extern int ncchkioi_vector_append (NC_chk_vector *, void *);
#endif