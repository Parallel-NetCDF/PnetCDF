#ifndef _nczipio_INTERNAL_H
#define _nczipio_INTERNAL_H

#include "nczipio_driver.h"

#define CHK_ERR_ALLREDUCE(V0,V1,V2,V3,V4,V5) \
        err = MPI_Allreduce(V0,V1,V2,V3,V4,V5); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Allreduce"); \
            DEBUG_RETURN_ERROR(err) \
        }

#define CHK_ERR_PACK(V0,V1,V2,V3,V4,V5,V6) \
        err = MPI_Pack(V0,V1,V2,V3,V4,V5,V6); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Pack"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_UNPACK(V0,V1,V2,V3,V4,V5,V6) \
        err = MPI_Unpack(V0,V1,V2,V3,V4,V5,V6); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Unpack"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_TYPE_COMMIT(V0) \
        err = MPI_Type_commit(V0); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Type_commit"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_TYPE_CREATE_SUBARRAY(V0,V1,V2,V3,V4,V5,V6) \
        err = MPI_Type_create_subarray(V0,V1,V2,V3,V4,V5,V6); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Type_create_subarray"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_WAITALL(V0,V1,V2) \
        err = MPI_Waitall(V0,V1,V2); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Waitall"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_MPROBE(V0,V1,V2,V3,V4) \
        err = MPI_Mprobe(V0,V1,V2,V3,V4); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Mprobe"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_GET_COUNT(V0,V1,V2) \
        err = MPI_Get_count(V0,V1,V2); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Get_count"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_IMRECV(V0,V1,V2,V3,V4) \
        err = MPI_Imrecv(V0,V1,V2,V3,V4); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Imrecv"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_ISEND(V0,V1,V2,V3,V4,V5,V6) \
        err = MPI_Isend(V0,V1,V2,V3,V4,V5,V6); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Isend"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_IRECV(V0,V1,V2,V3,V4,V5,V6) \
        err = MPI_Irecv(V0,V1,V2,V3,V4,V5,V6); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_Irecv"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_SET_VIEW(V0,V1,V2,V3,V4,V5) \
        err = MPI_File_set_view(V0,V1,V2,V3,V4,V5); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_File_set_view"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_READ_AT_ALL(V0,V1,V2,V3,V4,V5) \
        err = MPI_File_read_at_all(V0,V1,V2,V3,V4,V5); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_File_read_at_all"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ERR_WRITE_AT_ALL(V0,V1,V2,V3,V4,V5) \
        err = MPI_File_write_at_all(V0,V1,V2,V3,V4,V5); \
        if (err != MPI_SUCCESS){ \
            err = ncmpii_error_mpi2nc(err, "MPI_File_write_at_all"); \
            DEBUG_RETURN_ERROR(err) \
        }
#define CHK_ALLOC(V0) \
        if (V0 == NULL){ \
            DEBUG_RETURN_ERROR(NC_ENOMEM) \
        }

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
extern int nczipioi_print_profile(NC_zip*);

// Misc
extern int nczipioi_calc_chunk_owner(NC_zip*, NC_zip_var*, int, MPI_Offset**, MPI_Offset**);
extern int nczipioi_calc_chunk_size(NC_zip*, NC_zip_var*, int, MPI_Offset**, MPI_Offset**);

// Var
extern int nczipioi_var_init(NC_zip*, NC_zip_var*, int, MPI_Offset**, MPI_Offset**);
extern int nczipioi_load_var(NC_zip*, NC_zip_var*, int, int*);
extern int nczipioi_load_nvar(NC_zip*, int, int*);
extern int nczipioi_save_var(NC_zip*, NC_zip_var*);
extern int nczipioi_save_nvar(NC_zip*, int, int*);
extern void nczipioi_var_free(NC_zip_var*);
extern int nczipioi_var_resize(NC_zip*, NC_zip_var*);
extern int nczipioi_init_nvar(NC_zip*, int, int*, int, int*);

// Chunks
extern int nczipioi_chunk_itr_init(NC_zip_var*, const MPI_Offset*, const MPI_Offset*, MPI_Offset*, int*);
extern int nczipioi_chunk_itr_next(NC_zip_var*, const MPI_Offset*, const MPI_Offset*, MPI_Offset*, int*);
extern int get_chunk_overlap(NC_zip_var*, MPI_Offset*, const MPI_Offset*, const MPI_Offset*, MPI_Offset*, MPI_Offset*);
extern int get_chunk_id(NC_zip_var*, MPI_Offset*);
extern int get_chunk_itr(NC_zip_var*, int, MPI_Offset*);

// Get
//extern int nczipioi_get_var_old(NC_zip*, NC_zip_var*, MPI_Offset*, MPI_Offset*, MPI_Offset*, void*);
extern int nczipioi_get_var_cb_chunk(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);
extern int nczipioi_get_var_cb_proc(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);
extern int nczipioi_get_varn(NC_zip *, NC_zip_var*, int , MPI_Offset* const *, MPI_Offset* const *, const void*);
extern int nczipioi_get_varn_cb_chunk(NC_zip*, NC_zip_var*, int nreq, MPI_Offset* const*, MPI_Offset* const*, MPI_Offset* const*, void**);
extern int nczipioi_get_varn_cb_proc(NC_zip*, NC_zip_var*, int nreq, MPI_Offset* const*, MPI_Offset* const*, void**);
extern int nczipioi_iget_var(NC_zip*, int, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*, MPI_Offset, MPI_Datatype, int*);
extern int nczipioi_iget_varn(NC_zip*, int, int, MPI_Offset * const*, MPI_Offset * const*, void*, MPI_Offset, MPI_Datatype, int*);
extern int nczipioi_iget_cb_chunk(NC_zip*, int, int*, int*);
extern int nczipioi_iget_cb_proc(NC_zip*, int, int*, int*);

// Put
//extern int nczipioi_put_var_old(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);
extern int nczipioi_put_var(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);
extern int nczipioi_put_var_cb_chunk(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);
extern int nczipioi_put_var_cb_proc(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);
extern int nczipioi_put_varn(NC_zip*, NC_zip_var*, int, MPI_Offset* const*, MPI_Offset* const *, const void*);
extern int nczipioi_put_varn_cb_chunk(NC_zip*, NC_zip_var*, int, MPI_Offset* const*, MPI_Offset* const*, MPI_Offset* const*, void**);
extern int nczipioi_put_varn_cb_proc(NC_zip*, NC_zip_var*, int, MPI_Offset* const*, MPI_Offset* const*, void**);
extern int nczipioi_iput_var(NC_zip*, int, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, const void*, const void*, int*);
extern int nczipioi_iput_varn(NC_zip*, int, int, MPI_Offset * const*, MPI_Offset * const*, const void*, const void*, int*);
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