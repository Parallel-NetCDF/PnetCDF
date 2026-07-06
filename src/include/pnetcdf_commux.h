/*
 *  Copyright (C) 2026, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * A small C ABI for the installed Python commux package.
 */

#ifndef H_PNETCDF_COMMUX
#define H_PNETCDF_COMMUX

#include <stdint.h>

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #define PNC_COMMUX_API __declspec(dllexport)
  #else
    #define PNC_COMMUX_API __declspec(dllimport)
  #endif
#else
  #if __GNUC__ >= 4
    #define PNC_COMMUX_API __attribute__ ((visibility ("default")))
  #else
    #define PNC_COMMUX_API
  #endif
#endif

#define PNC_COMMUX_SUCCESS      0
#define PNC_COMMUX_ERR_RUNTIME -1
#define PNC_COMMUX_ERR_ARG     -2

typedef enum PNC_CommuxReduceOp {
    PNC_COMMUX_SUM = 0,
    PNC_COMMUX_MIN = 1,
    PNC_COMMUX_MAX = 2
} PNC_CommuxReduceOp;

#if defined(__cplusplus)
extern "C" {
#endif

PNC_COMMUX_API const char*
pnc_commux_last_error(void);

PNC_COMMUX_API int
pnc_commux_init(const char *backend, const char *init_method,
                int rank, int world_size);

PNC_COMMUX_API int
pnc_commux_init_env(const char *backend);

PNC_COMMUX_API int
pnc_commux_finalize(void);

PNC_COMMUX_API int
pnc_commux_rank(int *rankp);

PNC_COMMUX_API int
pnc_commux_size(int *sizep);

PNC_COMMUX_API int
pnc_commux_barrier(void);

PNC_COMMUX_API int
pnc_commux_allreduce_int64(int64_t value, PNC_CommuxReduceOp op,
                           int64_t *resultp);

PNC_COMMUX_API int
pnc_commux_send_int64(int64_t value, int dst, int tag);

PNC_COMMUX_API int
pnc_commux_recv_int64(int src, int tag, int64_t *valuep);

PNC_COMMUX_API int
pnc_commux_send_bytes(const void *data, int64_t nbytes, int dst, int tag);

PNC_COMMUX_API int
pnc_commux_recv_bytes(void *data, int64_t nbytes, int src, int tag);

PNC_COMMUX_API int
pnc_commux_file_open(const char *path, int flags, int mode, int64_t *handlep);

PNC_COMMUX_API int
pnc_commux_file_close(int64_t handle);

PNC_COMMUX_API int
pnc_commux_file_pwrite(int64_t handle, const void *data, int64_t nbytes,
                       int64_t offset);

PNC_COMMUX_API int
pnc_commux_file_sync(int64_t handle);

#if defined(__cplusplus)
}
#endif

#endif /* H_PNETCDF_COMMUX */
