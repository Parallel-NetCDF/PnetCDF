/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifndef H_PNCIO
#define H_PNCIO

#include <stdio.h>
#include <stdlib.h>
#include <sys/errno.h>
#include <unistd.h>   /* pwrite() */

#include <stdbool.h>
#include <string.h>     /* memcpy() */
#include <stddef.h>     /* size_t */
#include <sys/types.h>  /* off_t */
#include <assert.h>
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#define FDTYPE int

#include <pnc_debug.h>
#include <dispatch.h>
#include <common.h>

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
#define NMEASURES 8
#endif

#define PNCIO_LOCKS  300    /* file system supports fcntl()-style locking */
#define PNCIO_Feature(a, b) ((b == PNCIO_LOCKS) ? 1 : 0)

#if defined(F_SETLKW64)
#define PNCIO_UNLOCK(fd, offset, whence, len) \
        PNCIO_GEN_SetLock64(fd, F_SETLK, F_UNLCK, offset, whence, len)
#define PNCIO_WRITE_LOCK(fd, offset, whence, len) \
        PNCIO_GEN_SetLock64(fd, F_SETLKW, F_WRLCK, offset, whence, len)
#else
#define PNCIO_UNLOCK(fd, offset, whence, len) \
        PNCIO_GEN_SetLock(fd, F_SETLK, F_UNLCK, offset, whence, len)
#define PNCIO_WRITE_LOCK(fd, offset, whence, len) \
        PNCIO_GEN_SetLock(fd, F_SETLKW, F_WRLCK, offset, whence, len)
#endif


#define PNCIO_PERM          0666   /* file creation permission mask */

#define PNCIO_UFS           152    /* Unix file system */
#define PNCIO_LUSTRE        163    /* Lustre */
#define PNCIO_FSTYPE_MPIIO  -1     /* Use MPI-IO */
#define PNCIO_FSTYPE_CHECK  0      /* Use PnetCDF PNCIO drivers */

#define PNCIO_LUSTRE_MAX_OSTS 256  /* Maximum number of Lustre OSTs if hint
                                    * striping_factor is not set by user.
                                    */

#define PNCIO_CB_BUFFER_SIZE_DFLT     "16777216"
#define PNCIO_IND_RD_BUFFER_SIZE_DFLT "4194304"
#define PNCIO_IND_WR_BUFFER_SIZE_DFLT "524288"
#define PNCIO_CB_CONFIG_LIST_DFLT     "*:1"

/* PNCIO_DS_WR_NPAIRS_LB is the lower bound of the total number of
 *     offset-length pairs over the non-aggregator senders to be received by an
 *     I/O aggregator to skip the potentially expensive heap-merge sort that
 *     determines whether or not data sieving write is necessary.
 * PNCIO_DS_WR_NAGGRS_LB is the lower bound of the number of non-aggregators
 *     sending their offset-length pairs to an I/O aggregator.
 * Both conditions must be met to skip the heap-merge sort.
 *
 * When data sieving is enabled, read-modify-write will perform at each round
 * of two-phase I/O at each aggregator. The following describes whether
 * detecting "holes" in a write region is necessary, depending on the data
 * sieving hint, romio_ds_write, is set to enable/disable/automatic.
 *   + automatic - We need to check whether holes exist. If holes exist, the
 *       "read-modify" part must run. If not, "read-modify" can be skipped.
 *   + enable - "read-modify" part must perform, skip hole checking, and thus
 *       skip the heap-merge sort.
 *   + disable - "read-modify" part must skip, need not check holes, but must
 *       construct srt_off_len to merge all others_req[] into a single sorted
 *       list, which requires to call a heap-merge sort. This step is necessary
 *       because write data from all non-aggregators are received into the same
 *       write_buf, with a possibility of overlaps, and srt_off_len stores the
 *       coalesced offset-length pairs of individual non-contiguous write
 *       request and will be used to write them to the file.
 *
 * Heap-merge sort merges offset-length pairs received from all non-aggregators
 * into a single list, which can be expensive. Its cost can be even larger than
 * the cost of "read" in "read-modify-write". Below two constants are the lower
 * bounds used to determine whether or not to perform such sorting, when data
 * sieving is set to the automatic mode.
 */
#define PNCIO_DS_WR_NPAIRS_LB 8192
#define PNCIO_DS_WR_NAGGRS_LB 256
#define DO_HEAP_MERGE(nrecv, npairs) ((nrecv) > PNCIO_DS_WR_NAGGRS_LB || (npairs) > PNCIO_DS_WR_NPAIRS_LB)

#define PNCIO_TYPE_DECREASE 0x00000001  /* if not monotonic nondecreasing */
#define PNCIO_TYPE_OVERLAP  0x00000002  /* if contains overlapping regions */
#define PNCIO_TYPE_NEGATIVE 0x00000004  /* if one of displacements is negative */

#define PNCIO_HINT_AUTO -1
#define PNCIO_HINT_DISABLE 0
#define PNCIO_HINT_ENABLE 1

#define PNCIO_STRIPING_AUTO -1
#define PNCIO_STRIPING_INHERIT 0

typedef struct {
    int nc_striping;
    int striping_factor;
    int striping_unit;
    int start_iodevice;
    int cb_nodes;
    int cb_buffer_size;
    int ind_rd_buffer_size;
    int ind_wr_buffer_size;

    int romio_cb_read;
    int romio_cb_write;
    int romio_ds_read;
    int romio_ds_write;
    int romio_no_indep_rw;

    /* Hints for Lustre file system */
    int lustre_overstriping_ratio;

    /* Hints set by PnetCDF internally */
    int lustre_num_osts;
    int *ranklist;

} PNCIO_Hints;

typedef struct {
    MPI_Datatype type;      /* MPI derived datatype */
    MPI_Offset   size;      /* total size in bytes (sum of len[*]) */
    MPI_Count    count;     /* number of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset  *off;       /* [count] byte offsets */
    MPI_Offset  *len;       /* [count] block lengths in bytes */
#else
    MPI_Offset  *off;       /* [count] byte offsets */
    int         *len;       /* [count] block lengths in bytes */
#endif
    MPI_Count    idx;       /* index of off-len pairs consumed so far */
    MPI_Aint     rem;       /* remaining amount in the pair to be consumed */
    int          is_contig; /* whether view of file or buffer is contiguous */
} PNCIO_View;

typedef struct {
    MPI_Comm comm;          /* communicator indicating who called open */
    const char *filename;
    int file_system;        /* type of file system */

    int fd_sys;             /* system file descriptor */
    PNCIO_node_ids node_ids;/* node IDs of each rank */
    int access_mode;        /* Access mode (sequential, append, etc.),
                             * possibly modified to deal with
                             * data sieving or deferred open */

    int is_open;            /* no_indep_rw, 0: not open yet 1: is open */

    int skip_read;          /* whether to skip reads in read-modify-write */

    MPI_Offset disp;        /* file displacement */
    MPI_Datatype filetype;  /* file type set in fileview */
                            /* etype in fileview is always MPI_BYTE in PnetCDF */
    PNCIO_View flat_file;   /* flattern filetype */

    int atomicity;          /* true=atomic, false=nonatomic */
    char *io_buf;           /* two-phase buffer allocated out of i/o path */
    int is_agg;             /* bool: if I am an aggregator */
    int my_cb_nodes_index;  /* my index into fd->hints->ranklist[]. -1 if N/A */
    PNCIO_Hints *hints;     /* structure containing fs-indep. info values */
    MPI_Info info;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double write_timing[NMEASURES];
    double read_timing[NMEASURES];
    MPI_Count write_counter[NMEASURES];
    MPI_Count read_counter[NMEASURES];
#endif
} PNCIO_File;

typedef struct {
    MPI_Offset *offsets;  /* array of offsets */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset *lens;     /* array of lengths */
    MPI_Count  *mem_ptrs; /* array of pointers. used in the read/write phase to
                           * indicate where the data is stored in memory
                           * promoted to MPI_Count so we can construct types
                           * with _c versions
                           */
    MPI_Count   count;    /* size of above arrays */
#else
    int        *lens;
    MPI_Aint   *mem_ptrs;
    size_t      count;
#endif
    size_t curr; /* index of offsets/lens that is currently being processed */
} PNCIO_Access;

/*---- APIs -----------------------------------------------------------------*/
extern
int PNCIO_FileSysType(const char *filename);

extern
int PNCIO_File_open(MPI_Comm comm, const char *filename, int amode,
                MPI_Info info, PNCIO_File *fh);

extern
int PNCIO_File_close(PNCIO_File *fh);

extern
int PNCIO_File_set_view(PNCIO_File *fh, MPI_Offset disp, MPI_Datatype filetype,
                MPI_Aint npairs,
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count *offsets, MPI_Count *lengths
#else
                MPI_Offset *offsets, int *lengths
#endif
);

extern
int PNCIO_File_sync(PNCIO_File *fh);

extern
int PNCIO_File_delete(const char *filename);

extern
int PNCIO_File_set_size(PNCIO_File *fh, MPI_Offset size);

extern
int PNCIO_File_get_size(PNCIO_File *fh, MPI_Offset *size);

extern
int PNCIO_File_get_info(PNCIO_File *fh, MPI_Info *info_used);

extern
int PNCIO_File_SetInfo(PNCIO_File *fh, MPI_Info  users_info);

/* PNC I/O APIs */
extern
MPI_Offset PNCIO_File_write_at(PNCIO_File *fh, MPI_Offset offset,
                const void *buf, PNCIO_View buf_view);
extern
MPI_Offset PNCIO_File_write_at_all(PNCIO_File *fh, MPI_Offset offset,
                const void *buf, PNCIO_View buf_view);

extern
MPI_Offset PNCIO_File_read_at(PNCIO_File *fh, MPI_Offset offset, void *buf,
                PNCIO_View buf_view);
extern
MPI_Offset PNCIO_File_read_at_all(PNCIO_File *fh, MPI_Offset offset, void *buf,
                PNCIO_View buf_view);

extern
MPI_Offset PNCIO_WriteContig(PNCIO_File *fd, const void *buf,
                MPI_Offset w_size, MPI_Offset offset);

extern
MPI_Offset PNCIO_ReadContig(PNCIO_File *fd, void *buf, MPI_Offset r_size,
                MPI_Offset offset);

/* utility APIs */
extern
void PNCIO_Calc_file_domains(MPI_Offset * st_offsets,
                MPI_Offset *end_offsets, int nprocs, int nprocs_for_coll,
                MPI_Offset *min_st_offset_ptr, MPI_Offset **fd_start_ptr,
                MPI_Offset **fd_end_ptr, MPI_Offset *fd_size_ptr,
                int striping_unit);

extern
void PNCIO_Calc_my_req(PNCIO_File *fd, MPI_Offset min_st_offset,
                const MPI_Offset *fd_end, MPI_Offset fd_size,
                int nprocs, MPI_Count *count_my_req_procs_ptr,
                MPI_Count **count_my_req_per_proc_ptr,
                PNCIO_Access **my_req_ptr, MPI_Aint **buf_idx_ptr);

extern
void PNCIO_Calc_others_req(PNCIO_File *fd, MPI_Count count_my_req_procs,
                MPI_Count *count_my_req_per_proc, PNCIO_Access *my_req,
                int nprocs, int myrank, MPI_Count *count_others_req_procs_ptr,
                MPI_Count **count_others_req_per_proc_ptr,
                PNCIO_Access **others_req_ptr);

extern
void PNCIO_Free_my_req(MPI_Count *count_my_req_per_proc,
                PNCIO_Access *my_req, MPI_Aint *buf_idx);

extern
void PNCIO_Free_others_req(MPI_Count *count_others_req_per_proc,
                PNCIO_Access *others_req);


extern
int PNCIO_Calc_aggregator(const PNCIO_File *fd, MPI_Offset off, MPI_Offset min_off,
                MPI_Offset *len, MPI_Offset fd_size, const MPI_Offset *fd_end);

extern
void PNCIO_Heap_merge(PNCIO_Access *others_req, MPI_Count *count,
                MPI_Offset *srt_off, MPI_Count *srt_len, MPI_Count *start_pos,
                int nprocs, int nprocs_recv, MPI_Count total_elements);

/* Generic APIs */
extern
int PNCIO_GEN_SetLock(PNCIO_File *fd, int cmd, int type, MPI_Offset offset,
                int whence, MPI_Offset len);

extern
int PNCIO_GEN_SetLock64(PNCIO_File *fd, int cmd, int type, MPI_Offset offset,
                int whence, MPI_Offset len);

extern
MPI_Offset PNCIO_GEN_WriteStrided(PNCIO_File *fd, const void *buf,
                PNCIO_View buf_view, MPI_Offset offset);

extern
MPI_Offset PNCIO_GEN_ReadStrided_naive(PNCIO_File *fd, void *buf,
                PNCIO_View buf_view, MPI_Offset offset);

extern
MPI_Offset PNCIO_GEN_ReadStridedColl(PNCIO_File *fd, void *buf,
                PNCIO_View buf_view, MPI_Offset offset);

extern
MPI_Offset PNCIO_GEN_WriteStrided_naive(PNCIO_File *fd, const void *buf,
                PNCIO_View buf_view, MPI_Offset offset);

extern
MPI_Offset PNCIO_GEN_ReadStrided(PNCIO_File *fd, void *buf,
                PNCIO_View buf_view, MPI_Offset offset);

extern
MPI_Offset PNCIO_GEN_WriteStridedColl(PNCIO_File *fd, const void *buf,
                PNCIO_View buf_view, MPI_Offset offset);

/* Lustre */
extern
int PNCIO_Lustre_create(PNCIO_File *fd, int access_mode);

extern
int PNCIO_Lustre_open(PNCIO_File *fd);

extern
MPI_Offset PNCIO_LUSTRE_WriteStrided(PNCIO_File *fd, const void *buf,
                PNCIO_View buf_view, MPI_Offset offset);

extern
MPI_Offset PNCIO_LUSTRE_WriteStridedColl(PNCIO_File *fd, const void *buf,
                PNCIO_View buf_view, MPI_Offset offset);

#endif
