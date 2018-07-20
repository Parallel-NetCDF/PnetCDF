/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* $Id$ */

#ifndef _NCFOO_DRIVER_H
#define _NCFOO_DRIVER_H

#include <mpi.h>
#include <pnetcdf.h>
#include <dispatch.h>
#include <limits.h>
#include <unistd.h>

#define NC_LOG_TYPE_TEXT 1
#define NC_LOG_TYPE_SCHAR 2
#define NC_LOG_TYPE_UCHAR 3
#define NC_LOG_TYPE_SHORT 4
#define NC_LOG_TYPE_USHORT 5
#define NC_LOG_TYPE_INT 6
#define NC_LOG_TYPE_UINT 7
#define NC_LOG_TYPE_LONG 13
#define NC_LOG_TYPE_FLOAT 8
#define NC_LOG_TYPE_DOUBLE 9
#define NC_LOG_TYPE_LONGLONG 10
#define NC_LOG_TYPE_ULONGLONG 11
#define NC_LOG_TYPE_NATIVE 12

#define NC_LOG_API_KIND_VAR 1
#define NC_LOG_API_KIND_VAR1 2
#define NC_LOG_API_KIND_VARA 3
#define NC_LOG_API_KIND_VARS 4

#define NC_LOG_MAGIC_SIZE 8
#define NC_LOG_MAGIC "PnetCDF0"

#define NC_LOG_FORMAT_SIZE 8
#define NC_LOG_FORMAT_CDF_MAGIC "CDF0\0\0\0\0"
#define NC_LOG_FORMAT_HDF5_MAGIC "\211HDF\r\n\032\n"
#define NC_LOG_FORMAT_BP_MAGIC "BP\0\0\0\0\0\0"

#define NC_LOG_FALSE 0x00
#define NC_LOG_TRUE 0x01

#define NC_LOG_HINT_LOG_ENABLE 0x01
#define NC_LOG_HINT_DEL_ON_CLOSE 0x02
#define NC_LOG_HINT_FLUSH_ON_WAIT 0x04
#define NC_LOG_HINT_FLUSH_ON_SYNC 0x08
#define NC_LOG_HINT_FLUSH_ON_READ 0x10
#define NC_LOG_HINT_LOG_OVERWRITE 0x20
#define NC_LOG_HINT_LOG_CHECK 0x40
#define NC_LOG_HINT_LOG_SHARE 0x80

/* PATH_MAX after padding to 4 byte allignment */
#if PATH_MAX % 4 == 0
#define NC_LOG_PATH_MAX PATH_MAX
#elif PATH_MAX % 4 == 1
#define NC_LOG_PATH_MAX PATH_MAX + 3
#elif PATH_MAX % 4 == 2
#define NC_LOG_PATH_MAX PATH_MAX + 2
#elif PATH_MAX % 4 == 3
#define NC_LOG_PATH_MAX PATH_MAX + 1
#endif

/* Metadata header
 * Variable named according to the spec
 * ToDo: Replace int with 4 byte integer variable if int is not 4 byte
 */
typedef struct NC_bb_metadataheader {
    char magic[NC_LOG_MAGIC_SIZE];
    char format[NC_LOG_FORMAT_SIZE];
    int big_endian;
    int is_external;
    MPI_Offset num_ranks;
    MPI_Offset rank_id;
    MPI_Offset entry_begin;
    MPI_Offset max_ndims;
    MPI_Offset num_entries;
    int basenamelen;
    char basename[1];   /* The hack to keep basename inside the structure */
} NC_bb_metadataheader;

/* Metadata entry header
 * Variable named according to the spec
 * ToDo: Replace int with 4 byte integer variable if int is not 4 byte
 */
typedef struct NC_bb_metadataentry {
    MPI_Offset esize;
    int api_kind;
    int itype;
    int varid;
    int ndims;
    MPI_Offset data_off;
    MPI_Offset data_len;
} NC_bb_metadataentry;

typedef struct NC_bb_metadataptr {
    NC_bb_metadataentry *ptr;
    int valid;
    int reqid;
} NC_bb_metadataptr;

typedef struct NC_bb_metadataidx {
    NC_bb_metadataptr *entries;
    int nused;
    int nalloc;
} NC_bb_metadataidx;

/* Buffer structure */
typedef struct NC_bb_buffer {
    size_t nalloc;
    size_t nused;
    void *buffer;
} NC_bb_buffer;

/* Vector structure */
typedef struct NC_bb_sizevector {
    size_t nalloc;
    size_t nused;
    size_t *values;
} NC_bb_sizevector;

/* Vector structure */
typedef struct NC_bb_intvector {
    size_t nalloc;
    size_t nused;
    int *values;
} NC_bb_intvector;

/* Put_req structure */
typedef struct NC_bb_put_req {
    int valid;  // If this request object is in use (corresponding to some nonblocking request)
    int ready;  // If the corresponding log entry is flushed and the status is avaiable
    int status; // status of flushing the corresponding log entry
    int entrystart; // First log entry generated by this request
    int entryend;   // Last log entry generated by this request
} NC_bb_put_req;

/* Put_req structure */
typedef struct NC_bb_put_list {
    NC_bb_put_req *reqs;    // Array of request object
    int *ids;   // Array of request ids
    int nalloc; // Size of the pool
    int nused;  // Number of ids issued
} NC_bb_put_list;

/* Shared file object */
typedef struct NC_bb_sharedfile {
    int fd; // POSIX file descriptor
    int chanel; // Which chanel are we on (Usually according to the rank)
    int nchanel;    // How many chanel are there (how many process are sharing the file)
    size_t pos; // Logical file position within the fileview
    size_t bsize;   // Dividing blocksize
    size_t fsize;   // Current file size
} NC_bb_sharedfile;

/* Log structure */
typedef struct NC_bb {
    char metalogpath[PATH_MAX];    /* path of metadata log */
    char datalogpath[PATH_MAX];    /* path of data log */
    char logbase[PATH_MAX];        /* path of log files */
    int rank;
    int np;
    NC_bb_sharedfile *metalog_fd;    /* file handle of metadata log */
    NC_bb_sharedfile *datalog_fd;    /* file handle of data log */
    int recdimid;
    int inited;
    int hints;
    size_t datalogsize;
    NC_bb_buffer metadata; /* In memory metadata buffer that mirrors the metadata log */
    NC_bb_metadataidx metaidx;
    NC_bb_sizevector entrydatasize;    /* Array of metadata entries */
    int isflushing;   /* If log is flushing */
    MPI_Offset max_ndims;
    NC_bb_put_list putlist;
    MPI_Offset recdimsize;
    MPI_Offset flushbuffersize;
    MPI_Offset maxentrysize;
#ifdef PNETCDF_PROFILING
    /* Profiling information */
    MPI_Offset total_data;
    MPI_Offset total_meta;
    MPI_Offset max_buffer;
    double total_time;
    double create_time;
    double enddef_time;
    double put_time;
    double flush_time;
    double close_time;
    double flush_replay_time;
    double flush_data_rd_time;
    double flush_put_time;
    double flush_wait_time;
    double put_data_wr_time;
    double put_meta_wr_time;
    double put_num_wr_time;
#endif

    int                mode;        /* file _open/_create mode */
    int                flag;        /* define/data/collective/indep mode */
    int                ncid;
    char              *path;        /* path name */
    MPI_Comm           comm;        /* MPI communicator */
    MPI_Comm           logcomm;        /* MPI communicator */
    MPI_Info           info;
    void              *ncp;         /* pointer to driver's internal object */
    struct PNC_driver *ncmpio_driver;
} NC_bb;

int ncbbio_get_node_comm(MPI_Comm global_comm, MPI_Comm *node_comm);
int ncbbio_log_buffer_init(NC_bb_buffer * bp);
void ncbbio_log_buffer_free(NC_bb_buffer * bp);
char* ncbbio_log_buffer_alloc(NC_bb_buffer *bp, size_t size);
int ncbbio_log_sizearray_init(NC_bb_sizevector *sp);
void ncbbio_log_sizearray_free(NC_bb_sizevector *sp);
int ncbbio_log_sizearray_append(NC_bb_sizevector *sp, size_t size);
int log_flush(NC_bb *ncbbp);
int ncbbio_log_create(NC_bb *ncbbp);
int ncbbio_log_put_var(NC_bb *ncbbp, int varid, const MPI_Offset start[], const MPI_Offset count[], const MPI_Offset stride[], void *buf, MPI_Datatype buftype, MPI_Offset *putsize);
int ncbbio_log_close(NC_bb *ncbbp, int replay);
int ncbbio_log_flush(NC_bb *ncbbp);
int ncbbio_log_enddef(NC_bb *ncbbp);
int ncbbio_init(NC_bb *ncbbp);

int ncbbio_put_list_init(NC_bb *ncbbp);
int ncbbio_put_list_free(NC_bb *ncbbp);
int ncbbio_put_list_add(NC_bb *ncbbp, int *id);
int ncbbio_put_list_remove(NC_bb *ncbbp, int reqid);
int ncbbio_handle_put_req(NC_bb *ncbbp, int reqid, int *stat);
int ncbbio_handle_all_put_req(NC_bb *ncbbp);
int ncbbio_cancel_put_req(NC_bb *ncbbp, int reqid, int *stat);
int ncbbio_cancel_all_put_req(NC_bb *ncbbp);
int ncbbio_metaidx_init(NC_bb *ncbbp);
int ncbbio_metaidx_add(NC_bb *ncbbp, NC_bb_metadataentry *entry);
int ncbbio_metaidx_free(NC_bb *ncbbp);
int ncbbio_log_intvector_init(NC_bb_intvector *vp);
void ncbbio_log_intvector_free(NC_bb_intvector *vp);
int ncbbio_log_intvector_append(NC_bb_intvector *vp, int size);

int ncbbio_sharedfile_open(MPI_Comm comm, char *path, int flag, MPI_Info info, NC_bb_sharedfile **fh);
int ncbbio_sharedfile_close(NC_bb_sharedfile *f);
int ncbbio_sharedfile_pwrite(NC_bb_sharedfile *f, void *buf, size_t count, off_t offset);
int ncbbio_sharedfile_write(NC_bb_sharedfile *f, void *buf, size_t count);
int ncbbio_sharedfile_pread(NC_bb_sharedfile *f, void *buf, size_t count, off_t offset);
int ncbbio_sharedfile_read(NC_bb_sharedfile *f, void *buf, size_t count);
int ncbbio_sharedfile_seek(NC_bb_sharedfile *f, off_t offset, int whence);

void ncbbio_extract_hint(NC_bb *ncbbp, MPI_Info info);
void ncbbio_export_hint(NC_bb *ncbbp, MPI_Info info);

extern int
ncbbio_create(MPI_Comm comm, const char *path, int cmode, int ncid, MPI_Info info, void **ncdp);

extern int
ncbbio_open(MPI_Comm comm, const char *path, int omode, int ncid, MPI_Info info, void **ncdp);

extern int
ncbbio_close(void *ncdp);

extern int
ncbbio_enddef(void *ncdp);

extern int
ncbbio__enddef(void *ncdp, MPI_Offset h_minfree, MPI_Offset v_align, MPI_Offset v_minfree, MPI_Offset r_align);

extern int
ncbbio_redef(void *ncdp);

extern int
ncbbio_sync(void *ncdp);

extern int
ncbbio_flush(void *ncdp);

extern int
ncbbio_abort(void *ncdp);

extern int
ncbbio_set_fill(void *ncdp, int fill_mode, int *old_fill_mode);

extern int
ncbbio_fill_var_rec(void *ncdp, int varid, MPI_Offset recno);

extern int
ncbbio_inq(void *ncdp, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

extern int
ncbbio_inq_misc(void *ncdp, int *pathlen, char *path, int *num_fix_varsp,
               int *num_rec_varsp, int *striping_size, int *striping_count,
               MPI_Offset *header_size, MPI_Offset *header_extent,
               MPI_Offset *recsize, MPI_Offset *put_size, MPI_Offset *get_size,
               MPI_Info *info_used, int *nreqs, MPI_Offset *usage,
               MPI_Offset *buf_size);

extern int
ncbbio_sync_numrecs(void *ncdp);

extern int
ncbbio_begin_indep_data(void *ncdp);

extern int
ncbbio_end_indep_data(void *ncdp);

extern int
ncbbio_def_dim(void *ncdp, const char *name, MPI_Offset size, int *dimidp);

extern int
ncbbio_inq_dimid(void *ncdp, const char *name, int *dimidp);

extern int
ncbbio_inq_dim(void *ncdp, int dimid, char *name, MPI_Offset *lengthp);

extern int
ncbbio_rename_dim(void *ncdp, int dimid, const char *newname);

extern int
ncbbio_inq_att(void *ncdp, int varid, const char *name, nc_type *xtypep, MPI_Offset *lenp);

extern int
ncbbio_inq_attid(void *ncdp, int varid, const char *name, int *idp);

extern int
ncbbio_inq_attname(void *ncdp, int varid, int attnum, char *name);

extern int
ncbbio_copy_att(void *ncdp_in, int varid_in, const char *name, void *ncdp_out, int varid_out);

extern int
ncbbio_rename_att(void *ncdp, int varid, const char *name, const char *newname);

extern int
ncbbio_del_att(void *ncdp, int varid, const char *name);

extern int
ncbbio_get_att(void *ncdp, int varid, const char *name, void *value, MPI_Datatype itype);

extern int
ncbbio_put_att(void *ncdp, int varid, const char *name, nc_type xtype, MPI_Offset nelems, const void *value, MPI_Datatype itype);

extern int
ncbbio_def_var(void *ncdp, const char *name, nc_type type, int ndims, const int *dimids, int *varidp);

extern int
ncbbio_def_var_fill(void *ncdp, int varid, int nofill, const void *fill_value);

extern int
ncbbio_inq_var(void *ncdp, int varid, char *name, nc_type *xtypep, int *ndimsp,
               int *dimids, int *nattsp, MPI_Offset *offsetp, int *no_fill, void *fill_value);

extern int
ncbbio_inq_varid(void *ncdp, const char *name, int *varid);

extern int
ncbbio_rename_var(void *ncdp, int varid, const char *newname);

extern int
ncbbio_get_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncbbio_put_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncbbio_get_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncbbio_put_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncbbio_get_vard(void *ncdp, int varid, MPI_Datatype filetype, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncbbio_put_vard(void *ncdp, int varid, MPI_Datatype filetype, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncbbio_iget_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncbbio_iput_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncbbio_bput_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncbbio_iget_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncbbio_iput_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncbbio_bput_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncbbio_buffer_attach(void *ncdp, MPI_Offset bufsize);

extern int
ncbbio_buffer_detach(void *ncdp);

extern int
ncbbio_wait(void *ncdp, int num_reqs, int *req_ids, int *statuses, int reqMode);

extern int
ncbbio_cancel(void *ncdp, int num_reqs, int *req_ids, int *statuses);

#endif
