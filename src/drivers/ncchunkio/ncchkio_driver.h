/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _ncchkio_DRIVER_H
#define _ncchkio_DRIVER_H

#include <dispatch.h>
#include <mpi.h>
#include <pnetcdf.h>
#include <ncchk_filter_driver.h>

#include "ncchkioi_profile.h"

#define NC_CHK_VAR_RAW		  0
#define NC_CHK_VAR_COMPRESSED 1
#define NC_CHK_VAR_DATA		  2
#define NC_CHK_VAR_META		  3

#define NC_CHK_MAPPING_STATIC  0
#define NC_CHK_MAPPING_DYNAMIC 01

#define NC_CHK_COMM_CHUNK 0
#define NC_CHK_COMM_PROC  1

#define NC_CHK_ 1

/* Chunk cache structure */
typedef struct NC_chk_cache {
	char *buf;	   // Buffer
	size_t bsize;  // Size in byte
	int serial;	   // batch number to detect swap out of cache allocated in the same batch
	struct NC_chk_cache **ref;	// Ref to clr when it is swap out
	struct NC_chk_cache *prev;
	struct NC_chk_cache *next;
} NC_chk_cache;

/* Get_req structure */
typedef struct NC_chk_req {
	int varid;
	int nreq;
	MPI_Offset *start;
	MPI_Offset **starts;
	MPI_Offset *count;
	MPI_Offset **counts;
	MPI_Offset *stride;
	MPI_Offset bufcount;
	MPI_Datatype buftype;
	char *buf;
	char *xbuf;
	char **xbufs;
} NC_chk_req;

/* Get_req list structure */
typedef struct NC_chk_req_list {
	NC_chk_req *reqs;  // Array of request object
	int *ids;		   // Array of request ids
	int *pos;		   // Array of position of request ids in ids
	int nalloc;		   // Size of the pool
	int nused;		   // Number of ids issued
} NC_chk_req_list;

typedef struct NC_chk_var_chunk {
	MPI_Offset *start;
	MPI_Offset *xdata_offs;
	MPI_Offset *xdata_lens;
	int owner;
	char *data;
	char *xdata;
} NC_chk_var_chunk;

typedef struct NC_chk_chunk_index_entry {
	MPI_Offset off;
	int len;
} NC_chk_chunk_index_entry;

typedef struct NC_chk_var {
	int varkind;
	int isrec;
	int isnew;

	nc_type xtype;
	MPI_Datatype etype;
	int esize;

	int ndim;
	MPI_Offset *dimsize;
	int *dimids;

	int varid;

	int nchunk;
	int nchunkrec;
	int nchunkalloc;
	int nrec;
	int nrecalloc;
	int expanded;
	int chunksize;
	int *nchunks;
	int *cidsteps;
	int *chunk_owner;
	int *chunkdim;
	int *dirty;
	NC_chk_cache **chunk_cache;

	int nmychunk;
	int nmychunkrec;
	int *mychunks;

	MPI_Offset metaoff;
	NC_chk_chunk_index_entry *chunk_index;
	// MPI_Offset *data_offs;
	// int *data_lens;

	NCCHK_filter *filter_driver; /* Compression driver */
	int filter;

	int chunk_map_method;
} NC_chk_var;

typedef struct NC_chk_var_list {
	NC_chk_var *data;
	int cnt;
	int nalloc;
} NC_chk_var_list;

typedef struct NC_chk NC_chk; /* forward reference */
struct NC_chk {
	int mode; /* file _open/_create mode */
	int flag; /* define/data/collective/indep mode */
	int rank;
	int np;
	char *path;	   /* path name */
	MPI_Comm comm; /* MPI communicator */
	void *ncp;	   /* pointer to driver's internal object */
	struct PNC_driver *driver;
	int blockmapping;
	MPI_Offset recsize;	  /* record dim size */
	MPI_Offset recnalloc; /* record dim allocated */
	MPI_Offset default_recnalloc;
	int recdim; /* record dim id */
	NC_chk_var_list vars;
	NC_chk_req_list putlist, getlist;
	int comm_unit;
	int delay_init;
	int exact_cown;
	int max_ndim;
	int max_chunk_size;
	MPI_Offset nmychunks;  // Sum of nmychunk in everyvar
	int default_filter;
	int nwrite;
	MPI_Offset getsize;
	MPI_Offset putsize;
	size_t cache_limit;
	size_t cache_limit_hint;
	size_t cache_used;
	int cache_serial;
	NC_chk_cache *cache_head;
	NC_chk_cache *cache_tail;
	int ndim;				// Number of dim in file
	int *chunkdim;			// Default chunk dim for each dimension
	MPI_Offset cown_size;	// Size of all chunks owned
	MPI_Datatype overlaptype;
    MPI_Op max_cown_op;
	MPI_Offset assigned_chunks;
    double cown_ratio;
	size_t hdr_reserve;	// Additional reserve space in the file header
#ifdef PNETCDF_PROFILING
	NC_chk_timers profile;
	MPI_Offset sendsize;
	MPI_Offset recvsize;
	MPI_Offset var_size_sum;
	MPI_Offset var_zsize_sum;
	int nsend;
	int nrecv;
	int nremote;
	int nreq;
	int nlocal;
#endif
};

extern int ncchkio_create (
	MPI_Comm comm, const char *path, int cmode, int ncid, MPI_Info info, void **ncdp);

extern int ncchkio_open (
	MPI_Comm comm, const char *path, int omode, int ncid, MPI_Info info, void **ncdp);

extern int ncchkio_close (void *ncdp);

extern int ncchkio_enddef (void *ncdp);

extern int ncchkio__enddef (
	void *ncdp, MPI_Offset h_minfree, MPI_Offset v_align, MPI_Offset v_minfree, MPI_Offset r_align);

extern int ncchkio_redef (void *ncdp);

extern int ncchkio_sync (void *ncdp);

extern int ncchkio_flush (void *ncdp);

extern int ncchkio_abort (void *ncdp);

extern int ncchkio_set_fill (void *ncdp, int fill_mode, int *old_fill_mode);

extern int ncchkio_fill_var_rec (void *ncdp, int varid, MPI_Offset recno);

extern int ncchkio_inq (void *ncdp, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

extern int ncchkio_inq_misc (void *ncdp,
							 int *pathlen,
							 char *path,
							 int *num_fix_varsp,
							 int *num_rec_varsp,
							 int *striping_size,
							 int *striping_count,
							 MPI_Offset *header_size,
							 MPI_Offset *header_extent,
							 MPI_Offset *recsize,
							 MPI_Offset *put_size,
							 MPI_Offset *get_size,
							 MPI_Info *info_used,
							 int *nreqs,
							 MPI_Offset *usage,
							 MPI_Offset *buf_size);

extern int ncchkio_sync_numrecs (void *ncdp);

extern int ncchkio_begin_indep_data (void *ncdp);

extern int ncchkio_end_indep_data (void *ncdp);

extern int ncchkio_def_dim (void *ncdp, const char *name, MPI_Offset size, int *dimidp);

extern int ncchkio_inq_dimid (void *ncdp, const char *name, int *dimidp);

extern int ncchkio_inq_dim (void *ncdp, int dimid, char *name, MPI_Offset *lengthp);

extern int ncchkio_rename_dim (void *ncdp, int dimid, const char *newname);

extern int ncchkio_inq_att (
	void *ncdp, int varid, const char *name, nc_type *xtypep, MPI_Offset *lenp);

extern int ncchkio_inq_attid (void *ncdp, int varid, const char *name, int *idp);

extern int ncchkio_inq_attname (void *ncdp, int varid, int attnum, char *name);

extern int ncchkio_copy_att (
	void *ncdp_in, int varid_in, const char *name, void *ncdp_out, int varid_out);

extern int ncchkio_rename_att (void *ncdp, int varid, const char *name, const char *newname);

extern int ncchkio_del_att (void *ncdp, int varid, const char *name);

extern int ncchkio_get_att (
	void *ncdp, int varid, const char *name, void *value, MPI_Datatype itype);

extern int ncchkio_put_att (void *ncdp,
							int varid,
							const char *name,
							nc_type xtype,
							MPI_Offset nelems,
							const void *value,
							MPI_Datatype itype);

extern int ncchkio_def_var (
	void *ncdp, const char *name, nc_type type, int ndims, const int *dimids, int *varidp);

extern int ncchkio_def_var_fill (void *ncdp, int varid, int nofill, const void *fill_value);

extern int ncchkio_inq_var (void *ncdp,
							int varid,
							char *name,
							nc_type *xtypep,
							int *ndimsp,
							int *dimids,
							int *nattsp,
							MPI_Offset *offsetp,
							int *no_fill,
							void *fill_value);

extern int ncchkio_inq_varid (void *ncdp, const char *name, int *varid);

extern int ncchkio_rename_var (void *ncdp, int varid, const char *newname);

extern int ncchkio_get_var (void *ncdp,
							int varid,
							const MPI_Offset *start,
							const MPI_Offset *count,
							const MPI_Offset *stride,
							const MPI_Offset *imap,
							void *buf,
							MPI_Offset bufcount,
							MPI_Datatype buftype,
							int reqMode);

extern int ncchkio_put_var (void *ncdp,
							int varid,
							const MPI_Offset *start,
							const MPI_Offset *count,
							const MPI_Offset *stride,
							const MPI_Offset *imap,
							const void *buf,
							MPI_Offset bufcount,
							MPI_Datatype buftype,
							int reqMode);

extern int ncchkio_get_varn (void *ncdp,
							 int varid,
							 int num,
							 MPI_Offset *const *starts,
							 MPI_Offset *const *counts,
							 void *buf,
							 MPI_Offset bufcount,
							 MPI_Datatype buftype,
							 int reqMode);

extern int ncchkio_put_varn (void *ncdp,
							 int varid,
							 int num,
							 MPI_Offset *const *starts,
							 MPI_Offset *const *counts,
							 const void *buf,
							 MPI_Offset bufcount,
							 MPI_Datatype buftype,
							 int reqMode);

extern int ncchkio_get_vard (void *ncdp,
							 int varid,
							 MPI_Datatype filetype,
							 void *buf,
							 MPI_Offset bufcount,
							 MPI_Datatype buftype,
							 int reqMode);

extern int ncchkio_put_vard (void *ncdp,
							 int varid,
							 MPI_Datatype filetype,
							 const void *buf,
							 MPI_Offset bufcount,
							 MPI_Datatype buftype,
							 int reqMode);

extern int ncchkio_iget_var (void *ncdp,
							 int varid,
							 const MPI_Offset *start,
							 const MPI_Offset *count,
							 const MPI_Offset *stride,
							 const MPI_Offset *imap,
							 void *buf,
							 MPI_Offset bufcount,
							 MPI_Datatype buftype,
							 int *req,
							 int reqMode);

extern int ncchkio_iput_var (void *ncdp,
							 int varid,
							 const MPI_Offset *start,
							 const MPI_Offset *count,
							 const MPI_Offset *stride,
							 const MPI_Offset *imap,
							 const void *buf,
							 MPI_Offset bufcount,
							 MPI_Datatype buftype,
							 int *req,
							 int reqMode);

extern int ncchkio_bput_var (void *ncdp,
							 int varid,
							 const MPI_Offset *start,
							 const MPI_Offset *count,
							 const MPI_Offset *stride,
							 const MPI_Offset *imap,
							 const void *buf,
							 MPI_Offset bufcount,
							 MPI_Datatype buftype,
							 int *req,
							 int reqMode);

extern int ncchkio_iget_varn (void *ncdp,
							  int varid,
							  int num,
							  MPI_Offset *const *starts,
							  MPI_Offset *const *counts,
							  void *buf,
							  MPI_Offset bufcount,
							  MPI_Datatype buftype,
							  int *reqid,
							  int reqMode);

extern int ncchkio_iput_varn (void *ncdp,
							  int varid,
							  int num,
							  MPI_Offset *const *starts,
							  MPI_Offset *const *counts,
							  const void *buf,
							  MPI_Offset bufcount,
							  MPI_Datatype buftype,
							  int *reqid,
							  int reqMode);

extern int ncchkio_bput_varn (void *ncdp,
							  int varid,
							  int num,
							  MPI_Offset *const *starts,
							  MPI_Offset *const *counts,
							  const void *buf,
							  MPI_Offset bufcount,
							  MPI_Datatype buftype,
							  int *reqid,
							  int reqMode);

extern int ncchkio_buffer_attach (void *ncdp, MPI_Offset bufsize);

extern int ncchkio_buffer_detach (void *ncdp);

extern int ncchkio_wait (void *ncdp, int num_reqs, int *req_ids, int *statuses, int reqMode);

extern int ncchkio_cancel (void *ncdp, int num_reqs, int *req_ids, int *statuses);

#endif
