/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _nczipio_DRIVER_H
#define _nczipio_DRIVER_H

#include <mpi.h>
#include <pnetcdf.h>
#include <dispatch.h>
#include <zip_driver.h>

#define NC_ZIP_VAR_RAW 0
#define NC_ZIP_VAR_COMPRESSED 1
#define NC_ZIP_VAR_DATA 2

#define NC_ZIP_MAPPING_STATIC 0

typedef enum {
    STATIC = 0,           
    DYNAMIC = 1 
} mapping_strategy;

/* Get_req structure */
typedef struct NC_zip_req {
    int varid;
    int nreq;
    MPI_Offset *start;
    MPI_Offset **starts;
    MPI_Offset *count;
    MPI_Offset **counts;
    MPI_Offset *stride;
    char *buf;
    char *xbuf;
    char **xbufs;
} NC_zip_req;

/* Get_req list structure */
typedef struct NC_zip_req_list {
    NC_zip_req *reqs;    // Array of request object
    int *ids;   // Array of request ids
    int *pos;   // Array of position of request ids in ids
    int nalloc; // Size of the pool
    int nused;  // Number of ids issued
} NC_zip_req_list;

typedef struct NC_zip_var_chunk {
    MPI_Offset *start;
    MPI_Offset *xdata_offs;
    MPI_Offset *xdata_lens;
    int owner;
    char *data;
    char *xdata;
} NC_zip_var_chunk;

typedef struct NC_zip_var {
    int varkind;

    nc_type xtype;
    MPI_Datatype etype;
    int esize;

    int ndim;
    MPI_Offset *dimsize;
    int *dimids;
    MPI_Offset recsize;
    
    int varid;

    int nchunks;
    int chunksize;
    int *chunk_owner;
    int *chunkdim;
    char **chunk_cache;

    int nmychunks;
    int *mychunks;

    int datavarid;
    MPI_Offset *data_offs;
    MPI_Offset *data_lens;
    
    mapping_strategy chunk_map_method;
} NC_zip_var;

typedef struct NC_zip_var_list {
    NC_zip_var *data;
    int cnt;
    int nalloc;
} NC_zip_var_list;

typedef struct NC_zip NC_zip; /* forward reference */
struct NC_zip {
    int                mode;        /* file _open/_create mode */
    int                flag;        /* define/data/collective/indep mode */
    int                rank;
    int                  np;
    char              *path;        /* path name */
    MPI_Comm           comm;        /* MPI communicator */
    void              *ncp;         /* pointer to driver's internal object */
    struct PNC_driver *driver;
    NCZIP_driver      *zip;         /* Compression driver */
    int                zipdriver;
    int                blockmapping;
    NC_zip_var_list    vars;
};

extern int
nczipio_create(MPI_Comm comm, const char *path, int cmode, int ncid, MPI_Info info, void **ncdp);

extern int
nczipio_open(MPI_Comm comm, const char *path, int omode, int ncid, MPI_Info info, void **ncdp);

extern int
nczipio_close(void *ncdp);

extern int
nczipio_enddef(void *ncdp);

extern int
nczipio__enddef(void *ncdp, MPI_Offset h_minfree, MPI_Offset v_align, MPI_Offset v_minfree, MPI_Offset r_align);

extern int
nczipio_redef(void *ncdp);

extern int
nczipio_sync(void *ncdp);

extern int
nczipio_flush(void *ncdp);

extern int
nczipio_abort(void *ncdp);

extern int
nczipio_set_fill(void *ncdp, int fill_mode, int *old_fill_mode);

extern int
nczipio_fill_var_rec(void *ncdp, int varid, MPI_Offset recno);

extern int
nczipio_inq(void *ncdp, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

extern int
nczipio_inq_misc(void *ncdp, int *pathlen, char *path, int *num_fix_varsp,
               int *num_rec_varsp, int *striping_size, int *striping_count,
               MPI_Offset *header_size, MPI_Offset *header_extent,
               MPI_Offset *recsize, MPI_Offset *put_size, MPI_Offset *get_size,
               MPI_Info *info_used, int *nreqs, MPI_Offset *usage,
               MPI_Offset *buf_size);

extern int
nczipio_sync_numrecs(void *ncdp);

extern int
nczipio_begin_indep_data(void *ncdp);

extern int
nczipio_end_indep_data(void *ncdp);

extern int
nczipio_def_dim(void *ncdp, const char *name, MPI_Offset size, int *dimidp);

extern int
nczipio_inq_dimid(void *ncdp, const char *name, int *dimidp);

extern int
nczipio_inq_dim(void *ncdp, int dimid, char *name, MPI_Offset *lengthp);

extern int
nczipio_rename_dim(void *ncdp, int dimid, const char *newname);

extern int
nczipio_inq_att(void *ncdp, int varid, const char *name, nc_type *xtypep, MPI_Offset *lenp);

extern int
nczipio_inq_attid(void *ncdp, int varid, const char *name, int *idp);

extern int
nczipio_inq_attname(void *ncdp, int varid, int attnum, char *name);

extern int
nczipio_copy_att(void *ncdp_in, int varid_in, const char *name, void *ncdp_out, int varid_out);

extern int
nczipio_rename_att(void *ncdp, int varid, const char *name, const char *newname);

extern int
nczipio_del_att(void *ncdp, int varid, const char *name);

extern int
nczipio_get_att(void *ncdp, int varid, const char *name, void *value, MPI_Datatype itype);

extern int
nczipio_put_att(void *ncdp, int varid, const char *name, nc_type xtype, MPI_Offset nelems, const void *value, MPI_Datatype itype);

extern int
nczipio_def_var(void *ncdp, const char *name, nc_type type, int ndims, const int *dimids, int *varidp);

extern int
nczipio_def_var_fill(void *ncdp, int varid, int nofill, const void *fill_value);

extern int
nczipio_inq_var(void *ncdp, int varid, char *name, nc_type *xtypep, int *ndimsp,
               int *dimids, int *nattsp, MPI_Offset *offsetp, int *no_fill, void *fill_value);

extern int
nczipio_inq_varid(void *ncdp, const char *name, int *varid);

extern int
nczipio_rename_var(void *ncdp, int varid, const char *newname);

extern int
nczipio_get_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
nczipio_put_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
nczipio_get_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
nczipio_put_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
nczipio_get_vard(void *ncdp, int varid, MPI_Datatype filetype, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
nczipio_put_vard(void *ncdp, int varid, MPI_Datatype filetype, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
nczipio_iget_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
nczipio_iput_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
nczipio_bput_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
nczipio_iget_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
nczipio_iput_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
nczipio_bput_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
nczipio_buffer_attach(void *ncdp, MPI_Offset bufsize);

extern int
nczipio_buffer_detach(void *ncdp);

extern int
nczipio_wait(void *ncdp, int num_reqs, int *req_ids, int *statuses, int reqMode);

extern int
nczipio_cancel(void *ncdp, int num_reqs, int *req_ids, int *statuses);

#endif
