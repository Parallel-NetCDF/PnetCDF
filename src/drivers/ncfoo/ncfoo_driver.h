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

typedef struct NC_foo NC_foo; /* forward reference */
struct NC_foo {
    int                mode;        /* file _open/_create mode */
    int                flag;        /* define/data/collective/indep mode */
    char              *path;        /* path name */
    MPI_Comm           comm;        /* MPI communicator */
    void              *ncp;         /* pointer to driver's internal object */
    struct PNC_driver *driver;
};

extern int
ncfoo_create(MPI_Comm comm, const char *path, int cmode, int ncid, MPI_Info info, void **ncdp);

extern int
ncfoo_open(MPI_Comm comm, const char *path, int omode, int ncid, MPI_Info info, void **ncdp);

extern int
ncfoo_close(void *ncdp);

extern int
ncfoo_enddef(void *ncdp);

extern int
ncfoo__enddef(void *ncdp, MPI_Offset h_minfree, MPI_Offset v_align, MPI_Offset v_minfree, MPI_Offset r_align);

extern int
ncfoo_redef(void *ncdp);

extern int
ncfoo_sync(void *ncdp);

extern int
ncfoo_abort(void *ncdp);

extern int
ncfoo_set_fill(void *ncdp, int fill_mode, int *old_fill_mode);

extern int
ncfoo_fill_var_rec(void *ncdp, int varid, MPI_Offset recno);

extern int
ncfoo_inq(void *ncdp, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

extern int
ncfoo_inq_misc(void *ncdp, int *pathlen, char *path, int *num_fix_varsp,
               int *num_rec_varsp, int *striping_size, int *striping_count,
               MPI_Offset *header_size, MPI_Offset *header_extent,
               MPI_Offset *recsize, MPI_Offset *put_size, MPI_Offset *get_size,
               MPI_Info *info_used, int *nreqs, MPI_Offset *usage,
               MPI_Offset *buf_size);

extern int
ncfoo_sync_numrecs(void *ncdp);

extern int
ncfoo_begin_indep_data(void *ncdp);

extern int
ncfoo_end_indep_data(void *ncdp);

extern int
ncfoo_def_dim(void *ncdp, const char *name, MPI_Offset size, int *dimidp);

extern int
ncfoo_inq_dimid(void *ncdp, const char *name, int *dimidp);

extern int
ncfoo_inq_dim(void *ncdp, int dimid, char *name, MPI_Offset *lengthp);

extern int
ncfoo_rename_dim(void *ncdp, int dimid, const char *newname);

extern int
ncfoo_inq_att(void *ncdp, int varid, const char *name, nc_type *xtypep, MPI_Offset *lenp);

extern int
ncfoo_inq_attid(void *ncdp, int varid, const char *name, int *idp);

extern int
ncfoo_inq_attname(void *ncdp, int varid, int attnum, char *name);

extern int
ncfoo_copy_att(void *ncdp_in, int varid_in, const char *name, void *ncdp_out, int varid_out);

extern int
ncfoo_rename_att(void *ncdp, int varid, const char *name, const char *newname);

extern int
ncfoo_del_att(void *ncdp, int varid, const char *name);

extern int
ncfoo_get_att(void *ncdp, int varid, const char *name, void *value, MPI_Datatype itype);

extern int
ncfoo_put_att(void *ncdp, int varid, const char *name, nc_type xtype, MPI_Offset nelems, const void *value, MPI_Datatype itype);

extern int
ncfoo_def_var(void *ncdp, const char *name, nc_type type, int ndims, const int *dimids, int *varidp);

extern int
ncfoo_def_var_fill(void *ncdp, int varid, int nofill, const void *fill_value);

extern int
ncfoo_inq_var(void *ncdp, int varid, char *name, nc_type *xtypep, int *ndimsp,
               int *dimids, int *nattsp, MPI_Offset *offsetp, int *no_fill, void *fill_value);

extern int
ncfoo_inq_varid(void *ncdp, const char *name, int *varid);

extern int
ncfoo_rename_var(void *ncdp, int varid, const char *newname);

extern int
ncfoo_get_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncfoo_put_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncfoo_get_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncfoo_put_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncfoo_get_vard(void *ncdp, int varid, MPI_Datatype filetype, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncfoo_put_vard(void *ncdp, int varid, MPI_Datatype filetype, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncfoo_iget_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncfoo_iput_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncfoo_bput_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncfoo_iget_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncfoo_iput_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncfoo_bput_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncfoo_buffer_attach(void *ncdp, MPI_Offset bufsize);

extern int
ncfoo_buffer_detach(void *ncdp);

extern int
ncfoo_wait(void *ncdp, int num_reqs, int *req_ids, int *statuses, int reqMode);

extern int
ncfoo_cancel(void *ncdp, int num_reqs, int *req_ids, int *statuses);

#endif
