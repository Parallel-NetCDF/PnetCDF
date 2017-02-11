/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _NCMPI_DISPATCH_H
#define _NCMPI_DISPATCH_H

#include <mpi.h>
#include <pnetcdf.h>

extern int
ncmpii_create(MPI_Comm comm, const char *path, int cmode, MPI_Info info, void **ncdp);

extern int
ncmpii_open(MPI_Comm comm, const char *path, int omode, MPI_Info info, void **ncdp);

extern int
ncmpii_close(void *ncdp);

extern int
ncmpii_enddef(void *ncdp);

extern int
ncmpii__enddef(void *ncdp, MPI_Offset h_minfree, MPI_Offset v_align, MPI_Offset v_minfree, MPI_Offset r_align);

extern int
ncmpii_redef(void *ncdp);

extern int
ncmpii_sync(void *ncdp);

extern int
ncmpii_abort(void *ncdp);

extern int
ncmpii_set_fill(void *ncdp, int fill_mode, int *old_fill_mode);

extern int
ncmpii_inq(void *ncdp, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

extern int
ncmpii_inq_striping(void *ncdp, int *striping_size, int *striping_count);

extern int
ncmpii_inq_num_rec_vars(void *ncdp, int *num_rec_varsp);

extern int
ncmpii_inq_num_fix_vars(void *ncdp, int *num_fix_varsp);

extern int
ncmpii_inq_recsize(void *ncdp, MPI_Offset *recsize);

extern int
ncmpii_inq_put_size(void *ncdp, MPI_Offset *size);

extern int
ncmpii_inq_get_size(void *ncdp, MPI_Offset *size);

extern int
ncmpii_inq_header_size(void *ncdp, MPI_Offset *size);

extern int
ncmpii_inq_header_extent(void *ncdp, MPI_Offset *extent);

extern int
ncmpii_inq_file_info(void *ncdp, MPI_Info *info);

extern int
ncmpii_sync_numrecs(void *ncdp);

extern int
ncmpii_begin_indep_data(void *ncdp);

extern int
ncmpii_end_indep_data(void *ncdp);

extern int
ncmpii_def_dim(void *ncdp, const char *name, MPI_Offset size, int *dimidp);

extern int
ncmpii_inq_dimid(void *ncdp, const char *name, int *dimidp);

extern int
ncmpii_inq_dim(void *ncdp, int dimid, char *name, MPI_Offset *lengthp);

extern int
ncmpii_rename_dim(void *ncdp, int dimid, const char *newname);

extern int
ncmpii_inq_att(void *ncdp, int varid, const char *name, nc_type *xtypep, MPI_Offset *lenp);

extern int
ncmpii_inq_attid(void *ncdp, int varid, const char *name, int *idp); 

extern int
ncmpii_inq_attname(void *ncdp, int varid, int attnum, char *name);

extern int
ncmpii_copy_att(void *ncdp_in, int varid_in, const char *name, void *ncdp_out, int varid_out);

extern int
ncmpii_rename_att(void *ncdp, int varid, const char *name, const char *newname);

extern int
ncmpii_del_att(void *ncdp, int varid, const char *name);

extern int
ncmpii_def_var(void *ncdp, const char *name, nc_type type, int ndims, const int *dimids, int *varidp);

extern int
ncmpii_inq_var(void *ncdp, int varid, char *name, nc_type *xtypep, int *ndimsp,
               int *dimids, int *nattsp, MPI_Offset *offsetp, int *no_fill, void *fill_value);

extern int
ncmpii_inq_varid(void *ncdp, const char *name, int *varid);

extern int
ncmpii_rename_var(void *ncdp, int varid, const char *newname);

#endif
