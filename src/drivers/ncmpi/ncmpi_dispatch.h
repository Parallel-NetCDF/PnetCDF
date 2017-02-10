/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _NCMPI_DISPATCH_H
#define _NCMPI_DISPATCH_H

#include <mpi.h>

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

#endif
