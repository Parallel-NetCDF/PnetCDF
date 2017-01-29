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
ncmpii_begin_indep_data(void *ncdp);

extern int
ncmpii_end_indep_data(void *ncdp);

extern int
ncmpii_inq(void *ncdp, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

#endif
