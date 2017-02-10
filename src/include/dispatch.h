/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _DISPATCH_H_
#define _DISPATCH_H_

#include <pnetcdf.h>
#include <mpi.h>


typedef struct PNC_Dispatch PNC_Dispatch;

// extern PNC_Dispatch* ncmpi_dispatcher;
extern PNC_Dispatch* ncmpii_inq_dispatcher(void);

struct PNC_Dispatch {

// int model; /* one of the NC_FORMATX #'s */

/* APIs manipulate files */
int (*create)(MPI_Comm, const char*, int, MPI_Info, void**);
int (*open)(MPI_Comm, const char*, int, MPI_Info, void**);
int (*close)(void*);
int (*enddef)(void*);
int (*_enddef)(void*,MPI_Offset,MPI_Offset,MPI_Offset,MPI_Offset);
int (*redef)(void*);
int (*sync)(void*);
int (*abort)(void*);
int (*set_fill)(void*,int,int*);
int (*inq)(void*,int*,int*,int*,int*);
int (*inq_striping)(void*,int*,int*);
int (*inq_num_rec_vars)(void*,int*);
int (*inq_num_fix_vars)(void*,int*);
int (*inq_recsize)(void*,MPI_Offset*);
int (*inq_put_size)(void*,MPI_Offset*);
int (*inq_get_size)(void*,MPI_Offset*);
int (*inq_header_size)(void*,MPI_Offset*);
int (*inq_header_extent)(void*,MPI_Offset*);
int (*inq_file_info)(void*,MPI_Info*);
int (*sync_numrecs)(void*);
int (*begin_indep_data)(void*);
int (*end_indep_data)(void*);

/* APIs manipulate dimensions */
int (*def_dim)(void*,const char*,MPI_Offset,int*);
int (*inq_dimid)(void*,const char*,int*);
int (*inq_dim)(void*,int,char*,MPI_Offset*);
int (*rename_dim)(void*, int, const char*);

/* APIs read/write variables */

/* APIs read/write attributes */
int (*inq_att)(void*,int,const char*,nc_type*,MPI_Offset*);
int (*inq_attid)(void*,int,const char*,int*);
int (*inq_attname)(void*,int,int,char*);
int (*copy_att)(void*,int,const char*,void*,int);
int (*rename_att)(void*,int,const char*,const char*);
int (*del_att)(void*,int,const char*);

#ifdef NOT_YET
int (*get_att)(int, int, const char*, void*, nc_type);
int (*put_att)(int, int, const char*, nc_type, size_t, const void*, nc_type);

int (*def_var)(int, const char*, nc_type, int, const int*, int*);
int (*inq_varid)(int, const char*, int*);
int (*rename_var)(int, int, const char*);

int (*get_vara)(int, int, const size_t*, const size_t*, void*, nc_type);
int (*put_vara)(int, int, const size_t*, const size_t*, const void*, nc_type);

/* Added to solve Ferret performance problem with Opendap */
int (*get_vars)(int, int, const size_t*, const size_t*, const ptrdiff_t*, void*, nc_type);
int (*put_vars)(int, int, const size_t*, const size_t*, const ptrdiff_t*, const void*, nc_type);

int (*get_varm)(int, int, const size_t*, const size_t*, const ptrdiff_t*, const ptrdiff_t*, void*, nc_type);
int (*put_varm)(int, int, const size_t*, const size_t*, const ptrdiff_t*, const ptrdiff_t*, const void*, nc_type);

int (*inq_var_all)(int ncid, int varid, char *name, nc_type *xtypep,
               int *ndimsp, int *dimidsp, int *nattsp,
               int *shufflep, int *deflatep, int *deflate_levelp,
               int *fletcher32p, int *contiguousp, size_t *chunksizesp,
               int *no_fill, void *fill_valuep, int *endiannessp,
	       int *options_maskp, int *pixels_per_blockp);

int (*var_par_access)(int, int, int);

int (*def_var_fill)(int, int, int, const void*);

#endif
};

#ifdef NOT_YET
/* Following functions must be handled as non-dispatch */
const char* (*nc_inq_libvers)(void);
const char* (*nc_strerror)(int);
int (*nc_delete)(const char*path);

/* Define the common fields for NC and NC_FILE_INFO_T etc */
typedef struct NCcommon {
	int ext_ncid; /* uid << 16 */
	int int_ncid; /* unspecified other id */
	struct NC_Dispatch* dispatch;
	void* dispatchdata; /* per-protocol instance data */
	char* path; /* as specified at open or create */
} NCcommon;
#endif

/* Common Shared Structure for all Dispatched Objects */
typedef struct PNC {
    int   mode; /* as provided to _open/_create */
    int   format; /* file format */
    char *path;
    struct PNC_Dispatch *dispatch;
    void *ncp; /*per-'file' data; points to e.g. NC3_INFO data*/
} PNC;

int PNC_check_id(int ncid, PNC **pncp);

#endif /* _DISPATCH_H */
