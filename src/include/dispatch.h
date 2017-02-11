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
int (*inq_misc)(void*,int*,char*,int*,int*,int*,int*,MPI_Offset*,MPI_Offset*,MPI_Offset*,MPI_Offset*,MPI_Offset*,MPI_Info*);
int (*sync_numrecs)(void*);
int (*begin_indep_data)(void*);
int (*end_indep_data)(void*);

/* APIs manipulate dimensions */
int (*def_dim)(void*,const char*,MPI_Offset,int*);
int (*inq_dimid)(void*,const char*,int*);
int (*inq_dim)(void*,int,char*,MPI_Offset*);
int (*rename_dim)(void*, int, const char*);

/* APIs read/write attributes */
int (*inq_att)(void*,int,const char*,nc_type*,MPI_Offset*);
int (*inq_attid)(void*,int,const char*,int*);
int (*inq_attname)(void*,int,int,char*);
int (*copy_att)(void*,int,const char*,void*,int);
int (*rename_att)(void*,int,const char*,const char*);
int (*del_att)(void*,int,const char*);
int (*get_att)(void*,int,const char*,void*);
int (*get_att_text)(void*,int,const char*,char*);
int (*get_att_schar)(void*,int,const char*,signed char*);
int (*get_att_uchar)(void*,int,const char*,unsigned char*);
int (*get_att_short)(void*,int,const char*,short*);
int (*get_att_ushort)(void*,int,const char*,unsigned short*);
int (*get_att_int)(void*,int,const char*,int*);
int (*get_att_uint)(void*,int,const char*,unsigned int*);
int (*get_att_long)(void*,int,const char*,long*);
int (*get_att_float)(void*,int,const char*,float*);
int (*get_att_double)(void*,int,const char*,double*);
int (*get_att_longlong)(void*,int,const char*,long long*);
int (*get_att_ulonglong)(void*,int,const char*,unsigned long long*);
int (*put_att)(void*,int,const char*,nc_type,MPI_Offset,const void*);
int (*put_att_text)(void*,int,const char*,MPI_Offset,const char*);
int (*put_att_schar)(void*,int,const char*,nc_type,MPI_Offset,const signed char*);
int (*put_att_uchar)(void*,int,const char*,nc_type,MPI_Offset,const unsigned char*);
int (*put_att_short)(void*,int,const char*,nc_type,MPI_Offset,const short*);
int (*put_att_ushort)(void*,int,const char*,nc_type,MPI_Offset,const unsigned short*);
int (*put_att_int)(void*,int,const char*,nc_type,MPI_Offset,const int*);
int (*put_att_uint)(void*,int,const char*,nc_type,MPI_Offset,const unsigned int*);
int (*put_att_long)(void*,int,const char*,nc_type,MPI_Offset,const long*);
int (*put_att_float)(void*,int,const char*,nc_type,MPI_Offset,const float*);
int (*put_att_double)(void*,int,const char*,nc_type,MPI_Offset,const  double*);
int (*put_att_longlong)(void*,int,const char*,nc_type,MPI_Offset,const long long*);
int (*put_att_ulonglong)(void*,int,const char*,nc_type,MPI_Offset,const unsigned long long*);

/* APIs read/write variables */
int (*def_var)(void*,const char*,nc_type,int,const int*,int*);
int (*def_var_fill)(void*,int,int,const void*);
int (*inq_var)(void*,int,char*,nc_type*,int*,int*,int*,MPI_Offset*,int*,void*);
int (*inq_varid)(void*,const char*,int*);
int (*rename_var)(void*,int,const char*);

#ifdef NOT_YET
int (*get_att)(int, int, const char*, void*, nc_type);
int (*put_att)(int, int, const char*, nc_type, size_t, const void*, nc_type);


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
