/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _PNC_DISPATCH_H_
#define _PNC_DISPATCH_H_

#include <pnetcdf.h>
#include <mpi.h>

#define INDEP_IO 0
#define COLL_IO  1
#define NONBLOCKING_IO  -1

/* list of all API kinds */
typedef enum {
    API_VARD,
    API_VARN,
    API_VAR,
    API_VAR1,
    API_VARA,
    API_VARS,
    API_VARM
} api_kind;

struct PNC_Dispatch {
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
    int (*inq_misc)(void*,int*,char*,int*,int*,int*,int*,MPI_Offset*,MPI_Offset*,MPI_Offset*,MPI_Offset*,MPI_Offset*,MPI_Info*,int*,MPI_Offset*,MPI_Offset*);
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
    int (*get_att)(void*,int,const char*,void*,nc_type);
    int (*put_att)(void*,int,const char*,nc_type,MPI_Offset,const void*,nc_type);

    /* APIs read/write variables */
    int (*def_var)(void*,const char*,nc_type,int,const int*,int*);
    int (*def_var_fill)(void*,int,int,const void*);
    int (*fill_rec)(void*,int,MPI_Offset);
    int (*inq_var)(void*,int,char*,nc_type*,int*,int*,int*,MPI_Offset*,int*,void*);
    int (*inq_varid)(void*,const char*,int*);
    int (*rename_var)(void*,int,const char*);

    int (*get_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,void*,MPI_Offset,MPI_Datatype,api_kind,nc_type,int);
    int (*put_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const void*,MPI_Offset,MPI_Datatype,api_kind,nc_type,int);

    int (*get_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,void*,MPI_Offset,MPI_Datatype,nc_type,int);
    int (*put_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,const void*,MPI_Offset,MPI_Datatype,nc_type,int);

    int (*get_vard)(void*,int,MPI_Datatype,void*,MPI_Offset,MPI_Datatype,int);
    int (*put_vard)(void*,int,MPI_Datatype,const void*,MPI_Offset,MPI_Datatype,int);

    int (*iget_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,void*,MPI_Offset,MPI_Datatype,int*,api_kind,nc_type);
    int (*iput_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const void*,MPI_Offset,MPI_Datatype,int*,api_kind,nc_type);
    int (*bput_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const void*,MPI_Offset,MPI_Datatype,int*,api_kind,nc_type);

    int (*iget_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,void*,MPI_Offset,MPI_Datatype,int*,nc_type);
    int (*iput_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,const void*,MPI_Offset,MPI_Datatype,int*,nc_type);
    int (*bput_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,const void*,MPI_Offset,MPI_Datatype,int*,nc_type);

    int (*buffer_attach)(void*,MPI_Offset);
    int (*buffer_detach)(void*);
    int (*wait)(void*,int,int*,int*,int);
    int (*cancel)(void*,int,int*,int*);
};

typedef struct PNC_Dispatch PNC_Dispatch;

/* Common Shared Structure for all Dispatched Objects */
struct PNC {
    int                  mode;   /* file _open/_create mode */
    int                  format; /* file format */
    char                *path;   /* path name */
    struct PNC_Dispatch *dispatch;
    void                *ncp;    /* pointer to dispatcher data object */
};

typedef struct PNC PNC;

/* subroutine prototypes */

extern PNC_Dispatch* ncmpii_inq_dispatcher(void);

extern int PNC_check_id(int ncid, PNC **pncp);

#endif /* _PNC_DISPATCH_H_ */
