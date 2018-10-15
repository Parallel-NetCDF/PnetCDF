/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _PNC_DISPATCH_H
#define _PNC_DISPATCH_H

#include <pnetcdf.h>
#include <mpi.h>

/* various I/O modes or types */
#define NC_REQ_COLL    0x00000001  /* collective request */
#define NC_REQ_INDEP   0x00000002  /* independent request */
#define NC_REQ_WR      0x00000004  /* write request */
#define NC_REQ_RD      0x00000008  /* read request */
#define NC_REQ_ZERO    0x00000010  /* nonblocking API */
#define NC_REQ_HL      0x00000020  /* high-level API */
#define NC_REQ_FLEX    0x00000040  /* flexible API */
#define NC_REQ_BLK     0x00000080  /* blocking get/put API */
#define NC_REQ_NBI     0x00000100  /* nonblocking iget/iput API */
#define NC_REQ_NBB     0x00000200  /* nonblocking bput API */

#define NC_MODE_RDONLY 0x00001000  /* file is opned in read-only mode */
#define NC_MODE_DEF    0x00002000  /* in define mode */
#define NC_MODE_INDEP  0x00004000  /* in independent data mode */
#define NC_MODE_CREATE 0x00008000  /* just created and still in define mode */
#define NC_MODE_FILL   0x00010000  /* fill mode */
#define NC_MODE_SAFE   0x00020000  /* safe mode enabled */
#define NC_MODE_BB     0x00040000  /* burst buffering mode enabled */
#define NC_MODE_SWAP_ON            0x00080000  /* in-place byte swap enabled */
#define NC_MODE_SWAP_OFF           0x00100000  /* in-place byte swap disabled */
#define NC_MODE_STRICT_COORD_BOUND 0x00200000  /* strict coordinate bound check */

/* list of all API kinds */
typedef enum {
    API_VARD,
    API_VARN,
    API_VAR,
    API_VAR1,
    API_VARA,
    API_VARS,
    API_VARM
} NC_api;

struct PNC_driver {
    /* APIs manipulate files */
    int (*create)(MPI_Comm, const char*, int, int, MPI_Info, void**);
    int (*open)(MPI_Comm, const char*, int, int, MPI_Info, void**);
    int (*close)(void*);
    int (*enddef)(void*);
    int (*_enddef)(void*,MPI_Offset,MPI_Offset,MPI_Offset,MPI_Offset);
    int (*redef)(void*);
    int (*sync)(void*);
    int (*flush)(void*);
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
    int (*get_att)(void*,int,const char*,void*,MPI_Datatype);
    int (*put_att)(void*,int,const char*,nc_type,MPI_Offset,const void*,MPI_Datatype);

    /* APIs read/write variables */
    int (*def_var)(void*,const char*,nc_type,int,const int*,int*);
    int (*def_var_fill)(void*,int,int,const void*);
    int (*fill_var_rec)(void*,int,MPI_Offset);
    int (*inq_var)(void*,int,char*,nc_type*,int*,int*,int*,MPI_Offset*,int*,void*);
    int (*inq_varid)(void*,const char*,int*);
    int (*rename_var)(void*,int,const char*);

    int (*get_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,void*,MPI_Offset,MPI_Datatype,int);
    int (*put_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const void*,MPI_Offset,MPI_Datatype,int);

    int (*get_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,void*,MPI_Offset,MPI_Datatype,int);
    int (*put_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,const void*,MPI_Offset,MPI_Datatype,int);

    int (*get_vard)(void*,int,MPI_Datatype,void*,MPI_Offset,MPI_Datatype,int);
    int (*put_vard)(void*,int,MPI_Datatype,const void*,MPI_Offset,MPI_Datatype,int);

    int (*iget_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,void*,MPI_Offset,MPI_Datatype,int*,int);
    int (*iput_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const void*,MPI_Offset,MPI_Datatype,int*,int);
    int (*bput_var)(void*,int,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const MPI_Offset*,const void*,MPI_Offset,MPI_Datatype,int*,int);

    int (*iget_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,void*,MPI_Offset,MPI_Datatype,int*,int);
    int (*iput_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,const void*,MPI_Offset,MPI_Datatype,int*,int);
    int (*bput_varn)(void*,int,int,MPI_Offset* const*,MPI_Offset* const*,const void*,MPI_Offset,MPI_Datatype,int*,int);

    int (*buffer_attach)(void*,MPI_Offset);
    int (*buffer_detach)(void*);
    int (*wait)(void*,int,int*,int*,int);
    int (*cancel)(void*,int,int*,int*);
};

typedef struct PNC_driver PNC_driver;

#define PNC_VARS_CHUNK 64

struct PNC_var {
    int         ndims;
    int         recdim;   /* if >=0, this is a record variable */
    nc_type     xtype;    /* external NC data type */
    MPI_Offset *shape;    /* [ndims] */
};
typedef struct PNC_var PNC_var;

/* one dispatcher object per file: containing info independent from drivers,
 * and can be used for sanity checks, operations need not involve drivers
 */
struct PNC {
    int                mode;        /* file _open/_create mode */
    int                flag;        /* define/data/collective/indep mode */
    int                format;      /* file format */
    char              *path;        /* path name */
    MPI_Comm           comm;        /* MPI communicator */
    int                ndims;       /* number of dimensions defined */
    int                unlimdimid;  /* dim ID of NC_UNLIMITED */
    int                nvars;       /* number of variables */
    int                nrec_vars;   /* number of record variables */
    struct PNC_var    *vars;        /* array of variable objects */
    void              *ncp;         /* pointer to driver internal object */
    struct PNC_driver *driver;
};

typedef struct PNC PNC;

/* subroutine prototypes */

extern PNC_driver* ncmpio_inq_driver(void);

extern PNC_driver* nc4io_inq_driver(void);

extern PNC_driver* ncadios_inq_driver(void); 

extern PNC_driver* ncfoo_inq_driver(void);

extern PNC_driver* ncbbio_inq_driver(void);

extern int PNC_check_id(int ncid, PNC **pncp);

#endif /* _PNC_DISPATCH_H */
