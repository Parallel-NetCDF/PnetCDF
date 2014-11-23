/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */
#ifndef _NC_H_
#define _NC_H_

/*
 * netcdf library 'private' data structures, objects and interfaces
 */

#include <stddef.h>     /* size_t */
#include <sys/types.h>  /* off_t */

#include "ncio.h"       /* ncio */
#include "fbits.h"

/* C macros for TRACE MPI calls */
#ifdef PNETCDF_TRACE_MPI_COMM
#define TRACE_COMM(x) printf("TRACE-MPI-COMM: FILE %s FUNC %s() LINE %d calling %s()\n",__FILE__,__func__,__LINE__,#x),mpireturn=x
#else
#define TRACE_COMM(x) mpireturn=x
#endif

#ifdef PNETCDF_TRACE_MPI_IO
#define TRACE_IO(x) printf("TRACE-MPI-IO:   FILE %s FUNC %s() LINE %d calling %s()\n",__FILE__,__func__,__LINE__,#x),mpireturn=x
#else
#define TRACE_IO(x) mpireturn=x
#endif


/* define MPI_OFFSET if not defined */
#ifndef HAVE_MPI_OFFSET_DATATYPE
    #ifdef HAVE_MPI_LONG_LONG_INT
        #define MPI_OFFSET MPI_LONG_LONG_INT
    #else
        #define MPI_OFFSET MPI_INT
    #endif
#endif

/* XXX: this seems really low.  do we end up spending a ton of time mallocing?
 * could we reduce that by increasing this to something 21st century? */
#ifndef NC_ARRAY_GROWBY
#define NC_ARRAY_GROWBY 4
#endif

/* ncmpi_create/ncmpi_open set up header to be 'chunksize' big and to grow
 * by 'chunksize' as new items adde. This used to be 4k. 256k lets us read
 * in an entire climate header in one go */
#define NC_DEFAULT_CHUNKSIZE 256*1024

/* when variabe's nctype is NC_CHAR, I/O buffer's MPI type must be MPI_CHAR
 * and vice versa */
#define NCMPII_ECHAR(nctype, mpitype) ((((nctype) == NC_CHAR) == ((mpitype) != MPI_CHAR)) ? NC_ECHAR : NC_NOERR)

/*
 * The extern size of an empty
 * netcdf version 1 file.
 * The initial value of ncp->xsz.
 */
#define MIN_NC_XSZ 32

/* netcdf file format:
     netcdf_file  = header  data
     header       = magic  numrecs  dim_list  gatt_list  var_list
     magic        = 'C'  'D'  'F'  VERSION
     VERSION      = \x01 | \x02 | \x05
     numrecs      = NON_NEG | STREAMING
     dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     gatt_list    = att_list
     att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     var_list     = ABSENT | NC_VARIABLE   nelems  [var ...]
     ABSENT       = ZERO  ZERO                  // Means list is not present
     ZERO         = \x00 \x00 \x00 \x00         // 32-bit zero

  Minmum happens when nothing is defined, i.e.
     magic              -- 4 bytes
     numrecs            -- 4 bytes for CDF-1 and CDF-2, 8 bytes for CDF-5
     dim_list = ABSENT  -- 8 bytes
     gatt_list = ABSENT -- 8 bytes
     var_list = ABSENT  -- 8 bytes
*/

typedef struct NC NC; /* forward reference */

/*
 *  The internal data types
 */
typedef enum {
    NC_UNSPECIFIED =  0,
/*  NC_BITFIELD    =  7, */
/*  NC_STRING      =  8, */
    NC_DIMENSION   = 10,
    NC_VARIABLE    = 11,
    NC_ATTRIBUTE   = 12
} NCtype;


/*
 * Counted string for names and such
 */
typedef struct {
    /* all xdr'd */
    MPI_Offset  nchars;
    char       *cp;     /* [nchars+1] one additional char for '\0' */
} NC_string;

extern NC *
ncmpii_new_NC(const MPI_Offset *chunkp);

extern NC *
ncmpii_dup_NC(const NC *ref);

/* Begin defined in string.c */
extern void
ncmpii_free_NC_string(NC_string *ncstrp);

extern int
ncmpii_NC_check_name(const char *name, int file_ver);

extern NC_string *
ncmpii_new_NC_string(MPI_Offset slen, const char *str);

extern int
ncmpii_set_NC_string(NC_string *ncstrp, const char *str);

/* End defined in string.c */

/*
 * NC dimension stucture
 */
typedef struct {
    /* all xdr'd */
    NC_string *name;
    MPI_Offset size;
#ifdef ENABLE_SUBFILING
    int range[2]; /* subfile range {start, end} */
    MPI_Offset rcount; /* subfile range count */
    int num_subfiles;
#endif
} NC_dim;

/* the dimension ID returned from ncmpi_def_dim() is an integer pointer
 * which means the total number of defined dimension allowed in a file
 * is up to 2^31-1. Thus, the member ndefined below should be of type int.
 */
typedef struct NC_dimarray {
    int      nalloc;    /* number allocated >= ndefined */
    int      ndefined;  /* number of defined dimensions */
    NC_dim **value;
} NC_dimarray;

/* Begin defined in dim.c */

extern void
ncmpii_free_NC_dim(NC_dim *dimp);

extern NC_dim *
ncmpii_new_x_NC_dim(NC_string *name);

extern int
ncmpii_find_NC_Udim(const NC_dimarray *ncap, NC_dim **dimpp);

extern int
incr_NC_dimarray(NC_dimarray *ncap, NC_dim *newdimp);

extern NC_dim*
dup_NC_dim(const NC_dim *dimp);

/* dimarray */

extern void
ncmpii_free_NC_dimarray(NC_dimarray *ncap);

extern int
ncmpii_dup_NC_dimarray(NC_dimarray *ncap, const NC_dimarray *ref);

extern NC_dim *
ncmpii_elem_NC_dimarray(const NC_dimarray *ncap, int elem);

extern int
ncmpi_def_dim(int ncid, const char *name, MPI_Offset size, int *dimidp);

extern int
ncmpi_rename_dim( int ncid, int dimid, const char *newname);

extern int
ncmpi_inq_dimid(int ncid, const char *name, int *dimid_ptr);

extern int
ncmpi_inq_dim(int ncid, int dimid, char *name, MPI_Offset *sizep);

extern int
ncmpi_inq_dimname(int ncid, int dimid, char *name);

extern int
ncmpi_inq_dimlen(int ncid, int dimid, MPI_Offset *lenp);
/* End defined in dim.c */

/*
 * NC attribute
 */
typedef struct {
    MPI_Offset xsz;      /* amount of space at xvalue (4-byte aligned) */
    NC_string *name;     /* name of the attributes */
    nc_type    type;     /* the discriminant */
    MPI_Offset nelems;   /* number of attribute elements */
    void      *xvalue;   /* the actual data, in external representation */
} NC_attr;

typedef struct NC_attrarray {
    int       nalloc;    /* number allocated >= ndefined */
    int       ndefined;  /* number of defined attributes */
    NC_attr **value;
} NC_attrarray;

/* Begin defined in attr.c */

extern void
ncmpii_free_NC_attr(NC_attr *attrp);

extern NC_attr *
ncmpii_new_x_NC_attr(NC_string *strp, nc_type type, MPI_Offset nelems);

extern int
incr_NC_attrarray(NC_attrarray *ncap, NC_attr *newelemp);

extern NC_attr*
dup_NC_attr(const NC_attr *rattrp);

/* attrarray */

extern void
ncmpii_free_NC_attrarray(NC_attrarray *ncap);

extern int
ncmpii_dup_NC_attrarray(NC_attrarray *ncap, const NC_attrarray *ref);

extern NC_attr *
ncmpii_elem_NC_attrarray(const NC_attrarray *ncap, MPI_Offset elem);

extern int
ncmpi_put_att_text(int ncid, int varid, const char *name,
        MPI_Offset nelems, const char *value);

extern int
ncmpi_get_att_text(int ncid, int varid, const char *name, char *str);

extern int
ncmpi_put_att_schar(int ncid, int varid, const char *name,
        nc_type type, MPI_Offset nelems, const signed char *value);

extern int
ncmpi_get_att_schar(int ncid, int varid, const char *name, signed char *tp);

extern int
ncmpi_put_att_uchar(int ncid, int varid, const char *name,
        nc_type type, MPI_Offset nelems, const unsigned char *value);

extern int
ncmpi_get_att_uchar(int ncid, int varid, const char *name, unsigned char *tp);

extern int
ncmpi_put_att_short(int ncid, int varid, const char *name,
        nc_type type, MPI_Offset nelems, const short *value);

extern int
ncmpi_get_att_short(int ncid, int varid, const char *name, short *tp);

extern int
ncmpi_put_att_int(int ncid, int varid, const char *name,
        nc_type type, MPI_Offset nelems, const int *value);

extern int
ncmpi_get_att_int(int ncid, int varid, const char *name, int *tp);

extern int
ncmpi_put_att_long(int ncid, int varid, const char *name,
        nc_type type, MPI_Offset nelems, const long *value);

extern int
ncmpi_get_att_long(int ncid, int varid, const char *name, long *tp);

extern int
ncmpi_put_att_float(int ncid, int varid, const char *name,
        nc_type type, MPI_Offset nelems, const float *value);
extern int
ncmpi_get_att_float(int ncid, int varid, const char *name, float *tp);
extern int
ncmpi_put_att_double(int ncid, int varid, const char *name,
        nc_type type, MPI_Offset nelems, const double *value);
extern int
ncmpi_get_att_double(int ncid, int varid, const char *name, double *tp);

extern int
ncmpi_inq_attid(int ncid, int varid, const char *name, int *attnump);

extern int
ncmpi_inq_atttype(int ncid, int varid, const char *name, nc_type *datatypep);

extern int
ncmpi_inq_attlen(int ncid, int varid, const char *name, MPI_Offset *lenp);

extern int
ncmpi_inq_att(int ncid, int varid, const char *name,
        nc_type *datatypep, MPI_Offset *lenp);

extern int
ncmpi_copy_att(int ncid_in, int varid_in, const char *name,
        int ncid_out, int ovarid);

extern int
ncmpi_rename_att( int ncid, int varid, const char *name, const char *newname);

extern int
ncmpi_del_att(int ncid, int varid, const char *name);

extern int
ncmpi_inq_attname(int ncid, int varid, int attnum, char *name);
/* End defined in attr.c */

/*
 * NC variable: description and data
 */
typedef struct {
    MPI_Offset    xsz;    /* byte size of 1 array element */
    MPI_Offset   *shape;  /* dim->size of each dim */
    MPI_Offset   *dsizes; /* the right to left product of shape */
    NC_string    *name;   /* name of the variable */
    int           ndims;  /* number of dimensions */
#ifdef ENABLE_SUBFILING
    int           ndims_org; /* ndims before subfiling */
#endif
    int          *dimids; /* array of dimension IDs */
#ifdef ENABLE_SUBFILING
    int          *dimids_org; /* dimids before subfiling */
#endif
    NC_attrarray  attrs;  /* attribute array */
    nc_type       type;   /* variable's data type */
    MPI_Offset    len;    /* this is the "vsize" defined in header format, the
                             total size in bytes of the array variable.
                             For record variable, this is the record size */
    MPI_Offset    begin;  /* starting file offset of this variable */
#ifdef ENABLE_SUBFILING
    int           num_subfiles;
#endif
} NC_var;

/* note: we only allow less than 2^31-1 variables defined in a file */
typedef struct NC_vararray {
    int      nalloc;      /* number allocated >= ndefined */
    int      ndefined;    /* number of defined variables */
    int      num_rec_vars;/* number of defined record variables */
    NC_var **value;
} NC_vararray;

/* Begin defined in var.c */

extern void
ncmpii_free_NC_var(NC_var *varp);

extern NC_var *
ncmpii_new_x_NC_var(NC_string *strp, int ndims);

extern NC_var*
dup_NC_var(const NC_var *rvarp);

extern int
incr_NC_vararray(NC_vararray *ncap, NC_var *newvarp);

/* vararray */

extern void
ncmpii_free_NC_vararray(NC_vararray *ncap);

extern int
ncmpii_dup_NC_vararray(NC_vararray *ncap, const NC_vararray *ref);

extern int
ncmpii_NC_var_shape64(NC *ncp, NC_var *varp, const NC_dimarray *dims);

extern int
ncmpii_NC_check_vlen(NC_var *varp, MPI_Offset vlen_max);

extern NC_var *
ncmpii_NC_lookupvar(NC *ncp, int varid);

extern int
ncmpi_def_var(int ncid, const char *name, nc_type type,
        int ndims, const int *dimidsp, int *varidp);

extern int
ncmpi_rename_var(int ncid, int varid, const char *newname);

extern int
ncmpi_inq_var(int ncid, int varid, char *name, nc_type *typep,
        int *ndimsp, int *dimids, int *nattsp);

extern int
ncmpi_inq_varid(int ncid, const char *name, int *varid_ptr);

extern int
ncmpi_inq_varname(int ncid, int varid, char *name);

extern int
ncmpi_inq_vartype(int ncid, int varid, nc_type *typep);

extern int
ncmpi_inq_varndims(int ncid, int varid, int *ndimsp);

extern int
ncmpi_inq_vardimid(int ncid, int varid, int *dimids);

extern int
ncmpi_inq_varnatts(int ncid, int varid, int *nattsp);

extern int
ncmpi_rename_var(int ncid, int varid, const char *newname);
/* End defined in var.c */

#define IS_RECVAR(vp) \
        ((vp)->shape != NULL ? (*(vp)->shape == NC_UNLIMITED) : 0 )

/*
 *  *  The PnetCDF non-blocking I/O request type
 *   */
typedef struct NC_req {
    int            id;
    int            rw_flag;
    NC_var        *varp;
    void          *buf;    /* the original user buffer */
    void          *xbuf;   /* the buffer used to read/write, may point to
                              the same address as buf */
    int            iscontig_of_ptypes;
    int            need_swap_back_buf;
    int            use_abuf; /* whether use the attached buffer */
    int            is_imap;
    void          *tmpBuf;       /* callback tmp buffer to be freed */
    int            tmpBufSize;   /* size of tmp buffer */
    void          *userBuf;      /* callback user buffer to be unpacked from tmpBuf */
    int            userBufCount; /* user buffer count */
    MPI_Datatype   userBufType;  /* user buffer type */
    int            ndims;
    MPI_Offset    *start;  /* [ndims] */
    MPI_Offset    *count;  /* [ndims] */
    MPI_Offset    *stride; /* [ndims] */
    MPI_Offset     bnelems;
    MPI_Offset     fnelems; /* number of nc_type elements in file */
    MPI_Offset     lnelems;
    MPI_Offset     bufcount;
    MPI_Offset     offset_start;  /* starting of aggregate access region */
    MPI_Offset     offset_end;    /*   ending of aggregate access region */
    MPI_Datatype   buftype;
    MPI_Datatype   ptype;
    MPI_Datatype   imaptype;
    int           *status;
    int            num_subreqs;
    struct NC_req *subreqs;  /* [num_subreq] */
    struct NC_req *next;
} NC_req;

#define NC_ABUF_DEFAULT_TABLE_SIZE 128

typedef struct NC_buf_status {
    int        is_used;
    MPI_Offset req_size;
} NC_buf_status;

typedef struct NC_buf {
    MPI_Offset     size_allocated;
    MPI_Offset     size_used;
    int            table_size;
    NC_buf_status *occupy_table; /* [table_size] */
    int            tail;         /* index of last free entry */
    void          *buf;
} NC_buf;

struct NC {
    /* links to make list of open netcdf's */
    struct NC *next;
    struct NC *prev;
#ifdef ENABLE_SUBFILING
    int nc_num_subfiles; /* # of subfiles */
    int ncid_sf; /* ncid of subfile */
#endif
    /* contains the previous NC during redef. */
    struct NC *old;
    /* flags */
#define NC_INDEP  0x10000   /* in independent data mode, cleared by endindep */
#define NC_CREAT  0x20000   /* in create phase, cleared by ncenddef */
#define NC_INDEF  0x80000   /* in define mode, cleared by ncenddef */
#define NC_NSYNC  0x100000  /* synchronise numrecs on change */
#define NC_HSYNC  0x200000  /* synchronise whole header on change */
#define NC_NDIRTY 0x400000  /* numrecs has changed */
#define NC_HDIRTY 0x800000  /* header info has changed */
/* NC_NOFILL is defined in netcdf.h, historical interface */
    int           flags;
    int           safe_mode;    /* 0 or 1, for parameter consistentcy check */
    int           subfile_mode; /* 0 or 1, for disable/enable subfiling */
    ncio         *nciop;
    MPI_Offset    chunk;    /* largest extent this layer will request from
                               ncio->get() */
    MPI_Offset    xsz;      /* external size of this header, <= var[0].begin */
    MPI_Offset    begin_var;/* file offset of the first (non-record) var */
    MPI_Offset    begin_rec;/* file offset of the first 'record' */

    MPI_Offset    recsize;  /* length of 'record': sum of single record sizes
                               of all the record variables */
    MPI_Offset    numrecs;  /* number of 'records' allocated */
    NC_dimarray   dims;     /* dimensions defined */
    NC_attrarray  attrs;    /* global attributes defined */
    NC_vararray   vars;     /* variables defined */
    NC_req       *head;     /* linked list head of nonblocking requests */
    NC_req       *tail;     /* tail of the linked list */
    NC_buf       *abuf;     /* attached buffer, used by bput APIs */
};

#define NC_readonly(ncp) \
        (!fIsSet(ncp->nciop->ioflags, NC_WRITE))

#define NC_IsNew(ncp) \
        fIsSet((ncp)->flags, NC_CREAT)

#define NC_indep(ncp) \
        fIsSet((ncp)->flags, NC_INDEP)

#define NC_indef(ncp) \
        (NC_IsNew(ncp) || fIsSet((ncp)->flags, NC_INDEF))

#define set_NC_ndirty(ncp) \
        fSet((ncp)->flags, NC_NDIRTY)

#define NC_ndirty(ncp) \
        fIsSet((ncp)->flags, NC_NDIRTY)

#define set_NC_hdirty(ncp) \
        fSet((ncp)->flags, NC_HDIRTY)

#define NC_hdirty(ncp) \
        fIsSet((ncp)->flags, NC_HDIRTY)

#define NC_dofill(ncp) \
        (!fIsSet((ncp)->flags, NC_NOFILL))

#define NC_doHsync(ncp) \
        fIsSet((ncp)->flags, NC_HSYNC)

#define NC_doNsync(ncp) \
        fIsSet((ncp)->flags, NC_NSYNC)

#define NC_get_numrecs(ncp) \
        ((ncp)->numrecs)

#define NC_set_numrecs(ncp, nrecs) \
        {((ncp)->numrecs = (nrecs));}

#define NC_increase_numrecs(ncp, nrecs) \
        {if((nrecs) > (ncp)->numrecs) ((ncp)->numrecs = (nrecs));}

#define ErrIsHeaderDiff(err) \
        (NC_EMULTIDEFINE_FIRST >= err && err >= NC_EMULTIDEFINE_LAST)


/* Begin defined in nc.c */

extern int
ncmpii_NC_check_id(int ncid, NC **ncpp);

extern int
ncmpii_cktype(int cdf_ver, nc_type datatype);

extern MPI_Offset
ncmpix_howmany(nc_type type, MPI_Offset xbufsize);

extern int
ncmpii_dset_has_recvars(NC *ncp);

extern int
ncmpii_read_numrecs(NC *ncp);

extern int
ncmpii_write_numrecs(NC *ncp, MPI_Offset new_numrecs, int forceWrite);

extern int
ncmpii_write_header(NC *ncp);

extern int
ncmpii_NC_sync(NC *ncp, int doFsync);

extern void
ncmpii_free_NC(NC *ncp);

extern void
ncmpii_add_to_NCList(NC *ncp);

extern void
ncmpii_del_from_NCList(NC *ncp);

extern int
ncmpii_read_NC(NC *ncp);

extern int
ncmpii_enddef(NC *ncp);

extern int
ncmpii__enddef(NC *ncp, MPI_Offset h_minfree, MPI_Offset v_align,
               MPI_Offset v_minfree, MPI_Offset r_align);

extern int
ncmpii_NC_close(NC *ncp);

extern int
ncmpi_inq(int ncid, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

extern int
ncmpi_inq_ndims(int ncid, int *ndimsp);

extern int
ncmpi_inq_nvars(int ncid, int *nvarsp);

extern int
ncmpi_inq_natts(int ncid, int *nattsp);

extern int
ncmpi_inq_unlimdim(int ncid, int *xtendimp);

extern int
ncmpi_get_default_format(void);

extern int
ncmpi_inq_num_rec_vars(int ncid, int *nump);

extern int
ncmpi_inq_num_fix_vars(int ncid, int *nump);

/* End defined in nc.c */

#if 0
/* Begin defined in v1hpg.c */

extern size_t
ncx_len_NC(const NC *ncp, MPI_Offset sizeof_off_t);

extern int
ncx_put_NC(const NC *ncp, void **xpp, MPI_Offset offset, MPI_Offset extent);

extern int
nc_get_NC( NC *ncp);

/* End defined in v1hpg.c */

/* Begin defined in putget.c */

extern int
ncmpii_fill_NC_var(NC *ncp, const NC_var *varp, MPI_Offset recno);

extern int
ncmpii_inq_rec(int ncid, MPI_Offset *nrecvars, MPI_Offset *recvarids, MPI_Offset *recsizes);

extern int
ncmpii_get_rec(int ncid, MPI_Offset recnum, void **datap);

extern int
ncmpii_put_rec(int ncid, MPI_Offset recnum, void *const *datap);
#endif

/* End defined in putget.c */

/* Begin defined in header.c */
typedef struct bufferinfo {
    ncio       *nciop;
    MPI_Offset  offset;   /* current read/write offset in the file */
    int         version;  /* 1, 2, and 5 for CDF-1, 2, and 5 respectively */
    int         safe_mode;/* 0: disabled, 1: enabled */
    void       *base;     /* beginning of read/write buffer */
    void       *pos;      /* current position in buffer */
    MPI_Offset  size;     /* size of the buffer */
    MPI_Offset  index;    /* index of current position in buffer */
    MPI_Offset  put_size; /* amount of writes so far in bytes */
    MPI_Offset  get_size; /* amount of reads  so far in bytes */
} bufferinfo;

extern MPI_Offset
ncmpix_len_nctype(nc_type type);

#if 0
extern int
hdr_put_NC_attrarray(bufferinfo *pbp, const NC_attrarray *ncap);
#endif

extern MPI_Offset
ncmpii_hdr_len_NC(const NC *ncp);

extern int
ncmpii_hdr_get_NC(NC *ncp);

extern int
ncmpii_hdr_put_NC(NC *ncp, void *buf);

extern int
ncmpii_NC_computeshapes(NC *ncp);

extern int
ncmpii_hdr_check_NC(bufferinfo *getbuf, NC *ncp);
/* end defined in header.c */

/* begin defined in mpincio.c */
extern int
ncmpiio_create(MPI_Comm comm, const char *path, int ioflags, MPI_Info info,
               NC *ncp);

extern int
ncmpiio_open(MPI_Comm comm, const char *path, int ioflags, MPI_Info info,
             NC *ncp);
extern int
ncmpiio_sync(ncio *nciop);

extern int
ncmpiio_move(ncio *const nciop, MPI_Offset to, MPI_Offset from,
             MPI_Offset nbytes);

extern int
ncmpiio_move_fixed_vars(NC *ncp, NC *old);

extern int
ncmpiio_get_hint(NC *ncp, char *key, char *value, int *flag);

extern int
NC_computeshapes(NC *ncp);

/* end defined in mpincio.h */

/* begin defined in error.c */
const char* ncmpi_strerror(int err);

int ncmpii_handle_error(int mpi_errorcode, char *msg);

/* end defined in error.c */
/*
 * These functions are used to support
 * interface version 2 backward compatiblity.
 * N.B. these are tested in ../nc_test even though they are
 * not public. So, be careful to change the declarations in
 * ../nc_test/tests.h if you change these.
 */

int ncmpii_x_putn_schar (void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_uchar (void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_short (void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_ushort(void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_int   (void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_uint  (void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_float (void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_double(void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_int64 (void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_uint64(void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_schar (const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_uchar (const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_short (const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_ushort(const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_int   (const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_uint  (const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_float (const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_double(const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_int64 (const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_uint64(const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);

int NCedgeck(const NC *ncp, const NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *edges);

int NCstrideedgeck(const NC *ncp, const NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *edges, const MPI_Offset *stride);

int NCcoordck(NC *ncp, const NC_var *varp, const MPI_Offset *coord, const int rw_flag);


int ncmpii_need_convert(nc_type nctype,MPI_Datatype mpitype);

int ncmpii_need_swap(nc_type nctype,MPI_Datatype mpitype);

void ncmpii_swapn(void *dest_p, const void* src_p, MPI_Offset nelems, int esize);

void ncmpii_in_swapn(void *buf, MPI_Offset nelems, int esize);

int ncmpii_is_request_contiguous(NC *ncp, NC_var *varp,
                const MPI_Offset starts[], const MPI_Offset  counts[]);

int ncmpii_get_offset(NC *ncp, NC_var *varp, const MPI_Offset starts[],
                const MPI_Offset counts[], const MPI_Offset strides[],
                const int rw_flag, MPI_Offset *offset_ptr);

int ncmpii_check_mpifh(NC* ncp, MPI_File *mpifh, MPI_Comm comm,
                int collective);

int ncmpii_sync_numrecs(NC *ncp, MPI_Offset newnumrecs);

int ncmpii_vars_create_filetype(NC* ncp, NC_var* varp,
                const MPI_Offset start[], const MPI_Offset count[],
                const MPI_Offset stride[], int rw_flag, int *blocklen,
                MPI_Offset *offset, MPI_Datatype *filetype,
                int *is_filetype_contig);

int ncmpii_concatenate_datatypes(NC *ncp, int num, int *blocklens,
                MPI_Offset *displacements, MPI_Datatype *dtypes,
                MPI_Datatype *datatype);

extern int
ncmpii_getput_vars(NC *ncp, NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *count, const MPI_Offset *stride,
                void *buf, MPI_Offset bufcount, MPI_Datatype datatype,
                int rw_flag, int io_method, int isTempBuf);

extern int
ncmpii_getput_varm(NC *ncp, NC_var *varp, const MPI_Offset start[],
                const MPI_Offset count[], const MPI_Offset stride[],
                const MPI_Offset imap[], void *buf, MPI_Offset bufcount,
                MPI_Datatype datatype, int rw_flag, int io_method);
extern int
ncmpii_igetput_varm(NC *ncp, NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *stride, const MPI_Offset *imap,
                const MPI_Offset *count, void *buf, MPI_Offset bufcount,
                MPI_Datatype datatype, int *reqid, int rw_flag, int use_abuf,
                int isSameGroup);

extern int
ncmpii_wait(NC *ncp, int io_method, int num_reqs, int *req_ids,
                int *statuses);

extern int
ncmpii_cancel(NC *ncp, int num_req, int *req_ids, int *statuses);

extern int
ncmpii_getput_zero_req(NC *ncp, int rw_flag, int sync_numrecs);

extern void
ncmpii_inq_malloc_size(MPI_Offset *size);

extern int
ncmpii_inq_malloc_max_size(MPI_Offset *size);

extern int
ncmpii_inq_malloc_list(void);

extern int
ncmpii_inq_files_opened(int *num, int *ncids);

extern MPI_Datatype
ncmpii_nc2mpitype(nc_type type);

extern int
ncmpii_set_iget_callback(NC *ncp, int reqid, void *tmpBuf, int tmpBufSize,
                         void *userBuf, int userBufCount,
                         MPI_Datatype userBufType);

extern int
ncmpii_set_iput_callback(NC *ncp, int reqid, void *tmpPutBuf,
                         int need_swap_back_buf);

extern int
ncmpii_abuf_malloc(NC *ncp, MPI_Offset nbytes, void **buf);

#endif /* _NC_H_ */
