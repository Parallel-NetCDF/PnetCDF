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

/* for put request less than 4KB, copy it to a buffer and do byte swap there,
 * so if the user buffer is immutable (assuming smaller than 4KB), it will not
 * cause seg fault. Not a perfect solution, but should be sufficient for most
 * of the cases.
 */
#define NC_BYTE_SWAP_BUFFER_SIZE 4096

/* define MPI_OFFSET if not defined */
#ifndef HAVE_MPI_OFFSET_DATATYPE
    #ifdef HAVE_MPI_LONG_LONG_INT
        #define MPI_OFFSET MPI_LONG_LONG_INT
    #else
        #define MPI_OFFSET MPI_INT
    #endif
#endif

enum API_KIND {
    API_VARD, /* do not check start and count, no flexible APIs */
    API_VARN, /* do not check start and count */
    API_VAR,  /* do not check start and count */
    API_VAR1, /* check start */
    API_VARA, /* check start and count */
    API_VARS, /* check start and count */
    API_VARM  /* check start and count */
};

#define WRITE_REQ 0
#define READ_REQ  1

#define INDEP_IO 0
#define COLL_IO  1
#define NONBLOCKING_IO  -1

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


/* XXX: this seems really low.  do we end up spending a ton of time mallocing?
 * could we reduce that by increasing this to something 21st century? */
#ifndef NC_ARRAY_GROWBY
#define NC_ARRAY_GROWBY 4
#endif

/* ncmpi_create/ncmpi_open set up header to be 'chunksize' big and to grow
 * by 'chunksize' as new items added. This used to be 4k. 256k lets us read
 * in an entire climate header in one go */
#define NC_DEFAULT_CHUNKSIZE 262144

/* when variable's nctype is NC_CHAR, I/O buffer's MPI type must be MPI_CHAR
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

  Minimum happens when nothing is defined, i.e.
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
    NC_DIMENSION   = 10,  /* \x00 \x00 \x00 \x0A */
    NC_VARIABLE    = 11,  /* \x00 \x00 \x00 \x0B */
    NC_ATTRIBUTE   = 12   /* \x00 \x00 \x00 \x0C */
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
ncmpii_new_NC_string(size_t slen, const char *str);

extern int
ncmpii_set_NC_string(NC_string *ncstrp, const char *str);

/* End defined in string.c */

/*
 * NC dimension structure
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

#define NC_NAME_TABLE_CHUNK 16
#define HASH_TABLE_SIZE 256
/*
#define HASH_FUNC(x) ncmpii_jenkins_one_at_a_time_hash(x)
#define HASH_FUNC(x) ncmpii_additive_hash(x)
#define HASH_FUNC(x) ncmpii_rotating_hash(x)
#define HASH_FUNC(x) ncmpii_Pearson_hash(x)
*/
#define HASH_FUNC(x) ncmpii_Bernstein_hash(x)
/* Look like Bernstein's hashing function performs the best */

/* For the initial naive implementation of hashing:
 * #define HASH_FUNC(x) (unsigned char)x[0]
 * if used this simple hash function, HASH_TABLE_SIZE must be 256 which is the
 * number of possible keys can be stored in an unsigned char
 */

typedef struct NC_nametable {
    int  num;  /* number of elements added in the list array. Its value starts
                  with zero and incremented with new name is created. When its
                  value becomes a multiple of NC_NAME_TABLE_CHUNK, list will be
                  reallocated to a space of size (num + NC_NAME_TABLE_CHUNK) */
    int *list; /* dimension or variable IDs */
} NC_nametable;

/* the dimension ID returned from ncmpi_def_dim() is an integer pointer
 * which means the total number of defined dimension allowed in a file
 * is up to 2^31-1. Thus, the member ndefined below should be of type int.
 */
typedef struct NC_dimarray {
    int            nalloc;        /* number allocated >= ndefined */
    int            ndefined;      /* number of defined dimensions */
    int            unlimited_id;  /* -1 for not defined, otherwise >= 0 */
    NC_dim       **value;
    NC_nametable   nameT[HASH_TABLE_SIZE]; /* table for quick name lookup.
                    * indices 0, 1, ... HASH_TABLE_SIZE-1 are the hash keys.
                    * nameT[i].num is the number of hash collisions. The IDs of
                    * dimensions with names producing the same hash key i are
                    * stored in nameT[i].list[*]
                    */
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

/* End defined in dim.c */

/*
 * NC attribute
 *
 * Number of attributes is limited by 2^31-1 because the argument attnump in
 *  int nc_inq_attid(int ncid, int varid, const char *name, int *attnump);
 * is a signed 4-byte integer.
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

extern int
ncmpii_NC_findattr(const NC_attrarray *ncap, const char *uname);

/* attrarray */

extern void
ncmpii_free_NC_attrarray(NC_attrarray *ncap);

extern int
ncmpii_dup_NC_attrarray(NC_attrarray *ncap, const NC_attrarray *ref);

extern NC_attr *
ncmpii_elem_NC_attrarray(const NC_attrarray *ncap, MPI_Offset elem);

/* End defined in attr.c */

/*
 * NC variable: description and data
 */
typedef struct {
    int           xsz;    /* byte size of 1 array element */
    MPI_Offset   *shape;  /* dim->size of each dim */
    MPI_Offset   *dsizes; /* the right to left product of shape */
    NC_string    *name;   /* name of the variable */
    int           ndims;  /* number of dimensions */
    int          *dimids; /* array of dimension IDs */
    NC_attrarray  attrs;  /* attribute array */
    nc_type       type;   /* variable's data type */
    MPI_Offset    len;    /* this is the "vsize" defined in header format, the
                             total size in bytes of the array variable.
                             For record variable, this is the record size */
    MPI_Offset    begin;  /* starting file offset of this variable */
    char          no_fill;
#ifdef ENABLE_SUBFILING
    int           ndims_org;  /* ndims before subfiling */
    int          *dimids_org; /* dimids before subfiling */
    int           num_subfiles;
#endif
} NC_var;

/* note: we only allow less than 2^31-1 variables defined in a file */
typedef struct NC_vararray {
    int            nalloc;      /* number allocated >= ndefined */
    int            ndefined;    /* number of defined variables */
    int            num_rec_vars;/* number of defined record variables */
    NC_var       **value;
    NC_nametable   nameT[HASH_TABLE_SIZE]; /* table for quick name lookup.
                    * indices 0, 1, ... HASH_TABLE_SIZE-1 are the hash keys.
                    * nameT[i].num is the number of hash collisions. The IDs of
                    * variables with names producing the same hash key i are
                    * stored in nameT[i].list[*]
                    */
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
ncmpii_NC_var_shape64(NC_var *varp, const NC_dimarray *dims);

extern int
ncmpii_NC_check_vlen(NC_var *varp, MPI_Offset vlen_max);

extern int
ncmpii_NC_lookupvar(NC *ncp, int varid, NC_var **varp);

/* End defined in var.c */

#define IS_RECVAR(vp) \
        ((vp)->shape != NULL ? (*(vp)->shape == NC_UNLIMITED) : 0 )

/*
 *  The PnetCDF non-blocking I/O request type
 */
typedef struct NC_req {
    int            id;          /* even number for write, odd for read */
    void          *buf;         /* the original user buffer */
    void          *xbuf;        /* the buffer used to read/write, may point to
                                   the same address as buf */
    MPI_Offset     num_recs;    /* number of records requested (1 for
                                   fixed-size variable) */
    int            buftype_is_contig;
    int            need_swap_back_buf;
    int            abuf_index;  /* index in the abuf occupy_table
                                   -1 means not using attached buffer */

    void          *tmpBuf;      /* tmp buffer to be freed, used only by
                                   nonblocking varn when buftype is noncontig */
    void          *userBuf;     /* user buffer to be unpacked from tmpBuf. used
                                   only by by nonblocking varn when buftype is
                                   noncontig */
    NC_var        *varp;
    MPI_Offset    *start;        /* [varp->ndims] */
    MPI_Offset    *count;        /* [varp->ndims] */
    MPI_Offset    *stride;       /* [varp->ndims] */
    MPI_Offset     bnelems;      /* number of elements in user buffer */
    MPI_Offset     offset_start; /* starting of aggregate access region */
    MPI_Offset     offset_end;   /*   ending of aggregate access region */
    MPI_Offset     bufcount;     /* the number of buftype in this request */
    MPI_Datatype   buftype;      /* user defined derived data type */
    MPI_Datatype   ptype;        /* element data type in buftype */
    MPI_Datatype   imaptype;     /* derived data type constructed from imap */
    int           *status;       /* pointer to user's status */
} NC_req;

#define NC_ABUF_DEFAULT_TABLE_SIZE 128

typedef struct NC_buf_status {
    int        is_used;
    MPI_Aint   buf_addr;
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

/* chunk size for allocating read/write nonblocking request lists */
#define NC_REQUEST_CHUNK 1024

struct NC {
    /* linked list of currently opened netcdf files */
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
#define NC_CREAT  0x20000   /* in create phase, cleared by enddef */
#define NC_INDEF  0x80000   /* in define mode, cleared by enddef */
#define NC_NSYNC  0x100000  /* synchronise numrecs on change */
#define NC_HSYNC  0x200000  /* synchronise whole header on change */
#define NC_NDIRTY 0x400000  /* numrecs has changed */
#define NC_HDIRTY 0x800000  /* header info has changed */
/* NC_NOFILL is defined in netcdf.h, historical interface */
    int           format;   /* 1, 2, or 5 corresponding to CDF-1, 2, or 5 */
    int           flags;
    int           safe_mode;    /* 0 or 1, for parameter consistency check */
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
    int           numGetReqs;  /* number of pending nonblocking get requests */
    int           numPutReqs;  /* number of pending nonblocking put requests */
    NC_req       *get_list;    /* list of nonblocking read requests */
    NC_req       *put_list;    /* list of nonblocking write requests */
    NC_buf       *abuf;        /* attached buffer, used by bput APIs */
};

#define NC_readonly(ncp) \
        (!fIsSet((ncp)->nciop->ioflags, NC_WRITE))

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

#define NC_doFsync(ncp) \
        fIsSet((ncp)->nciop->ioflags, NC_SHARE)

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
        (NC_EMULTIDEFINE_FIRST >= (err) && (err) >= NC_EMULTIDEFINE_LAST)

#define IsPrimityMPIType(buftype) (buftype == MPI_FLOAT          || \
                                   buftype == MPI_DOUBLE         || \
                                   buftype == MPI_INT            || \
                                   buftype == MPI_CHAR           || \
                                   buftype == MPI_SIGNED_CHAR    || \
                                   buftype == MPI_UNSIGNED_CHAR  || \
                                   buftype == MPI_SHORT          || \
                                   buftype == MPI_UNSIGNED_SHORT || \
                                   buftype == MPI_UNSIGNED       || \
                                   buftype == MPI_LONG           || \
                                   buftype == MPI_LONG_LONG_INT  || \
                                   buftype == MPI_UNSIGNED_LONG_LONG)

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
ncmpii_write_header(NC *ncp);

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
ncmpii_close(NC *ncp);

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
} bufferinfo;

extern int
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

int ncmpii_handle_error(int mpi_errorcode, char *msg);

/* end defined in error.c */
/*
 * These functions are used to support
 * interface version 2 backward compatibility.
 * N.B. these are tested in ../nc_test even though they are
 * not public. So, be careful to change the declarations in
 * ../nc_test/tests.h if you change these.
 */

int ncmpii_x_putn_NC_CHAR  (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_putn_NC_BYTE  (int cdf_ver,
                           void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_UBYTE (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_SHORT (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_USHORT(void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_INT   (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_UINT  (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_FLOAT (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_DOUBLE(void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_INT64 (void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);
int ncmpii_x_putn_NC_UINT64(void *xbuf, const void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype, void *fillp);

int ncmpii_x_getn_NC_CHAR  (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_BYTE  (int cdf_ver,
                           const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_UBYTE (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_SHORT (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_USHORT(const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_INT   (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_UINT  (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_FLOAT (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_DOUBLE(const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_INT64 (const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);
int ncmpii_x_getn_NC_UINT64(const void *xbuf, void *buf, MPI_Offset nelems,
                           MPI_Datatype datatype);

int NC_start_count_stride_ck(const NC *ncp, const NC_var *varp,
                const MPI_Offset *start, const MPI_Offset *count,
                const MPI_Offset *stride, const int rw_flag);

int ncmpii_need_convert(int format, nc_type nctype, MPI_Datatype mpitype);

int ncmpii_need_swap(nc_type nctype,MPI_Datatype mpitype);

void ncmpii_swapn(void *dest_p, const void* src_p, MPI_Offset nelems, int esize);

void ncmpii_in_swapn(void *buf, MPI_Offset nelems, int esize);

int ncmpii_is_request_contiguous(NC *ncp, NC_var *varp,
                const MPI_Offset starts[], const MPI_Offset  counts[]);

int ncmpii_get_offset(NC *ncp, NC_var *varp, const MPI_Offset starts[],
                const MPI_Offset counts[], const MPI_Offset strides[],
                const int rw_flag, MPI_Offset *offset_ptr);

int ncmpii_check_mpifh(NC* ncp, int collective);

int ncmpii_write_numrecs(NC *ncp, MPI_Offset new_numrecs);

int ncmpii_sync_numrecs(NC *ncp, MPI_Offset newnumrecs);

int ncmpii_vars_create_filetype(NC* ncp, NC_var* varp,
                const MPI_Offset start[], const MPI_Offset count[],
                const MPI_Offset stride[], int rw_flag, int *blocklen,
                MPI_Offset *offset, MPI_Datatype *filetype,
                int *is_filetype_contig);

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
ncmpii_inq_malloc_size(MPI_Offset *size);

extern int
ncmpii_inq_malloc_max_size(MPI_Offset *size);

extern int
ncmpii_inq_malloc_list(void);

#ifdef PNC_MALLOC_TRACE
void ncmpii_init_malloc_tracing(void);
#endif

extern void *
NCI_Malloc_fn(size_t size, const int lineno, const char *func,
              const char *filename);

extern void *
NCI_Calloc_fn(size_t nelem, size_t elsize, const int lineno, const char *func,
              const char *filename);

extern void *
NCI_Realloc_fn(void *ptr, size_t size, const int lineno, const char *func,
               const char *filename);

extern void
NCI_Free_fn(void *ptr, const int lineno, const char *func,
            const char *filename);

extern int
ncmpii_inq_files_opened(int *num, int *ncids);

extern MPI_Datatype
ncmpii_nc2mpitype(nc_type type);

extern nc_type
ncmpii_mpi2nctype(MPI_Datatype itype);

extern int
ncmpii_end_indep_data(NC *ncp); 

extern int                
ncmpii_file_set_view(NC *ncp, MPI_File fh, MPI_Offset *offset, MPI_Datatype filetype);

extern int
ncmpii_create_imaptype(NC_var *varp, const MPI_Offset *count,
                       const MPI_Offset *imap, const MPI_Offset  bnelems,
                       const int el_size, MPI_Datatype ptype,
                       MPI_Datatype *imaptype);

extern int
ncmpii_calc_datatype_elems(NC *ncp, NC_var *varp, const MPI_Offset *start,
                           const MPI_Offset *count, const MPI_Offset *stride,
                           int rw_flag, MPI_Datatype buftype,
                           MPI_Datatype *ptype, MPI_Offset *bufcount,
                           MPI_Offset *bnelems, MPI_Offset *nbytes,
                           int *el_size, int *buftype_is_contig);

extern int
ncmpii_fill_vars(NC *ncp);

extern int
ncmpii_sanity_check(int ncid, int varid, const MPI_Offset *start,
                    const MPI_Offset *count, const MPI_Offset *stride,
                    MPI_Offset bufcount, MPI_Datatype buftype,
                    enum API_KIND api, int isFlexibleAPI, int mustInDataMode,
                    int rw_flag, int io_method, NC **ncp, NC_var **varp);

extern char*
ncmpii_err_code_name(int err);

extern int
ncmpii_jenkins_one_at_a_time_hash(const char *str_name);

extern int
ncmpii_additive_hash(const char *str_name);

extern int
ncmpii_rotating_hash(const char *str_name);

extern int
ncmpii_Bernstein_hash(const char *str_name);

extern int
ncmpii_Pearson_hash(const char *str_name);

extern int
ncmpii_update_name_lookup_table(NC_nametable *nameT, const int id,
                                const char *oldname, const char *newname);

extern int
ncmpii_inq_var_fill(NC_var *varp, void *fill_value);

extern int
ncmpii_inq_default_fill_value(int type, void *fill_value);

extern int
ncmpii_getput_zero_req(NC *ncp, int rw_flag);

#endif /* _NC_H_ */
