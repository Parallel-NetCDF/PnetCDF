/*
 *	Copyright 1996, University Corporation for Atmospheric Research
 *      See netcdf/COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id$ */
#ifndef _NC_H_
#define _NC_H_

/*
 *	netcdf library 'private' data structures, objects and interfaces
 */

#include "ncconfig.h"

#include <stddef.h>	/* size_t */
#include <sys/types.h>	/* off_t */
#include "pnetcdf.h"
#include "ncio.h"	/* ncio */
#include "fbits.h"


#ifndef NC_ARRAY_GROWBY
#define NC_ARRAY_GROWBY 4
#endif

/*
 * The extern size of an empty
 * netcdf version 1 file.
 * The initial value of ncp->xsz.
 */
#define MIN_NC_XSZ 32

typedef struct NC NC; /* forward reference */

/*
 *  The internal data types
 */
typedef enum {
	NC_UNSPECIFIED = 0,
/* future	NC_BITFIELD = 7, */
/*	NC_STRING =	8,	*/
	NC_DIMENSION =	10,
	NC_VARIABLE =	11,
	NC_ATTRIBUTE =	12
} NCtype;


/*
 * Counted string for names and such
 */
typedef struct {
	/* all xdr'd */
	size_t nchars;
	char *cp;
} NC_string;

extern NC *
ncmpii_new_NC(const size_t *chunkp);

extern NC *
ncmpii_dup_NC(const NC *ref);

/* Begin defined in string.c */
extern void
ncmpii_free_NC_string(NC_string *ncstrp);

extern int
ncmpii_NC_check_name(const char *name);

extern NC_string *
ncmpii_new_NC_string(size_t slen, const char *str);

extern int
ncmpii_set_NC_string(NC_string *ncstrp, const char *str);

/* End defined in string.c */

/*
 * NC dimension stucture
 */
typedef struct {
	/* all xdr'd */
	NC_string *name;
	size_t size;
} NC_dim;

typedef struct NC_dimarray {
	size_t nalloc;		/* number allocated >= nelems */
	/* below gets xdr'd */
	/* NCtype type = NC_DIMENSION */
	size_t nelems;		/* length of the array */
	NC_dim **value;
} NC_dimarray;

/* Begin defined in dim.c */

extern void
ncmpii_free_NC_dim(NC_dim *dimp);

extern NC_dim *
ncmpii_new_x_NC_dim(NC_string *name);

extern int
ncmpii_find_NC_Udim(const NC_dimarray *ncap, NC_dim **dimpp);

/* dimarray */

extern void
ncmpii_free_NC_dimarrayV0(NC_dimarray *ncap);

extern void
ncmpii_free_NC_dimarrayV(NC_dimarray *ncap);

extern int
ncmpii_dup_NC_dimarrayV(NC_dimarray *ncap, const NC_dimarray *ref);

extern NC_dim *
ncmpii_elem_NC_dimarray(const NC_dimarray *ncap, size_t elem);

extern int
ncmpi_def_dim(int ncid, const char *name, MPI_Offset size, int *dimidp);

extern int
ncmpi_rename_dim( int ncid, int dimid, const char *newname);

extern int
ncmpi_inq_dimid(int ncid, const char *name, int *dimid_ptr);

extern int
ncmpi_inq_dim(int ncid, int dimid, char *name, int *sizep);

extern int 
ncmpi_inq_dimname(int ncid, int dimid, char *name);

extern int 
ncmpi_inq_dimlen(int ncid, int dimid, int *lenp);
/* End defined in dim.c */

/*
 * NC attribute
 */
typedef struct {
	size_t xsz;		/* amount of space at xvalue */
	/* below gets xdr'd */
	NC_string *name;
	nc_type type;		/* the discriminant */
	size_t nelems;		/* length of the array */
	void *xvalue;		/* the actual data, in external representation */
} NC_attr;

typedef struct NC_attrarray {
	size_t nalloc;		/* number allocated >= nelems */
	/* below gets xdr'd */
	/* NCtype type = NC_ATTRIBUTE */
	size_t nelems;		/* length of the array */
	NC_attr **value;
} NC_attrarray;

/* Begin defined in attr.c */

extern void
ncmpii_free_NC_attr(NC_attr *attrp);

extern NC_attr *
ncmpii_new_x_NC_attr(
	NC_string *strp,
	nc_type type,
	size_t nelems);

extern NC_attr **
ncmpii_NC_findattr(const NC_attrarray *ncap, const char *name);

/* attrarray */

extern void
ncmpii_free_NC_attrarrayV0(NC_attrarray *ncap);

extern void
ncmpii_free_NC_attrarrayV(NC_attrarray *ncap);

extern int
ncmpii_dup_NC_attrarrayV(NC_attrarray *ncap, const NC_attrarray *ref);

extern NC_attr *
ncmpii_elem_NC_attrarray(const NC_attrarray *ncap, size_t elem);

extern int
ncmpi_put_att_text(int ncid, int varid, const char *name,
	int nelems, const char *value);

extern int
ncmpi_get_att_text(int ncid, int varid, const char *name, char *str);

extern int
ncmpi_put_att_schar(int ncid, int varid, const char *name,
	nc_type type, int nelems, const signed char *value);

extern int
ncmpi_get_att_schar(int ncid, int varid, const char *name, signed char *tp);

extern int
ncmpi_put_att_uchar(int ncid, int varid, const char *name,
	nc_type type, int nelems, const unsigned char *value);

extern int
ncmpi_get_att_uchar(int ncid, int varid, const char *name, unsigned char *tp);

extern int
ncmpi_put_att_short(int ncid, int varid, const char *name,
	nc_type type, int nelems, const short *value);

extern int
ncmpi_get_att_short(int ncid, int varid, const char *name, short *tp);

extern int
ncmpi_put_att_int(int ncid, int varid, const char *name,
	nc_type type, int nelems, const int *value);

extern int
ncmpi_get_att_int(int ncid, int varid, const char *name, int *tp);

extern int
ncmpi_put_att_long(int ncid, int varid, const char *name,
	nc_type type, int nelems, const long *value);

extern int
ncmpi_get_att_long(int ncid, int varid, const char *name, long *tp);

extern int
ncmpi_put_att_float(int ncid, int varid, const char *name,
	nc_type type, int nelems, const float *value);
extern int
ncmpi_get_att_float(int ncid, int varid, const char *name, float *tp);
extern int
ncmpi_put_att_double(int ncid, int varid, const char *name,
	nc_type type, int nelems, const double *value);
extern int
ncmpi_get_att_double(int ncid, int varid, const char *name, double *tp);

extern int 
ncmpi_inq_attid(int ncid, int varid, const char *name, int *attnump);

extern int 
ncmpi_inq_atttype(int ncid, int varid, const char *name, nc_type *datatypep);

extern int 
ncmpi_inq_attlen(int ncid, int varid, const char *name, int *lenp);

extern int
ncmpi_inq_att(int ncid, int varid, const char *name, 
	nc_type *datatypep, int *lenp);

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
	size_t xsz;		/* xszof 1 element */
	size_t *shape; /* compiled info: dim->size of each dim */
	size_t *dsizes; /* compiled info: the right to left product of shape */
	/* below gets xdr'd */
	NC_string *name;
	/* next two: formerly NC_iarray *assoc */ /* user definition */
	size_t ndims;	/* assoc->count */
	int *dimids;	/* assoc->value */
	NC_attrarray attrs;
	nc_type type;		/* the discriminant */
	size_t len;		/* the total length originally allocated */
	MPI_Offset begin;
} NC_var;

typedef struct NC_vararray {
	size_t nalloc;		/* number allocated >= nelems */
	/* below gets xdr'd */
	/* NCtype type = NC_VARIABLE */
	size_t nelems;		/* length of the array */
	NC_var **value;
} NC_vararray;

/* Begin defined in var.c */

extern void
ncmpii_free_NC_var(NC_var *varp);

extern NC_var *
ncmpii_new_x_NC_var(
	NC_string *strp,
	size_t ndims);

/* vararray */

extern void
ncmpii_free_NC_vararrayV0(NC_vararray *ncap);

extern void
ncmpii_free_NC_vararrayV(NC_vararray *ncap);

extern int
ncmpii_dup_NC_vararrayV(NC_vararray *ncap, const NC_vararray *ref);

extern int
ncmpii_NC_var_shape(NC_var *varp, const NC_dimarray *dims);

extern int
ncmpii_NC_findvar(const NC_vararray *ncap, const char *name, NC_var **varpp);

extern int
ncmpii_NC_check_vlen(NC_var *varp, size_t vlen_max);

extern NC_var *
ncmpii_NC_lookupvar(NC *ncp, int varid);

extern int
ncmpi_def_var( int ncid, const char *name, nc_type type,
	 int ndims, const int *dimids, int *varidp);

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

struct NC {
	/* links to make list of open netcdf's */
	struct NC *next;
	struct NC *prev;
	/* contains the previous NC during redef. */
	struct NC *old;
	/* flags */
#define NC_INDEP 1	/* in independent data mode, cleared by endindep */
#define NC_CREAT 2	/* in create phase, cleared by ncenddef */
#define NC_INDEF 8	/* in define mode, cleared by ncenddef */
#define NC_NSYNC 0x10	/* synchronise numrecs on change */
#define NC_HSYNC 0x20	/* synchronise whole header on change */
#define NC_NDIRTY 0x40	/* numrecs has changed */
#define NC_HDIRTY 0x80  /* header info has changed */
/*	NC_NOFILL in netcdf.h, historical interface */
	int flags;
	ncio *nciop;
	size_t chunk;	/* largest extent this layer will request from ncio->get() */
	MPI_Offset xsz;	/* external size of this header, <= var[0].begin */
	MPI_Offset begin_var; /* position of the first (non-record) var */
	MPI_Offset begin_rec; /* position of the first 'record' */
	/* don't constrain maximu sinze of record unnecessarily */
	MPI_Offset recsize;	/* length of 'record' */	
	/* below gets xdr'd */
	size_t numrecs; /* number of 'records' allocated */
	NC_dimarray dims;
	NC_attrarray attrs;
	NC_vararray vars;
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
/* Begin defined in nc.c */

extern int
ncmpii_NC_check_id(int ncid, NC **ncpp);

extern int
ncmpii_cktype(nc_type datatype);

extern size_t
ncmpix_howmany(nc_type type, size_t xbufsize);

extern int
ncmpii_read_numrecs(NC *ncp);

extern int
ncmpii_write_numrecs(NC *ncp);

extern int
ncmpii_NC_sync(NC *ncp);

extern void
ncmpii_free_NC(NC *ncp);

extern void
ncmpii_add_to_NCList(NC *ncp);

extern void
ncmpii_del_from_NCList(NC *ncp);

extern int
ncmpii_read_NC(NC *ncp);

extern int 
ncmpii_NC_enddef(NC *ncp);

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

/* End defined in nc.c */
/* Begin defined in v1hpg.c */

extern size_t
ncx_len_NC(const NC *ncp, size_t sizeof_off_t);

extern int
ncx_put_NC(const NC *ncp, void **xpp, off_t offset, size_t extent);

extern int
nc_get_NC( NC *ncp);

/* End defined in v1hpg.c */

#if 0
/* Begin defined in putget.c */

extern int
ncmpii_fill_NC_var(NC *ncp, const NC_var *varp, size_t recno);

extern int
ncmpii_inq_rec(int ncid, size_t *nrecvars, int *recvarids, size_t *recsizes);

extern int
ncmpii_get_rec(int ncid, size_t recnum, void **datap);

extern int
ncmpii_put_rec(int ncid, size_t recnum, void *const *datap);
#endif

/* End defined in putget.c */

/* Begin defined in header.c */
typedef struct bufferinfo {
  ncio *nciop;		
  MPI_Offset offset;	/* current read/write offset in the file */
  int version;		/* either 1 for normal netcdf or 
			   2 for 8-byte offset version */
  void *base;     	/* beginning of read/write buffer */
  void *pos;      	/* current position in buffer */
  size_t size;		/* size of the buffer */
  size_t index;		/* index of current position in buffer */
} bufferinfo;  

extern size_t 
ncmpix_len_nctype(nc_type type);

#if 0
extern int
hdr_put_NC_attrarray(bufferinfo *pbp, const NC_attrarray *ncap);
#endif

extern size_t
ncmpii_hdr_len_NC(const NC *ncp, size_t sizeof_off_t);

extern int
ncmpii_hdr_get_NC(NC *ncp);

extern int 
ncmpii_hdr_put_NC(NC *ncp, void *buf);

extern int
ncmpii_NC_computeshapes(NC *ncp);

/* end defined in header.c */

/* begin defined in mpincio.c */
extern int
ncmpiio_create(MPI_Comm comm, const char *path, int ioflags, MPI_Info info,
            ncio **nciopp);

extern int
ncmpiio_open(MPI_Comm comm, const char *path, int ioflags, MPI_Info info,
          ncio **nciopp);
extern int
ncmpiio_sync(ncio *nciop);

extern int
ncmpiio_move(ncio *const nciop, off_t to, off_t from, size_t nbytes);

extern int
NC_computeshapes(NC *ncp);

/* end defined in mpincio.h */

/* begin defined in error.c */
const char * nc_strerror(int err);

void ncmpii_handle_error(int rank, int mpi_status, char *msg);
/* end defined in error.c */
/*
 * These functions are used to support
 * interface version 2 backward compatiblity.
 * N.B. these are tested in ../nc_test even though they are
 * not public. So, be careful to change the declarations in
 * ../nc_test/tests.h if you change these.
 */

extern int
ncmpii_put_att(int ncid, int varid, const char *name, nc_type datatype,
	size_t len, const void *value);

extern int
ncmpii_get_att(int ncid, int varid, const char *name, void *value);

extern int
ncmpii_put_var1(int ncid, int varid, const MPI_Offset *index, const void *value);

extern int
ncmpii_get_var1(int ncid, int varid, const MPI_Offset *index, void *value);

extern int
ncmpii_put_vara(int ncid, int varid,
	 const MPI_Offset *start, const MPI_Offset *count, const void *value);

extern int
ncmpii_get_vara(int ncid, int varid,
	 const MPI_Offset *start, const MPI_Offset *count, void *value);

extern int
ncmpii_put_vars(int ncid, int varid,
	 const MPI_Offset *start, const MPI_Offset *count, const ptrdiff_t *stride,
	 const void * value);

extern int
ncmpii_get_vars(int ncid, int varid,
	 const MPI_Offset *start, const MPI_Offset *count, const ptrdiff_t *stride,
	 void * value);

extern int
ncmpii_put_varm(int ncid, int varid,
	 const MPI_Offset *start, const MPI_Offset *count, const ptrdiff_t *stride,
	 const ptrdiff_t * map, const void *value);

extern int
ncmpii_get_varm(int ncid, int varid,
	 const MPI_Offset *start, const MPI_Offset *count, const ptrdiff_t *stride,
	 const ptrdiff_t * map, void *value);


#endif /* _NC_H_ */
