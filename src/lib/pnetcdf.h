/*
 * Copyright 1993-1996 University Corporation for Atmospheric Research/Unidata
 * 
 * Portions of this software were developed by the Unidata Program at the 
 * University Corporation for Atmospheric Research.
 * 
 * Access and use of this software shall impose the following obligations
 * and understandings on the user. The user is granted the right, without
 * any fee or cost, to use, copy, modify, alter, enhance and distribute
 * this software, and any derivative works thereof, and its supporting
 * documentation for any purpose whatsoever, provided that this entire
 * notice appears in all copies of the software, derivative works and
 * supporting documentation.  Further, UCAR requests that the user credit
 * UCAR/Unidata in any publications that result from the use of this
 * software or in any product that includes this software. The names UCAR
 * and/or Unidata, however, may not be used in any advertising or publicity
 * to endorse or promote any products or commercial entity unless specific
 * written permission is obtained from UCAR/Unidata. The user also
 * understands that UCAR/Unidata is not obligated to provide the user with
 * any support, consulting, training or assistance of any kind with regard
 * to the use, operation and performance of this software nor to provide
 * the user with any updates, revisions, new versions or "bug fixes."
 * 
 * THIS SOFTWARE IS PROVIDED BY UCAR/UNIDATA "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL UCAR/UNIDATA BE LIABLE FOR ANY SPECIAL,
 * INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
 * FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
 * WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
 */
/* "$Id$" */

#ifndef _NETCDF_
#define _NETCDF_

#include <stddef.h> /* size_t, ptrdiff_t */
#include <errno.h>  /* netcdf functions sometimes return system errors */

#if defined(__cplusplus)
extern "C" {
#endif

/*
 *  The netcdf external data types
 */
typedef enum {
	NC_NAT =	0,	/* NAT = 'Not A Type' (c.f. NaN) */
	NC_BYTE =	1,	/* signed 1 byte integer */
	NC_CHAR =	2,	/* ISO/ASCII character */
	NC_SHORT =	3,	/* signed 2 byte integer */
	NC_INT =	4,	/* signed 4 byte integer */
	NC_FLOAT =	5,	/* single precision floating point number */
	NC_DOUBLE =	6	/* double precision floating point number */
} nc_type;

typedef nc_type ncmpi_type;


/*
 * 	Default fill values, used unless _FillValue attribute is set.
 * These values are stuffed into newly allocated space as appropriate.
 * The hope is that one might use these to notice that a particular datum
 * has not been set.
 */
#define NC_FILL_BYTE	((signed char)-127)
#define NC_FILL_CHAR	((char)0)
#define NC_FILL_SHORT	((short)-32767)
#define NC_FILL_INT	(-2147483647L)
#define NC_FILL_FLOAT	(9.9692099683868690e+36f) /* near 15 * 2^119 */
#define NC_FILL_DOUBLE	(9.9692099683868690e+36)


/*
 * The above values are defaults.
 * If you wish a variable to use a different value than the above
 * defaults, create an attribute with the same type as the variable
 * and the following reserved name. The value you give the attribute
 * will be used as the fill value for that variable.
 */
#define _FillValue	"_FillValue"


/*
 * 'mode' flags for nccreate and ncopen
 */
#define NC_NOWRITE	0	/* default is read only */
#define NC_WRITE    	0x1	/* read & write */
#define NC_CLOBBER	0
#define NC_NOCLOBBER	0x4	/* Don't destroy existing file on create */
#define NC_FILL		0	/* argument to ncsetfill to clear NC_NOFILL */
#define NC_NOFILL	0x100	/* Don't fill data section an records */
#define NC_LOCK		0x0400	/* Use locking if available */
#define NC_SHARE	0x0800	/* Share updates, limit cacheing */

/*
 * Let nc__create() or nc__open() figure out
 * as suitable chunk size.
 */
#define NC_SIZEHINT_DEFAULT 0

/*
 * In nc__enddef(), align to the chunk size.
 */
#define NC_ALIGN_CHUNK ((size_t)(-1))

/*
 * 'size' argument to ncdimdef for an unlimited dimension
 */
#define NC_UNLIMITED 0L

/*
 * attribute id to put/get a global attribute
 */
#define NC_GLOBAL -1


/*
 * These maximums are enforced by the interface, to facilitate writing
 * applications and utilities.  However, nothing is statically allocated to
 * these sizes internally.
 */
#define NC_MAX_DIMS	100	 /* max dimensions per file */
#define NC_MAX_ATTRS	2000	 /* max global or per variable attributes */
#define NC_MAX_VARS	2000	 /* max variables per file */
#define NC_MAX_NAME	128	 /* max length of a name */
#define NC_MAX_VAR_DIMS	NC_MAX_DIMS /* max per variable dimensions */


/*
 * The netcdf version 3 functions all return integer error status.
 * These are the possible values, in addition to certain
 * values from the system errno.h.
 */

#define NC_ISSYSERR(err)	((err) > 0)

#define	NC_NOERR	0	/* No Error */

#define NC_ENOTINDEP	(-26)	/* Operation not allowed in collective data mode */
#define NC_EINDEP	(-27)	/* Operation not allowed in independent data mode */
#define NC_EFILE	(-28)	/* Unknown error in file operation */
#define NC_EREAD	(-29)	/* Unknown error in reading file */
#define NC_EWRITE	(-30)	/* Unknown error in writting to file */
#define NC_EMULTIDEFINE (-31)	/* NC definations on multiprocesses conflict */
#define NC_EOFILE	(-32)	/* file open/creation failed */
#define	NC_EBADID	(-33)	/* Not a netcdf id */
#define	NC_ENFILE	(-34)	/* Too many netcdfs open */
#define	NC_EEXIST	(-35)	/* netcdf file exists && NC_NOCLOBBER */
#define	NC_EINVAL	(-36)	/* Invalid Argument */
#define	NC_EPERM	(-37)	/* Write to read only */
#define	NC_ENOTINDEFINE	(-38)	/* Operation not allowed in data mode */
#define	NC_EINDEFINE	(-39)	/* Operation not allowed in define mode */
#define	NC_EINVALCOORDS	(-40)	/* Index exceeds dimension bound */
#define	NC_EMAXDIMS	(-41)	/* NC_MAX_DIMS exceeded */
#define	NC_ENAMEINUSE	(-42)	/* String match to name in use */
#define NC_ENOTATT	(-43)	/* Attribute not found */
#define	NC_EMAXATTS	(-44)	/* NC_MAX_ATTRS exceeded */
#define NC_EBADTYPE	(-45)	/* Not a netcdf data type */
#define NC_EBADDIM	(-46)	/* Invalid dimension id or name */
#define NC_EUNLIMPOS	(-47)	/* NC_UNLIMITED in the wrong index */
#define	NC_EMAXVARS	(-48)	/* NC_MAX_VARS exceeded */
#define NC_ENOTVAR	(-49)	/* Variable not found */
#define NC_EGLOBAL	(-50)	/* Action prohibited on NC_GLOBAL varid */
#define NC_ENOTNC	(-51)	/* Not a netcdf file */
#define NC_ESTS        	(-52)	/* In Fortran, string too short */
#define NC_EMAXNAME    	(-53)	/* NC_MAX_NAME exceeded */
#define NC_EUNLIMIT    	(-54)	/* NC_UNLIMITED size already in use */
#define NC_ENORECVARS  	(-55)	/* nc_rec op when there are no record vars */
#define NC_ECHAR	(-56)	/* Attempt to convert between text & numbers */
#define NC_EEDGE	(-57)	/* Edge+start exceeds dimension bound */
#define NC_ESTRIDE	(-58)	/* Illegal stride */
#define NC_EBADNAME	(-59)	/* Attribute or variable name
                                         contains illegal characters */
/* N.B. following must match value in ncx.h */
#define NC_ERANGE	(-60)	/* Math result not representable */
#define NC_ENOMEM	(-61)	/* Memory allocation (malloc) failure */

/*
 * The Interface
 */

/* Declaration modifiers for DLL support (MSC et al) */

#if defined(DLL_NETCDF) /* define when library is a DLL */
#  if defined(DLL_EXPORT) /* define when building the library */
#   define MSC_EXTRA __declspec(dllexport)
#  else
#   define MSC_EXTRA __declspec(dllimport)
#  endif
#else
#define MSC_EXTRA
#endif	/* defined(DLL_NETCDF) */

# define EXTERNL extern MSC_EXTRA

#if defined(__cplusplus)
}
#endif

#include <mpi.h>

/* Begin Prototypes */


const char *
ncmpi_strerror(int err);

/* Begin Dataset Functions */


int ncmpi_create(MPI_Comm comm, const char *path, int cmode, MPI_Info info, int *ncidp); 

int ncmpi_open(MPI_Comm comm, const char *path, int omode, MPI_Info info, int *ncidp);


int ncmpi_enddef(int ncid);


int ncmpi_redef(int ncid);


int ncmpi_sync(int ncid);


int ncmpi_abort(int ncid);


int ncmpi_begin_indep_data(int ncid);


int ncmpi_end_indep_data(int ncid);


int ncmpi_close(int ncid);

/* End Dataset Functions */

/* Begin Define Mode Functions */


int ncmpi_def_dim(int ncid, const char *name, size_t len, int *idp);


int ncmpi_def_var(int ncid, const char *name, nc_type xtype, 
              int ndims, const int *dimidsp, int *varidp);


int ncmpi_rename_dim(int ncid, int dimid, const char *name);


int ncmpi_rename_var(int ncid, int varid, const char *name);

/* End Define Mode Functions */

/* Begin Inquiry Functions */


int ncmpi_inq(int ncid, int *ndimsp, int *nvarsp,
          int *ngattsp, int *unlimdimidp); 


int ncmpi_inq_ndims(int ncid, int *ndimsp);


int ncmpi_inq_nvars(int ncid, int *nvarsp);


int ncmpi_inq_natts(int ncid, int *ngattsp);


int ncmpi_inq_unlimdim(int ncid, int *unlimdimidp);


int ncmpi_inq_dimid(int ncid, const char *name, int *idp);


int ncmpi_inq_dim(int ncid, int dimid, char *name, size_t *lenp);


int ncmpi_inq_dimname(int ncid, int dimid, char *name);


int ncmpi_inq_dimlen(int ncid, int dimid, size_t *lenp);


int ncmpi_inq_var(int ncid, int varid, char *name,
              nc_type *xtypep, int *ndimsp, int *dimidsp,
              int *nattsp);


int ncmpi_inq_varid(int ncid, const char *name, int *varidp);


int ncmpi_inq_varname(int ncid, int varid, char *name);


int ncmpi_inq_vartype(int ncid, int varid, nc_type *xtypep);


int ncmpi_inq_varndims(int ncid, int varid, int *ndimsp);


int ncmpi_inq_vardimid(int ncid, int varid, int *dimidsp);


int ncmpi_inq_varnatts(int ncid, int varid, int *nattsp);

/* End Inquiry Functions */

/* Begin _att */


int ncmpi_inq_att(int ncid, int varid, const char *name,
              nc_type *xtypep, size_t *lenp);


int ncmpi_inq_attid(int ncid, int varid, const char *name, int *idp);


int ncmpi_inq_atttype(int ncid, int varid, const char *name,
                  nc_type *xtypep);


int ncmpi_inq_attlen(int ncid, int varid, const char *name,
                 size_t *lenp);


int ncmpi_inq_attname(int ncid, int varid, int attnum, char *name);


int ncmpi_copy_att(int ncid_in, int varid_in, const char *name,
               int ncid_out, int varid_out);


int ncmpi_rename_att(int ncid, int varid, const char *name,
                 const char *newname);


int ncmpi_del_att(int ncid, int varid, const char *name);


int ncmpi_put_att_text(int ncid, int varid, const char *name, size_t len,
                   const char *op);


int ncmpi_get_att_text(int ncid, int varid, const char *name, char *ip);


int ncmpi_put_att_uchar(int ncid, int varid, const char *name,
                    nc_type xtype, size_t len, const unsigned char *op);


int ncmpi_get_att_uchar(int ncid, int varid, const char *name,
                    unsigned char *ip);


int ncmpi_put_att_schar(int ncid, int varid, const char *name,
                    nc_type xtype, size_t len, const signed char *op);


int ncmpi_get_att_schar(int ncid, int varid, const char *name,
                    signed char *ip);


int ncmpi_put_att_short(int ncid, int varid, const char *name,
                    nc_type xtype, size_t len, const short *op);


int ncmpi_get_att_short(int ncid, int varid, const char *name, short *ip);


int ncmpi_put_att_int(int ncid, int varid, const char *name,
                  nc_type xtype, size_t len, const int *op);


int ncmpi_get_att_int(int ncid, int varid, const char *name, int *ip);


int ncmpi_put_att_long(int ncid, int varid, const char *name,
                   nc_type xtype, size_t len, const long *op);


int ncmpi_get_att_long(int ncid, int varid, const char *name, long *ip);


int ncmpi_put_att_float(int ncid, int varid, const char *name,
                    nc_type xtype, size_t len, const float *op);


int ncmpi_get_att_float(int ncid, int varid, const char *name, float *ip);


int ncmpi_put_att_double(int ncid, int varid, const char *name,
                     nc_type xtype, size_t len, const double *op);


int ncmpi_get_att_double(int ncid, int varid, const char *name,
                     double *ip);

/* End _att */

/* Begin {put,get}_var1 */


int ncmpi_put_var1(int ncid, int varid,
               const size_t index[],
               const void *buf, int bufcount,
               MPI_Datatype datatype);


int ncmpi_get_var1(int ncid, int varid,
               const size_t index[],
               void *buf, int bufcount,
               MPI_Datatype datatype);


int ncmpi_put_var1_text(int ncid, int varid,
                    const size_t index[],
                    const char *op);


int ncmpi_put_var1_short(int ncid, int varid,
                     const size_t index[],
                     const short *op);


int ncmpi_put_var1_int(int ncid, int varid,
                   const size_t index[],
                   const int *op);


int ncmpi_put_var1_float(int ncid, int varid,
                     const size_t index[],
                     const float *op);


int ncmpi_put_var1_double(int ncid, int varid,
                      const size_t index[],
                      const double *op);


int ncmpi_get_var1_text(int ncid, int varid,
                    const size_t index[],
                    char *ip);


int ncmpi_get_var1_short(int ncid, int varid,
                     const size_t index[],
                     short *ip);


int ncmpi_get_var1_int(int ncid, int varid,
                   const size_t index[],
                   int *ip);


int ncmpi_get_var1_float(int ncid, int varid,
                     const size_t index[],
                     float *ip);


int ncmpi_get_var1_double(int ncid, int varid,
                      const size_t index[],
                      double *ip);

/* End {put,get}_var1 */

/* Begin {put,get}_var */  


int ncmpi_put_var(int ncid, int varid, const void *buf, int bufcount, MPI_Datatype datatype);


int ncmpi_get_var(int ncid, int varid, void *buf, int bufcount, MPI_Datatype datatype);


int ncmpi_get_var_all(int ncid, int varid, void *buf, int bufcount, MPI_Datatype datatype);


int ncmpi_put_var_text(int ncid, int varid, const char *op);


int ncmpi_put_var_short(int ncid, int varid, const short *op);


int ncmpi_put_var_int(int ncid, int varid, const int *op);


int ncmpi_put_var_float(int ncid, int varid, const float *op);


int ncmpi_put_var_double(int ncid, int varid, const double *op);


int ncmpi_get_var_text(int ncid, int varid, char *ip);


int ncmpi_get_var_short(int ncid, int varid, short *ip);


int ncmpi_get_var_int(int ncid, int varid, int *ip);


int ncmpi_get_var_float(int ncid, int varid, float *ip);


int ncmpi_get_var_double(int ncid, int varid, double *ip);


int ncmpi_get_var_text_all(int ncid, int varid, char *ip);


int ncmpi_get_var_short_all(int ncid, int varid, short *ip);


int ncmpi_get_var_int_all(int ncid, int varid, int *ip);


int ncmpi_get_var_float_all(int ncid, int varid, float *ip);


int ncmpi_get_var_double_all(int ncid, int varid, double *ip);

/* End {put,get}_var */

/* Begin {put,get}_vara */


int ncmpi_put_vara_all(int ncid, int varid,
                   const size_t start[], const size_t count[],
                   const void *buf, int bufcount,
                   MPI_Datatype datatype);


int ncmpi_get_vara_all(int ncid, int varid,
                   const size_t start[], const size_t count[],
                   void *buf, int bufcount,
                   MPI_Datatype datatype);


int ncmpi_put_vara(int ncid, int varid,
               const size_t start[], const size_t count[],
               const void *buf, int bufcount,
               MPI_Datatype datatype);


int ncmpi_get_vara(int ncid, int varid,
               const size_t start[], const size_t count[],
               void *buf, int bufcount,
               MPI_Datatype datatype);


int ncmpi_put_vara_text_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const char *op);


int ncmpi_put_vara_text(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const char *op);


int ncmpi_put_vara_short_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const short *op);


int ncmpi_put_vara_short(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const short *op);


int ncmpi_put_vara_int_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const int *op);


int ncmpi_put_vara_int(int ncid, int varid,
                const size_t start[], const size_t count[],
                const int *op);


int ncmpi_put_vara_float_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const float *op);


int ncmpi_put_vara_float(int ncid, int varid,
                const size_t start[], const size_t count[],
                const float *op);


int ncmpi_put_vara_double_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const double *op);


int ncmpi_put_vara_double(int ncid, int varid,
                const size_t start[], const size_t count[],
                const double *op);


int ncmpi_get_vara_text_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    char *ip);


int ncmpi_get_vara_text(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    char *ip);


int ncmpi_get_vara_short_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    short *ip);


int ncmpi_get_vara_short(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    short *ip);


int ncmpi_get_vara_int_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    int *ip);


int ncmpi_get_vara_int(int ncid, int varid,
                const size_t start[], const size_t count[],
                int *ip);


int ncmpi_get_vara_float_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    float *ip);


int ncmpi_get_vara_float(int ncid, int varid,
                const size_t start[], const size_t count[],
                float *ip);


int ncmpi_get_vara_double_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    double *ip);


int ncmpi_get_vara_double(int ncid, int varid,
                const size_t start[], const size_t count[],
                double *ip);

/* End {put,get}_vara */

/* Begin {put,get}_vars */


int ncmpi_put_vars_all(int ncid, int varid,
                   const size_t start[],
                   const size_t count[],
                   const size_t stride[],
                   const void *buf, int bufcount,
                   MPI_Datatype datatype);


int ncmpi_get_vars_all(int ncid, int varid,
                   const size_t start[],
                   const size_t count[],
                   const size_t stride[],
                   void *buf, int bufcount,
                   MPI_Datatype datatype);


int ncmpi_put_vars(int ncid, int varid,
               const size_t start[],
               const size_t count[],
               const size_t stride[],
               const void *buf, int bufcount,
               MPI_Datatype datatype);


int ncmpi_get_vars(int ncid, int varid,
               const size_t start[],
               const size_t count[],
               const size_t stride[],
               void *buf, int bufcount,
               MPI_Datatype datatype);


int ncmpi_put_vars_text_all(int ncid, int varid,
                        const size_t start[],
                        const size_t count[],
                        const size_t stride[],
                        const char *op);


int ncmpi_put_vars_text(int ncid, int varid,
                    const size_t start[],
                    const size_t count[],
                    const size_t stride[],
                    const char *op);


int ncmpi_put_vars_short_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         const short *op);


int ncmpi_put_vars_short(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     const short *op);


int ncmpi_put_vars_int_all(int ncid, int varid,
                       const size_t start[],
                       const size_t count[],
                       const size_t stride[],
                       const int *op);


int ncmpi_put_vars_int(int ncid, int varid,
                   const size_t start[],
                   const size_t count[],
                   const size_t stride[],
                   const int *op);


int ncmpi_put_vars_float_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         const float *op);


int ncmpi_put_vars_float(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     const float *op);


int ncmpi_put_vars_double_all(int ncid, int varid,
                          const size_t start[],
                          const size_t count[],
                          const size_t stride[],
                          const double *op);


int ncmpi_put_vars_double(int ncid, int varid,
                      const size_t start[],
                      const size_t count[],
                      const size_t stride[],
                      const double *op);


int ncmpi_get_vars_text_all(int ncid, int varid,
                        const size_t start[],
                        const size_t count[],
                        const size_t stride[],
                        char *ip);


int ncmpi_get_vars_text(int ncid, int varid,
                    const size_t start[],
                    const size_t count[],
                    const size_t stride[],
                    char *ip);


int ncmpi_get_vars_short_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         short *ip);


int ncmpi_get_vars_short(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     short *ip);


int ncmpi_get_vars_int_all(int ncid, int varid,
                       const size_t start[],
                       const size_t count[],
                       const size_t stride[],
                       int *ip);


int ncmpi_get_vars_int(int ncid, int varid,
                   const size_t start[],
                   const size_t count[],
                   const size_t stride[],
                   int *ip);


int ncmpi_get_vars_float_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         float *ip);


int ncmpi_get_vars_float(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     float *ip);


int ncmpi_get_vars_double_all(int ncid, int varid,
                          const size_t start[],
                          const size_t count[],
                          const size_t stride[],
                          double *ip);


int ncmpi_get_vars_double(int ncid, int varid,
                      const size_t start[],
                      const size_t count[],
                      const size_t stride[],
                      double *ip);

/* End Prototypes */

/* End {put,get}_vars */

/* These macros are defined in serial netcdf (3.5.0) for backwards
 * compatibility with older netcdf code.   We aren't concerned with backwards
 * compatibility, so if your code doesn't compile with parallel-netcdf, maybe
 * this is why: 
 *
 *
 *  OLD NAME                 NEW NAME
 *  ----------------------------------
 *  FILL_BYTE       NC_FILL_BYTE
 *  FILL_CHAR       NC_FILL_CHAR
 *  FILL_SHORT      NC_FILL_SHORT
 *  FILL_LONG       NC_FILL_INT
 *  FILL_FLOAT      NC_FILL_FLOAT
 *  FILL_DOUBLE     NC_FILL_DOUBLE
 *
 *  MAX_NC_DIMS     NC_MAX_DIMS
 *  MAX_NC_ATTRS    NC_MAX_ATTRS
 *  MAX_NC_VARS     NC_MAX_VARS
 *  MAX_NC_NAME     NC_MAX_NAME
 *  MAX_VAR_DIMS    NC_MAX_VAR_DIMS
 */
#endif /* _NETCDF_ */
