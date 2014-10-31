/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* "$Id$" */

#ifndef _NCX_H_
#define _NCX_H_

/*
 * An external data representation interface.
 *
 * This started out as a general replacement for ONC XDR,
 * specifically, the xdrmem family of functions.
 * 
 * We eventually realized that we could write more portable
 * code if we decoupled any association between the 'C' types
 * and the external types. (XDR has this association between the 'C'
 * types and the external representations, like xdr_int() takes
 * an int argument and goes to an external int representation.) 
 * So, now there is a matrix of functions.
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <assert.h>
#include <errno.h>
#include <sys/types.h> /* off_t */
#include <limits.h>
#include <float.h>
#include <stdint.h>

#include <mpi.h>

#include "rnd.h"
#include "nctypes.h"


/*
 * External sizes of the primitive elements.
 */
#define X_SIZEOF_CHAR		1
#define X_SIZEOF_SHORT		2
#define X_SIZEOF_INT		4	/* xdr_int */
#if 0
#define X_SIZEOF_LONG		8	/* xdr_long_long */
#endif
#define X_SIZEOF_FLOAT		4
#define X_SIZEOF_DOUBLE		8

/* wkliao: not sure I am doing this right ... */
#define X_SIZEOF_UBYTE		1
#define X_SIZEOF_USHORT		2
#define X_SIZEOF_UINT		4
#define X_SIZEOF_INT64		8
#define X_SIZEOF_UINT64		8
#define X_SIZEOF_STRING		8        /* char pointer ? */

/*
 * For now, netcdf is limited to 32 bit sizes. 
 * If compiled with support for "large files" then 
 * netcdf will use a 64 bit off_t and it can then write a file 
 * using 64 bit offsets.
 *  see also X_SIZE_MAX, X_OFF_MAX below
 */
#define X_SIZEOF_OFF_T		(sizeof(off_t))
#define X_SIZEOF_SIZE_T		X_SIZEOF_INT

#ifndef OLD_SIZES
/* copied from netCDF libsrc/ncx.h */
/*
 * limits of the external representation
 */
#define X_SCHAR_MIN     (-128)
#define X_SCHAR_MAX     127
#define X_UCHAR_MAX     255U
#define X_SHORT_MIN     (-32768)
#define X_SHRT_MIN      X_SHORT_MIN     /* alias compatible with limits.h */
#define X_SHORT_MAX     32767
#define X_SHRT_MAX      X_SHORT_MAX     /* alias compatible with limits.h */
#define X_USHORT_MAX    65535U
#define X_USHRT_MAX     X_USHORT_MAX    /* alias compatible with limits.h */
#define X_INT_MIN       (-2147483647-1)
#define X_INT_MAX       2147483647
#define X_UINT_MAX      4294967295U
#define X_LONGLONG_MIN  (-9223372036854775807LL-1LL)
#define X_LONGLONG_MAX  9223372036854775807LL
#define X_ULONGLONG_MAX 18446744073709551615ULL
#define X_FLOAT_MAX     3.402823466e+38f
#define X_FLOAT_MIN     (-X_FLOAT_MAX)
#define X_FLT_MAX       X_FLOAT_MAX     /* alias compatible with limits.h */
#if CRAYFLOAT
/* ldexp(1. - ldexp(.5 , -46), 1024) */
#define X_DOUBLE_MAX    1.79769313486230e+308
#else
/* scalb(1. - scalb(.5 , -52), 1024) */
#define X_DOUBLE_MAX    1.7976931348623157e+308 
#endif
#define X_DOUBLE_MIN    (-X_DOUBLE_MAX)
#define X_DBL_MAX       X_DOUBLE_MAX    /* alias compatible with limits.h */

#define X_SIZE_MAX      X_UINT_MAX
#define X_OFF_MAX       X_INT_MAX

/* copied from netCDF libsrc/ncx.c */
/* alias poorly named limits.h macros */
#define  SHORT_MAX  SHRT_MAX
#define  SHORT_MIN  SHRT_MIN
#define USHORT_MAX USHRT_MAX
#ifndef LLONG_MAX
#   define LLONG_MAX    9223372036854775807LL
#   define LLONG_MIN    (-LLONG_MAX - 1LL)
#   define ULLONG_MAX   18446744073709551615ULL
#endif
#ifndef LONG_LONG_MAX
#define LONG_LONG_MAX LLONG_MAX
#endif
#ifndef LONG_LONG_MIN
#define LONG_LONG_MIN LLONG_MIN
#endif
#ifndef ULONG_LONG_MAX
#define ULONG_LONG_MAX ULLONG_MAX
#endif
#include <float.h>
#ifndef FLT_MAX /* This POSIX macro missing on some systems */
# ifndef NO_IEEE_FLOAT
# define FLT_MAX 3.40282347e+38f
# else
# error "You will need to define FLT_MAX"
# endif
#endif
/* alias poorly named float.h macros */
#define FLOAT_MAX FLT_MAX
#define FLOAT_MIN (-FLT_MAX)
#define DOUBLE_MAX DBL_MAX
#define DOUBLE_MIN (-DBL_MAX)
#define FLOAT_MAX_EXP FLT_MAX_EXP
#define DOUBLE_MAX_EXP DBL_MAX_EXP
#include <assert.h>
#define UCHAR_MIN 0
#define Min(a,b) ((a) < (b) ? (a) : (b))
#define Max(a,b) ((a) > (b) ? (a) : (b))

/*
 *  * If the machine's float domain is "smaller" than the external one
 *   * use the machine domain
 *    */
#if defined(FLT_MAX_EXP) && FLT_MAX_EXP < 128 /* 128 is X_FLT_MAX_EXP */
#undef X_FLOAT_MAX
# define X_FLOAT_MAX FLT_MAX
#undef X_FLOAT_MIN
# define X_FLOAT_MIN (-X_FLOAT_MAX)
#endif

#if _SX /* NEC SUPER UX */
#define LOOPCNT 256    /* must be no longer than hardware vector length */
#if _INT64
#undef  INT_MAX /* workaround cpp bug */
#define INT_MAX  X_INT_MAX
#undef  INT_MIN /* workaround cpp bug */
#define INT_MIN  X_INT_MIN
#undef  LONG_MAX /* workaround cpp bug */
#define LONG_MAX  X_INT_MAX
#undef  LONG_MIN /* workaround cpp bug */
#define LONG_MIN  X_INT_MIN
#elif _LONG64
#undef  LONG_MAX /* workaround cpp bug */
#define LONG_MAX  4294967295L
#undef  LONG_MIN /* workaround cpp bug */
#define LONG_MIN -4294967295L
#endif
#if !_FLOAT0
#error "FLOAT1 and FLOAT2 not supported"
#endif
#endif /* _SX */

#else
/*
 * limits of the external representation
 * we rely on ANSI-C defined constants in limits.h. Do any modern environments
 * not have these?
 */
#define X_SCHAR_MIN	SCHAR_MIN
#define X_SCHAR_MAX	SCHAR_MAX
#define X_UCHAR_MAX	UCHAR_MAX
#define X_SHORT_MIN	SHRT_MIN
#define X_SHRT_MIN	X_SHORT_MIN	/* alias compatible with limits.h */
#define X_SHORT_MAX	SHRT_MAX
#define X_SHRT_MAX	X_SHORT_MAX	/* alias compatible with limits.h */
#define X_USHORT_MAX	USHRT_MAX
#define X_USHRT_MAX	X_USHORT_MAX	/* alias compatible with limits.h */
#define X_INT_MIN	INT_MIN
#define X_INT_MAX	INT_MAX
#define X_UINT_MAX	UINT_MAX
#if 0
#define X_LONG_MIN	(-2147483647-1)
#define X_LONG_MAX	2147483647
#define X_ULONG_MAX	4294967295U
#endif
#define X_FLOAT_MAX	FLT_MAX
#define X_FLOAT_MIN	-(FLT_MAX)
#define X_FLT_MAX	X_FLOAT_MAX	/* alias compatible with limits.h */
/* scalb(1. - scalb(.5 , -52), 1024) */
#define X_DOUBLE_MAX	DBL_MAX
#define X_DOUBLE_MIN	-(DBL_MAX)
#define X_DBL_MAX	X_DOUBLE_MAX	/* alias compatible with limits.h */

#define X_SIZE_MAX	X_UINT_MAX	
#define X_OFF_MAX	X_INT_MAX

#define X_INT64_T_MAX	9223372036854775807LL
#define X_UINT64_T_MAX	18446744073709551616U

#ifndef LLONG_MAX
#define LLONG_MAX  0x7fffffffffffffffLL
#endif
#ifndef LLONG_MIN
#define LLONG_MIN (-0x7fffffffffffffffLL-1)
#endif
#ifndef ULLONG_MAX
#define ULLONG_MAX  0xffffffffffffffffULL
#endif

#ifndef X_INT64_MAX
#define X_INT64_MAX    LLONG_MAX
#endif

#ifndef X_INT64_MIN
#define X_INT64_MIN    LLONG_MIN
#endif

#ifndef X_UINT64_MAX
#define X_UINT64_MAX  ULLONG_MAX
#endif

/* alias poorly named limits.h macros */
#define  SHORT_MAX  SHRT_MAX
#define  SHORT_MIN  SHRT_MIN
#define USHORT_MAX USHRT_MAX

#ifndef FLT_MAX /* This POSIX macro missing on some systems */
# ifndef NO_IEEE_FLOAT
# define FLT_MAX 3.40282347e+38f
# else
# error "You will need to define FLT_MAX"
# endif
#endif

/*
 * If the machine's float domain is "smaller" than the external one
 * use the machine domain
 */
#if defined(FLT_MAX_EXP) && FLT_MAX_EXP < 128 /* 128 is X_FLT_MAX_EXP */
#undef X_FLOAT_MAX
# define X_FLOAT_MAX FLT_MAX
#undef X_FLOAT_MIN
# define X_FLOAT_MIN (-X_FLOAT_MAX)
#endif

#if defined(_SX) && _SX != 0 /* NEC SUPER UX */
#if _INT64
#undef  INT_MAX /* workaround cpp bug */
#define INT_MAX  X_INT_MAX
#undef  INT_MIN /* workaround cpp bug */
#define INT_MIN  X_INT_MIN
#undef  LONG_MAX /* workaround cpp bug */
#define LONG_MAX  X_INT_MAX
#undef  LONG_MIN /* workaround cpp bug */
#define LONG_MIN  X_INT_MIN
#elif _LONG64
#undef  LONG_MAX /* workaround cpp bug */
#define LONG_MAX  4294967295L
#undef  LONG_MIN /* workaround cpp bug */
#define LONG_MIN -4294967295L
#endif
#endif /* _SX */

#endif

/* Begin ncx_len */

/*
 * ncx_len_xxx() interfaces are defined as macros below, 
 * These give the length of an array of nelems of the type.
 * N.B. The 'char' and 'short' interfaces give the X_ALIGNED length.
 */
#define X_ALIGN			4	/* a.k.a. BYTES_PER_XDR_UNIT */

static const char nada[X_ALIGN] = {0, 0, 0, 0};

#define ncmpix_len_char(nelems) \
	_RNDUP((nelems), X_ALIGN)

#define ncmpix_len_short(nelems) \
	(((nelems) + (nelems)%2)  * X_SIZEOF_SHORT)

#define ncmpix_len_int(nelems) \
	((nelems) * X_SIZEOF_INT)

#if 0
#define ncmpix_len_long(nelems) \
	((nelems) * X_SIZEOF_LONG)
#endif

#define ncmpix_len_float(nelems) \
	((nelems) * X_SIZEOF_FLOAT)

#define ncmpix_len_double(nelems) \
	((nelems) * X_SIZEOF_DOUBLE)

#define ncmpix_len_int64(nelems) \
	((nelems) * X_SIZEOF_INT64)

/* End ncx_len */

#ifdef __CHAR_UNSIGNED__
	/* 'char' is unsigned, declare ncbyte as 'signed char' */
typedef signed char schar;

#else
	/* 'char' is signed */
typedef signed char schar;

#endif	/* __CHAR_UNSIGNED__ */


#ifndef WORDS_BIGENDIAN
/* LITTLE_ENDIAN: DEC and intel */
/*
 * Routines to convert to BIGENDIAN.
 * Optimize the swapn?b() and swap?b() routines aggressivly.
 */

#define SWAP2(a) ( (((a) & 0xff) << 8) | \
                   (((a) >> 8) & 0xff) )

#define SWAP4(a) ( ((a) << 24) | \
                  (((a) <<  8) & 0x00ff0000) | \
                  (((a) >>  8) & 0x0000ff00) | \
                  (((a) >> 24) & 0x000000ff) )

void
ncmpii_swapn2b(void *dst, const void *src, MPI_Offset nn);
void
ncmpii_swapn4b(void *dst, const void *src, MPI_Offset nn);

# ifndef vax
void
ncmpii_swap4b(void *dst, const void *src);
void
ncmpii_swap8b(void *dst, const void *src);
void
ncmpii_swapn8b(void *dst, const void *src, MPI_Offset nn);
# endif /* !vax */

#endif /* LITTLE_ENDIAN */


/*
 * Primitive numeric conversion functions.
 * The `put' functions convert from native internal
 * type to the external type, while the `get' functions
 * convert from the external to the internal.
 *
 * These take the form
 *	int ncx_get_{external_type}_{internal_type}(
 *		const void *xp,
 *		internal_type *ip
 *	);
 *	int ncx_put_{external_type}_{internal_type}(
 *		void *xp,
 *		const internal_type *ip
 *	);
 * where
 *	`external_type' and `internal_type' chosen from
		schar
		uchar
		short
		ushort
		int
		uint
		long
		ulong
		float
		double
 *
 * Not all combinations make sense.
 * We may not implement all combinations that make sense.
 * The netcdf functions that use this ncx interface don't
 * use these primitive conversion functions. They use the
 * aggregate conversion functions declared below.
 *
 * Storage for a single element of external type is at the `void * xp'
 * argument.
 *
 * Storage for a single element of internal type is at `ip' argument.
 *
 * These functions return 0 (NC_NOERR) when no error occured,
 * or NC_ERANGE when the value being converted is too large.
 * When NC_ERANGE occurs, an undefined (implementation dependent)
 * conversion may have occured.
 *
 * Note that loss of precision may occur silently.
 *
 */

/*
 * Other primitive conversion functions
 * N.B. slightly different interface
 * Used by netcdf.
 */

/* ncx_get_int_size_t */
extern int
ncmpix_get_size_t(const void **xpp, MPI_Offset *ulp, int sizeof_t);
/* ncx_get_int_off_t */
extern int
ncmpix_get_off_t(const void **xpp, MPI_Offset *lp, int sizeof_off_t);

/* ncx_put_int_size_t */
extern int
ncmpix_put_size_t(void **xpp, const MPI_Offset ulp, int sizeof_off_t);
extern int
ncmpix_put_size_t1(void **xpp, const MPI_Offset *ulp);
/* ncx_put_int_off_t */
extern int
ncmpix_put_off_t(void **xpp, const MPI_Offset *lp, MPI_Offset sizeof_off_t);

extern int
ncmpix_put_int32(void **xpp, const int ip);
extern int              
ncmpix_put_int64(void **xpp, const MPI_Offset ip);
extern int 
ncmpix_get_int32(void **xpp, int *ip);
extern int
ncmpix_get_int64(void **xpp, MPI_Offset *llp);

/*
 * Aggregate numeric conversion functions.
 * Convert an array.  Replaces xdr_array(...).
 * These functions are used by netcdf. Unlike the xdr
 * interface, we optimize for aggregate conversions.
 * This functions should be implemented to take advantage
 * of multiple processor / parallel hardware where available.
 *
 * These take the form
 *	int ncx_getn_{external_type}_{internal_type}(
 *		const void *xpp,
 *		size_t nelems,
 *		internal_type *ip
 *	);
 *	int ncx_putn_{external_type}_{internal_type}(
 *		void **xpp,
 *		size_t nelems,
 *		const internal_type *ip
 *	);
 * Where the types are as in the primitive numeric conversion functions.
 *
 * The value of the pointer to pointer argument, *xpp, is
 * expected to reference storage for `nelems' of the external
 * type.  On return, it modified to reference just past the last
 * converted external element.
 *
 * The types whose external size is less than X_ALIGN also have `pad'
 * interfaces. These round (and zero fill on put) *xpp up to X_ALIGN
 * boundaries. (This is the usual xdr behavior.)
 *
 * The `ip' argument should point to an array of `nelems' of
 * internal_type.
 *
 * Range errors (NC_ERANGE) for a individual values in the array 
 * DO NOT terminate the array conversion. All elements are converted,
 * with some having undefined values.
 * If any range error occurs, the function returns NC_ERANGE.
 *
 */

/*---- schar ----------------------------------------------------------------*/
extern int
ncmpix_getn_schar_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_schar_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_schar_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_schar_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_getn_schar_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_schar_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_schar_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_schar_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_schar_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_schar_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_schar_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_pad_getn_schar_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_pad_getn_schar_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_pad_getn_schar_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_pad_getn_schar_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_pad_getn_schar_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_pad_getn_schar_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_pad_getn_schar_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_pad_getn_schar_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_pad_getn_schar_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_pad_getn_schar_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_pad_getn_schar_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_schar_schar (void **xpp, MPI_Offset nelems, const schar   *ip);
extern int
ncmpix_putn_schar_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_schar_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_schar_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_schar_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_schar_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_schar_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_schar_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_schar_double(void **xpp, MPI_Offset nelems, const double  *ip);
extern int
ncmpix_putn_schar_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_schar_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);
 
extern int
ncmpix_pad_putn_schar_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_pad_putn_schar_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_pad_putn_schar_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_pad_putn_schar_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_pad_putn_schar_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_pad_putn_schar_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_pad_putn_schar_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_pad_putn_schar_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_pad_putn_schar_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_pad_putn_schar_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_pad_putn_schar_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);

/*---- uchar ----------------------------------------------------------------*/
extern int
ncmpix_getn_uchar_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_uchar_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_uchar_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_uchar_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_getn_uchar_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_uchar_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_uchar_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_uchar_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_uchar_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_uchar_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_uchar_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_pad_getn_uchar_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_pad_getn_uchar_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_pad_getn_uchar_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_pad_getn_uchar_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_pad_getn_uchar_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_pad_getn_uchar_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_pad_getn_uchar_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_pad_getn_uchar_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_pad_getn_uchar_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_pad_getn_uchar_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_pad_getn_uchar_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_uchar_schar (void **xpp, MPI_Offset nelems, const schar   *ip);
extern int
ncmpix_putn_uchar_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_uchar_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_uchar_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_uchar_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_uchar_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_uchar_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_uchar_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_uchar_double(void **xpp, MPI_Offset nelems, const double  *ip);
extern int
ncmpix_putn_uchar_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_uchar_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);
 
extern int
ncmpix_pad_putn_uchar_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_pad_putn_uchar_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_pad_putn_uchar_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_pad_putn_uchar_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_pad_putn_uchar_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_pad_putn_uchar_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_pad_putn_uchar_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_pad_putn_uchar_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_pad_putn_uchar_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_pad_putn_uchar_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_pad_putn_uchar_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);

/*---- short ----------------------------------------------------------------*/
extern int
ncmpix_getn_short_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_short_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_short_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_short_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_getn_short_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_short_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_short_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_short_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_short_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_short_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_short_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_pad_getn_short_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_pad_getn_short_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_pad_getn_short_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_pad_getn_short_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_pad_getn_short_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_pad_getn_short_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_pad_getn_short_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_pad_getn_short_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_pad_getn_short_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_pad_getn_short_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_pad_getn_short_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_short_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_putn_short_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_short_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_short_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_short_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_short_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_short_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_short_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_short_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_putn_short_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_short_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);
 
extern int
ncmpix_pad_putn_short_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_pad_putn_short_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_pad_putn_short_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_pad_putn_short_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_pad_putn_short_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_pad_putn_short_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_pad_putn_short_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_pad_putn_short_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_pad_putn_short_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_pad_putn_short_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_pad_putn_short_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);

/*---- ushort ---------------------------------------------------------------*/
extern int
ncmpix_getn_ushort_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_ushort_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_ushort_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_ushort_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_getn_ushort_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_ushort_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_ushort_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_ushort_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_ushort_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_ushort_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_ushort_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_pad_getn_ushort_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_pad_getn_ushort_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_pad_getn_ushort_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_pad_getn_ushort_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_pad_getn_ushort_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_pad_getn_ushort_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_pad_getn_ushort_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_pad_getn_ushort_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_pad_getn_ushort_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_pad_getn_ushort_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_pad_getn_ushort_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_ushort_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_putn_ushort_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_ushort_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_ushort_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_ushort_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_ushort_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_ushort_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_ushort_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_ushort_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_putn_ushort_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_ushort_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);
 
extern int
ncmpix_pad_putn_ushort_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_pad_putn_ushort_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_pad_putn_ushort_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_pad_putn_ushort_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_pad_putn_ushort_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_pad_putn_ushort_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_pad_putn_ushort_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_pad_putn_ushort_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_pad_putn_ushort_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_pad_putn_ushort_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_pad_putn_ushort_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);

/*---- int ------------------------------------------------------------------*/
extern int
ncmpix_getn_int_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_int_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_int_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_int_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_getn_int_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_int_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_int_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_long_long (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_int_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_int_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_int_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_int_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_int_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_putn_int_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_int_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_int_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_int_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_int_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_int_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_long_long (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_int_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_int_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_putn_int_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_int_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);
 
/*---- uint -----------------------------------------------------------------*/
extern int
ncmpix_getn_uint_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_uint_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_uint_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_uint_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_getn_uint_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_uint_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_uint_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_long_long (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_uint_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_uint_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_uint_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_uint_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_uint_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_putn_uint_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_uint_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_uint_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_uint_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_uint_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_uint_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_long_long (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_uint_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_uint_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_putn_uint_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_uint_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);
 
/*---- float ----------------------------------------------------------------*/
extern int
ncmpix_getn_float_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_float_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_float_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_float_ushort(const void **xpp, MPI_Offset nelems, ushort *ip);
extern int
ncmpix_getn_float_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_float_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_float_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_float_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_float_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_float_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_float_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_float_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_putn_float_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_float_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_float_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_float_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_float_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_float_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_float_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_float_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_putn_float_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_float_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);
 
/*---- double ---------------------------------------------------------------*/
extern int
ncmpix_getn_double_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_double_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_double_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_double_ushort(const void **xpp, MPI_Offset nelems, ushort  *ip);
extern int
ncmpix_getn_double_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_double_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_double_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_double_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_double_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_double_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_double_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_double_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_putn_double_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_double_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_double_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_double_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_double_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_double_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_double_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_double_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_putn_double_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_double_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);

/*---- int64 ----------------------------------------------------------------*/
extern int
ncmpix_getn_int64_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_int64_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_int64_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_int64_ushort(const void **xpp, MPI_Offset nelems, ushort  *ip);
extern int
ncmpix_getn_int64_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_int64_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_int64_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_int64_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_int64_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_int64_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_int64_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_int64_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_putn_int64_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_int64_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_int64_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_int64_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_int64_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_int64_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_int64_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_int64_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_putn_int64_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_int64_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);

/*---- uint64 ---------------------------------------------------------------*/
extern int
ncmpix_getn_uint64_schar (const void **xpp, MPI_Offset nelems, schar  *ip);
extern int
ncmpix_getn_uint64_uchar (const void **xpp, MPI_Offset nelems, uchar  *ip);
extern int
ncmpix_getn_uint64_short (const void **xpp, MPI_Offset nelems, short  *ip);
extern int
ncmpix_getn_uint64_ushort(const void **xpp, MPI_Offset nelems, ushort  *ip);
extern int
ncmpix_getn_uint64_int   (const void **xpp, MPI_Offset nelems, int    *ip);
extern int
ncmpix_getn_uint64_uint  (const void **xpp, MPI_Offset nelems, uint   *ip);
extern int
ncmpix_getn_uint64_long  (const void **xpp, MPI_Offset nelems, long   *ip);
extern int
ncmpix_getn_uint64_float (const void **xpp, MPI_Offset nelems, float  *ip);
extern int
ncmpix_getn_uint64_double(const void **xpp, MPI_Offset nelems, double *ip);
extern int
ncmpix_getn_uint64_int64 (const void **xpp, MPI_Offset nelems, int64  *ip);
extern int
ncmpix_getn_uint64_uint64(const void **xpp, MPI_Offset nelems, uint64 *ip);

extern int
ncmpix_putn_uint64_schar (void **xpp, MPI_Offset nelems, const schar  *ip);
extern int
ncmpix_putn_uint64_uchar (void **xpp, MPI_Offset nelems, const uchar  *ip);
extern int
ncmpix_putn_uint64_short (void **xpp, MPI_Offset nelems, const short  *ip);
extern int
ncmpix_putn_uint64_ushort(void **xpp, MPI_Offset nelems, const ushort *ip);
extern int
ncmpix_putn_uint64_int   (void **xpp, MPI_Offset nelems, const int    *ip);
extern int
ncmpix_putn_uint64_uint  (void **xpp, MPI_Offset nelems, const uint   *ip);
extern int
ncmpix_putn_uint64_long  (void **xpp, MPI_Offset nelems, const long   *ip);
extern int
ncmpix_putn_uint64_float (void **xpp, MPI_Offset nelems, const float  *ip);
extern int
ncmpix_putn_uint64_double(void **xpp, MPI_Offset nelems, const double *ip);
extern int
ncmpix_putn_uint64_int64 (void **xpp, MPI_Offset nelems, const int64  *ip);
extern int
ncmpix_putn_uint64_uint64(void **xpp, MPI_Offset nelems, const uint64 *ip);

 

/*
 * Other aggregate conversion functions.
 */

/* read ASCII characters */
extern int
ncmpix_getn_text(const void **xpp, MPI_Offset nchars, char *cp);
extern int
ncmpix_pad_getn_text(const void **xpp, MPI_Offset nchars, char *cp, nc_type dummy);

/* write ASCII characters */
extern int
ncmpix_putn_text(void **xpp, MPI_Offset nchars, const char *cp);
extern int
ncmpix_pad_putn_text(void **xpp, MPI_Offset nchars, const char *cp, nc_type dummy);

/* for symmetry */
#define ncmpix_getn_char_char(xpp, nelems, fillp) ncmpix_getn_text(xpp, nelems, fillp)
#define ncmpix_putn_char_char(xpp, nelems, fillp) ncmpix_putn_text(xpp, nelems, fillp)

/* read opaque data */
extern int
ncmpix_getn_void(const void **xpp, MPI_Offset nchars, void *vp);
extern int
ncmpix_pad_getn_void(const void **xpp, MPI_Offset nchars, void *vp);

/* write opaque data */
extern int
ncmpix_putn_void(void **xpp, MPI_Offset nchars, const void *vp);
extern int
ncmpix_pad_putn_void(void **xpp, MPI_Offset nchars, const void *vp);

#endif /* _NCX_H_ */
