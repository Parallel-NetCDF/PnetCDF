/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *	Copyright 1996, University Corporation for Atmospheric Research
 *	See netcdf/COPYRIGHT file for copying and redistribution conditions.
 * 	
 * 	This file contains some routines derived from code
 *	which is copyrighted by Sun Microsystems, Inc.
 *	The "#ifdef vax" versions of
 *		 ncmpix_put_float_float()
 *		 ncmpix_get_float_float()
 *		 ncmpix_put_double_double()
 *		 ncmpix_get_double_double()
 *		 ncmpix_putn_float_float()
 *		 ncmpix_getn_float_float()
 *		 ncmpix_putn_double_double()
 *		 ncmpix_getn_double_double()
 * 	are derived from xdr_float() and xdr_double() routines
 *	in the freely available, copyrighted Sun RPCSRC 3.9
 *	distribution, xdr_float.c.
 * 	Our "value added" is that these are always memory to memory,
 *	they handle IEEE subnormals properly, and their "n" versions
 *	operate speedily on arrays.
 */
/* $Id$ */

/*
 * An external data representation interface.
 */

#include "nc.h"
#include "ncx.h"
#include <string.h>
#include <limits.h>
#include <sys/types.h>


/*
 * Primitive numeric conversion functions.
 */

/* We don't implement x_schar and x_uchar primitives. */

/*---- x_size_t -------------------------------------------------------------*/

#if SIZEOF_SIZE_T < X_SIZEOF_SIZE_T
#error "x_size_t implementation"
/* netcdf requires MPI_Offset which can hold a values from 0 to 2^31 -1 */
#endif

int
ncmpix_put_size_t1(void **xpp, const MPI_Offset *ulp)
{
	/* similar to put_ix_int() */
	uchar *cp = (uchar *) *xpp;
	assert(*ulp <= X_SIZE_MAX);

	*cp++ = (uchar)((*ulp) >> 24);
	*cp++ = (uchar)(((*ulp) & 0x00ff0000) >> 16);
	*cp++ = (uchar)(((*ulp) & 0x0000ff00) >>  8);
	*cp   = (uchar)((*ulp) & 0x000000ff);

	*xpp = (void *)((char *)(*xpp) + X_SIZEOF_SIZE_T);
	return NC_NOERR;
}

/*---< ncmpix_put_int32() >---------------------------------------------------*/
/* copy the contents of lp (a signed 32-bit integer) to xpp in Big Endian
 * form and advance *xpp 4 bytes
 */
int
ncmpix_put_int32(void      **xpp,
                 const int   ip)
{
#ifdef WORDS_BIGENDIAN
    int *ptr = (int*) (*xpp); /* typecast to 4-byte integer */
    *ptr = ip;
#else
    /* bitwise shifts below are to produce an integer in Big Endian */
    /* cannot call swap4b(), as lp is 8-byte */
    uchar *cp = (uchar *) *xpp;
    *cp++ = (uchar)((ip & 0xff000000) >> 24);
    *cp++ = (uchar)((ip & 0x00ff0000) >> 16);
    *cp++ = (uchar)((ip & 0x0000ff00) >>  8);
    *cp   = (uchar)( ip & 0x000000ff);
#endif
    /* advance *xpp 4 bytes */
    *xpp  = (void *)((char *)(*xpp) + 4);

    return NC_NOERR;
}

/*---< ncmpix_put_int64() >---------------------------------------------------*/
/* copy the contents of lp (a signed 64-bit integer) to xpp in Big Endian
 * form and advance *xpp 8 bytes
 */
int
ncmpix_put_int64(void             **xpp,
                 const MPI_Offset   ip)
{
#ifdef WORDS_BIGENDIAN
    MPI_Offset *ptr = (MPI_Offset*) (*xpp); /* typecast to 8-byte integer */
    *ptr = ip;
#else
    uchar *cp = (uchar *) *xpp;
    /* below is the same as calling swap8b(*xpp, &ip) */
    *cp++ = (uchar)((ip & 0xff00000000000000ULL) >> 56);
    *cp++ = (uchar)((ip & 0x00ff000000000000ULL) >> 48);
    *cp++ = (uchar)((ip & 0x0000ff0000000000ULL) >> 40);
    *cp++ = (uchar)((ip & 0x000000ff00000000ULL) >> 32);
    *cp++ = (uchar)((ip & 0x00000000ff000000ULL) >> 24);
    *cp++ = (uchar)((ip & 0x0000000000ff0000ULL) >> 16);
    *cp++ = (uchar)((ip & 0x000000000000ff00ULL) >>  8);
    *cp   = (uchar)( ip & 0x00000000000000ffULL);
#endif
    /* advance *xpp 8 bytes */
    *xpp  = (void *)((char *)(*xpp) + 8);

    return NC_NOERR;
}

/*---< ncmpix_put_size_t() >-------------------------------------------------*/
/* copy the contents of lp (a signed 64-bit integer) to xpp in Big Endian form
 * and advance *xpp to next size_t
 */
int
ncmpix_put_size_t(void             **xpp,
                  const MPI_Offset   lp,
                  int                sizeof_t) /* 4 or 8 */
{
    /* similar to put_ix_int() */
    uchar *cp = (uchar *) *xpp;

    /* No negative offsets stored in netcdf */
    /* if (lp < 0) return ERANGE; */

    assert(sizeof_t == 4 || sizeof_t == 8);

    /* bitwise shifts below are to produce an integer in Big Endian */
    if (sizeof_t == 4 ) {
#ifdef WORDS_BIGENDIAN
        int *ptr = (int*) (*xpp); /* typecast to 4-byte integer */
        *ptr = (int)lp;
#else
        /* cannot call swap4b(), as lp is 8-byte */
        *cp++ = (uchar)((lp & 0xff000000) >> 24);
        *cp++ = (uchar)((lp & 0x00ff0000) >> 16);
        *cp++ = (uchar)((lp & 0x0000ff00) >>  8);
        *cp   = (uchar)( lp & 0x000000ff);
#endif
    }
    else {
#ifdef WORDS_BIGENDIAN
        MPI_Offset *ptr = (MPI_Offset*) (*xpp); /* typecast to 8-byte integer */
        *ptr = lp;
#else
        /* below is the same as calling swap8b(*xpp, &lp) */
        *cp++ = (uchar)((lp & 0xff00000000000000ULL) >> 56);
        *cp++ = (uchar)((lp & 0x00ff000000000000ULL) >> 48);
        *cp++ = (uchar)((lp & 0x0000ff0000000000ULL) >> 40);
        *cp++ = (uchar)((lp & 0x000000ff00000000ULL) >> 32);
        *cp++ = (uchar)((lp & 0x00000000ff000000ULL) >> 24);
        *cp++ = (uchar)((lp & 0x0000000000ff0000ULL) >> 16);
        *cp++ = (uchar)((lp & 0x000000000000ff00ULL) >>  8);
        *cp   = (uchar)( lp & 0x00000000000000ffULL);
#endif
    }

    /* advance *xpp to next size_t */
    *xpp  = (void *)((char *)(*xpp) + sizeof_t);

    return NC_NOERR;
}

/*----< ncmpix_get_int32() >--------------------------------------------------*/
int
ncmpix_get_int32(void **xpp,
                 int   *ip)
{
    const uchar *cp = (const uchar *) *xpp;

    /* cannot call swap4b(), as lp is 8-byte */
    *ip  = (*cp++ << 24);
    *ip |= (*cp++ << 16);
    *ip |= (*cp++ <<  8);
    *ip |=  *cp; 

    /* advance *xpp 4 bytes */
    *xpp = (void *)((const char *)(*xpp) + 4);

    return NC_NOERR;
}

/*----< ncmpix_get_int64() >-------------------------------------------------*/
int
ncmpix_get_int64(void       **xpp,
                 MPI_Offset  *llp)
{
    const uchar *cp = (const uchar *) *xpp;

    /* below is the same as calling swap8b(llp, *xpp) */
    *llp  = ((MPI_Offset)(*cp++) << 56);
    *llp |= ((MPI_Offset)(*cp++) << 48);
    *llp |= ((MPI_Offset)(*cp++) << 40);
    *llp |= ((MPI_Offset)(*cp++) << 32);
    *llp |= ((MPI_Offset)(*cp++) << 24);
    *llp |= ((MPI_Offset)(*cp++) << 16);
    *llp |= ((MPI_Offset)(*cp++) <<  8);
    *llp |=  (MPI_Offset)*cp;

    /* advance *xpp 8 bytes */
    *xpp = (void *)((const char *)(*xpp) + 8);

    return NC_NOERR;
}

/*----< ncmpix_get_size_t() >-------------------------------------------------*/
int
ncmpix_get_size_t(const void **xpp,
                  MPI_Offset  *lp,
                  int          sizeof_t)  /* 4 or 8 */
{
    /* similar to get_ix_int() */
    const uchar *cp = (const uchar *) *xpp;

    assert(sizeof_t == 4 || sizeof_t == 8);

    if (sizeof_t == 4) {
        /* cannot call swap4b(), as lp is 8-byte */
        *lp  = (*cp++ << 24);
        *lp |= (*cp++ << 16);
        *lp |= (*cp++ <<  8);
        *lp |=  *cp; 
    }
    else {
        /* below is the same as calling swap8b(lp, *xpp) */
        *lp  = ((off_t)(*cp++) << 56);
        *lp |= ((off_t)(*cp++) << 48);
        *lp |= ((off_t)(*cp++) << 40);
        *lp |= ((off_t)(*cp++) << 32);
        *lp |= ((off_t)(*cp++) << 24);
        *lp |= ((off_t)(*cp++) << 16);
        *lp |= ((off_t)(*cp++) <<  8);
        *lp |=  (off_t)*cp;
    }

    /* advance *xpp to next size_t */
    *xpp = (const void *)((const char *)(*xpp) + sizeof_t);

    return NC_NOERR;
}

/* below ncmpix_put_off_t() and ncmpix_get_off_t() are the same as
 * ncmpix_put_size_t() and ncmpix_get_size_t()
 */
#if 0
/* x_off_t */

/* A previous version of this function would try to let systems with a 4 byte
 * off_t read and write the 8-byte offsets of a CDF-2 file.  For simplicity, we
 * now ensure that platforms with a 4-byte off_t will never open a CDF-2 file
 * no matter what the size. 
 */
int
ncmpix_put_off_t(void **xpp, const MPI_Offset *lp, MPI_Offset sizeof_off_t)
{
	/* similar to put_ix_int() */
	uchar *cp = (uchar *) *xpp;
		/* No negative offsets stored in netcdf */
	if (*lp < 0) {
		/* assume this is an overflow of a 32-bit int */
		return ERANGE;
	}
	assert(sizeof_off_t == 4 || sizeof_off_t == 8);
	if (sizeof_off_t == 4 ) {
		*cp++ = (uchar)(((*lp) & 0xff000000) >> 24);
		*cp++ = (uchar)(((*lp) & 0x00ff0000) >> 16);
		*cp++ = (uchar)(((*lp) & 0x0000ff00) >>  8);
		*cp   = (uchar)( (*lp) & 0x000000ff);
	} else {
		*cp++ = (uchar)(((*lp) & 0xff00000000000000ULL) >> 56);
		*cp++ = (uchar)(((*lp) & 0x00ff000000000000ULL) >> 48);
		*cp++ = (uchar)(((*lp) & 0x0000ff0000000000ULL) >> 40);
		*cp++ = (uchar)(((*lp) & 0x000000ff00000000ULL) >> 32);
		*cp++ = (uchar)(((*lp) & 0x00000000ff000000ULL) >> 24);
		*cp++ = (uchar)(((*lp) & 0x0000000000ff0000ULL) >> 16);
		*cp++ = (uchar)(((*lp) & 0x000000000000ff00ULL) >>  8);
		*cp   = (uchar)( (*lp) & 0x00000000000000ffULL);
	}
	*xpp = (void *)((char *)(*xpp) + sizeof_off_t);
	return NC_NOERR;
}

/* see comments for ncmpix_put_off_t */
int
ncmpix_get_off_t(const void **xpp, MPI_Offset *lp, MPI_Offset sizeof_off_t)
{
	/* similar to get_ix_int() */
	const uchar *cp = (const uchar *) *xpp;
	assert(sizeof_off_t == 4 || sizeof_off_t == 8);
	
       if (sizeof_off_t == 4) {
               *lp = *cp++ << 24;
               *lp |= (*cp++ << 16);
               *lp |= (*cp++ <<  8);
               *lp |= *cp; 
       } else {
               *lp =  ((off_t)(*cp++) << 56);
               *lp |= ((off_t)(*cp++) << 48);
               *lp |= ((off_t)(*cp++) << 40);
               *lp |= ((off_t)(*cp++) << 32);
               *lp |= ((off_t)(*cp++) << 24);
               *lp |= ((off_t)(*cp++) << 16);
               *lp |= ((off_t)(*cp++) <<  8);
               *lp |=  (off_t)*cp;
       }
       *xpp = (const void *)((const char *)(*xpp) + sizeof_off_t);


	return NC_NOERR;
}
#endif

/*
 * Aggregate numeric conversion functions.
 *
 *     getputn_schar.c
 *     getputn_uchar.c
 *     getputn_short.c
 *     getputn_ushort.c
 *     getputn_int.c
 *     getputn_uint.c
 *     getputn_float.c
 *     getputn_double.c
 *     getputn_int64.c
 *     getputn_uint64.c
*/

/*
 * Other aggregate conversion functions.
 */

/* text */

int
ncmpix_getn_text(const void **xpp, MPI_Offset nelems, char *tp)
{
	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);
	return NC_NOERR;

}

/*----< ncmpix_pad_getn_text() >----------------------------------------------*/
int
ncmpix_pad_getn_text(const void **xpp,
                     MPI_Offset   nelems,
                     char        *tp,
                     nc_type      dummy)                     
{
    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    memcpy(tp, *xpp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems + rndup);

    return NC_NOERR;
}

/*----< ncmpix_putn_text() >--------------------------------------------------*/
int
ncmpix_putn_text(void       **xpp,
                 MPI_Offset   nelems,
                 const char  *tp)
{
    memcpy(*xpp, tp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    return NC_NOERR;
}

/*----< ncmpix_pad_putn_text() >----------------------------------------------*/
int
ncmpix_pad_putn_text(void       **xpp,
                     MPI_Offset   nelems,
                     const char  *tp,
                     nc_type      dummy)                     
{
    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    memcpy(*xpp, tp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    if (rndup) { /* zero padding and advance the pointer */
        memcpy(*xpp, nada, rndup);
        *xpp = (void *)((char *)(*xpp) + rndup);
    }

    return NC_NOERR;
}


/* opaque */

int
ncmpix_getn_void(const void **xpp, MPI_Offset nelems, void *tp)
{
	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);
	return NC_NOERR;

}

int
ncmpix_pad_getn_void(const void **xpp, MPI_Offset nelems, void *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems + rndup);

	return NC_NOERR;

}

int
ncmpix_putn_void(void **xpp, MPI_Offset nelems, const void *tp)
{
	(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	return NC_NOERR;

}

int
ncmpix_pad_putn_void(void **xpp, MPI_Offset nelems, const void *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	if(rndup)
	{
		(void) memcpy(*xpp, nada, rndup);
		*xpp = (void *)((char *)(*xpp) + rndup);
	}
	
	return NC_NOERR;

}

