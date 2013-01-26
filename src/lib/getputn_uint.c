/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include "ncx.h"

/* ftype is the variable's nc_type defined in file, eg. int64
 * btype is the I/O buffer's C data type, eg. long long
 * buftype is I/O bufer's MPI data type, eg. MPI_UNSIGNED_LONG_LONG
 * apitype is data type appeared in the API names, eg. ncmpi_get_vara_longlong
 */

/*---- x_uint ---------------------------------------------------------------*/

#if USHORT_MAX == X_UINT_MAX
    typedef unsigned short ix_uint;
    #define SIZEOF_IX_UINT SIZEOF_USHORT
    #define IX_UINT_MAX USHORT_MAX
#elif UINT_MAX  >= X_UINT_MAX
    typedef unsigned int ix_uint;
    #define SIZEOF_IX_UINT SIZEOF_UINT
    #define IX_UINT_MAX UINT_MAX
#elif ULONG_MAX  >= X_UINT_MAX
    typedef unsigned long ix_uint;
    #define SIZEOF_IX_UINT SIZEOF_ULONG
    #define IX_UINT_MAX ULONG_MAX
#else
    #error "ix_uint implementation"
#endif


static void
get_ix_uint(const void *xp, ix_uint *ip)
{
    const uchar *cp = (const uchar *) xp;

    *ip  =  *cp++ << 24;
    *ip |= (*cp++ << 16);
    *ip |= (*cp++ << 8);
    *ip |= *cp; 
}

static void
put_ix_uint(void *xp, const ix_uint *ip)
{
    uchar *cp = (uchar *) xp;

    *cp++ = (*ip) >> 24;
    *cp++ = ((*ip) & 0x00ff0000) >> 16;
    *cp++ = ((*ip) & 0x0000ff00) >>  8;
    *cp   = ((*ip) & 0x000000ff);
}

/*----< ncmpix_get_uint_short() >--------------------------------------------*/
static int
ncmpix_get_uint_short(const void *xp, short *ip)
{
    ix_uint xx=0;
    get_ix_uint(xp, &xx);
    *ip = xx;
#if IX_UINT_MAX > SHORT_MAX
    if (xx > SHORT_MAX)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}

/*----< ncmpix_get_uint_ushort() >-------------------------------------------*/
static int
ncmpix_get_uint_ushort(const void *xp, ushort *ip)
{
    ix_uint xx=0;
    get_ix_uint(xp, &xx);
    *ip = xx;
#if IX_UINT_MAX > USHORT_MAX
    if (xx > USHORT_MAX)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}

/*----< ncmpix_get_uint_int() >----------------------------------------------*/
static int
ncmpix_get_uint_int(const void *xp, int *ip)
{
    ix_uint xx;
    get_ix_uint(xp, &xx);
    *ip = xx;
#if IX_UINT_MAX > INT_MAX
    if (xx > INT_MAX)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}

#if SIZEOF_IX_UINT != SIZEOF_UINT
/*----< ncmpix_get_uint_uint() >---------------------------------------------*/
static int
ncmpix_get_uint_uint(const void *xp, uint *ip)
{
#if SIZEOF_IX_UINT == SIZEOF_UINT && IX_UINT_MAX == UINT_MAX
    get_ix_uint(xp, (ix_uint *)ip);
    return NC_NOERR;
#else
    ix_uint xx;
    get_ix_uint(xp, &xx);
    *ip = xx;
#  if IX_UINT_MAX > UINT_MAX
    if (xx > UINT_MAX)
        return NC_ERANGE;
#  endif
    return NC_NOERR;
#endif
}
#endif

/*----< ncmpix_get_uint_long() >---------------------------------------------*/
static int
ncmpix_get_uint_long(const void *xp, long *ip)
{
    ix_uint xx;
    get_ix_uint(xp, &xx);
    *ip = xx;
#if IX_UINT_MAX > LONG_MAX
    if (xx > LONG_MAX)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}


#define GET_UINT(btype, range_check)                                          \
static int                                                                    \
ncmpix_get_uint_##btype(const void *xp, btype *ip)                            \
{                                                                             \
    ix_uint xx;                                                               \
    get_ix_uint(xp, &xx);                                                     \
    *ip = xx;                                                                 \
    range_check                                                               \
    return NC_NOERR;                                                          \
}
/*----< ncmpix_get_uint_schar() >--------------------------------------------*/
/*----< ncmpix_get_uint_uchar() >--------------------------------------------*/
GET_UINT(schar, if (xx > SCHAR_MAX) return NC_ERANGE;)
GET_UINT(uchar, if (xx > UCHAR_MAX) return NC_ERANGE;)
/*----< ncmpix_get_uint_float() >--------------------------------------------*/
/*----< ncmpix_get_uint_double() >-------------------------------------------*/
/*----< ncmpix_get_uint_int64() >--------------------------------------------*/
/*----< ncmpix_get_uint_uint64() >-------------------------------------------*/
GET_UINT(float,)
GET_UINT(double,)
GET_UINT(int64,)
GET_UINT(uint64,)


/*----< ncmpix_put_uint_schar() >--------------------------------------------*/
static int
ncmpix_put_uint_schar(void *xp, const schar *ip)
{
    uchar *cp = (uchar *) xp;
    *cp++ = 0x00;
    *cp++ = 0x00;
    *cp++ = 0x00;
    *cp = (uchar)*ip;

    if (*ip < 0)
        return NC_ERANGE;

    return NC_NOERR;
}

/*----< ncmpix_put_uint_uchar() >--------------------------------------------*/
static int
ncmpix_put_uint_uchar(void *xp, const uchar *ip)
{
    uchar *cp = (uchar *) xp;
    *cp++ = 0x00;
    *cp++ = 0x00;
    *cp++ = 0x00;
    *cp   = *ip;
    return NC_NOERR;
}

/*----< ncmpix_put_uint_short() >--------------------------------------------*/
static int
ncmpix_put_uint_short(void *xp, const short *ip)
{
    ix_uint xx = (ix_uint)(*ip);
    put_ix_uint(xp, &xx);
#if IX_UINT_MAX < SHORT_MAX
    if (*ip > X_UINT_MAX || *ip < 0)
        return NC_ERANGE;
#else
    if (*ip < 0)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}

/*----< ncmpix_put_uint_ushort() >-------------------------------------------*/
static int
ncmpix_put_uint_ushort(void *xp, const ushort *ip)
{
#if SIZEOF_IX_UINT == SIZEOF_USHORT && IX_UINT_MAX == USHORT_MAX
    put_ix_uint(xp, (ix_uint *)ip);
    return NC_NOERR;
#else
    ix_uint xx = (ix_uint)(*ip);
    put_ix_uint(xp, &xx);
#   if IX_UINT_MAX < USHORT_MAX
    if (*ip > X_UINT_MAX)
        return NC_ERANGE;
#   endif
    return NC_NOERR;
#endif
}

/*----< ncmpix_put_uint_int() >----------------------------------------------*/
static int
ncmpix_put_uint_int(void *xp, const int *ip)
{
    ix_uint xx = (ix_uint)(*ip);
    put_ix_uint(xp, &xx);
#if IX_UINT_MAX < INT_MAX
    if (*ip > X_UINT_MAX || *ip < 0)
        return NC_ERANGE;
#else
    if (*ip < 0)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}

#if SIZEOF_IX_UINT != SIZEOF_UINT
/*----< ncmpix_put_uint_uint() >---------------------------------------------*/
static int
ncmpix_put_uint_uint(void *xp, const uint *ip)
{
#if SIZEOF_IX_UINT == SIZEOF_UINT && IX_UINT_MAX == UINT_MAX
    put_ix_uint(xp, (ix_uint *)ip);
    return NC_NOERR;
#else
    ix_uint xx = (ix_uint)(*ip);
    put_ix_uint(xp, &xx);
#   if IX_UINT_MAX < UINT_MAX
    if (*ip > X_UINT_MAX)
        return NC_ERANGE;
#   endif
    return NC_NOERR;
#endif
}
#endif

/*----< ncmpix_put_uint_long() >---------------------------------------------*/
static int
ncmpix_put_uint_long(void *xp, const long *ip)
{
    ix_uint xx = (ix_uint)(*ip);
    put_ix_uint(xp, &xx);
#if IX_UINT_MAX < LONG_MAX
    if (*ip > X_UINT_MAX || *ip < 0)
        return NC_ERANGE;
#else
    if (*ip < 0)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}

#define PUT_UINT(btype)                                                       \
static int                                                                    \
ncmpix_put_uint_##btype(void *xp, const btype *ip)                            \
{                                                                             \
    ix_uint xx = (ix_uint)(*ip);                                              \
    put_ix_uint(xp, &xx);                                                     \
    if ((double)*ip > (double)X_UINT_MAX || *ip < 0)                          \
    /* typecast to double before the comparison. this is because rounding     \
     * error when typecasting to float, eg. X_UINT_MAX 4294967295 and         \
     * (float) 4294967295 will become 4294967296.0000                         \
     */                                                                       \
        return NC_ERANGE;                                                     \
    return NC_NOERR;                                                          \
}
/*----< ncmpix_put_uint_float() >--------------------------------------------*/
/*----< ncmpix_put_uint_double() >-------------------------------------------*/
/*----< ncmpix_put_uint_int64() >--------------------------------------------*/
PUT_UINT(float)
PUT_UINT(double)
PUT_UINT(int64)

/*----< ncmpix_put_uint_uint64() >-------------------------------------------*/
static int
ncmpix_put_uint_uint64(void *xp, const uint64 *ip)
{
    ix_uint xx = (ix_uint)(*ip);
    put_ix_uint(xp, &xx);
    if (*ip > X_UINT_MAX)
        return NC_ERANGE;
    return NC_NOERR;
}


/*---- uint -----------------------------------------------------------------*/

#define GETN_UINT(btype)                                                      \
int                                                                           \
ncmpix_getn_uint_##btype(const void **xpp, MPI_Offset nelems, btype *tp)      \
{                                                                             \
    const char *xp = (const char *) *xpp;                                     \
    int status = NC_NOERR;                                                    \
                                                                              \
    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_UINT, tp++) {               \
        const int lstatus = ncmpix_get_uint_##btype(xp, tp);                  \
        if (lstatus != NC_NOERR) status = lstatus;                            \
    }                                                                         \
                                                                              \
    *xpp = (void *)xp;                                                        \
    return status;                                                            \
}
/*----< ncmpix_getn_uint_schar() >-------------------------------------------*/
/*----< ncmpix_getn_uint_uchar() >-------------------------------------------*/
/*----< ncmpix_getn_uint_short() >-------------------------------------------*/
/*----< ncmpix_getn_uint_int() >---------------------------------------------*/
/*----< ncmpix_getn_uint_float() >-------------------------------------------*/
/*----< ncmpix_getn_uint_double() >------------------------------------------*/
/*----< ncmpix_getn_uint_ushort() >------------------------------------------*/
/*----< ncmpix_getn_uint_int64() >-------------------------------------------*/
/*----< ncmpix_getn_uint_uint64() >------------------------------------------*/
GETN_UINT(schar)
GETN_UINT(uchar)
GETN_UINT(short)
GETN_UINT(int)
GETN_UINT(float)
GETN_UINT(double)
GETN_UINT(ushort)
GETN_UINT(int64)
GETN_UINT(uint64)

/*----< ncmpix_getn_uint_uint() >--------------------------------------------*/
#if X_SIZEOF_UINT == SIZEOF_UINT
/* optimized version */
int
ncmpix_getn_uint_uint(const void **xpp, MPI_Offset nelems, uint *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(tp, *xpp, nelems * sizeof(uint));
# else
    ncmpii_swapn4b(tp, *xpp, nelems);
# endif
    *xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_UINT);
    return NC_NOERR;
}
#else
GETN_UINT(uint)
#endif

/*----< ncmpix_getn_uint_long() >--------------------------------------------*/
#if X_SIZEOF_UINT == SIZEOF_LONG
/* optimized version */
int
ncmpix_getn_uint_long(const void **xpp, MPI_Offset nelems, long *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(tp, *xpp, nelems * sizeof(long));
# else
    ncmpii_swapn4b(tp, *xpp, nelems);
# endif
    *xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_UINT);
    return NC_NOERR;
}
#else
GETN_UINT(long)
#endif

#define PUTN_UINT(btype)                                                      \
int                                                                           \
ncmpix_putn_uint_##btype(void **xpp, MPI_Offset nelems, const btype *tp)      \
{                                                                             \
    char *xp = (char *) *xpp;                                                 \
    int status = NC_NOERR;                                                    \
                                                                              \
    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_UINT, tp++) {               \
        int lstatus = ncmpix_put_uint_##btype(xp, tp);                        \
        if (lstatus != NC_NOERR) status = lstatus;                            \
    }                                                                         \
                                                                              \
    *xpp = (void *)xp;                                                        \
    return status;                                                            \
}
/*----< ncmpix_putn_uint_schar() >-------------------------------------------*/
/*----< ncmpix_putn_uint_uchar() >-------------------------------------------*/
/*----< ncmpix_putn_uint_short() >-------------------------------------------*/
/*----< ncmpix_putn_uint_int() >---------------------------------------------*/
/*----< ncmpix_putn_uint_float() >-------------------------------------------*/
/*----< ncmpix_putn_uint_double() >------------------------------------------*/
/*----< ncmpix_putn_uint_ushort() >------------------------------------------*/
/*----< ncmpix_putn_uint_int64() >-------------------------------------------*/
/*----< ncmpix_putn_uint_uint64() >------------------------------------------*/
PUTN_UINT(schar)
PUTN_UINT(uchar)
PUTN_UINT(short)
PUTN_UINT(int)
PUTN_UINT(float)
PUTN_UINT(double)
PUTN_UINT(ushort)
PUTN_UINT(int64)
PUTN_UINT(uint64)

/*----< ncmpix_putn_uint_uint() >--------------------------------------------*/
#if X_SIZEOF_UINT == SIZEOF_UINT
/* optimized version */
int
ncmpix_putn_uint_uint(void **xpp, MPI_Offset nelems, const uint *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(*xpp, tp, nelems * X_SIZEOF_UINT);
# else
    ncmpii_swapn4b(*xpp, tp, nelems);
# endif
    *xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_UINT);
    return NC_NOERR;
}
#else
PUTN_UINT(uint)
#endif

/*----< ncmpix_putn_uint_long() >--------------------------------------------*/
#if X_SIZEOF_UINT == SIZEOF_LONG
/* optimized version */
int
ncmpix_putn_uint_long(void **xpp, MPI_Offset nelems, const long *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(*xpp, tp, nelems * X_SIZEOF_UINT);
# else
    ncmpii_swapn4b(*xpp, tp, nelems);
# endif
    *xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_UINT);
    return NC_NOERR;
}
#else
PUTN_UINT(long)
#endif

