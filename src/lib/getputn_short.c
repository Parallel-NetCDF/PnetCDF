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

/*---- x_short --------------------------------------------------------------*/

#if SHORT_MAX == X_SHORT_MAX
    typedef short ix_short;
    #define SIZEOF_IX_SHORT  SIZEOF_SHORT
    #define IX_SHORT_MAX     SHORT_MAX
#elif INT_MAX >= X_SHORT_MAX
    typedef int ix_short;
    #define SIZEOF_IX_SHORT  SIZEOF_INT
    #define IX_SHORT_MAX     INT_MAX
#elif LONG_MAX >= X_SHORT_MAX
    typedef long ix_short;
    #define SIZEOF_IX_SHORT  SIZEOF_LONG
    #define IX_SHORT_MAX     LONG_MAX
#else
    #error "ix_short ix_ushort implementation"
#endif

static void
get_ix_short(const void *xp, ix_short *ip)
{
    const uchar *cp = (const uchar *) xp;
    *ip = *cp++ << 8;
#if SIZEOF_IX_SHORT > X_SIZEOF_SHORT
    if (*ip & 0x8000) {
        /* extern is negative */
        *ip |= (~(0xffff)); /* N.B. Assumes "twos complement" */
    }
#endif
    *ip |= *cp; 
}

static void
put_ix_short(void *xp, const ix_short *ip)
{
    uchar *cp = (uchar *) xp;
    *cp++ = (*ip) >> 8;
    *cp = (*ip) & 0xff;
}


#if X_SIZEOF_SHORT != SIZEOF_SHORT
/*----< ncmpix_get_short_short() >-------------------------------------------*/
static int
ncmpix_get_short_short(const void *xp, short *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_SHORT && IX_SHORT_MAX == SHORT_MAX
    get_ix_short(xp, (ix_short *)ip);
    return NC_NOERR;
#else
    ix_short xx;
    get_ix_short(xp, &xx);
    *ip = xx;
#   if IX_SHORT_MAX > SHORT_MAX
    if (xx > SHORT_MAX || xx < SHORT_MIN)
        return NC_ERANGE;
#   endif
    return NC_NOERR;
#endif
}
#endif

/*----< ncmpix_get_short_int() >---------------------------------------------*/
int
ncmpix_get_short_int(const void *xp, int *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_INT && IX_SHORT_MAX == INT_MAX
    get_ix_short(xp, (ix_short *)ip);
    return NC_NOERR;
#else
    ix_short xx;
    get_ix_short(xp, &xx);
    *ip = xx;
#   if IX_SHORT_MAX > INT_MAX
    if (xx > INT_MAX || xx < INT_MIN)
        return NC_ERANGE;
#   endif
    return NC_NOERR;
#endif
}

/*----< ncmpix_get_short_long() >--------------------------------------------*/
static int
ncmpix_get_short_long(const void *xp, long *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_LONG && IX_SHORT_MAX == LONG_MAX
    get_ix_short(xp, (ix_short *)ip);
    return NC_NOERR;
#else
    /* assert(LONG_MAX >= X_SHORT_MAX); */
    ix_short xx;
    get_ix_short(xp, &xx);
    *ip = xx;
    return NC_NOERR;
#endif
}

#define GET_SHORT(btype, range_check)                                         \
static int                                                                    \
ncmpix_get_short_##btype(const void *xp, btype *ip)                           \
{                                                                             \
    ix_short xx;                                                              \
    get_ix_short(xp, &xx); /* get a short in the form of local Endianness */  \
    *ip = xx;              /* typecast to btype */                            \
    range_check                                                               \
    return NC_NOERR;                                                            \
}
/*----< ncmpix_get_short_schar() >-------------------------------------------*/
/*----< ncmpix_get_short_uchar() >-------------------------------------------*/
GET_SHORT(schar,  if (xx > SCHAR_MAX || xx < SCHAR_MIN) return NC_ERANGE;)
GET_SHORT(uchar,  if (xx > UCHAR_MAX || xx < 0)         return NC_ERANGE;)
/*----< ncmpix_get_short_float() >-------------------------------------------*/
/*----< ncmpix_get_short_double() >------------------------------------------*/
/*----< ncmpix_get_short_int64() >-------------------------------------------*/
GET_SHORT(float,  )
GET_SHORT(double, )
GET_SHORT(int64,  )
/*----< ncmpix_get_short_ushort() >------------------------------------------*/
/*----< ncmpix_get_short_uint() >--------------------------------------------*/
/*----< ncmpix_get_short_uint64() >------------------------------------------*/
GET_SHORT(ushort, if (xx < 0) return NC_ERANGE;)
GET_SHORT(uint,   if (xx < 0) return NC_ERANGE;)
GET_SHORT(uint64, if (xx < 0) return NC_ERANGE;)


/*----< ncmpix_put_short_schar() >-------------------------------------------*/
static int
ncmpix_put_short_schar(void *xp, const schar *ip)
{
    uchar *cp = (uchar *) xp;
    
    /* copy the signed bit from schar to short */
    if (*ip & 0x80)    /* 0x80 = 10000000(bin) = -127(dec) */
        /* ip is negative */
        *cp++ = 0xff;  /* 0xff = 11111111(bin) = -0(dec) */
        /* now the higher 8 bits are all 1s */
    else
        /* ip is positive */
        *cp++ = 0;

    *cp = (uchar)*ip;  /* the lower 8-bits */
    return NC_NOERR;
}

static int
ncmpix_put_short_uchar(void *xp, const uchar *ip)
{
    uchar *cp = (uchar *) xp;
    *cp++ = 0;
    *cp = *ip;
    return NC_NOERR;
}

#if X_SIZEOF_SHORT != SIZEOF_SHORT
static int
ncmpix_put_short_short(void *xp, const short *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_SHORT && X_SHORT_MAX == SHORT_MAX
    put_ix_short(xp, (const ix_short *)ip);
    return NC_NOERR;
#else
    ix_short xx = (ix_short)*ip;
    put_ix_short(xp, &xx);
# if X_SHORT_MAX < SHORT_MAX
    if (*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
        return NC_ERANGE;
# endif
    return NC_NOERR;
#endif
}
#endif

static int
ncmpix_put_short_int(void *xp, const int *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_INT && X_SHORT_MAX == INT_MAX
    put_ix_short(xp, (const ix_short *)ip);
    return NC_NOERR;
#else
    ix_short xx = (ix_short)*ip;  /* typecasting int to short */
    put_ix_short(xp, &xx);
# if X_SHORT_MAX < INT_MAX
    if (*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
        return NC_ERANGE;
# endif
    return NC_NOERR;
#endif
}

static int
ncmpix_put_short_long(void *xp, const long *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_LONG && X_SHORT_MAX == LONG_MAX
    put_ix_short(xp, (const ix_short *)ip);
    return NC_NOERR;
#else
    ix_short xx = (ix_short)*ip;
    put_ix_short(xp, &xx);
# if X_SHORT_MAX < LONG_MAX
    if (*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
        return NC_ERANGE;
# endif
    return NC_NOERR;
#endif
}

#define PUT_SHORT(btype, range_check)                                         \
static int                                                                    \
ncmpix_put_short_##btype(void *xp, const btype *ip)                           \
{                                                                             \
    ix_short xx = (ix_short) *ip;                                             \
    put_ix_short(xp, &xx);                                                    \
    range_check                                                               \
    return NC_NOERR;                                                          \
}
/*----< ncmpix_put_short_float() >-------------------------------------------*/
/*----< ncmpix_put_short_double() >------------------------------------------*/
/*----< ncmpix_put_short_int64() >-------------------------------------------*/
PUT_SHORT(float,  if (*ip > X_SHORT_MAX || *ip < X_SHORT_MIN) return NC_ERANGE;)
PUT_SHORT(double, if (*ip > X_SHORT_MAX || *ip < X_SHORT_MIN) return NC_ERANGE;)
PUT_SHORT(int64,  if (*ip > X_SHORT_MAX || *ip < X_SHORT_MIN) return NC_ERANGE;)
/*----< ncmpix_put_short_ushort() >------------------------------------------*/
/*----< ncmpix_put_short_uint() >--------------------------------------------*/
/*----< ncmpix_put_short_uint64() >------------------------------------------*/
PUT_SHORT(ushort, if (*ip > X_SHORT_MAX) return NC_ERANGE;)
PUT_SHORT(uint,   if (*ip > X_SHORT_MAX) return NC_ERANGE;)
PUT_SHORT(uint64, if (*ip > X_SHORT_MAX) return NC_ERANGE;)

/*---- short ----------------------------------------------------------------*/

#define GETN_SHORT(fn, btype, pad)                                            \
int                                                                           \
ncmpix_##fn(const void **xpp, MPI_Offset nelems, btype *tp)                   \
{                                                                             \
    const char *xp = (const char *) *xpp;                                     \
    int status = NC_NOERR;                                                    \
                                                                              \
    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++) {              \
        const int lstatus = ncmpix_get_short_##btype(xp, tp);                 \
        if (lstatus != NC_NOERR) status = lstatus;                            \
    }                                                                         \
                                                                              \
    if (pad && nelems % 2 != 0)                                               \
        xp += X_SIZEOF_SHORT;                                                 \
                                                                              \
    *xpp = (void *)xp;                                                        \
    return status;                                                            \
}

/*----< ncmpix_getn_short_schar() >------------------------------------------*/
/*----< ncmpix_getn_short_uchar() >------------------------------------------*/
/*----< ncmpix_getn_short_int() >--------------------------------------------*/
/*----< ncmpix_getn_short_long() >-------------------------------------------*/
/*----< ncmpix_getn_short_float() >------------------------------------------*/
/*----< ncmpix_getn_short_double() >-----------------------------------------*/
/*----< ncmpix_getn_short_ushort() >-----------------------------------------*/
/*----< ncmpix_getn_short_uint() >-------------------------------------------*/
/*----< ncmpix_getn_short_int64() >------------------------------------------*/
/*----< ncmpix_getn_short_uint64() >-----------------------------------------*/
GETN_SHORT(    getn_short_schar,  schar,  0)
GETN_SHORT(    getn_short_uchar,  uchar,  0)
GETN_SHORT(    getn_short_int,    int,    0)
GETN_SHORT(    getn_short_long,   long,   0)
GETN_SHORT(    getn_short_float,  float,  0)
GETN_SHORT(    getn_short_double, double, 0)
GETN_SHORT(    getn_short_ushort, ushort, 0)
GETN_SHORT(    getn_short_uint,   uint,   0)
GETN_SHORT(    getn_short_int64,  int64,  0)
GETN_SHORT(    getn_short_uint64, uint64, 0)

/*----< ncmpix_pad_getn_short_schar() >--------------------------------------*/
/*----< ncmpix_pad_getn_short_uchar() >--------------------------------------*/
/*----< ncmpix_pad_getn_short_int() >----------------------------------------*/
/*----< ncmpix_pad_getn_short_long() >---------------------------------------*/
/*----< ncmpix_pad_getn_short_float() >--------------------------------------*/
/*----< ncmpix_pad_getn_short_double() >-------------------------------------*/
/*----< ncmpix_pad_getn_short_ushort() >-------------------------------------*/
/*----< ncmpix_pad_getn_short_uint() >---------------------------------------*/
/*----< ncmpix_pad_getn_short_int64() >--------------------------------------*/
/*----< ncmpix_pad_getn_short_uint64() >-------------------------------------*/
GETN_SHORT(pad_getn_short_schar,  schar,  1)
GETN_SHORT(pad_getn_short_uchar,  uchar,  1)
GETN_SHORT(pad_getn_short_int,    int,    1)
GETN_SHORT(pad_getn_short_long,   long,   1)
GETN_SHORT(pad_getn_short_float,  float,  1)
GETN_SHORT(pad_getn_short_double, double, 1)
GETN_SHORT(pad_getn_short_ushort, ushort, 1)
GETN_SHORT(pad_getn_short_uint,   uint,   1)
GETN_SHORT(pad_getn_short_int64,  int64,  1)
GETN_SHORT(pad_getn_short_uint64, uint64, 1)

/*----< ncmpix_getn_short_short() >------------------------------------------*/
/*----< ncmpix_pad_getn_short_short() >--------------------------------------*/
#if X_SIZEOF_SHORT == SIZEOF_SHORT
/* optimized version */
int
ncmpix_getn_short_short(const void **xpp, MPI_Offset nelems, short *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(tp, *xpp, nelems * sizeof(short));
# else
    ncmpii_swapn2b(tp, *xpp, nelems);
# endif
    *xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_SHORT);
    return NC_NOERR;
}
int
ncmpix_pad_getn_short_short(const void **xpp, MPI_Offset nelems, short *tp)
{
    const MPI_Offset rndup = nelems % 2;
# ifdef WORDS_BIGENDIAN
    memcpy(tp, *xpp, nelems * sizeof(short));
# else
    ncmpii_swapn2b(tp, *xpp, nelems);
# endif
    *xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_SHORT + rndup);
    return NC_NOERR;
}
#else
GETN_SHORT(    getn_short_short,  short,  0)
GETN_SHORT(pad_getn_short_short,  short,  1)
#endif

#define PUTN_SHORT(fn, btype, pad)                                            \
int                                                                           \
ncmpix_##fn(void **xpp, MPI_Offset nelems, const btype *tp)                   \
{   /* put tp[nelems] (type btype) to xpp[nelems] (type short) */             \
    char *xp = (char *) *xpp;                                                 \
    int status = NC_NOERR;                                                    \
                                                                              \
    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++) {              \
        int lstatus = ncmpix_put_short_##btype(xp, tp);                       \
        if (lstatus != NC_NOERR) status = lstatus;                            \
    }                                                                         \
                                                                              \
    if (pad && nelems % 2 != 0) {                                             \
        memcpy(xp, nada, X_SIZEOF_SHORT);                                     \
        xp += X_SIZEOF_SHORT;                                                 \
    }                                                                         \
                                                                              \
    *xpp = (void *)xp;                                                        \
    return status;                                                            \
}

/*----< ncmpix_putn_short_schar() >------------------------------------------*/
/*----< ncmpix_putn_short_uchar() >------------------------------------------*/
/*----< ncmpix_putn_short_int() >--------------------------------------------*/
/*----< ncmpix_putn_short_long() >-------------------------------------------*/
/*----< ncmpix_putn_short_float() >------------------------------------------*/
/*----< ncmpix_putn_short_double() >-----------------------------------------*/
/*----< ncmpix_putn_short_ushort() >-----------------------------------------*/
/*----< ncmpix_putn_short_uint() >-------------------------------------------*/
/*----< ncmpix_putn_short_int64() >------------------------------------------*/
/*----< ncmpix_putn_short_uint64() >-----------------------------------------*/
PUTN_SHORT(    putn_short_schar,  schar,  0)
PUTN_SHORT(    putn_short_uchar,  uchar,  0)
PUTN_SHORT(    putn_short_int,    int,    0)
PUTN_SHORT(    putn_short_long,   long,   0)
PUTN_SHORT(    putn_short_float,  float,  0)
PUTN_SHORT(    putn_short_double, double, 0)
PUTN_SHORT(    putn_short_ushort, ushort, 0)
PUTN_SHORT(    putn_short_uint,   uint,   0)
PUTN_SHORT(    putn_short_int64,  int64,  0)
PUTN_SHORT(    putn_short_uint64, uint64, 0)

/*----< ncmpix_pad_putn_short_schar() >--------------------------------------*/
/*----< ncmpix_pad_putn_short_uchar() >--------------------------------------*/
/*----< ncmpix_pad_putn_short_int() >----------------------------------------*/
/*----< ncmpix_pad_putn_short_long() >---------------------------------------*/
/*----< ncmpix_pad_putn_short_float() >--------------------------------------*/
/*----< ncmpix_pad_putn_short_double() >-------------------------------------*/
/*----< ncmpix_pad_putn_short_ushort() >-------------------------------------*/
/*----< ncmpix_pad_putn_short_uint() >---------------------------------------*/
/*----< ncmpix_pad_putn_short_int64() >--------------------------------------*/
/*----< ncmpix_pad_putn_short_uint64() >-------------------------------------*/
PUTN_SHORT(pad_putn_short_schar,  schar,  1)
PUTN_SHORT(pad_putn_short_uchar,  uchar,  1)
PUTN_SHORT(pad_putn_short_int,    int,    1)
PUTN_SHORT(pad_putn_short_long,   long,   1)
PUTN_SHORT(pad_putn_short_float,  float,  1)
PUTN_SHORT(pad_putn_short_double, double, 1)
PUTN_SHORT(pad_putn_short_ushort, ushort, 1)
PUTN_SHORT(pad_putn_short_uint,   uint,   1)
PUTN_SHORT(pad_putn_short_int64,  int64,  1)
PUTN_SHORT(pad_putn_short_uint64, uint64, 1)

/*----< ncmpix_putn_short_short() >------------------------------------------*/
/*----< ncmpix_pad_putn_short_short() >--------------------------------------*/
#if X_SIZEOF_SHORT == SIZEOF_SHORT
/* optimized version */
int
ncmpix_putn_short_short(void **xpp, MPI_Offset nelems, const short *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(*xpp, tp, nelems * X_SIZEOF_SHORT);
# else
    ncmpii_swapn2b(*xpp, tp, nelems);
# endif
    *xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_SHORT);
    return NC_NOERR;
}
int
ncmpix_pad_putn_short_short(void **xpp, MPI_Offset nelems, const short *tp)
{
    const MPI_Offset rndup = nelems % 2;
# ifdef WORDS_BIGENDIAN
    memcpy(*xpp, tp, nelems * X_SIZEOF_SHORT);
# else
    ncmpii_swapn2b(*xpp, tp, nelems);
# endif
    *xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_SHORT + rndup);
    return NC_NOERR;
}
#else
PUTN_SHORT(    putn_short_short,  short,  0)
PUTN_SHORT(pad_putn_short_short,  short,  1)
#endif

