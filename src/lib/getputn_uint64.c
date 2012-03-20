/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#include "ncx.h"

/* ftype is the variable's nc_type defined in file, eg. int64
 * btype is the I/O buffer's C data type, eg. long long
 * buftype is I/O bufer's MPI data type, eg. MPI_UNSIGNED_LONG_LONG
 * apitype is data type appeared in the API names, eg. ncmpi_get_vara_longlong
 */

/*---- x_uint64 -------------------------------------------------------------*/

/*
static void
get_ix_uint64(const void *xp, uint64 *ip)
{
    *ip = *((uint64*)xp);
#ifndef WORDS_BIGENDIAN
    SWAP8B(ip);
#endif
}

static void
put_ix_uint64(void *xp, const uint64 *ip)
{
    *((uint64*) xp) = *ip;
#ifndef WORDS_BIGENDIAN
    SWAP8B(xp);
#endif
}
*/

static void
get_ix_uint64(const void *xp, uint64 *ip)
{
    /* are these bit shifting faster than byte swap? */
    const uchar *cp = (const uchar *) xp;

    *ip  = ((uint64)(*cp++) << 56);
    *ip |= ((uint64)(*cp++) << 48);
    *ip |= ((uint64)(*cp++) << 40);
    *ip |= ((uint64)(*cp++) << 32);
    *ip |= ((uint64)(*cp++) << 24);
    *ip |= ((uint64)(*cp++) << 16);
    *ip |= ((uint64)(*cp++) <<  8);
    *ip |=  (uint64)*cp;
}

static void
put_ix_uint64(void *xp, const uint64 *ip)
{
    uchar *cp = (uchar *) xp;

    *cp++ = (*ip) >> 56;
    *cp++ = ((*ip) & 0x00ff000000000000ULL) >> 48;
    *cp++ = ((*ip) & 0x0000ff0000000000ULL) >> 40;
    *cp++ = ((*ip) & 0x000000ff00000000ULL) >> 32;
    *cp++ = ((*ip) & 0x00000000ff000000ULL) >> 24;
    *cp++ = ((*ip) & 0x0000000000ff0000ULL) >> 16;
    *cp++ = ((*ip) & 0x000000000000ff00ULL) >>  8;
    *cp   = ((*ip) & 0x00000000000000ffULL);
}

#define GET_UINT64(btype, range_check)                                        \
static int                                                                    \
ncmpix_get_uint64_##btype(const void *xp, btype *ip)                          \
{                                                                             \
    uint64 xx;                                                                \
    get_ix_uint64(xp, &xx);                                                   \
    *ip = xx;                                                                 \
    range_check       /* check if can fit into btype */                       \
    return NC_NOERR;                                                          \
}
/* for smaller-sized   signed types, check if the got uint64 is too big (schar, short, int, long, float)
 * for smaller-sized unsigned types, check if the got uint64 is too big (uchar, ushort, uint)
 * for equal-sized     signed types, check if the got uint64 is too big (double, int64)
 * for equal-sized   unsigned types, no check is needed (uint64)
 */
/*----< ncmpix_get_uint64_schar() >------------------------------------------*/
/*----< ncmpix_get_uint64_uchar() >------------------------------------------*/
/*----< ncmpix_get_uint64_short() >------------------------------------------*/
/*----< ncmpix_get_uint64_int() >--------------------------------------------*/
/*----< ncmpix_get_uint64_long() >-------------------------------------------*/
/*----< ncmpix_get_uint64_float() >------------------------------------------*/
/*----< ncmpix_get_uint64_double() >-----------------------------------------*/
/*----< ncmpix_get_uint64_ushort() >-----------------------------------------*/
/*----< ncmpix_get_uint64_uint() >-------------------------------------------*/
/*----< ncmpix_get_uint64_int64() >------------------------------------------*/
/*----< ncmpix_get_uint64_uint64() >-----------------------------------------*/
GET_UINT64(schar,  if (xx > SCHAR_MAX) return NC_ERANGE;)
GET_UINT64(uchar,  if (xx > UCHAR_MAX) return NC_ERANGE;)
GET_UINT64(short,  if (xx > SHRT_MAX)  return NC_ERANGE;)
GET_UINT64(int,    if (xx > INT_MAX)   return NC_ERANGE;)
#if SIZEOF_LONG == X_SIZEOF_INT
static int
ncmpix_get_uint64_long(const void *xp, long *ip)
{
    return ncmpix_get_uint64_int(xp, (int*)ip);
}
#else
GET_UINT64(long,   if (xx > LONG_MAX)  return NC_ERANGE;)
#endif
GET_UINT64(float,  if (xx > FLT_MAX)   return NC_ERANGE;)
GET_UINT64(double, if (xx > DBL_MAX)   return NC_ERANGE;)
GET_UINT64(ushort, if (xx > USHRT_MAX) return NC_ERANGE;)
GET_UINT64(uint,   if (xx > UINT_MAX)  return NC_ERANGE;)
GET_UINT64(int64,  if (xx > LLONG_MAX) return NC_ERANGE;)
GET_UINT64(uint64,)

#define PUT_UINT64(btype, range_check)                                        \
static int                                                                    \
ncmpix_put_uint64_##btype(void *xp, const btype *ip)                          \
{                                                                             \
    uint64 xx = (uint64) *ip;                                                 \
    put_ix_uint64(xp, &xx);                                                   \
    range_check       /* check if can fit into uint64 */                      \
    return NC_NOERR;                                                          \
}
/* for smaller-sized   signed types, check if the put value is negative (schar, short, int, long, float)
 * for smaller-sized unsigned types, no check is needed (uchar, ushort, uint)
 * for equal-sized     signed types, check if the put value is negative (double, int64)
 * for equal-sized   unsigned types, no check is needed (uint64)
 */
/*----< ncmpix_put_uint64_uchar() >------------------------------------------*/
/*----< ncmpix_put_uint64_ushort() >-----------------------------------------*/
/*----< ncmpix_put_uint64_uint() >-------------------------------------------*/
/*----< ncmpix_put_uint64_uint64() >-----------------------------------------*/
PUT_UINT64(uchar,)
PUT_UINT64(uint,)
PUT_UINT64(ushort,)
PUT_UINT64(uint64,)
/*----< ncmpix_put_uint64_schar() >------------------------------------------*/
/*----< ncmpix_put_uint64_short() >------------------------------------------*/
/*----< ncmpix_put_uint64_int() >--------------------------------------------*/
/*----< ncmpix_put_uint64_long() >-------------------------------------------*/
/*----< ncmpix_put_uint64_float() >------------------------------------------*/
/*----< ncmpix_put_uint64_double() >-----------------------------------------*/
/*----< ncmpix_put_uint64_int64() >------------------------------------------*/
PUT_UINT64(schar,  if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(short,  if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(int,    if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(long,   if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(float,  if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(double, if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(int64,  if (*ip < 0) return NC_ERANGE;)


/*---- uint64 ----------------------------------------------------------------*/

#define GETN_UINT64(btype)                                                     \
int                                                                           \
ncmpix_getn_uint64_##btype(const void **xpp, MPI_Offset nelems, btype *tp)      \
{                                                                             \
    const char *xp = (const char *) *xpp;                                     \
    int status = NC_NOERR;                                                    \
                                                                              \
    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_UINT64, tp++) {             \
        const int lstatus = ncmpix_get_uint64_##btype(xp, tp);                 \
        if (lstatus != NC_NOERR) status = lstatus;                            \
    }                                                                         \
                                                                              \
    *xpp = (void *)xp;                                                        \
    return status;                                                            \
}

/*----< ncmpix_getn_uint64_schar() >-----------------------------------------*/
/*----< ncmpix_getn_uint64_uchar() >-----------------------------------------*/
/*----< ncmpix_getn_uint64_short() >-----------------------------------------*/
/*----< ncmpix_getn_uint64_int() >-------------------------------------------*/
/*----< ncmpix_getn_uint64_long() >------------------------------------------*/
/*----< ncmpix_getn_uint64_float() >-----------------------------------------*/
/*----< ncmpix_getn_uint64_double() >----------------------------------------*/
/*----< ncmpix_getn_uint64_ushort() >----------------------------------------*/
/*----< ncmpix_getn_uint64_uint() >------------------------------------------*/
/*----< ncmpix_getn_uint64_int64() >-----------------------------------------*/
GETN_UINT64(schar)
GETN_UINT64(uchar)
GETN_UINT64(short)
GETN_UINT64(int)
GETN_UINT64(long)
GETN_UINT64(float)
GETN_UINT64(double)
GETN_UINT64(ushort)
GETN_UINT64(uint)
GETN_UINT64(int64)

/*----< ncmpix_getn_uint64_uint64() >----------------------------------------*/
/* optimized version */
int
ncmpix_getn_uint64_uint64(const void **xpp, MPI_Offset nelems, uint64 *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(tp, *xpp, nelems * sizeof(uint64));
# else
    ncmpii_swapn8b(tp, *xpp, nelems);
# endif
    *xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_UINT64);
    return NC_NOERR;
}

#define PUTN_UINT64(btype)                                                    \
int                                                                           \
ncmpix_putn_uint64_##btype(void **xpp, MPI_Offset nelems, const btype *tp)    \
{                                                                             \
    char *xp = (char *) *xpp;                                                 \
    int status = NC_NOERR;                                                    \
                                                                              \
    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_UINT64, tp++) {             \
        int lstatus = ncmpix_put_uint64_##btype(xp, tp);                      \
        if (lstatus != NC_NOERR) status = lstatus;                            \
    }                                                                         \
                                                                              \
    *xpp = (void *)xp;                                                        \
    return status;                                                            \
}

/*----< ncmpix_putn_uint64_schar() >-----------------------------------------*/
/*----< ncmpix_putn_uint64_uchar() >-----------------------------------------*/
/*----< ncmpix_putn_uint64_short() >-----------------------------------------*/
/*----< ncmpix_putn_uint64_int() >-------------------------------------------*/
/*----< ncmpix_putn_uint64_long() >------------------------------------------*/
/*----< ncmpix_putn_uint64_float() >-----------------------------------------*/
/*----< ncmpix_putn_uint64_double() >----------------------------------------*/
/*----< ncmpix_putn_uint64_ushort() >----------------------------------------*/
/*----< ncmpix_putn_uint64_uint() >------------------------------------------*/
/*----< ncmpix_putn_uint64_int64() >-----------------------------------------*/
PUTN_UINT64(schar)
PUTN_UINT64(uchar)
PUTN_UINT64(short)
PUTN_UINT64(int)
PUTN_UINT64(long)
PUTN_UINT64(float)
PUTN_UINT64(double)
PUTN_UINT64(ushort)
PUTN_UINT64(uint)
PUTN_UINT64(int64)

/*----< ncmpix_putn_uint64_uint64() >----------------------------------------*/
/* optimized version */
int
ncmpix_putn_uint64_uint64(void **xpp, MPI_Offset nelems, const uint64 *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(*xpp, tp, nelems * X_SIZEOF_UINT64);
# else
    ncmpii_swapn8b(*xpp, tp, nelems);
# endif
    *xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_UINT64);
    return NC_NOERR;
}


