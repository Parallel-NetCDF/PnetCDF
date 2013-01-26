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

/*---- x_double -------------------------------------------------------------*/

#if X_SIZEOF_DOUBLE == SIZEOF_DOUBLE  && !defined(NO_IEEE_FLOAT)

static void
get_ix_double(const void *xp, double *ip)
{
#ifdef WORDS_BIGENDIAN
    memcpy(ip, xp, sizeof(double));
#else
    ncmpii_swap8b(ip, xp);
#endif
}

static void
put_ix_double(void *xp, const double *ip)
{
#ifdef WORDS_BIGENDIAN
    memcpy(xp, ip, X_SIZEOF_DOUBLE);
#else
    ncmpii_swap8b(xp, ip);
#endif
}

#elif vax

/* What IEEE double precision floating point looks like on a Vax */
struct ieee_double {
    unsigned int exp_hi   : 7;
    unsigned int sign     : 1;
    unsigned int mant_6   : 4;
    unsigned int exp_lo   : 4;
    unsigned int mant_5   : 8;
    unsigned int mant_4   : 8;
    unsigned int mant_lo  : 32;
};

/* Vax double precision floating point */
struct vax_double {
    unsigned int mantissa1 : 7;
    unsigned int exp       : 8;
    unsigned int sign      : 1;
    unsigned int mantissa2 : 16;
    unsigned int mantissa3 : 16;
    unsigned int mantissa4 : 16;
};

#define VAX_DBL_BIAS   0x81
#define IEEE_DBL_BIAS  0x3ff
#define MASK(nbits)    ((1 << nbits) - 1)

static const struct dbl_limits {
    struct    vax_double d;
    struct    ieee_double ieee;
} dbl_limits[2] = {
    {{ 0x7f, 0xff, 0x0, 0xffff, 0xffff, 0xffff },   /* Max Vax */
     { 0x7f, 0x0,  0x0, 0xf,    0x0,    0x0, 0x0}}, /* Max IEEE */
    {{ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},               /* Min Vax */
     { 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}},         /* Min IEEE */
};


static void
get_ix_double(const void *xp, double *ip)
{
    struct vax_double *const vdp = (struct vax_double *)ip;
    const struct ieee_double *const idp = (const struct ieee_double *) xp;
    const struct dbl_limits *lim;
    int ii;
    for (ii = 0, lim = dbl_limits;
         ii < sizeof(dbl_limits)/sizeof(struct dbl_limits);
         ii++, lim++) {
         if ((idp->mant_lo == lim->ieee.mant_lo) &&
             (idp->mant_4 == lim->ieee.mant_4) &&
             (idp->mant_5 == lim->ieee.mant_5) &&
             (idp->mant_6 == lim->ieee.mant_6) &&
             (idp->exp_lo == lim->ieee.exp_lo) &&
             (idp->exp_hi == lim->ieee.exp_hi)) {
             *vdp = lim->d;
             goto doneit;
        }
    }
    {
        unsigned exp = idp->exp_hi << 4 | idp->exp_lo;
        vdp->exp = exp - IEEE_DBL_BIAS + VAX_DBL_BIAS;
    }
    {
        unsigned mant_hi = (idp->mant_6 << 16)
                         | (idp->mant_5 << 8)
                         | idp->mant_4;
        unsigned mant_lo = SWAP4(idp->mant_lo);
        vdp->mantissa1 = (mant_hi >> 13);
        vdp->mantissa2 = ((mant_hi & MASK(13)) << 3) | (mant_lo >> 29);
        vdp->mantissa3 = (mant_lo >> 13);
        vdp->mantissa4 = (mant_lo << 3);
    }
    doneit:
        vdp->sign = idp->sign;
}

static void
put_ix_double(void *xp, const double *ip)
{
    const struct vax_double *const vdp = (const struct vax_double *)ip;
    struct ieee_double *const idp = (struct ieee_double *) xp;

    if ((vdp->mantissa4 > (dbl_limits[0].d.mantissa4 - 3)) &&
        (vdp->mantissa3 == dbl_limits[0].d.mantissa3) &&
        (vdp->mantissa2 == dbl_limits[0].d.mantissa2) &&
        (vdp->mantissa1 == dbl_limits[0].d.mantissa1) &&
        (vdp->exp == dbl_limits[0].d.exp)) {
        *idp = dbl_limits[0].ieee;
        goto shipit;
    }
    if ((vdp->mantissa4 == dbl_limits[1].d.mantissa4) &&
        (vdp->mantissa3 == dbl_limits[1].d.mantissa3) &&
        (vdp->mantissa2 == dbl_limits[1].d.mantissa2) &&
        (vdp->mantissa1 == dbl_limits[1].d.mantissa1) &&
        (vdp->exp == dbl_limits[1].d.exp)) {
        *idp = dbl_limits[1].ieee;
        goto shipit;
    }

    {
        unsigned exp = vdp->exp - VAX_DBL_BIAS + IEEE_DBL_BIAS;

        unsigned mant_lo = ((vdp->mantissa2 & MASK(3)) << 29) |
                            (vdp->mantissa3 << 13) |
                            ((vdp->mantissa4 >> 3) & MASK(13));

        unsigned mant_hi = (vdp->mantissa1 << 13)
                 | (vdp->mantissa2 >> 3);

        if ((vdp->mantissa4 & 7) > 4) {
            /* round up */
            mant_lo++;
            if (mant_lo == 0) {
                mant_hi++;
                if (mant_hi > 0xffffff) {
                    mant_hi = 0;
                    exp++;
                }
            }
        }

        idp->mant_lo = SWAP4(mant_lo);
        idp->mant_6  = mant_hi >> 16;
        idp->mant_5  = (mant_hi & 0xff00) >> 8;
        idp->mant_4  = mant_hi;
        idp->exp_hi  = exp >> 4;
        idp->exp_lo  = exp;
    }
        
    shipit:
        idp->sign = vdp->sign;

}

    /* vax */
#else
#error "ix_double implementation"
#endif

#define GET_DOUBLE(btype, range_check)                                        \
static int                                                                    \
ncmpix_get_double_##btype(const void *xp, btype *ip)                          \
{                                                                             \
    double xx=0.0;                                                            \
    get_ix_double(xp, &xx);                                                   \
    *ip = xx;                                                                 \
    range_check             /* check if can fit into btype */                 \
    return NC_NOERR;                                                          \
}
/*----< ncmpix_get_double_schar() >------------------------------------------*/
/*----< ncmpix_get_double_short() >------------------------------------------*/
/*----< ncmpix_get_double_int() >--------------------------------------------*/
/*----< ncmpix_get_double_long() >-------------------------------------------*/
/*----< ncmpix_get_double_float() >------------------------------------------*/
/*----< ncmpix_get_double_int64() >------------------------------------------*/
GET_DOUBLE(schar,  if (xx > SCHAR_MAX  || xx < SCHAR_MIN) return NC_ERANGE;)
GET_DOUBLE(short,  if (xx > SHRT_MAX   || xx < SHRT_MIN)  return NC_ERANGE;)
GET_DOUBLE(int,    if (xx > INT_MAX    || xx < INT_MIN)   return NC_ERANGE;)
GET_DOUBLE(long,   if (xx > LONG_MAX   || xx < LONG_MIN)  return NC_ERANGE;)
GET_DOUBLE(float,  if (xx > FLT_MAX    || xx < -FLT_MAX)  return NC_ERANGE;)
GET_DOUBLE(int64,  if (xx > LLONG_MAX  || xx < LLONG_MIN) return NC_ERANGE;)
/*----< ncmpix_get_double_uchar() >------------------------------------------*/
/*----< ncmpix_get_double_ushort() >-----------------------------------------*/
/*----< ncmpix_get_double_uint() >-------------------------------------------*/
/*----< ncmpix_get_double_uint64() >-----------------------------------------*/
GET_DOUBLE(uchar,  if (xx > UCHAR_MAX  || xx < 0) return NC_ERANGE;)
GET_DOUBLE(ushort, if (xx > USHRT_MAX  || xx < 0) return NC_ERANGE;)
GET_DOUBLE(uint,   if (xx > UINT_MAX   || xx < 0) return NC_ERANGE;)
GET_DOUBLE(uint64, if (xx > ULLONG_MAX || xx < 0) return NC_ERANGE;)

static int
ncmpix_get_double_double(const void *xp, double *ip)
{
    get_ix_double(xp, ip);
    return NC_NOERR;
}


#define PUT_DOUBLE(btype)                                                     \
static int                                                                    \
ncmpix_put_double_##btype(void *xp, const btype *ip)                          \
{                                                                             \
    double xx = (double) *ip;                                                 \
    put_ix_double(xp, &xx);                                                   \
    return NC_NOERR;                                                          \
}
/*----< ncmpix_put_double_schar() >------------------------------------------*/
/*----< ncmpix_put_double_uchar() >------------------------------------------*/
/*----< ncmpix_put_double_short() >------------------------------------------*/
/*----< ncmpix_put_double_int() >--------------------------------------------*/
/*----< ncmpix_put_double_long() >-------------------------------------------*/
/*----< ncmpix_put_double_float() >------------------------------------------*/
/*----< ncmpix_put_double_ushort() >-----------------------------------------*/
/*----< ncmpix_put_double_uint() >-------------------------------------------*/
/*----< ncmpix_put_double_int64() >------------------------------------------*/
/*----< ncmpix_put_double_uint64() >-----------------------------------------*/
PUT_DOUBLE(schar)
PUT_DOUBLE(uchar)
PUT_DOUBLE(short)
PUT_DOUBLE(int)
PUT_DOUBLE(long)
PUT_DOUBLE(float)
PUT_DOUBLE(ushort)
PUT_DOUBLE(uint)
PUT_DOUBLE(int64)
PUT_DOUBLE(uint64)

static int
ncmpix_put_double_double(void *xp, const double *ip)
{
    put_ix_double(xp, ip);
#ifdef NO_IEEE_FLOAT
    if(*ip > X_DOUBLE_MAX || *ip < X_DOUBLE_MIN)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}

/*---- double ---------------------------------------------------------------*/

#define GETN_DOUBLE(btype)                                                    \
int                                                                           \
ncmpix_getn_double_##btype(const void **xpp, MPI_Offset nelems, btype *tp)    \
{                                                                             \
    const char *xp = (const char *) *xpp;                                     \
    int status = NC_NOERR;                                                    \
                                                                              \
    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++) {             \
        const int lstatus = ncmpix_get_double_##btype(xp, tp);                \
        if (lstatus != NC_NOERR) status = lstatus;                            \
    }                                                                         \
                                                                              \
    *xpp = (void *)xp;                                                        \
    return status;                                                            \
}
/*----< ncmpix_getn_double_schar() >-----------------------------------------*/
/*----< ncmpix_getn_double_uchar() >-----------------------------------------*/
/*----< ncmpix_getn_double_short() >-----------------------------------------*/
/*----< ncmpix_getn_double_int() >-------------------------------------------*/
/*----< ncmpix_getn_double_long() >------------------------------------------*/
/*----< ncmpix_getn_double_float() >-----------------------------------------*/
/*----< ncmpix_getn_double_ushort() >----------------------------------------*/
/*----< ncmpix_getn_double_uint() >------------------------------------------*/
/*----< ncmpix_getn_double_int64() >-----------------------------------------*/
/*----< ncmpix_getn_double_uint64() >----------------------------------------*/
GETN_DOUBLE(schar)
GETN_DOUBLE(uchar)
GETN_DOUBLE(short)
GETN_DOUBLE(int)
GETN_DOUBLE(long)
GETN_DOUBLE(float)
GETN_DOUBLE(ushort)
GETN_DOUBLE(uint)
GETN_DOUBLE(int64)
GETN_DOUBLE(uint64)

/*----< ncmpix_getn_double_double() >----------------------------------------*/
#if X_SIZEOF_DOUBLE == SIZEOF_DOUBLE && !defined(NO_IEEE_FLOAT)
/* optimized version */
int
ncmpix_getn_double_double(const void **xpp, MPI_Offset nelems, double *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(tp, *xpp, nelems * sizeof(double));
# else
    ncmpii_swapn8b(tp, *xpp, nelems);
# endif
    *xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_DOUBLE);
    return NC_NOERR;
}
#elif vax
int
ncmpix_getn_double_double(const void **xpp, MPI_Offset ndoubles, double *ip)
{
    double *const end = ip + ndoubles;

    while(ip < end) {
        struct vax_double *const vdp = (struct vax_double *)ip;
        const struct ieee_double *const idp = (const struct ieee_double *) (*xpp);
        const struct dbl_limits *lim;
        int ii;
        for (ii = 0, lim = dbl_limits;
             ii < sizeof(dbl_limits)/sizeof(struct dbl_limits);
             ii++, lim++) {
             if ((idp->mant_lo == lim->ieee.mant_lo) &&
                 (idp->mant_4  == lim->ieee.mant_4)  &&
                 (idp->mant_5  == lim->ieee.mant_5)  &&
                 (idp->mant_6  == lim->ieee.mant_6)  &&
                 (idp->exp_lo  == lim->ieee.exp_lo)  &&
                 (idp->exp_hi  == lim->ieee.exp_hi)) {
                 *vdp = lim->d;
                 goto doneit;
             }
        }
        {
            unsigned exp = idp->exp_hi << 4 | idp->exp_lo;
            vdp->exp = exp - IEEE_DBL_BIAS + VAX_DBL_BIAS;
        }
        {
            unsigned mant_hi = ((idp->mant_6 << 16)
                             | (idp->mant_5 << 8)
                             | idp->mant_4);
            unsigned mant_lo = SWAP4(idp->mant_lo);
            vdp->mantissa1 = (mant_hi >> 13);
            vdp->mantissa2 = ((mant_hi & MASK(13)) << 3)
                           | (mant_lo >> 29);
            vdp->mantissa3 = (mant_lo >> 13);
            vdp->mantissa4 = (mant_lo << 3);
        }
    doneit:
        vdp->sign = idp->sign;

        ip++;
        *xpp = (char *)(*xpp) + X_SIZEOF_DOUBLE;
    }
    return NC_NOERR;
}
    /* vax */

#else
GETN_DOUBLE(double)
#endif

#define PUTN_DOUBLE(btype)                                                    \
int                                                                           \
ncmpix_putn_double_##btype(void **xpp, MPI_Offset nelems, const btype *tp)    \
{                                                                             \
    char *xp = (char *) *xpp;                                                 \
    int status = NC_NOERR;                                                    \
                                                                              \
    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++) {             \
        int lstatus = ncmpix_put_double_##btype(xp, tp);                      \
        if (lstatus != NC_NOERR) status = lstatus;                            \
    }                                                                         \
                                                                              \
    *xpp = (void *)xp;                                                        \
    return status;                                                            \
}
/*----< ncmpix_putn_double_schar() >-----------------------------------------*/
/*----< ncmpix_putn_double_uchar() >-----------------------------------------*/
/*----< ncmpix_putn_double_short() >-----------------------------------------*/
/*----< ncmpix_putn_double_int() >-------------------------------------------*/
/*----< ncmpix_putn_double_long() >------------------------------------------*/
/*----< ncmpix_putn_double_float() >-----------------------------------------*/
/*----< ncmpix_putn_double_ushort() >----------------------------------------*/
/*----< ncmpix_putn_double_uint() >------------------------------------------*/
/*----< ncmpix_putn_double_int64() >-----------------------------------------*/
/*----< ncmpix_putn_double_uint64() >----------------------------------------*/
PUTN_DOUBLE(schar)
PUTN_DOUBLE(uchar)
PUTN_DOUBLE(short)
PUTN_DOUBLE(int)
PUTN_DOUBLE(long)
PUTN_DOUBLE(float)
PUTN_DOUBLE(ushort)
PUTN_DOUBLE(uint)
PUTN_DOUBLE(int64)
PUTN_DOUBLE(uint64)

/*----< ncmpix_putn_double_double() >----------------------------------------*/
#if X_SIZEOF_DOUBLE == SIZEOF_DOUBLE && !defined(NO_IEEE_FLOAT)
/* optimized version */
int
ncmpix_putn_double_double(void **xpp, MPI_Offset nelems, const double *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(*xpp, tp, nelems * X_SIZEOF_DOUBLE);
# else
    ncmpii_swapn8b(*xpp, tp, nelems);
# endif
    *xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_DOUBLE);
    return NC_NOERR;
}
#elif vax
int
ncmpix_putn_double_double(void **xpp, MPI_Offset ndoubles, const double *ip)
{
    const double *const end = ip + ndoubles;

    while(ip < end) {
        const struct vax_double *const vdp = (const struct vax_double *)ip;
        struct ieee_double *const idp = (struct ieee_double *) (*xpp);

        if ((vdp->mantissa4 > (dbl_limits[0].d.mantissa4 - 3)) &&
            (vdp->mantissa3 == dbl_limits[0].d.mantissa3) &&
            (vdp->mantissa2 == dbl_limits[0].d.mantissa2) &&
            (vdp->mantissa1 == dbl_limits[0].d.mantissa1) &&
            (vdp->exp == dbl_limits[0].d.exp)) {
            *idp = dbl_limits[0].ieee;
            goto shipit;
        }
        if ((vdp->mantissa4 == dbl_limits[1].d.mantissa4) &&
            (vdp->mantissa3 == dbl_limits[1].d.mantissa3) &&
            (vdp->mantissa2 == dbl_limits[1].d.mantissa2) &&
            (vdp->mantissa1 == dbl_limits[1].d.mantissa1) &&
            (vdp->exp == dbl_limits[1].d.exp)) {
            *idp = dbl_limits[1].ieee;
            goto shipit;
        }

        {
            unsigned exp = vdp->exp - VAX_DBL_BIAS + IEEE_DBL_BIAS;
            unsigned mant_lo = ((vdp->mantissa2 & MASK(3)) << 29) |
                                (vdp->mantissa3 << 13) |
                               ((vdp->mantissa4 >> 3) & MASK(13));
            unsigned mant_hi = (vdp->mantissa1 << 13)
                             | (vdp->mantissa2 >> 3);

            if ((vdp->mantissa4 & 7) > 4) {
                /* round up */
                mant_lo++;
                if (mant_lo == 0) {
                    mant_hi++;
                    if (mant_hi > 0xffffff) {
                        mant_hi = 0;
                        exp++;
                    }
                }
            }

            idp->mant_lo = SWAP4(mant_lo);
            idp->mant_6  = mant_hi >> 16;
            idp->mant_5  = (mant_hi & 0xff00) >> 8;
            idp->mant_4  = mant_hi;
            idp->exp_hi  = exp >> 4;
            idp->exp_lo  = exp;
        }
        
    shipit:
        idp->sign = vdp->sign;

        ip++;
        *xpp = (char *)(*xpp) + X_SIZEOF_DOUBLE;
    }
    return NC_NOERR;
}
    /* vax */

#else
PUTN_DOUBLE(double)
#endif

