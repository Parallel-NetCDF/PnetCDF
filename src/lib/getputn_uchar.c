#include "ncx.h"

/*---- uchar ----------------------------------------------------------------*/

/*----< ncmpix_getn_uchar_schar() >------------------------------------------*/
int
ncmpix_getn_uchar_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
    /* file type is uchar, buffer type is schar */
    int status = NC_NOERR;
    uchar *xp = (uchar *) *xpp;

    while (nelems-- > 0) {
        if (*xp > X_SCHAR_MAX) /* uchar might be too big for schar */
            status = NC_ERANGE;
        *tp++ = *xp++;
    }

    return status;
}

/*----< ncmpix_pad_getn_uchar_schar() >--------------------------------------*/
int
ncmpix_pad_getn_uchar_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
    int status = NC_NOERR;
    uchar *xp = (uchar *) *xpp;

    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    while (nelems-- > 0) {
        if (*xp > X_SCHAR_MAX) /* uchar might be too big for schar */
            status = NC_ERANGE;
        *tp++ = *xp++;
    }
    *xpp = (void *)(xp + rndup);

    return status;
}

/*----< ncmpix_getn_uchar_uchar() >------------------------------------------*/
int
ncmpix_getn_uchar_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
    /* file type is uchar, buffer type is uchar */
    memcpy(tp, *xpp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    return NC_NOERR;
}

#define GETN_UCHAR(btype)                                                     \
int                                                                           \
ncmpix_getn_uchar_##btype(const void **xpp, MPI_Offset nelems, btype *tp)     \
{                                                                             \
    /* there is no ENDIANness issue, as uchar is 1 byte */                    \
    uchar *xp = (uchar *) *xpp;                                               \
    while (nelems-- != 0)                                                     \
        *tp++ = *xp++;                                                        \
    *xpp = (const void *)xp;                                                  \
                                                                              \
    return NC_NOERR;                                                          \
}

/*----< ncmpix_getn_uchar_short() >------------------------------------------*/
/*----< ncmpix_getn_uchar_int() >--------------------------------------------*/
/*----< ncmpix_getn_uchar_long() >-------------------------------------------*/
/*----< ncmpix_getn_uchar_float() >------------------------------------------*/
/*----< ncmpix_getn_uchar_double() >-----------------------------------------*/
/*----< ncmpix_getn_uchar_ushort() >-----------------------------------------*/
/*----< ncmpix_getn_uchar_uint() >-------------------------------------------*/
/*----< ncmpix_getn_uchar_int64() >------------------------------------------*/
/*----< ncmpix_getn_uchar_uint64() >-----------------------------------------*/
GETN_UCHAR(short)
GETN_UCHAR(int)
GETN_UCHAR(long)
GETN_UCHAR(float)
GETN_UCHAR(double)
GETN_UCHAR(ushort)
GETN_UCHAR(uint)
GETN_UCHAR(int64)
GETN_UCHAR(uint64)

/*----< ncmpix_pad_getn_uchar_uchar() >--------------------------------------*/
int
ncmpix_pad_getn_uchar_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    memcpy(tp, *xpp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems + rndup);

    return NC_NOERR;
}

#define PAD_GETN_UCHAR(btype)                                                 \
int                                                                           \
ncmpix_pad_getn_uchar_##btype(const void **xpp, MPI_Offset nelems, btype *tp) \
{                                                                             \
    MPI_Offset rndup = nelems % X_ALIGN;                                      \
    if (rndup) rndup = X_ALIGN - rndup;                                       \
                                                                              \
    /* there is no ENDIANness issue, as uchar is 1 byte */                    \
    uchar *xp = (uchar *) *xpp;                                               \
    while (nelems-- != 0)                                                     \
        *tp++ = *xp++;                                                        \
    *xpp = (void *)((char *)(*xpp) + rndup);                                  \
                                                                              \
    return NC_NOERR;                                                          \
}

/*----< ncmpix_pad_getn_uchar_short() >--------------------------------------*/
/*----< ncmpix_pad_getn_uchar_int() >----------------------------------------*/
/*----< ncmpix_pad_getn_uchar_long() >---------------------------------------*/
/*----< ncmpix_pad_getn_uchar_float() >--------------------------------------*/
/*----< ncmpix_pad_getn_uchar_double() >-------------------------------------*/
/*----< ncmpix_pad_getn_uchar_ushort() >-------------------------------------*/
/*----< ncmpix_pad_getn_uchar_uint() >---------------------------------------*/
/*----< ncmpix_pad_getn_uchar_int64() >--------------------------------------*/
/*----< ncmpix_pad_getn_uchar_uint64() >-------------------------------------*/
PAD_GETN_UCHAR(short)
PAD_GETN_UCHAR(int)
PAD_GETN_UCHAR(long)
PAD_GETN_UCHAR(float)
PAD_GETN_UCHAR(double)
PAD_GETN_UCHAR(ushort)
PAD_GETN_UCHAR(uint)
PAD_GETN_UCHAR(int64)
PAD_GETN_UCHAR(uint64)


/*----< ncmpix_putn_uchar_uchar() >------------------------------------------*/
int
ncmpix_putn_uchar_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
    memcpy(*xpp, tp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    return NC_NOERR;
}

/*----< ncmpix_pad_putn_uchar_uchar() >--------------------------------------*/
int
ncmpix_pad_putn_uchar_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    memcpy(*xpp, tp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    if (rndup) {
        memcpy(*xpp, nada, rndup);
        *xpp = (void *)((char *)(*xpp) + rndup);
    }

    return NC_NOERR;
}

#define PUTN_UCHAR(btype, range_check)                                        \
int                                                                           \
ncmpix_putn_uchar_##btype(void **xpp, MPI_Offset nelems, const btype *tp)     \
{                                                                             \
    int status=NC_NOERR;                                                      \
    uchar *xp = (uchar *) *xpp;                                               \
                                                                              \
    while (nelems-- != 0) {                                                   \
        range_check           /* check if can fit into uchar */               \
        *xp++ = (uchar) *tp++;                                                \
    }                                                                         \
    *xpp = (void *)xp;                                                        \
                                                                              \
    return status;                                                            \
}
/*----< ncmpix_putn_uchar_schar() >------------------------------------------*/
/*----< ncmpix_putn_uchar_short() >------------------------------------------*/
/*----< ncmpix_putn_uchar_int() >--------------------------------------------*/
/*----< ncmpix_putn_uchar_long() >-------------------------------------------*/
/*----< ncmpix_putn_uchar_float() >------------------------------------------*/
/*----< ncmpix_putn_uchar_double() >-----------------------------------------*/
/*----< ncmpix_putn_uchar_int64() >------------------------------------------*/
PUTN_UCHAR(schar,  if (*tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(short,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(int,    if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(long,   if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(float,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(double, if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(int64,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
/*----< ncmpix_putn_uchar_ushort() >-----------------------------------------*/
/*----< ncmpix_putn_uchar_uint() >-------------------------------------------*/
/*----< ncmpix_putn_uchar_uint64() >-----------------------------------------*/
PUTN_UCHAR(ushort, if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)
PUTN_UCHAR(uint,   if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)
PUTN_UCHAR(uint64, if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)

#define PAD_PUTN_UCHAR(btype, range_check)                                    \
int                                                                           \
ncmpix_pad_putn_uchar_##btype(void **xpp, MPI_Offset nelems, const btype *tp) \
{                                                                             \
    int status=NC_NOERR;                                                      \
    uchar *xp = (uchar *) *xpp;                                               \
                                                                              \
    MPI_Offset rndup = nelems % X_ALIGN;                                      \
    if (rndup) rndup = X_ALIGN - rndup;                                       \
                                                                              \
    while (nelems-- != 0) {                                                   \
        range_check         /* check if can fit into uchar */                 \
        *xp++ = (uchar) *tp++;                                                \
    }                                                                         \
                                                                              \
    if (rndup) {                                                              \
        memcpy(xp, nada, rndup);                                              \
        xp += rndup;                                                          \
    }                                                                         \
    *xpp = (void *)xp;                                                        \
                                                                              \
    return status;                                                            \
}
/*----< ncmpix_pad_putn_uchar_schar() >--------------------------------------*/
/*----< ncmpix_pad_putn_uchar_short() >--------------------------------------*/
/*----< ncmpix_pad_putn_uchar_int() >----------------------------------------*/
/*----< ncmpix_pad_putn_uchar_long() >---------------------------------------*/
/*----< ncmpix_pad_putn_uchar_float() >--------------------------------------*/
/*----< ncmpix_pad_putn_uchar_double() >-------------------------------------*/
/*----< ncmpix_pad_putn_uchar_int64() >--------------------------------------*/
PAD_PUTN_UCHAR(schar,  if (*tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(short,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(int,    if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(long,   if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(float,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(double, if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(int64,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
/*----< ncmpix_pad_putn_uchar_ushort() >-------------------------------------*/
/*----< ncmpix_pad_putn_uchar_uint() >---------------------------------------*/
/*----< ncmpix_pad_putn_uchar_uint64() >-------------------------------------*/
PAD_PUTN_UCHAR(ushort, if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)
PAD_PUTN_UCHAR(uint,   if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)
PAD_PUTN_UCHAR(uint64, if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)


