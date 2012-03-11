#include "nc.h"
#include "ncx.h"

/*---- schar ----------------------------------------------------------------*/

/*----< ncmpix_getn_schar_schar() >------------------------------------------*/
ncmpix_getn_schar_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
    memcpy(tp, *xpp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);
    return NC_NOERR;
}

#define GETN_SCHAR(btype, range_check)                                        \
int                                                                           \
ncmpix_getn_schar_##btype(const void **xpp, MPI_Offset nelems, btype *tp)     \
{                                                                             \
    /* file type (variable defined in file) is schar, buffer type is btype */ \
    int    status = NC_NOERR;                                                 \
    schar *xp = (schar *) *xpp;                                               \
                                                                              \
    /* there is no ENDIANness issue, as schar is 1 byte */                    \
    while (nelems-- != 0) {                                                   \
        range_check                                                           \
        *tp++ = *xp++;                                                        \
    }                                                                         \
    *xpp = (void *)xp;                                                        \
                                                                              \
    return status;                                                            \
}

/*----< ncmpix_getn_schar_uchar() >------------------------------------------*/
/*----< ncmpix_getn_schar_ushort() >-----------------------------------------*/
/*----< ncmpix_getn_schar_uint() >-------------------------------------------*/
/*----< ncmpix_getn_schar_uint64() >-----------------------------------------*/
GETN_SCHAR(uchar,  if (*xp < 0) status = NC_ERANGE;)
GETN_SCHAR(ushort, if (*xp < 0) status = NC_ERANGE;)
GETN_SCHAR(uint,   if (*xp < 0) status = NC_ERANGE;)
GETN_SCHAR(uint64, if (*xp < 0) status = NC_ERANGE;)
/*----< ncmpix_getn_schar_short() >------------------------------------------*/
/*----< ncmpix_getn_schar_int() >--------------------------------------------*/
/*----< ncmpix_getn_schar_long() >-------------------------------------------*/
/*----< ncmpix_getn_schar_float() >------------------------------------------*/
/*----< ncmpix_getn_schar_double() >-----------------------------------------*/
/*----< ncmpix_getn_schar_int64() >------------------------------------------*/
GETN_SCHAR(short,)
GETN_SCHAR(int,)
GETN_SCHAR(long,)
GETN_SCHAR(float,)
GETN_SCHAR(double,)
GETN_SCHAR(int64,)


/*----< ncmpix_pad_getn_schar_schar() >--------------------------------------*/
int
ncmpix_pad_getn_schar_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
    /* get n elements of schar from the variable defined as schar (NC_BYTE,
     * or NC_CHAR) in the file */
    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    memcpy(tp, *xpp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems + rndup);

    return NC_NOERR;
}

#define PAD_GETN_SCHAR(btype, range_check)                                    \
int                                                                           \
ncmpix_pad_getn_schar_##btype(const void **xpp, MPI_Offset nelems, btype *tp) \
{                                                                             \
    /* get n elements of "btype" from the variable defined as schar (NC_BYTE, \
     * or NC_CHAR) in the file */                                             \
    int    status = NC_NOERR;                                                 \
    schar *xp = (schar *) *xpp;                                               \
                                                                              \
    MPI_Offset rndup = nelems % X_ALIGN;                                      \
    if (rndup) rndup = X_ALIGN - rndup;                                       \
                                                                              \
    /* there is no ENDIANness issue, as schar is 1 byte */                    \
    while (nelems-- != 0) {                                                   \
        range_check         /* check if can fit into btype */                 \
        *tp++ = *xp++;                                                        \
    }                                                                         \
    *xpp = (void *)(xp + rndup);                                              \
                                                                              \
    return status;                                                            \
}
/*----< ncmpix_pad_getn_schar_uchar() >--------------------------------------*/
/*----< ncmpix_pad_getn_schar_ushort() >-------------------------------------*/
/*----< ncmpix_pad_getn_schar_uint() >---------------------------------------*/
/*----< ncmpix_pad_getn_schar_uint64() >-------------------------------------*/
PAD_GETN_SCHAR(uchar,  if (*xp < 0) status = NC_ERANGE;)
PAD_GETN_SCHAR(ushort, if (*xp < 0) status = NC_ERANGE;)
PAD_GETN_SCHAR(uint,   if (*xp < 0) status = NC_ERANGE;)
PAD_GETN_SCHAR(uint64, if (*xp < 0) status = NC_ERANGE;)
/*----< ncmpix_pad_getn_schar_short() >--------------------------------------*/
/*----< ncmpix_pad_getn_schar_int() >----------------------------------------*/
/*----< ncmpix_pad_getn_schar_long() >---------------------------------------*/
/*----< ncmpix_pad_getn_schar_float() >--------------------------------------*/
/*----< ncmpix_pad_getn_schar_double() >-------------------------------------*/
/*----< ncmpix_pad_getn_schar_int64() >--------------------------------------*/
PAD_GETN_SCHAR(short,)
PAD_GETN_SCHAR(int,)
PAD_GETN_SCHAR(long,)
PAD_GETN_SCHAR(float,)
PAD_GETN_SCHAR(double,)
PAD_GETN_SCHAR(int64,)

int
/*----< ncmpix_putn_schar_schar() >------------------------------------------*/
ncmpix_putn_schar_schar(void **xpp, MPI_Offset nelems, const schar *tp)
{
    memcpy(*xpp, tp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    return NC_NOERR;
}

/*----< ncmpix_pad_putn_schar_schar() >--------------------------------------*/
int
ncmpix_pad_putn_schar_schar(void **xpp, MPI_Offset nelems, const schar *tp)
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


#define PUTN_SCHAR(btype, range_check)                                        \
int                                                                           \
ncmpix_putn_schar_##btype(void **xpp, MPI_Offset nelems, const btype *tp)     \
{                                                                             \
    /* put n elements of data in btype to file type schar */                  \
    int status=NC_NOERR;                                                      \
    schar *xp = (schar *) *xpp;                                               \
                                                                              \
    while (nelems-- != 0) {                                                   \
        range_check         /* check if can fit into schar */                 \
        *xp++ = (schar) *tp++;                                                \
    }                                                                         \
    *xpp = (void *)xp;                                                        \
                                                                              \
    return status;                                                            \
}
/*----< ncmpix_putn_schar_uchar() >------------------------------------------*/
/*----< ncmpix_putn_schar_ushort() >-----------------------------------------*/
/*----< ncmpix_putn_schar_uint() >-------------------------------------------*/
/*----< ncmpix_putn_schar_uint64() >-----------------------------------------*/
PUTN_SCHAR(uchar,  if (*tp > X_SCHAR_MAX) status = NC_ERANGE;)
PUTN_SCHAR(ushort, if (*tp > X_SCHAR_MAX) status = NC_ERANGE;)
PUTN_SCHAR(uint,   if (*tp > X_SCHAR_MAX) status = NC_ERANGE;)
PUTN_SCHAR(uint64, if (*tp > X_SCHAR_MAX) status = NC_ERANGE;)
/*----< ncmpix_putn_schar_short() >------------------------------------------*/
/*----< ncmpix_putn_schar_int() >--------------------------------------------*/
/*----< ncmpix_putn_schar_long() >-------------------------------------------*/
/*----< ncmpix_putn_schar_float() >------------------------------------------*/
/*----< ncmpix_putn_schar_double() >-----------------------------------------*/
/*----< ncmpix_putn_schar_int64() >------------------------------------------*/
PUTN_SCHAR(short,  if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PUTN_SCHAR(int,    if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PUTN_SCHAR(long,   if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PUTN_SCHAR(float,  if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PUTN_SCHAR(double, if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PUTN_SCHAR(int64,  if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)

#define PAD_PUTN_SCHAR(btype, range_check)                                    \
int                                                                           \
ncmpix_pad_putn_schar_##btype(void **xpp, MPI_Offset nelems, const btype *tp) \
{                                                                             \
    /* put n elements of btype data type to the variable defined as schar     \
       in the file */                                                         \
    int status=NC_NOERR;                                                      \
    schar *xp = (schar *) *xpp;                                               \
                                                                              \
    MPI_Offset rndup = nelems % X_ALIGN;                                      \
    if (rndup) rndup = X_ALIGN - rndup;                                       \
                                                                              \
    while (nelems-- != 0) {                                                   \
        range_check             /* check if can fit into schar */             \
        *xp++ = (schar) *tp++;                                                \
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
/*----< ncmpix_pad_putn_schar_uchar() >--------------------------------------*/
/*----< ncmpix_pad_putn_schar_ushort() >-------------------------------------*/
/*----< ncmpix_pad_putn_schar_uint() >---------------------------------------*/
/*----< ncmpix_pad_putn_schar_uint64() >-------------------------------------*/
PAD_PUTN_SCHAR(uchar,  if (*tp > X_SCHAR_MAX) status = NC_ERANGE;)
PAD_PUTN_SCHAR(ushort, if (*tp > X_SCHAR_MAX) status = NC_ERANGE;)
PAD_PUTN_SCHAR(uint,   if (*tp > X_SCHAR_MAX) status = NC_ERANGE;)
PAD_PUTN_SCHAR(uint64, if (*tp > X_SCHAR_MAX) status = NC_ERANGE;)
/*----< ncmpix_pad_putn_schar_short() >--------------------------------------*/
/*----< ncmpix_pad_putn_schar_int() >----------------------------------------*/
/*----< ncmpix_pad_putn_schar_long() >---------------------------------------*/
/*----< ncmpix_pad_putn_schar_float() >--------------------------------------*/
/*----< ncmpix_pad_putn_schar_double() >-------------------------------------*/
/*----< ncmpix_pad_putn_schar_int64() >--------------------------------------*/
PAD_PUTN_SCHAR(short,  if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PAD_PUTN_SCHAR(int,    if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PAD_PUTN_SCHAR(long,   if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PAD_PUTN_SCHAR(float,  if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PAD_PUTN_SCHAR(double, if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)
PAD_PUTN_SCHAR(int64,  if (*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN) status = NC_ERANGE;)

// PAD_PUTN_SCHAR(uchar,  NO_CHECK_SCHAR_RANGE)
/* wkliao: In netcdf3, there is no no range check for schar_uchar case. This
 * does not seem right ...
 */

// PUTN_SCHAR(uchar,  (0))
/* wkliao: In netcdf3, there is no no range check for schar_uchar case. That
 * does not seem right ...
 */
