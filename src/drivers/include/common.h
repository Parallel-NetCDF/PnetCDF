/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _COMMON_H
#define _COMMON_H

#include <mpi.h>
#include <pnetcdf.h>

/*
 * Macros for dealing with flag bits.
 */
#define fSet(t, f)      ((t) |=  (f))
#define fClr(t, f)      ((t) &= ~(f))
#define fIsSet(t, f)    ((t) &   (f))
#define fMask(t, f)     ((t) & ~ (f))

#ifndef MAX
#define MAX(mm,nn) (((mm) > (nn)) ? (mm) : (nn))
#endif
#ifndef MIN
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))
#endif

/* useful for aligning memory */
#define _RNDUP(x, unit)      ((((x) + (unit) - 1) / (unit)) * (unit))
#define _RNDDOWN(x, unit)    ((x) - ((x)%(unit)))

/* #define M_RND_UNIT   (sizeof(double))
 * SIZEOF_DOUBLE is defined in config.h
 */
#define M_RND_UNIT        SIZEOF_DOUBLE
#define M_RNDUP(x)        _RNDUP(x, M_RND_UNIT)
#define M_RNDDOWN(x)      _RNDDOWN(x, M_RND_UNIT)
#define D_RNDUP(x, align) _RNDUP(x, (off_t)(align))

/* for put request less than 4KB, copy it to a buffer and do byte swap there,
 * so if the user buffer is immutable (assuming smaller than 4KB), it will not
 * cause seg fault. Not a perfect solution, but should be sufficient for most
 * of the cases.
 */
#define NC_BYTE_SWAP_BUFFER_SIZE 4096

#ifdef WORDS_BIGENDIAN
#define NEED_BYTE_SWAP(xtype,itype) 0
#else
#define NEED_BYTE_SWAP(xtype,itype)                              \
    ((xtype == NC_CHAR  && itype == MPI_CHAR)           ||       \
     (xtype == NC_BYTE  && itype == MPI_SIGNED_CHAR)    ||       \
     (xtype == NC_UBYTE && itype == MPI_UNSIGNED_CHAR)) ? 0 : 1
#endif

extern void *
NCI_Malloc_fn(size_t size, const int lineno, const char *func,
              const char *filename);

extern void *
NCI_Calloc_fn(size_t nelem, size_t elsize, const int lineno, const char *func,
              const char *filename);

extern void *
NCI_Realloc_fn(void *ptr, size_t size, const int lineno, const char *func,
               const char *filename);

extern void
NCI_Free_fn(void *ptr, const int lineno, const char *func,
            const char *filename);

#define NCI_Malloc(a)    NCI_Malloc_fn(a,__LINE__,__func__,__FILE__)
#define NCI_Calloc(a,b)  NCI_Calloc_fn(a,b,__LINE__,__func__,__FILE__)
#define NCI_Realloc(a,b) NCI_Realloc_fn(a,b,__LINE__,__func__,__FILE__)
#define NCI_Free(a)      NCI_Free_fn(a,__LINE__,__func__,__FILE__)

extern int
ncmpii_inq_malloc_size(size_t *size);

extern int
ncmpii_inq_malloc_max_size(size_t *size);

extern int
ncmpii_inq_malloc_list(void);

extern int
ncmpii_dtype_decode(MPI_Datatype dtype, MPI_Datatype *ptype, int *el_size,
                    MPI_Offset *nelems, int *isderived,
                    int *iscontig_of_ptypes);

extern int
ncmpii_create_imaptype(int ndims, const MPI_Offset *count,
                       const MPI_Offset *imap, MPI_Datatype ptype,
                       MPI_Datatype *imaptype);

extern int
ncmpii_error_mpi2nc(int mpi_errorcode, char *msg);

extern int
ncmpii_error_posix2nc(char *err_msg);

#ifdef ENABLE_ADIOS
extern int
ncmpii_error_adios2nc(int adios_err, char *err_msg);
#endif

extern int
ncmpii_check_name(const char *name, int file_ver);

extern MPI_Datatype
ncmpii_nc2mpitype(nc_type xtype);

extern int
ncmpii_xlen_nc_type(nc_type xtype, int *size);

extern int
ncmpii_buftype_decode(int ndims, nc_type xtype, const MPI_Offset *count,
                      MPI_Offset bufcount, MPI_Datatype buftype,
                      MPI_Datatype *etype, int *esize, MPI_Offset *nelems,
                      MPI_Offset *xnbytes, int *isContig);

extern int
ncmpii_pack(int ndims, const MPI_Offset *count, const MPI_Offset *imap,
            void *buf, MPI_Offset bufcount, MPI_Datatype buftype,
            MPI_Offset *bnelems, MPI_Datatype *ptype, void **cbuf);

#if 0
extern int
ncmpii_put_cast_swap(int format, MPI_Offset nelems, nc_type xtype,
                     MPI_Datatype itype, void *fillp, void *ibuf,
                     int isNewBuf, void **xbuf);

extern int
ncmpii_get_cast_swap(int format, MPI_Offset nelems, nc_type xtype,
                     MPI_Datatype etype, void *buf, void *xbuf, void **ibuf);
#endif

extern int
ncmpii_need_convert(int format, nc_type xtype, MPI_Datatype mpitype);

extern void
ncmpii_in_swapn(void *buf, MPI_Offset nelems, int esize);

extern int
ncmpii_putn_NC_CHAR  (void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_putn_NC_BYTE  (int cdf_ver,
                      void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_UBYTE (void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_SHORT (void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_USHORT(void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_INT   (void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_UINT  (void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_FLOAT (void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_DOUBLE(void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_INT64 (void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);
extern int
ncmpii_putn_NC_UINT64(void *xbuf, const void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype, void *fillp);

extern int
ncmpii_getn_NC_CHAR  (const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_BYTE  (int cdf_ver,
                      const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_UBYTE (const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_SHORT (const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_USHORT(const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_INT   (const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_UINT  (const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_FLOAT (const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_DOUBLE(const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_INT64 (const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);
extern int
ncmpii_getn_NC_UINT64(const void *xbuf, void *buf, MPI_Offset nelems,
                      MPI_Datatype datatype);

extern int
ncmpii_utf8_normalize(const char *str, char **normalp);

extern int
ncmpii_utf8_validate(const char* name);

typedef struct hash_map_node {
    char *key;
    int val;
    struct hash_map_node *next;
} hash_map_node;


typedef struct hash_map {
    unsigned int (*hash)(const char* key);
    hash_map_node **table;
    unsigned int size;
} hash_map;

int hash_map_init(hash_map *map, int size, unsigned int (*hash)(const char* key));
int hash_map_free(hash_map *map);
int hash_map_add(hash_map *map, char *key, int val);
int hash_map_find(hash_map *map, char *key, int *val);

#ifndef HAVE_STRDUP
extern char *strdup(const char *s);
#endif
#ifndef HAVE_STRCASECMP
extern int strcasecmp(const char *s1, const char *s2);
#endif

#endif

