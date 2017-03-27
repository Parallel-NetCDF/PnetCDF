/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <assert.h>
#include <sys/types.h>  /* open() */
#include <sys/stat.h>   /* open() */
#include <fcntl.h>      /* open() */
#include <unistd.h>     /* read() */
#include <string.h>
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <errno.h>

#include <mpi.h>

#include <ncx.h>
#include <macro.h>

/*
 * "magic number" at beginning of file: 0x43444601 (big endian) 
 */
static const schar ncmagic[] = {'C', 'D', 'F', 0x01}; 

/* Prototypes for functions used only in this file */
static int val_get_NCtype(bufferinfo *gbp, NCtype *typep);
static int val_get_size_t(bufferinfo *gbp, MPI_Offset *sp);
static int val_get_NC_string(bufferinfo *gbp, NC_string **ncstrpp);
static int val_get_NC_dim(bufferinfo *gbp, NC_dim **dimpp);
static int val_get_NC_dimarray(bufferinfo *gbp, NC_dimarray *ncap);
static int val_get_nc_type(bufferinfo *gbp, nc_type *typep);
static int val_get_NC_attrV(bufferinfo *gbp, NC_attr *attrp);
static int val_get_NC_attr(bufferinfo *gbp, NC_attr **attrpp);
static int val_get_NC_attrarray(bufferinfo *gbp, NC_attrarray *ncap);
static int val_get_NC_var(bufferinfo *gbp, NC_var **varpp);
static int val_get_NC_vararray(bufferinfo *gbp, NC_vararray *ncap);
static int val_get_NC(NC *ncp);

#define ABSENT 0
#define X_SIZEOF_INT 4
static int x_sizeof_NON_NEG;

#define ABORT(e) { \
    printf("Abort at line %d: %s\n",__LINE__,ncmpi_strerror(e)); \
    MPI_Finalize(); \
    return 1; \
}


/*
 * Fetch the next header chunk.
 */
static int
val_fetch(bufferinfo *gbp) {
    ssize_t nn = 0;
    MPI_Offset slack;        /* any leftover data in the buffer */
    MPI_Aint pos_addr, base_addr;

    assert(gbp->base != NULL);
  
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    slack = gbp->size - (pos_addr - base_addr);
    /* if gbp->pos and gbp->base are the same, there is no leftover buffer data
     * to worry about.  
     * In the other extreme, where gbp->size == (gbp->pos - gbp->base), then all
     * data in the buffer has been consumed */
    if (slack == gbp->size) slack = 0;

    memset(gbp->base, 0, gbp->size);
    gbp->pos = gbp->base;
    gbp->index = 0;

    if (-1 == lseek(gbp->nciop->fd, gbp->offset-slack, SEEK_SET)) {
        printf("Error at line %d: lseek %s\n",__LINE__,strerror(errno));
        return -1;
    }
    nn = read(gbp->nciop->fd, gbp->base, gbp->size);
    if (nn == -1) {
        printf("Error at line %d: lseek %s\n",__LINE__,strerror(errno));
        return -1;
    }
/*
    if (nn < gbp->size) {
        printf("Error: file header size is less than expected\n");
printf("Error: pos_addr=%ld base_addr=%ld gbp->size=%lld nn=%zd\n",pos_addr,base_addr,gbp->size,nn);
        return NC_ENOTNC;
    }
*/
    gbp->offset += (gbp->size - slack);

    return NC_NOERR;
}

/*
 * Ensure that 'nextread' bytes are available.
 */
static int
val_check_buffer(bufferinfo *gbp,
                 MPI_Offset  nextread)
{
    MPI_Aint pos_addr, base_addr;

#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    if (pos_addr + nextread <= base_addr + gbp->size)
        return NC_NOERR;

    return val_fetch(gbp);
} 

static int
val_get_NCtype(bufferinfo *gbp, NCtype *typep)
{
    unsigned int type = 0;
    int status = val_check_buffer(gbp, x_sizeof_NON_NEG);
    if (status != NC_NOERR) {
        printf("NC component type is expected for ");
        return status;
    }

    status = ncmpix_get_uint32((const void**)(&gbp->pos), &type);
    if (status != NC_NOERR) return status;
    *typep = (NCtype) type;
    return NC_NOERR;
}

static int
val_get_size_t(bufferinfo *gbp, MPI_Offset *sp) {
  int sizeof_t = (gbp->version < 5) ? 4 : 8; 
  int status = val_check_buffer(gbp, sizeof_t);
  if (status != NC_NOERR) {
    printf("size is expected for ");
    return status; 
  }
  if (gbp->version < 5) {
      unsigned int tmp=0;
      status = ncmpix_get_uint32((const void **)(&gbp->pos), &tmp);
      *sp = (MPI_Offset)tmp;
  }
  else {
      unsigned long long tmp=0;
      status = ncmpix_get_uint64((const void **)(&gbp->pos), &tmp);
      *sp = (MPI_Offset)tmp;
  }
  return status;
}

static int
val_get_NC_string(bufferinfo *gbp, NC_string **ncstrpp) {
  int status;
  MPI_Offset  nchars = 0, padding, bufremain, strcount; 
  NC_string *ncstrp;
  char *cpos, pad[X_ALIGN-1];
  MPI_Aint pos_addr, base_addr;

  status = val_get_size_t(gbp, &nchars);
  if (status != NC_NOERR) {
    printf("the name string of ");
    return status;
  }

  ncstrp = ncmpii_new_NC_string(nchars, NULL);
  if (ncstrp == NULL)
    return NC_ENOMEM;

  padding = _RNDUP(X_SIZEOF_CHAR * ncstrp->nchars, X_ALIGN)
            - X_SIZEOF_CHAR * ncstrp->nchars;
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
  bufremain = gbp->size - (pos_addr - base_addr);
  cpos = ncstrp->cp;

  while (nchars > 0) {
    if (bufremain > 0) {
      strcount = MIN(bufremain, X_SIZEOF_CHAR * nchars); 
      (void) memcpy(cpos, gbp->pos, strcount);
      nchars -= strcount/X_SIZEOF_CHAR;
      gbp->pos = (void *)((char *)gbp->pos + strcount);
      cpos += strcount; 
      bufremain -= strcount;
    } else {
      status = val_fetch(gbp);
      if(status != NC_NOERR) {
	printf("fetching the name string of ");
        ncmpii_free_NC_string(ncstrp);
        return status;
      } 
      bufremain = gbp->size;
    }
  }

  memset(pad, 0, X_ALIGN-1);
  status = val_check_buffer(gbp, padding);
  if(status != NC_NOERR) {
    printf("fetching padding for the name string of ");
    ncmpii_free_NC_string(ncstrp);
    return status;
  } 
  if (memcmp(gbp->pos, pad, padding) != 0) {
    printf("Error @ [0x%8.8Lx]: \n\tPadding should be 0x00 for the name string alignment of ", (long long unsigned)
	   (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size));
    ncmpii_free_NC_string(ncstrp);
    return NC_EINVAL;
  }
  gbp->pos = (void *)((char *)gbp->pos + padding);
  
  *ncstrpp = ncstrp;
  
  return NC_NOERR;  
}

static int
val_get_NC_dim(bufferinfo *gbp, NC_dim **dimpp) {
    int status;
    NC_string *ncstrp;
    NC_dim *dimp;

    status = val_get_NC_string(gbp, &ncstrp);
    if (status != NC_NOERR) return status;

    dimp = ncmpii_new_x_NC_dim(ncstrp);
    if(dimp == NULL) return NC_ENOMEM;

    status = val_get_size_t(gbp, &dimp->size);
    if (status != NC_NOERR) {
        printf("\"%s\" - ", ncstrp->cp);
        ncmpii_free_NC_dim(dimp); /* frees name */
        return status;
    }

    *dimpp = dimp;

    return NC_NOERR;
}

static int
val_get_NC_dimarray(bufferinfo *gbp, NC_dimarray *ncap)
{
    /* netCDF file format:
     *  ...
     * dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_DIMENSION = \x00 \x00 \x00 \x0A         // tag for list of dimensions
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int status;
    NCtype type = NC_UNSPECIFIED; 
    NC_dim **dpp, **end;
    int dim;
    unsigned long long err_addr;
    MPI_Offset tmp;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    status = val_get_NCtype(gbp, &type);
    if (status != NC_NOERR) {
        printf("preamble of \n");
        return status; 
    }

    /* get nelems */
    status = val_get_size_t(gbp, &tmp);
    if (status != NC_NOERR) {
        printf("the length of ");
        return status;
    }
    ncap->ndefined = tmp; /* number of allowable defined variables < 2^32 */

    err_addr = ((size_t)gbp->pos - (size_t)gbp->base) + (gbp->offset - gbp->size) -
                (X_SIZEOF_INT + x_sizeof_NON_NEG); 

    if (ncap->ndefined == 0) {
        /* no dimension defined */
        if (type != ABSENT) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component type, while ABSENT is expected for ");
            return NC_ENOTNC;
        }
    } else {
        if (type != NC_DIMENSION) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component type, while NC_DIMENSION is expected as number of dimensions is %d for ", ncap->ndefined);
            return NC_ENOTNC;
        }

        /* check each dimension */
        ncap->value = (NC_dim **) NCI_Malloc(ncap->ndefined * sizeof(NC_dim *));
        if (ncap->value == NULL) return NC_ENOMEM;
        ncap->nalloc = ncap->ndefined;

        dpp = ncap->value;
        end = &dpp[ncap->ndefined];
        for (/*NADA*/ dim = 0; dpp < end; dpp++, dim++) {
            status = val_get_NC_dim(gbp, dpp);
            if (status != NC_NOERR) {
	        printf("dimension[%d] in ", dim);
                ncap->ndefined = dpp - ncap->value;
                ncmpii_free_NC_dimarray(ncap);
                return status;
            }
        }
    }

    return NC_NOERR;
}

static int
val_get_nc_type(bufferinfo *gbp, nc_type *typep) {
    /* NCtype is 4-byte integer */
    unsigned int type = 0;
    int status = val_check_buffer(gbp, 4);
    if (status != NC_NOERR) return status;

    /* get a 4-byte integer */
    status = ncmpix_get_uint32((const void**)(&gbp->pos), &type);
    gbp->index += X_SIZEOF_INT;
    if (status != NC_NOERR) return status;

  if (   type != NC_BYTE
      && type != NC_UBYTE
      && type != NC_CHAR
      && type != NC_SHORT
      && type != NC_USHORT
      && type != NC_INT
      && type != NC_UINT
      && type != NC_FLOAT
      && type != NC_DOUBLE
      && type != NC_INT64
      && type != NC_UINT64) {
    printf("Error @ [0x%8.8Lx]: \n\tUnknown data type for the values of ",
	   (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - X_SIZEOF_INT));
    return NC_EINVAL; 
  }
 
  *typep = (nc_type) type;

  return NC_NOERR;
}

/*
 * Get the values of an attribute  
 */
static int
val_get_NC_attrV(bufferinfo *gbp, NC_attr *attrp) {
    int status;
    void *value = attrp->xvalue;
    char pad[X_ALIGN-1]; 
    MPI_Offset nvalues = attrp->nelems, esz, padding, bufremain, attcount;
    MPI_Aint pos_addr, base_addr;

    esz = ncmpix_len_nctype(attrp->type);
    padding = attrp->xsz - esz * nvalues;
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    bufremain = gbp->size - (pos_addr - base_addr);

  while (nvalues > 0) {
    if (bufremain > 0) {
      attcount = MIN(bufremain, esz * nvalues);
      (void) memcpy(value, gbp->pos, attcount);
      nvalues -= attcount/esz;
      gbp->pos = (void *)((char *)gbp->pos + attcount);
      value = (void *)((char *)value + attcount);
      bufremain -= attcount;
    } else {
      status = val_fetch(gbp);
      if(status != NC_NOERR) {
	printf("fetching the values of ");
        return status;
      }
      bufremain = gbp->size;
    }
  }
 
  memset(pad, 0, X_ALIGN-1);
  if (memcmp(gbp->pos, pad, padding) != 0) {
    printf("Error @ [0x%8.8Lx]: \n\tPadding should be 0x00 for the values alignment of ",
           (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size)); 
    return NC_EINVAL;
  }
  gbp->pos = (void *)((char *)gbp->pos + padding);

  return NC_NOERR;
}

static int
val_get_NC_attr(bufferinfo *gbp, NC_attr **attrpp) {
  NC_string *strp;
  int status;
  nc_type type; 
  MPI_Offset nelems;
  NC_attr *attrp;

  status = val_get_NC_string(gbp, &strp);
  if(status != NC_NOERR)
    return status;

  status = val_get_nc_type(gbp, &type);
  if(status != NC_NOERR) {
    printf("\"%s\" - ", strp->cp);
    ncmpii_free_NC_string(strp);
    return status;
  }

  status = val_get_size_t(gbp, &nelems); 
  if(status != NC_NOERR) {
    printf("the values of \"%s\" - ", strp->cp);
    ncmpii_free_NC_string(strp);
    return status;
  }

  attrp = ncmpii_new_x_NC_attr(strp, type, nelems);
  if(attrp == NULL) {
    ncmpii_free_NC_string(strp);
    return status;
  }

  status = val_get_NC_attrV(gbp, attrp);
  if(status != NC_NOERR) {
    printf("\"%s\" - ", strp->cp);
    ncmpii_free_NC_attr(attrp); /* frees strp */ 
    return status;
  }

  *attrpp = attrp; 
  
  return NC_NOERR; 
}

static int
val_get_NC_attrarray(bufferinfo *gbp, NC_attrarray *ncap)
{
    /* netCDF file format:
     *  ...
     * att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_ATTRIBUTE = \x00 \x00 \x00 \x0C         // tag for list of attributes
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int status;
    NCtype type = NC_UNSPECIFIED;
    NC_attr **app, **end;
    int att;
    MPI_Offset tmp;
    unsigned long long err_addr;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    status = val_get_NCtype(gbp, &type);
    if (status != NC_NOERR) {
        printf("preamble of ");
        return status; 
    }

    /* get nelems */
    status = val_get_size_t(gbp, &tmp);
    if (status != NC_NOERR) {
        printf("the length of ");
        return status;
    }
    ncap->ndefined = tmp; /* number of allowable defined variables < 2^32 */

    err_addr = ((size_t)gbp->pos - (size_t)gbp->base) + (gbp->offset - gbp->size) -
                (X_SIZEOF_INT + x_sizeof_NON_NEG); 

    if (ncap->ndefined == 0) {
        /* no attribute defined */
        if (type != ABSENT) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component type, while ABSENT is expected for ");
            return NC_ENOTNC;
        }
    } else {
        if (type != NC_ATTRIBUTE) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component type, while NC_ATTRIBUTE is expected as number of dimensions is %d for ", ncap->ndefined);
            return NC_ENOTNC;
        }

        ncap->value = (NC_attr **) NCI_Malloc(ncap->ndefined * sizeof(NC_attr *));
        if (ncap->value == NULL)
            return NC_ENOMEM;
        ncap->nalloc = ncap->ndefined; 

        app = ncap->value;
        end = &app[ncap->ndefined];
        for( /*NADA*/ att = 0; app < end; app++, att++) {
            status = val_get_NC_attr(gbp, app);
            if (status != NC_NOERR) {
	        printf("attribute[%d] of ", att);
                ncap->ndefined = app - ncap->value;
                ncmpii_free_NC_attrarray(ncap);
                return status;
            }
        }
    }
  
    return NC_NOERR;
}

static int
val_get_NC_var(bufferinfo *gbp, NC_var **varpp)
{
    /* netCDF file format:
     * netcdf_file = header data
     * header      = magic numrecs dim_list gatt_list var_list
     *  ...
     * var         = name nelems [dimid ...] vatt_list nc_type vsize begin
     * nelems      = NON_NEG
     * dimid       = NON_NEG
     * vatt_list   = att_list
     * nc_type     = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * vsize       = NON_NEG
     * begin       = OFFSET        // Variable start location.
     * OFFSET      = <non-negative INT> |  // CDF-1
     *               <non-negative INT64>  // CDF-2 and CDF-5
     * NON_NEG     = <non-negative INT> |  // CDF-1 and CDF-2
     *               <non-negative INT64>  // CDF-5
     */
    NC_string *strp;
    int dim, status;
    MPI_Offset ndims, dimid;
    NC_var *varp;

    status = val_get_NC_string(gbp, &strp);
    if (status != NC_NOERR) return status;

    status = val_get_size_t(gbp, &ndims);
    if (status != NC_NOERR) {
        printf("the dimid list of \"%s\" - ", strp->cp);
        ncmpii_free_NC_string(strp); 
        return status;
    }

    varp = ncmpii_new_x_NC_var(strp, ndims);
    if (varp == NULL) {
        ncmpii_free_NC_string(strp);
        return NC_ENOMEM;
    }

    for (dim=0; dim<ndims; dim++) {
        status = val_check_buffer(gbp, (gbp->version < 5 ? 4 : 8));
        if (status != NC_NOERR) {
            printf("the dimid[%d] is expected for \"%s\" - ", dim, strp->cp);
            ncmpii_free_NC_var(varp);
            return status;
        }
        if (gbp->version < 5) {
            unsigned int tmp=0;
            status = ncmpix_get_uint32((const void **)(&gbp->pos), &tmp);
            dimid = (int)tmp;
        }
        else {
            unsigned long long tmp=0;
            status = ncmpix_get_uint64((const void **)(&gbp->pos), &tmp);
            dimid = (int)tmp;
        }
        varp->dimids[dim] = dimid;
        if (status != NC_NOERR) {
            ncmpii_free_NC_var(varp);
            return status;
        }
    }

    status = val_get_NC_attrarray(gbp, &varp->attrs);
    if (status != NC_NOERR) {
        printf("ATTRIBUTE list of \"%s\" - ", strp->cp);
        ncmpii_free_NC_var(varp);
        return status;
    }

    status = val_get_nc_type(gbp, &varp->type);
    if (status != NC_NOERR) {
        printf("\"%s\" - ", strp->cp);
        ncmpii_free_NC_var(varp);
        return status;
    } 

    /* TODO: instead of getting vsize from file, we recalculate it */
    status = val_get_size_t(gbp, &varp->len);
    if (status != NC_NOERR) {
        printf("the data of  \"%s\" - ", strp->cp);
        ncmpii_free_NC_var(varp);
        return status;
    }

    status = val_check_buffer(gbp, (gbp->version < 5 ? 4 : 8));
    if (status != NC_NOERR) {
        printf("offset is expected for the data of \"%s\" - ", strp->cp);
        ncmpii_free_NC_var(varp);
        return status;
    }
    if (gbp->version == 1) {
        unsigned int tmp=0;
        status = ncmpix_get_uint32((const void **)(&gbp->pos), &tmp);
        varp->begin = (MPI_Offset)tmp;
    }
    else {
        unsigned long long tmp=0;
        status = ncmpix_get_uint64((const void **)(&gbp->pos), &tmp);
        varp->begin = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) {
        ncmpii_free_NC_var(varp);
        return status;
    }

    *varpp = varp;
    return NC_NOERR;
}

static int
val_get_NC_vararray(bufferinfo *gbp, NC_vararray *ncap)
{
    /* netCDF file format:
     * netcdf_file = header  data
     * header      = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var_list    = ABSENT | NC_VARIABLE   nelems  [var ...]
     * ABSENT      = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *               ZERO  ZERO64  // for CDF-5
     * ZERO        = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64      = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_VARIABLE = \x00 \x00 \x00 \x0B         // tag for list of variables
     * nelems      = NON_NEG       // number of elements in following sequence
     * NON_NEG     = <non-negative INT> |        // CDF-1 and CDF-2
     *               <non-negative INT64>        // CDF-5
     */
    int status;
    NCtype type = NC_UNSPECIFIED;
    NC_var **vpp, **end;
    int var;
    MPI_Offset tmp;
    unsigned long long err_addr;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL); 

    status = val_get_NCtype(gbp, &type);
    if (status != NC_NOERR) {
        printf("preamble of ");
        return status;
    }
 
    status = val_get_size_t(gbp, &tmp);
    if(status != NC_NOERR) {
        printf("the length of ");
        return status;
    }
    ncap->ndefined = tmp; /* number of allowable defined variables < 2^32 */
 
    err_addr = ((size_t)gbp->pos - (size_t)gbp->base) + (gbp->offset - gbp->size) -
                (X_SIZEOF_INT + x_sizeof_NON_NEG);

    if(ncap->ndefined == 0) {
        if (type != ABSENT) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component type, while ABSENT is expected for ");
            return NC_ENOTNC;
        }
    } else {
        if (type != NC_VARIABLE) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component type, while NC_VARIABLE is expected as number of dimensions is %d for ", ncap->ndefined);
            return NC_ENOTNC;
        }
 
        ncap->value = (NC_var **) NCI_Malloc(ncap->ndefined * sizeof(NC_var *));
        if (ncap->value == NULL) return NC_ENOMEM; 
        ncap->nalloc = ncap->ndefined;

        vpp = ncap->value;
        end = &vpp[ncap->ndefined];
        for (/*NADA*/ var = 0; vpp < end; vpp++, var++) {
            status = val_get_NC_var(gbp, vpp);
            if (status != NC_NOERR) {
                printf("variable[%d] in ", var);
                ncap->ndefined = vpp - ncap->value;
                ncmpii_free_NC_vararray(ncap);
                return status;
            }
        }
    }

    return NC_NOERR;
}

/*
 * Given a valid ncp, check all variables for their sizes against the maximal
 * allowable sizes. Different CDF formation versions have different maximal
 * sizes. This function returns NC_EVARSIZE if any variable has a bad len
 * (product of non-rec dim sizes too large), else return NC_NOERR.
 */
static int
val_NC_check_vlens(NC *ncp)
{
    NC_var **vpp;
    /* maximum permitted variable size (or size of one record's worth
       of a record variable) in bytes.  This is different for format 1
       and format 2. */
    MPI_Offset ii, vlen_max, rec_vars_count;
    MPI_Offset large_fix_vars_count, large_rec_vars_count;
    int last = 0;

    if (ncp->vars.ndefined == 0)
        return NC_NOERR;

    if (ncp->format >= 5) /* CDF-5 */
        return NC_NOERR;

    /* only CDF-1 and CDF-2 need to continue */

    if (ncp->flags & NC_64BIT_OFFSET) /* CDF2 format */
        vlen_max = X_UINT_MAX - 3; /* "- 3" handles rounded-up size */
    else
        vlen_max = X_INT_MAX - 3; /* CDF1 format */

    /* Loop through vars, first pass is for non-record variables */
    large_fix_vars_count = 0;
    rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.ndefined; ii++, vpp++) {
        if (!IS_RECVAR(*vpp)) {
            last = 0;
            if (ncmpii_NC_check_vlen(*vpp, vlen_max) == 0) {
                /* check this variable's shape product against vlen_max */
                large_fix_vars_count++;
                last = 1;
            }
        } else {
            rec_vars_count++;
        }
    }
    /* OK if last non-record variable size too large, since not used to
       compute an offset */
    if (large_fix_vars_count > 1) {  /* only one "too-large" variable allowed */
        printf("CDF-%d format allows only one large fixed-size variable\n",ncp->format);
        return NC_EVARSIZE;
    }

    /* The only "too-large" variable must be the last one defined */
    if (large_fix_vars_count == 1 && last == 0) {
        printf("CDF-%d format allows only one large fixed-size variable and it must be the last one defined\n",ncp->format);
        return NC_EVARSIZE;
    }

    if (rec_vars_count == 0) return NC_NOERR;

    /* if there is a "too-large" fixed-size variable, no record variable is
     * allowed */
    if (large_fix_vars_count == 1) {
        printf("CDF-%d format allows only one large fixed-size variable when there is no record variable defined\n",ncp->format);
        return NC_EVARSIZE;
    }

    /* Loop through vars, second pass is for record variables.   */
    large_rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.ndefined; ii++, vpp++) {
        if (IS_RECVAR(*vpp)) {
            last = 0;
            if (ncmpii_NC_check_vlen(*vpp, vlen_max) == 0) {
                /* check this variable's shape product against vlen_max */
                large_rec_vars_count++;
                last = 1;
            }
        }
    }

    /* For CDF-2, no record variable can require more than 2^32 - 4 bytes of
     * storage for each record's worth of data, unless it is the last record
     * variable. See
     * http://www.unidata.ucar.edu/software/netcdf/docs/file_structure_and_performance.html#offset_format_limitations
     */
    if (large_rec_vars_count > 1) { /* only one "too-large" variable allowed */
        printf("CDF-%d format allows only one large record variable\n",ncp->format);
        return NC_EVARSIZE;
    }

    /* and it has to be the last one */
    if (large_rec_vars_count == 1 && last == 0) {
        printf("CDF-%d format allows only one large record variable and it must be the last one defined\n",ncp->format);
        return NC_EVARSIZE;
    }

    return NC_NOERR;
}

static int
val_get_NC(NC *ncp)
{
    int status;
    bufferinfo getbuf;
    char magic[sizeof(ncmagic)];
    MPI_Offset nrecs = 0;
    MPI_Aint pos_addr, base_addr;

    /* Initialize the get buffer that stores the header read from the file */
    getbuf.nciop = ncp->nciop;
    getbuf.offset = 0;     /* read from start of the file */

    /* CDF-5's minimum header size is 4 bytes more than CDF-1 and CDF-2's */
    getbuf.size = _RNDUP( MAX(MIN_NC_XSZ+4, ncp->chunk), X_ALIGN );
    if (getbuf.size > NC_DEFAULT_CHUNKSIZE)
        getbuf.size = NC_DEFAULT_CHUNKSIZE;

    getbuf.pos = getbuf.base = (void *)NCI_Malloc(getbuf.size);
    getbuf.index = 0;

    /* Fetch the next header chunk. The chunk is 'gbp->size' bytes big
     * netcdf_file = header data
     * header      = magic numrecs dim_list gatt_list var_list
     */
    status = val_fetch(&getbuf);
    if (status != NC_NOERR) goto fn_exit;
  
    /* First get the file format information, magic */
    memset(magic, 0, sizeof(magic));
    status = ncmpix_getn_text((const void **)(&getbuf.pos), sizeof(magic), magic);
    getbuf.index += sizeof(magic);

    if (memcmp(magic, ncmagic, sizeof(ncmagic)-1) != 0) {
        printf("Error: Unknow file signature, (C D F \\001, \\002, or \\005) is expected!\n");
        status = NC_ENOTNC;
        goto fn_exit;
    }

    /* check version number in last byte of magic */
    if (magic[sizeof(ncmagic)-1] == 0x1) {
        getbuf.version = 1;
        ncp->format = 1;
        fSet(ncp->flags, NC_32BIT);
    } else if (magic[sizeof(ncmagic)-1] == 0x2) {
        getbuf.version = 2;
        ncp->format = 2;
        fSet(ncp->flags, NC_64BIT_OFFSET);
        if (sizeof(MPI_Offset) != 8) {
            /* take the easy way out: if we can't support all CDF-2
             * files, return immediately */
            status = NC_ESMALL;
            goto fn_exit;
        }
    } else if (magic[sizeof(ncmagic)-1] == 0x5) {
        getbuf.version = 5;
        ncp->format = 5;
        fSet(ncp->flags, NC_64BIT_DATA);
        if (sizeof(MPI_Offset) != 8) {
            status = NC_ESMALL;
            goto fn_exit;
        }
    } else {
        status = NC_ENOTNC;
        goto fn_exit;
    }

    /* header = magic numrecs dim_list gatt_list var_list
     * Check numrecs
     */
    status = val_check_buffer(&getbuf, (getbuf.version < 5) ? 4 : 8);
    if (status != NC_NOERR) {
        printf("Error: number of records is expected!\n");
        status = NC_ENOTNC;
        goto fn_exit;
    }

    /* get numrecs from getbuf into ncp */
    if (getbuf.version < 5) {
        unsigned int tmp=0;
        status = ncmpix_get_uint32((const void **)(&getbuf.pos), &tmp);
        nrecs = (MPI_Offset)tmp;
    }
    else {
        unsigned long long tmp=0;
        status = ncmpix_get_uint64((const void **)(&getbuf.pos), &tmp);
        nrecs = (MPI_Offset)tmp;
    }
    if (status != NC_NOERR) goto fn_exit;

    if (getbuf.version < 5)
        x_sizeof_NON_NEG = 4;
    else
        x_sizeof_NON_NEG = 8;

    /* advance buffer for 4 or 8 bytes */
    getbuf.index += x_sizeof_NON_NEG;
    ncp->numrecs = nrecs;

#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(getbuf.pos,  &pos_addr);
    MPI_Get_address(getbuf.base, &base_addr);
#else
    MPI_Address(getbuf.pos,  &pos_addr);
    MPI_Address(getbuf.base, &base_addr);
#endif
    assert(pos_addr < base_addr + getbuf.size);

    /* header = magic numrecs dim_list gatt_list var_list
     * dim_list = ABSENT | NC_DIMENSION  nelems  [dim ...]
     * Check dim_list
     */
    status = val_get_NC_dimarray(&getbuf, &ncp->dims);
    if (status != NC_NOERR) {
        printf("DIMENSION list!\n");
        goto fn_exit;
    }

    /* header = magic numrecs dim_list gatt_list var_list
     * att_list = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     * Check att_list
     */
    status = val_get_NC_attrarray(&getbuf, &ncp->attrs); 
    if (status != NC_NOERR) {
        printf("GLOBAL ATTRIBUTE list!\n");
        goto fn_exit;
    }

    /* header = magic numrecs dim_list gatt_list var_list
     * var_list    = ABSENT | NC_VARIABLE   nelems  [var ...]
     * Check var_list
     */
    status = val_get_NC_vararray(&getbuf, &ncp->vars);
    if (status != NC_NOERR) {
        printf("VARIABLE list!\n");
        goto fn_exit;
    }

    ncp->xsz = ncmpii_hdr_len_NC(ncp);
    status = ncmpii_NC_computeshapes(ncp);
    if (status != NC_NOERR) goto fn_exit;

    status = val_NC_check_vlens(ncp);
    if (status != NC_NOERR) goto fn_exit;

fn_exit:
    NCI_Free(getbuf.base);

    return status;
}

/* End Of get NC */

int main(int argc, char **argv)
{
    char *filename;
    int fd, status=NC_NOERR;
    NC *ncp=NULL;
    struct stat ncfilestat;

    MPI_Init(&argc, &argv);

    if (argc != 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        MPI_Finalize();
        return 1;
    } 
    filename = argv[1];

    fd = open(filename, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Error on open file %s: %s\n", filename,strerror(errno));
        MPI_Finalize();
        return 1;
    }

    /* Allocate NC object */
    ncp = ncmpii_new_NC(NULL);
    if (ncp == NULL) {
        status = NC_ENOMEM;
        ABORT(status)
        goto prog_exit;
    }

    ncp->nciop = ncmpiio_new(filename, NC_NOWRITE);
    if (ncp->nciop == NULL) {
        status = NC_ENOMEM;
        ABORT(status)
        goto prog_exit;
    }

    ncp->nciop->fd = fd;

    /* read and validate the header */
    status = val_get_NC(ncp);
    if (status != NC_NOERR && status != -1) goto prog_exit;

    /* check data size */
    if (-1 == fstat(ncp->nciop->fd, &ncfilestat)) {
        printf("Error at line %d fstat (%s)\n",__LINE__,strerror(errno));
        status = NC_EFILE;
        goto prog_exit;
    }
    if ( ncp->begin_rec + ncp->recsize * ncp->numrecs < ncfilestat.st_size ) {
        printf("Error: file size is larger than defined!\n");
        status = NC_EFILE;
        goto prog_exit;
    }
    else if ( ncp->numrecs > 0 &&
              ncp->begin_rec + ncp->recsize * (ncp->numrecs - 1) > ncfilestat.st_size ) {
        printf("Error: file size is less than expected!\n");
        printf("\tbegin_rec=%lld recsize=%lld numrecs=%lld ncfilestat.st_size=%lld\n",ncp->begin_rec, ncp->recsize, ncp->numrecs, (long long) ncfilestat.st_size);
        status = NC_EFILE;
        goto prog_exit;
    }

prog_exit:
    /* free allocated buffers and close the file */
    if (ncp != NULL) {
        if (ncp->nciop != NULL) ncmpiio_free(ncp->nciop);
        ncmpii_free_NC(ncp);
    }
    close(fd);

    if (status == NC_NOERR)
        printf("The netCDF file is validated!\n");
    else
        printf("Error: %s\n",ncmpi_strerror(status));
    MPI_Finalize();
    return 0;
}
