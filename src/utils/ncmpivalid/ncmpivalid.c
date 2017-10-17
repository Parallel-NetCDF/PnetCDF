/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>  /* open() */
#include <sys/stat.h>   /* open() */
#include <fcntl.h>      /* open() */
#include <unistd.h>     /* read() getopt() */
#include <string.h>     /* strcpy(), strncpy() */
#include <assert.h>
#include <errno.h>

#include <mpi.h>

/* TODO: should not use any PnetCDF source codes, as this CDF format validate
 * utility should run independently from PnetCDF
 */
#include <ncmpio_NC.h>
#include <ncx.h>
#include <common.h>

#ifndef EXIT_FAILURE
#ifndef vms
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#else
/* In OpenVMS, success is indicated by odd values and failure by even values. */
#define EXIT_SUCCESS 1
#define EXIT_FAILURE 0
#endif
#endif

#define DEBUG
#ifdef DEBUG
#define DEBUG_RETURN(e) {                                             \
    printf("(Error %s at line %d in file %s)\n",#e,__LINE__,__FILE__);\
    return e;                                                         \
}
#else
#define DEBUG_RETURN(e) return e;
#endif

/*
 * "magic number" at beginning of file: 0x43444601 (big Edian) 
 */
static const schar ncmagic[] = {'C', 'D', 'F', 0x01}; 

#define ABSENT 0
#define X_SIZEOF_INT 4
static int x_sizeof_NON_NEG;

/*----< ncmpio_xlen_nc_type() >----------------------------------------------*/
/* return the length of external NC data type */
static int
xlen_nc_type(nc_type xtype) {
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  return 1;
        case NC_SHORT:
        case NC_USHORT: return 2;
        case NC_INT:
        case NC_UINT:
        case NC_FLOAT:  return 4;
        case NC_DOUBLE:
        case NC_INT64:
        case NC_UINT64: return 8;
        default: DEBUG_RETURN(NC_EBADTYPE)
    }
}

/* calculate the followings
 *   ncp->begin_var           first variable's offset, file header extent
 *   ncp->begin_rec           first record variable's offset
 *   ncp->recsize             sum of all single record size of all variables
 *   ncp->vars.value[*]->len  individual variable size (record size)
 */
static int
compute_var_shape(NC *ncp)
{
    int i, j, err;
    NC_var *first_var = NULL;       /* first "non-record" var */
    NC_var *first_rec = NULL;       /* first "record" var */

    if (ncp->vars.ndefined == 0) return NC_NOERR;

    ncp->begin_var = ncp->xsz;
    ncp->begin_rec = ncp->xsz;
    ncp->recsize   = 0;

    for (i=0; i<ncp->vars.ndefined; i++) {
        /* check if dimids are valid */
        for (j=0; j<ncp->vars.value[i]->ndims; j++) {
            if (ncp->vars.value[i]->dimids[j] < 0 ||
                ncp->vars.value[i]->dimids[j] >= ncp->dims.ndefined)
                DEBUG_RETURN(NC_EBADDIM) /* dimid is not defined */
        }
        /* ncp->vars.value[i]->len will be recomputed from dimensions in
         * ncmpio_NC_var_shape64() */
        err = ncmpio_NC_var_shape64(ncp->vars.value[i], &ncp->dims);
        if (err != NC_NOERR) return err;

        if (IS_RECVAR(ncp->vars.value[i])) {
            if (first_rec == NULL) first_rec = ncp->vars.value[i];
            ncp->recsize += ncp->vars.value[i]->len;
        }
        else { /* fixed-size variable */
            if (first_var == NULL) first_var = ncp->vars.value[i];
            ncp->begin_rec = ncp->vars.value[i]->begin
                           + ncp->vars.value[i]->len;
        }
    }

    if (first_rec != NULL) {
        if (ncp->begin_rec > first_rec->begin)
            DEBUG_RETURN(NC_ENOTNC) /* not a netCDF file or corrupted */

        ncp->begin_rec = first_rec->begin;
        /*
         * for special case of exactly one record variable, pack value
         */
        if (ncp->recsize == first_rec->len)
            ncp->recsize = *first_rec->dsizes * first_rec->xsz;
    }

    if (first_var != NULL)
        ncp->begin_var = first_var->begin;
    else
        ncp->begin_var = ncp->begin_rec;

    if (ncp->begin_var <= 0 || ncp->xsz > ncp->begin_var ||
        ncp->begin_rec <= 0 || ncp->begin_var > ncp->begin_rec)
        DEBUG_RETURN(NC_ENOTNC) /* not a netCDF file or corrupted */

    return NC_NOERR;
}

/*
 * Fetch the next header chunk.
 */
static int
val_fetch(int fd, bufferinfo *gbp) {
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

    if (-1 == lseek(fd, gbp->offset-slack, SEEK_SET)) {
        printf("Error at line %d: lseek %s\n",__LINE__,strerror(errno));
        return -1;
    }
    nn = read(fd, gbp->base, gbp->size);
    if (nn == -1) {
        printf("Error at line %d: lseek %s\n",__LINE__,strerror(errno));
        return -1;
    }
/*
    if (nn < gbp->size) {
        printf("Error: file header size is less than expected\n");
printf("Error: pos_addr=%ld base_addr=%ld gbp->size=%lld nn=%zd\n",pos_addr,base_addr,gbp->size,nn);
        DEBUG_RETURN(NC_ENOTNC)
    }
*/
    gbp->offset += (gbp->size - slack);

    return NC_NOERR;
}

/*
 * Ensure that 'nextread' bytes are available.
 */
static int
val_check_buffer(int         fd,
                 bufferinfo *gbp,
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

    return val_fetch(fd, gbp);
} 

static int
val_get_NC_tag(int fd, bufferinfo *gbp, NC_tag *tagp)
{
    unsigned int tag = 0;
    int status = val_check_buffer(fd, gbp, x_sizeof_NON_NEG);
    if (status != NC_NOERR) {
        printf("NC component tag is expected for ");
        return status;
    }

    status = ncmpix_get_uint32((const void**)(&gbp->pos), &tag);
    if (status != NC_NOERR) return status;
    *tagp = (NC_tag) tag;
    return NC_NOERR;
}

static int
val_get_size_t(int fd, bufferinfo *gbp, MPI_Offset *sp) {
  int sizeof_t = (gbp->version < 5) ? 4 : 8; 
  int status = val_check_buffer(fd, gbp, sizeof_t);
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
val_get_NC_string(int fd, bufferinfo *gbp, char **namep) {
    int status;
    char *cpos, pad[X_ALIGN-1];
    MPI_Offset nchars=0, padding, bufremain, strcount;
    MPI_Aint pos_addr, base_addr;

    *namep = NULL;
    status = val_get_size_t(fd, gbp, &nchars);
    if (status != NC_NOERR) {
        printf("the name string of ");
        return status;
    }

    *namep = (char*) NCI_Malloc((size_t)nchars + 1);
    if (*namep == NULL) DEBUG_RETURN(NC_ENOMEM)
    (*namep)[nchars] = '\0'; /* add terminal character */

    padding = _RNDUP(nchars, X_ALIGN) - nchars;
#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(gbp->pos,  &pos_addr);
    MPI_Get_address(gbp->base, &base_addr);
#else
    MPI_Address(gbp->pos,  &pos_addr);
    MPI_Address(gbp->base, &base_addr);
#endif
    bufremain = gbp->size - (pos_addr - base_addr);
    cpos = *namep;

    while (nchars > 0) {
        if (bufremain > 0) {
            strcount = MIN(bufremain, nchars); 
            (void) memcpy(cpos, gbp->pos, strcount);
            nchars -= strcount;
            gbp->pos = (void *)((char *)gbp->pos + strcount);
            cpos += strcount; 
            bufremain -= strcount;
        } else {
            status = val_fetch(fd, gbp);
            if (status != NC_NOERR) {
                printf("fetching the name string of ");
                NCI_Free(*namep);
                *namep = NULL;
                return status;
            } 
            bufremain = gbp->size;
        }
    }

    if (padding > 0) {
        memset(pad, 0, X_ALIGN-1);
        status = val_check_buffer(fd, gbp, padding);
        if (status != NC_NOERR) {
            printf("fetching padding for the name string of ");
            return status;
        } 
        if (memcmp(gbp->pos, pad, padding) != 0) {
            printf("Error @ [0x%8.8Lx]: \n\tPadding should be 0x00 for the name string alignment of ", (long long unsigned)
	           (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size));
            NCI_Free(*namep);
            *namep = NULL;
            DEBUG_RETURN(NC_ENOTNC)
        }
        gbp->pos = (void *)((char *)gbp->pos + padding);
    }

    return NC_NOERR;  
}

static int
val_get_NC_dim(int fd, bufferinfo *gbp, NC_dim **dimpp) {
    int status;
    char *name=NULL;
    NC_dim *dimp;

    *dimpp = NULL;

    status = val_get_NC_string(fd, gbp, &name);
    if (status != NC_NOERR) {
        if (name != NULL) NCI_Free(name);
        return status;
    }

    dimp = (NC_dim*) NCI_Malloc(sizeof(NC_dim));
    if (dimp == NULL) {
        if (name != NULL) NCI_Free(name);
        DEBUG_RETURN(NC_ENOMEM)
    }
    dimp->name     = name;
    dimp->name_len = strlen(name);

    status = val_get_size_t(fd, gbp, &dimp->size);
    if (status != NC_NOERR) { /* frees dimp */
        printf("\"%s\" - ", name);
        NCI_Free(dimp->name);
        NCI_Free(dimp);
        return status;
    }

    *dimpp = dimp;

    return NC_NOERR;
}

static int
val_get_NC_dimarray(int fd, bufferinfo *gbp, NC_dimarray *ncap)
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
    NC_tag tag = NC_UNSPECIFIED; 
    int dim;
    unsigned long long err_addr;
    MPI_Offset tmp;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* read NC_tag (NC_DIMENSION or ZERO) from gbp buffer */
    status = val_get_NC_tag(fd, gbp, &tag);
    if (status != NC_NOERR) {
        printf("preamble of \n");
        return status; 
    }

    /* read nelems (number of dimensions) from gbp buffer */
    status = val_get_size_t(fd, gbp, &tmp);
    if (status != NC_NOERR) {
        printf("the length of ");
        return status;
    }
    if (tmp > NC_MAX_DIMS) {
        /* number of allowable defined dimensions NC_MAX_DIMS */
        printf("the length of ");
        return NC_EMAXDIMS;
    }
    ncap->ndefined = (int)tmp;

    err_addr = (size_t)gbp->pos - (size_t)gbp->base + gbp->offset - gbp->size -
               (X_SIZEOF_INT + x_sizeof_NON_NEG); 

    if (ncap->ndefined == 0) {
        /* no dimension defined */
        /* From the CDF file format specification, the tag is either
         * NC_DIMENSION or ABSENT (ZERO), but we follow NetCDF library to skip
         * checking the tag when ndefined is zero.
         */
        return NC_NOERR;
#if 0
        if (tag != ABSENT) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component tag, while ABSENT is expected for ");
            DEBUG_RETURN(NC_ENOTNC)
        }
#endif
    } else {
        if (tag != NC_DIMENSION) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component tag, while NC_DIMENSION is expected as number of dimensions is %d for ", ncap->ndefined);
            DEBUG_RETURN(NC_ENOTNC)
        }

        /* check each dimension */
        size_t alloc_size = (size_t)ncap->ndefined + NC_ARRAY_GROWBY;
        ncap->value = (NC_dim **) NCI_Calloc(alloc_size, sizeof(NC_dim *));
        if (ncap->value == NULL) DEBUG_RETURN(NC_ENOMEM)

        for (dim=0; dim<ncap->ndefined; dim++) {
            status = val_get_NC_dim(fd, gbp, &ncap->value[dim]);
            if (status != NC_NOERR) {
	        printf("dimension[%d] in ", dim);
                ncap->ndefined = dim;
                ncmpio_free_NC_dimarray(ncap);
                return status;
            }
        }
    }

    return NC_NOERR;
}

static int
val_get_nc_type(int fd, bufferinfo *gbp, nc_type *xtypep) {
    /* nc_type is 4-byte integer */
    unsigned int xtype = 0;
    int status = val_check_buffer(fd, gbp, 4);
    if (status != NC_NOERR) return status;

    /* get a 4-byte integer */
    status = ncmpix_get_uint32((const void**)(&gbp->pos), &xtype);
    if (status != NC_NOERR) return status;

    if (xtype < NC_BYTE) goto err_exit;

    if (gbp->version < 5) {
        if (xtype > NC_DOUBLE) goto err_exit;
    }
    else if (xtype > NC_UINT64) goto err_exit;

    *xtypep = (nc_type) xtype;

    return NC_NOERR;

err_exit:
    printf("Error @ [0x%8.8Lx]: \n\tUnknown NC data xtype for the values of ",
	   (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - X_SIZEOF_INT));
    DEBUG_RETURN(NC_EBADTYPE)
}

/*
 * Get the values of an attribute  
 */
static int
val_get_NC_attrV(int fd, bufferinfo *gbp, NC_attr *attrp) {
    int status;
    void *value = attrp->xvalue;
    char pad[X_ALIGN-1]; 
    MPI_Offset nvalues, padding, bufremain, attcount;
    MPI_Aint pos_addr, base_addr;

    nvalues = attrp->nelems * xlen_nc_type(attrp->xtype);
    padding = attrp->xsz - nvalues;
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
            attcount = MIN(bufremain, nvalues);
            (void) memcpy(value, gbp->pos, attcount);
            nvalues -= attcount;
            gbp->pos = (void *)((char *)gbp->pos + attcount);
            value = (void *)((char *)value + attcount);
            bufremain -= attcount;
        } else {
            status = val_fetch(fd, gbp);
            if(status != NC_NOERR) {
	        printf("fetching the values of ");
                return status;
            }
            bufremain = gbp->size;
        }
    }
 
    if (padding > 0) {
        memset(pad, 0, X_ALIGN-1);
        if (memcmp(gbp->pos, pad, padding) != 0) {
            printf("Error @ [0x%8.8Lx]: \n\tPadding should be 0x00 for the values alignment of ",
                   (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size)); 
            DEBUG_RETURN(NC_ENOTNC)
        }
        gbp->pos = (void *)((char *)gbp->pos + padding);
    }

    return NC_NOERR;
}

static MPI_Offset
x_len_NC_attrV(nc_type    xtype,
               MPI_Offset nelems)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  return _RNDUP(nelems, 4);
        case NC_SHORT:
        case NC_USHORT: return ((nelems + (nelems)%2) * 2);
        case NC_INT:    return (nelems * 4);
        case NC_UINT:   return (nelems * 4);
        case NC_FLOAT:  return (nelems * 4);
        case NC_DOUBLE: return (nelems * 8);
        case NC_INT64:  return (nelems * 8);
        case NC_UINT64: return (nelems * 8);
        default: fprintf(stderr, "Error: bad xtype(%d) in %s\n",xtype,__func__);
    }
    return 0;
}

static int
new_NC_attr(char        *name,
            nc_type      xtype,
            MPI_Offset   nelems,
            NC_attr    **attrp)
{
    *attrp = (NC_attr*) NCI_Malloc(sizeof(NC_attr));
    if (*attrp == NULL ) DEBUG_RETURN(NC_ENOMEM)

    (*attrp)->xtype    = xtype;
    (*attrp)->xsz      = 0;
    (*attrp)->nelems   = nelems;
    (*attrp)->xvalue   = NULL;
    (*attrp)->name     = name;
    (*attrp)->name_len = strlen(name);

    if (nelems > 0) {
        MPI_Offset xsz = x_len_NC_attrV(xtype, nelems);
        (*attrp)->xsz    = xsz;
        (*attrp)->xvalue = NCI_Malloc((size_t)xsz);
        if ((*attrp)->xvalue == NULL) {
            NCI_Free(*attrp);
            *attrp = NULL;
            DEBUG_RETURN(NC_ENOMEM)
        }
    }
    return NC_NOERR;
}

static int
val_get_NC_attr(int fd, bufferinfo *gbp, NC_attr **attrpp) {
  char *name=NULL;
  int status;
  nc_type xtype; 
  MPI_Offset nelems;
  NC_attr *attrp;

  status = val_get_NC_string(fd, gbp, &name);
  if (status != NC_NOERR) {
      if (name != NULL) NCI_Free(name);
      return status;
  }

  status = val_get_nc_type(fd, gbp, &xtype);
  if(status != NC_NOERR) {
    printf("\"%s\" - ", name);
    if (name != NULL) NCI_Free(name);
    return status;
  }

  status = val_get_size_t(fd, gbp, &nelems); 
  if(status != NC_NOERR) {
    printf("the values of \"%s\" - ", name);
    if (name != NULL) NCI_Free(name);
    return status;
  }

  status = new_NC_attr(name, xtype, nelems, &attrp);
  if(status != NC_NOERR) {
    if (name != NULL) NCI_Free(name);
    return status;
  }

  status = val_get_NC_attrV(fd, gbp, attrp);
  if(status != NC_NOERR) {
    printf("\"%s\" - ", name);
    NCI_Free(attrp->name);
    NCI_Free(attrp->xvalue);
    NCI_Free(attrp);
    return status;
  }

  *attrpp = attrp; 
  
  return NC_NOERR; 
}

static int
val_get_NC_attrarray(int fd, bufferinfo *gbp, NC_attrarray *ncap)
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
    NC_tag tag = NC_UNSPECIFIED;
    int att;
    MPI_Offset tmp;
    unsigned long long err_addr;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL);

    /* read NC_tag (NC_ATTRIBUTE or ZERO) from gbp buffer */
    status = val_get_NC_tag(fd, gbp, &tag);
    if (status != NC_NOERR) {
        printf("preamble of ");
        return status; 
    }

    /* read nelems (number of attributes) from gbp buffer */
    status = val_get_size_t(fd, gbp, &tmp);
    if (status != NC_NOERR) {
        printf("the length of ");
        return status;
    }
    if (tmp > NC_MAX_ATTRS) {
        /* number of allowable defined attributes NC_MAX_ATTRS */
        printf("the length of ");
        return NC_EMAXATTS;
    }
    ncap->ndefined = (int)tmp;

    err_addr = (size_t)gbp->pos - (size_t)gbp->base + gbp->offset - gbp->size -
               (X_SIZEOF_INT + x_sizeof_NON_NEG); 

    if (ncap->ndefined == 0) {
        /* no attribute defined */
        /* From the CDF file format specification, the tag is either
         * NC_ATTRIBUTE or ABSENT (ZERO), but we follow NetCDF library to skip
         * checking the tag when ndefined is zero.
         */
        return NC_NOERR;
#if 0
        if (tag != ABSENT) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component tag, while ABSENT is expected for ");
            DEBUG_RETURN(NC_ENOTNC)
        }
#endif
    } else {
        if (tag != NC_ATTRIBUTE) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component tag, while NC_ATTRIBUTE is expected as number of dimensions is %d for ", ncap->ndefined);
            DEBUG_RETURN(NC_ENOTNC)
        }

        size_t alloc_size = (size_t)ncap->ndefined + NC_ARRAY_GROWBY;
        ncap->value = (NC_attr **) NCI_Calloc(alloc_size, sizeof(NC_attr *));
        if (ncap->value == NULL) DEBUG_RETURN(NC_ENOMEM)

        for (att=0; att<ncap->ndefined; att++) {
            status = val_get_NC_attr(fd, gbp, &ncap->value[att]);
            if (status != NC_NOERR) {
	        printf("attribute[%d] of ", att);
                ncap->ndefined = att;
                ncmpio_free_NC_attrarray(ncap);
                return status;
            }
        }
    }
  
    return NC_NOERR;
}

/*----< ncmpio_new_NC_var() >------------------------------------------------*/
static NC_var *
val_new_NC_var(char *name, int ndims)
{
    NC_var *varp;

    varp = (NC_var *) NCI_Calloc(1, sizeof(NC_var));
    if (varp == NULL) return NULL;

    if (ndims > 0) {
        varp->shape  = (MPI_Offset*)NCI_Calloc(ndims, SIZEOF_MPI_OFFSET);
        varp->dsizes = (MPI_Offset*)NCI_Calloc(ndims, SIZEOF_MPI_OFFSET);
        varp->dimids = (int *)      NCI_Calloc(ndims, SIZEOF_INT);
    }

    varp->name     = name;
    varp->name_len = strlen(name);
    varp->ndims    = ndims;
    varp->xsz      = 0;
    varp->len      = 0;
    varp->begin    = 0;

    return varp;
}

static int
val_get_NC_var(int          fd,
               bufferinfo  *gbp,
               NC_var     **varpp,
               int          f_ndims) /* no. dimensions defined in file */

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
    char *name=NULL;
    int dim, dimid, status;
    MPI_Offset ndims;
    NC_var *varp;

    status = val_get_NC_string(fd, gbp, &name);
    if (status != NC_NOERR) {
        if (name != NULL) NCI_Free(name);
        return status;
    }

    status = val_get_size_t(fd, gbp, &ndims);
    if (status != NC_NOERR) {
        printf("the dimid list of \"%s\" - ", name);
        if (name != NULL) NCI_Free(name);
        return status;
    }
    /* cannot be more than NC_MAX_VAR_DIMS */
    if (ndims > NC_MAX_VAR_DIMS) {
        if (name != NULL) NCI_Free(name);
        DEBUG_RETURN(NC_EMAXDIMS)
    }

    varp = val_new_NC_var(name, ndims);
    if (varp == NULL) {
        if (name != NULL) NCI_Free(name);
        DEBUG_RETURN(NC_ENOMEM)
    }

    for (dim=0; dim<ndims; dim++) {
        status = val_check_buffer(fd, gbp, (gbp->version < 5 ? 4 : 8));
        if (status != NC_NOERR) {
            printf("the dimid[%d] is expected for \"%s\" - ", dim, name);
            ncmpio_free_NC_var(varp);
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
        if (status != NC_NOERR) {
            ncmpio_free_NC_var(varp);
            return status;
        }
        /* dimid should be < f_ndims (num of dimensions defined in file) */
        if (dimid >= f_ndims) {
            ncmpio_free_NC_var(varp);
            DEBUG_RETURN(NC_EBADDIM)
        }
        varp->dimids[dim] = dimid;
    }

    /* var = name nelems [dimid ...] vatt_list nc_type vsize begin
     *                               ^^^^^^^^^                     */
    status = val_get_NC_attrarray(fd, gbp, &varp->attrs);
    if (status != NC_NOERR) {
        printf("ATTRIBUTE list of \"%s\" - ", name);
        ncmpio_free_NC_var(varp);
        return status;
    }

    /* var = name nelems [dimid ...] vatt_list nc_type vsize begin
     *                                         ^^^^^^^             */
    status = val_get_nc_type(fd, gbp, &varp->xtype);
    if (status != NC_NOERR) {
        printf("\"%s\" - ", name);
        ncmpio_free_NC_var(varp);
        return status;
    } 

    status = ncmpii_xlen_nc_type(varp->xtype, &varp->xsz);
    if (status != NC_NOERR) {
        printf("\"%s\" - ", name);
        ncmpio_free_NC_var(varp);
        return status;
    } 

    /* var = name nelems [dimid ...] vatt_list nc_type vsize begin
     *                                                 ^^^^^      
     * instead of use vsize from file, we recalculate it in
     * compute_var_shape() */
    status = val_get_size_t(fd, gbp, &varp->len);
    if (status != NC_NOERR) {
        printf("the data of  \"%s\" - ", name);
        ncmpio_free_NC_var(varp);
        return status;
    }

    status = val_check_buffer(fd, gbp, (gbp->version == 1 ? 4 : 8));
    if (status != NC_NOERR) {
        printf("offset is expected for the data of \"%s\" - ", name);
        ncmpio_free_NC_var(varp);
        return status;
    }
    /* var = name nelems [dimid ...] vatt_list nc_type vsize begin
     *                                                       ^^^^^ */
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
        ncmpio_free_NC_var(varp);
        return status;
    }

    *varpp = varp;
    return NC_NOERR;
}

static int
val_get_NC_vararray(int          fd,
                    bufferinfo  *gbp,
                    NC_vararray *ncap,
                    int          f_ndims) /* no. dimensions defined in file */

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
    NC_tag tag = NC_UNSPECIFIED;
    int var;
    MPI_Offset tmp;
    unsigned long long err_addr;

    assert(gbp != NULL && gbp->pos != NULL);
    assert(ncap != NULL);
    assert(ncap->value == NULL); 

    /* read NC_tag (NC_VARIABLE or ZERO) from gbp buffer */
    status = val_get_NC_tag(fd, gbp, &tag);
    if (status != NC_NOERR) {
        printf("preamble of ");
        return status;
    }
 
    /* read nelems (number of variables) from gbp buffer */
    status = val_get_size_t(fd, gbp, &tmp);
    if (status != NC_NOERR) {
        printf("the length of ");
        return status;
    }
    if (tmp > NC_MAX_VARS) {
        /* number of allowable defined variables NC_MAX_VARS */
        printf("the length of ");
        return NC_EMAXVARS;
    }
    ncap->ndefined = (int)tmp;

    err_addr = (size_t)gbp->pos - (size_t)gbp->base + gbp->offset - gbp->size -
               (X_SIZEOF_INT + x_sizeof_NON_NEG);

    if (ncap->ndefined == 0) {
        /* From the CDF file format specification, the tag is either
         * NC_VARIABLE or ABSENT (ZERO), but we follow NetCDF library to skip
         * checking the tag when ndefined is zero.
         */
        return NC_NOERR;
#if 0
        if (tag != ABSENT) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component tag, while ABSENT is expected for ");
            DEBUG_RETURN(NC_ENOTNC)
        }
#endif
    } else {
        if (tag != NC_VARIABLE) {
            printf("Error @ [0x%8.8Lx]:\n", err_addr);
            printf("\tInvalid NC component tag, while NC_VARIABLE is expected as number of dimensions is %d for ", ncap->ndefined);
            DEBUG_RETURN(NC_ENOTNC)
        }
 
        size_t alloc_size = (size_t)ncap->ndefined + NC_ARRAY_GROWBY;
        ncap->value = (NC_var **) NCI_Calloc(alloc_size, sizeof(NC_var *));
        if (ncap->value == NULL) DEBUG_RETURN(NC_ENOMEM) 

        for (var=0; var<ncap->ndefined; var++) {
            status = val_get_NC_var(fd, gbp, &ncap->value[var], f_ndims);
            if (status != NC_NOERR) {
                printf("variable[%d] in ", var);
                ncap->ndefined = var;
                ncmpio_free_NC_vararray(ncap);
                return status;
            }
        }
    }

    return NC_NOERR;
}

/*----< NC_check_vlen() >----------------------------------------------------*/
/* Check whether variable size is less than or equal to vlen_max,
 * without overflowing in arithmetic calculations.  If OK, return 1,
 * else, return 0.  For CDF1 format or for CDF2 format on non-LFS
 * platforms, vlen_max should be 2^31 - 4, but for CDF2 format on
 * systems with LFS it should be 2^32 - 4.
 */
static int
NC_check_vlen(NC_var     *varp,
              MPI_Offset  vlen_max)
{
    int i;
    MPI_Offset prod=varp->xsz;     /* product of xsz and dimensions so far */

    for (i = IS_RECVAR(varp) ? 1 : 0; i < varp->ndims; i++) {
        if (varp->shape[i] > vlen_max / prod) {
            return 0;           /* size in bytes won't fit in a 32-bit int */
        }
        prod *= varp->shape[i];
    }
    return 1;
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
        vlen_max = X_INT64_MAX - 3; /* "- 3" handles rounded-up size */
    else if (ncp->flags & NC_64BIT_OFFSET) /* CDF2 format */
        vlen_max = X_UINT_MAX  - 3; /* "- 3" handles rounded-up size */
    else
        vlen_max = X_INT_MAX   - 3; /* CDF1 format */

    /* Loop through vars, first pass is for non-record variables */
    large_fix_vars_count = 0;
    rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.ndefined; ii++, vpp++) {
        if (!IS_RECVAR(*vpp)) {
            last = 0;
            if (NC_check_vlen(*vpp, vlen_max) == 0) {
                /* check this variable's shape product against vlen_max */
                if (ncp->format >= 5) /* CDF-5 */
                    DEBUG_RETURN(NC_EVARSIZE)
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
        DEBUG_RETURN(NC_EVARSIZE)
    }

    /* The only "too-large" variable must be the last one defined */
    if (large_fix_vars_count == 1 && last == 0) {
        printf("CDF-%d format allows only one large fixed-size variable and it must be the last one defined\n",ncp->format);
        DEBUG_RETURN(NC_EVARSIZE)
    }

    if (rec_vars_count == 0) return NC_NOERR;

    /* if there is a "too-large" fixed-size variable, no record variable is
     * allowed */
    if (large_fix_vars_count == 1) {
        printf("CDF-%d format allows only one large fixed-size variable when there is no record variable defined\n",ncp->format);
        DEBUG_RETURN(NC_EVARSIZE)
    }

    /* Loop through vars, second pass is for record variables.   */
    large_rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.ndefined; ii++, vpp++) {
        if (IS_RECVAR(*vpp)) {
            last = 0;
            if (NC_check_vlen(*vpp, vlen_max) == 0) {
                /* check this variable's shape product against vlen_max */
                if (ncp->format >= 5) /* CDF-5 */
                    DEBUG_RETURN(NC_EVARSIZE)
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
        DEBUG_RETURN(NC_EVARSIZE)
    }

    /* and it has to be the last one */
    if (large_rec_vars_count == 1 && last == 0) {
        printf("CDF-%d format allows only one large record variable and it must be the last one defined\n",ncp->format);
        DEBUG_RETURN(NC_EVARSIZE)
    }

    return NC_NOERR;
}

/*
 * Given a valid ncp, check all variables for their begins whether in an
 * increasing order.
 */
static int
val_NC_check_voff(NC *ncp)
{
    NC_var *varp;
    MPI_Offset i, prev, prev_off;

    if (ncp->vars.ndefined == 0) return NC_NOERR;

    /* Loop through vars, first pass is for non-record variables */
    prev_off = ncp->begin_var;
    prev     = 0;
    for (i=0; i<ncp->vars.ndefined; i++) {
        varp = ncp->vars.value[i];
        if (IS_RECVAR(varp)) continue;

        if (varp->begin < prev_off) {
            if (i == 0)
                printf("Variable \"%s\" begin offset (%lld) is less than header extent (%lld)\n", varp->name, varp->begin, prev_off);
            else
                printf("Variable \"%s\" begin offset (%lld) is less than previous variable \"%s\" end offset (%lld)\n", varp->name, varp->begin, ncp->vars.value[prev]->name, prev_off);
            DEBUG_RETURN(NC_ENOTNC)
        }
        prev_off = varp->begin + varp->len;
        prev = i;
    }

    if (ncp->begin_rec < prev_off) {
        printf("Record variable section begin offset (%lld) is less than fix-sized variable section end offset (%lld)\n", varp->begin, prev_off);
        DEBUG_RETURN(NC_ENOTNC)
    }

    /* Loop through vars, second pass is for record variables */
    prev_off = ncp->begin_rec;
    prev     = 0;
    for (i=0; i<ncp->vars.ndefined; i++) {
        varp = ncp->vars.value[i];
        if (!IS_RECVAR(varp)) continue;

        if (varp->begin < prev_off) {
            if (i == 0)
                printf("Variable \"%s\" begin offset (%lld) is less than record variable section begin offset (%lld)\n", varp->name, varp->begin, prev_off);
            else
                printf("Variable \"%s\" begin offset (%lld) is less than previous variable \"%s\" end offset (%lld)\n", varp->name, varp->begin, ncp->vars.value[prev]->name, prev_off);
            DEBUG_RETURN(NC_ENOTNC)
        }
        prev_off = varp->begin + varp->len;
        prev = i;
    }

    return NC_NOERR;
}

static int
val_get_NC(int fd, NC *ncp)
{
    int status;
    bufferinfo getbuf;
    char magic[sizeof(ncmagic)];
    MPI_Offset nrecs = 0;
    MPI_Aint pos_addr, base_addr;

    /* Initialize the get buffer that stores the header read from the file */
    getbuf.comm          = ncp->comm;
    getbuf.collective_fh = ncp->collective_fh;
    getbuf.offset        = 0;     /* read from start of the file */

    /* CDF-5's minimum header size is 4 bytes more than CDF-1 and CDF-2's */
    getbuf.size = _RNDUP( MAX(MIN_NC_XSZ+4, ncp->chunk), X_ALIGN );

    getbuf.pos = getbuf.base = (void *)NCI_Malloc(getbuf.size);

    /* Fetch the next header chunk. The chunk is 'gbp->size' bytes big
     * netcdf_file = header data
     * header      = magic numrecs dim_list gatt_list var_list
     */
    status = val_fetch(fd, &getbuf);
    if (status != NC_NOERR) goto fn_exit;
  
    /* First get the file format information, magic */
    memset(magic, 0, sizeof(magic));
    status = ncmpix_getn_text((const void **)(&getbuf.pos), sizeof(magic), magic);

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
    status = val_check_buffer(fd, &getbuf, (getbuf.version < 5) ? 4 : 8);
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
    status = val_get_NC_dimarray(fd, &getbuf, &ncp->dims);
    if (status != NC_NOERR) {
        printf("DIMENSION list!\n");
        goto fn_exit;
    }

    /* header = magic numrecs dim_list gatt_list var_list
     * att_list = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     * Check att_list
     */
    status = val_get_NC_attrarray(fd, &getbuf, &ncp->attrs); 
    if (status != NC_NOERR) {
        printf("GLOBAL ATTRIBUTE list!\n");
        goto fn_exit;
    }

    /* header = magic numrecs dim_list gatt_list var_list
     * var_list    = ABSENT | NC_VARIABLE   nelems  [var ...]
     * Check var_list
     */
    status = val_get_NC_vararray(fd, &getbuf, &ncp->vars, ncp->dims.ndefined);
    if (status != NC_NOERR) {
        printf("VARIABLE list!\n");
        goto fn_exit;
    }

    ncp->xsz = ncmpio_hdr_len_NC(ncp);

    status = compute_var_shape(ncp);
    if (status != NC_NOERR) goto fn_exit;

    status = val_NC_check_vlens(ncp);
    if (status != NC_NOERR) goto fn_exit;

    status = val_NC_check_voff(ncp);
    if (status != NC_NOERR) goto fn_exit;

fn_exit:
    NCI_Free(getbuf.base);

    return status;
}

/* End Of get NC */

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] file_name\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (exit 1 when fail, 0 success)\n"
    "       filename: input netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char **argv)
{
    extern int optind;
    char filename[512], *cmd;
    int i, rank, nprocs, verbose=1, fd, status=NC_NOERR;
    NC *ncp=NULL;
    struct stat ncfilestat;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    cmd = (char*) malloc(strlen(argv[0])+1);
    strcpy(cmd,argv[0]);

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'h':
            default:  if (rank==0) usage(cmd);
                      free(cmd);
                      MPI_Finalize();
                      return 1;
        }
    argc -= optind;
    argv += optind;
    if (argc != 1) {
        if (rank==0) usage(cmd);
        MPI_Finalize();
        free(cmd);
        return 1;
    }
    free(cmd);
    snprintf(filename, 512, "%s", argv[0]);

    fd = open(filename, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Error on open file %s (%s)\n",
                filename,strerror(errno));
        MPI_Finalize();
        return 1;
    }

    /* Allocate NC object */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) {
        status = NC_ENOMEM;
        printf("Error at line %d when calling ncmpio_new_NC()\n",__LINE__);
        goto prog_exit;
    }

    /* read and validate the header */
    status = val_get_NC(fd, ncp);
    if (status != NC_NOERR && status != -1) goto prog_exit;

    /* check data size */
    if (-1 == fstat(fd, &ncfilestat)) {
        printf("Error at line %d fstat (%s)\n",__LINE__,strerror(errno));
        status = NC_EFILE;
        goto prog_exit;
    }
    if (ncp->numrecs > 0) {
        MPI_Offset expect_fsize;
        expect_fsize = ncp->begin_rec + ncp->recsize * ncp->numrecs;
        if (expect_fsize < ncfilestat.st_size) {
            printf("Error: file size (%lld) is larger than expected (%lld)!\n",(long long)ncfilestat.st_size, expect_fsize);
            printf("\tbegin_rec=%lld recsize=%lld numrecs=%lld ncfilestat.st_size=%lld\n",ncp->begin_rec, ncp->recsize, ncp->numrecs, (long long) ncfilestat.st_size);
            status = NC_EFILE;
            goto prog_exit;
        }
        else if (expect_fsize > ncfilestat.st_size) {
            /* if file header are valid and the only error is the file size
             * less than expected, then this is due to partial data written
             * to the variable while the file is in no fill mode */
            if (verbose)
                printf("Warning: file size (%lld) is less than expected (%lld)!\n",(long long)ncfilestat.st_size, expect_fsize);
        }
    }
    else {
        MPI_Offset expect_fsize;
        if (ncp->vars.ndefined == 0)
            expect_fsize = ncp->xsz;
        else
            /* find the size of last fix-sized variable */
            expect_fsize = ncp->vars.value[ncp->vars.ndefined-1]->begin +
                           ncp->vars.value[ncp->vars.ndefined-1]->len;
        if (expect_fsize < ncfilestat.st_size) {
            printf("Error: file size (%lld) is larger than expected (%lld)!\n",(long long)ncfilestat.st_size, expect_fsize);
            status = NC_EFILE;
            goto prog_exit;
        }
        else if (expect_fsize > ncfilestat.st_size) {
            /* if file header are valid and the only error is the file size
             * less than expected, then this is due to partial data written
             * to the variable while the file is in no fill mode */
            if (verbose)
                printf("Warning: file size (%lld) is less than expected (%lld)!\n",(long long)ncfilestat.st_size, expect_fsize);
        }
    }

prog_exit:
    if (ncp != NULL) {
        if (ncp->dims.value  != NULL) NCI_Free(ncp->dims.value);
        if (ncp->attrs.value != NULL) NCI_Free(ncp->attrs.value);
        if (ncp->vars.value  != NULL) NCI_Free(ncp->vars.value);
        NCI_Free(ncp);
    }
    close(fd);

    if (status == NC_NOERR) {
        if (verbose)
            printf("File \"%s\" is a valid NetCDF file.\n",filename);
    }
    else
        printf("Error: %s\n",ncmpi_strerror(status));

    MPI_Finalize();

    exit((status == NC_NOERR) ? EXIT_SUCCESS : EXIT_FAILURE);
}
