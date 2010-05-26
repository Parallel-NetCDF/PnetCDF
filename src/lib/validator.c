/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#include "nc.h"
#include <mpi.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include "ncx.h"
#include "macro.h"

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

static int val_fetch(bufferinfo *gbp, MPI_Offset fsize);
static int val_check_buffer(bufferinfo *gbp, MPI_Offset nextread);

/* Begin Of get NC */

/*
 * Fetch the next header chunk.
 */
static int
val_fetch(bufferinfo *gbp, MPI_Offset fsize) {
  char *buf;
  ssize_t nn = 0, bufsize = 0;

  assert(gbp->base != NULL);
  
  fsize = _RNDUP(fsize, X_ALIGN);
  (void) memset(gbp->base, 0, gbp->size);
  buf = gbp->pos = gbp->base;

  lseek(gbp->nciop->fd, gbp->offset, SEEK_SET);
  while ( (bufsize < gbp->size) && (nn = read(gbp->nciop->fd, buf, gbp->size-bufsize)) > 0 ) {
    buf += nn;
    bufsize += nn;
  }
  gbp->offset += bufsize; 

  if (bufsize < fsize) {
    printf("Error @ [0x%8.8Lx]: \n\tUnexpected EOF, while ",
	   (long long unsigned) gbp->offset);
    return -1;
  }
    
  gbp->size = bufsize;

  return NC_NOERR;
}

/*
 * Ensure that 'nextread' bytes are available.
 */
static int
val_check_buffer(bufferinfo *gbp, MPI_Offset nextread) {
  if ((char *)gbp->pos + nextread <= (char *)gbp->base + gbp->size)
    return NC_NOERR;
  return val_fetch(gbp, MIN(gbp->size, nextread));
} 

static int
val_get_NCtype(bufferinfo *gbp, NCtype *typep) {
  int type = 0;
  int status = val_check_buffer(gbp, X_SIZEOF_INT);
  if (status != NC_NOERR) {
    printf("NC component type is expected for ");
    return status;
  }

  status =  ncmpix_get_int_int(gbp->pos, &type);
  gbp->pos = (void *)((char *)gbp->pos + X_SIZEOF_INT);
  if (status != NC_NOERR)
    return status;
  *typep = (NCtype) type;
  return NC_NOERR;
}

static int
val_get_size_t(bufferinfo *gbp, MPI_Offset *sp) {
  int sizeof_t = gbp->version == 5 ? 8 : 4; 
  int status = val_check_buffer(gbp, sizeof_t);
  if (status != NC_NOERR) {
    printf("size is expected for ");
    return status; 
  }
  return ncmpix_get_size_t((const void **)(&gbp->pos), sp, sizeof_t );
}

static int
val_get_NC_string(bufferinfo *gbp, NC_string **ncstrpp) {
  int status;
  MPI_Offset  nchars = 0, padding, bufremain, strcount; 
  NC_string *ncstrp;
  char *cpos;
  char pad[X_ALIGN-1];

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
  bufremain = gbp->size - (size_t)((char *)gbp->pos - (char *)gbp->base);
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
      status = val_fetch(gbp, MIN(gbp->size, X_SIZEOF_CHAR * nchars));
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
    return EINVAL;
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
  if (status != NC_NOERR) 
    return status;

  dimp = ncmpii_new_x_NC_dim(ncstrp);
  if(dimp == NULL)
    return NC_ENOMEM;

  status = val_get_size_t(gbp, &dimp->size);
  if(status != NC_NOERR) {
    printf("\"%s\" - ", ncstrp->cp);
    ncmpii_free_NC_dim(dimp); /* frees name */
    return status;
  }

  *dimpp = dimp;

  return NC_NOERR;
}

static int
val_get_NC_dimarray(bufferinfo *gbp, NC_dimarray *ncap) {
  int status;
  NCtype type = NC_UNSPECIFIED; 
  NC_dim **dpp, **end;
  int dim;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL);

  status = val_get_NCtype(gbp, &type);
  if(status != NC_NOERR) {
    printf("preamble of ");
    return status; 
  }

  status = val_get_size_t(gbp, &ncap->nelems);
  if(status != NC_NOERR) {
    printf("the length of ");
    return status;
  }

  if(ncap->nelems == 0) {
    if (type != NC_DIMENSION && type != NC_UNSPECIFIED) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
	      (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T));
      printf("NC_DIMENSION or NC_UNSPECIFIED is expected for ");
      return EINVAL;
    }
  } else {
    if(type != NC_DIMENSION) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
	      (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T));
      printf("NC_DIMENSION is expected since number of dimensions is %d for ", (int)ncap->nelems);
      return EINVAL;
    }

    ncap->value = (NC_dim **) NCI_Malloc(ncap->nelems * sizeof(NC_dim *));
    if(ncap->value == NULL)
      return NC_ENOMEM;
    ncap->nalloc = ncap->nelems;

    dpp = ncap->value;
    end = &dpp[ncap->nelems];
    for( /*NADA*/ dim = 0; dpp < end; dpp++, dim++) {
      status = val_get_NC_dim(gbp, dpp);
      if (status != NC_NOERR) {
	printf("dimension[%d] in ", dim);
        ncap->nelems = dpp - ncap->value;
        ncmpii_free_NC_dimarrayV(ncap);
        return status;
      }
    }
  }

  return NC_NOERR;
}

static int
val_get_nc_type(bufferinfo *gbp, nc_type *typep) {
  int type = 0;
  int status = val_check_buffer(gbp, X_SIZEOF_INT);
  if(status != NC_NOERR) {
    printf("data type is expected for the values of ");
    return status;
  }

  status =  ncmpix_get_int_int(gbp->pos, &type);
  if(status != NC_NOERR) 
    return status;
  gbp->pos = (void *)((char *)gbp->pos + X_SIZEOF_INT); 

  if (   type != NC_BYTE
      && type != NC_CHAR
      && type != NC_SHORT
      && type != NC_INT
      && type != NC_FLOAT
      && type != NC_DOUBLE) {
    printf("Error @ [0x%8.8Lx]: \n\tUnknown data type for the values of ",
	   (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - X_SIZEOF_INT));
    return EINVAL; 
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

  esz = ncmpix_len_nctype(attrp->type);
  padding = attrp->xsz - esz * nvalues;
  bufremain = gbp->size - (size_t)((char *)gbp->pos - (char *)gbp->base);

  while (nvalues > 0) {
    if (bufremain > 0) {
      attcount = MIN(bufremain, esz * nvalues);
      (void) memcpy(value, gbp->pos, attcount);
      nvalues -= attcount/esz;
      gbp->pos = (void *)((char *)gbp->pos + attcount);
      value = (void *)((char *)value + attcount);
      bufremain -= attcount;
    } else {
      status = val_fetch(gbp, MIN(gbp->size, esz * nvalues));
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
    return EINVAL;
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
val_get_NC_attrarray(bufferinfo *gbp, NC_attrarray *ncap){
  int status;
  NCtype type = NC_UNSPECIFIED;
  NC_attr **app, **end;
  int att;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL);

  status = val_get_NCtype(gbp, &type);
  if(status != NC_NOERR) {
    printf("preamble of ");
    return status; 
  }

  status = val_get_size_t(gbp, &ncap->nelems);
  if(status != NC_NOERR) {
    printf("the length of ");
    return status;
  }

  if(ncap->nelems == 0) {
    if (type != NC_ATTRIBUTE && type != NC_UNSPECIFIED) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
              (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T));
      printf("NC_ATTRIBUTE or NC_UNSPECIFIED is expected for "); 
      return EINVAL;
    }
  } else {
    if(type != NC_ATTRIBUTE) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
              (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T));
      printf("NC_ATTRIBUTE is expected since number of attributes is %d for ", (int)ncap->nelems);  
      return EINVAL;
    }

    ncap->value = (NC_attr **) NCI_Malloc(ncap->nelems * sizeof(NC_attr *));
    if(ncap->value == NULL)
      return NC_ENOMEM;
    ncap->nalloc = ncap->nelems; 

    app = ncap->value;
    end = &app[ncap->nelems];
    for( /*NADA*/ att = 0; app < end; app++, att++) {
      status = val_get_NC_attr(gbp, app);
      if (status != NC_NOERR) {
	printf("attribute[%d] of ", att);
        ncap->nelems = app - ncap->value;
        ncmpii_free_NC_attrarrayV(ncap);
        return status;
      }
    }
  }
  
  return NC_NOERR;
}

static int
val_get_NC_var(bufferinfo *gbp, NC_var **varpp) {
  NC_string *strp;
  int status;
  MPI_Offset ndims, *tmp_dim;
  size_t dim;
  NC_var *varp;

  status = val_get_NC_string(gbp, &strp);
  if(status != NC_NOERR)
    return status;

  status = val_get_size_t(gbp, &ndims);
  if(status != NC_NOERR) {
     printf("the dimid list of \"%s\" - ", strp->cp);
     ncmpii_free_NC_string(strp); 
     return status;
  }

  varp = ncmpii_new_x_NC_var(strp, ndims);
  if(varp == NULL) {
    ncmpii_free_NC_string(strp);
    return NC_ENOMEM;
  }

  for (dim = 0; dim < ndims; dim++ ) {
    status = val_check_buffer(gbp, (gbp->version == 5 ? 8 : 4));
    if(status != NC_NOERR) {
      printf("the dimid[%d] is expected for \"%s\" - ", (int)dim, strp->cp);
      ncmpii_free_NC_var(varp);
      return status;
    }
    tmp_dim = (MPI_Offset*) (varp->dimids + dim);
    status = ncmpix_getn_long_long((const void **)(&gbp->pos), 
                              1, tmp_dim);
    if(status != NC_NOERR) {
      ncmpii_free_NC_var(varp);
      return status;
    }
  }

  status = val_get_NC_attrarray(gbp, &varp->attrs);
  if(status != NC_NOERR) {
    printf("ATTRIBUTE list of \"%s\" - ", strp->cp);
    ncmpii_free_NC_var(varp);
    return status;
  }

  status = val_get_nc_type(gbp, &varp->type);
  if(status != NC_NOERR) {
    printf("\"%s\" - ", strp->cp);
    ncmpii_free_NC_var(varp);
    return status;
  } 

  status = val_get_size_t(gbp, &varp->len);
  if(status != NC_NOERR) {
    printf("the data of  \"%s\" - ", strp->cp);
    ncmpii_free_NC_var(varp);
    return status;
  }

  status = val_check_buffer(gbp, (gbp->version == 5 ? 8 : 4));
  if(status != NC_NOERR) {
    printf("offset is expected for the data of \"%s\" - ", strp->cp);
    ncmpii_free_NC_var(varp);
    return status;
  }
  status = ncmpix_get_off_t((const void **)&gbp->pos,
                         &varp->begin, (gbp->version == 1 ? 4 : 8));
  if(status != NC_NOERR) {
    ncmpii_free_NC_var(varp);
    return status;
  }

  *varpp = varp;
  return NC_NOERR;
}

static int
val_get_NC_vararray(bufferinfo *gbp, NC_vararray *ncap) {
  int status;
  NCtype type = NC_UNSPECIFIED;
  NC_var **vpp, **end;
  int var;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL); 

  status = val_get_NCtype(gbp, &type);
  if(status != NC_NOERR) {
    printf("preamble of ");
    return status;
  }
 
  status = val_get_size_t(gbp, &ncap->nelems);
  if(status != NC_NOERR) {
    printf("the length of ");
    return status;
  }
 
  if(ncap->nelems == 0) {
    if (type != NC_VARIABLE && type != NC_UNSPECIFIED) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
              (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T));
      printf("NC_VARIABLE or NC_UNSPECIFIED is expected for ");
      return EINVAL;
    }
  } else {
    if(type != NC_VARIABLE) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
              (long long unsigned) (((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T));
      printf("NC_VARIABLE is expected since number of variables is %d for ", (int)ncap->nelems);        
      return EINVAL;
    }
 
    ncap->value = (NC_var **) NCI_Malloc(ncap->nelems * sizeof(NC_var *));
    if(ncap->value == NULL)
      return NC_ENOMEM; 
    ncap->nalloc = ncap->nelems;

    vpp = ncap->value;
    end = &vpp[ncap->nelems];
    for( /*NADA*/ var = 0; vpp < end; vpp++, var++) {
      status = val_get_NC_var(gbp, vpp);
      if (status != NC_NOERR) {
        printf("variable[%d] in ", var);
        ncap->nelems = vpp - ncap->value;
        ncmpii_free_NC_vararrayV(ncap);
        return status;
      }
    }
  }

  return NC_NOERR;
}

static int
val_get_NC(NC *ncp) {
  int status;
  bufferinfo getbuf;
  schar magic[sizeof(ncmagic)];
  MPI_Offset nrecs = 0;

  assert(ncp != NULL);

  /* Initialize the get buffer */

  getbuf.nciop = ncp->nciop;
  getbuf.offset = 0; 	/* read from start of the file */
  getbuf.size = _RNDUP( MAX(MIN_NC_XSZ, ncp->chunk), X_ALIGN );
  if (getbuf.size > 4096)
    getbuf.size = 4096;
  getbuf.pos = getbuf.base = (void *)NCI_Malloc(getbuf.size);

  status = val_fetch(&getbuf, sizeof(magic));
  if(status != NC_NOERR) {
    printf("magic number (C D F \\001) is expected!\n");
    return status;
  }
  
  /* Get the header from get buffer */

  (void) memset(magic, 0, sizeof(magic));
  status = ncmpix_getn_schar_schar(
          (const void **)(&getbuf.pos), sizeof(magic), magic);
  if(memcmp(magic, ncmagic, sizeof(ncmagic)) != 0) {
    printf("Error @ [0x%8.8x]: \n\tUnknow magic number, while (C D F \\001) is expected!\n", (unsigned) 0);
    NCI_Free(getbuf.base);
    return NC_ENOTNC;
  }

  status = val_check_buffer(&getbuf, X_SIZEOF_SIZE_T);
  if(status != NC_NOERR) {
    printf("number of records is expected!\n");
    NCI_Free(getbuf.base);
    return status;
  }
  status = ncmpix_get_size_t((const void **)(&getbuf.pos), &nrecs, getbuf.version == 5 ? 8 : 4);
  if(status != NC_NOERR) {
    NCI_Free(getbuf.base);
    return status;
  }
  ncp->numrecs = nrecs;

  assert((char *)getbuf.pos < (char *)getbuf.base + getbuf.size);

  status = val_get_NC_dimarray(&getbuf, &ncp->dims);
  if(status != NC_NOERR) {
    printf("DIMENSION list!\n");
    NCI_Free(getbuf.base);
    return status;
  }

  status = val_get_NC_attrarray(&getbuf, &ncp->attrs); 
  if(status != NC_NOERR) {
    printf("GLOBAL ATTRIBUTE list!\n");
    NCI_Free(getbuf.base);
    return status;
  }

  status = val_get_NC_vararray(&getbuf, &ncp->vars);
  if(status != NC_NOERR) {
    printf("VARIABLE list!\n");
    NCI_Free(getbuf.base);
    return status; 
  }

  ncp->xsz = ncmpii_hdr_len_NC(ncp, (getbuf.version == 1 ? 4 : 8)); 
  status = ncmpii_NC_computeshapes(ncp);
  NCI_Free(getbuf.base);

  return status;
}

/* End Of get NC */

int 
main(int argc, char **argv) {

  char *ncfile;
  int status;
  NC *ncp;
  struct stat ncfilestat;

  if (argc < 2) {
    printf("Missing ncfile name. Usage:\n\t ncvalid <ncfile>\n");
    exit(1);
  } 

  if (argc > 2) {
    printf("Too many arguments. Usage:\n\t ncvalid <ncfile>\n");
    exit(1);
  }

  ncfile = argv[1];

  /* open the netCDF file */

  ncp = ncmpii_new_NC(NULL);
  if(ncp == NULL) {
    printf("Not enough memory!\n");
    return 0; 
  }

  ncp->nciop = ncmpiio_new(ncfile, NC_NOWRITE);
  if(ncp->nciop == NULL) {
    ncmpii_free_NC(ncp);
    printf("Not enough memory!\n");
    return 0; 
  }

  if ( (*((int *)&ncp->nciop->fd) = open(ncfile, O_RDONLY)) < 0 ) {
    printf("Can not open file: %s\n", ncfile);
    ncmpiio_free(ncp->nciop);
    ncmpii_free_NC(ncp);
    return 0;
  }

  /* read to validate the header */

  status = val_get_NC(ncp);
  if (status !=  0) {
    close(ncp->nciop->fd);
    ncmpiio_free(ncp->nciop);
    ncmpii_free_NC(ncp);
    return 0;
  }

  /* check data size */
  
  fstat(ncp->nciop->fd, &ncfilestat);
  if ( ncp->begin_rec + ncp->recsize * ncp->numrecs < ncfilestat.st_size ) {
    printf("Error: \n\tData size is larger than defined!\n");
    close(ncp->nciop->fd);
    ncmpiio_free(ncp->nciop);
    ncmpii_free_NC(ncp);
    return 0;  
  } else if ( ncp->begin_rec + ncp->recsize * (ncp->numrecs - 1) >= ncfilestat.st_size ) {
    printf("Error: \n\tData size is less than expected!\n");
    close(ncp->nciop->fd);
    ncmpiio_free(ncp->nciop);
    ncmpii_free_NC(ncp);
    return 0;
  }


  /* close the file */

  close(ncp->nciop->fd);
  ncmpiio_free(ncp->nciop);
  ncmpii_free_NC(ncp);

  printf("The netCDF file is validated!\n");

  return 0;
}
