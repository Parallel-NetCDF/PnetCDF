/*********************************************************************************
 *
 * This file is written by Northwestern University and Argonne National Laboratory
 *
 ********************************************************************************/
#include <mpi.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "nc.h"
#include "ncx.h"

#undef MAX  /* system may define MAX somewhere and complain */
#undef MIN  /* system may define MIN somewhere and complain */
#define MAX(mm,nn) (((mm) > (nn)) ? (mm) : (nn))
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))

/*
 * "magic number" at beginning of file: 0x43444601 (big endian) 
 */
static const schar ncmagic[] = {'C', 'D', 'F', 0x01}; 

/* Begin Of get NC */

/*
 * Fetch the next header chunk.
 */
int
val_fetch(bufferinfo *gbp, size_t fsize) {
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
    printf("Error @ [0x%8.8Lx]: \n\tUnexpected EOF, while ", gbp->offset);
    return -1;
  }
    
  gbp->size = bufsize;

  return ENOERR;
}

/*
 * Ensure that 'nextread' bytes are available.
 */
int
val_check_buffer(bufferinfo *gbp, size_t nextread) {
  if ((char *)gbp->pos + nextread <= (char *)gbp->base + gbp->size)
    return ENOERR;
  return val_fetch(gbp, MIN(gbp->size, nextread));
} 

int
val_get_NCtype(bufferinfo *gbp, NCtype *typep) {
  int type = 0;
  int status = val_check_buffer(gbp, X_SIZEOF_INT);
  if (status != ENOERR) {
    printf("NC component type is expected for ");
    return status;
  }

  status =  ncx_get_int_int(gbp->pos, &type);
  gbp->pos = (void *)((char *)gbp->pos + X_SIZEOF_INT);
  if (status != ENOERR)
    return status;
  *typep = (NCtype) type;
  return ENOERR;
}

int
val_get_size_t(bufferinfo *gbp, size_t *sp) {
  int status = val_check_buffer(gbp, X_SIZEOF_SIZE_T);
  if (status != ENOERR) {
    printf("size is expected for ");
    return status; 
  }
  return ncx_get_size_t((const void **)(&gbp->pos), sp);
}

int
val_get_NC_string(bufferinfo *gbp, NC_string **ncstrpp) {
  int status;
  size_t  nchars = 0, padding, bufremain, strcount; 
  NC_string *ncstrp;
  char *cpos;
  char pad[X_ALIGN-1];

  status = val_get_size_t(gbp, &nchars);
  if (status != ENOERR) {
    printf("the name string of ");
    return status;
  }

  ncstrp = new_NC_string(nchars, NULL);
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
      if(status != ENOERR) {
	printf("fetching the name string of ");
        free_NC_string(ncstrp);
        return status;
      } 
      bufremain = gbp->size;
    }
  }

  memset(pad, 0, X_ALIGN-1);
  status = val_check_buffer(gbp, padding);
  if(status != ENOERR) {
    printf("fetching padding for the name string of ");
    free_NC_string(ncstrp);
    return status;
  } 
  if (memcmp(gbp->pos, pad, padding) != 0) {
    printf("Error @ [0x%8.8Lx]: \n\tPadding should be 0x00 for the name string alignment of ",
	   ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size);
    free_NC_string(ncstrp);
    return EINVAL;
  }
  gbp->pos = (void *)((char *)gbp->pos + padding);
  
  *ncstrpp = ncstrp;
  
  return ENOERR;  
}

int
val_get_NC_dim(bufferinfo *gbp, NC_dim **dimpp) {
  int status;
  NC_string *ncstrp;
  NC_dim *dimp;

  status = val_get_NC_string(gbp, &ncstrp);
  if (status != ENOERR) 
    return status;

  dimp = new_x_NC_dim(ncstrp);
  if(dimp == NULL)
    return NC_ENOMEM;

  status = val_get_size_t(gbp, &dimp->size);
  if(status != ENOERR) {
    printf("\"%s\" - ", ncstrp->cp);
    free_NC_dim(dimp); /* frees name */
    return status;
  }

  *dimpp = dimp;

  return ENOERR;
}

int
val_get_NC_dimarray(bufferinfo *gbp, NC_dimarray *ncap) {
  int status;
  NCtype type = NC_UNSPECIFIED; 
  NC_dim **dpp, **end;
  int dim;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL);

  status = val_get_NCtype(gbp, &type);
  if(status != ENOERR) {
    printf("preamble of ");
    return status; 
  }

  status = val_get_size_t(gbp, &ncap->nelems);
  if(status != ENOERR) {
    printf("the length of ");
    return status;
  }

  if(ncap->nelems == 0) {
    if (type != NC_DIMENSION && type != NC_UNSPECIFIED) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
	      ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T);
      printf("NC_DIMENSION or NC_UNSPECIFIED is expected for ");
      return EINVAL;
    }
  } else {
    if(type != NC_DIMENSION) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
	      ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T);
      printf("NC_DIMENSION is expected since number of dimensions is %d for ", ncap->nelems);
      return EINVAL;
    }

    ncap->value = (NC_dim **) malloc(ncap->nelems * sizeof(NC_dim *));
    if(ncap->value == NULL)
      return NC_ENOMEM;
    ncap->nalloc = ncap->nelems;

    dpp = ncap->value;
    end = &dpp[ncap->nelems];
    for( /*NADA*/ dim = 0; dpp < end; dpp++, dim++) {
      status = val_get_NC_dim(gbp, dpp);
      if (status != ENOERR) {
	printf("dimension[%d] in ", dim);
        ncap->nelems = dpp - ncap->value;
        free_NC_dimarrayV(ncap);
        return status;
      }
    }
  }

  return ENOERR;
}

int
val_get_nc_type(bufferinfo *gbp, nc_type *typep) {
  int type = 0;
  int status = val_check_buffer(gbp, X_SIZEOF_INT);
  if(status != ENOERR) {
    printf("data type is expected for the values of ");
    return status;
  }

  status =  ncx_get_int_int(gbp->pos, &type);
  if(status != ENOERR) 
    return status;
  gbp->pos = (void *)((char *)gbp->pos + X_SIZEOF_INT); 

  if (   type != NC_BYTE
      && type != NC_CHAR
      && type != NC_SHORT
      && type != NC_INT
      && type != NC_FLOAT
      && type != NC_DOUBLE) {
    printf("Error @ [0x%8.8Lx]: \n\tUnknown data type for the values of ",
	   ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - X_SIZEOF_INT);
    return EINVAL; 
  }
 
  *typep = (nc_type) type;

  return ENOERR;
}

/*
 * Get the values of an attribute  
 */
int
val_get_NC_attrV(bufferinfo *gbp, NC_attr *attrp) {
  int status;
  void *value = attrp->xvalue;
  char pad[X_ALIGN-1]; 
  size_t nvalues = attrp->nelems, esz, padding, bufremain, attcount;

  esz = ncx_len_nctype(attrp->type);
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
      if(status != ENOERR) {
	printf("fetching the values of ");
        return status;
      }
      bufremain = gbp->size;
    }
  }
 
  memset(pad, 0, X_ALIGN-1);
  if (memcmp(gbp->pos, pad, padding) != 0) {
    printf("Error @ [0x%8.8Lx]: \n\tPadding should be 0x00 for the values alignment of ",
           ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size); 
    return EINVAL;
  }
  gbp->pos = (void *)((char *)gbp->pos + padding);

  return ENOERR;
}

int
val_get_NC_attr(bufferinfo *gbp, NC_attr **attrpp) {
  NC_string *strp;
  int status;
  nc_type type; 
  size_t nelems;
  NC_attr *attrp;

  status = val_get_NC_string(gbp, &strp);
  if(status != ENOERR)
    return status;

  status = val_get_nc_type(gbp, &type);
  if(status != ENOERR) {
    printf("\"%s\" - ", strp->cp);
    free_NC_string(strp);
    return status;
  }

  status = val_get_size_t(gbp, &nelems); 
  if(status != ENOERR) {
    printf("the values of \"%s\" - ", strp->cp);
    free_NC_string(strp);
    return status;
  }

  attrp = new_x_NC_attr(strp, type, nelems);
  if(attrp == NULL) {
    free_NC_string(strp);
    return status;
  }

  status = val_get_NC_attrV(gbp, attrp);
  if(status != ENOERR) {
    printf("\"%s\" - ", strp->cp);
    free_NC_attr(attrp); /* frees strp */ 
    return status;
  }

  *attrpp = attrp; 
  
  return ENOERR; 
}

int
val_get_NC_attrarray(bufferinfo *gbp, NC_attrarray *ncap){
  int status;
  NCtype type = NC_UNSPECIFIED;
  NC_attr **app, **end;
  int att;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL);

  status = val_get_NCtype(gbp, &type);
  if(status != ENOERR) {
    printf("preamble of ");
    return status; 
  }

  status = val_get_size_t(gbp, &ncap->nelems);
  if(status != ENOERR) {
    printf("the length of ");
    return status;
  }

  if(ncap->nelems == 0) {
    if (type != NC_ATTRIBUTE && type != NC_UNSPECIFIED) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
              ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T);
      printf("NC_ATTRIBUTE or NC_UNSPECIFIED is expected for "); 
      return EINVAL;
    }
  } else {
    if(type != NC_ATTRIBUTE) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
              ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T);
      printf("NC_ATTRIBUTE is expected since number of attributes is %d for ", ncap->nelems);  
      return EINVAL;
    }

    ncap->value = (NC_attr **) malloc(ncap->nelems * sizeof(NC_attr *));
    if(ncap->value == NULL)
      return NC_ENOMEM;
    ncap->nalloc = ncap->nelems; 

    app = ncap->value;
    end = &app[ncap->nelems];
    for( /*NADA*/ att = 0; app < end; app++, att++) {
      status = val_get_NC_attr(gbp, app);
      if (status != ENOERR) {
	printf("attribute[%d] of ", att);
        ncap->nelems = app - ncap->value;
        free_NC_attrarrayV(ncap);
        return status;
      }
    }
  }
  
  return ENOERR;
}

int
val_get_NC_var(bufferinfo *gbp, NC_var **varpp) {
  NC_string *strp;
  int status;
  size_t ndims, dim;
  NC_var *varp;

  status = val_get_NC_string(gbp, &strp);
  if(status != ENOERR)
    return status;

  status = val_get_size_t(gbp, &ndims);
  if(status != ENOERR) {
     printf("the dimid list of \"%s\" - ", strp->cp);
     free_NC_string(strp); 
     return status;
  }

  varp = new_x_NC_var(strp, ndims);
  if(varp == NULL) {
    free_NC_string(strp);
    return NC_ENOMEM;
  }

  for (dim = 0; dim < ndims; dim++ ) {
    status = val_check_buffer(gbp, X_SIZEOF_INT);
    if(status != ENOERR) {
      printf("the dimid[%d] is expected for \"%s\" - ", dim, strp->cp);
      free_NC_var(varp);
      return status;
    }
    status = ncx_getn_int_int((const void **)(&gbp->pos), 
                              1, varp->dimids + dim);
    if(status != ENOERR) {
      free_NC_var(varp);
      return status;
    }
  }

  status = val_get_NC_attrarray(gbp, &varp->attrs);
  if(status != ENOERR) {
    printf("ATTRIBUTE list of \"%s\" - ", strp->cp);
    free_NC_var(varp);
    return status;
  }

  status = val_get_nc_type(gbp, &varp->type);
  if(status != ENOERR) {
    printf("\"%s\" - ", strp->cp);
    free_NC_var(varp);
    return status;
  } 

  status = val_get_size_t(gbp, &varp->len);
  if(status != ENOERR) {
    printf("the data of  \"%s\" - ", strp->cp);
    free_NC_var(varp);
    return status;
  }

  status = val_check_buffer(gbp, X_SIZEOF_OFF_T);
  if(status != ENOERR) {
    printf("offset is expected for the data of \"%s\" - ", strp->cp);
    free_NC_var(varp);
    return status;
  }
  status = ncx_get_off_t((const void **)&gbp->pos,
                         &varp->begin);
  if(status != ENOERR) {
    free_NC_var(varp);
    return status;
  }

  *varpp = varp;
  return ENOERR;
}

int
val_get_NC_vararray(bufferinfo *gbp, NC_vararray *ncap) {
  int status;
  NCtype type = NC_UNSPECIFIED;
  NC_var **vpp, **end;
  int var;

  assert(gbp != NULL && gbp->pos != NULL);
  assert(ncap != NULL);
  assert(ncap->value == NULL); 

  status = val_get_NCtype(gbp, &type);
  if(status != ENOERR) {
    printf("preamble of ");
    return status;
  }
 
  status = val_get_size_t(gbp, &ncap->nelems);
  if(status != ENOERR) {
    printf("the length of ");
    return status;
  }
 
  if(ncap->nelems == 0) {
    if (type != NC_VARIABLE && type != NC_UNSPECIFIED) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
              ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T);
      printf("NC_VARIABLE or NC_UNSPECIFIED is expected for ");
      return EINVAL;
    }
  } else {
    if(type != NC_VARIABLE) {
      printf("Error @ [0x%8.8Lx]: \n\tInvalid NC component type, while ",
              ((size_t) gbp->pos - (size_t) gbp->base) + gbp->offset - gbp->size - 2 * X_SIZEOF_SIZE_T);
      printf("NC_VARIABLE is expected since number of variables is %d for ", ncap->nelems);        
      return EINVAL;
    }
 
    ncap->value = (NC_var **) malloc(ncap->nelems * sizeof(NC_var *));
    if(ncap->value == NULL)
      return NC_ENOMEM; 
    ncap->nalloc = ncap->nelems;

    vpp = ncap->value;
    end = &vpp[ncap->nelems];
    for( /*NADA*/ var = 0; vpp < end; vpp++, var++) {
      status = val_get_NC_var(gbp, vpp);
      if (status != ENOERR) {
        printf("variable[%d] in ", var);
        ncap->nelems = vpp - ncap->value;
        free_NC_vararrayV(ncap);
        return status;
      }
    }
  }

  return ENOERR;
}

int
val_get_NC(NC *ncp) {
  int status;
  bufferinfo getbuf;
  schar magic[sizeof(ncmagic)];
  size_t nrecs = 0;

  assert(ncp != NULL);

  /* Initialize the get buffer */

  getbuf.nciop = ncp->nciop;
  getbuf.offset = 0; 	/* read from start of the file */
  getbuf.size = _RNDUP( MAX(MIN_NC_XSZ, ncp->chunk), X_ALIGN );
  if (getbuf.size > 4096)
    getbuf.size = 4096;
  getbuf.pos = getbuf.base = (void *)malloc(getbuf.size);

  status = val_fetch(&getbuf, sizeof(magic));
  if(status != ENOERR) {
    printf("magic number (C D F \\001) is expected!\n");
    return status;
  }
  
  /* Get the header from get buffer */

  (void) memset(magic, 0, sizeof(magic));
  status = ncx_getn_schar_schar(
          (const void **)(&getbuf.pos), sizeof(magic), magic);
  if(memcmp(magic, ncmagic, sizeof(ncmagic)) != 0) {
    printf("Error @ [0x%8.8x]: \n\tUnknow magic number, while (C D F \\001) is expected!\n", 0);
    free(getbuf.base);
    return NC_ENOTNC;
  }

  status = val_check_buffer(&getbuf, X_SIZEOF_SIZE_T);
  if(status != ENOERR) {
    printf("number of records is expected!\n");
    free(getbuf.base);
    return status;
  }
  status = ncx_get_size_t((const void **)(&getbuf.pos), &nrecs);
  if(status != ENOERR) {
    free(getbuf.base);
    return status;
  }
  ncp->numrecs = nrecs;

  assert((char *)getbuf.pos < (char *)getbuf.base + getbuf.size);

  status = val_get_NC_dimarray(&getbuf, &ncp->dims);
  if(status != ENOERR) {
    printf("DIMENSION list!\n");
    free(getbuf.base);
    return status;
  }

  status = val_get_NC_attrarray(&getbuf, &ncp->attrs); 
  if(status != ENOERR) {
    printf("GLOBAL ATTRIBUTE list!\n");
    free(getbuf.base);
    return status;
  }

  status = val_get_NC_vararray(&getbuf, &ncp->vars);
  if(status != ENOERR) {
    printf("VARIABLE list!\n");
    free(getbuf.base);
    return status; 
  }

  ncp->xsz = hdr_len_NC(ncp); 
  status = NC_computeshapes(ncp);
  free(getbuf.base);

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

  ncp = new_NC(NULL);
  if(ncp == NULL) {
    printf("Not enough memory!\n");
    return 0; 
  }

  ncp->nciop = ncmpiio_new(ncfile, NC_NOWRITE);
  if(ncp->nciop == NULL) {
    free_NC(ncp);
    printf("Not enough memory!\n");
    return 0; 
  }

  if ( (*((int *)&ncp->nciop->fd) = open(ncfile, O_RDONLY)) < 0 ) {
    printf("Can not open file: %s\n", ncfile);
    ncmpiio_free(ncp->nciop);
    free_NC(ncp);
    return 0;
  }

  /* read to validate the header */

  status = val_get_NC(ncp);
  if (status !=  0) {
    close(ncp->nciop->fd);
    ncmpiio_free(ncp->nciop);
    free_NC(ncp);
    return 0;
  }

  /* check data size */
  
  fstat(ncp->nciop->fd, &ncfilestat);
  if ( ncp->begin_rec + ncp->recsize * ncp->numrecs < ncfilestat.st_size ) {
    printf("Error: \n\tData size is larger than defined!\n");
    close(ncp->nciop->fd);
    ncmpiio_free(ncp->nciop);
    free_NC(ncp);
    return 0;  
  } else if ( ncp->begin_rec + ncp->recsize * (ncp->numrecs - 1) >= ncfilestat.st_size ) {
    printf("Error: \n\tData size is less than expected!\n");
    close(ncp->nciop->fd);
    ncmpiio_free(ncp->nciop);
    free_NC(ncp);
    return 0;
  }


  /* close the file */

  close(ncp->nciop->fd);
  ncmpiio_free(ncp->nciop);
  free_NC(ncp);

  printf("The netCDF file is validated!\n");

  return 0;
}
