/*
 *	Copyright 1996, University Corporation for Atmospheric Research
 *      See netcdf/COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id$ */

#include "nc.h"
#include "rnd.h"
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "ncx.h"

/* list of open netcdf's */
static NC *NClist = NULL;

/* This is the default create format for ncmpi_create and nc__create. */
static int default_create_format = NC_FORMAT_CLASSIC;

/* These have to do with version numbers. */
#define MAGIC_NUM_LEN 4
#define VER_CLASSIC 1
#define VER_64BIT_OFFSET 2
#define VER_HDF5 3




/* Prototypes for functions used only in this file */
static int move_data_r(NC *ncp, NC *old);
static int move_recs_r(NC *ncp, NC *old);
static int move_vars_r(NC *ncp, NC *old);
static int write_NC(NC *ncp);
static int NC_begins(NC *ncp, size_t h_minfree, size_t v_align,
		      size_t v_minfree, size_t r_align);
static int NC_check_def(MPI_Comm comm, void *buf, size_t nn);

#if 0
static int enddef(NC *ncp);
static int nc_sync(int ncid);
static int nc_set_fill(int ncid, int fillmode, int *old_mode_ptr);
#endif

void
ncmpii_add_to_NCList(NC *ncp)
{
	assert(ncp != NULL);

	ncp->prev = NULL;
	if(NClist != NULL)
		NClist->prev = ncp;
	ncp->next = NClist;
	NClist = ncp;
}

void
ncmpii_del_from_NCList(NC *ncp)
{
	assert(ncp != NULL);

	if(NClist == ncp)
	{
		assert(ncp->prev == NULL);
		NClist = ncp->next;
	}
	else
	{
		assert(ncp->prev != NULL);
		ncp->prev->next = ncp->next;
	}

	if(ncp->next != NULL)
		ncp->next->prev = ncp->prev;

	ncp->next = NULL;
	ncp->prev = NULL;
}

/*
 * Check the data set definitions across all processes by
 * comparing the header buffer streams of all processes.
 */
static int
NC_check_def(MPI_Comm comm, void *buf, size_t nn) {
  int rank;
  int errcheck, compare = 0;
  void *cmpbuf;
  int max_size;

  MPI_Comm_rank(comm, &rank);

  if (rank == 0)
    max_size = nn;
  MPI_Bcast(&max_size, 1, MPI_INT, 0, comm);

  compare = max_size - nn;
  MPI_Allreduce(&compare, &errcheck, 1, MPI_INT, MPI_LOR, comm);
  if (errcheck)
    return NC_EMULTIDEFINE;

  if (rank == 0) 
    cmpbuf = buf;
  else
    cmpbuf = (void *)malloc(nn);

  MPI_Bcast(cmpbuf, nn, MPI_BYTE, 0, comm);

  if (rank != 0) {
    compare = memcmp(buf, cmpbuf, nn);
    free(cmpbuf);
  }

  MPI_Allreduce(&compare, &errcheck, 1, MPI_INT, MPI_LOR, comm);

  if (errcheck)
    return NC_EMULTIDEFINE;
  else 
    return NC_NOERR;
}

int
ncmpii_NC_check_id(int ncid, NC **ncpp)
{
	NC *ncp;

	if(ncid >= 0)
	{
		for(ncp = NClist; ncp != NULL; ncp = ncp->next)
		{
			if(ncp->nciop->fd == ncid)
			{
				*ncpp = ncp;
				return NC_NOERR; /* normal return */
			}
		}
	}

	/* else, not found */
	return NC_EBADID;
}


/* static */
void
ncmpii_free_NC(NC *ncp)
{
	if(ncp == NULL)
		return;
	ncmpii_free_NC_dimarrayV(&ncp->dims);
	ncmpii_free_NC_attrarrayV(&ncp->attrs);
	ncmpii_free_NC_vararrayV(&ncp->vars);
	free(ncp);
}


/* static */
NC *
ncmpii_new_NC(const size_t *chunkp)
{
	NC *ncp;

	ncp = (NC *) malloc(sizeof(NC));
	if(ncp == NULL)
		return NULL;
	(void) memset(ncp, 0, sizeof(NC));

	ncp->xsz = MIN_NC_XSZ;
	assert(ncp->xsz == ncmpii_hdr_len_NC(ncp, 0)); 
	
	ncp->chunk = chunkp != NULL ? *chunkp : NC_SIZEHINT_DEFAULT;

	return ncp;
}

/* This function sets a default create flag that will be logically
   or'd to whatever flags are passed into nc_create for all future
   calls to nc_create.
   Valid default create flags are NC_64BIT_OFFSET, NC_CLOBBER,
   NC_LOCK, NC_SHARE. */
int
ncmpi_set_default_format(int format, int *old_formatp)
{
    /* Return existing format if desired. */
    if (old_formatp)
      *old_formatp = default_create_format;

    /* Make sure only valid format is set. */
    if (format != NC_FORMAT_CLASSIC && format != NC_FORMAT_64BIT)
      return NC_EINVAL;
    default_create_format = format;
    return NC_NOERR;
}

/* returns a value suituable for a create flag.  Will return one or more of the
 * following values ORed together:
 * NC_64BIT_OFFSET, NC_CLOBBER, NC_LOCK, NC_SHARE */
int
ncmpii_get_default_format(void)
{
	return default_create_format;
}

/* static */
NC *
ncmpii_dup_NC(const NC *ref)
{
	NC *ncp;

	ncp = (NC *) malloc(sizeof(NC));
	if(ncp == NULL)
		return NULL;
	(void) memset(ncp, 0, sizeof(NC));

	if(ncmpii_dup_NC_dimarrayV(&ncp->dims, &ref->dims) != NC_NOERR)
		goto err;
	if(ncmpii_dup_NC_attrarrayV(&ncp->attrs, &ref->attrs) != NC_NOERR)
		goto err;
	if(ncmpii_dup_NC_vararrayV(&ncp->vars, &ref->vars) != NC_NOERR)
		goto err;

	ncp->xsz = ref->xsz;
	ncp->begin_var = ref->begin_var;
	ncp->begin_rec = ref->begin_rec;
	ncp->recsize = ref->recsize;
	NC_set_numrecs(ncp, NC_get_numrecs(ref));
	return ncp;
err:
	ncmpii_free_NC(ncp);
	return NULL;
}


/*
 *  Verify that this is a user nc_type
 * Formerly
NCcktype()
 * Sense of the return is changed.
 */
int
ncmpii_cktype(nc_type type)
{
	switch((int)type){
	case NC_BYTE:
	case NC_CHAR:
	case NC_SHORT:
	case NC_INT:
	case NC_FLOAT:
	case NC_DOUBLE:
		return(NC_NOERR);
	}
	return(NC_EBADTYPE);
}


/*
 * How many objects of 'type'
 * will fit into xbufsize?
 */
size_t
ncmpix_howmany(nc_type type, size_t xbufsize)
{
	switch(type){
	case NC_BYTE:
	case NC_CHAR:
		return xbufsize;
	case NC_SHORT:
		return xbufsize/X_SIZEOF_SHORT;
	case NC_INT:
		return xbufsize/X_SIZEOF_INT;
	case NC_FLOAT:
		return xbufsize/X_SIZEOF_FLOAT;
	case NC_DOUBLE:
		return xbufsize/X_SIZEOF_DOUBLE;
	default:
		assert("ncmpix_howmany: Bad type" == 0);
		return(0);
	}
}

#define	D_RNDUP(x, align) _RNDUP(x, (off_t)(align))

/*
 * Compute each variable's 'begin' offset,
 * update 'begin_rec' as well.
 */
static int
NC_begins(NC *ncp,
	size_t h_minfree, size_t v_align,
	size_t v_minfree, size_t r_align)
{
	size_t ii;
	size_t sizeof_off_t;
	off_t index = 0;
	NC_var **vpp;
	NC_var *last = NULL;

	if(v_align == NC_ALIGN_CHUNK)
		v_align = ncp->chunk;
	if(r_align == NC_ALIGN_CHUNK)
		r_align = ncp->chunk;

	if (fIsSet(ncp->flags, NC_64BIT_OFFSET)) {
		sizeof_off_t = 8;
	} else {
		sizeof_off_t = 4;
	}

	ncp->xsz = ncmpii_hdr_len_NC(ncp, sizeof_off_t);

	if(ncp->vars.nelems == 0) 
		return NC_NOERR;

	/* only (re)calculate begin_var if there is not sufficient space in header
	   or start of non-record variables is not aligned as requested by valign */
	if (ncp->begin_var < ncp->xsz + h_minfree ||
	    ncp->begin_var != D_RNDUP(ncp->begin_var, v_align) ) 
	{
	  index = (off_t) ncp->xsz;
	  ncp->begin_var = D_RNDUP(index, v_align);
	  if(ncp->begin_var < index + h_minfree)
	  {
	    ncp->begin_var = D_RNDUP(index + (off_t)h_minfree, v_align);
	  }
	}
	index = ncp->begin_var;

	/* loop thru vars, first pass is for the 'non-record' vars */
	vpp = ncp->vars.value;
	for(ii = 0; ii < ncp->vars.nelems ; ii++, vpp++)
	{
		if( IS_RECVAR(*vpp) )
		{
			/* skip record variables on this pass */
			continue;
		}
#if 0
fprintf(stderr, "    VAR %d %s: %ld\n", ii, (*vpp)->name->cp, (long)index);
#endif
		if (sizeof_off_t == 4 && (index > X_OFF_MAX || index < 0))
		{
			return NC_EVARSIZE;
		}

		(*vpp)->begin = index;
		index += (*vpp)->len;
	}

	/* only (re)calculate begin_rec if there is not sufficient
	   space at end of non-record variables or if start of record
	   variables is not aligned as requested by r_align */
	if (ncp->begin_rec < index + v_minfree ||
	    ncp->begin_rec != D_RNDUP(ncp->begin_rec, r_align) )
	{
	  ncp->begin_rec = D_RNDUP(index, r_align);
	  if(ncp->begin_rec < index + v_minfree)
	  {
	    ncp->begin_rec = D_RNDUP(index + (off_t)v_minfree, r_align);
	  }
	}
	index = ncp->begin_rec;

	ncp->recsize = 0;

	/* loop thru vars, second pass is for the 'record' vars */
	vpp = (NC_var **)ncp->vars.value;
	for(ii = 0; ii < ncp->vars.nelems; ii++, vpp++)
	{
		if( !IS_RECVAR(*vpp) )
		{
			/* skip non-record variables on this pass */
			continue;
		}

#if 0
fprintf(stderr, "    REC %d %s: %ld\n", ii, (*vpp)->name->cp, (long)index);
#endif
		if (sizeof_off_t == 4 && (index > X_OFF_MAX || index < 0))
		{
			return NC_EVARSIZE;
		}
		(*vpp)->begin = index;
		index += (*vpp)->len;
		/* check if record size must fit in 32-bits */
#if SIZEOF_OFF_T == SIZEOF_SIZE_T && SIZEOF_SIZE_T == 4
		if (ncp->recsize > X_UINT_MAX - (*vpp)->len)
		{
			return NC_EVARSIZE;
		}
#endif
		ncp->recsize += (*vpp)->len;
		last = (*vpp);
	}

	/*
	 * for special case of exactly one record variable, pack value
	 */
	if(last != NULL && ncp->recsize == last->len)
		ncp->recsize = *last->dsizes * last->xsz;

	if(NC_IsNew(ncp))
		NC_set_numrecs(ncp, 0);

	return NC_NOERR;
}

#define NC_NUMRECS_OFFSET 4
#define NC_NUMRECS_EXTENT 4

/*
 * Read just the numrecs member.
 * (A relatively expensive way to do things.)
 */

 
int
ncmpii_read_numrecs(NC *ncp) {
  int status = NC_NOERR, mpireturn;
  size_t nrecs;
  void *buf, *pos;
  MPI_Status mpistatus;
  int rank;

  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  assert(!NC_indef(ncp));
 
  pos = buf = (void *)malloc(X_SIZEOF_SIZE_T);

  /* reset the file view */
  mpireturn = MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE,
                    MPI_BYTE, "native", ncp->nciop->mpiinfo);
  if (mpireturn != MPI_SUCCESS) {
    char errorString[512];
    int  errorStringLen;
    MPI_Error_string(mpireturn, errorString, &errorStringLen);
    printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
    return NC_EREAD;
  }
 
  mpireturn = MPI_File_read_at(ncp->nciop->collective_fh, NC_NUMRECS_OFFSET,
                               buf, X_SIZEOF_SIZE_T, MPI_BYTE, &mpistatus);
 
  if (mpireturn != MPI_SUCCESS) {
    char errorString[512];
    int  errorStringLen;
    MPI_Error_string(mpireturn, errorString, &errorStringLen);
    printf("%2d: MPI_File_read_at error = %s\n", rank, errorString);
    return NC_EREAD;
  } 

  status = ncmpix_get_size_t((const void **)&pos, &nrecs);
  ncp->numrecs = nrecs;
 
  free(buf);
 
  return status;
}
 
/*
 * Write out just the numrecs member.
 * (A relatively expensive way to do things.)
 */

/*
 * Collective operation implicit
 */

int
ncmpii_write_numrecs(NC *ncp) {
  int status = NC_NOERR, mpireturn;
  size_t nrecs;
  void *buf, *pos; 
  MPI_Status mpistatus;
  MPI_Comm comm;
  int rank;
  
  assert(!NC_readonly(ncp));
  assert(!NC_indef(ncp));

  comm = ncp->nciop->comm;
  MPI_Comm_rank(comm, &rank);

  nrecs = ncp->numrecs;
  pos = buf = (void *)malloc(X_SIZEOF_SIZE_T);
  status = ncmpix_put_size_t(&pos, &nrecs);

  if(NC_indep(ncp) && NC_independentFhOpened(ncp->nciop)) {
    mpireturn = MPI_File_sync(ncp->nciop->independent_fh);
    if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_sync error = %s\n", rank, errorString);
        return NC_EWRITE;
    }
    MPI_Barrier(comm);
  }

  /* reset the file view */
  mpireturn = MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, 
		    MPI_BYTE, "native", ncp->nciop->mpiinfo);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

  if (rank == 0) {
    mpireturn = MPI_File_write_at(ncp->nciop->collective_fh, NC_NUMRECS_OFFSET,
				  buf, X_SIZEOF_SIZE_T, MPI_BYTE, &mpistatus); 
  } 

  MPI_Bcast(&mpireturn, 1, MPI_INT, 0, comm);

  if (mpireturn != MPI_SUCCESS) {
    char errorString[512];
    int  errorStringLen;
    MPI_Error_string(mpireturn, errorString, &errorStringLen);
    printf("%2d: MPI_File_write_at error = %s\n", rank, errorString);
    status = NC_EWRITE;
  } else {
    mpireturn = MPI_File_sync(ncp->nciop->collective_fh); 
    if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_sync error = %s\n", rank, errorString);
        status = NC_EWRITE;
    }

    MPI_Barrier(comm);
    fClr(ncp->flags, NC_NDIRTY);  
  }

  free(buf);

  return status;
}

/*
 * Read in the header
 * It is expensive.
 */

int
ncmpii_read_NC(NC *ncp) {
  int status = NC_NOERR;

  ncmpii_free_NC_dimarrayV(&ncp->dims);
  ncmpii_free_NC_attrarrayV(&ncp->attrs);
  ncmpii_free_NC_vararrayV(&ncp->vars); 

  status = ncmpii_hdr_get_NC(ncp);

  if(status == NC_NOERR)
    fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);

  return status;
}

/*
 * Write out the header
 */

static int
write_NC(NC *ncp)
{
  int status = NC_NOERR, mpireturn;
  void *buf;
  MPI_Status mpistatus;
  int rank;
 
  assert(!NC_readonly(ncp));
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

  buf = (void *)malloc(ncp->xsz);
  status = ncmpii_hdr_put_NC(ncp, buf);
  if(status != NC_NOERR) {
    free(buf);
    return status;
  }
  status = NC_check_def(ncp->nciop->comm, buf, ncp->xsz);
  if (status != NC_NOERR) {
    free(buf);
    return status;
  }

  /* reset the file view */
  mpireturn = MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE,
                    MPI_BYTE, "native", ncp->nciop->mpiinfo);
  if (mpireturn != MPI_SUCCESS) {
    char errorString[512];
    int  errorStringLen;
    MPI_Error_string(mpireturn, errorString, &errorStringLen);
    printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
    return NC_EWRITE;
  }

  if (rank == 0) {
    mpireturn = MPI_File_write_at(ncp->nciop->collective_fh, 0, buf, 
			          ncp->xsz, MPI_BYTE, &mpistatus);
  }

  MPI_Bcast(&mpireturn, 1, MPI_INT, 0, ncp->nciop->comm);

  if (mpireturn != MPI_SUCCESS) {
    char errorString[512];
    int  errorStringLen;
    MPI_Error_string(mpireturn, errorString, &errorStringLen);
    printf("%2d: MPI_File_write_at error = %s\n", rank, errorString);
    return NC_EWRITE;
  }

  fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);
  free(buf);
 
  return status;
} 

/*
 * Write the header or the numrecs if necessary.
 */
int
ncmpii_NC_sync(NC *ncp)
{
  int mynumrecs, numrecs;

  assert(!NC_readonly(ncp));

  /* collect and set the max numrecs due to difference by independent write */

  mynumrecs = ncp->numrecs;
  MPI_Allreduce(&mynumrecs, &numrecs, 1, MPI_INT, MPI_MAX, ncp->nciop->comm);
  if (numrecs > ncp->numrecs) {
    ncp->numrecs = numrecs;
    set_NC_ndirty(ncp);
  }

  if(NC_hdirty(ncp)) {
    return write_NC(ncp);
  }
  /* else */

  if(NC_ndirty(ncp)) {
    return ncmpii_write_numrecs(ncp);
  }
  /* else */

  return NC_NOERR;
}


/*
 * Move the records "out". 
 * Fill as needed.
 */

static int
move_data_r(NC *ncp, NC *old) {
  /* no new variable inserted, move the whole contiguous data part */
  ncp->numrecs = old->numrecs;
  return ncmpiio_move(ncp->nciop, ncp->begin_var, old->begin_var, 
              old->begin_rec - old->begin_var + old->recsize * old->numrecs);
}

static int
move_recs_r(NC *ncp, NC *old) {
  int status;
  int recno;
  const size_t old_nrecs = old->numrecs;
  const size_t ncp_recsize = ncp->recsize;
  const size_t old_recsize = old->recsize;
  const off_t ncp_off = ncp->begin_rec;
  const off_t old_off = old->begin_rec;

  assert(ncp_recsize >= old_recsize);

  if (ncp_recsize == old_recsize) {
    
    /* No new rec var inserted, move all rec vars as a whole */

    status = ncmpiio_move(ncp->nciop, ncp_off, old_off, 
                       old_recsize * old_nrecs);
    if(status != NC_NOERR)
      return status;
  } else {

    /* else, new rec var inserted, to be moved one record at a time */

    for (recno = (int)old_nrecs -1; recno >= 0; recno--) {
      status = ncmpiio_move(ncp->nciop, 
                         ncp_off+recno*ncp_recsize, 
                         old_off+recno*old_recsize, 
                         old_recsize);
      if(status != NC_NOERR)
        return status;
    } 
  }

  ncp->numrecs = old_nrecs;
  
  return NC_NOERR;
}


/*
 * Move the "non record" variables "out". 
 * Fill as needed.
 */

static int
move_vars_r(NC *ncp, NC *old) {
  return ncmpiio_move(ncp->nciop, ncp->begin_var, old->begin_var, 
                   old->begin_rec - old->begin_var); 
}

#if 0
 /*
 * Given a valid ncp, return NC_EVARSIZE if any variable has a bad len 
 * (product of non-rec dim sizes too large), else return NC_NOERR.
 */
static int
ncmpii_NC_check_vlens(NC *ncp)
{
    NC_var **vpp;
    /* maximum permitted variable size (or size of one record's worth
       of a record variable) in bytes.  This is different for format 1
       and format 2. */
    size_t vlen_max;
    size_t ii;
    size_t large_vars_count;
    size_t rec_vars_count;
    int last=-1;

    if(ncp->vars.nelems == 0) 
       return NC_NOERR;

    if ((ncp->flags & NC_64BIT_OFFSET) && sizeof(off_t) > 4) {
       /* CDF2 format and LFS */
       vlen_max = X_UINT_MAX - 3; /* "- 3" handles rounded-up size */
    } else {
       /* CDF1 format */
       vlen_max = X_INT_MAX - 3;
    }
    /* Loop through vars, first pass is for non-record variables.   */
    large_vars_count = 0;
    rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.nelems; ii++, vpp++) {
       if( !IS_RECVAR(*vpp) ) {
           last = 0;
           if( ncmpii_NC_check_vlen(*vpp, vlen_max) == 0 ) {
               large_vars_count++;
               last = 1;
           }
       } else {
         rec_vars_count++;
       }
    }
    /* OK if last non-record variable size too large, since not used to 
       compute an offset */
    if( large_vars_count > 1) { /* only one "too-large" variable allowed */
      return NC_EVARSIZE;
    }
    /* and it has to be the last one */ 
    if( large_vars_count == 1 && last == 0) { 
      return NC_EVARSIZE;
    }
    if( rec_vars_count > 0 ) {
       /* and if it's the last one, there can't be any record variables */
       if( large_vars_count == 1 && last == 1) {
           return NC_EVARSIZE;
       }
       /* Loop through vars, second pass is for record variables.   */
       large_vars_count = 0;
       vpp = ncp->vars.value;
       for (ii = 0; ii < ncp->vars.nelems; ii++, vpp++) {
           if( IS_RECVAR(*vpp) ) {
               last = 0;
               if( ncmpii_NC_check_vlen(*vpp, vlen_max) == 0 ) {
                   large_vars_count++;
                   last = 1;
               }
           }
       }
       /* OK if last record variable size too large, since not used to 
          compute an offset */
       if( large_vars_count > 1) { /* only one "too-large" variable allowed */
           return NC_EVARSIZE;
       }
       /* and it has to be the last one */ 
       if( large_vars_count == 1 && last == 0) { 
           return NC_EVARSIZE;
       }
    }
    return NC_NOERR;
}
#endif
 
int 
ncmpii_NC_enddef(NC *ncp) {
  int status = NC_NOERR;
  MPI_Comm comm;
  int mpireturn;
  int rank;

  assert(!NC_readonly(ncp));
  assert(NC_indef(ncp)); 

  comm = ncp->nciop->comm;

  MPI_Comm_rank(comm, &rank);

  NC_begins(ncp, 0, 1, 0, 1);

  /* serial netcdf calls a check on dimension lenghths here */

  /* To be updated */
  if(ncp->old != NULL) {
    /* a plain redef, not a create */
    assert(!NC_IsNew(ncp));
    assert(fIsSet(ncp->flags, NC_INDEF));
    assert(ncp->begin_rec >= ncp->old->begin_rec);
    assert(ncp->begin_var >= ncp->old->begin_var);
    assert(ncp->vars.nelems >= ncp->old->vars.nelems);
 
    mpireturn = MPI_File_sync(ncp->nciop->collective_fh);
    if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_sync error = %s\n", rank, errorString);
        return NC_EWRITE;
    }
    /*
     * Barrier needed switching between read and write 
     * Important, MPI_File_sync doesn't ensure barrier 
     */
    MPI_Barrier(comm); 

    if(ncp->vars.nelems != 0) {
      if(ncp->begin_rec > ncp->old->begin_rec) {
        if (ncp->vars.nelems == ncp->old->vars.nelems) {
          status = move_data_r(ncp, ncp->old);
          if(status != NC_NOERR)
            return status;
        } else {
          status = move_recs_r(ncp, ncp->old);
          if(status != NC_NOERR)
            return status;
          if(ncp->begin_var > ncp->old->begin_var) {
            status = move_vars_r(ncp, ncp->old);
            if(status != NC_NOERR)
              return status;
          }
        }
      } else { /* ... ncp->begin_rec > ncp->old->begin_rec */
        /* Even if (ncp->begin_rec == ncp->old->begin_rec)
         * and     (ncp->begin_var == ncp->old->begin_var)
         * might still have added a new record variable
         */
        if(ncp->recsize > ncp->old->recsize) {
          status = move_recs_r(ncp, ncp->old);
          if(status != NC_NOERR)
            return status;
        }
      }
    }
  } /* ... ncp->old != NULL */
 
  status = write_NC(ncp);
  if (status != NC_NOERR)
    return status;
 
  if(ncp->old != NULL) {
    ncmpii_free_NC(ncp->old);
    ncp->old = NULL;
  }
 
  fClr(ncp->flags, NC_CREAT | NC_INDEF);
  mpireturn = MPI_File_sync(ncp->nciop->collective_fh);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_sync error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

  MPI_Barrier(comm);

  return NC_NOERR;
}

#if 0
static int 
enddef(NC *ncp)
{
  assert(!NC_readonly(ncp));
  if(!NC_indef(ncp))
    return(NC_ENOTINDEFINE);

  NC_begins(ncp, 0, 1, 0, 1);

  if(ncp->old != NULL)
  {
    ncmpii_free_NC(ncp->old);
    ncp->old = NULL;
  }

  fClr(ncp->flags, NC_CREAT | NC_INDEF);

  return NC_NOERR;
}
#endif


/* Public */

int 
ncmpii_NC_close(NC *ncp) {
  int status = NC_NOERR;

  if(NC_indef(ncp)) {
    status = ncmpii_NC_enddef(ncp); /* TODO: defaults */
    if(status != NC_NOERR ) {
      /* To do: Abort new definition, if any */
      if (ncp->old != NULL) {
        ncmpii_free_NC(ncp->old);
        ncp->old = NULL;
        fClr(ncp->flags, NC_INDEF);
      }
    }
  }
  else if(!NC_readonly(ncp)) {
    status = ncmpii_NC_sync(ncp);
    if (status != NC_NOERR)
      return status;
  }
 
  (void) ncmpiio_close(ncp->nciop, 0);
  ncp->nciop = NULL;
 
  ncmpii_del_from_NCList(ncp);
 
  ncmpii_free_NC(ncp);
 
  return status;
}

int
ncmpi_inq(int ncid,
	int *ndimsp,
	int *nvarsp,
	int *nattsp,
	int *xtendimp)
{
	int status;
	NC *ncp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(ndimsp != NULL)
		*ndimsp = (int) ncp->dims.nelems;
	if(nvarsp != NULL)
		*nvarsp = (int) ncp->vars.nelems;
	if(nattsp != NULL)
		*nattsp = (int) ncp->attrs.nelems;
	if(xtendimp != NULL)
		*xtendimp = ncmpii_find_NC_Udim(&ncp->dims, NULL);

	return NC_NOERR;
}

int 
ncmpi_inq_ndims(int ncid, int *ndimsp)
{
	int status;
	NC *ncp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(ndimsp != NULL)
		*ndimsp = (int) ncp->dims.nelems;

	return NC_NOERR;
}

int 
ncmpi_inq_nvars(int ncid, int *nvarsp)
{
	int status;
	NC *ncp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(nvarsp != NULL)
		*nvarsp = (int) ncp->vars.nelems;

	return NC_NOERR;
}

int 
ncmpi_inq_natts(int ncid, int *nattsp)
{
	int status;
	NC *ncp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(nattsp != NULL)
		*nattsp = (int) ncp->attrs.nelems;

	return NC_NOERR;
}

int 
ncmpi_inq_unlimdim(int ncid, int *xtendimp)
{
	int status;
	NC *ncp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(xtendimp != NULL)
		*xtendimp = ncmpii_find_NC_Udim(&ncp->dims, NULL);

	return NC_NOERR;
}

#if 0
static int
nc_sync(int ncid)
{
	int status;
	NC *ncp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(NC_indef(ncp))
		return NC_EINDEFINE;

	if(NC_readonly(ncp))
	{
		return ncmpii_read_NC(ncp);
	}
	/* else, read/write */

	status = ncmpii_NC_sync(ncp);
	if(status != NC_NOERR)
		return status;

	return ncp->nciop->sync(ncp->nciop);
}


static int
nc_set_fill(int ncid,
	int fillmode, int *old_mode_ptr)
{
	int status;
	NC *ncp;
	int oldmode;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	oldmode = fIsSet(ncp->flags, NC_NOFILL) ? NC_NOFILL : NC_FILL;

	if(fillmode == NC_NOFILL)
	{
		fSet(ncp->flags, NC_NOFILL);
	}
	else if(fillmode == NC_FILL)
	{
		if(fIsSet(ncp->flags, NC_NOFILL))
		{
			/*
			 * We are changing back to fill mode
			 * so do a sync
			 */
			status = ncmpii_NC_sync(ncp);
			if(status != NC_NOERR)
				return status;
		}
		fClr(ncp->flags, NC_NOFILL);
	}
	else
	{
		return NC_EINVAL; /* Invalid fillmode */
	}

	if(old_mode_ptr != NULL)
		*old_mode_ptr = oldmode;

	return NC_NOERR;
}
#endif

/*ARGSUSED*/

