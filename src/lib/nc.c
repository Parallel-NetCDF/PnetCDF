/*
 *	Copyright 1996, University Corporation for Atmospheric Research
 *      See netcdf/COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id$ */

#include "nc.h"
#include "rnd.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ncx.h"
#if defined(LOCKNUMREC) /* && _CRAYMPP */
#  include <mpp/shmem.h>
#  include <intrinsics.h>
#endif

/* list of open netcdf's */
static NC *NClist = NULL;

/* static */
void
add_to_NCList(NC *ncp)
{
	assert(ncp != NULL);

	ncp->prev = NULL;
	if(NClist != NULL)
		NClist->prev = ncp;
	ncp->next = NClist;
	NClist = ncp;
}

/* static */
void
del_from_NCList(NC *ncp)
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

#ifdef USE_MPIO /* Following interface added by Jianwei Li */

/*
 * Check the data set definitions across all processes by
 * comparing the header buffer streams of all processes.
 */
int
NC_check_def(MPI_Comm comm, void *buf, size_t nn) {
  int rank;
  int errcheck, compare = 0;
  void *cmpbuf;

  MPI_Comm_rank(comm, &rank);

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

#endif /* USE_MPIO */

int
NC_check_id(int ncid, NC **ncpp)
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
free_NC(NC *ncp)
{
	if(ncp == NULL)
		return;
	free_NC_dimarrayV(&ncp->dims);
	free_NC_attrarrayV(&ncp->attrs);
	free_NC_vararrayV(&ncp->vars);
#if _CRAYMPP && defined(LOCKNUMREC)
	shfree(ncp);
#else
	free(ncp);
#endif /* _CRAYMPP && LOCKNUMREC */
}


/* static */
NC *
new_NC(const size_t *chunkp)
{
	NC *ncp;

#if _CRAYMPP && defined(LOCKNUMREC)
	ncp = (NC *) shmalloc(sizeof(NC));
#else
	ncp = (NC *) malloc(sizeof(NC));
#endif /* _CRAYMPP && LOCKNUMREC */
	if(ncp == NULL)
		return NULL;
	(void) memset(ncp, 0, sizeof(NC));

	ncp->xsz = MIN_NC_XSZ;
	assert(ncp->xsz == hdr_len_NC(ncp));
	
	ncp->chunk = chunkp != NULL ? *chunkp : NC_SIZEHINT_DEFAULT;

	return ncp;
}


/* static */
NC *
dup_NC(const NC *ref)
{
	NC *ncp;

#if _CRAYMPP && defined(LOCKNUMREC)
	ncp = (NC *) shmalloc(sizeof(NC));
#else
	ncp = (NC *) malloc(sizeof(NC));
#endif /* _CRAYMPP && LOCKNUMREC */
	if(ncp == NULL)
		return NULL;
	(void) memset(ncp, 0, sizeof(NC));

	if(dup_NC_dimarrayV(&ncp->dims, &ref->dims) != NC_NOERR)
		goto err;
	if(dup_NC_attrarrayV(&ncp->attrs, &ref->attrs) != NC_NOERR)
		goto err;
	if(dup_NC_vararrayV(&ncp->vars, &ref->vars) != NC_NOERR)
		goto err;

	ncp->xsz = ref->xsz;
	ncp->begin_var = ref->begin_var;
	ncp->begin_rec = ref->begin_rec;
	ncp->recsize = ref->recsize;
	NC_set_numrecs(ncp, NC_get_numrecs(ref));
	return ncp;
err:
	free_NC(ncp);
	return NULL;
}


/*
 *  Verify that this is a user nc_type
 * Formerly
NCcktype()
 * Sense of the return is changed.
 */
int
nc_cktype(nc_type type)
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
ncx_howmany(nc_type type, size_t xbufsize)
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
	}
	assert("ncx_howmany: Bad type" == 0);
	return(0);
}

#define	D_RNDUP(x, align) _RNDUP(x, (off_t)(align))

/*
 * Compute each variable's 'begin' offset,
 * update 'begin_rec' as well.
 */
/* static */
void
NC_begins(NC *ncp,
	size_t h_minfree, size_t v_align,
	size_t v_minfree, size_t r_align)
{
	size_t ii;
	off_t index = 0;
	NC_var **vpp;
	NC_var *last = NULL;

	if(v_align == NC_ALIGN_CHUNK)
		v_align = ncp->chunk;
	if(r_align == NC_ALIGN_CHUNK)
		r_align = ncp->chunk;

	ncp->xsz = hdr_len_NC(ncp);

	if(ncp->vars.nelems == 0) 
		return;

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
		(*vpp)->begin = index;
		index += (*vpp)->len;
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

}

#define NC_NUMRECS_OFFSET 4
#define NC_NUMRECS_EXTENT 4

/*
 * Read just the numrecs member.
 * (A relatively expensive way to do things.)
 */

#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */
 
int
read_numrecs(NC *ncp) {
  int status = NC_NOERR, mpireturn;
  size_t nrecs;
  void *buf, *pos;
  MPI_Status mpistatus;
 
  assert(!NC_indef(ncp));
 
  pos = buf = (void *)malloc(X_SIZEOF_SIZE_T);

  /* reset the file view */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE,
                    MPI_BYTE, "native", ncp->nciop->mpiinfo);
 
  mpireturn = MPI_File_read_at(ncp->nciop->collective_fh, NC_NUMRECS_OFFSET,
                               buf, X_SIZEOF_SIZE_T, MPI_BYTE, &mpistatus);
 
  if (mpireturn != MPI_SUCCESS) {
    status = NC_EWRITE;
  } 

  status = ncx_get_size_t((const void **)&pos, &nrecs);
  ncp->numrecs = nrecs;
 
  free(buf);
 
  return status;
}
 
#else /* Following interface modified as above by Jianwei Li */

int
read_numrecs(NC *ncp)
{
	int status = NC_NOERR;
	const void *xp;
	size_t nrecs = NC_get_numrecs(ncp);

	assert(!NC_indef(ncp));

	status = ncp->nciop->get(ncp->nciop,
		 NC_NUMRECS_OFFSET, NC_NUMRECS_EXTENT, 0, (void **)&xp);
					/* cast away const */
	if(status != NC_NOERR)
		return status;

	status = ncx_get_size_t(&xp, &nrecs);

	(void) ncp->nciop->rel(ncp->nciop, NC_NUMRECS_OFFSET, 0);

	if(status == NC_NOERR)
	{
		NC_set_numrecs(ncp, nrecs);
		fClr(ncp->flags, NC_NDIRTY);
	}

	return status;
}

#endif

/*
 * Write out just the numrecs member.
 * (A relatively expensive way to do things.)
 */
#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */ 

/*
 * Collective operation implicit
 */

int
write_numrecs(NC *ncp) {
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
  status = ncx_put_size_t(&pos, &nrecs);

  if(NC_indep(ncp) && NC_independentFhOpened(ncp->nciop)) {
    MPI_File_sync(ncp->nciop->independent_fh);
    MPI_Barrier(comm);
  }

  /* reset the file view */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, 
		    MPI_BYTE, "native", ncp->nciop->mpiinfo);

  if (rank == 0) {
    mpireturn = MPI_File_write_at(ncp->nciop->collective_fh, NC_NUMRECS_OFFSET,
				  buf, X_SIZEOF_SIZE_T, MPI_BYTE, &mpistatus); 
  } 

  MPI_Bcast(&mpireturn, 1, MPI_INT, 0, comm);

  if (mpireturn != MPI_SUCCESS) {
    status = NC_EWRITE;
  } else {
    MPI_File_sync(ncp->nciop->collective_fh); 
    MPI_Barrier(comm);
    fClr(ncp->flags, NC_NDIRTY);  
  }

  free(buf);

  return status;
}

#else /* Following interface modified as above by Jianwei Li */ 

int
write_numrecs(NC *ncp)
{
	int status = NC_NOERR;
	void *xp;

	assert(!NC_readonly(ncp));
	assert(!NC_indef(ncp));

	status = ncp->nciop->get(ncp->nciop,
		 NC_NUMRECS_OFFSET, NC_NUMRECS_EXTENT, RGN_WRITE, &xp);
	if(status != NC_NOERR)
		return status;

	{
		const size_t nrecs = NC_get_numrecs(ncp);
		status = ncx_put_size_t(&xp, &nrecs);
	}

	(void) ncp->nciop->rel(ncp->nciop, NC_NUMRECS_OFFSET, RGN_MODIFIED);

	if(status == NC_NOERR)
		fClr(ncp->flags, NC_NDIRTY);

	return status;
}

#endif /* USE_MPIO */

/*
 * Read in the header
 * It is expensive.
 */
#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */

int
read_NC(NC *ncp) {
  int status = NC_NOERR;

  free_NC_dimarrayV(&ncp->dims);
  free_NC_attrarrayV(&ncp->attrs);
  free_NC_vararrayV(&ncp->vars); 

  status = hdr_get_NC(ncp);

  if(status == NC_NOERR)
    fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);

  return status;
}

#else /* Following interface modified as above by Jianwei Li */

/* static */
int
read_NC(NC *ncp)
{
	int status = NC_NOERR;

	free_NC_dimarrayV(&ncp->dims);
	free_NC_attrarrayV(&ncp->attrs);
	free_NC_vararrayV(&ncp->vars);

	status = nc_get_NC(ncp);

	if(status == NC_NOERR)
		fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);

	return status;
}

#endif /* USE_MPIO */

/*
 * Write out the header
 */
#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */

int
write_NC(NC *ncp)
{
  int status = NC_NOERR, mpireturn;
  void *buf;
  MPI_Status mpistatus;
  int rank;
 
  assert(!NC_readonly(ncp));
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

  buf = (void *)malloc(ncp->xsz);
  status = hdr_put_NC(ncp, buf);
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
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE,
                    MPI_BYTE, "native", ncp->nciop->mpiinfo);

  if (rank == 0) {
    mpireturn = MPI_File_write_at(ncp->nciop->collective_fh, 0, buf, 
			          ncp->xsz, MPI_BYTE, &mpistatus);
  }

  fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);
  free(buf);
 
  return status;
} 

#else /* Following interface modified as above by Jianwei Li */

static int
write_NC(NC *ncp)
{
	int status = NC_NOERR;

	assert(!NC_readonly(ncp));

	status = ncx_put_NC(ncp, NULL, 0, 0);

	if(status == NC_NOERR)
		fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);

	return status;
}

#endif /* USE_MPIO */

/*
 * Write the header or the numrecs if necessary.
 */
int
NC_sync(NC *ncp)
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
    return write_numrecs(ncp);
  }
  /* else */

  return NC_NOERR;
}


#ifndef USE_MPIO /* out-of-date code removed by Jianwei Li */

/*
 * Initialize the 'non-record' variables.
 */
static int
fillerup(NC *ncp)
{
	int status = NC_NOERR;
	size_t ii;
	NC_var **varpp;

	assert(!NC_readonly(ncp));
	assert(NC_dofill(ncp));

	/* loop thru vars */
	varpp = ncp->vars.value;
	for(ii = 0; ii < ncp->vars.nelems; ii++, varpp++)
	{
		if(IS_RECVAR(*varpp))
		{
			/* skip record variables */
			continue;
		}

		status = fill_NC_var(ncp, *varpp, 0);
		if(status != NC_NOERR)
			break;
	}
	return status;
}

/* Begin endef */

/*
 */
static int
fill_added_recs(NC *gnu, NC *old)
{
	NC_var ** const gnu_varpp = (NC_var **)gnu->vars.value;

	const int old_nrecs = (int) NC_get_numrecs(old);
	int recno = 0;
	for(; recno < old_nrecs; recno++)
	{
		int varid = (int)old->vars.nelems;
		for(; varid < (int)gnu->vars.nelems; varid++)
		{
			const NC_var *const gnu_varp = *(gnu_varpp + varid);
			if(!IS_RECVAR(gnu_varp))
			{
				/* skip non-record variables */
				continue;
			}
			/* else */
			{
			const int status = fill_NC_var(gnu, gnu_varp, recno);
			if(status != NC_NOERR)
				return status;
			}
		}
	}

	return NC_NOERR;
}

/*
 */
static int
fill_added(NC *gnu, NC *old)
{
	NC_var ** const gnu_varpp = (NC_var **)gnu->vars.value;
	int varid = (int)old->vars.nelems;

	for(; varid < (int)gnu->vars.nelems; varid++)
	{
		const NC_var *const gnu_varp = *(gnu_varpp + varid);
		if(IS_RECVAR(gnu_varp))
		{
			/* skip record variables */
			continue;
		}
		/* else */
		{
		const int status = fill_NC_var(gnu, gnu_varp, 0);
		if(status != NC_NOERR)
			return status;
		}
	}

	return NC_NOERR;
}

#endif /* USE_MPIO */

/*
 * Move the records "out". 
 * Fill as needed.
 */
#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */

int
move_data_r(NC *ncp, NC *old) {
  /* no new variable inserted, move the whole contiguous data part */
  ncp->numrecs = old->numrecs;
  return ncio_move(ncp->nciop, ncp->begin_var, old->begin_var, 
              old->begin_rec - old->begin_var + old->recsize * old->numrecs);
}

int
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

    status = ncio_move(ncp->nciop, ncp_off, old_off, 
                       old_recsize * old_nrecs);
    if(status != NC_NOERR)
      return status;
  } else {

    /* else, new rec var inserted, to be moved one record at a time */

    for (recno = (int)old_nrecs -1; recno >= 0; recno--) {
      status = ncio_move(ncp->nciop, 
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

#else /* Following interface modified as above by Jianwei Li */

/* static */
int
move_recs_r(NC *gnu, NC *old)
{
	int status;
	int recno;
	int varid;
	NC_var **gnu_varpp = (NC_var **)gnu->vars.value;
	NC_var **old_varpp = (NC_var **)old->vars.value;
	NC_var *gnu_varp;
	NC_var *old_varp;
	off_t gnu_off;
	off_t old_off;
	const size_t old_nrecs = NC_get_numrecs(old);
	
	/* Don't parallelize this loop */
	for(recno = (int)old_nrecs -1; recno >= 0; recno--)
	{
	/* Don't parallelize this loop */
	for(varid = (int)old->vars.nelems -1; varid >= 0; varid--)
	{
		gnu_varp = *(gnu_varpp + varid);
		if(!IS_RECVAR(gnu_varp))
		{
			/* skip non-record variables on this pass */
			continue;
		}
		/* else */

		/* else, a pre-existing variable */
		old_varp = *(old_varpp + varid);
		gnu_off = gnu_varp->begin + (off_t)(gnu->recsize * recno);
		old_off = old_varp->begin + (off_t)(old->recsize * recno);

		if(gnu_off == old_off)
			continue; 	/* nothing to do */

		assert(gnu_off > old_off);
	
		status = gnu->nciop->move(gnu->nciop, gnu_off, old_off,
			 old_varp->len, 0);

		if(status != NC_NOERR)
			return status;
		
	}
	}

	NC_set_numrecs(gnu, old_nrecs);

	return NC_NOERR;
}

#endif /* USE_MPIO */

/*
 * Move the "non record" variables "out". 
 * Fill as needed.
 */
#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */

int
move_vars_r(NC *ncp, NC *old) {
  return ncio_move(ncp->nciop, ncp->begin_var, old->begin_var, 
                   old->begin_rec - old->begin_var); 
}
 
#else /* Following interface modified as above by Jianwei Li */   

/* static */
int
move_vars_r(NC *gnu, NC *old)
{
	int status;
	int varid;
	NC_var **gnu_varpp = (NC_var **)gnu->vars.value;
	NC_var **old_varpp = (NC_var **)old->vars.value;
	NC_var *gnu_varp;
	NC_var *old_varp;
	off_t gnu_off;
	off_t old_off;
	
	/* Don't parallelize this loop */
	for(varid = (int)old->vars.nelems -1;
		 varid >= 0; varid--)
	{
		gnu_varp = *(gnu_varpp + varid);
		if(IS_RECVAR(gnu_varp))
		{
			/* skip record variables on this pass */
			continue;
		}
		/* else */

		old_varp = *(old_varpp + varid);
		gnu_off = gnu_varp->begin;
		old_off = old_varp->begin;
	
		if(gnu_off == old_off)
			continue; 	/* nothing to do */

		assert(gnu_off > old_off);

		status = gnu->nciop->move(gnu->nciop, gnu_off, old_off,
			 old_varp->len, 0);

		if(status != NC_NOERR)
			return status;
		
	}

	return NC_NOERR;
}

#endif /* USE_MPIO */

#ifdef USE_MPIO /* function below added by Jianwei Li */

int 
NC_enddef(NC *ncp) {
  int status = NC_NOERR;
  MPI_Comm comm;

  assert(!NC_readonly(ncp));
  assert(NC_indef(ncp)); 

  comm = ncp->nciop->comm;

  NC_begins(ncp, 0, 1, 0, 1);
 
  /* To be updated */
  if(ncp->old != NULL) {
    /* a plain redef, not a create */
    assert(!NC_IsNew(ncp));
    assert(fIsSet(ncp->flags, NC_INDEF));
    assert(ncp->begin_rec >= ncp->old->begin_rec);
    assert(ncp->begin_var >= ncp->old->begin_var);
    assert(ncp->vars.nelems >= ncp->old->vars.nelems);
 
    MPI_File_sync(ncp->nciop->collective_fh);
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
      } else {
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
  }
 
  status = write_NC(ncp);
  if (status != NC_NOERR)
    return status;
 
  if(ncp->old != NULL) {
    free_NC(ncp->old);
    ncp->old = NULL;
  }
 
  fClr(ncp->flags, NC_CREAT | NC_INDEF);
  MPI_File_sync(ncp->nciop->collective_fh);
  MPI_Barrier(comm);

  return NC_NOERR;
}


int 
enddef(NC *ncp)
{
  assert(!NC_readonly(ncp));
  if(!NC_indef(ncp))
    return(NC_ENOTINDEFINE);

  NC_begins(ncp, 0, 1, 0, 1);

  if(ncp->old != NULL)
  {
    free_NC(ncp->old);
    ncp->old = NULL;
  }

  fClr(ncp->flags, NC_CREAT | NC_INDEF);

  return NC_NOERR;
}

#endif

/*
 *  End define mode.
 *  Common code for ncendef, ncclose(endef)
 *  Flushes I/O buffers.
 */
static int
NC_endef(NC *ncp,
	size_t h_minfree, size_t v_align,
	size_t v_minfree, size_t r_align)
{
	int status = NC_NOERR;

	assert(!NC_readonly(ncp));
	assert(NC_indef(ncp));

	NC_begins(ncp, h_minfree, v_align, v_minfree, r_align);

	if(ncp->old != NULL)
	{
		/* a plain redef, not a create */
		assert(!NC_IsNew(ncp));
		assert(fIsSet(ncp->flags, NC_INDEF));
		assert(ncp->begin_rec >= ncp->old->begin_rec);
		assert(ncp->begin_var >= ncp->old->begin_var);

		if(ncp->vars.nelems != 0)
		{
		if(ncp->begin_rec > ncp->old->begin_rec)
		{
			status = move_recs_r(ncp, ncp->old);
			if(status != NC_NOERR)
				return status;
			if(ncp->begin_var > ncp->old->begin_var)
			{
				status = move_vars_r(ncp, ncp->old);
				if(status != NC_NOERR)
					return status;
			} 
			/* else if (ncp->begin_var == ncp->old->begin_var) { NOOP } */
		}
		else
		{	/* Even if (ncp->begin_rec == ncp->old->begin_rec)
			   and     (ncp->begin_var == ncp->old->begin_var)
			   might still have added a new record variable */
		        if(ncp->recsize > ncp->old->recsize)
			{
			        status = move_recs_r(ncp, ncp->old);
				if(status != NC_NOERR)
				      return status;
			}
		}
		}
	}

	status = write_NC(ncp);
	if(status != NC_NOERR)
		return status;

#ifndef USE_MPIO /* out-of-date code removed by Jianwei Li */

	if(NC_dofill(ncp))
	{
		if(NC_IsNew(ncp))
		{
			status = fillerup(ncp);
			if(status != NC_NOERR)
				return status;
			
		}
		else if(ncp->vars.nelems > ncp->old->vars.nelems)
		{
			status = fill_added(ncp, ncp->old);
			if(status != NC_NOERR)
				return status;
			status = fill_added_recs(ncp, ncp->old);
			if(status != NC_NOERR)
				return status;
		}
	}

#endif /* USE_MPIO */

	if(ncp->old != NULL)
	{
		free_NC(ncp->old);
		ncp->old = NULL;
	}

	fClr(ncp->flags, NC_CREAT | NC_INDEF);

	return ncp->nciop->sync(ncp->nciop);
}

#ifdef LOCKNUMREC
static int
NC_init_pe(NC *ncp, int basepe) {
	if (basepe < 0 || basepe >= _num_pes()) {
		return NC_EINVAL; /* invalid base pe */
	}
	/* initialize common values */
	ncp->lock[LOCKNUMREC_VALUE] = 0;
	ncp->lock[LOCKNUMREC_LOCK] = 0;
	ncp->lock[LOCKNUMREC_SERVING] = 0;
	ncp->lock[LOCKNUMREC_BASEPE] =  basepe;
	return NC_NOERR;
}
#endif

/* Public */

#ifdef USE_MPIO /* This section is added by Jianwei Li */

int 
NC_close(NC *ncp) {
  int status = NC_NOERR;

  if(NC_indef(ncp)) {
    status = NC_enddef(ncp); /* TODO: defaults */
    if(status != NC_NOERR ) {
      /* To do: Abort new definition, if any */
      if (ncp->old != NULL) {
        free_NC(ncp->old);
        ncp->old = NULL;
        fClr(ncp->flags, NC_INDEF);
      }
    }
  }
  else if(!NC_readonly(ncp)) {
    status = NC_sync(ncp);
    if (status != NC_NOERR)
      return status;
  }
 
  (void) ncio_close(ncp->nciop, 0);
  ncp->nciop = NULL;
 
  del_from_NCList(ncp);
 
  free_NC(ncp);
 
  return status;
}

#else

int
nc_close(int ncid)
{
	int status = NC_NOERR;
	NC *ncp; 

	status = NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(NC_indef(ncp))
	{
		status = NC_endef(ncp, 0, 1, 0, 1); /* TODO: defaults */
		if(status != NC_NOERR )
		{
			(void) nc_abort(ncid);
			return status;
		}
	}
	else if(!NC_readonly(ncp))
	{
		status = NC_sync(ncp);
	}

	(void) ncio_close(ncp->nciop, 0);
	ncp->nciop = NULL;

	del_from_NCList(ncp);

	free_NC(ncp);

	return status;
}
#endif /* #ifdef USE_MPI */

int
nc_inq(int ncid,
	int *ndimsp,
	int *nvarsp,
	int *nattsp,
	int *xtendimp)
{
	int status;
	NC *ncp;

	status = NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(ndimsp != NULL)
		*ndimsp = (int) ncp->dims.nelems;
	if(nvarsp != NULL)
		*nvarsp = (int) ncp->vars.nelems;
	if(nattsp != NULL)
		*nattsp = (int) ncp->attrs.nelems;
	if(xtendimp != NULL)
		*xtendimp = find_NC_Udim(&ncp->dims, NULL);

	return NC_NOERR;
}

int 
nc_inq_ndims(int ncid, int *ndimsp)
{
	int status;
	NC *ncp;

	status = NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(ndimsp != NULL)
		*ndimsp = (int) ncp->dims.nelems;

	return NC_NOERR;
}

int 
nc_inq_nvars(int ncid, int *nvarsp)
{
	int status;
	NC *ncp;

	status = NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(nvarsp != NULL)
		*nvarsp = (int) ncp->vars.nelems;

	return NC_NOERR;
}

int 
nc_inq_natts(int ncid, int *nattsp)
{
	int status;
	NC *ncp;

	status = NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(nattsp != NULL)
		*nattsp = (int) ncp->attrs.nelems;

	return NC_NOERR;
}

int 
nc_inq_unlimdim(int ncid, int *xtendimp)
{
	int status;
	NC *ncp;

	status = NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(xtendimp != NULL)
		*xtendimp = find_NC_Udim(&ncp->dims, NULL);

	return NC_NOERR;
}


int
nc_sync(int ncid)
{
	int status;
	NC *ncp;

	status = NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(NC_indef(ncp))
		return NC_EINDEFINE;

	if(NC_readonly(ncp))
	{
		return read_NC(ncp);
	}
	/* else, read/write */

	status = NC_sync(ncp);
	if(status != NC_NOERR)
		return status;

	return ncp->nciop->sync(ncp->nciop);
}


int
nc_set_fill(int ncid,
	int fillmode, int *old_mode_ptr)
{
	int status;
	NC *ncp;
	int oldmode;

	status = NC_check_id(ncid, &ncp); 
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
			status = NC_sync(ncp);
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

#ifdef LOCKNUMREC

/* create function versions of the NC_*_numrecs macros */
size_t NC_get_numrecs(const NC *ncp) {
	shmem_t numrec;
	shmem_short_get(&numrec, (shmem_t *) ncp->lock + LOCKNUMREC_VALUE, 1,
		ncp->lock[LOCKNUMREC_BASEPE]);
	return (size_t) numrec;
}

void NC_set_numrecs(NC *ncp, size_t nrecs) {
	shmem_t numrec = (shmem_t) nrecs;
	/* update local value too */
	ncp->lock[LOCKNUMREC_VALUE] = (ushmem_t) numrec;
	shmem_short_put((shmem_t *) ncp->lock + LOCKNUMREC_VALUE, &numrec, 1,
		ncp->lock[LOCKNUMREC_BASEPE]);
}

void NC_increase_numrecs(NC *ncp, size_t nrecs) {
	/* this is only called in one place that's already protected
	 * by a lock ... so don't worry about it */
	if (nrecs > NC_get_numrecs(ncp))
		NC_set_numrecs(ncp, nrecs);
}

#endif /* LOCKNUMREC */

/* everyone in communicator group will be executing this */
/*ARGSUSED*/
int
nc_set_base_pe(int ncid, int pe)
{
#if _CRAYMPP && defined(LOCKNUMREC)
	int status;
	NC *ncp;
	shmem_t numrecs;

	if ((status = NC_check_id(ncid, &ncp)) != NC_NOERR) {
		return status;
	}
	if (pe < 0 || pe >= _num_pes()) {
		return NC_EINVAL; /* invalid base pe */
	}

	numrecs = (shmem_t) NC_get_numrecs(ncp);

	ncp->lock[LOCKNUMREC_VALUE] = (ushmem_t) numrecs;

	/* update serving & lock values for a "smooth" transition */
	/* note that the "real" server will being doing this as well */
	/* as all the rest in the group */
	/* must have syncronization before & after this step */
	shmem_short_get(
		(shmem_t *) ncp->lock + LOCKNUMREC_SERVING,
		(shmem_t *) ncp->lock + LOCKNUMREC_SERVING,
		1, ncp->lock[LOCKNUMREC_BASEPE]);

	shmem_short_get(
		(shmem_t *) ncp->lock + LOCKNUMREC_LOCK,
		(shmem_t *) ncp->lock + LOCKNUMREC_LOCK,
		1, ncp->lock[LOCKNUMREC_BASEPE]);

	/* complete transition */
	ncp->lock[LOCKNUMREC_BASEPE] = (ushmem_t) pe;

#endif /* _CRAYMPP && LOCKNUMREC */
	return NC_NOERR;
}

/*ARGSUSED*/
int
nc_inq_base_pe(int ncid, int *pe)
{
#if _CRAYMPP && defined(LOCKNUMREC)
	int status;
	NC *ncp;

	if ((status = NC_check_id(ncid, &ncp)) != NC_NOERR) {
		return status;
	}

	*pe = (int) ncp->lock[LOCKNUMREC_BASEPE];
#else
	/*
	 * !_CRAYMPP, only pe 0 is valid
	 */
	*pe = 0;
#endif /* _CRAYMPP && LOCKNUMREC */
	return NC_NOERR;
}
