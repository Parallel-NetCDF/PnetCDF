/*
 *	Copyright 1996, University Corporation for Atmospheric Research
 *      See netcdf/COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id$ */

#include "nc.h"
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>
#include <assert.h>
#include "ncx.h"
#include "fbits.h"

/*
 * Free dim
 * Formerly
NC_free_dim(dim)
 */
void
ncmpii_free_NC_dim(NC_dim *dimp)
{
	if(dimp == NULL)
		return;
	ncmpii_free_NC_string(dimp->name);
	free(dimp);
}


NC_dim *
ncmpii_new_x_NC_dim(NC_string *name)
{
	NC_dim *dimp;

	dimp = (NC_dim *) malloc(sizeof(NC_dim));
	if(dimp == NULL)
		return NULL;

	dimp->name = name;
	dimp->size = 0;

	return(dimp);
}


/*
 * Formerly
NC_new_dim(const char *name, long size)
 */
static NC_dim *
ncmpii_new_NC_dim(const char *name, MPI_Offset size)
{
	NC_string *strp;
	NC_dim *dimp;

	strp = ncmpii_new_NC_string(strlen(name), name);
	if(strp == NULL)
		return NULL;

	dimp = ncmpii_new_x_NC_dim(strp);
	if(dimp == NULL)
	{
		ncmpii_free_NC_string(strp);
		return NULL;
	}

	dimp->size = size;

	return(dimp);
}


static NC_dim *
dup_NC_dim(const NC_dim *dimp)
{
	return ncmpii_new_NC_dim(dimp->name->cp, dimp->size);
}

/*
 * Step thru NC_DIMENSION array, seeking the UNLIMITED dimension.
 * Return dimid or -1 on not found.
 * *dimpp is set to the appropriate NC_dim.
 * The loop structure is odd. In order to parallelize,
 * we moved a clearer 'break' inside the loop body to the loop test.
 */
int
ncmpii_find_NC_Udim(const NC_dimarray *ncap, NC_dim **dimpp)
{
	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return -1;

	{
	int dimid = 0;
	NC_dim **loc = ncap->value;

	for(; (MPI_Offset) dimid < ncap->nelems
			 && (*loc)->size != NC_UNLIMITED; dimid++, loc++)
	{
		/*EMPTY*/
	}
	if(dimid >= ncap->nelems)
		return(-1); /* not found */
	/* else, normal return */
	if(dimpp != NULL)
		*dimpp = *loc;
	return dimid;
	}
}


/*
 * Step thru NC_DIMENSION array, seeking match on name.
 * Return dimid or -1 on not found.
 * *dimpp is set to the appropriate NC_dim.
 * The loop structure is odd. In order to parallelize,
 * we moved a clearer 'break' inside the loop body to the loop test.
 */
static int
NC_finddim(const NC_dimarray *ncap, const char *name, NC_dim **dimpp)
{

	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return -1;

	{
	MPI_Offset slen = strlen(name);
	int dimid = 0;
	NC_dim **loc = (NC_dim **) ncap->value;

	for(; (MPI_Offset) dimid < ncap->nelems
			&& (strlen((*loc)->name->cp) != slen
				|| strncmp((*loc)->name->cp, name, slen) != 0);
		 dimid++, loc++)
	{
		/*EMPTY*/
	}
	if(dimid >= ncap->nelems)
		return(-1); /* not found */
	/* else, normal return */
	if(dimpp != NULL)
			*dimpp = *loc;
	return(dimid);
	}
}


/* dimarray */


/*
 * Free the stuff "in" (referred to by) an NC_dimarray.
 * Leaves the array itself allocated.
 */
void
ncmpii_free_NC_dimarrayV0(NC_dimarray *ncap)
{
	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return;

	assert(ncap->value != NULL);

	{
		NC_dim **dpp = ncap->value;
		NC_dim *const *const end = &dpp[ncap->nelems];
		for( /*NADA*/; dpp < end; dpp++)
		{
			ncmpii_free_NC_dim(*dpp);
			*dpp = NULL;
		}
	}
	ncap->nelems = 0;
}


/*
 * Free NC_dimarray values.
 * formerly
NC_free_array()
 */
void
ncmpii_free_NC_dimarrayV(NC_dimarray *ncap)
{
	assert(ncap != NULL);
	
	if(ncap->nalloc == 0)
		return;

	assert(ncap->value != NULL);

	ncmpii_free_NC_dimarrayV0(ncap);

	free(ncap->value);
	ncap->value = NULL;
	ncap->nalloc = 0;
}


int
ncmpii_dup_NC_dimarrayV(NC_dimarray *ncap, const NC_dimarray *ref)
{
	int status = NC_NOERR;

	assert(ref != NULL);
	assert(ncap != NULL);

	if(ref->nelems != 0)
	{
		const MPI_Offset sz = ref->nelems * sizeof(NC_dim *);
		ncap->value = (NC_dim **) malloc(sz);
		if(ncap->value == NULL)
			return NC_ENOMEM;
		(void) memset(ncap->value, 0, sz);
		ncap->nalloc = ref->nelems;
	}

	ncap->nelems = 0;
	{
		NC_dim **dpp = ncap->value;
		const NC_dim **drpp = (const NC_dim **)ref->value;
		NC_dim *const *const end = &dpp[ref->nelems];
		for( /*NADA*/; dpp < end; drpp++, dpp++, ncap->nelems++)
		{
			*dpp = dup_NC_dim(*drpp);
			if(*dpp == NULL)
			{
				status = NC_ENOMEM;
				break;
			}
		}
	}

	if(status != NC_NOERR)
	{
		ncmpii_free_NC_dimarrayV(ncap);
		return status;
	}

	assert(ncap->nelems == ref->nelems);

	return NC_NOERR;
}


/*
 * Add a new handle on the end of an array of handles
 * Formerly
NC_incr_array(array, tail)
 */
static int
incr_NC_dimarray(NC_dimarray *ncap, NC_dim *newelemp)
{
	NC_dim **vp;

	assert(ncap != NULL);

	if(ncap->nalloc == 0)
	{
		assert(ncap->nelems == 0);
		vp = (NC_dim **) malloc(NC_ARRAY_GROWBY * sizeof(NC_dim *));
		if(vp == NULL)
			return NC_ENOMEM;
		ncap->value = vp;
		ncap->nalloc = NC_ARRAY_GROWBY;
	}
	else if(ncap->nelems +1 > ncap->nalloc)
	{
		vp = (NC_dim **) realloc(ncap->value,
			(ncap->nalloc + NC_ARRAY_GROWBY) * sizeof(NC_dim *));
		if(vp == NULL)
			return NC_ENOMEM;
		ncap->value = vp;
		ncap->nalloc += NC_ARRAY_GROWBY;
	}

	if(newelemp != NULL)
	{
		ncap->value[ncap->nelems] = newelemp;
		ncap->nelems++;
	}
	return NC_NOERR;
}


NC_dim *
ncmpii_elem_NC_dimarray(const NC_dimarray *ncap, size_t elem)
{
	assert(ncap != NULL);
		/* cast needed for braindead systems with signed MPI_Offset */
	if(ncap->nelems == 0 || (unsigned long long) elem >= ncap->nelems)
		return NULL;

	assert(ncap->value != NULL);

	return ncap->value[elem];
}


/* Public */

int
ncmpi_def_dim(int ncid, const char *name, MPI_Offset size, int *dimidp)
{
	int status;
	NC *ncp;
	int dimid;
	NC_dim *dimp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	if(!NC_indef(ncp))
		return NC_ENOTINDEFINE;

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	/* MPI_Offset is usually a signed value, but serial netcdf uses 
	 * MPI_Offset -- normally unsigned */
	if ((ncp->flags & NC_64BIT_OFFSET) && sizeof(off_t) > 4) {
		/* CDF2 format and LFS */
		if (size > X_UINT_MAX - 3 || (size < 0)) 
			/* "-3" handles rounded-up size */
			return NC_EDIMSIZE;
	} else if ((ncp->flags & NC_64BIT_DATA)) {
		/* CDF5 format*/
		if (size < 0) 
			/* "-3" handles rounded-up size */
			return NC_EDIMSIZE;
	} else {
		/* CDF1 format */
		if (size > X_INT_MAX - 3 || (size < 0))
			return NC_EDIMSIZE;
	}

	if(size == NC_UNLIMITED)
	{
		dimid = ncmpii_find_NC_Udim(&ncp->dims, &dimp);
		if(dimid != -1)
		{
			assert(dimid != -1);
			return NC_EUNLIMIT;
		}
	}

	if(ncp->dims.nelems >= NC_MAX_DIMS)
		return NC_EMAXDIMS;

	dimid = NC_finddim(&ncp->dims, name, &dimp);
	if(dimid != -1)
		return NC_ENAMEINUSE;
	
	dimp = ncmpii_new_NC_dim(name, size);
	if(dimp == NULL)
		return NC_ENOMEM;
	status = incr_NC_dimarray(&ncp->dims, dimp);
	if(status != NC_NOERR)
	{
		ncmpii_free_NC_dim(dimp);
		return status;
	}

	if(dimidp != NULL)
		*dimidp = (int)ncp->dims.nelems -1;
	return NC_NOERR;
}


int
ncmpi_inq_dimid(int ncid, const char *name, int *dimid_ptr)
{
	int status;
	NC *ncp;
	int dimid;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	dimid = NC_finddim(&ncp->dims, name, NULL);

	if(dimid == -1)
		return NC_EBADDIM;

	*dimid_ptr = dimid;
	return NC_NOERR;
}


int
ncmpi_inq_dim(int ncid, int dimid, char *name, MPI_Offset *sizep)
{
	int status;
	NC *ncp;
	NC_dim *dimp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	dimp = ncmpii_elem_NC_dimarray(&ncp->dims, (size_t) dimid);
	if(dimp == NULL)
		return NC_EBADDIM;

	if(name != NULL)
	{
		(void)strncpy(name, dimp->name->cp, 
			dimp->name->nchars);
		name[dimp->name->nchars] = 0;
	}
	if(sizep != 0)
	{
		if(dimp->size == NC_UNLIMITED)
			*sizep = NC_get_numrecs(ncp);
		else
			*sizep = dimp->size;	
	}
	return NC_NOERR;
}


int 
ncmpi_inq_dimname(int ncid, int dimid, char *name)
{
	int status;
	NC *ncp;
	NC_dim *dimp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	dimp = ncmpii_elem_NC_dimarray(&ncp->dims, (size_t) dimid);
	if(dimp == NULL)
		return NC_EBADDIM;

	if(name != NULL)
	{
		(void)strncpy(name, dimp->name->cp, 
			dimp->name->nchars);
		name[dimp->name->nchars] = 0;
	}

	return NC_NOERR;
}


int 
ncmpi_inq_dimlen(int ncid, int dimid, MPI_Offset *lenp)
{
	int status;
	NC *ncp;
	NC_dim *dimp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	dimp = ncmpii_elem_NC_dimarray(&ncp->dims, (size_t) dimid);
	if(dimp == NULL)
		return NC_EBADDIM;

	if(lenp != 0)
	{
		if(dimp->size == NC_UNLIMITED)
			*lenp = NC_get_numrecs(ncp);
		else
			*lenp = dimp->size;	
	}
	return NC_NOERR;
}


int
ncmpi_rename_dim( int ncid, int dimid, const char *newname)
{
    int status, existid;
    NC *ncp;
    NC_dim *dimp;

    status = ncmpii_NC_check_id(ncid, &ncp); 
    if (status != NC_NOERR)
        return status;

    if (NC_readonly(ncp))
        return NC_EPERM;

    status = ncmpii_NC_check_name(newname);
    if (status != NC_NOERR)
        return status;

    existid = NC_finddim(&ncp->dims, newname, &dimp);
    if (existid != -1)
        return NC_ENAMEINUSE;

    dimp = ncmpii_elem_NC_dimarray(&ncp->dims, (size_t) dimid);
    if (dimp == NULL)
        return NC_EBADDIM;

    if (NC_indef(ncp)) {
        NC_string *old = dimp->name;
        NC_string *newStr = ncmpii_new_NC_string(strlen(newname), newname);
        if (newStr == NULL)
            return NC_ENOMEM;
        dimp->name = newStr;
        ncmpii_free_NC_string(old);
        return NC_NOERR;
    }
    /* else, not in define mode */

    status = ncmpii_set_NC_string(dimp->name, newname);
    if (status != NC_NOERR)
        return status;

    set_NC_hdirty(ncp);

    if (NC_doHsync(ncp)) {
        status = ncmpii_NC_sync(ncp, 1);
        if (status != NC_NOERR)
            return status;
    }

    return NC_NOERR;
}
