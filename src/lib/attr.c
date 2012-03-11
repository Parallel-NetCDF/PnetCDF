/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include <stdio.h>
#include "nc.h"
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>
#include <assert.h>
#include "nc.h"
#include "ncx.h"
#include "fbits.h"
#include "rnd.h"
#include "macro.h"


/*
 * Free attr
 * Formerly
NC_free_attr()
 */
void
ncmpii_free_NC_attr(NC_attr *attrp)
{

	if(attrp == NULL)
		return;
	ncmpii_free_NC_string(attrp->name);
	NCI_Free(attrp);
}


/*----< ncmpix_len_NC_attrV() >----------------------------------------------*/
/*
 * How much space will 'nelems' of 'type' take in
 *  external representation (as the values of an attribute)?
 */
static size_t
ncmpix_len_NC_attrV(nc_type type,
                    size_t  nelems)
{
    switch(type) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  return ncmpix_len_char(nelems);
        case NC_SHORT:
        case NC_USHORT: return ncmpix_len_short(nelems);
        case NC_INT:
        case NC_UINT:   return ncmpix_len_int(nelems);
        case NC_FLOAT:  return ncmpix_len_float(nelems);
        case NC_DOUBLE: return ncmpix_len_double(nelems);
        case NC_INT64:
        case NC_UINT64: return ncmpix_len_int64(nelems);
        case NC_STRING: return ncmpix_len_string(nelems);
        default: assert("ncmpix_len_NC_attr bad type" == 0);
    }
    return 0;
}


NC_attr *
ncmpii_new_x_NC_attr(
	NC_string *strp,
	nc_type type,
	MPI_Offset nelems)
{
	NC_attr *attrp;
	const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
	size_t sz = M_RNDUP(sizeof(NC_attr));

	assert(!(xsz == 0 && nelems != 0));

	sz += xsz;

	attrp = (NC_attr *) NCI_Malloc(sz);
	if(attrp == NULL )
		return NULL;

	attrp->xsz = xsz;

	attrp->name = strp;
	attrp->type = type;
	attrp->nelems = nelems;
	if(xsz != 0)
		attrp->xvalue = (char *)attrp + M_RNDUP(sizeof(NC_attr));
	else
		attrp->xvalue = NULL;

	return(attrp);
}


/*
 * Formerly
NC_new_attr(name,type,count,value)
 */
static NC_attr *
ncmpii_new_NC_attr(
	const char *name,
	nc_type type,
	MPI_Offset nelems)
{
	NC_string *strp;
	NC_attr *attrp;

	assert(name != NULL && *name != 0);

	strp = ncmpii_new_NC_string(strlen(name), name);
	if(strp == NULL)
		return NULL;
	
	attrp = ncmpii_new_x_NC_attr(strp, type, nelems);
	if(attrp == NULL)
	{
		ncmpii_free_NC_string(strp);
		return NULL;
	}

	return(attrp);
}


static NC_attr *
dup_NC_attr(const NC_attr *rattrp)
{
	NC_attr *attrp = ncmpii_new_NC_attr(rattrp->name->cp,
		 rattrp->type, rattrp->nelems);
	if(attrp == NULL)
		return NULL;
	(void) memcpy(attrp->xvalue, rattrp->xvalue, rattrp->xsz);
	return attrp;
}

/* attrarray */

/*
 * Free the stuff "in" (referred to by) an NC_attrarray.
 * Leaves the array itself allocated.
 */
void
ncmpii_free_NC_attrarrayV0(NC_attrarray *ncap)
{
	assert(ncap != NULL);

	if(ncap->ndefined == 0)
		return;

	assert(ncap->value != NULL);

	{
		NC_attr **app = ncap->value;
		NC_attr *const *const end = &app[ncap->ndefined];
		for( /*NADA*/; app < end; app++)
		{
			ncmpii_free_NC_attr(*app);
			*app = NULL;
		}
	}
	ncap->ndefined = 0;
}


/*
 * Free NC_attrarray values.
 * formerly
NC_free_array()
 */
void
ncmpii_free_NC_attrarrayV(NC_attrarray *ncap)
{
	assert(ncap != NULL);
	
	if(ncap->nalloc == 0)
		return;

	assert(ncap->value != NULL);

	ncmpii_free_NC_attrarrayV0(ncap);

	NCI_Free(ncap->value);
	ncap->value = NULL;
	ncap->nalloc = 0;
}


int
ncmpii_dup_NC_attrarrayV(NC_attrarray *ncap, const NC_attrarray *ref)
{
	int status = NC_NOERR;

	assert(ref != NULL);
	assert(ncap != NULL);

	if(ref->ndefined != 0)
	{
		const size_t sz = ref->ndefined * sizeof(NC_attr *);
		ncap->value = (NC_attr **) NCI_Malloc(sz);
		if(ncap->value == NULL)
			return NC_ENOMEM;

		(void) memset(ncap->value, 0, sz);
		ncap->nalloc = ref->ndefined;
	}

	ncap->ndefined = 0;
	{
		NC_attr **app = ncap->value;
		const NC_attr **drpp = (const NC_attr **)ref->value;
		NC_attr *const *const end = &app[ref->ndefined];
		for( /*NADA*/; app < end; drpp++, app++, ncap->ndefined++)
		{
			*app = dup_NC_attr(*drpp);
			if(*app == NULL)
			{
				status = NC_ENOMEM;
				break;
			}
		}
	}

	if(status != NC_NOERR)
	{
		ncmpii_free_NC_attrarrayV(ncap);
		return status;
	}

	assert(ncap->ndefined == ref->ndefined);

	return NC_NOERR;
}


/*
 * Add a new handle on the end of an array of handles
 * Formerly
NC_incr_array(array, tail)
 */
static int
incr_NC_attrarray(NC_attrarray *ncap, NC_attr *newelemp)
{
	NC_attr **vp;

	assert(ncap != NULL);

	if(ncap->nalloc == 0)
	{
		assert(ncap->ndefined == 0);
		vp = (NC_attr **) NCI_Malloc(NC_ARRAY_GROWBY * sizeof(NC_attr *));
		if(vp == NULL)
			return NC_ENOMEM;

		ncap->value = vp;
		ncap->nalloc = NC_ARRAY_GROWBY;
	}
	else if(ncap->ndefined +1 > ncap->nalloc)
	{
		vp = (NC_attr **) NCI_Realloc(ncap->value,
			(ncap->nalloc + NC_ARRAY_GROWBY) * sizeof(NC_attr *));
		if(vp == NULL)
			return NC_ENOMEM;
	
		ncap->value = vp;
		ncap->nalloc += NC_ARRAY_GROWBY;
	}

	if(newelemp != NULL)
	{
		ncap->value[ncap->ndefined] = newelemp;
		ncap->ndefined++;
	}
	return NC_NOERR;
}


static NC_attr *
elem_NC_attrarray(const NC_attrarray *ncap, MPI_Offset elem)
{
	assert(ncap != NULL);
	if((elem < 0) || ncap->ndefined == 0 || elem >= ncap->ndefined)
		return NULL;

	assert(ncap->value != NULL);

	return ncap->value[elem];
}

/* End attarray per se */

/*----< NC_attrarray0() >----------------------------------------------------*/
/*
 * Given ncp and varid, return ptr to array of attributes
 * else NULL on error. This is equivalent to validate varid.
 */
static NC_attrarray *
NC_attrarray0(NC  *ncp,
              int  varid)
{
    if (varid == NC_GLOBAL) /* Global attribute, attach to cdf */
        return &ncp->attrs;

    if (varid >= 0 && varid < ncp->vars.ndefined)
        return &ncp->vars.value[varid]->attrs;

    return NULL;
}


/*----< ncmpii_NC_findattr() >------------------------------------------------*/
/*
 * Step thru NC_ATTRIBUTE array, seeking match on name.
 *  return match or NULL if Not Found.
 */
static int
ncmpii_NC_findattr(const NC_attrarray *ncap,
                   const char         *name)
{
    int i;

    assert(ncap != NULL);

    // if (ncap->ndefined == 0) return NULL; /* none created yet */

    for (i=0; i<ncap->ndefined; i++) {
        if (strcmp(ncap->value[i]->name->cp, name) == 0)
            return i;
    }
    return -1;
}


/*
 * Look up by ncid, varid and name, return NULL if not found
 */
static int 
NC_lookupattr(int ncid,
	int varid,
	const char *name, /* attribute name */
	NC_attr **attrpp) /* modified on return */
{
	int indx, status;
	NC *ncp;
	NC_attrarray *ncap;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	indx = ncmpii_NC_findattr(ncap, name);
	if(indx == -1)
		return NC_ENOTATT;

	if(attrpp != NULL)
		*attrpp = ncap->value[indx];

	return NC_NOERR;
}

/* Public */

int
ncmpi_inq_attname(int ncid, int varid, int attnum, char *name)

{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	attrp = elem_NC_attrarray(ncap, attnum);
	if(attrp == NULL)
		return NC_ENOTATT;

	(void) strncpy(name, attrp->name->cp, attrp->name->nchars);
	name[attrp->name->nchars] = 0;

	return NC_NOERR;
}


int 
ncmpi_inq_attid(int ncid, int varid, const char *name, int *attnump)
{
	int indx, status;
	NC *ncp;
	NC_attrarray *ncap;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;
	

	indx = ncmpii_NC_findattr(ncap, name);
	if(indx == -1)
		return NC_ENOTATT;

	if(attnump != NULL)
		*attnump = indx;

	return NC_NOERR;
}

int 
ncmpi_inq_atttype(int ncid, int varid, const char *name, nc_type *datatypep)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(datatypep != NULL)
		*datatypep = attrp->type;

	return NC_NOERR;
}

int 
ncmpi_inq_attlen(int ncid, int varid, const char *name, MPI_Offset *lenp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(lenp != NULL)
		*lenp = attrp->nelems;

	return NC_NOERR;
}

int
ncmpi_inq_att(int ncid,
	int varid,
	const char *name, /* input, attribute name */
	nc_type *datatypep,
	MPI_Offset *lenp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(datatypep != NULL)
		*datatypep = attrp->type;
	if(lenp != NULL)
		*lenp = attrp->nelems;

	return NC_NOERR;
}


int
ncmpi_rename_att( int ncid, int varid, const char *name, const char *newname)
{
    int indx, file_ver, status;
    NC *ncp;
    NC_attrarray *ncap;
    NC_attr *attrp;
    NC_string *newStr, *old;

    /* sortof inline clone of NC_lookupattr() */
    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    if (NC_readonly(ncp))
        return NC_EPERM;

    ncap = NC_attrarray0(ncp, varid);
    if (ncap == NULL)
        return NC_ENOTVAR;

    file_ver = 1;
    if (fIsSet(ncp->flags, NC_64BIT_OFFSET))
        file_ver = 2;
    else if (fIsSet(ncp->flags, NC_64BIT_DATA))
        file_ver = 5;

    status = ncmpii_NC_check_name(newname, file_ver);
    if (status != NC_NOERR)
        return status;

    indx = ncmpii_NC_findattr(ncap, name);
    if (indx < 0)
        return NC_ENOTATT;
    attrp = ncap->value[indx];
    /* end inline clone NC_lookupattr() */

    if (ncmpii_NC_findattr(ncap, newname) >= 0)
        /* name in use */
        return NC_ENAMEINUSE;

    old = attrp->name;
    if (NC_indef(ncp)) {
        newStr = ncmpii_new_NC_string(strlen(newname), newname);
        if (newStr == NULL)
            return NC_ENOMEM;
        attrp->name = newStr;
        ncmpii_free_NC_string(old);
        return NC_NOERR;
    }
    /* else, no in define mode */
    status = ncmpii_set_NC_string(old, newname);
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


int
ncmpi_copy_att(int ncid_in, int varid_in, const char *name, int ncid_out, int ovarid)
{
	int indx, status;
	NC_attr *iattrp;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = NC_lookupattr(ncid_in, varid_in, name, &iattrp);
	if(status != NC_NOERR)
		return status;

	status = ncmpii_NC_check_id(ncid_out, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, ovarid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	indx = ncmpii_NC_findattr(ncap, name);
	if (indx >= 0) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			attrp = ncap->value[indx]; /* convenience */
	
			if(iattrp->xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = iattrp->xsz;
			attrp->type = iattrp->type;
			attrp->nelems = iattrp->nelems;

			(void) memcpy(attrp->xvalue, iattrp->xvalue,
				iattrp->xsz);
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				status = ncmpii_NC_sync(ncp, 1);
				if(status != NC_NOERR)
					return status;
			}

			return NC_NOERR;
		}
		/* else, redefine using existing array slot */
		old = ncap->value[indx];
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->ndefined >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	attrp = ncmpii_new_NC_attr(name, iattrp->type, iattrp->nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	(void) memcpy(attrp->xvalue, iattrp->xvalue,
		iattrp->xsz);

	if(indx >= 0)
	{
		assert(old != NULL);
		ncap->value[indx] = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		status = incr_NC_attrarray(ncap, attrp);
		if(status != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return status;
		}
	}

	return NC_NOERR;
}


int
ncmpi_del_att(int ncid, int varid, const char *name)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	int attrid;
	MPI_Offset slen;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(!NC_indef(ncp))
		return NC_ENOTINDEFINE;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

			/* sortof inline ncmpii_NC_findattr() */
	slen = strlen(name);

	attrpp = (NC_attr **) ncap->value;
	for(attrid = 0; (size_t) attrid < ncap->ndefined; attrid++, attrpp++)
	{
		if( slen == (*attrpp)->name->nchars &&
			strncmp(name, (*attrpp)->name->cp, slen) == 0)
		{
			old = *attrpp;
			break;
		}
	}
	if( (size_t) attrid == ncap->ndefined )
		return NC_ENOTATT;
			/* end inline NC_findattr() */

	/* shuffle down */
	for(attrid++; (size_t) attrid < ncap->ndefined; attrid++)
	{
		*attrpp = *(attrpp + 1);
		attrpp++;
	}
	*attrpp = NULL;
	/* decrement count */
	ncap->ndefined--;

	ncmpii_free_NC_attr(old);

	return NC_NOERR;
}

static nc_type longtype = (sizeof(long) == sizeof(int) ? NC_INT : NC_INT64);

/* ncmpi_pad_xxx APIs are only applicable for filetypes of size smaller
 * than 4 bytes
 */
#define PAD_GETN_FILETYPE(filetype)                                           \
static int                                                                    \
ncmpix_pad_getn_##filetype(const void **xpp,                                  \
                           MPI_Offset   nelems,                               \
                           void        *tp,                                   \
                           nc_type      buftype)                              \
{                                                                             \
    switch(buftype) {                                                         \
        case NC_CHAR:                                                         \
        case NC_BYTE:                                                         \
            return ncmpix_pad_getn_##filetype##_schar (xpp, nelems, tp);      \
        case NC_UBYTE:                                                        \
            return ncmpix_pad_getn_##filetype##_uchar (xpp, nelems, tp);      \
        case NC_SHORT:                                                        \
            return ncmpix_pad_getn_##filetype##_short (xpp, nelems, tp);      \
        case NC_USHORT:                                                       \
            return ncmpix_pad_getn_##filetype##_ushort(xpp, nelems, tp);      \
        case NC_INT:                                                          \
            return ncmpix_pad_getn_##filetype##_int   (xpp, nelems, tp);      \
        case NC_UINT:                                                         \
            return ncmpix_pad_getn_##filetype##_uint  (xpp, nelems, tp);      \
        case NC_FLOAT:                                                        \
            return ncmpix_pad_getn_##filetype##_float (xpp, nelems, tp);      \
        case NC_DOUBLE:                                                       \
            return ncmpix_pad_getn_##filetype##_double(xpp, nelems, tp);      \
        case NC_INT64:                                                        \
            return ncmpix_pad_getn_##filetype##_int64 (xpp, nelems, tp);      \
        case NC_UINT64:                                                       \
            return ncmpix_pad_getn_##filetype##_uint64(xpp, nelems, tp);      \
        default:                                                              \
            assert("ncmpix_pad_getn_" #filetype "invalid buffer type" == 0);  \
            return NC_EBADTYPE;                                               \
    }                                                                         \
}
/*----< ncmpix_pad_getn_schar() >--------------------------------------------*/
/*----< ncmpix_pad_getn_uchar() >--------------------------------------------*/
/*----< ncmpix_pad_getn_short() >--------------------------------------------*/
/*----< ncmpix_pad_getn_ushort() >-------------------------------------------*/
PAD_GETN_FILETYPE(schar)
PAD_GETN_FILETYPE(uchar)
PAD_GETN_FILETYPE(short)
PAD_GETN_FILETYPE(ushort)

#define GETN_FILETYPE(filetype)                                               \
static int                                                                    \
ncmpix_getn_##filetype(const void **xpp,                                      \
                       MPI_Offset   nelems,                                   \
                       void        *tp,                                       \
                       nc_type      buftype)                                  \
{                                                                             \
    switch(buftype) {                                                         \
        case NC_CHAR:                                                         \
        case NC_BYTE:                                                         \
            return ncmpix_getn_##filetype##_schar (xpp, nelems, tp);          \
        case NC_UBYTE:                                                        \
            return ncmpix_getn_##filetype##_uchar (xpp, nelems, tp);          \
        case NC_SHORT:                                                        \
            return ncmpix_getn_##filetype##_short (xpp, nelems, tp);          \
        case NC_USHORT:                                                       \
            return ncmpix_getn_##filetype##_ushort(xpp, nelems, tp);          \
        case NC_INT:                                                          \
            return ncmpix_getn_##filetype##_int   (xpp, nelems, tp);          \
        case NC_UINT:                                                         \
            return ncmpix_getn_##filetype##_uint  (xpp, nelems, tp);          \
        case NC_FLOAT:                                                        \
            return ncmpix_getn_##filetype##_float (xpp, nelems, tp);          \
        case NC_DOUBLE:                                                       \
            return ncmpix_getn_##filetype##_double(xpp, nelems, tp);          \
        case NC_INT64:                                                        \
            return ncmpix_getn_##filetype##_int64 (xpp, nelems, tp);          \
        case NC_UINT64:                                                       \
            return ncmpix_getn_##filetype##_uint64(xpp, nelems, tp);          \
        default:                                                              \
            assert("ncmpix_pad_getn_" #filetype "invalid buffer type" == 0);  \
            return NC_EBADTYPE;                                               \
    }                                                                         \
}

/*----< ncmpix_getn_int() >--------------------------------------------------*/
/*----< ncmpix_getn_uint() >-------------------------------------------------*/
/*----< ncmpix_getn_float() >------------------------------------------------*/
/*----< ncmpix_getn_double() >-----------------------------------------------*/
/*----< ncmpix_getn_int64() >------------------------------------------------*/
/*----< ncmpix_getn_uint64() >-----------------------------------------------*/
GETN_FILETYPE(int)
GETN_FILETYPE(uint)
GETN_FILETYPE(float)
GETN_FILETYPE(double)
GETN_FILETYPE(int64)
GETN_FILETYPE(uint64)

/*----< ncmpix_pad_getn() >--------------------------------------------------*/
/* padding only applicable to file types of size smaller than 4 bytes */
static int
ncmpix_pad_getn(const void **xpp,
                MPI_Offset   nelems,
                void        *tp,
                nc_type      filetype,
                nc_type      buftype)
{
    /* get n elements from (filetype*)*xpp to (buftype*)tp */
    /* Checking for character-number conversion should have already been done */

    switch(filetype) {
        case NC_CHAR:
        case NC_BYTE:
            return ncmpix_pad_getn_schar (xpp, nelems, tp, buftype);
        case NC_UBYTE:
            return ncmpix_pad_getn_uchar (xpp, nelems, tp, buftype);
        case NC_SHORT:
            return ncmpix_pad_getn_short (xpp, nelems, tp, buftype);
        case NC_USHORT:
            return ncmpix_pad_getn_ushort(xpp, nelems, tp, buftype);
        case NC_INT:
            return ncmpix_getn_int       (xpp, nelems, tp, buftype);
        case NC_UINT:
            return ncmpix_getn_uint      (xpp, nelems, tp, buftype);
        case NC_FLOAT:
            return ncmpix_getn_float     (xpp, nelems, tp, buftype);
        case NC_DOUBLE:
            return ncmpix_getn_double    (xpp, nelems, tp, buftype);
        case NC_INT64:
            return ncmpix_getn_int64     (xpp, nelems, tp, buftype);
        case NC_UINT64:
            return ncmpix_getn_uint64    (xpp, nelems, tp, buftype);
        default: 
            assert("ncmpix_pad_getn invalid filetype" == 0);
            return NC_EBADTYPE;
    }
}

/*----< ncmpii_get_att() >---------------------------------------------------*/
static int
ncmpii_get_att(int         ncid,
               int         varid,
               const char *name,
               void       *tp,       /* I/O buffer */
               nc_type     buftype)  /* I/O buffer data type */
{
    int status;
    NC_attr *attrp;

    status = NC_lookupattr(ncid, varid, name, &attrp);
    if (status != NC_NOERR) return status;

    if (attrp->nelems == 0) return NC_NOERR;

    /* No character conversions are allowed. */
    if (attrp->type != buftype &&
        (attrp->type == NC_CHAR   || buftype == NC_CHAR ||
         attrp->type == NC_STRING || buftype == NC_STRING))
        return NC_ECHAR;

    const void *xp = attrp->xvalue;
    return ncmpix_pad_getn(&xp, attrp->nelems, tp, attrp->type, buftype);
}

/*----< ncmpi_get_att_text() >-----------------------------------------------*/
int
ncmpi_get_att_text(int ncid, int varid, const char *name, char *value)
{
    return ncmpii_get_att(ncid, varid, name, value, NC_CHAR);
}

#define GET_ATT_TYPE(fntype, ext_buftype, nc_buftype)                         \
int                                                                           \
ncmpi_get_att_##fntype(int ncid, int varid, const char  *name,                \
                       ext_buftype *value)                                    \
{                                                                             \
    return ncmpii_get_att(ncid, varid, name, value, nc_buftype);              \
}
/*----< ncmpi_get_att_schar() >----------------------------------------------*/
/*----< ncmpi_get_att_uchar() >----------------------------------------------*/
/*----< ncmpi_get_att_ubyte() >----------------------------------------------*/
/*----< ncmpi_get_att_short() >----------------------------------------------*/
/*----< ncmpi_get_att_ushort() >---------------------------------------------*/
/*----< ncmpi_get_att_int() >------------------------------------------------*/
/*----< ncmpi_get_att_uint() >-----------------------------------------------*/
/*----< ncmpi_get_att_long() >-----------------------------------------------*/
/*----< ncmpi_get_att_float() >----------------------------------------------*/
/*----< ncmpi_get_att_double() >---------------------------------------------*/
/*----< ncmpi_get_att_longlong() >-------------------------------------------*/
/*----< ncmpi_get_att_ulonglong() >------------------------------------------*/
/*----< ncmpi_get_att_string() >---------------------------------------------*/
GET_ATT_TYPE(schar,     signed char,        NC_BYTE)
GET_ATT_TYPE(uchar,     unsigned char,      NC_UBYTE)
GET_ATT_TYPE(ubyte,     unsigned char,      NC_UBYTE)
GET_ATT_TYPE(short,     short,              NC_SHORT)
GET_ATT_TYPE(ushort,    unsigned short,     NC_USHORT)
GET_ATT_TYPE(int,       int,                NC_INT)
GET_ATT_TYPE(uint,      unsigned int,       NC_UINT)
GET_ATT_TYPE(long,      long,               longtype)
GET_ATT_TYPE(float,     float,              NC_FLOAT)
GET_ATT_TYPE(double,    double,             NC_DOUBLE)
GET_ATT_TYPE(longlong,  long long,          NC_INT64)
GET_ATT_TYPE(ulonglong, unsigned long long, NC_UINT64)
// GET_ATT_TYPE(string, char*,              NC_STRING)

int
ncmpi_get_att_string(int ncid, int varid, const char  *name, char **value)
{
    printf("Error: string type is not yet supported\n");
    return NC_ENOTSUPPORT;
}


#define PAD_PUTN_FILETYPE(ftype)                                              \
static int                                                                    \
ncmpix_pad_putn_##ftype(void       **xpp,                                     \
                        MPI_Offset   nelems,                                  \
                        const void  *tp,                                      \
                        nc_type      btype)                                   \
{                                                                             \
    switch(btype) {                                                           \
        case NC_CHAR:                                                         \
        case NC_BYTE:                                                         \
            return ncmpix_pad_putn_##ftype##_schar (xpp, nelems, tp);         \
        case NC_UBYTE:                                                        \
            return ncmpix_pad_putn_##ftype##_uchar (xpp, nelems, tp);         \
        case NC_SHORT:                                                        \
            return ncmpix_pad_putn_##ftype##_short (xpp, nelems, tp);         \
        case NC_USHORT:                                                       \
            return ncmpix_pad_putn_##ftype##_ushort(xpp, nelems, tp);         \
        case NC_INT:                                                          \
            return ncmpix_pad_putn_##ftype##_int   (xpp, nelems, tp);         \
        case NC_UINT:                                                         \
            return ncmpix_pad_putn_##ftype##_uint  (xpp, nelems, tp);         \
        case NC_FLOAT:                                                        \
            return ncmpix_pad_putn_##ftype##_float (xpp, nelems, tp);         \
        case NC_DOUBLE:                                                       \
            return ncmpix_pad_putn_##ftype##_double(xpp, nelems, tp);         \
        case NC_INT64:                                                        \
            return ncmpix_pad_putn_##ftype##_int64 (xpp, nelems, tp);         \
        case NC_UINT64:                                                       \
            return ncmpix_pad_putn_##ftype##_uint64(xpp, nelems, tp);         \
        default:                                                              \
            assert("ncmpix_pad_putn_" #ftype "invalid type" == 0);            \
            return NC_EBADTYPE;                                               \
    }                                                                         \
}
/*----< ncmpix_pad_putn_schar() >--------------------------------------------*/
/*----< ncmpix_pad_putn_uchar() >--------------------------------------------*/
/*----< ncmpix_pad_putn_short() >--------------------------------------------*/
/*----< ncmpix_pad_putn_ushort() >-------------------------------------------*/
PAD_PUTN_FILETYPE(schar)
PAD_PUTN_FILETYPE(uchar)
PAD_PUTN_FILETYPE(short)
PAD_PUTN_FILETYPE(ushort)

#define PUTN_FILETYPE(ftype)                                                  \
static int                                                                    \
ncmpix_putn_##ftype(void           **xpp,                                     \
                    MPI_Offset       nelems,                                  \
                    const void *tp,                                           \
                    nc_type          btype)                                   \
{                                                                             \
    switch(btype) {                                                           \
        case NC_CHAR:                                                         \
        case NC_BYTE:                                                         \
            return ncmpix_putn_##ftype##_schar (xpp, nelems, tp);             \
        case NC_UBYTE:                                                        \
            return ncmpix_putn_##ftype##_uchar (xpp, nelems, tp);             \
        case NC_SHORT:                                                        \
            return ncmpix_putn_##ftype##_short (xpp, nelems, tp);             \
        case NC_USHORT:                                                       \
            return ncmpix_putn_##ftype##_ushort(xpp, nelems, tp);             \
        case NC_INT:                                                          \
            return ncmpix_putn_##ftype##_int   (xpp, nelems, tp);             \
        case NC_UINT:                                                         \
            return ncmpix_putn_##ftype##_uint  (xpp, nelems, tp);             \
        case NC_FLOAT:                                                        \
            return ncmpix_putn_##ftype##_float (xpp, nelems, tp);             \
        case NC_DOUBLE:                                                       \
            return ncmpix_putn_##ftype##_double(xpp, nelems, tp);             \
        case NC_INT64:                                                        \
            return ncmpix_putn_##ftype##_int64 (xpp, nelems, tp);             \
        case NC_UINT64:                                                       \
            return ncmpix_putn_##ftype##_uint64(xpp, nelems, tp);             \
        default:                                                              \
            assert("ncmpix_putn_" #ftype "invalid type" == 0);                \
            return NC_EBADTYPE;                                               \
    }                                                                         \
}
/*----< ncmpix_putn_int() >--------------------------------------------------*/
/*----< ncmpix_putn_uint() >-------------------------------------------------*/
/*----< ncmpix_putn_float() >------------------------------------------------*/
/*----< ncmpix_putn_double() >-----------------------------------------------*/
/*----< ncmpix_putn_int64() >------------------------------------------------*/
/*----< ncmpix_putn_uint64() >-----------------------------------------------*/
PUTN_FILETYPE(int)
PUTN_FILETYPE(uint)
PUTN_FILETYPE(float)
PUTN_FILETYPE(double)
PUTN_FILETYPE(int64)
PUTN_FILETYPE(uint64)

/*----< ncmpix_pad_putn() >--------------------------------------------------*/
/* padding only applicable to file types of size smaller than 4 bytes */
static int
ncmpix_pad_putn(void       **xpp,
                MPI_Offset   nelems,
                const void  *tp,
                nc_type      filetype,
                nc_type      buftype)
{
    /* put n elements from (buftype*)tp to (filetype*)*xpp */
    /* Checking for character-number conversion should have been done */

    switch(filetype) {
        case NC_CHAR:
        case NC_BYTE:
            return ncmpix_pad_putn_schar (xpp, nelems, tp, buftype);
        case NC_UBYTE:
            return ncmpix_pad_putn_uchar (xpp, nelems, tp, buftype);
        case NC_SHORT:
            return ncmpix_pad_putn_short (xpp, nelems, tp, buftype);
        case NC_USHORT:
            return ncmpix_pad_putn_ushort(xpp, nelems, tp, buftype);
        case NC_INT:
            return ncmpix_putn_int       (xpp, nelems, tp, buftype);
        case NC_UINT:
            return ncmpix_putn_uint      (xpp, nelems, tp, buftype);
        case NC_FLOAT:
            return ncmpix_putn_float     (xpp, nelems, tp, buftype);
        case NC_DOUBLE:
            return ncmpix_putn_double    (xpp, nelems, tp, buftype);
        case NC_INT64:
            return ncmpix_putn_int64     (xpp, nelems, tp, buftype);
        case NC_UINT64:
            return ncmpix_putn_uint64    (xpp, nelems, tp, buftype);
        default: 
            assert("ncmpix_pad_putn invalid filetype" == 0);
            return NC_EBADTYPE;
    }
}

/*----< ncmpii_put_att() >---------------------------------------------------*/
/* Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters
 */
static int
ncmpii_put_att(int         ncid,
               int         varid,
               const char *name,     /* attribute name */
               nc_type     filetype, /* type defined in file header */
               MPI_Offset  nelems,   /* number of elements of type buftype */
               const void *buf,      /* I/O buffer */
               nc_type     buftype)  /* I/O buffer type */
{
    int indx, file_ver, status;
    NC *ncp;
    NC_attrarray *ncap;
    NC_attr *attrp, *old=NULL;

   if (!name || strlen(name) > NC_MAX_NAME)
      return NC_EBADNAME;

    /* Should CDF-5 allow very large file header? */
    // if (len > X_INT_MAX) return NC_EINVAL;

    /* get the file ID */
    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    /* file should be opened with writable permission */
    if (NC_readonly(ncp)) return NC_EPERM;

    /* nelems can be zero, i.e. an attribute with only its name */
    if (nelems > 0 && buf == NULL)
        return NC_EINVAL; /* Null arg */

    /* get the file format version */
    file_ver = 1;
    if (fIsSet(ncp->flags, NC_64BIT_OFFSET))
        file_ver = 2;
    else if (fIsSet(ncp->flags, NC_64BIT_DATA))
        file_ver = 5;

    if (nelems < 0 || (nelems > X_INT_MAX && file_ver <= 2))
        return NC_EINVAL; /* Invalid nelems */

    /* No character conversions are allowed. */
    if (filetype != buftype &&
        (filetype == NC_CHAR   || buftype == NC_CHAR ||
         filetype == NC_STRING || buftype == NC_STRING))
        return NC_ECHAR;

    /* check if the attribute name is legal */
    status = ncmpii_NC_check_name(name, file_ver);
    if (status != NC_NOERR) return status;

    /* check if filetype is valid, as filetype is given by user
     * no need to check buftype, as buftype is set internally
     */
    status = ncmpii_cktype(filetype);
    if (status != NC_NOERR) return status;

    /* get the pointer to the attribute array */
    ncap = NC_attrarray0(ncp, varid);
    if (ncap == NULL) return NC_ENOTVAR;

    indx = ncmpii_NC_findattr(ncap, name);
    if (indx >= 0) { /* name in use */
        if (!NC_indef(ncp)) {
            /* not in define mode, meaning to over-write attribute value */

            const size_t xsz = ncmpix_len_NC_attrV(filetype, nelems);
            /* xsz is the total size of this attribute */

            attrp = ncap->value[indx]; /* convenience */

            if (xsz > attrp->xsz) /* new attribute requires a larger space */
                return NC_ENOTINDEFINE;
            /* else, we can reuse existing without redef */

            attrp->xsz = xsz;
            attrp->type = filetype;
            attrp->nelems = nelems;

            if (nelems != 0) {
                void *xp = attrp->xvalue;
                /* using xp to prevent change the pointer attr->xvalue,
                 * as ncmpix_pad_putn() advances the first argument
                 * with nelems elements
                 */
                status = ncmpix_pad_putn(&xp, nelems, buf, filetype, buftype);
                if (status != NC_NOERR) return status;
            }

            set_NC_hdirty(ncp); /* indicate the header is dirty */

            if (NC_doHsync(ncp)) {
                /* if NC_SHARE is used, we need to flush to file now */
                status = ncmpii_NC_sync(ncp, 1);
                if (status != NC_NOERR) return status;
            }

            return NC_NOERR;
        }
        /* else, redefine using existing array slot */
        old = ncap->value[indx];
    }
    else { /* name never been used */
        /* creating new attributes must be done in define mode */
        if (!NC_indef(ncp)) return NC_ENOTINDEFINE;

        if (ncap->ndefined >= NC_MAX_ATTRS)
            return NC_EMAXATTS;
    }

    /* create a new attribute object */
    attrp = ncmpii_new_NC_attr(name, filetype, nelems);
    if (attrp == NULL) return NC_ENOMEM;

    if (nelems != 0) { /* non-zero length attribute */
        void *xp = attrp->xvalue;
        status = ncmpix_pad_putn(&xp, nelems, buf, filetype, buftype);
        /* wkliao: no immediately return error code here? Strange ...
         *         instead, continue calling incr_NC_attrarray() then
         *         return the error code at the end
         */
    }

    if (indx >= 0) { /* modify the existing attribute */
        assert(old != NULL);
        ncap->value[indx] = attrp;
        ncmpii_free_NC_attr(old);
    }
    else { /* creating a new attribute */
        int lstatus = incr_NC_attrarray(ncap, attrp);
        if (lstatus != NC_NOERR) {
            ncmpii_free_NC_attr(attrp);
            return lstatus;
        }
    }

    return status;
}

/*----< ncmpi_put_att_text() >-----------------------------------------------*/
int
ncmpi_put_att_text(int ncid, int varid, const char  *name,
                   MPI_Offset nelems, const char *value)
{
    return ncmpii_put_att(ncid, varid, name, NC_CHAR,
                          nelems, value, NC_CHAR);
}

#define PUT_ATT_TYPE(fntype, ext_buftype, nc_buftype)                         \
int                                                                           \
ncmpi_put_att_##fntype(int ncid, int varid, const char  *name, nc_type xtype, \
                       MPI_Offset nelems, const ext_buftype *value)           \
{                                                                             \
    return ncmpii_put_att(ncid, varid, name, xtype, nelems, value,            \
                          nc_buftype);                                        \
}
/*----< ncmpi_put_att_schar() >----------------------------------------------*/
/*----< ncmpi_put_att_uchar() >----------------------------------------------*/
/*----< ncmpi_put_att_ubyte() >----------------------------------------------*/
/*----< ncmpi_put_att_short() >----------------------------------------------*/
/*----< ncmpi_put_att_ushort() >---------------------------------------------*/
/*----< ncmpi_put_att_int() >------------------------------------------------*/
/*----< ncmpi_put_att_uint() >-----------------------------------------------*/
/*----< ncmpi_put_att_long() >-----------------------------------------------*/
/*----< ncmpi_put_att_float() >----------------------------------------------*/
/*----< ncmpi_put_att_double() >---------------------------------------------*/
/*----< ncmpi_put_att_longlong() >-------------------------------------------*/
/*----< ncmpi_put_att_ulonglong() >------------------------------------------*/
/*----< ncmpi_put_att_string() >---------------------------------------------*/
PUT_ATT_TYPE(schar,     signed char,        NC_BYTE)
PUT_ATT_TYPE(uchar,     unsigned char,      NC_UBYTE)
PUT_ATT_TYPE(ubyte,     unsigned char,      NC_UBYTE)
PUT_ATT_TYPE(short,     short,              NC_SHORT)
PUT_ATT_TYPE(ushort,    unsigned short,     NC_USHORT)
PUT_ATT_TYPE(int,       int,                NC_INT)
PUT_ATT_TYPE(uint,      unsigned int,       NC_UINT)
PUT_ATT_TYPE(long,      long,               longtype)
PUT_ATT_TYPE(float,     float,              NC_FLOAT)
PUT_ATT_TYPE(double,    double,             NC_DOUBLE)
PUT_ATT_TYPE(longlong,  long long,          NC_INT64)
PUT_ATT_TYPE(ulonglong, unsigned long long, NC_UINT64)
// PUT_ATT_TYPE(string, char*,              NC_STRING)

int
ncmpi_put_att_string(int ncid, int varid, const char  *name,
                     MPI_Offset nelems, const char **value)
{
    printf("Error: string type is not yet supported\n");
    return NC_ENOTSUPPORT;
    // return ncmpii_put_att(ncid, varid, name, NC_STRING,
    //                       nelems, value, NC_STRING);
}

/* For netCDF, the type mapping between file types and buffer types
 * are based on netcdf4. Check APIs of nc_put_att_xxx from source files
 *     netCDF/netcdf-4.1.3/libdispatch/att.c
 *     netCDF/netcdf-4.1.3/libsrc4/nc4attr.c
 *
 * Note that schar means signed 1-byte integers in attributes. Hence the call
 * below is illegal. NC_ECHAR will return, indicating the error on trying
 * type conversion between characters and numbers.
 *
 * ncmpi_put_att_schar(ncid, varid, "attr name", NC_CHAR, strlen(attrp), attrp);
 *
 * This rule and mapping apply for variables as well. See APIs of
 * nc_put_vara_xxx from source files
 *     netCDF/netcdf-4.1.3/libdispatch/var.c
 *     netCDF/netcdf-4.1.3/libsrc4/nc4var.c
 *
 */
