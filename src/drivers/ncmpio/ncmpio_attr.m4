dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in
 * src/dispatchers/attribute.c and src/dispatchers/attr_getput.m4
 *
 * ncmpi_inq_att()     : dispatcher->inq_att()
 * ncmpi_inq_attid()   : dispatcher->inq_attid()
 * ncmpi_inq_attname() : dispatcher->inq_attname()
 * ncmpi_copy_att()    : dispatcher->copy_att()
 * ncmpi_rename_att()  : dispatcher->rename_att()
 * ncmpi_del_att()     : dispatcher->del_att()
 * ncmpi_put_att()     : dispatcher->put_att()
 * ncmpi_get_att()     : dispatcher->get_att()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncx.h>
#include "ncmpio_NC.h"

include(`foreach.m4')dnl
include(`utils.m4')dnl

/*----< x_len_NC_attrV() >---------------------------------------------------*/
/* How much space will 'nelems' of 'xtype' take in external representation.
 * Note the space is aligned in 4-byte boundary.
 */
static MPI_Offset
x_len_NC_attrV(nc_type    xtype,
               MPI_Offset nelems)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  return _RNDUP(nelems, 4);
        case NC_SHORT:
        case NC_USHORT: return ((nelems + nelems%2) * 2);
        case NC_INT:
        case NC_UINT:
        case NC_FLOAT:  return (nelems * 4);
        case NC_DOUBLE:
        case NC_INT64:
        case NC_UINT64: return (nelems * 8);
        default: fprintf(stderr, "Error: bad type(%d) in %s\n",xtype,__func__);
    }
    return 0;
}

/*----< ncmpio_new_NC_attr() >-----------------------------------------------*/
/*
 * IN:  name is an already normalized attribute name (NULL terminated)
 * OUT: (*attrp)->xvalue is malloc-ed
 */
int
ncmpio_new_NC_attr(char        *name,
                   size_t       name_len,
                   nc_type      xtype,
                   MPI_Offset   nelems,
                   NC_attr    **attrp)
{
    *attrp = (NC_attr*) NCI_Malloc(sizeof(NC_attr));
    if (*attrp == NULL ) DEBUG_RETURN_ERROR(NC_ENOMEM)

    (*attrp)->xtype    = xtype;
    (*attrp)->xsz      = 0;
    (*attrp)->nelems   = nelems;
    (*attrp)->xvalue   = NULL;
    (*attrp)->name     = name;
    (*attrp)->name_len = name_len;

    if (nelems > 0) {
        /* obtain 4-byte aligned size of space to store the values */
        MPI_Offset xsz = x_len_NC_attrV(xtype, nelems);
        (*attrp)->xsz    = xsz;
        (*attrp)->xvalue = NCI_Malloc((size_t)xsz);
        if ((*attrp)->xvalue == NULL) {
            NCI_Free(*attrp);
            *attrp = NULL;
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
    }
    return NC_NOERR;
}

/*----< ncmpio_free_NC_attr() >----------------------------------------------*/
void
ncmpio_free_NC_attr(NC_attr *attrp)
{
    assert(attrp != NULL);

    if (attrp->xvalue != NULL)
        NCI_Free(attrp->xvalue);
    NCI_Free(attrp->name);
}

/*----< dup_NC_attr() >------------------------------------------------------*/
/* duplicate an NC_attr object */
static int
dup_NC_attr(const NC_attr *rattrp, NC_attr **attrp)
{
    char *name;

    /* rattrp->name has already been normalized */
    name = (char*) NCI_Malloc(rattrp->name_len + 1);
    if (name == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    strncpy(name, rattrp->name, rattrp->name_len);
    name[rattrp->name_len] = '\0';

    return ncmpio_new_NC_attr(name, rattrp->name_len, rattrp->xtype,
                              rattrp->nelems, attrp);
}

/* attrarray */

/*----< ncmpio_free_NC_attrarray() >-----------------------------------------*/
/* Free NC_attrarray values. */
void
ncmpio_free_NC_attrarray(NC_attrarray *ncap)
{
    int i;

    assert(ncap != NULL);

    if (ncap->value != NULL) {
        /* when error is detected reading NC_ATTRIBUTE tag, ncap->ndefined can
         * be > 0 and ncap->value is still NULL
         */
        for (i=0; i<ncap->ndefined; i++) {
            if (ncap->value[i] == NULL) continue;
            ncmpio_free_NC_attr(ncap->value[i]);
            NCI_Free(ncap->value[i]);
        }

        /* attributes can be deleted, thus ncap->value can still be allocated
         * while ncap->ndefined == 0 */
        NCI_Free(ncap->value);
        ncap->value = NULL;
    }
    ncap->ndefined = 0;

#ifndef SEARCH_NAME_LINEARLY
    /* free space allocated for attribute name lookup table */
    if (ncap->nameT != NULL) {
        ncmpio_hash_table_free(ncap->nameT, ncap->hash_size);
        NCI_Free(ncap->nameT);
        ncap->nameT = NULL;
        ncap->hash_size = 0;
    }
#endif
}

/*----< ncmpio_dup_NC_attrarray() >------------------------------------------*/
int
ncmpio_dup_NC_attrarray(NC_attrarray *ncap, const NC_attrarray *ref)
{
    int i, status=NC_NOERR;

    assert(ref != NULL);
    assert(ncap != NULL);

    if (ref->ndefined == 0) { /* return now, if no attribute is defined */
        ncap->ndefined = 0;
        ncap->value    = NULL;
        return NC_NOERR;
    }

    if (ref->ndefined > 0) {
        size_t alloc_size = _RNDUP(ref->ndefined, PNC_ARRAY_GROWBY);
        ncap->value = (NC_attr **) NCI_Calloc(alloc_size, sizeof(NC_attr*));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    ncap->ndefined = 0;
    for (i=0; i<ref->ndefined; i++) {
        status = dup_NC_attr(ref->value[i], &ncap->value[i]);
        if (status != NC_NOERR) {
            ncmpio_free_NC_attrarray(ncap);
            return status;
        }
        ncap->ndefined++;
    }

    assert(ncap->ndefined == ref->ndefined);

#ifndef SEARCH_NAME_LINEARLY
    /* allocate hashing lookup table, if not allocated yet */
    if (ncap->nameT == NULL)
        ncap->nameT = NCI_Calloc(ncap->hash_size, sizeof(NC_nametable));

    /* duplicate attribute name lookup table */
    ncmpio_hash_table_copy(ncap->nameT, ref->nameT, ncap->hash_size);
#endif

    return NC_NOERR;
}

/*----< incr_NC_attrarray() >------------------------------------------------*/
/* Add a new handle at the end of an array of handles */
static int
incr_NC_attrarray(int isGlobal, NC_attrarray *ncap, NC_attr *new_attr)
{
    size_t growby = (isGlobal) ? PNC_ARRAY_GROWBY : PNC_VATTR_ARRAY_GROWBY;

    assert(ncap != NULL);
    assert(new_attr != NULL);

    if (ncap->ndefined % growby == 0) {
        /* grow the array to accommodate the new handle */
        size_t alloc_size = (size_t)ncap->ndefined + growby;

        ncap->value = (NC_attr **) NCI_Realloc(ncap->value,
                                   alloc_size * sizeof(NC_attr*));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    ncap->value[ncap->ndefined++] = new_attr;

    return NC_NOERR;
}

/* End attrarray per se */

/*----< NC_attrarray0() >----------------------------------------------------*/
/* Given ncp and varid, return pointer to array of attributes
 * else NULL on error. This can also be used to validate varid.
 */
static NC_attrarray *
NC_attrarray0(NC *ncp, int varid)
{
    if (varid == NC_GLOBAL) /* Global attribute */
        return &ncp->attrs;

    if (varid >= 0 && varid < ncp->vars.ndefined)
        return &ncp->vars.value[varid]->attrs;

    return NULL;
}

/*----< ncmpio_NC_findattr() >------------------------------------------------*/
/* Step thru NC_ATTRIBUTE array, seeking match on name.
 * return match or -1 if Not Found.
 */
int
ncmpio_NC_findattr(const NC_attrarray *ncap,
                   const char         *name) /* normalized string */
{
#ifndef SEARCH_NAME_LINEARLY
    int key;
#endif
    int i;
    size_t nchars;

    assert(ncap != NULL);

    if (ncap->ndefined == 0) return -1; /* none created yet */

#ifdef SEARCH_NAME_LINEARLY
    /* we assume the number of attributes is small and use the
     * following linear search. If the number is expected to be large, then
     * it is better to use the hashing as for variables and dimensions.
     */
    nchars = strlen(name);
    for (i=0; i<ncap->ndefined; i++) {
        if (ncap->value[i]->name_len == nchars &&
            strcmp(ncap->value[i]->name, name) == 0)
            return i;
    }
#else
    /* hash name into a key for name lookup */
    key = HASH_FUNC(name, ncap->hash_size);

    /* check the list (all names sharing the same key) using linear search */
    nchars = strlen(name);
    for (i=0; i<ncap->nameT[key].num; i++) {
        int attr_id = ncap->nameT[key].list[i];
        if (ncap->value[attr_id]->name_len == nchars &&
            strcmp(name, ncap->value[attr_id]->name) == 0) {
            return attr_id; /* the name already exists */
        }
    }
#endif
    return -1; /* the name has never been used */
}

/*----< NC_lookupattr() >----------------------------------------------------*/
/* Look up by ncid, ncap, and name */
static int
NC_lookupattr(const NC_attrarray  *ncap,
              const char          *name,   /* normalized attribute name */
              NC_attr            **attrpp) /* modified on return */
{
    int indx;

    /* requires validity of ncid and ncap already been checked */

    indx = ncmpio_NC_findattr(ncap, name);
    if (indx == -1) DEBUG_RETURN_ERROR(NC_ENOTATT)

    if (attrpp != NULL)
        *attrpp = ncap->value[indx];

    return NC_NOERR;
}

/* Public */

/*----< ncmpio_inq_attname() >-----------------------------------------------*/
/* This is an independent subroutine */
int
ncmpio_inq_attname(void *ncdp,
                   int   varid,
                   int   attid,
                   char *name)   /* out */
{
    NC *ncp=(NC*)ncdp;
    NC_attrarray *ncap;
    NC_attr *attrp;

    /* check varid and get pointer to the NC_attrarray */
    ncap = NC_attrarray0(ncp, varid);
    if (ncap == NULL) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* check attribute ID */
    if ((attid < 0) || ncap->ndefined == 0 || attid >= ncap->ndefined)
        DEBUG_RETURN_ERROR(NC_ENOTATT)

    assert(ncap->value != NULL);

    attrp = ncap->value[attid];

    if (name == NULL) DEBUG_RETURN_ERROR(NC_EINVAL)

    /* in PnetCDF, attrp->name is always NULL character terminated */
    strcpy(name, attrp->name);

    return NC_NOERR;
}

/*----< ncmpio_inq_attid() >-------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpio_inq_attid(void       *ncdp,
                 int         varid,
                 const char *name,
                 int        *attidp)  /* out */
{
    int indx, err;
    char *nname=NULL; /* normalized name */
    NC *ncp=(NC*)ncdp;
    NC_attrarray *ncap;

    ncap = NC_attrarray0(ncp, varid);
    if (ncap == NULL) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) return err;

    indx = ncmpio_NC_findattr(ncap, nname);
    NCI_Free(nname);
    if (indx == -1) DEBUG_RETURN_ERROR(NC_ENOTATT)

    if (attidp != NULL) *attidp = indx;

    return NC_NOERR;
}

/*----< ncmpio_inq_att() >---------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpio_inq_att(void       *ncdp,
               int         varid,
               const char *name, /* input, attribute name */
               nc_type    *datatypep,
               MPI_Offset *lenp)
{
    int err;
    char *nname=NULL;    /* normalized name */
    NC *ncp=(NC*)ncdp;
    NC_attr *attrp;
    NC_attrarray *ncap;

    ncap = NC_attrarray0(ncp, varid);
    if (ncap == NULL) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) return err;

    err = NC_lookupattr(ncap, nname, &attrp);
    NCI_Free(nname);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

    if (datatypep != NULL) *datatypep = attrp->xtype;

    if (lenp != NULL) *lenp = attrp->nelems;

    return NC_NOERR;
}

/*----< ncmpio_rename_att() >-------------------------------------------------*/
/* This API is collective. If the new name is longer than the old name, this
 * API must be called in define mode.
 */
int
ncmpio_rename_att(void       *ncdp,
                  int         varid,
                  const char *name,
                  const char *newname)
{
    int attr_id=-1, err;
    char *nname=NULL;    /* normalized name */
    char *nnewname=NULL; /* normalized newname */
    size_t nnewname_len=0;
    NC *ncp=(NC*)ncdp;
    NC_attrarray *ncap=NULL;
    NC_attr *attrp=NULL;

    ncap = NC_attrarray0(ncp, varid);
    if (ncap == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) goto err_check;

    attr_id = ncmpio_NC_findattr(ncap, nname);
    if (attr_id < 0) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTATT)
        goto err_check;
    }

    attrp = ncap->value[attr_id];

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(newname, &nnewname);
    if (err != NC_NOERR) goto err_check;

    nnewname_len = strlen(nnewname);

    if (ncmpio_NC_findattr(ncap, nnewname) >= 0) {
        /* name in use */
        DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
        goto err_check;
    }

    if (! NC_indef(ncp) && attrp->name_len < nnewname_len) {
        /* when data mode, nnewname cannot be longer than the old one */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

err_check:
    if (nname != NULL) NCI_Free(nname);
    if (ncp->safe_mode) {
        int minE, mpireturn;

        /* check error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nnewname != NULL) NCI_Free(nnewname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        }
        if (minE != NC_NOERR) {
            if (nnewname != NULL) NCI_Free(nnewname);
            return minE;
        }
    }

    if (err != NC_NOERR) {
        if (nnewname != NULL) NCI_Free(nnewname);
        return err;
    }

    assert(attrp != NULL);

#ifndef SEARCH_NAME_LINEARLY
    ncmpio_hash_replace(ncap->nameT, ncap->hash_size, attrp->name, nnewname, attr_id);
#endif

    /* replace the old name with new name */
    NCI_Free(attrp->name);
    attrp->name     = nnewname;
    attrp->name_len = nnewname_len;

    if (! NC_indef(ncp)) { /* when file is in data mode */
        /* Let root write the entire header to the file. Note that we cannot
         * just update the name in-place in file header, because if the file
         * space occupied by the name shrinks, all the metadata following it
         * must be moved ahead.
         */
        err = ncmpio_write_header(ncp);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
    }

    return err;
}


/*----< ncmpio_copy_att() >---------------------------------------------------*/
/* This API is collective for processes that opened ncdp_out.
 * If the attribute does not exist in ncdp_out, then this API must be called
 * when ncdp_out is in define mode.
 * If the attribute does exist in ncdp_out and the attribute in ncdp_in is
 * larger than the one in ncdp_out, then this API must be called when ncdp_out
 * is in define mode.
 */
int
ncmpio_copy_att(void       *ncdp_in,
                int         varid_in,
                const char *name,
                void       *ncdp_out,
                int         varid_out)
{
    int indx=0, err;
    char *nname=NULL;    /* normalized name */
    NC *ncp_in=(NC*)ncdp_in, *ncp_out=(NC*)ncdp_out;
    NC_attrarray *ncap_out=NULL, *ncap_in;
    NC_attr *iattrp=NULL, *attrp=NULL;

    ncap_in = NC_attrarray0(ncp_in, varid_in);
    if (ncap_in == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    ncap_out = NC_attrarray0(ncp_out, varid_out);
    if (ncap_out == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) goto err_check;

    err = NC_lookupattr(ncap_in, nname, &iattrp);
    if (err != NC_NOERR) {
        assert(iattrp == NULL);
        DEBUG_TRACE_ERROR(err)
        goto err_check;
    }

    indx = ncmpio_NC_findattr(ncap_out, nname);

    if (indx >= 0) { /* name in use in ncap_out */
        if (ncdp_in == ncdp_out && varid_in == varid_out)
            /* self copy is not considered an error */
            goto err_check;

        if (!NC_indef(ncp_out) &&  /* not allowed in data mode */
            iattrp->xsz > ncap_out->value[indx]->xsz) {
            DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
            goto err_check;
        }
    }
    else { /* attribute does not exist in ncdp_out */
        if (!NC_indef(ncp_out)) {
            /* add new attribute is not allowed in data mode */
            DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
            goto err_check;
        }
        /* Note we no longer limit the number of attributes, as CDF file
         * formats impose no such limit. Thus, the value of NC_MAX_ATTRS has
         * been changed NC_MAX_INT, as NC_attrarray.ndefined is of type signed
         * int and so is natts argument in ncmpi_inq_varnatts()
         */
        if (ncap_out->ndefined == NC_MAX_ATTRS) {
            DEBUG_ASSIGN_ERROR(err, NC_EMAXATTS)
            goto err_check;
        }
    }

err_check:
    if (ncp_out->safe_mode) {
        int minE, mpireturn;

        /* check the error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN,
                                  ncp_out->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nname != NULL) NCI_Free(nname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        }
        if (minE != NC_NOERR) {
            if (nname != NULL) NCI_Free(nname);
            return minE;
        }
    }

    if (err != NC_NOERR) {
        if (nname != NULL) NCI_Free(nname);
        return err;
    }
    assert(ncap_out != NULL);
    assert(nname != NULL);

    if (indx >= 0) { /* name in use in ncdp_out */
        NCI_Free(nname);
        if (ncdp_in == ncdp_out && varid_in == varid_out) {
            /* self copy is not considered an error */
            return NC_NOERR;
        }

        /* reuse existing attribute array slot without redef */
        attrp = ncap_out->value[indx];

        if (iattrp->xsz > attrp->xsz) {
            if (attrp->xvalue != NULL) NCI_Free(attrp->xvalue);
            attrp->xvalue = NCI_Malloc((size_t)iattrp->xsz);
            if (attrp->xvalue == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
        attrp->xsz    = iattrp->xsz;
        attrp->xtype  = iattrp->xtype;
        attrp->nelems = iattrp->nelems;
    }
    else { /* attribute does not exist in ncdp_out */
        err = ncmpio_new_NC_attr(nname, iattrp->name_len, iattrp->xtype,
                                 iattrp->nelems, &attrp);
        if (err != NC_NOERR) return err;

#ifndef SEARCH_NAME_LINEARLY
        /* allocate hashing lookup table, if not allocated yet */
        if (ncap_out->nameT == NULL)
            ncap_out->nameT = NCI_Calloc(ncap_out->hash_size, sizeof(NC_nametable));

        /* insert nname to the lookup table */
        ncmpio_hash_insert(ncap_out->nameT, ncap_out->hash_size, nname, ncap_out->ndefined);
#endif

        err = incr_NC_attrarray((varid_out == NC_GLOBAL), ncap_out, attrp);
        if (err != NC_NOERR) return err;
    }

    if (iattrp->xsz > 0)
        memcpy(attrp->xvalue, iattrp->xvalue, (size_t)iattrp->xsz);

    if (!NC_indef(ncp_out)) { /* called in data mode */
        /* Let root write the entire header to the file. Note that we
         * cannot just update the variable name in its space occupied in
         * the file header, because if the file space occupied by the name
         * shrinks, all the metadata following it must be moved ahead.
         */
        err = ncmpio_write_header(ncp_out); /* update file header */
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
    }

    return err;
}

/*----< ncmpio_del_att() >---------------------------------------------------*/
/* This is a collective subroutine and must be called in define mode */
int
ncmpio_del_att(void       *ncdp,
               int         varid,
               const char *name)
{
    int err, attrid=-1;
    char *nname=NULL; /* normalized name */
    NC *ncp=(NC*)ncdp;
    NC_attrarray *ncap=NULL;

    /* check NC_ENOTVAR */
    ncap = NC_attrarray0(ncp, varid);
    if (ncap == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTVAR)
        goto err_check;
    }

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) goto err_check;

    attrid = ncmpio_NC_findattr(ncap, nname);
    if (attrid == -1) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTATT)
        goto err_check;
    }

#ifndef SEARCH_NAME_LINEARLY
    /* delete name entry from hash table */
    err = ncmpio_hash_delete(ncap->nameT, ncap->hash_size, nname, attrid);
    if (err != NC_NOERR) goto err_check;
#endif

err_check:
    if (nname != NULL) NCI_Free(nname);
    if (ncp->safe_mode) {
        int minE, mpireturn;

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN,ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }

    if (err != NC_NOERR) return err;
    assert(ncap != NULL);

    /* delete attribute */
    if (ncap->value[attrid]->xvalue != NULL)
        NCI_Free(ncap->value[attrid]->xvalue);
    NCI_Free(ncap->value[attrid]->name);
    NCI_Free(ncap->value[attrid]);

    /* shuffle down */
    for (; attrid < ncap->ndefined-1; attrid++)
        ncap->value[attrid] = ncap->value[attrid+1];

    /* decrement count */
    ncap->ndefined--;

    return NC_NOERR;
}

dnl
dnl GET_ATT(fntype)
dnl
define(`GET_ATT',dnl
`dnl
/*----< get_att_$1() >-------------------------------------------------------*/
/* This is an independent subroutine */
static int
get_att_$1(nc_type       xtype,
           const void  **xp,
           MPI_Offset    nelems,
           NC2ITYPE($1) *buf)
{
    switch(xtype) {
        /* possible error returned of this switch block is NC_ERANGE */
        case NC_BYTE:
            return ncmpix_pad_getn_NC_BYTE_$1  (xp, nelems, buf);
        case NC_UBYTE:
            return ncmpix_pad_getn_NC_UBYTE_$1 (xp, nelems, buf);
        case NC_SHORT:
            return ncmpix_pad_getn_NC_SHORT_$1 (xp, nelems, buf);
        case NC_USHORT:
            return ncmpix_pad_getn_NC_USHORT_$1(xp, nelems, buf);
        case NC_INT:
            return ncmpix_getn_NC_INT_$1       (xp, nelems, buf);
        case NC_UINT:
            return ncmpix_getn_NC_UINT_$1      (xp, nelems, buf);
        case NC_FLOAT:
            return ncmpix_getn_NC_FLOAT_$1     (xp, nelems, buf);
        case NC_DOUBLE:
            return ncmpix_getn_NC_DOUBLE_$1    (xp, nelems, buf);
        case NC_INT64:
            return ncmpix_getn_NC_INT64_$1     (xp, nelems, buf);
        case NC_UINT64:
            return ncmpix_getn_NC_UINT64_$1    (xp, nelems, buf);
        default:
            /* this error is unlikely, but an internal error if happened */
            fprintf(stderr, "Error: bad attrp->xtype(%d) in %s\n",
                    xtype,__func__);
            DEBUG_RETURN_ERROR(NC_EBADTYPE)
    }
}
')dnl
dnl
foreach(`iType', (schar,uchar,short,ushort,int,uint,float,double,longlong,ulonglong),
        `GET_ATT(iType)dnl
')dnl
dnl
/*----< ncmpio_get_att() >---------------------------------------------------*/
/* This is an independent subroutine.
 * itype: the internal data type of buf. Its value can be MPI_DATATYPE_NULL or
 *        MPI promitive data type. For MPI_DATATYPE_NULL, it means the data
 *        type of buf matches the external type of attribute defined in file.
 *        For other MPI promitive data type, it corresponds to the type names
 *        shown in the API name.
 */
int
ncmpio_get_att(void         *ncdp,
               int           varid,
               const char   *name,
               void         *buf,
               MPI_Datatype  itype)
{
    int           err;
    char         *nname=NULL; /* normalized name */
    NC           *ncp=(NC*)ncdp;
    NC_attr      *attrp;
    NC_attrarray *ncap=NULL;
    const void   *xp;
    nc_type       xtype;
    MPI_Offset    nelems;

    /* sanity checks for varid and name has been done in dispatcher */

    /* obtain NC_attrarray object pointer */
    if (varid == NC_GLOBAL) ncap = &ncp->attrs;
    else                    ncap = &ncp->vars.value[varid]->attrs;

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) return err;

    /* whether the attr exists (check NC_ENOTATT) */
    err = NC_lookupattr(ncap, nname, &attrp);
    NCI_Free(nname);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

    if (attrp->nelems == 0) return NC_NOERR;

    if (itype == MPI_DATATYPE_NULL)
        /* this is called from API ncmpi_get_att() where the internal and
         * external data types match */
        itype = ncmpii_nc2mpitype(attrp->xtype);

    /* No character conversions are allowed. */
    if ((attrp->xtype == NC_CHAR && itype != MPI_CHAR) ||
        (attrp->xtype != NC_CHAR && itype == MPI_CHAR))
        DEBUG_RETURN_ERROR(NC_ECHAR)

    if (buf == NULL) DEBUG_RETURN_ERROR(NC_EINVAL)

    /* must use xp, as ncmpix_pad_getn_text moves xp ahead */
    xp     = attrp->xvalue;
    xtype  = attrp->xtype;
    nelems = attrp->nelems;

    /* cannot use switch statement, as MPI_Datatype may not be int,
     * for instance, it is a C struct in OpenMPI
     */
    if (itype == MPI_CHAR)
        return ncmpix_pad_getn_text    (&xp, nelems,      (char*)buf);
    else if (itype == MPI_SIGNED_CHAR)
        return get_att_schar    (xtype, &xp, nelems,     (schar*)buf);
    else if (itype == MPI_UNSIGNED_CHAR) {
        if (xtype == NC_BYTE && ncp->format < NC_FORMAT_CDF5)
            xtype = NC_UBYTE; /* special case: no NC_ERANGE check */
        return get_att_uchar    (xtype, &xp, nelems,     (uchar*)buf);
    }
    else if (itype == MPI_SHORT)
        return get_att_short    (xtype, &xp, nelems,     (short*)buf);
    else if (itype == MPI_UNSIGNED_SHORT)
        return get_att_ushort   (xtype, &xp, nelems,    (ushort*)buf);
    else if (itype == MPI_INT)
        return get_att_int      (xtype, &xp, nelems,       (int*)buf);
    else if (itype == MPI_UNSIGNED)
        return get_att_uint     (xtype, &xp, nelems,      (uint*)buf);
    else if (itype == MPI_FLOAT)
        return get_att_float    (xtype, &xp, nelems,     (float*)buf);
    else if (itype == MPI_DOUBLE)
        return get_att_double   (xtype, &xp, nelems,    (double*)buf);
    else if (itype == MPI_LONG_LONG_INT)
        return get_att_longlong (xtype, &xp, nelems,  (longlong*)buf);
    else if (itype == MPI_UNSIGNED_LONG_LONG)
        return get_att_ulonglong(xtype, &xp, nelems, (ulonglong*)buf);
    DEBUG_RETURN_ERROR(NC_EBADTYPE)
}

dnl
dnl PUTN_ITYPE(_pad, itype)
dnl
define(`PUTN_ITYPE',dnl
`dnl
/*----< putn_$1() >----------------------------------------------------------*/
/* This is a collective subroutine */
static int
putn_$1(void       **xpp,    /* buffer to be written to file */
        MPI_Offset   nelems, /* no. elements in user buffer */
        const $1    *buf,    /* user buffer of type $1 */
        nc_type      xtype,  /* external NC type */
        void        *fillp)  /* fill value in internal representation */
{
    switch(xtype) {
        case NC_BYTE:
            return ncmpix_pad_putn_NC_BYTE_$1  (xpp, nelems, buf, fillp);
        case NC_UBYTE:
            return ncmpix_pad_putn_NC_UBYTE_$1 (xpp, nelems, buf, fillp);
        case NC_SHORT:
            return ncmpix_pad_putn_NC_SHORT_$1 (xpp, nelems, buf, fillp);
        case NC_USHORT:
            return ncmpix_pad_putn_NC_USHORT_$1(xpp, nelems, buf, fillp);
        case NC_INT:
            return ncmpix_putn_NC_INT_$1       (xpp, nelems, buf, fillp);
        case NC_UINT:
            return ncmpix_putn_NC_UINT_$1      (xpp, nelems, buf, fillp);
        case NC_FLOAT:
            return ncmpix_putn_NC_FLOAT_$1     (xpp, nelems, buf, fillp);
        case NC_DOUBLE:
            return ncmpix_putn_NC_DOUBLE_$1    (xpp, nelems, buf, fillp);
        case NC_INT64:
            return ncmpix_putn_NC_INT64_$1     (xpp, nelems, buf, fillp);
        case NC_UINT64:
            return ncmpix_putn_NC_UINT64_$1    (xpp, nelems, buf, fillp);
        case NC_CHAR:
            DEBUG_RETURN_ERROR(NC_ECHAR) /* NC_ECHAR check is done earlier */
        default: fprintf(stderr, "Error: bad xtype(%d) in %s\n",xtype,__func__);
            DEBUG_RETURN_ERROR(NC_EBADTYPE)
    }
}
')dnl
dnl
foreach(`iType', (schar,uchar,short,ushort,int,uint,float,double,longlong,ulonglong),
        `PUTN_ITYPE(iType)dnl
')dnl

/* For netCDF, the type mapping between file types and buffer types
 * are based on netcdf4. Check APIs of nc_put_att_xxx from source files
 *     netCDF/netcdf-x.x.x/libdispatch/att.c
 *     netCDF/netcdf-x.x.x/libsrc4/nc4attr.c
 *
 * Note that schar means signed 1-byte integers in attributes. Hence the call
 * below is illegal (NC_ECHAR will return), indicating the error on trying
 * type conversion between characters and numbers.
 *
 * ncmpi_put_att_schar(ncid, varid, "attr name", NC_CHAR, strlen(attrp), attrp);
 *
 * This rule and mapping apply for variables as well. See APIs of
 * nc_put_vara_xxx from source files
 *     netCDF/netcdf-x.x.x/libdispatch/var.c
 *     netCDF/netcdf-x.x.x/libsrc4/nc4var.c
 *
 */

/*----< ncmpio_put_att() >---------------------------------------------------*/
/* This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, then this
 * API must be called when the file is in define mode.
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters
 *
 * Note ncmpi_put_att_text will never return NC_ERANGE error, as text is not
 * convertible to numerical types.
 *
 * itype: the internal data type of buf. Its value can be MPI_DATATYPE_NULL or
 *        MPI promitive data type. For MPI_DATATYPE_NULL, it means the data
 *        type of buf matches the external type of attribute defined in file.
 *        For other MPI promitive data type, it corresponds to the type names
 *        shown in the API name.
 */
int
ncmpio_put_att(void         *ncdp,
               int           varid,
               const char   *name,
               nc_type       xtype,  /* external (file/NC) data type */
               MPI_Offset    nelems,
               const void   *buf,
               MPI_Datatype  itype)  /* internal (memory) data type */
{
    int indx=0, err;
    char *nname=NULL; /* normalized name */
    MPI_Offset xsz=0;
    NC *ncp=(NC*)ncdp;
    NC_attrarray *ncap=NULL;
    NC_attr *attrp=NULL;

    /* sanity checks for varid, name, xtype has been done in dispatcher */

    /* If this is the _FillValue attribute, then let PnetCDF return the
     * same error codes as netCDF
     */
    if (varid != NC_GLOBAL && !strcmp(name, _FillValue)) {
        /* Fill value must be of the same data type */
        if (xtype != ncp->vars.value[varid]->xtype) {
            DEBUG_ASSIGN_ERROR(err, NC_EBADTYPE)
            goto err_check;
        }

        /* Fill value must have exactly one value */
        if (nelems != 1) {
            DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
            goto err_check;
        }

        /* Only allowed in initial define mode (i.e. variable is newly defined) */
        if (ncp->old != NULL && varid < ncp->old->vars.ndefined) {
            DEBUG_ASSIGN_ERROR(err, NC_ELATEFILL)
            goto err_check;
        }
    }

    xsz = x_len_NC_attrV(xtype, nelems);
    /* xsz is the total aligned size of this attribute */

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) goto err_check;

    /* obtain NC_attrarray object pointer, ncap */
    if (varid == NC_GLOBAL) ncap = &ncp->attrs;
    else                    ncap = &ncp->vars.value[varid]->attrs;

    /* check whether attribute already exists */
    indx = ncmpio_NC_findattr(ncap, nname);

    if (indx >= 0) { /* name in use */
        /* xsz is the total size of this attribute */
        if (!NC_indef(ncp) && xsz > ncap->value[indx]->xsz) {
            /* The new attribute requires a larger space, which is not allowed
             * in data mode */
            DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
            goto err_check;
        }
    }
    else { /* attribute does not exist in ncid */
        if (!NC_indef(ncp)) {
            /* add new attribute is not allowed in data mode */
            DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
            goto err_check;
        }
        /* Note we no longer limit the number of attributes, as CDF file formats
         * impose no such limit. Thus, the value of NC_MAX_ATTRS has been
         * changed NC_MAX_INT, as NC_attrarray.ndefined is of type signed int
         * and so is natts argument in ncmpi_inq_varnatts()
         */
        if (ncap->ndefined == NC_MAX_ATTRS) {
            DEBUG_ASSIGN_ERROR(err, NC_EMAXATTS)
            goto err_check;
        }
    }

err_check:
    if (ncp->safe_mode) { /* check the error code across processes */
        int minE, mpireturn;

        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nname != NULL) NCI_Free(nname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        }
        if (minE != NC_NOERR) {
            if (nname != NULL) NCI_Free(nname);
            return minE;
        }
        /* argument consistency check has been done at the dispatchers */
    }

    if (err != NC_NOERR) {
        if (nname != NULL) NCI_Free(nname);
        return err;
    }
    assert(ncap != NULL);
    assert(nname != NULL);

    if (indx >= 0) { /* name in use */
        NCI_Free(nname);
        attrp = ncap->value[indx]; /* convenience */

        if (xsz > attrp->xsz) { /* new attribute requires a larger space */
            if (attrp->xvalue != NULL) NCI_Free(attrp->xvalue);
            attrp->xvalue = NCI_Malloc((size_t)xsz);
            if (attrp->xvalue == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
        attrp->xsz    = xsz;
        attrp->xtype  = xtype;
        attrp->nelems = nelems;
    }
    else { /* attribute does not exist in ncid */
        err = ncmpio_new_NC_attr(nname, strlen(nname), xtype, nelems, &attrp);
        if (err != NC_NOERR) return err;

#ifndef SEARCH_NAME_LINEARLY
        /* allocate hashing lookup table, if not allocated yet */
        if (ncap->nameT == NULL)
            ncap->nameT = NCI_Calloc(ncap->hash_size, sizeof(NC_nametable));

        /* insert nname to the lookup table */
        ncmpio_hash_insert(ncap->nameT, ncap->hash_size, nname, ncap->ndefined);
#endif

        err = incr_NC_attrarray((varid == NC_GLOBAL), ncap, attrp);
        if (err != NC_NOERR) return err;
    }

    if (nelems != 0 && buf != NULL) { /* non-zero length attribute */
        /* using xp below to prevent change the pointer attr->xvalue, as
         * ncmpix_pad_putn_<type>() advances the first argument with nelems
         * elements. Note that attrp->xvalue is malloc-ed with a buffer of
         * size that is aligned with a 4-byte boundary.
         */
        unsigned char fill[8]; /* fill value in internal representation */
        void *xp = attrp->xvalue;

        if (itype == MPI_CHAR) {
            err = ncmpix_pad_putn_text(&xp, nelems, buf);
        }
        else if (ncp->format < 5 && itype == MPI_UNSIGNED_CHAR &&
                 xtype == NC_BYTE) { /* special case: no NC_ERANGE check */
            /* find the fill value of NC_UBYTE */
            err = ncmpio_inq_default_fill_value(NC_UBYTE, &fill);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
            err = putn_uchar(&xp, nelems, buf, NC_UBYTE, &fill);
        }
        else {
            /* find the fill value */
            err = ncmpio_inq_default_fill_value(xtype, &fill);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

            /* cannot use switch statement, as MPI_Datatype may not be int,
             * for instance, it is a C struct in OpenMPI
             */
            /* MPI_CHAR: _text API has been handled above */
            if (itype == MPI_SIGNED_CHAR)
                err = putn_schar    (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_UNSIGNED_CHAR)
                err = putn_uchar    (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_SHORT)
                err = putn_short    (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_UNSIGNED_SHORT)
                err = putn_ushort   (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_INT)
                err = putn_int      (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_UNSIGNED)
                err = putn_uint     (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_FLOAT)
                err = putn_float    (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_DOUBLE)
                err = putn_double   (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_LONG_LONG_INT)
                err = putn_longlong (&xp, nelems, buf, xtype, &fill);
            else if (itype == MPI_UNSIGNED_LONG_LONG)
                err = putn_ulonglong(&xp, nelems, buf, xtype, &fill);
            else DEBUG_ASSIGN_ERROR(err, NC_EBADTYPE)
        }

        /* no immediately return error code here? Strange ...
         * Instead, we continue and call incr_NC_attrarray() to add
         * this attribute (for create case) as it is legal. But if
         * we return error and reject this attribute, then nc_test will
         * fail with this error message below:
         * FAILURE at line 252 of test_read.c: ncmpi_inq: wrong number
         * of global atts returned, 3
         * Check netCDF-4, it is doing the same thing!
         *
         * One of the error codes returned from ncmpix_pad_putn_<type>() is
         * NC_ERANGE, meaning one or more elements are type overflow.
         * Should we reject the entire attribute array if only part of
         * the array overflow? For netCDF4, the answer is NO.
         */
/*
        if (err != NC_NOERR) {
            if (attrp->xvalue != NULL) NCI_Free(attrp->xvalue);
            NCI_Free(attrp->name);
            NCI_Free(attrp);
            DEBUG_RETURN_ERROR(err)
        }
*/
    }

    if (!NC_indef(ncp)) { /* called in data mode */
        /* Let root write the entire header to the file. Note that we
         * cannot just update the attribute in its space occupied in the
         * file header, because if the file space occupied by the attribute
         * shrinks, all the metadata following it must be moved ahead.
         */
        int status;
        status = ncmpio_write_header(ncp); /* update file header */
        if (err == NC_NOERR) err = status;
    }

    return err;
}
