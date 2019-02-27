/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in
 * src/dispatchers/variable.c
 *
 * ncmpi_def_var()    : dispatcher->def_var()
 * ncmpi_inq_varid()  : dispatcher->inq_varid()
 * ncmpi_inq_var()    : dispatcher->inq_var()
 * ncmpi_rename_var() : dispatcher->rename_var()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h> /* memset() */
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncx.h>
#include "ncmpio_NC.h"

/*----< ncmpio_free_NC_var() >-----------------------------------------------*/
/* Free NC_var object */
void
ncmpio_free_NC_var(NC_var *varp)
{
    if (varp == NULL) return;

    ncmpio_free_NC_attrarray(&varp->attrs);
    NCI_Free(varp->name);
#ifdef ENABLE_SUBFILING
    if (varp->num_subfiles > 1) /* deallocate it */
        NCI_Free(varp->dimids_org);
#endif
    if (varp->shape  != NULL) NCI_Free(varp->shape);
    if (varp->dsizes != NULL) NCI_Free(varp->dsizes);
    if (varp->dimids != NULL) NCI_Free(varp->dimids);

    NCI_Free(varp);
}

/*----< ncmpio_new_NC_var() >------------------------------------------------*/
NC_var *
ncmpio_new_NC_var(char *name, int ndims)
{
    NC_var *varp;

    varp = (NC_var *) NCI_Calloc(1, sizeof(NC_var));
    if (varp == NULL) return NULL;

    if (ndims > 0) {
        varp->shape  = (MPI_Offset*)NCI_Calloc(ndims, SIZEOF_MPI_OFFSET);
        varp->dsizes = (MPI_Offset*)NCI_Calloc(ndims, SIZEOF_MPI_OFFSET);
        varp->dimids = (int *)      NCI_Calloc(ndims, SIZEOF_INT);
    }

    varp->name     = name;         /* name has been malloc-ed */
    varp->name_len = strlen(name); /* name has been NULL checked */
    varp->ndims    = ndims;

    return varp;
}

/*----< dup_NC_var() >-------------------------------------------------------*/
static NC_var *
dup_NC_var(const NC_var *rvarp)
{
    char *name;
    NC_var *varp;

    /* note that name in rvarp->name is already normalized */
    name = (char*) NCI_Malloc(strlen(rvarp->name)+1);
    if (name == NULL) return NULL;
    strcpy(name, rvarp->name);

    /* allocate a NC_var object */
    varp = ncmpio_new_NC_var(name, rvarp->ndims);
    if (varp == NULL ) return NULL;

    varp->xtype = rvarp->xtype;

    /* copy dimids[] */
    if (rvarp->ndims != 0 && rvarp->dimids != NULL)
        memcpy(varp->dimids, rvarp->dimids, (size_t)rvarp->ndims * SIZEOF_INT);

    /* copy attributes */
    if (ncmpio_dup_NC_attrarray(&varp->attrs, &rvarp->attrs) != NC_NOERR) {
        ncmpio_free_NC_var(varp);
        return NULL;
    }

    /* copy the contents of shape may not be necessary, as one must call
     * compute_var_shape() to recompute it after a new variable is created
     */
    memcpy(varp->shape,  rvarp->shape,  (size_t)rvarp->ndims * SIZEOF_MPI_OFFSET);
    memcpy(varp->dsizes, rvarp->dsizes, (size_t)rvarp->ndims * SIZEOF_MPI_OFFSET);
    varp->xsz   = rvarp->xsz;
    varp->len   = rvarp->len;
    varp->begin = rvarp->begin;

    return varp;
}

/* vararray */

/*----< ncmpio_free_NC_vararray() >------------------------------------------*/
void
ncmpio_free_NC_vararray(NC_vararray *ncap)
{
    int i;

    assert(ncap != NULL);
    if (ncap->ndefined == 0) return;

    if (ncap->value != NULL) {
        /* when error is detected reading NC_VARIABLE tag, ncap->ndefined can
         * be > 0 and ncap->value is still NULL
         */
        for (i=0; i<ncap->ndefined; i++) {
            if (ncap->value[i] != NULL)
                ncmpio_free_NC_var(ncap->value[i]);
        }
        NCI_Free(ncap->value);
        ncap->value    = NULL;
    }
    ncap->ndefined = 0;

#ifndef SEARCH_NAME_LINEARLY
    /* free space allocated for var name lookup table */
    ncmpio_hash_table_free(ncap->nameT);
#endif
}

/*----< ncmpio_dup_NC_vararray() >-------------------------------------------*/
int
ncmpio_dup_NC_vararray(NC_vararray       *ncap,
                       const NC_vararray *ref)
{
    int i, status=NC_NOERR;

    assert(ref != NULL);
    assert(ncap != NULL);

    if (ref->ndefined == 0) {
        ncap->ndefined = 0;
        ncap->value = NULL;
        return NC_NOERR;
    }

    if (ref->ndefined > 0) {
        size_t alloc_size = _RNDUP(ref->ndefined, NC_ARRAY_GROWBY);
        ncap->value = (NC_var **) NCI_Calloc(alloc_size, sizeof(NC_var*));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    /* duplicate one NC_var object at a time */
    ncap->ndefined = 0;
    for (i=0; i<ref->ndefined; i++) {
        ncap->value[i] = dup_NC_var(ref->value[i]);
        if (ncap->value[i] == NULL) {
            DEBUG_ASSIGN_ERROR(status, NC_ENOMEM)
            break;
        }
        ncap->ndefined++;
    }
    if (status != NC_NOERR) {
        ncmpio_free_NC_vararray(ncap);
        return status;
    }
    assert(ncap->ndefined == ref->ndefined);

#ifndef SEARCH_NAME_LINEARLY
    /* duplicate var name lookup table */
    ncmpio_hash_table_copy(ncap->nameT, ref->nameT);
#endif

    return NC_NOERR;
}


/* End vararray per se */


#ifdef SEARCH_NAME_LINEARLY
/*----< NC_findvar() >-------------------------------------------------------*/
/*
 * Step thru NC_VARIABLE array, seeking match on name.
 * If found, set the variable ID pointed by vardip, otherwise return NC_ENOTVAR
 */
static int
NC_findvar(const NC_vararray *ncap,
           const char        *name,  /* normalized name */
           int               *varidp)
{
    int varid;
    size_t nchars;

    assert (ncap != NULL);

    if (ncap->ndefined == 0) return NC_ENOTVAR;

    nchars = strlen(name);
    for (varid=0; varid<ncap->ndefined; varid++) {
        if (ncap->value[varid]->name_len == nchars &&
            strcmp(ncap->value[varid]->name, name) == 0) {
            if (varidp != NULL) *varidp = varid;
            return NC_NOERR; /* found it */
        }
    }

    return NC_ENOTVAR; /* not found */
}
#else
/*----< NC_findvar() >-------------------------------------------------------*/
/* Check if the name has been used.
 * If yes, set the variable ID pointed by vardip, otherwise return NC_ENOTVAR
 */
static int
NC_findvar(const NC_vararray  *ncap,
           const char         *name,  /* normalized name */
           int                *varidp)
{
    int i, key, varid;
    size_t nchars;

    assert (ncap != NULL);

    if (ncap->ndefined == 0) return NC_ENOTVAR;

    /* hash the var name into a key for name lookup */
    key = HASH_FUNC(name);

    /* check the list using linear search */
    nchars = strlen(name);
    for (i=0; i<ncap->nameT[key].num; i++) {
        varid = ncap->nameT[key].list[i];
        if (ncap->value[varid]->name_len == nchars &&
            strcmp(ncap->value[varid]->name, name) == 0) {
            if (varidp != NULL) *varidp = varid;
            return NC_NOERR; /* the name already exists */
        }
    }

    return NC_ENOTVAR; /* the name has never been used */
}
#endif

/*----< ncmpio_NC_var_shape64() >--------------------------------------------*/
/* set varp->xsz, varp->shape and varp->len of a variable */
int
ncmpio_NC_var_shape64(NC_var            *varp,
                      const NC_dimarray *dims)
{
    int i;
    MPI_Offset product = 1;

    if (varp->ndims == 0) goto out;

    /* determine shape[] of the variable */
    for (i=0; i<varp->ndims; i++) {
        /* For file create, varp->dimids[i] has been checked in ncmpi_def_var()
         * in dispatchers/variable.c. For file open, it has been checked in
         * hdr_get_NC_var() in ncmpio_header_get.c */

        varp->shape[i] = dims->value[varp->dimids[i]]->size;

        /* check for record variable, only the highest dimension can
         * be unlimited */
        if (varp->shape[i] == NC_UNLIMITED && i != 0)
            DEBUG_RETURN_ERROR(NC_EUNLIMPOS)
    }

    /*
     * compute dsizes[], from right to left product of shape
     * For example, a 3D array of size 5x4x3 in C order,
     * For fixed-size variable: dsizes[0]=60 dsizes[1]=12 dsizes[2]=3
     * For record     variable: dsizes[0]=12 dsizes[1]=12 dsizes[2]=3
     */
    product = 1;
    if (varp->ndims == 1) {
        if (varp->shape[0] == NC_UNLIMITED)
            varp->dsizes[0] = 1;
        else {
            varp->dsizes[0] = varp->shape[0];
            product = varp->shape[0];
        }
    }
    else { /* varp->ndims > 1 */
        varp->dsizes[varp->ndims-1] = varp->shape[varp->ndims-1];
        product = varp->shape[varp->ndims-1];
        for (i=varp->ndims-2; i>=0; i--) {
            if (varp->shape[i] != NC_UNLIMITED)
                product *= varp->shape[i];
            varp->dsizes[i] = product;
        }
    }

out :
    /* No variable size can be > X_INT64_MAX - 3 */
    if (0 == ncmpio_NC_check_vlen(varp, X_INT64_MAX-3)) return NC_EVARSIZE;

    /*
     * For CDF-1 and CDF-2 formats, the total number of array elements
     * cannot exceed 2^32, unless this variable is the last fixed-size
     * variable, there is no record variable, and the file starting
     * offset of this variable is less than 2GiB.
     * We will check this in ncmpi_enddef() which calls ncmpio_NC_enddef()
     * which calls ncmpio_NC_check_vlens()
     *
    if (ncp->format < 5 && product >= X_UINT_MAX)
        DEBUG_RETURN_ERROR(NC_EVARSIZE)
     */

    /*
     * align variable size to 4 byte boundary, required by all netcdf
     * file formats
     */
    varp->len = product * varp->xsz;
    if (varp->len % 4 > 0)
        varp->len += 4 - varp->len % 4; /* round up */

    return NC_NOERR;
}

/*----< ncmpio_def_var() >---------------------------------------------------*/
int
ncmpio_def_var(void       *ncdp,
               const char *name,
               nc_type     xtype,
               int         ndims,
               const int  *dimids,
               int        *varidp)
{
    int err=NC_NOERR;
    char *nname=NULL; /* normalized name */
    NC *ncp=(NC*)ncdp;
    NC_var *varp=NULL;

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) goto err_check;

    /* allocate a new NC_var object */
    varp = ncmpio_new_NC_var(nname, ndims);
    if (varp == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto err_check;
    }
    /* sanity check for xtype has been done at dispatchers */
    varp->xtype = xtype;
    ncmpii_xlen_nc_type(xtype, &varp->xsz);

    /* copy dimids[] */
    if (ndims != 0 && dimids != NULL)
        memcpy(varp->dimids, dimids, (size_t)ndims * SIZEOF_INT);

    /* set up array dimensional structures */
    err = ncmpio_NC_var_shape64(varp, &ncp->dims);
    if (err != NC_NOERR) {
        ncmpio_free_NC_var(varp);
        nname = NULL; /* already freed in ncmpio_free_NC_var() */
        goto err_check;
    }

    /* allocate/expand ncp->vars.value array */
    if (ncp->vars.ndefined % NC_ARRAY_GROWBY == 0) {
        size_t alloc_size = (size_t)ncp->vars.ndefined + NC_ARRAY_GROWBY;
        ncp->vars.value = (NC_var **) NCI_Realloc(ncp->vars.value,
                                      alloc_size * sizeof(NC_var*));
        if (ncp->vars.value == NULL) {
            ncmpio_free_NC_var(varp);
            nname = NULL; /* already freed in ncmpio_free_NC_var() */
            err = NC_ENOMEM;
            goto err_check;
        }
    }

    varp->varid = ncp->vars.ndefined; /* varid */

    /* Add a new handle to the end of an array of handles */
    ncp->vars.value[ncp->vars.ndefined] = varp;

    ncp->vars.ndefined++;

err_check:
    if (ncp->safe_mode) {
        int minE, mpireturn;

        /* first check the error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, ncp->comm);
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

    assert(nname != NULL);

#ifndef SEARCH_NAME_LINEARLY
    /* insert nname to the lookup table */
    ncmpio_hash_insert(ncp->vars.nameT, nname, varp->varid);
#endif

    if (varidp != NULL) *varidp = varp->varid;

    /* default is NOFILL */
    varp->no_fill = 1;

    /* change to FILL only if the entire dataset fill mode is FILL */
    if (NC_dofill(ncp)) varp->no_fill = 0;

    return NC_NOERR;
}


/*----< ncmpio_inq_varid() >-------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpio_inq_varid(void       *ncdp,
                 const char *name,
                 int        *varid)
{
    int err=NC_NOERR;
    char *nname=NULL; /* normalized name */
    NC *ncp=(NC*)ncdp;

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

    err = NC_findvar(&ncp->vars, nname, varid);
    NCI_Free(nname);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

    return NC_NOERR;
}

/*----< ncmpio_inq_var() >---------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpio_inq_var(void       *ncdp,
               int         varid,
               char       *name,
               nc_type    *xtypep,
               int        *ndimsp,
               int        *dimids,
               int        *nattsp,
               MPI_Offset *offsetp,
               int        *no_fillp,    /* OUT: 1 not fill mode, 0 fill mode */
               void       *fill_valuep) /* OUT: user-defined or default fill value */
{
    int err=NC_NOERR;
    NC *ncp=(NC*)ncdp;
    NC_var *varp=NULL;

    /* using NC_GLOBAL in varid is illegal for this API, except for inquiring
     * natts. See
     * http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2015/msg00196.html
     * Checking NC_GLOBAL has been done at the calling routines at top level.
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)
     */

    if (varid == NC_GLOBAL) {
        /* in this case, all other pointer arguments must be NULLs */
        if (nattsp != NULL)
            *nattsp = ncp->attrs.ndefined;
        return NC_NOERR;
    }

    varp = ncp->vars.value[varid];

    if (name != NULL)
        /* in PnetCDF, name is always NULL character terminated */
        strcpy(name, varp->name);

    if (xtypep != NULL)
        *xtypep = varp->xtype;

    if (ndimsp != NULL) {
#ifdef ENABLE_SUBFILING
        /* varp->num_subfiles is already set during open or enddef */
        if (varp->num_subfiles > 1)
            *ndimsp = varp->ndims_org;
        else
#endif
            *ndimsp = varp->ndims;
    }
    if (dimids != NULL) {
#ifdef ENABLE_SUBFILING
        /* varp->dimids_org is already set during open or enddef */
        if (varp->num_subfiles > 1)
            memcpy(dimids, varp->dimids_org, (size_t)varp->ndims_org * SIZEOF_INT);
        else
#endif
            memcpy(dimids, varp->dimids, (size_t)varp->ndims * SIZEOF_INT);
    }
    if (nattsp != NULL) *nattsp = varp->attrs.ndefined;

    if (offsetp != NULL) *offsetp = varp->begin;

    if (no_fillp != NULL) *no_fillp = varp->no_fill;

    if (fill_valuep != NULL) {
        err = ncmpio_inq_var_fill(varp, fill_valuep);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
    }

    return NC_NOERR;
}


/*----< ncmpio_rename_var() >------------------------------------------------*/
/* This API is collective.
 * If the new name is longer than the old name, the netCDF file must be in
 * define mode. Otherwise, it can be called in either define or data mode.
 */
int
ncmpio_rename_var(void       *ncdp,
                  int         varid,
                  const char *newname)
{
    int err=NC_NOERR;
    char *nnewname=NULL; /* normalized name */
    size_t nnewname_len=0;
    NC *ncp=(NC*)ncdp;
    NC_var *varp=NULL;

    /* check whether variable ID is valid */
    /* sanity check for ncdp and varid has been done in dispatchers */
    varp = ncp->vars.value[varid];

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(newname, &nnewname);
    if (err != NC_NOERR) goto err_check;

    nnewname_len = strlen(nnewname);

    if (! NC_indef(ncp) && varp->name_len < nnewname_len) {
        /* when in data mode, newname cannot be longer thne the old one */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

#ifndef SEARCH_NAME_LINEARLY
    /* update var name lookup table */
    err = ncmpio_update_name_lookup_table(ncp->vars.nameT, varid,
                        ncp->vars.value[varid]->name, nnewname);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR(err)
        goto err_check;
    }
#endif

err_check:
    if (ncp->safe_mode) {
        int minE, mpireturn;

        /* First check error code so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN,ncp->comm);
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

    assert(varp != NULL);

    /* replace the old name with new name */
    NCI_Free(varp->name);
    varp->name     = nnewname;
    varp->name_len = nnewname_len;

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

