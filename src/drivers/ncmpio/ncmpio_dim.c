/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in
 * src/dispatchers/dimension.c
 *
 * ncmpi_def_dim()    : dispatcher->def_dim()
 * ncmpi_inq_dimid()  : dispatcher->inq_dimid()
 * ncmpi_inq_dim()    : dispatcher->inq_dim()
 * ncmpi_rename_dim() : dispatcher->rename_dim()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"
#include <ncx.h>

/*----< dup_NC_dim() >-------------------------------------------------------*/
static int
dup_NC_dim(const NC_dim *rdimp, NC_dim **dimp)
{
    *dimp = (NC_dim*) NCI_Malloc(sizeof(NC_dim));
    if (*dimp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    (*dimp)->size     = rdimp->size;
    (*dimp)->name_len = strlen(rdimp->name)+1;
    (*dimp)->name     = (char*) NCI_Malloc((*dimp)->name_len);
    if ((*dimp)->name == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    strcpy((*dimp)->name, rdimp->name);

    return NC_NOERR;
}

#ifdef SEARCH_NAME_LINEARLY
/*----< NC_finddim() >-------------------------------------------------------*/
/*
 * Step thru NC_DIMENSION array, seeking match on name.
 * If found, set the dim ID pointed by dimidp, otherwise return NC_EBADDIM
 */
static int
NC_finddim(const NC_dimarray *ncap,
           const char        *name,  /* normalized dim name */
           int               *dimidp)
{
    int dimid;
    size_t nchars;

    if (ncap->ndefined == 0) return NC_EBADDIM;

    /* note that the number of dimensions allowed is < 2^32 */
    nchars = strlen(name);
    for (dimid=0; dimid<ncap->ndefined; dimid++) {
        if (ncap->value[dimid]->name_len == nchars &&
            strcmp(ncap->value[dimid]->name, name) == 0) {
            /* found the matched name */
            if (dimidp != NULL) *dimidp = dimid;
            return NC_NOERR; /* found it */
        }
    }
    return NC_EBADDIM; /* the name is not found */
}
#else
/*----< NC_finddim() >-------------------------------------------------------*/
/*
 * Search name from hash table ncap->nameT.
 * If found, set the dim ID pointed by dimidp, otherwise return NC_EBADDIM
 */
static int
NC_finddim(const NC_dimarray *ncap,
           const char        *name,  /* normalized dim name */
           int               *dimidp)
{
    int i, key, dimid;
    size_t nchars;

    if (ncap->ndefined == 0) return NC_EBADDIM;

    /* hash the dim name into a key for name lookup */
    key = HASH_FUNC(name);

    /* check the list using linear search */
    nchars = strlen(name);
    for (i=0; i<ncap->nameT[key].num; i++) {
        dimid = ncap->nameT[key].list[i];
        if (ncap->value[dimid]->name_len == nchars &&
            strcmp(name, ncap->value[dimid]->name) == 0) {
            if (dimidp != NULL) *dimidp = dimid;
            return NC_NOERR; /* the name already exists */
        }
    }
    return NC_EBADDIM; /* the name has never been used */
}
#endif

/* dimarray */

/*----< ncmpio_free_NC_dimarray() >------------------------------------------*/
void
ncmpio_free_NC_dimarray(NC_dimarray *ncap)
{
    int i;

    assert(ncap != NULL);
    if (ncap->ndefined == 0) return;

    if (ncap->value != NULL) {
        /* when error is detected reading NC_DIMENSION tag, ncap->ndefined can
         * be > 0 and ncap->value is still NULL
         */
        for (i=0; i<ncap->ndefined; i++) {
            /* when error is detected reading dimension i, ncap->value[i] can
             * still be NULL
             */
            if (ncap->value[i] == NULL) break;
            NCI_Free(ncap->value[i]->name);
            NCI_Free(ncap->value[i]);
        }
        NCI_Free(ncap->value);
        ncap->value = NULL;
    }
    ncap->ndefined = 0;

#ifndef SEARCH_NAME_LINEARLY
    /* free space allocated for dim name lookup table */
    ncmpio_hash_table_free(ncap->nameT);
#endif
}

/*----< ncmpio_dup_NC_dimarray() >-------------------------------------------*/
int
ncmpio_dup_NC_dimarray(NC_dimarray *ncap, const NC_dimarray *ref)
{
    int i, status=NC_NOERR;

    assert(ref != NULL);
    assert(ncap != NULL);

    if (ref->ndefined == 0) {
        ncap->ndefined = 0;
        ncap->value    = NULL;
        return NC_NOERR;
    }

    /* allocate array of NC_dim objects */
    if (ref->ndefined > 0) {
        size_t alloc_size = _RNDUP(ref->ndefined, NC_ARRAY_GROWBY);
        ncap->value = (NC_dim**) NCI_Calloc(alloc_size, sizeof(NC_dim*));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    }

    /* duplicate each NC_dim objects */
    ncap->ndefined = 0;
    for (i=0; i<ref->ndefined; i++) {
        status = dup_NC_dim(ref->value[i], &ncap->value[i]);
        if (status != NC_NOERR) {
            ncmpio_free_NC_dimarray(ncap);
            return status;
        }
        ncap->ndefined++;
    }
    assert(ncap->ndefined == ref->ndefined);

#ifndef SEARCH_NAME_LINEARLY
    /* duplicate dim name lookup table */
    ncmpio_hash_table_copy(ncap->nameT, ref->nameT);
#endif

    return NC_NOERR;
}

/*----< ncmpio_def_dim() >---------------------------------------------------*/
int
ncmpio_def_dim(void       *ncdp,    /* IN:  NC object */
               const char *name,    /* IN:  name of dimension */
               MPI_Offset  size,    /* IN:  dimension size */
               int        *dimidp)  /* OUT: dimension ID */
{
    int dimid, err=NC_NOERR;
    char *nname=NULL;  /* normalized name */
    NC *ncp=(NC*)ncdp;
    NC_dim *dimp=NULL;

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) return err;

    /* create a new dimension object (dimp->name points to nname) */
    dimp = (NC_dim*) NCI_Malloc(sizeof(NC_dim));
    if (dimp == NULL) {
        NCI_Free(nname);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    dimp->size     = size;
    dimp->name     = nname;
    dimp->name_len = strlen(nname);

    /* allocate/expand ncp->dims.value array */
    if (ncp->dims.ndefined % NC_ARRAY_GROWBY == 0) {
        size_t alloc_size = (size_t)ncp->dims.ndefined + NC_ARRAY_GROWBY;

        ncp->dims.value = (NC_dim **) NCI_Realloc(ncp->dims.value,
                                      alloc_size * sizeof(NC_dim*));
        if (ncp->dims.value == NULL) {
            NCI_Free(nname);
            NCI_Free(dimp);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
    }

    dimid = ncp->dims.ndefined;

    /* Add a new dim handle to the end of handle array */
    ncp->dims.value[dimid] = dimp;

    if (size == NC_UNLIMITED) ncp->dims.unlimited_id = dimid;

    ncp->dims.ndefined++;

#ifndef SEARCH_NAME_LINEARLY
    ncmpio_hash_insert(ncp->dims.nameT, nname, dimid);
#endif

    if (dimidp != NULL) *dimidp = dimid;

    return err;
}

/*----< ncmpio_inq_dimid() >-------------------------------------------------*/
int
ncmpio_inq_dimid(void       *ncdp,
                 const char *name,
                 int        *dimid)
{
    int err=NC_NOERR;
    char *nname=NULL; /* normalized name */
    NC *ncp=(NC*)ncdp;

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(name, &nname);
    if (err != NC_NOERR) return err;

    err = NC_finddim(&ncp->dims, nname, dimid);
    NCI_Free(nname);

    return err;
}

/*----< ncmpio_inq_dim() >---------------------------------------------------*/
int
ncmpio_inq_dim(void       *ncdp,
               int         dimid,
               char       *name,
               MPI_Offset *sizep)
{
    NC_dim *dimp;
    NC *ncp=(NC*)ncdp;

    /* sanity check for dimid has been done at dispatchers */
    dimp = ncp->dims.value[dimid];

    if (name != NULL)
        /* in PnetCDF, name is always NULL character terminated */
        strcpy(name, dimp->name);

    if (sizep != NULL) {
        if (dimp->size == NC_UNLIMITED)
            *sizep = ncp->numrecs;
        else
            *sizep = dimp->size;
    }
    return NC_NOERR;
}

/*----< ncmpio_rename_dim() >-------------------------------------------------*/
/* This API is collective and can be called in either define or data mode..
 * If the new name is longer than the old name, the netCDF dataset must be in
 * the define mode.
 */
int
ncmpio_rename_dim(void       *ncdp,
                  int         dimid,
                  const char *newname)
{
    int err=NC_NOERR;
    char *nnewname=NULL; /* normalized newname */
    size_t nnewname_len=0;
    NC *ncp=(NC*)ncdp;
    NC_dim *dimp=NULL;

    /* create a normalized character string */
    err = ncmpii_utf8_normalize(newname, &nnewname);
    if (err != NC_NOERR) goto err_check;

    nnewname_len = strlen(nnewname);

    /* sanity check for dimid has been done at dispatchers */
    dimp = ncp->dims.value[dimid];

    if (! NC_indef(ncp) && dimp->name_len < nnewname_len) {
        /* when in data mode, newname cannot be longer than the old one */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

#ifndef SEARCH_NAME_LINEARLY
    /* update dim name lookup table, by removing the old name and add
     * the new name */
    err = ncmpio_update_name_lookup_table(ncp->dims.nameT, dimid,
                             ncp->dims.value[dimid]->name, nnewname);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR(err)
        goto err_check;
    }
#endif

err_check:
    if (ncp->safe_mode) {
        int status, mpireturn;

        /* check the error so far across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            NCI_Free(nnewname);
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        }
        if (status != NC_NOERR) {
            NCI_Free(nnewname);
            return status;
        }
    }

    if (err != NC_NOERR) {
        if (nnewname != NULL) NCI_Free(nnewname);
        return err;
    }

    /* replace the old name with new name */
    assert(dimp != NULL);
    NCI_Free(dimp->name);
    dimp->name     = nnewname;
    dimp->name_len = nnewname_len;

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
