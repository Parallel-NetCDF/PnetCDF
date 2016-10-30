/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include "nc.h"
#include "ncx.h"
#include "fbits.h"
#include "macro.h"
#include "utf8proc.h"

/*
 * Free dim
 * Formerly
NC_free_dim(dim)
 */
inline void
ncmpii_free_NC_dim(NC_dim *dimp)
{
    if (dimp == NULL) return;
    ncmpii_free_NC_string(dimp->name);
    NCI_Free(dimp);
}


/* allocate and return a new NC_dim object */
inline NC_dim *
ncmpii_new_x_NC_dim(NC_string *name)
{
    NC_dim *dimp;

    dimp = (NC_dim *) NCI_Malloc(sizeof(NC_dim));
    if (dimp == NULL) return NULL;

    dimp->name = name;
    dimp->size = 0;

    return(dimp);
}

/*----< ncmpii_new_NC_dim() >------------------------------------------------*/
/*
 * Formerly, NC_new_dim(const char *name, long size)
 */
static int
ncmpii_new_NC_dim(NC_dimarray  *ncap,
                  const char   *name, /* normalized dim name */
                  MPI_Offset    size,
                  NC_dim      **dimp)
{
    NC_string *strp;

    if (strlen(name) == 0) DEBUG_RETURN_ERROR(NC_EBADNAME)

    strp = ncmpii_new_NC_string(strlen(name), name);
    if (strp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    *dimp = ncmpii_new_x_NC_dim(strp);
    if (*dimp == NULL) {
        ncmpii_free_NC_string(strp);
        DEBUG_RETURN_ERROR(NC_ENOMEM)
    }
    (*dimp)->size = size;

#ifndef SEARCH_NAME_LINEARLY
    if (ncap != NULL) { /* insert new dim to hash table */
        int key;
        NC_nametable *nameT = ncap->nameT; /* dim name lookup table */

        /* We use the first char as key for name lookup */
        key = HASH_FUNC(name);

        /* allocate or expand the space for nameT[key].list */
        if (nameT[key].num % NC_NAME_TABLE_CHUNK == 0)
            nameT[key].list = (int*) NCI_Realloc(nameT[key].list,
                              (size_t)(nameT[key].num+NC_NAME_TABLE_CHUNK) * sizeof(int));

        /* add the new variable ID to the name lookup table
         * the new varid will be ncap->ndefined
         */
        nameT[key].list[nameT[key].num] = ncap->ndefined;
        nameT[key].num++;
    }
    /* else case is for dimension duplication called from dup_NC_dim() */
#endif

    return NC_NOERR;
}

/*----< dup_NC_dim() >-------------------------------------------------------*/
NC_dim*
dup_NC_dim(const NC_dim *rdimp)
{
    int err;
    NC_dim *dimp;

    /* rdimp->name->cp is a normalized string */
    err = ncmpii_new_NC_dim(NULL, rdimp->name->cp, rdimp->size, &dimp);
    if (err != NC_NOERR) return NULL;
    return dimp;
}

/*----< ncmpii_find_NC_Udim() >----------------------------------------------*/
/*
 * Step thru NC_DIMENSION array, seeking the UNLIMITED dimension.
 * Return dimid or -1 on not found.
 * *dimpp is set to the appropriate NC_dim.
 */
int
ncmpii_find_NC_Udim(const NC_dimarray  *ncap,
                    NC_dim            **dimpp)
{
    int dimid;

    assert(ncap != NULL);

    if (ncap->ndefined == 0) return -1;

    /* note that the number of dimensions allowed is < 2^32 */
    for (dimid=0; dimid<ncap->ndefined; dimid++)
        if (ncap->value[dimid]->size == NC_UNLIMITED) {
            /* found the matched name */
            if (dimpp != NULL)
                *dimpp = ncap->value[dimid];
            return dimid;
        }

    /* not found */
    return -1;
}

#ifdef SEARCH_NAME_LINEARLY
/*----< ncmpii_NC_finddim() >------------------------------------------------*/
/*
 * Step thru NC_DIMENSION array, seeking match on name.
 * If found, set the dim ID pointed by dimidp, otherwise return NC_EBADDIM
 */
static int
ncmpii_NC_finddim(const NC_dimarray *ncap,
                  const char        *name,  /* normalized dim name */
                  int               *dimidp)
{
    int dimid;
    size_t nchars=strlen(name);

    assert(ncap != NULL);

    if (ncap->ndefined == 0) return NC_EBADDIM;

    /* note that the number of dimensions allowed is < 2^32 */
    for (dimid=0; dimid<ncap->ndefined; dimid++) {
        if (ncap->value[dimid]->name->nchars == nchars &&
            strncmp(ncap->value[dimid]->name->cp, name, nchars) == 0) {
            /* found the matched name */
            if (dimidp != NULL) *dimidp = dimid;
            return NC_NOERR; /* found it */
        }
    }
    return NC_EBADDIM; /* the name is not found */
}
#else
/*----< ncmpii_NC_finddim() >------------------------------------------------*/
/*
 * Search name from hash table ncap->nameT.
 * If found, set the dim ID pointed by dimidp, otherwise return NC_EBADDIM
 */
static int
ncmpii_NC_finddim(const NC_dimarray *ncap,
                  const char        *name,  /* normalized dim name */
                  int               *dimidp)
{
    int i, key, dimid;

    assert(ncap != NULL);

    if (ncap->ndefined == 0) return NC_EBADDIM;

    /* hash the dim name into a key for name lookup */
    key = HASH_FUNC(name);

    /* check the list using linear search */
    for (i=0; i<ncap->nameT[key].num; i++) {
        dimid = ncap->nameT[key].list[i];
        if (strcmp(name, ncap->value[dimid]->name->cp) == 0) {
            if (dimidp != NULL) *dimidp = dimid;
            return NC_NOERR; /* the name already exists */
        }
    }

    return NC_EBADDIM; /* the name has never been used */
}
#endif


/* dimarray */


/*----< ncmpii_free_NC_dimarray() >------------------------------------------*/
/*
 * Free NC_dimarray values.
 * formerly
NC_free_array()
 */
inline void
ncmpii_free_NC_dimarray(NC_dimarray *ncap)
{
    int i;

    assert(ncap != NULL);
    if (ncap->nalloc == 0) return;

    assert(ncap->value != NULL);
    for (i=0; i<ncap->ndefined; i++)
        ncmpii_free_NC_dim(ncap->value[i]);

    NCI_Free(ncap->value);
    ncap->value    = NULL;
    ncap->nalloc   = 0;
    ncap->ndefined = 0;

    /* free space allocated for dim name lookup table */
    for (i=0; i<HASH_TABLE_SIZE; i++) {
        if (ncap->nameT[i].num > 0)
            NCI_Free(ncap->nameT[i].list);
        ncap->nameT[i].num = 0;
    }
}


/*----< ncmpii_dup_NC_dimarray() >-------------------------------------------*/
int
ncmpii_dup_NC_dimarray(NC_dimarray *ncap, const NC_dimarray *ref)
{
    int i, status=NC_NOERR;

    assert(ref != NULL);
    assert(ncap != NULL);

    if (ref->nalloc == 0) {
        ncap->nalloc   = 0;
        ncap->ndefined = 0;
        ncap->value    = NULL;
        return NC_NOERR;
    }

    if (ref->nalloc > 0) {
        ncap->value = (NC_dim **) NCI_Calloc((size_t)ref->nalloc, sizeof(NC_dim *));
        if (ncap->value == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
        ncap->nalloc = ref->nalloc;
    }

    ncap->ndefined = 0;
    for (i=0; i<ref->ndefined; i++) {
        ncap->value[i] = dup_NC_dim(ref->value[i]);
        if (ncap->value[i] == NULL) {
            DEBUG_ASSIGN_ERROR(status, NC_ENOMEM)
            break;
        }
    }

    if (status != NC_NOERR) {
        ncmpii_free_NC_dimarray(ncap);
        return status;
    }

    ncap->ndefined = ref->ndefined;

    /* duplicate dim name lookup table */
    for (i=0; i<HASH_TABLE_SIZE; i++) {
        ncap->nameT[i].num = ref->nameT[i].num;
        ncap->nameT[i].list = NULL;
        if (ncap->nameT[i].num > 0) {
            ncap->nameT[i].list = NCI_Malloc((size_t)ncap->nameT[i].num * sizeof(int));
            memcpy(ncap->nameT[i].list, ref->nameT[i].list,
                   (size_t)ncap->nameT[i].num * sizeof(int));
        }
    }

    return NC_NOERR;
}


/*----< incr_NC_dimarray() >---------------------------------------------- --*/
/*
 * Add a new handle to the end of an array of handles
 * Formerly, NC_incr_array(array, tail)
 */
int
incr_NC_dimarray(NC_dimarray *ncap,
                 NC_dim      *newdimp)
{
    NC_dim **vp;

    assert(ncap != NULL);

    if (ncap->nalloc == 0) {
        assert(ncap->ndefined == 0);
        vp = (NC_dim **) NCI_Malloc(NC_ARRAY_GROWBY * sizeof(NC_dim *));
        if (vp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

        ncap->value = vp;
        ncap->nalloc = NC_ARRAY_GROWBY;
    }
    else if (ncap->ndefined + 1 > ncap->nalloc) {
        vp = (NC_dim **) NCI_Realloc(ncap->value,
             (size_t)(ncap->nalloc + NC_ARRAY_GROWBY) * sizeof(NC_dim *));
        if (vp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

        ncap->value = vp;
        ncap->nalloc += NC_ARRAY_GROWBY;
    }
    /* else here means some space still available */

    if (newdimp != NULL) {
        ncap->value[ncap->ndefined] = newdimp;
        ncap->ndefined++;
    }

    return NC_NOERR;
}


/*----< ncmpii_elem_NC_dimarray() >------------------------------------------*/
inline NC_dim *
ncmpii_elem_NC_dimarray(const NC_dimarray *ncap,
                        int                dimid)
{
    /* returns the dimension ID defined earlier */
    assert(ncap != NULL);

    if (dimid < 0 || ncap->ndefined == 0 || dimid >= ncap->ndefined)
        return NULL;

    assert(ncap->value != NULL);

    return ncap->value[dimid];
}


/* Public */

/*----< ncmpi_def_dim() >----------------------------------------------------*/
int
ncmpi_def_dim(int         ncid,    /* IN:  file ID */
              const char *name,    /* IN:  name of dimension */
              MPI_Offset  size,    /* IN:  dimension size */
              int        *dimidp)  /* OUT: dimension ID */
{
    int dimid, err;
    char *nname=NULL;  /* normalized name */
    NC *ncp=NULL;
    NC_dim *dimp=NULL;

    /* check if ncid is valid */
    err = ncmpii_NC_check_id(ncid, &ncp);
    if (err != NC_NOERR || ncp == NULL) DEBUG_RETURN_ERROR(err)

    /* must be called in define mode */
    if (!NC_indef(ncp)) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

    if (name == NULL || *name == 0 || strlen(name) > NC_MAX_NAME) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    /* check if the name string is legal for the netcdf format */
    err = ncmpii_NC_check_name(name, ncp->format);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR
        goto err_check;
    }

    /* MPI_Offset is usually a signed value, but serial netcdf uses
     * size_t -- normally unsigned
     * In 1999 ISO C standard, size_t is an unsigned integer type of at least
     * 16 bit. */
    if (ncp->format == 2) { /* CDF-2 format, max is 2^32-4 */
        if (size > X_UINT_MAX - 3 || (size < 0))
            /* "-3" handles rounded-up size */
            err = NC_EDIMSIZE;
    } else if (ncp->format == 5) { /* CDF-5 format*/
        if (size < 0)
            err = NC_EDIMSIZE;
    } else { /* CDF-1 format, max is 2^31-4 */
        if (size > X_INT_MAX - 3 || (size < 0))
            /* "-3" handles rounded-up size */
            err = NC_EDIMSIZE;
    }
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR
        goto err_check;
    }

    if (size == NC_UNLIMITED) {
#if 0
        /* check for any existing unlimited dimension, netcdf allows
         * one per file
         */
        dimid = ncmpii_find_NC_Udim(&ncp->dims, &dimp);
        if (dimid != -1) {
            DEBUG_ASSIGN_ERROR(err, NC_EUNLIMIT) /* found an existing one */
            goto err_check;
        }
#endif
        if (ncp->dims.unlimited_id != -1) {
            DEBUG_ASSIGN_ERROR(err, NC_EUNLIMIT) /* already defined */
            goto err_check;
        }
    }

    /* check if number of dimensions exceeds the upper bound */
    if (ncp->dims.ndefined >= NC_MAX_DIMS) {
        DEBUG_ASSIGN_ERROR(err, NC_EMAXDIMS)
        goto err_check;
    }

    /* create a normalized character string */
    nname = (char *)ncmpii_utf8proc_NFC((const unsigned char *)name);
    if (nname == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto err_check;
    }

    /* check if the name string is previously used */
    err = ncmpii_NC_finddim(&ncp->dims, nname, NULL);
    if (err != NC_EBADDIM) {
        DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
        goto err_check;
    }
    else
        err = NC_NOERR;

err_check:
    if (ncp->safe_mode) {
        int status, mpireturn;
        char root_name[NC_MAX_NAME];
        MPI_Offset root_size;

        /* check if name is consistent among all processes */
        if (name == NULL || *name == 0)
            root_name[0] = 0;
        else
            strncpy(root_name, name, NC_MAX_NAME);
        TRACE_COMM(MPI_Bcast)(root_name, NC_MAX_NAME, MPI_CHAR, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nname != NULL) free(nname);
            return ncmpii_handle_error(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && strcmp(root_name, name))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_NAME)

        /* check if sizes are consistent across all processes */
        root_size = size;
        TRACE_COMM(MPI_Bcast)(&root_size, 1, MPI_OFFSET, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nname != NULL) free(nname);
            return ncmpii_handle_error(mpireturn, "MPI_Bcast");
        }
        if (err == NC_NOERR && root_size != size)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_SIZE)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS) {
            if (nname != NULL) free(nname);
            return ncmpii_handle_error(mpireturn, "MPI_Allreduce");
        }
        if (err == NC_NOERR) err = status;
    }

    if (err != NC_NOERR) {
        if (nname != NULL) free(nname);
        return err;
    }

    assert(nname != NULL);

    /* create a new dimension object */
    err = ncmpii_new_NC_dim(&ncp->dims, nname, size, &dimp);
    free(nname);
    if (err != NC_NOERR) {
        if (dimp != NULL) ncmpii_free_NC_dim(dimp);
        DEBUG_RETURN_ERROR(err)
    }

    /* Add a new dim handle to the end of handle array */
    err = incr_NC_dimarray(&ncp->dims, dimp);
    if (err != NC_NOERR) {
        if (dimp != NULL) ncmpii_free_NC_dim(dimp);
        DEBUG_RETURN_ERROR(err)
    }

    /* ncp->dims.ndefined has been increased in incr_NC_dimarray() */
    dimid = (int)ncp->dims.ndefined -1;

    if (size == NC_UNLIMITED) ncp->dims.unlimited_id = dimid;

    if (dimidp != NULL) *dimidp = dimid;


    return err;
}


/*----< ncmpi_inq_dimid() >--------------------------------------------------*/
int
ncmpi_inq_dimid(int         ncid,
                const char *name,
                int        *dimid)
{
    int err;
    char *nname=NULL; /* normalized name */
    NC *ncp=NULL;

    err = ncmpii_NC_check_id(ncid, &ncp);
    if (err != NC_NOERR || ncp == NULL) DEBUG_RETURN_ERROR(err)

    if (name == NULL || *name == 0 || strlen(name) > NC_MAX_NAME)
        DEBUG_RETURN_ERROR(NC_EBADNAME)

    /* create a normalized character string */
    nname = (char *)ncmpii_utf8proc_NFC((const unsigned char *)name);
    if (nname == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    err = ncmpii_NC_finddim(&ncp->dims, nname, dimid);
    free(nname);

    return err;
}


/*----< ncmpi_inq_dim() >----------------------------------------------------*/
int
ncmpi_inq_dim(int         ncid,
              int         dimid,
              char       *name,
              MPI_Offset *sizep)
{
    int status;
    NC *ncp=NULL;
    NC_dim *dimp=NULL;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR || ncp == NULL) DEBUG_RETURN_ERROR(status)

    dimp = ncmpii_elem_NC_dimarray(&ncp->dims, dimid);
    if (dimp == NULL) DEBUG_RETURN_ERROR(NC_EBADDIM)

    if (name != NULL)
        /* in PnetCDF, name->cp is always NULL character terminated */
        strcpy(name, dimp->name->cp);

    if (sizep != NULL) {
        if (dimp->size == NC_UNLIMITED)
            *sizep = NC_get_numrecs(ncp);
        else
            *sizep = dimp->size;
    }
    return NC_NOERR;
}


/*----< ncmpi_inq_dimname() >------------------------------------------------*/
int
ncmpi_inq_dimname(int   ncid,
                  int   dimid,
                  char *name)
{
    return ncmpi_inq_dim(ncid, dimid, name, NULL);
}


/*----< ncmpi_inq_dimlen() >-------------------------------------------------*/
int
ncmpi_inq_dimlen(int         ncid,
                 int         dimid,
                 MPI_Offset *lenp)
{
    return ncmpi_inq_dim(ncid, dimid, NULL, lenp);
}


/*----< ncmpi_rename_dim() >--------------------------------------------------*/
/* This API is collective and can be called in either define or data mode..
 * If the new name is longer than the old name, the netCDF dataset must be in
 * the define mode.
 */
int
ncmpi_rename_dim(int         ncid,
                 int         dimid,
                 const char *newname)
{
    int err;
    char *nnewname=NULL; /* normalized newname */
    NC *ncp=NULL;
    NC_dim *dimp=NULL;
    NC_string *newStr=NULL;

    /* check whether ncid is valid */
    err = ncmpii_NC_check_id(ncid, &ncp);
    if (err != NC_NOERR || ncp == NULL) DEBUG_RETURN_ERROR(err)

    /* check file's write permission */
    if (NC_readonly(ncp)) {
        DEBUG_ASSIGN_ERROR(err, NC_EPERM)
        goto err_check;
    }

    if (newname == NULL || *newname == 0 || strlen(newname) > NC_MAX_NAME) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADNAME)
        goto err_check;
    }

    /* check whether newname is legal */
    err = ncmpii_NC_check_name(newname, ncp->format);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR
        goto err_check;
    }

    /* create a normalized character string */
    nnewname = (char *)ncmpii_utf8proc_NFC((const unsigned char *)newname);
    if (nnewname == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto err_check;
    }

    /* check whether newname is already in use */
    err = ncmpii_NC_finddim(&ncp->dims, nnewname, NULL);
    if (err != NC_EBADDIM) {
        DEBUG_ASSIGN_ERROR(err, NC_ENAMEINUSE)
        goto err_check;
    }

    /* retrieve dim object */
    dimp = ncmpii_elem_NC_dimarray(&ncp->dims, dimid);
    if (dimp == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_EBADDIM)
        goto err_check;
    }

    if (! NC_indef(ncp) && /* when file is in data mode */
        dimp->name->nchars < (MPI_Offset)strlen(nnewname)) {
        /* must in define mode when newname is longer */
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

    newStr = ncmpii_new_NC_string(strlen(nnewname), nnewname);
    if (newStr == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
        goto err_check;
    }

#ifndef SEARCH_NAME_LINEARLY
    /* update dim name lookup table, by removing the old name and add
     * the new name */
    err = ncmpii_update_name_lookup_table(ncp->dims.nameT, dimid,
          ncp->dims.value[dimid]->name->cp, nnewname);
    if (err != NC_NOERR) {
        DEBUG_TRACE_ERROR
        goto err_check;
    }
#endif

err_check:
    if (nnewname != NULL) free(nnewname);

    if (ncp->safe_mode) {
        int root_dimid, status, mpireturn;
        char root_name[NC_MAX_NAME];

        /* check if newname is consistent among all processes */
        if (newname == NULL || *newname == 0)
            root_name[0] = 0;
        else
            strncpy(root_name, newname, NC_MAX_NAME);
        TRACE_COMM(MPI_Bcast)(root_name, NC_MAX_NAME, MPI_CHAR, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && strcmp(root_name, newname))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_NAME)

        /* check if dimid is consistent across all processes */
        root_dimid = dimid;
        TRACE_COMM(MPI_Bcast)(&root_dimid, 1, MPI_INT, 0, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_dimid != dimid)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_handle_error(mpireturn, "MPI_Allreduce");
        if (err == NC_NOERR) err = status;
    }

    if (err != NC_NOERR) {
        if (newStr != NULL) ncmpii_free_NC_string(newStr);
        return err;
    }

    /* replace the old name with new name */
    assert(dimp != NULL);
    ncmpii_free_NC_string(dimp->name);
    dimp->name = newStr;

    if (! NC_indef(ncp)) { /* when file is in data mode */
        /* Let root write the entire header to the file. Note that we cannot
         * just update the variable name in its space occupied in the file
         * header, because if the file space occupied by the name shrinks, all
         * the metadata following it must be moved ahead.
         */
        err = ncmpii_write_header(ncp);
        if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
    }

    return err;
}
