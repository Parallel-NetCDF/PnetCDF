/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 * and src/dispatchers/variable.c
 *
 * ncmpi_set_fill()     : dispatcher->set_fill()
 * ncmpi_fill_var_rec() : dispatcher->fill_rec()
 * ncmpi_def_var_fill() : dispatcher->def_var_fill()
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
#include <ncx.h>
#include "ncmpio_NC.h"

#define CHECK_ERROR(status) {                                                \
    if (ncp->safe_mode == 1) {                                               \
        int g_status;                                                        \
        TRACE_COMM(MPI_Allreduce)(&status, &g_status, 1, MPI_INT, MPI_MIN,   \
                                  ncp->comm);                                \
        if (g_status != NC_NOERR) return g_status;                           \
    }                                                                        \
    else if (status != NC_NOERR)                                             \
        return status;                                                       \
}

/* The default fill values defined in pnetcdf.h.inc must be the same as the
 * ones defined in netCDF-4 and match the hexadecimal values set below
 *
#define NC_FILL_BYTE    ((signed char)-127)
#define NC_FILL_CHAR    ((char)0)
#define NC_FILL_SHORT   ((short)-32767)
#define NC_FILL_INT     (-2147483647)
#define NC_FILL_FLOAT   (9.9692099683868690e+36f)
#define NC_FILL_DOUBLE  (9.9692099683868690e+36)
#define NC_FILL_UBYTE   (255)
#define NC_FILL_USHORT  (65535)
#define NC_FILL_UINT    (4294967295U)
#define NC_FILL_INT64   ((long long)-9223372036854775806LL)
#define NC_FILL_UINT64  ((unsigned long long)18446744073709551614ULL)
*/
static unsigned char FILL_CHAR[1]   = {0x00};
static unsigned char FILL_BYTE[1]   = {0x81};
static unsigned char FILL_SHORT[2]  = {0x80, 0x01};
static unsigned char FILL_INT[4]    = {0x80, 0x00, 0x00, 0x01};
static unsigned char FILL_FLOAT[4]  = {0x7C, 0xF0, 0x00, 0x00};
static unsigned char FILL_DOUBLE[8] = {0x47, 0x9E, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
static unsigned char FILL_UBYTE[1]  = {0xFF};
static unsigned char FILL_USHORT[2] = {0xFF, 0xFF};
static unsigned char FILL_UINT[4]   = {0xFF, 0xFF, 0xFF, 0xFF};
static unsigned char FILL_INT64[8]  = {0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02};
static unsigned char FILL_UINT64[8] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFE};

/*----< ncmpio_inq_default_fill_value() >------------------------------------*/
/* copy the default fill value to the memory space pointed by fill_value */
int
ncmpio_inq_default_fill_value(int xtype, void *fillp)
{
    if (fillp == NULL) return NC_NOERR;

    switch(xtype) {
        case NC_CHAR   :               *(char*)fillp = NC_FILL_CHAR;   break;
        case NC_BYTE   :        *(signed char*)fillp = NC_FILL_BYTE;   break;
        case NC_SHORT  :              *(short*)fillp = NC_FILL_SHORT;  break;
        case NC_INT    :                *(int*)fillp = NC_FILL_INT;    break;
        case NC_FLOAT  :              *(float*)fillp = NC_FILL_FLOAT;  break;
        case NC_DOUBLE :             *(double*)fillp = NC_FILL_DOUBLE; break;
        case NC_UBYTE  :      *(unsigned char*)fillp = NC_FILL_UBYTE;  break;
        case NC_USHORT :     *(unsigned short*)fillp = NC_FILL_USHORT; break;
        case NC_UINT   :       *(unsigned int*)fillp = NC_FILL_UINT;   break;
        case NC_INT64  :          *(long long*)fillp = NC_FILL_INT64;  break;
        case NC_UINT64 : *(unsigned long long*)fillp = NC_FILL_UINT64; break;
        default : DEBUG_RETURN_ERROR(NC_EBADTYPE)
    }
    return NC_NOERR;
}

/*----< fill_var_buf() >----------------------------------------------------*/
/* fill the buffer, buf, with either user-defined fill values or default
 * values */
static int
fill_var_buf(const NC_var *varp,
             MPI_Offset    bnelems, /* number of elements in buf */
             void         *buf)
{
    int i, indx;

    indx = ncmpio_NC_findattr(&varp->attrs, _FillValue);
    if (indx >= 0) {
        /* User defined fill value */
        NC_attr *attrp = varp->attrs.value[indx];
        if (attrp->xtype != varp->xtype || attrp->nelems != 1)
            DEBUG_RETURN_ERROR(NC_EBADTYPE)

        /* Use the user defined value */
        char *bufp = buf;
        for (i=0; i<bnelems; i++) {
            memcpy(bufp, attrp->xvalue, (size_t)varp->xsz);
            bufp += varp->xsz;
        }
    }
    else { /* use the default */
        void *xvalue;
        switch(varp->xtype) {
            case NC_CHAR   : xvalue = &FILL_CHAR[0];   break;
            case NC_BYTE   : xvalue = &FILL_BYTE[0];   break;
            case NC_SHORT  : xvalue = &FILL_SHORT[0];  break;
            case NC_INT    : xvalue = &FILL_INT[0];    break;
            case NC_FLOAT  : xvalue = &FILL_FLOAT[0];  break;
            case NC_DOUBLE : xvalue = &FILL_DOUBLE[0]; break;
            case NC_UBYTE  : xvalue = &FILL_UBYTE[0];  break;
            case NC_USHORT : xvalue = &FILL_USHORT[0]; break;
            case NC_UINT   : xvalue = &FILL_UINT[0];   break;
            case NC_INT64  : xvalue = &FILL_INT64[0];  break;
            case NC_UINT64 : xvalue = &FILL_UINT64[0]; break;
            default : DEBUG_RETURN_ERROR(NC_EBADTYPE)
        }

        char *bufp = buf;
        for (i=0; i<bnelems; i++) {
            memcpy(bufp, xvalue, (size_t)varp->xsz);
            bufp += varp->xsz;
        }
    }
    return NC_NOERR;
}

/*----< fill_var_rec() >-----------------------------------------------------*/
/* This function is a collective.
 * For fixed-size variables, write the entire variable with fill values and
 * ignore argument recno. For record variables, write one record of that
 * variable with pre-defined/supplied fill value.
 */
static int
fill_var_rec(NC         *ncp,
             NC_var     *varp,
             MPI_Offset  recno) /* record number */
{
    int err, mpireturn, rank, nprocs;
    void *buf;
    MPI_Offset var_len, start, count, offset;
    MPI_File fh;
    MPI_Status mpistatus;

    MPI_Comm_rank(ncp->comm, &rank);
    MPI_Comm_size(ncp->comm, &nprocs);

    if (varp->ndims == 0) /* scalar variable */
        var_len = 1;
    else if (varp->ndims == 1 && IS_RECVAR(varp))
        var_len = 1;
    else if (IS_RECVAR(varp))
        var_len = varp->dsizes[1];
    else
        var_len = varp->dsizes[0];

    /* divide total number of elements of this variable among all processes */
    count = var_len / nprocs;
    start = count * rank;
    if (rank < var_len % nprocs) {
        start += rank;
        count++;
    }
    else {
        start += var_len % nprocs;
    }

    /* allocate buffer space */
    buf = NCI_Malloc((size_t)(count * varp->xsz));

    /* fill buffer with fill values */
    err = fill_var_buf(varp, count, buf);
    if (err != NC_NOERR) {
        NCI_Free(buf);
        return err;
    }

    /* calculate the starting file offset for each process */
    offset = varp->begin;
    if (IS_RECVAR(varp))
        offset += ncp->recsize * recno;
    offset += start * varp->xsz;

    fh = ncp->collective_fh;

    /* make the entire file visible */
    TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                MPI_INFO_NULL);

    count *= varp->xsz;
    if (count != (int)count) DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
    if (err != NC_NOERR)
        count = 0; /* participate collective write with 0-length request */

    /* write to variable collectively */
    TRACE_IO(MPI_File_write_at_all)(fh, offset, buf, (int)count, MPI_BYTE,
                                    &mpistatus);
    NCI_Free(buf);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_File_write_at_all");

    if (err != NC_NOERR) return err;

    if (IS_RECVAR(varp)) { /* update header's number of records in memory */
        /* recno may be differenet among, if safe mode is disabled. In
         * addition, recno can be > or < then ncp->numrecs. We need to sync
         * first before update numrecs field in file.
         *
         * First, find the max numrecs among all processes.
         */
        MPI_Offset max_numrecs, numrecs=recno+1;
        TRACE_COMM(MPI_Allreduce)(&numrecs, &max_numrecs, 1, MPI_OFFSET,
                                  MPI_MAX, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");

        /* In collective mode, ncp->numrecs is always sync-ed among processes */
        if (ncp->numrecs < max_numrecs) {
            err = ncmpio_write_numrecs(ncp, max_numrecs);
            if (err != NC_NOERR) return err;
            ncp->numrecs = max_numrecs;
        }
    }

    return NC_NOERR;
}

#ifdef FILL_ONE_VAR_AT_A_TIME
/*----< fillerup() >---------------------------------------------------------*/
/* fill the newly created fixed-size variables */
static int
fillerup(NC *ncp)
{
    int i, indx, err=NC_NOERR;

    /* loop thru all variables */
    for (i=0; i<ncp->vars.ndefined; i++) {
        if (IS_RECVAR(ncp->vars.value[i]))
            /* skip record variables. In PnetCDF record variables mus be
             * explicitly filled by calling ncmpi_fill_var_rec() */
            continue;

        /* check if _FillValue attribute is defined */
        indx = ncmpio_NC_findattr(&ncp->vars.value[i]->attrs, _FillValue);

        /* only if filling this variable is requested. Fill mode can be
         * enabled by 2 ways: explictly call to ncmpi_def_var_fill() or put
         * the attribute named _FillValue */
        if (ncp->vars.value[i]->no_fill && indx == -1) continue;

        /* collectively fill the entire variable */
        err = fill_var_rec(ncp, ncp->vars.value[i], 0);
        if (err != NC_NOERR) break;
    }
    return err;
}

/*----< fill_added() >-------------------------------------------------------*/
/* fill the newly added variables */
static int
fill_added(NC *ncp, NC *old_ncp)
{
    int indx, err=NC_NOERR, varid;

    /* loop thru all new variables */
    varid = old_ncp->vars.ndefined;
    for (; varid<ncp->vars.ndefined; varid++) {
        if (IS_RECVAR(ncp->vars.value[varid]))
            /* skip record variables */
            continue;

        /* check if _FillValue attribute is defined */
        indx = ncmpio_NC_findattr(&ncp->vars.value[varid]->attrs, _FillValue);

        /* only if filling this variable is requested */
        if (ncp->vars.value[varid]->no_fill && indx == -1) continue;

        /* collectively fill the entire variable */
        err = fill_var_rec(ncp, ncp->vars.value[varid], 0);
        if (err != NC_NOERR) break;
    }
    return err;
}

/*----< fill_added_recs() >--------------------------------------------------*/
/* for each newly added record variable, we fill the records one at a time */
static int
fill_added_recs(NC *ncp, NC *old_ncp)
{
    MPI_Offset old_nrecs = old_ncp->numrecs;
    int indx, err=NC_NOERR, recno, varid;

    /* loop thru all old records */
    for (recno=0; recno<old_nrecs; recno++) {
        /* check newly added variables only */
        for (varid=old_ncp->vars.ndefined; varid<ncp->vars.ndefined; varid++) {
            if (!IS_RECVAR(ncp->vars.value[varid]))
                /* skip non-record variables */
                continue;

            /* check if _FillValue attribute is defined */
            indx = ncmpio_NC_findattr(&ncp->vars.value[varid]->attrs, _FillValue);

            /* only if filling this variable is requested */
            if (ncp->vars.value[varid]->no_fill && indx == -1) continue;

            /* collectively fill the record */
            err = fill_var_rec(ncp, ncp->vars.value[varid], recno);
            if (err != NC_NOERR) return err;
        }
    }
    return NC_NOERR;
}
#endif

/*----< fillerup_aggregate() >-----------------------------------------------*/
/* fill the newly created/added variables
 * This version aggregates all writes into a single one
 */
static int
fillerup_aggregate(NC *ncp, NC *old_ncp)
{
    int i, j, k, rank, nprocs, start_vid, recno;
    int nVarsFill, *blocklengths;
    int mpireturn, err, status=NC_NOERR;
    char *buf_ptr, *noFill;
    void *buf;
    size_t nsegs;
    MPI_Offset buf_len, var_len, nrecs, start, *count;
    MPI_Aint *offset;
    MPI_Datatype filetype;
    MPI_File fh;
    MPI_Status mpistatus;
    NC_var *varp;

    MPI_Comm_rank(ncp->comm, &rank);
    MPI_Comm_size(ncp->comm, &nprocs);

    /* find the starting vid for newly added variables */
    start_vid = 0;
    nrecs = 0;  /* the current number of records */
    if (old_ncp != NULL) {
        start_vid = old_ncp->vars.ndefined;
        nrecs = old_ncp->numrecs;
    }

    noFill = (char*) NCI_Malloc((size_t)(ncp->vars.ndefined - start_vid));
    nVarsFill = 0;

#ifdef _CHECK_FILL_MODE_CONSISTENCY
    /* Because fill mode is not part of file header, we must broadcast root's
     * variables' fill modes and overwrite local's if an inconsistency is found
     * Note ncp->vars.ndefined is already made consistent by this point.
     */
    for (i=start_vid; i<ncp->vars.ndefined; i++)
        noFill[i-start_vid] = (char)(ncp->vars.value[i]->no_fill);
        TRACE_COMM(MPI_Bcast)(noFill, (ncp->vars.ndefined - start_vid),
                              MPI_BYTE, 0, ncp->comm);
    for (i=start_vid; i<ncp->vars.ndefined; i++) {
        /* overwrite local's mode */
        ncp->vars.value[i]->no_fill = noFill[i-start_vid];
        if (!noFill[i-start_vid]) nVarsFill++;
    }
#else
    for (i=start_vid; i<ncp->vars.ndefined; i++) {
        noFill[i-start_vid] = ncp->vars.value[i]->no_fill;
        if (!noFill[i-start_vid]) nVarsFill++;
    }
#endif
    if (nVarsFill == 0) { /* no variables in fill mode */
        NCI_Free(noFill);
        return NC_NOERR;
    }

    /* find the number of write segments (upper bound) */
    nsegs = (size_t)(ncp->vars.ndefined + ncp->vars.num_rec_vars * nrecs);
    count  = (MPI_Offset*) NCI_Malloc(nsegs * SIZEOF_MPI_OFFSET);
    offset = (MPI_Aint*)   NCI_Malloc(nsegs * SIZEOF_MPI_AINT);

    /* calculate each segment's offset and count */
    buf_len = 0; /* total write amount, used to allocate buffer */
    j = 0;
    for (i=start_vid; i<ncp->vars.ndefined; i++) {
        if (noFill[i-start_vid]) continue;

        varp = ncp->vars.value[i];

        if (IS_RECVAR(varp)) continue; /* first, fixed-size variables only */

        if (varp->ndims == 0) var_len = 1; /* scalar */
        else                  var_len = varp->dsizes[0];

        /* divide evenly total number of variable's elements among processes */
        count[j] = var_len / nprocs;
        start = count[j] * rank;
        if (rank < var_len % nprocs) {
            start += rank;
            count[j]++;
        }
        else
            start += var_len % nprocs;

        /* calculate the starting file offset */
        start *= varp->xsz;
        start += varp->begin;
        offset[j] = (MPI_Aint)start;
        if (start != offset[j]) {
            DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
            if (status == NC_NOERR) status = err;
            noFill[i-start_vid] = 1; /* skip this variable */
            continue;
        }
        /* add up the buffer size */
        buf_len += count[j] * varp->xsz;

        j++; /* increase j even when count[j] is zero */
    }

    /* loop thru all record variables to find the aggregated write amount */
    for (recno=0; recno<nrecs; recno++) {
        for (i=start_vid; i<ncp->vars.ndefined; i++) {
            if (noFill[i-start_vid]) continue;

            varp = ncp->vars.value[i];
            if (!IS_RECVAR(varp)) continue; /* record variables only */

            if (varp->ndims <= 1) var_len = 1;
            else                  var_len = varp->dsizes[1];

            /* divide total number of variable's elements among all processes */
            count[j] = var_len / nprocs;
            start = count[j] * rank;
            if (rank < var_len % nprocs) {
                start += rank;
                count[j]++;
            }
            else
                start += var_len % nprocs;

            /* calculate the starting file offset */
            start *= varp->xsz;
            start += varp->begin + ncp->recsize * recno;
            offset[j] = (MPI_Aint)start;
            if (start != offset[j]) {
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
                if (status == NC_NOERR) status = err;
                noFill[i-start_vid] = 1; /* skip this variable */
                continue;
            }
            /* add up the buffer size */
            buf_len += count[j] * varp->xsz;

            j++; /* increase j even when count[j] is zero */
        }
    }
    /* j is now the number of valid write segments */

    if (status == NC_NOERR && j == 0) {
        NCI_Free(noFill);
        NCI_Free(count);
        NCI_Free(offset);
        return NC_NOERR;
    }

    /* allocate one contiguous buffer space for all writes */
    blocklengths = (int*) NCI_Malloc((size_t)j * SIZEOF_INT);
    buf = NCI_Malloc((size_t)buf_len);
    buf_ptr = (char*)buf;

    /* fill write buffers for fixed-size variables first */
    j = k = 0;
    for (i=start_vid; i<ncp->vars.ndefined; i++) {
        if (noFill[i-start_vid]) continue;

        varp = ncp->vars.value[i];
        if (IS_RECVAR(varp)) continue;

        if (k < j) {  /* coalesce count[] and offset[] */
            count[k]  = count[j];
            offset[k] = offset[j];
        }
        j++;
        err = fill_var_buf(varp, count[k], buf_ptr);
        if (err != NC_NOERR) {
            if (status == NC_NOERR) status = err;
            continue; /* skip this request */
        }

        count[k] *= varp->xsz;
        if (count[k] != (int)count[k]) {
            DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
            if (status == NC_NOERR) status = err;
            continue; /* skip this request */
        }
        buf_ptr += count[k];
        blocklengths[k] = (int)count[k];
        k++;
    }
    /* k is the number of valid write requests thus far */

    /* loop thru all record variables to fill write buffers */
    for (recno=0; recno<nrecs; recno++) {
        for (i=start_vid; i<ncp->vars.ndefined; i++) {
            if (noFill[i-start_vid]) continue;

            varp = ncp->vars.value[i];
            if (!IS_RECVAR(varp)) continue;

            if (k < j) {  /* coalesce count[] and offset[] */
                count[k]  = count[j];
                offset[k] = offset[j];
            }
            j++;
            err = fill_var_buf(varp, count[k], buf_ptr);
            if (err != NC_NOERR) {
                if (status == NC_NOERR) status = err;
                continue; /* skip this request */
            }

            count[k] *= varp->xsz;
            if (count[k] != (int)count[k]) {
                DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
                if (status == NC_NOERR) status = err;
                continue; /* skip this request */
            }
            buf_ptr += count[k];
            blocklengths[k] = (int)count[k];
            k++;
        }
    }
    /* k is the number of valid write requests */
    NCI_Free(noFill);

    if (k == 0) {
        filetype = MPI_BYTE;
    }
    else {
        /* create fileview: a list of contiguous segment for each variable */
        mpireturn = MPI_Type_create_hindexed(k, blocklengths, offset, MPI_BYTE,
                                             &filetype);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else
            MPI_Type_commit(&filetype);
    }

    NCI_Free(blocklengths);
    NCI_Free(count);
    NCI_Free(offset);

    fh = ncp->collective_fh;

    TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, filetype, "native",
                                MPI_INFO_NULL);
    if (k > 0) MPI_Type_free(&filetype);

    if (buf_len != (int)buf_len) {
        if (status == NC_NOERR) status = NC_EINTOVERFLOW;
        buf_len = 0; /* skip this write */
    }

    /* write to variable collectively */
    TRACE_IO(MPI_File_write_at_all)(fh, 0, buf, (int)buf_len, MPI_BYTE, &mpistatus);

    NCI_Free(buf);

    TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                MPI_INFO_NULL);
    if (mpireturn != MPI_SUCCESS)
        if (status == NC_NOERR)
            status = ncmpii_error_mpi2nc(mpireturn, "MPI_File_set_view");

    return status;
}

/*----< ncmpio_fill_vars() >-------------------------------------------------*/
int
ncmpio_fill_vars(NC *ncp)
{
    int status=NC_NOERR;

#ifdef FILL_ONE_VAR_AT_A_TIME
    int err;

    /* fill variables according to their fill mode settings */
    if (NC_IsNew(ncp))
        /* file is just created */
        status = fillerup(ncp);
    else if (ncp->vars.ndefined > ncp->old->vars.ndefined) {
        /* old file, but new variables have been added */
        status = fill_added(ncp, ncp->old);

        err = fill_added_recs(ncp, ncp->old);
        if (status == NC_NOERR) status = err;
    }
#else
    /* fill variables according to their fill mode settings */
    if (NC_IsNew(ncp))
        /* file is just created and never exits define mode */
        status = fillerup_aggregate(ncp, NULL);
    else
        /* old file, but new variables have been added */
        status = fillerup_aggregate(ncp, ncp->old);
#endif
    return status;
}

/*----< ncmpio_fill_var_rec() >----------------------------------------------*/
/* fill an entire record of a record variable
 * this API is collective, must be called in data mode */
int
ncmpio_fill_var_rec(void      *ncdp,
                    int        varid,
                    MPI_Offset recno) /* record index, ignored if non-record var */
{
    int     indx, err=NC_NOERR;
    NC     *ncp=(NC*)ncdp;
    NC_var *varp=NULL;

    /* check whether variable ID is valid */
    /* sanity check for ncdp and varid has been done in dispatchers */
    varp = ncp->vars.value[varid];

    /* error if this is not a record variable */
    if (!IS_RECVAR(varp)) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTRECVAR)
        goto err_check;
    }

    /* check if _FillValue attribute is defined */
    indx = ncmpio_NC_findattr(&varp->attrs, _FillValue);

    /* error if the fill mode of this variable is not on */
    if (varp->no_fill && indx == -1) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTFILL)
        goto err_check;
    }

err_check:
    if (ncp->safe_mode) { /* consistency check */
        int root_varid, status, mpireturn;
        MPI_Offset root_recno;

        /* check if varid is consistent across all processes */
        root_varid = varid;
        TRACE_COMM(MPI_Bcast)(&root_varid, 1, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_varid != varid)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* check if recno is consistent across all processes */
        root_recno = recno;
        TRACE_COMM(MPI_Bcast)(&root_recno, 1, MPI_OFFSET, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && root_recno != recno)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");

        if (err == NC_NOERR) err = status;
    }
    if (err != NC_NOERR) return err;

    assert(varp != NULL);

    return fill_var_rec(ncp, varp, recno);
}

/*----< ncmpio_set_fill() >--------------------------------------------------*/
/* this API is collective, must be called in define mode, contrary to netCDF
 * where nc_set_fill() can also be called in data mode. The reason of PnetCDF
 * enforcing this requirement is because PnetCDF only fills fixed-size
 * variables at ncmpi_enddef() and record variables in ncmpi_fill_var_rec().
 */
int
ncmpio_set_fill(void *ncdp,
                int   fill_mode,
                int  *old_fill_mode)
{
    int i, mpireturn, oldmode;
    NC *ncp = (NC*)ncdp;

    if (ncp->safe_mode) {
        int err, status, root_fill_mode=fill_mode;

        TRACE_COMM(MPI_Bcast)(&root_fill_mode, 1, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return  ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (fill_mode != root_fill_mode)
            /* dataset's fill mode is inconsistent with root's */
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FILL_MODE)
        else
            err = NC_NOERR;

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (err == NC_NOERR) err = status;
        if (err != NC_NOERR) return err;
    }

    oldmode = fIsSet(ncp->flags, NC_MODE_FILL) ? NC_FILL : NC_NOFILL;

    if (fill_mode == NC_NOFILL)
        fClr(ncp->flags, NC_MODE_FILL);
    else if (fill_mode == NC_FILL)
        fSet(ncp->flags, NC_MODE_FILL);
    else
        DEBUG_RETURN_ERROR(NC_EINVAL) /* Invalid fill_mode */

    if (old_fill_mode != NULL) *old_fill_mode = oldmode;

    /* loop thru all variables defined so far to set/overwrite its fill mode */
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.value[i]->no_fill = (fill_mode == NC_NOFILL);

    /* once the file's fill mode is set, any new variables defined after this
     * call will check NC_dofill(ncp) and set their no_fill accordingly. See
     * ncmpi_def_var() */

    return NC_NOERR;
}

/*----< ncmpio_def_var_fill() >----------------------------------------------*/
/* this API is collective, and must be called in define mode */
int
ncmpio_def_var_fill(void       *ncdp,
                    int         varid,
                    int         no_fill,    /* 1: no fill, 0: fill */
                    const void *fill_value) /* when NULL, use default fill value */
{
    int err=NC_NOERR;
    NC *ncp=(NC*)ncdp;
    NC_var *varp=NULL;

    /* sanity check for ncdp and varid has been done in dispatchers */
    varp = ncp->vars.value[varid];

    if (ncp->safe_mode) {
        int root_ids[3], my_fill_null, minE, mpireturn;

        /* check if varid, no_fill, fill_value, are consistent */
        my_fill_null = (fill_value == NULL) ? 1 : 0;;
        root_ids[0] = varid;
        root_ids[1] = no_fill;
        root_ids[2] = my_fill_null;
        TRACE_COMM(MPI_Bcast)(&root_ids, 3, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR && (root_ids[0] != varid ||
            root_ids[1] != no_fill || root_ids[2] != my_fill_null))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* check fill_value, if not NULL, is consistent among processes */
        if (varp!= NULL && root_ids[2] == 0) {
            double root_fill_value; /* max nc_type space: 8 bytes */
            if (fill_value != NULL)
                memcpy(&root_fill_value, fill_value, (size_t)varp->xsz);
            TRACE_COMM(MPI_Bcast)(&root_fill_value, varp->xsz, MPI_BYTE, 0,
                                  ncp->comm);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
            if (err == NC_NOERR && fill_value != NULL &&
                memcmp(fill_value, &root_fill_value, (size_t)varp->xsz))
                /* variable's fill value is inconsistent with root's */
                DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_FILL_VALUE)
        }

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (err == NC_NOERR) err = minE;
    }
    if (err != NC_NOERR) return err;

    assert(varp != NULL);

    if (no_fill)
        varp->no_fill = 1;
    else
        varp->no_fill = 0;

    /* Are we setting a fill value? */
    if (fill_value != NULL && !varp->no_fill) {

        /* create/overwrite attribute _FillValue */
        err = ncmpio_put_att(ncdp, varid, _FillValue, varp->xtype,
                             1, fill_value, ncmpii_nc2mpitype(varp->xtype));
        if (err != NC_NOERR) return err;
    }

    return NC_NOERR;
}

/*----< ncmpio_inq_var_fill() >----------------------------------------------*/
int
ncmpio_inq_var_fill(NC_var *varp,
                    void   *fill_value) /* OUT: user-defined or
                                                default fill value */
{
    int i, err=NC_NOERR;
    const void *xp;
    NC_attrarray *ncap=&varp->attrs;

    assert(varp != NULL); /* NC_GLOBAL varid is illegal in this context */

    /* Check if _FillValue is defined for this variable */
    for (i=0; i<ncap->ndefined; i++) {
        if (strcmp(ncap->value[i]->name, _FillValue) == 0)
            break;
    }
    if (i == ncap->ndefined) { /* attribute _FillValue is not set */
        /* NetCDF 4.4.1 and prior does not use global attribute _FillValue if
         * it is not defined for the variable. Default fill values are used.
         * See fill_NC_var() in putget.m4.
         */
        err = ncmpio_inq_default_fill_value(varp->xtype, fill_value);
        return err;
    }

    /* retrieve user-defined attribute _FillValue */
    xp = ncap->value[i]->xvalue;

    /* value stored in xvalue is in external representation, may need byte-swap */
    switch(varp->xtype) {
        case NC_CHAR:   return ncmpix_getn_text               (&xp, 1,               (char*)fill_value);
        case NC_BYTE:   return ncmpix_getn_NC_BYTE_schar      (&xp, 1,        (signed char*)fill_value);
        case NC_UBYTE:  return ncmpix_getn_NC_UBYTE_uchar     (&xp, 1,      (unsigned char*)fill_value);
        case NC_SHORT:  return ncmpix_getn_NC_SHORT_short     (&xp, 1,              (short*)fill_value);
        case NC_USHORT: return ncmpix_getn_NC_USHORT_ushort   (&xp, 1,     (unsigned short*)fill_value);
        case NC_INT:    return ncmpix_getn_NC_INT_int         (&xp, 1,                (int*)fill_value);
        case NC_UINT:   return ncmpix_getn_NC_UINT_uint       (&xp, 1,       (unsigned int*)fill_value);
        case NC_FLOAT:  return ncmpix_getn_NC_FLOAT_float     (&xp, 1,              (float*)fill_value);
        case NC_DOUBLE: return ncmpix_getn_NC_DOUBLE_double   (&xp, 1,             (double*)fill_value);
        case NC_INT64:  return ncmpix_getn_NC_INT64_longlong  (&xp, 1,          (long long*)fill_value);
        case NC_UINT64: return ncmpix_getn_NC_UINT64_ulonglong(&xp, 1, (unsigned long long*)fill_value);
        default: return NC_EBADTYPE;
    }
}

