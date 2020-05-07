/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_enddef()  : dispatcher->enddef()
 * ncmpi__enddef() : dispatcher->_enddef()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>  /* strtol() */
#include <string.h>  /* memset() */
#include <assert.h>
#include <errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncx.h>
#include "ncmpio_NC.h"
#ifdef ENABLE_SUBFILING
#include "ncmpio_subfile.h"
#endif


/*----< move_file_block() >--------------------------------------------------*/
static int
move_file_block(NC         *ncp,
                MPI_Offset  to,
                MPI_Offset  from,
                MPI_Offset  nbytes)
{
    int rank, nprocs, bufcount, mpireturn, err, status=NC_NOERR, min_st;
    void *buf;
    int chunk_size=1048576; /* move 1 MB per process at a time */
    MPI_Status mpistatus;

    MPI_Comm_size(ncp->comm, &nprocs);
    MPI_Comm_rank(ncp->comm, &rank);

    /* if the file striping unit size is known (obtained from MPI-IO), then
     * we use that instead of 1 MB */
    if (ncp->striping_unit > 0) chunk_size = ncp->striping_unit;

    /* buf will be used as a temporal buffer to move data in chunks, i.e.
     * read a chunk and later write to the new location */
    buf = NCI_Malloc((size_t)chunk_size);
    if (buf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* make fileview entire file visible */
    TRACE_IO(MPI_File_set_view)(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE,
                                "native", MPI_INFO_NULL);

    /* move the variable starting from its tail toward its beginning */
    while (nbytes > 0) {
        int get_size=0;

        /* calculate how much to move at each time */
        bufcount = chunk_size;
        if (nbytes < (MPI_Offset)nprocs * chunk_size) {
            /* handle the last group of chunks */
            MPI_Offset rem_chunks = nbytes / chunk_size;
            if (rank > rem_chunks) /* these processes do not read/write */
                bufcount = 0;
            else if (rank == rem_chunks) /* this process reads/writes less */
                bufcount = (int)(nbytes % chunk_size);
            nbytes = 0;
        }
        else {
            nbytes -= chunk_size*nprocs;
        }

        /* explicitly initialize mpistatus object to 0. For zero-length read,
         * MPI_Get_count may report incorrect result for some MPICH version,
         * due to the uninitialized MPI_Status object passed to MPI-IO calls.
         * Thus we initialize it above to work around.
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));

        /* read the original data @ from+nbytes+rank*chunk_size */
        TRACE_IO(MPI_File_read_at_all)(ncp->collective_fh,
                                       from+nbytes+rank*chunk_size,
                                       buf, bufcount, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_read_at_all");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(status, NC_EREAD)
        }
        else {
            /* for zero-length read, MPI_Get_count may report incorrect result
             * for some MPICH version, due to the uninitialized MPI_Status
             * object passed to MPI-IO calls. Thus we initialize it above to
             * work around. See MPICH ticket:
             * https://trac.mpich.org/projects/mpich/ticket/2332
             *
             * Note we cannot set bufcount to get_size, as the actual size
             * read from a file may be less than bufcount. Because we are
             * moving whatever read to a new file offset, we must use the
             * amount actually read to call MPI_File_write_at_all below.
             */
            MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
            ncp->get_size += get_size;
        }

        /* MPI_Barrier(ncp->comm); */
        /* important, in case new region overlaps old region */
        TRACE_COMM(MPI_Allreduce)(&status, &min_st, 1, MPI_INT, MPI_MIN,
                                  ncp->comm);
        status = min_st;
        if (status != NC_NOERR) break;

        /* write to new location @ to+nbytes+rank*chunk_size
         *
         * Ideally, we should write the amount of get_size returned from a call
         * to MPI_Get_count in the below MPI write. This is in case some
         * variables are defined but never been written. The value returned by
         * MPI_Get_count is supposed to be the actual amount read by the MPI
         * read call. If partial data (or none) is available for read, then we
         * should just write that amount. Note this MPI write is collective,
         * and thus all processes must participate the call even if get_size
         * is 0. However, in some MPICH versions MPI_Get_count fails to report
         * the correct value due to an internal error that fails to initialize
         * the MPI_Status object. Therefore, the solution can be either to
         * explicitly initialize the status object to zeros, or to just use
         * bufcount for write. Note that the latter will write the variables
         * that have not been written before. Below uses the former option.
         */
#ifdef _USE_MPI_GET_COUNT
        /* explicitly initialize mpistatus object to 0. For zero-length read,
         * MPI_Get_count may report incorrect result for some MPICH version,
         * due to the uninitialized MPI_Status object passed to MPI-IO calls.
         * Thus we initialize it above to work around.
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));
#endif
        TRACE_IO(MPI_File_write_at_all)(ncp->collective_fh,
                                        to+nbytes+rank*chunk_size,
                                        buf, get_size /* bufcount */,
                                        MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_write_at_all");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
        }
        else {
#ifdef _USE_MPI_GET_COUNT
            int put_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
            ncp->put_size += put_size;
#else
            ncp->put_size += get_size; /* or bufcount */
#endif
        }
        TRACE_COMM(MPI_Allreduce)(&status, &min_st, 1, MPI_INT, MPI_MIN, ncp->comm);
        status = min_st;
        if (status != NC_NOERR) break;
    }
    NCI_Free(buf);
    return status;
}

/*----< move_fixed_vars() >--------------------------------------------------*/
/* move one fixed variable at a time, only when the new begin > old begin */
static int
move_fixed_vars(NC *ncp, NC *old)
{
    int i, err, status=NC_NOERR;

    /* move starting from the last fixed variable */
    for (i=old->vars.ndefined-1; i>=0; i--) {
        if (IS_RECVAR(old->vars.value[i])) continue;

        MPI_Offset from = old->vars.value[i]->begin;
        MPI_Offset to   = ncp->vars.value[i]->begin;
        if (to > from) {
            err = move_file_block(ncp, to, from, ncp->vars.value[i]->len);
            if (status == NC_NOERR) status = err;
        }
    }
    return status;
}

/*----< move_record_vars() >-------------------------------------------------*/
/* Move the record variables from lower offsets (old) to higher offsets. */
static int
move_record_vars(NC *ncp, NC *old) {
    int err;
    MPI_Offset recno;
    MPI_Offset nrecs = ncp->numrecs;
    MPI_Offset ncp_recsize = ncp->recsize;
    MPI_Offset old_recsize = old->recsize;
    MPI_Offset ncp_off = ncp->begin_rec;
    MPI_Offset old_off = old->begin_rec;

    assert(ncp_recsize >= old_recsize);

    if (ncp_recsize == old_recsize) {
        if (ncp_recsize == 0) /* no record variable defined yet */
            return NC_NOERR;

        /* No new record variable inserted, move the entire record variables
         * as a whole */
        err = move_file_block(ncp, ncp_off, old_off, ncp_recsize * nrecs);
        if (err != NC_NOERR) return err;
    } else {
        /* new record variables inserted, move one whole record at a time */
        for (recno = nrecs-1; recno >= 0; recno--) {
            err = move_file_block(ncp, ncp_off+recno*ncp_recsize,
                                       old_off+recno*old_recsize, old_recsize);
            if (err != NC_NOERR) return err;
        }
    }
    return NC_NOERR;
}

/*----< NC_begins() >--------------------------------------------------------*/
/*
 * This function is only called at enddef().
 * It computes each variable's 'begin' offset, and sets/updates the followings:
 *    ncp->xsz                   ---- header size
 *    ncp->vars.value[*]->begin  ---- each variable's 'begin' offset
 *    ncp->begin_var             ---- offset of first non-record variable
 *    ncp->begin_rec             ---- offset of first     record variable
 *    ncp->recsize               ---- sum of single records
 *    ncp->numrecs               ---- number of records (set only if new file)
 */
static int
NC_begins(NC *ncp)
{
    int i, j, rank, mpireturn;
    MPI_Offset end_var=0;
    NC_var *last = NULL;
    NC_var *first_var = NULL;       /* first "non-record" var */

    /* For CDF-1 and 2 formats, a variable's "begin" in the header is 4 bytes.
     * For CDF-5, it is 8 bytes.
     */

    /* get the true header size (not header extent) */
    MPI_Comm_rank(ncp->comm, &rank);
    ncp->xsz = ncmpio_hdr_len_NC(ncp);

    if (ncp->safe_mode) { /* this consistency check is redundant as metadata is
                             kept consistent at all time when safe mode is on */
        int err, status;
        MPI_Offset root_xsz = ncp->xsz;

        /* only root's header size matters */
        TRACE_COMM(MPI_Bcast)(&root_xsz, 1, MPI_OFFSET, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
            DEBUG_RETURN_ERROR(err)
        }

        err = NC_NOERR;
        if (root_xsz != ncp->xsz) err = NC_EMULTIDEFINE;

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
            DEBUG_RETURN_ERROR(err)
        }
        if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
    }

    /* This function is called in ncmpi_enddef(), which can happen either when
     * creating a new file or opening an existing file with metadata modified.
     * For the former case, ncp->begin_var == 0 here.
     * For the latter case, we set begin_var a new value only if the new header
     * grows out of its extent or the start of non-record variables is not
     * aligned as requested by ncp->h_align.
     * Note ncp->xsz is header size and ncp->begin_var is header extent.
     * Add the minimum header free space requested by user.
     */
    if (ncp->vars.ndefined > 0)
        ncp->begin_var = D_RNDUP(ncp->xsz + ncp->h_minfree, ncp->h_align);
    else /* no variable defined, ignore alignment and set header extent to
          * header size */
        ncp->begin_var = ncp->xsz;

    if (ncp->old != NULL) {
        /* If this define mode was entered from a redef(), we check whether
         * the new begin_var against the old begin_var. We do not shrink
         * the header extent.
         */
        if (ncp->begin_var < ncp->old->begin_var)
            ncp->begin_var = ncp->old->begin_var;
    }

    /* ncp->begin_var is the aligned starting file offset of the first
       variable, also the extent of file header */

    /* Now calculate the starting file offsets for all variables.
       loop thru vars, first pass is for the 'non-record' vars */
    end_var = ncp->begin_var;
    for (j=0, i=0; i<ncp->vars.ndefined; i++) {
        /* skip record variables on this pass */
        if (IS_RECVAR(ncp->vars.value[i])) continue;

        if (first_var == NULL) first_var = ncp->vars.value[i];

        /* for CDF-1 check if over the file size limit 32-bit integer */
        if (ncp->format == 1 && end_var > X_OFF_MAX)
            DEBUG_RETURN_ERROR(NC_EVARSIZE)

        /* this will pad out non-record variables with zero to the
         * requested alignment.  record variables are a bit trickier.
         * we don't do anything special with them */
        ncp->vars.value[i]->begin = D_RNDUP(end_var, ncp->fx_v_align);

        if (ncp->old != NULL) {
            /* move to the next fixed variable */
            for (; j<ncp->old->vars.ndefined; j++)
                if (!IS_RECVAR(ncp->old->vars.value[j]))
                    break;
            if (j < ncp->old->vars.ndefined) {
                if (ncp->vars.value[i]->begin < ncp->old->vars.value[j]->begin)
                    /* the first ncp->vars.ndefined non-record variables should
                       be the same. If the new begin is smaller, reuse the old
                       begin */
                    ncp->vars.value[i]->begin = ncp->old->vars.value[j]->begin;
                j++;
            }
        }
        /* end_var is the end offset of variable i */
        end_var = ncp->vars.value[i]->begin + ncp->vars.value[i]->len;
    }

    /* end_var now is pointing to the end of last non-record variable */

    /* only (re)calculate begin_rec if there is not sufficient
     * space at end of non-record variables or if start of record
     * variables is not aligned as requested by ncp->r_align.
     * If the existing begin_rec is already >= index, then leave the
     * begin_rec as is (in case some non-record variables are deleted)
     */
    if (ncp->begin_rec < end_var ||
        ncp->begin_rec != D_RNDUP(ncp->begin_rec, ncp->fx_v_align))
        ncp->begin_rec = D_RNDUP(end_var, ncp->fx_v_align);

    /* expand free space for fixed variable section */
    if (ncp->begin_rec < end_var + ncp->v_minfree)
        ncp->begin_rec = D_RNDUP(end_var + ncp->v_minfree, ncp->fx_v_align);

    /* align the starting offset for record variable section */
    if (ncp->r_align > 1)
        ncp->begin_rec = D_RNDUP(ncp->begin_rec, ncp->r_align);

    if (ncp->old != NULL) {
        /* check whether the new begin_rec is smaller */
        if (ncp->begin_rec < ncp->old->begin_rec)
            ncp->begin_rec = ncp->old->begin_rec;
    }

    if (first_var != NULL) ncp->begin_var = first_var->begin;
    else                   ncp->begin_var = ncp->begin_rec;

    end_var = ncp->begin_rec;
    /* end_var now is pointing to the beginning of record variables
     * note that this can be larger than the end of last non-record variable
     */

    ncp->recsize = 0;

    /* TODO: alignment for record variables (maybe using a new hint) */

    /* loop thru vars, second pass is for the 'record' vars,
     * re-calculate the starting offset for each record variable */
    for (j=0, i=0; i<ncp->vars.ndefined; i++) {
        if (!IS_RECVAR(ncp->vars.value[i]))
            /* skip non-record variables on this pass */
            continue;

        /* X_OFF_MAX is the max of 32-bit integer */
        if (ncp->format == 1 && end_var > X_OFF_MAX)
            DEBUG_RETURN_ERROR(NC_EVARSIZE)

        /* A few attempts at aligning record variables have failed
         * (either with range error or 'value read not that expected',
         * or with an error in ncmpi_redef()).  Not sufficient to align
         * 'begin', but haven't figured out what else to adjust */
        ncp->vars.value[i]->begin = end_var;

        if (ncp->old != NULL) {
            /* move to the next record variable */
            for (; j<ncp->old->vars.ndefined; j++)
                if (IS_RECVAR(ncp->old->vars.value[j]))
                    break;
            if (j < ncp->old->vars.ndefined) {
                if (ncp->vars.value[i]->begin < ncp->old->vars.value[j]->begin)
                    /* if the new begin is smaller, use the old begin */
                    ncp->vars.value[i]->begin = ncp->old->vars.value[j]->begin;
                j++;
            }
        }
        end_var += ncp->vars.value[i]->len;
        /* end_var is the end offset of record variable i */

        /* check if record size must fit in 32-bits (for CDF-1) */
#if SIZEOF_OFF_T == SIZEOF_SIZE_T && SIZEOF_SIZE_T == 4
        if (ncp->recsize > X_UINT_MAX - ncp->vars.value[i]->len)
            DEBUG_RETURN_ERROR(NC_EVARSIZE)
#endif
        ncp->recsize += ncp->vars.value[i]->len;
        last = ncp->vars.value[i];
    }

    /*
     * for special case (Check CDF-1 and CDF-2 file format specifications.)
     * "A special case: Where there is exactly one record variable, we drop the
     * requirement that each record be four-byte aligned, so in this case there
     * is no record padding."
     */
    if (last != NULL) {
        if (ncp->recsize == last->len) {
            /* exactly one record variable, pack value */
            ncp->recsize = *last->dsizes * last->xsz;
        }
#if 0
        else if (last->len == UINT32_MAX) { /* huge last record variable */
            ncp->recsize += *last->dsizes * last->xsz;
        }
#endif
    }

/* below is only needed if alignment is performed on record variables */
#if 0
    /*
     * for special case of exactly one record variable, pack value
     */
    /* if there is exactly one record variable, then there is no need to
     * pad for alignment -- there's nothing after it */
    if (last != NULL && ncp->recsize == last->len)
        ncp->recsize = *last->dsizes * last->xsz;
#endif

    if (NC_IsNew(ncp)) ncp->numrecs = 0;

    return NC_NOERR;
}

/*----< write_NC() >---------------------------------------------------------*/
/*
 * This function is collective and only called by enddef().
 * Write out the header
 * 1. Call ncmpio_hdr_put_NC() to copy the header object, ncp, to a buffer.
 * 2. Process rank 0 writes the header to file.
 */
static int
write_NC(NC *ncp)
{
    void *buf=NULL;
    int status=NC_NOERR, mpireturn, err, rank, header_wlen;
    size_t bufLen;
    MPI_Status mpistatus;

    assert(!NC_readonly(ncp));

    if (ncp->begin_var > X_INT_MAX) /* a fatal error */
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    MPI_Comm_rank(ncp->comm, &rank);

    /* only rank 0's header gets written to the file */
    if (rank == 0) {

        /* In NC_begins(), root's ncp->xsz and ncp->begin_var, root's header
         * size and extent, have been broadcast (sync-ed) among processes.
         */
#ifdef ENABLE_NULL_BYTE_HEADER_PADDING
        /* NetCDF classic file formats require the file header null-byte
         * padded. Thus we must calloc a buffer of size equal to file header
         * extent.
         */
        header_wlen = (int) ncp->begin_var;
        bufLen = _RNDUP(header_wlen, X_ALIGN);
        buf = NCI_Calloc(bufLen, 1);
#else
        /* Do not write padding area (between ncp->xsz and ncp->begin_var) */
        header_wlen = (int) ncp->xsz;
        bufLen = _RNDUP(header_wlen, X_ALIGN);
        buf = NCI_Malloc(bufLen);
#endif

        /* copy the entire local header object to buf */
        status = ncmpio_hdr_put_NC(ncp, buf);
        if (status != NC_NOERR) /* a fatal error */
            goto fn_exit;

        /* For non-fatal error, we continue to write header to the file, as now
         * the header object in memory has been sync-ed across all processes.
         */

        /* rank 0's fileview already includes the file header */

#ifdef _USE_MPI_GET_COUNT
        /* explicitly initialize mpistatus object to 0. For zero-length read,
         * MPI_Get_count may report incorrect result for some MPICH version,
         * due to the uninitialized MPI_Status object passed to MPI-IO calls.
         * Thus we initialize it above to work around.
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));
#endif
        /* write the header extent, not just header size, because NetCDF file
         * format specification requires null byte padding for header.
         */
        TRACE_IO(MPI_File_write_at)(ncp->collective_fh, 0, buf,
                                    header_wlen, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_write_at");
            /* write has failed, which is more serious than inconsistency */
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
        }
        else {
#ifdef _USE_MPI_GET_COUNT
            int put_size;
            MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
            ncp->put_size += put_size;
#else
            ncp->put_size += header_wlen;
#endif
        }
    }

fn_exit:
    if (buf != NULL) NCI_Free(buf);

    if (ncp->safe_mode == 1) {
        /* broadcast root's status, because only root writes to the file */
        int root_status = status;
        TRACE_COMM(MPI_Bcast)(&root_status, 1, MPI_INT, 0, ncp->comm);
        /* root's write has failed, which is more serious than inconsistency */
        if (root_status == NC_EWRITE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
    }

    fClr(ncp->flags, NC_NDIRTY);

    return status;
}

/* Many subroutines called in ncmpio__enddef() are collective. We check the
 * error codes of all processes only in safe mode, so the program can stop
 * collectively, if any one process got an error. However, when safe mode is
 * off, we simply return the error and program may hang if some processes
 * do not get error and proceed to the next subroutine call.
 */
#define CHECK_ERROR(err) {                                              \
    if (ncp->safe_mode == 1) {                                          \
        int status;                                                     \
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,   \
                                  ncp->comm);                           \
        if (mpireturn != MPI_SUCCESS) {                                 \
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");      \
            DEBUG_RETURN_ERROR(err)                                     \
        }                                                               \
        if (status != NC_NOERR) return status;                          \
    }                                                                   \
    else if (err != NC_NOERR)                                           \
        return err;                                                     \
}

/*----< ncmpio_NC_check_vlen() >---------------------------------------------*/
/* Check whether variable size is less than or equal to vlen_max,
 * without overflowing in arithmetic calculations.  If OK, return 1,
 * else, return 0.  For CDF1 format or for CDF2 format on non-LFS
 * platforms, vlen_max should be 2^31 - 4, but for CDF2 format on
 * systems with LFS it should be 2^32 - 4.
 */
int
ncmpio_NC_check_vlen(NC_var     *varp,
                     MPI_Offset  vlen_max)
{
    int i;
    MPI_Offset prod=varp->xsz;     /* product of xsz and dimensions so far */

    for (i = IS_RECVAR(varp) ? 1 : 0; i < varp->ndims; i++) {
        if (varp->shape[i] > vlen_max / prod) {
            return 0;           /* size in bytes > vlen_max */
        }
        prod *= varp->shape[i];
    }
    return 1;
}

/*----< ncmpio_NC_check_vlens() >--------------------------------------------*/
/* Given a valid ncp, check all variables for their sizes against the maximal
 * allowable sizes. Different CDF formation versions have different maximal
 * sizes. This function returns NC_EVARSIZE if any variable has a bad len
 * (product of non-rec dim sizes too large), else return NC_NOERR.
 */
int
ncmpio_NC_check_vlens(NC *ncp)
{
    int last = 0;
    MPI_Offset i, vlen_max, rec_vars_count;
    MPI_Offset large_fix_vars_count, large_rec_vars_count;
    NC_var *varp;

    if (ncp->vars.ndefined == 0) /* no variable defined */
        return NC_NOERR;

    /* maximum permitted variable size (or size of one record's worth
       of a record variable) in bytes. It is different between format 1
       2 and 5. */

    if (ncp->format >= 5) /* CDF-5 format max */
        vlen_max = X_INT64_MAX - 3; /* "- 3" handles rounded-up size */
    else if (ncp->format == 2) /* CDF2 format */
        vlen_max = X_UINT_MAX  - 3; /* "- 3" handles rounded-up size */
    else
        vlen_max = X_INT_MAX   - 3; /* CDF1 format */

    /* Loop through vars, first pass is for non-record variables */
    large_fix_vars_count = 0;
    rec_vars_count = 0;
    for (i=0; i<ncp->vars.ndefined; i++) {
        varp = ncp->vars.value[i];
        if (IS_RECVAR(varp)) {
            rec_vars_count++;
            continue;
        }

        last = 0;
        if (ncmpio_NC_check_vlen(varp, vlen_max) == 0) {
            /* check this variable's shape product against vlen_max */

            if (ncp->format >= 5) /* variable too big for CDF-5 */
                DEBUG_RETURN_ERROR(NC_EVARSIZE)

            large_fix_vars_count++;
            last = 1;
        }
    }
    /* OK if last non-record variable size too large, since not used to
       compute an offset */
    if (large_fix_vars_count > 1)  /* only one "too-large" variable allowed */
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    /* The only "too-large" variable must be the last one defined */
    if (large_fix_vars_count == 1 && last == 0)
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    if (rec_vars_count == 0) return NC_NOERR;

    /* if there is a "too-large" fixed-size variable, no record variable is
     * allowed */
    if (large_fix_vars_count == 1)
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    /* Loop through vars, second pass is for record variables.   */
    large_rec_vars_count = 0;
    for (i=0; i<ncp->vars.ndefined; i++) {
        varp = ncp->vars.value[i];
        if (!IS_RECVAR(varp)) continue;

        last = 0;
        if (ncmpio_NC_check_vlen(varp, vlen_max) == 0) {
            /* check this variable's shape product against vlen_max */

            if (ncp->format >= 5) /* variable too big for CDF-5 */
                DEBUG_RETURN_ERROR(NC_EVARSIZE)

            large_rec_vars_count++;
            last = 1;
        }
    }

    /* For CDF-2, no record variable can require more than 2^32 - 4 bytes of
     * storage for each record's worth of data, unless it is the last record
     * variable. See
     * http://www.unidata.ucar.edu/software/netcdf/docs/file_structure_and_performance.html#offset_format_limitations
     */
    if (large_rec_vars_count > 1)  /* only one "too-large" variable allowed */
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    /* and it has to be the last one */
    if (large_rec_vars_count == 1 && last == 0)
        DEBUG_RETURN_ERROR(NC_EVARSIZE)

    return NC_NOERR;
}

#ifdef VAR_BEGIN_IN_ARBITRARY_ORDER
typedef struct {
    MPI_Offset off;      /* starting file offset of a variable */
    MPI_Offset len;      /* length in bytes of a variable */
    int        ID;       /* variable index ID */
} off_len;

/*----< off_compare() >------------------------------------------------------*/
/* used for sorting the offsets of the off_len array */
static int
off_compare(const void *a, const void *b)
{
    if (((off_len*)a)->off > ((off_len*)b)->off) return  1;
    if (((off_len*)a)->off < ((off_len*)b)->off) return -1;
    return 0;
}
#endif

/*----< ncmpio_NC_check_voffs() >--------------------------------------------*/
/*
 * Given a valid ncp, check whether the file starting offsets (begin) of all
 * variables follows the same increasing order as they were defined.
 *
 * In NetCDF User's Guide, Chapter "File Structure and Performance", Section
 * "Parts of a NetCDF Class File", the following statement implies such
 * checking. "The order in which the variable data appears in each data section
 * is the same as the order in which the variables were defined, in increasing
 * numerical order by netCDF variable ID." URLs are given below.
 * https://www.unidata.ucar.edu/software/netcdf/documentation/historic/netcdf/Classic-File-Parts.html
 * https://www.unidata.ucar.edu/software/netcdf/docs/file_structure_and_performance.html#classic_file_parts
 *
 * However, the CDF file format specification does not require such order.
 * NetCDF version 4.6.0 and priors do not enforce this checking, but all
 * assume this requirement. See subroutine NC_computeshapes() in libsrc/v1hpg.c.
 * Similarly for PnetCDF, this check was not enforced until 1.9.0. Therefore,
 * it is important to keep this check to avoid potential problems.
 *
 * It appears that python scipy.netcdf does not follow this. An example can be
 * found in the NetCDF discussion thread:
 * https://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2018/msg00050.html
 * A scipy.netcdf program opens a NetCDF file with a few variables already
 * defined, enters define mode through redef call, and adds a new variable. The
 * scipy.netcdf implementation probably chooses to insert the new variable
 * entry in the front of "var_list" in file header. Technically speaking, this
 * does not violate the classic file format specification, but may result in
 * the file starting offsets ("begin" entry) of all variables defined in the
 * file header failed to appear in an increasing order. To obtain the file
 * header extent, one must scan the "begin" entry of all variables and find the
 * minimum as the extent.
 */
int
ncmpio_NC_check_voffs(NC *ncp)
{
    int i, num_fix_vars;

    if (ncp->vars.ndefined == 0) return NC_NOERR;

    num_fix_vars = ncp->vars.ndefined - ncp->vars.num_rec_vars;

#ifdef VAR_BEGIN_IN_ARBITRARY_ORDER
    int j;
    off_len *var_off_len;
    MPI_Offset var_end, max_var_end;

    if (num_fix_vars == 0) goto check_rec_var;

    /* check non-record variables first */
    var_off_len = (off_len*) NCI_Malloc(num_fix_vars * sizeof(off_len));
    for (i=0, j=0; i<ncp->vars.ndefined; i++) {
        NC_var *varp = ncp->vars.value[i];
        if (varp->begin < ncp->xsz) {
            if (ncp->safe_mode) {
                printf("Variable %s begin offset (%lld) is less than file header extent (%lld)\n",
                       varp->name, varp->begin, ncp->xsz);
            }
            NCI_Free(var_off_len);
            DEBUG_RETURN_ERROR(NC_ENOTNC)
        }
        if (IS_RECVAR(varp)) continue;
        var_off_len[j].off = varp->begin;
        var_off_len[j].len = varp->len;
        var_off_len[j].ID  = i;
        j++;
    }
    assert(j == num_fix_vars);

    for (i=1; i<num_fix_vars; i++) {
        if (var_off_len[i].off < var_off_len[i-1].off)
            break;
    }

    if (i < num_fix_vars)
        /* sort the off-len array into an increasing order */
        qsort(var_off_len, num_fix_vars, sizeof(off_len), off_compare);

    max_var_end = var_off_len[0].off + var_off_len[0].len;
    for (i=1; i<num_fix_vars; i++) {
        if (var_off_len[i].off < var_off_len[i-1].off + var_off_len[i-1].len) {
            if (ncp->safe_mode) {
                NC_var *var_cur = ncp->vars.value[var_off_len[i].ID];
                NC_var *var_prv = ncp->vars.value[var_off_len[i-1].ID];
                printf("Variable %s begin offset (%lld) overlaps variable %s (begin=%lld, length=%lld)\n",
                       var_cur->name, var_cur->begin, var_prv->name, var_prv->begin, var_prv->len);
            }
            NCI_Free(var_off_len);
            DEBUG_RETURN_ERROR(NC_ENOTNC)
        }
        var_end = var_off_len[i].off + var_off_len[i].len;
        max_var_end = MAX(max_var_end, var_end);
    }

    if (ncp->begin_rec < max_var_end) {
        if (ncp->safe_mode)
            printf("Record variable section begin (%lld) is less than fixed-size variable section end (%lld)\n",
                   ncp->begin_rec, max_var_end);
        NCI_Free(var_off_len);
        DEBUG_RETURN_ERROR(NC_ENOTNC)
    }
    NCI_Free(var_off_len);

check_rec_var:
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* check record variables */
    var_off_len = (off_len*) NCI_Malloc(ncp->vars.num_rec_vars * sizeof(off_len));
    for (i=0, j=0; i<ncp->vars.ndefined; i++) {
        NC_var *varp = ncp->vars.value[i];
        if (!IS_RECVAR(varp)) continue;
        var_off_len[j].off = varp->begin;
        var_off_len[j].len = varp->len;
        var_off_len[j].ID  = i;
        j++;
    }
    assert(j == ncp->vars.num_rec_vars);

    for (i=1; i<ncp->vars.num_rec_vars; i++) {
        if (var_off_len[i].off < var_off_len[i-1].off)
            break;
    }

    if (i < ncp->vars.num_rec_vars)
        /* sort the off-len array into an increasing order */
        qsort(var_off_len, ncp->vars.num_rec_vars, sizeof(off_len), off_compare);

    for (i=1; i<ncp->vars.num_rec_vars; i++) {
        if (var_off_len[i].off < var_off_len[i-1].off + var_off_len[i-1].len) {
            if (ncp->safe_mode) {
                NC_var *var_cur = ncp->vars.value[var_off_len[i].ID];
                NC_var *var_prv = ncp->vars.value[var_off_len[i-1].ID];
                printf("Variable %s begin offset (%lld) overlaps variable %s (begin=%lld, length=%lld)\n",
                       var_cur->name, var_cur->begin, var_prv->name, var_prv->begin, var_prv->len);
            }
            NCI_Free(var_off_len);
            DEBUG_RETURN_ERROR(NC_ENOTNC)
        }
    }
    NCI_Free(var_off_len);
#else
    /* Loop through vars, first pass is for non-record variables */
    if (num_fix_vars == 0) goto check_rec_var;

    int prev = 0;
    MPI_Offset prev_off = ncp->begin_var;

    for (i=0; i<ncp->vars.ndefined; i++) {
        NC_var *varp = ncp->vars.value[i];
        if (IS_RECVAR(varp)) continue;

        if (varp->begin < prev_off) {
            if (ncp->safe_mode) {
                if (i == 0)
                    printf("Variable \"%s\" begin offset (%lld) is less than header extent (%lld)\n",
                           varp->name, varp->begin, prev_off);
                else
                    printf("Variable \"%s\" begin offset (%lld) is less than previous variable \"%s\" end offset (%lld)\n",
                           varp->name, varp->begin, ncp->vars.value[prev]->name, prev_off);
            }
            DEBUG_RETURN_ERROR(NC_ENOTNC)
        }
        prev_off = varp->begin + varp->len;
        prev     = i;
    }

    if (ncp->begin_rec < prev_off) {
        if (ncp->safe_mode)
            printf("Record variable section begin offset (%lld) is less than fixed-size variable section end offset (%lld)\n",
                   ncp->begin_rec, prev_off);
        DEBUG_RETURN_ERROR(NC_ENOTNC)
    }

check_rec_var:
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* Loop through vars, second pass is for record variables */
    prev_off = ncp->begin_rec;
    prev     = 0;
    for (i=0; i<ncp->vars.ndefined; i++) {
        NC_var *varp = ncp->vars.value[i];
        if (!IS_RECVAR(varp)) continue;

        if (varp->begin < prev_off) {
            if (ncp->safe_mode) {
                printf("Variable \"%s\" begin offset (%lld) is less than previous variable end offset (%lld)\n",
                           varp->name, varp->begin, prev_off);
                if (i == 0)
                    printf("Variable \"%s\" begin offset (%lld) is less than record variable section begin offset (%lld)\n",
                           varp->name, varp->begin, prev_off);
                else
                    printf("Variable \"%s\" begin offset (%lld) is less than previous variable \"%s\" end offset (%lld)\n",
                           varp->name, varp->begin, ncp->vars.value[prev]->name, prev_off);
            }
            DEBUG_RETURN_ERROR(NC_ENOTNC)
        }
        prev_off = varp->begin + varp->len;
        prev     = i;
    }
#endif

    return NC_NOERR;
}

/*----< ncmpio__enddef() >---------------------------------------------------*/
/* This is a collective subroutine.
 * h_minfree  Sets the pad at the end of the "header" section.
 * v_align    Controls the alignment of the beginning of the data section for
 *            fixed size variables.
 * v_minfree  Sets the pad at the end of the data section for fixed size
 *            variables.
 * r_align    Controls the alignment of the beginning of the data section for
 *            variables which have an unlimited dimension (record variables).
 */
int
ncmpio__enddef(void       *ncdp,
               MPI_Offset  h_minfree,
               MPI_Offset  v_align,
               MPI_Offset  v_minfree,
               MPI_Offset  r_align)
{
    int i, flag, striping_unit, mpireturn, err=NC_NOERR, status=NC_NOERR;
    char value[MPI_MAX_INFO_VAL];
    MPI_Offset all_fix_var_size;
    NC *ncp = (NC*)ncdp;

    /* sanity check for NC_ENOTINDEFINE, NC_EINVAL, NC_EMULTIDEFINE_FNC_ARGS
     * has been done at dispatchers */
    ncp->h_minfree = h_minfree;
    ncp->v_minfree = v_minfree;

    /* calculate a good align size for PnetCDF level hints:
     * header_align_size and var_align_size based on the MPI-IO hint
     * striping_unit. This hint can be either supplied by the user or obtained
     * from MPI-IO (for example, ROMIO's Lustre driver makes a system call to
     * get the striping parameters of a file).
     */
    MPI_Info_get(ncp->mpiinfo, "striping_unit", MPI_MAX_INFO_VAL-1, value,
                 &flag);
    striping_unit = 0;
    if (flag) {
        errno = 0;
        striping_unit = (int)strtol(value,NULL,10);
        if (errno != 0) striping_unit = 0;
    }
    ncp->striping_unit = striping_unit;

    all_fix_var_size = 0;  /* sum of all defined fixed-size variables */
    for (i=0; i<ncp->vars.ndefined; i++) {
        if (IS_RECVAR(ncp->vars.value[i])) continue;
        all_fix_var_size += ncp->vars.value[i]->len;
    }

    /* ncp->h_align, ncp->v_align, ncp->r_align, and ncp->chunk have been
     * set during file create/open */

    if (ncp->h_align == 0) { /* user info does not set nc_header_align_size */
        if (v_align > 0)
            ncp->h_align = v_align;
        else if (all_fix_var_size > 0 && r_align > 0) /* no fixed-size variable */
            ncp->h_align = r_align;
        else if (striping_unit &&
            all_fix_var_size > FILE_ALIGNMENT_LB * striping_unit)
            /* if striping_unit is available and file size sufficiently large */
            ncp->h_align = striping_unit;
        else
            ncp->h_align = FILE_ALIGNMENT_DEFAULT;
    }
    /* else respect user hint */

    if (ncp->v_align == 0) { /* user info does not set nc_var_align_size */
        if (v_align > 0) /* else respect user hint */
            ncp->v_align = v_align;
#ifdef USE_AGGRESSIVE_ALIGN
        else if (striping_unit &&
                 all_fix_var_size > FILE_ALIGNMENT_LB * striping_unit)
            /* if striping_unit is available and file size sufficiently large */
            ncp->v_align = striping_unit;
        else
            ncp->v_align = FILE_ALIGNMENT_DEFAULT;
#endif
    }

    if (ncp->r_align == 0) { /* user info/env does not set nc_record_align_size */
        if (r_align > 0) /* else respect user hint */
            ncp->r_align = r_align;
#ifdef USE_AGGRESSIVE_ALIGN
        if (striping_unit)
            ncp->r_align = striping_unit;
        else
            ncp->r_align = FILE_ALIGNMENT_DEFAULT;
#endif
    }
    /* else ncp->r_align is already set by user/env, ignore the one passed by
     * the argument r_align of this subroutine.
     */

    /* all CDF formats require 4-bytes alignment */
    if (ncp->h_align == 0)    ncp->h_align = 4;
    else                      ncp->h_align = D_RNDUP(ncp->h_align, 4);
    if (ncp->v_align == 0)    ncp->v_align = 4;
    else                      ncp->v_align = D_RNDUP(ncp->v_align, 4);
    if (ncp->fx_v_align == 0) ncp->fx_v_align = 4;
    else                      ncp->fx_v_align = D_RNDUP(ncp->fx_v_align, 4);
    if (ncp->r_align == 0)    ncp->r_align = 4;
    else                      ncp->r_align = D_RNDUP(ncp->r_align, 4);

    /* reflect the hint changes to the MPI info object, so the user can inquire
     * what the true hint values are being used
     */
    sprintf(value, "%lld", ncp->h_align);
    MPI_Info_set(ncp->mpiinfo, "nc_header_align_size", value);
    sprintf(value, "%lld", ncp->fx_v_align);
    MPI_Info_set(ncp->mpiinfo, "nc_var_align_size", value);
    sprintf(value, "%lld", ncp->r_align);
    MPI_Info_set(ncp->mpiinfo, "nc_record_align_size", value);

#ifdef ENABLE_SUBFILING
    sprintf(value, "%d", ncp->num_subfiles);
    MPI_Info_set(ncp->mpiinfo, "nc_num_subfiles", value);
    if (ncp->num_subfiles > 1) {
        /* TODO: should return subfile-related msg when there's an error */
        err = ncmpio_subfile_partition(ncp);
        CHECK_ERROR(err)
    }
#else
    MPI_Info_set(ncp->mpiinfo, "pnetcdf_subfiling", "disable");
    MPI_Info_set(ncp->mpiinfo, "nc_num_subfiles", "0");
#endif

    /* check whether sizes of all variables are legal */
    err = ncmpio_NC_check_vlens(ncp);
    CHECK_ERROR(err)

    /* When ncp->old == NULL, this enddef is called the first time after file
     * create call. In this case, we compute each variable's 'begin', starting
     * file offset as well as the offsets of record variables.
     * When ncp->old != NULL, this enddef is called after a redef. In this
     * case, we re-used all variable offsets as many as possible.
     *
     * Note in NC_begins, root broadcasts ncp->xsz, the file header size, to
     * all processes.
     */
    err = NC_begins(ncp);
    CHECK_ERROR(err)

    /* update the total number of record variables */
    ncp->vars.num_rec_vars = 0;
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.num_rec_vars += IS_RECVAR(ncp->vars.value[i]);

    if (ncp->safe_mode) {
        /* check whether variable begins are in an increasing order.
         * This check is for debugging purpose. */
        err = ncmpio_NC_check_voffs(ncp);
        CHECK_ERROR(err)
    }

#ifdef ENABLE_SUBFILING
    if (ncp->num_subfiles > 1) {
        /* get ncp info for the subfile */
        err = NC_begins(ncp->ncp_sf);
        CHECK_ERROR(err)

        if (ncp->safe_mode) {
            /* check whether variable begins are in an increasing order.
             * This check is for debugging purpose. */
            err = ncmpio_NC_check_voffs(ncp->ncp_sf);
            CHECK_ERROR(err)
        }
    }
#endif

    if (ncp->old != NULL) {
        /* The current define mode was entered from ncmpi_redef, not from
         * ncmpi_create. We must check if header has been expanded.
         */

        assert(!NC_IsNew(ncp));
        assert(fIsSet(ncp->flags, NC_MODE_DEF));
        assert(ncp->begin_rec >= ncp->old->begin_rec);
        assert(ncp->begin_var >= ncp->old->begin_var);
        assert(ncp->vars.ndefined >= ncp->old->vars.ndefined);
        /* ncp->numrecs has already sync-ed in ncmpi_redef */

        if (ncp->vars.ndefined > 0) { /* no. record and non-record variables */
            if (ncp->begin_var > ncp->old->begin_var) {
                /* header size increases, shift the entire data part down */
                /* shift record variables first */
                err = move_record_vars(ncp, ncp->old);
                CHECK_ERROR(err)

                /* shift non-record variables */
                /* err = move_vars_r(ncp, ncp->old); */
                err = move_fixed_vars(ncp, ncp->old);
                CHECK_ERROR(err)
            }
            else if (ncp->begin_rec > ncp->old->begin_rec ||
                     ncp->recsize   > ncp->old->recsize) {
                /* number of non-record variables increases, or
                   number of records of record variables increases,
                   shift and move all record variables down */
                err = move_record_vars(ncp, ncp->old);
                CHECK_ERROR(err)
            }
        }
    } /* ... ncp->old != NULL */

    /* first sync header objects in memory across all processes, and then root
     * writes the header to file. Note safe_mode error check will be done in
     * write_NC() */
    status = write_NC(ncp);

    /* we should continue to exit define mode, even if header is inconsistent
     * among processes, so the program can proceed, say to close file properly.
     * However, if ErrIsHeaderDiff(status) is true, this error should
     * be considered fatal, as inconsistency is about the data structure,
     * rather then contents (such as attribute values) */

#ifdef ENABLE_SUBFILING
    /* write header to subfile */
    if (ncp->num_subfiles > 1) {
        err = write_NC(ncp->ncp_sf);
        if (status == NC_NOERR) status = err;
    }
#endif

    /* fill variables according to their fill mode settings */
    if (ncp->vars.ndefined > 0) {
        err = ncmpio_fill_vars(ncp);
        if (status == NC_NOERR) status = err;
    }

    if (ncp->old != NULL) {
        ncmpio_free_NC(ncp->old);
        ncp->old = NULL;
    }
    fClr(ncp->flags, NC_MODE_CREATE | NC_MODE_DEF);

#ifdef ENABLE_SUBFILING
    if (ncp->num_subfiles > 1)
        fClr(ncp->ncp_sf->flags, NC_MODE_CREATE | NC_MODE_DEF);
#endif

    /* If the user sets NC_SHARE, we enforce a stronger data consistency */
    if (NC_doFsync(ncp))
        ncmpio_file_sync(ncp); /* calling MPI_File_sync() */

    return status;
}

/*----< ncmpio_enddef() >----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpio_enddef(void *ncdp)
{
    return ncmpio__enddef(ncdp, 0, 0, 0, 0);
}

