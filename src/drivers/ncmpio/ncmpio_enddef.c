/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

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

/* Divide the amount of data to be moved into chunks of size MOVE_UNIT each,
 * and assign chunks to all processes. If the number of chunks is larger than
 * the number of processes, carry out the data movement in multiple rounds.
 */
#define MOVE_UNIT 16777216

#ifdef USE_POSIX_IO_TO_MOVE
/*----< move_file_block() >-------------------------------------------------*/
/* Call POSIX I/O subroutines to move data */
#include <fcntl.h>      /* open() */
#include <sys/types.h>  /* open() */
#include <sys/stat.h>   /* open() */
#include <unistd.h>     /* pread(), pwrite(), close() */

static int
move_file_block(NC         *ncp,
                MPI_Offset  to,     /* destination starting file offset */
                MPI_Offset  from,   /* source      starting file offset */
                MPI_Offset  nbytes) /* amount to be moved */
{
    int fd, rank, nprocs, status=NC_NOERR, do_open;
    void *buf;
    size_t num_moves, mv_amnt, p_units;
    off_t off_last, off_from, off_to;
    char *path = ncmpii_remove_file_system_type_prefix(ncp->path);

    /* check if this is a valid move request */
    if (to == from || nbytes == 0) return NC_NOERR;

    rank = ncp->rank;
    nprocs = ncp->nprocs;

    /* buf will be used as a temporal buffer to move data in chunks, i.e.
     * read a chunk and later write to the new location
     */
    buf = NCI_Malloc(MOVE_UNIT);
    if (buf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    p_units = MOVE_UNIT * nprocs;
    num_moves = nbytes / p_units;
    if (nbytes % p_units) num_moves++;
    off_last = (num_moves - 1) * p_units + rank * MOVE_UNIT;
    off_from = from + off_last;
    off_to   = to   + off_last;
    mv_amnt  = nbytes % p_units;
    if (mv_amnt == 0 && nbytes > 0) mv_amnt = p_units;

    /* determine the subset of processes that have data to move */
    do_open = 0;
    if (nbytes >= p_units)
        do_open = 1;
    else {
        MPI_Offset n_units = nbytes / MOVE_UNIT;
        if (nbytes % MOVE_UNIT) n_units++;
        if (rank < n_units) do_open = 1;
    }

    if (do_open && (fd = open(path, O_RDWR)) == -1) {
        fprintf(stderr,"Error at %s line %d: open file %s (%s)\n",
                __func__,__LINE__,path,strerror(errno));
        DEBUG_RETURN_ERROR(NC_EFILE)
    }

    /* move the data section starting from its tail toward its beginning */
    while (nbytes > 0) {
        size_t chunk_size;
        ssize_t get_size, put_size;

        if (mv_amnt == p_units) {
            /* each rank moves amount of chunk_size */
            chunk_size = MOVE_UNIT;
        }
        else {
            /* when total move amount is less than p_units */
            size_t num_chunks = mv_amnt / MOVE_UNIT;
            if (mv_amnt % MOVE_UNIT) num_chunks++;
            if (rank < num_chunks) {
                chunk_size = MOVE_UNIT;
                if (rank == num_chunks - 1 && mv_amnt % MOVE_UNIT > 0)
                    chunk_size = mv_amnt % MOVE_UNIT;
                assert(chunk_size > 0);
            }
            else
                chunk_size = 0;
        }

        if (chunk_size > 0) {
            /* read from file at off_from for amount of chunk_size */
            get_size = pread(fd, buf, chunk_size, off_from);
            if (get_size < 0) {
                fprintf(stderr,
                "Error at %s line %d: pread file %s offset %lld size %zd (%s)\n",
                __func__,__LINE__,path,(long long)off_from,chunk_size,strerror(errno));
                DEBUG_RETURN_ERROR(NC_EREAD)
            }
            ncp->get_size += get_size;
        }
        else
            get_size = 0;

        /* to prevent from one rank's write run faster than other's read */
        if (ncp->nprocs > 1) MPI_Barrier(ncp->comm);

        if (get_size > 0) {
            /* Write to new location at off_to for amount of get_size. Assuming
             * the call to MPI_Get_count() above returns the actual amount of
             * data read from the file, i.e. get_size.
             */
            put_size = pwrite(fd, buf, get_size, off_to);
            if (put_size < 0) {
                fprintf(stderr,
                "Error at %s line %d: pwrite file %s offset %lld size %zd (%s)\n",
                __func__,__LINE__,path,(long long)off_to,get_size,strerror(errno));
                DEBUG_RETURN_ERROR(NC_EREAD)
            }
            ncp->put_size += put_size;
        }

        /* move on to the next round */
        mv_amnt   = p_units;
        off_from -= mv_amnt;
        off_to   -= mv_amnt;
        nbytes   -= mv_amnt;
    }

    if (do_open && close(fd) == -1)
        DEBUG_RETURN_ERROR(NC_EFILE)

    NCI_Free(buf);
    return status;
}
#else
/*----< move_file_block() >-------------------------------------------------*/
/* Call MPI I/O subroutines to move data */
static int
move_file_block(NC         *ncp,
                MPI_Offset  to,     /* destination starting file offset */
                MPI_Offset  from,   /* source      starting file offset */
                MPI_Offset  nbytes) /* amount to be moved */
{
    int rank, nprocs, status=NC_NOERR, do_coll;
    void *buf;
    size_t num_moves, mv_amnt, p_units;
    MPI_Offset off_last, off_from, off_to, rlen, wlen;
    MPI_Comm comm;

    if (to == from || nbytes == 0) return NC_NOERR;

    /* If intra-node aggregation is enabled, then only the aggregators perform
     * the movement.
     */
    if (ncp->num_aggrs_per_node > 0 && ncp->ina_comm == MPI_COMM_NULL)
        return NC_NOERR;

    comm = (ncp->ina_comm == MPI_COMM_NULL) ? ncp->comm : ncp->ina_comm;
    rank = (ncp->ina_comm == MPI_COMM_NULL) ? ncp->rank : ncp->ina_rank;
    nprocs = (ncp->ina_comm == MPI_COMM_NULL) ? ncp->nprocs : ncp->ina_nprocs;

    /* MPI-IO fileview has been reset in ncmpi_redef() to make the entire file
     * visible
     */

    /* Use MPI collective I/O subroutines to move data, only if nproc > 1 and
     * MPI-IO hint "romio_no_indep_rw" is set to true. Otherwise, use MPI
     * independent I/O subroutines, as the data partitioned among processes are
     * not interleaved and thus need no collective I/O.
     */
    do_coll = (nprocs > 1 && fIsSet(ncp->flags, NC_HCOLL));

    /* buf will be used as a temporal buffer to move data in chunks, i.e.
     * read a chunk and later write to the new location
     */
    buf = NCI_Malloc(MOVE_UNIT);
    if (buf == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    p_units = MOVE_UNIT * nprocs;
    num_moves = nbytes / p_units;
    if (nbytes % p_units) num_moves++;
    off_last = (num_moves - 1) * p_units + rank * MOVE_UNIT;
    off_from = from + off_last;
    off_to   = to   + off_last;
    mv_amnt  = nbytes % p_units;
    if (mv_amnt == 0 && nbytes > 0) mv_amnt = p_units;

    /* move the data section starting from its tail toward its beginning */
    while (nbytes > 0) {
        int chunk_size;

        if (mv_amnt == p_units) {
            /* each rank moves amount of chunk_size */
            chunk_size = MOVE_UNIT;
        }
        else {
            /* when total move amount is less than p_units */
            size_t num_chunks = mv_amnt / MOVE_UNIT;
            if (mv_amnt % MOVE_UNIT) num_chunks++;
            if (rank < num_chunks) {
                chunk_size = MOVE_UNIT;
                if (rank == num_chunks - 1 && mv_amnt % MOVE_UNIT > 0)
                    chunk_size = mv_amnt % MOVE_UNIT;
                assert(chunk_size > 0);
            }
            else
                chunk_size = 0;
        }

        PNCIO_View buf_view;
        buf_view.type = MPI_BYTE;
        buf_view.size = chunk_size;
        buf_view.count = 1;
        buf_view.is_contig = 1;

        /* read from file at off_from for amount of chunk_size */
        rlen = 0;
        if (do_coll)
            rlen = ncmpio_file_read_at_all(ncp, off_from, buf, buf_view);
        else if (chunk_size > 0)
            rlen = ncmpio_file_read_at(ncp, off_from, buf, buf_view);
        if (status == NC_NOERR && rlen < 0) status = (int)rlen;

        /* to prevent from one rank's write run faster than other's read */
        if (nprocs > 1) MPI_Barrier(comm);

        /* Write to new location at off_to for amount of rlen, the actual read
         * amount is rlen.
         */
        buf_view.size = rlen;
        wlen = 0;
        if (do_coll && rlen > 0)
            wlen = ncmpio_file_write_at_all(ncp, off_to, buf, buf_view);
        else if (rlen > 0)
            wlen = ncmpio_file_write_at(ncp, off_to, buf, buf_view);
        if (status == NC_NOERR && wlen < 0) status = (int)wlen;

        /* move on to the next round */
        mv_amnt   = p_units;
        off_from -= mv_amnt;
        off_to   -= mv_amnt;
        nbytes   -= mv_amnt;
    }

    NCI_Free(buf);
    return status;
}
#endif

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
NC_begins(NC         *ncp,
          MPI_Offset  v_align,
          MPI_Offset  r_align)
{
    int i, j, mpireturn;
    MPI_Offset end_var=0;
    NC_var *last = NULL;

    /* For CDF-1 and 2 formats, a variable's "begin" in the header is 4 bytes.
     * For CDF-5, it is 8 bytes.
     */

    /* get the true header size (not header extent) */
    ncp->xsz = ncmpio_hdr_len_NC(ncp);

    if (ncp->safe_mode && ncp->nprocs > 1) {
        /* this consistency check is redundant as metadata is kept consistent
         * at all time when safe mode is on
         */
        int err, status;
        MPI_Offset root_xsz = ncp->xsz;

        /* only root's header size matters */
        TRACE_COMM(MPI_Bcast)(&root_xsz, 1, MPI_OFFSET, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
            DEBUG_RETURN_ERROR(err)
        }

        err = NC_NOERR;
        if (root_xsz != ncp->xsz) DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,ncp->comm);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
            DEBUG_RETURN_ERROR(err)
        }
        if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
    }

    if (ncp->vars.ndefined == 0) {
        /* There is no variable defined, ignore alignment and set header extent
         * to header size.
         */
        ncp->begin_var = MAX(ncp->begin_var, ncp->xsz);
        ncp->begin_rec = ncp->begin_var;
        ncp->recsize   = 0;
        ncp->numrecs   = 0;
        ncp->v_align   = 0;
        ncp->r_align   = 0;
        return NC_NOERR;
    }

    /* calculate a good file extent alignment size based on user hints.
     * The precedence of hints:
     * + 1st priority: hints set in the environment variable PNETCDF_HINTS,
     *                 i.e. nc_var_align_size and nc_record_align_size
     *                 e.g. PNETCDF_HINTS="nc_var_align_size=1024"
     * + 2nd priority: hints set in the MPI info objects passed into calls to
     *                 ncmpi_create() and ncmpi_open()
     *                 e.g. MPI_Info_set("nc_var_align_size", "1024");
     * + 3rd priority: hints passed from arguments of ncmpi__enddef()
     *                 i.e. v_align and r_align
     *                 e.g. ncmpi__enddef(..., v_align=1024,...)
     *
     * Default values
     *       NC_DEFAULT_H_MINFREE for h_minfree
     *       NC_DEFAULT_V_ALIGN   for v_align
     *       NC_DEFAULT_V_MINFREE for v_minfree
     *       NC_DEFAULT_R_ALIGN   for r_align
     */

    /* determine header extent (alignment for the data section) */
    if (ncp->info_v_align == -1) {
        /* hint nc_var_align_size is not set */

        /* argument v_align is set by user */
        if (v_align > 0)
            ncp->v_align = D_RNDUP(v_align, 4);
        else
            ncp->v_align = NC_DEFAULT_V_ALIGN;
    }
    else
        ncp->v_align = D_RNDUP(ncp->info_v_align, 4);

    /* determine alignment for record variable section */
    if (ncp->info_r_align == -1) {
        /* hint nc_record_align_size is not set */

        /* argument r_align is set by user */
        if (r_align > 0)
            ncp->r_align = D_RNDUP(r_align, 4);
        else if (ncp->vars.ndefined > ncp->vars.num_rec_vars)
            ncp->r_align = NC_DEFAULT_R_ALIGN;
        else
            ncp->r_align = NC_DEFAULT_V_ALIGN;
    }
    else
        ncp->r_align = D_RNDUP(ncp->info_r_align, 4);

    /* This function is called in ncmpi_enddef(), which can happen either when
     * creating a new file and first time call to ncmpi_enddef(), or other
     * case, e.g. opening an existing file, calling ncmpi_redef(), and then
     * ncmpi_enddef(). For the former case, ncp->begin_var == 0. For the latter
     * case, ncp->begin_var must be > 0, as it is the orignial header extent.
     * We increase begin_var only if the new header size grows out of its
     * original extent, or the start of variable section is not aligned as
     * requested by ncp->v_align. Note ncp->xsz is header size and
     * ncp->begin_var is header extent. Growth of header extent must also
     * respect the minimum header free space requested by user.
     */

    /* warrant a free space at the end of header section */
    ncp->begin_var = MAX(ncp->begin_var, ncp->xsz + ncp->h_minfree);

    /* Previously begin_var may be calculated using a different h_minfree and
     * v_align. Thus it can be larger than this round's calculation.
     */
    if (ncp->old != NULL)
        ncp->begin_var = MAX(ncp->begin_var, ncp->old->begin_var);

    /* align header extent if there are fix-sized variables */
    if (ncp->vars.ndefined > ncp->vars.num_rec_vars)
        ncp->begin_var = D_RNDUP(ncp->begin_var, ncp->v_align);

    /* ncp->begin_var is the aligned starting file offset of the first
     * variable (also data section), which is the extent of file header
     * (header section). File extent may contain free space for header to grow.
     */

    /* Now calculate the starting file offsets for all variables.
     * loop thru vars, first pass is for the 'non-record' vars
     */
    end_var = ncp->begin_var;
    for (j=0, i=0; i<ncp->vars.ndefined; i++) {
        /* skip record variables on this pass */
        if (IS_RECVAR(ncp->vars.value[i])) continue;

        /* for CDF-1 check if over the file size limit 32-bit integer */
        if (ncp->format == 1 && end_var > NC_MAX_INT)
            DEBUG_RETURN_ERROR(NC_EVARSIZE)

        /* this will pad out non-record variables with the 4-byte alignment */
        ncp->vars.value[i]->begin = D_RNDUP(end_var, 4);

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
    /* end_var now is pointing to the end of last fix-sized variable */

    ncp->fix_end = D_RNDUP(end_var, 4);

    /* warrant a free space at the end of fix-sized variable section */
    if (ncp->vars.ndefined > ncp->vars.num_rec_vars)
        ncp->begin_rec = ncp->fix_end + ncp->v_minfree;
    else /* Ignore v_minfree when there is no fix-sized variable. */
        ncp->begin_rec = ncp->fix_end;

    /* Previously begin_rec may be calculated using a different v_minfree and
     * r_align. Thus it can be larger than this round's calculation.
     */
    if (ncp->old != NULL)
        ncp->begin_rec = MAX(ncp->begin_rec, ncp->old->begin_rec);

    /* align the starting offset of record variable section */
    ncp->begin_rec = D_RNDUP(ncp->begin_rec, ncp->r_align);

    /* When there is no fix_sized variable, set begin_var == begin_rec */
    if (ncp->vars.ndefined == ncp->vars.num_rec_vars)
        ncp->begin_var = ncp->begin_rec;

    if (ncp->old != NULL) {
        assert(ncp->begin_var >= ncp->old->begin_var);
        assert(ncp->begin_rec >= ncp->old->begin_rec);
    }

    /* Alignment r_align is only applicable to the record variable section,
     * not individual record variables.
     */

    /* Loop through record variables and calculate the starting offset of each
     * record variable.
     */
    end_var = ncp->begin_rec;
    ncp->recsize = 0;
    for (j=0, i=0; i<ncp->vars.ndefined; i++) {
        if (!IS_RECVAR(ncp->vars.value[i]))
            /* skip non-record variables on this pass */
            continue;

        /* NC_MAX_INT is the max of 32-bit integer */
        if (ncp->format == 1 && end_var > NC_MAX_INT)
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
        if (ncp->recsize > NC_MAX_UINT - ncp->vars.value[i]->len)
            DEBUG_RETURN_ERROR(NC_EVARSIZE)
#endif
        ncp->recsize += ncp->vars.value[i]->len;
        last = ncp->vars.value[i];
    }

    /* For special case (Check CDF-1 and CDF-2 file format specifications.)
     * "A special case: Where there is exactly one record variable, we drop the
     * requirement that each record be four-byte aligned, so in this case there
     * is no record padding."
     */
    if (last != NULL) { /* i.e. at least one record variable */
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

#if 0
    /* This code block is to align individual record variable, which is no
     * longer needed.
     *
     * for special case of exactly one record variable, pack value
     *
     * if there is exactly one record variable, then there is no need to
     * pad for alignment -- there's nothing after it.
     */
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
    int status=NC_NOERR, is_coll=0;
    MPI_Offset i, header_wlen, ntimes;
    PNCIO_View buf_view;

    assert(!NC_readonly(ncp));

    buf_view.is_contig = 1;

    /* Depending on whether NC_HCOLL is set, writing file header can be done
     * through either MPI collective or independent write call.
     * When * ncp->nprocs == 1, ncp->collective_fh == ncp->independent_fh
     * For those ranks participating the collective MPI write call, their
     * is_coll is set to 1, otherwise 0.
     */
    if (fIsSet(ncp->flags, NC_HCOLL)) {
        if (ncp->num_aggrs_per_node > 0)
            is_coll = (ncp->ina_nprocs > 1 && ncp->rank == ncp->my_aggr);
        else
            is_coll = (ncp->nprocs > 1);
    }

    /* In NC_begins(), root's ncp->xsz and ncp->begin_var, root's header
     * size and extent, have been broadcast (sync-ed) among processes.
     */

#ifdef ENABLE_NULL_BYTE_HEADER_PADDING
    /* NetCDF classic file formats require the file header null-byte padded.
     * PnetCDF's default is not to write the padding area (between ncp->xsz and
     * ncp->begin_var). When this padding feature is enabled, we write the
     * padding area only when writing the header the first time, i.e. creating
     * a new file, or the new header extent becomes larger than the old one.
     */
    if (ncp->old == NULL || ncp->begin_var > ncp->old->begin_var)
        header_wlen = ncp->begin_var;
    else
        header_wlen = ncp->xsz;
#else
    /* Do not write padding area (between ncp->xsz and ncp->begin_var) */
    header_wlen = ncp->xsz;
#endif

    header_wlen = PNETCDF_RNDUP(header_wlen, X_ALIGN);

    /* if header_wlen is > NC_MAX_INT, then write the header in chunks.
     * Note reading file header is already done in chunks. See
     * ncmpio_hdr_get_NC().
     */
    ntimes = header_wlen / NC_MAX_INT;
    if (header_wlen % NC_MAX_INT) ntimes++;

    /* only rank 0's header gets written to the file */
    if (ncp->rank == 0) {
        char *buf=NULL, *buf_ptr;
        MPI_Offset offset, remain;

#ifdef ENABLE_NULL_BYTE_HEADER_PADDING
        /* NetCDF classic file formats require the file header null-byte
         * padded. Thus we must calloc a buffer of size equal to file header
         * extent.
         */
        buf = (char*)NCI_Calloc(header_wlen, 1);
#else
        /* Do not write padding area (between ncp->xsz and ncp->begin_var) */
        buf = (char*)NCI_Malloc(header_wlen);
#endif

        /* copy the entire local header object to buf */
        status = ncmpio_hdr_put_NC(ncp, buf);
        if (status != NC_NOERR) /* a fatal error */
            goto fn_exit;

        /* For non-fatal error, we continue to write header to the file, as now
         * the header object in memory has been sync-ed across all processes.
         */

        /* rank 0's fileview already includes the file header */

        /* write the header in chunks */
        offset = 0;
        remain = header_wlen;
        buf_ptr = buf;
        buf_view.type = MPI_BYTE;
        buf_view.count = 1;
        for (i=0; i<ntimes; i++) {
            MPI_Offset wlen;
            buf_view.size = MIN(remain, NC_MAX_INT);
            if (is_coll)
                wlen = ncmpio_file_write_at_all(ncp, offset, buf_ptr, buf_view);
            else
                wlen = ncmpio_file_write_at(ncp, offset, buf_ptr, buf_view);
            if (status == NC_NOERR && wlen < 0) status = (int)wlen;

            offset  += buf_view.size;
            buf_ptr += buf_view.size;
            remain  -= buf_view.size;
        }
        NCI_Free(buf);
    }
    else if (is_coll) {
        /* other processes participate the collective call */
        buf_view.size = 0;
        for (i=0; i<ntimes; i++)
            ncmpio_file_write_at_all(ncp, 0, NULL, buf_view);
    }

fn_exit:
    if (ncp->safe_mode == 1 && ncp->nprocs > 1) {
        /* broadcast root's status, because only root writes to the file */
        int mpireturn, root_status = status;
        TRACE_COMM(MPI_Bcast)(&root_status, 1, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            status = ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        else if (root_status == NC_EWRITE)
            /* root's write has failed, more serious than inconsistency */
            DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
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
    if (ncp->safe_mode == 1 && ncp->nprocs > 1) {                       \
        int min_err;                                                    \
        TRACE_COMM(MPI_Allreduce)(&err, &min_err, 1, MPI_INT, MPI_MIN,  \
                                  ncp->comm);                           \
        if (mpireturn != MPI_SUCCESS) {                                 \
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");      \
            DEBUG_RETURN_ERROR(err)                                     \
        }                                                               \
        if (min_err != NC_NOERR) return min_err;                        \
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
        vlen_max = NC_MAX_INT64 - 3; /* "- 3" handles rounded-up size */
    else if (ncp->format == 2) /* CDF2 format */
        vlen_max = NC_MAX_UINT  - 3; /* "- 3" handles rounded-up size */
    else
        vlen_max = NC_MAX_INT   - 3; /* CDF1 format */

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
    int i, num_fix_vars, prev;
    MPI_Offset prev_off;

    if (ncp->vars.ndefined == 0) return NC_NOERR;

    num_fix_vars = ncp->vars.ndefined - ncp->vars.num_rec_vars;

#ifdef VAR_BEGIN_IN_ARBITRARY_ORDER
    int j;
    off_len *var_off_len;
    MPI_Offset var_end, max_var_end;

    if (num_fix_vars == 0) goto check_rec_var;

    /* check non-record variables first */
    var_off_len = (off_len*) NCI_Malloc(sizeof(off_len) * num_fix_vars);
    for (i=0, j=0; i<ncp->vars.ndefined; i++) {
        NC_var *varp = ncp->vars.value[i];
        if (varp->begin < ncp->xsz) {
            if (ncp->safe_mode) {
                printf("Variable %s begin offset ("OFFFMT") is less than file header extent ("OFFFMT")\n",
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
                printf("Variable %s begin offset ("OFFFMT") overlaps variable %s (begin="OFFFMT", length="OFFFMT")\n",
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
            printf("Record variable section begin ("OFFFMT") is less than fixed-size variable section end ("OFFFMT")\n",
                   ncp->begin_rec, max_var_end);
        NCI_Free(var_off_len);
        DEBUG_RETURN_ERROR(NC_ENOTNC)
    }
    NCI_Free(var_off_len);

check_rec_var:
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* check record variables */
    var_off_len = (off_len*) NCI_Malloc(sizeof(off_len) * ncp->vars.num_rec_vars);
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
                printf("Variable %s begin offset ("OFFFMT") overlaps variable %s (begin="OFFFMT", length="OFFFMT")\n",
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

    prev = 0;
    prev_off = ncp->begin_var;

    for (i=0; i<ncp->vars.ndefined; i++) {
        NC_var *varp = ncp->vars.value[i];
        if (IS_RECVAR(varp)) continue;

        if (varp->begin < prev_off) {
            if (ncp->safe_mode) {
                if (i == 0)
                    printf("Variable \"%s\" begin offset ("OFFFMT") is less than header extent ("OFFFMT")\n",
                           varp->name, varp->begin, prev_off);
                else
                    printf("Variable \"%s\" begin offset ("OFFFMT") is less than previous variable \"%s\" end offset ("OFFFMT")\n",
                           varp->name, varp->begin, ncp->vars.value[prev]->name, prev_off);
            }
            DEBUG_RETURN_ERROR(NC_ENOTNC)
        }
        prev_off = varp->begin + varp->len;
        prev     = i;
    }

    if (ncp->begin_rec < prev_off) {
        if (ncp->safe_mode)
            printf("Record variable section begin offset ("OFFFMT") is less than fixed-size variable section end offset ("OFFFMT")\n",
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
                printf("Variable \"%s\" begin offset ("OFFFMT") is less than previous variable end offset ("OFFFMT")\n",
                           varp->name, varp->begin, prev_off);
                if (i == 0)
                    printf("Variable \"%s\" begin offset ("OFFFMT") is less than record variable section begin offset ("OFFFMT")\n",
                           varp->name, varp->begin, prev_off);
                else
                    printf("Variable \"%s\" begin offset ("OFFFMT") is less than previous variable \"%s\" end offset ("OFFFMT")\n",
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

#if 0
/*----< read_hints() >-------------------------------------------------------*/
/* check only the following hints set in environment variable PNETCDF_HINTS or
 * MPI_Info object passed to ncmpi_create() and ncmpi_open().
 * nc_header_align_size, nc_var_align_size, and nc_record_align_size
 */
static void
read_hints(NC *ncp)
{
    char *warn_str="Warning: skip ill-formed hint set in PNETCDF_HINTS";
    char *env_str, *env_str_cpy, *hint, *next_hint, *key, *val, *deli;
    char *hint_saved=NULL;

    /* reset hints from environment variable PNETCDF_HINTS */
    ncp->env_v_align = -1;
    ncp->env_r_align = -1;

    /* get hints from the environment variable PNETCDF_HINTS, a string of
     * hints separated by ";" and each hint is in the form of hint=value. E.g.
     * "cb_nodes=16;romio_ds_write=true". If this environment variable is set,
     * it overrides the same hints that were set by MPI_Info_set() called in
     * the application program.
     */
    env_str = getenv("PNETCDF_HINTS");
    if (env_str == NULL) return;

    env_str_cpy = strdup(env_str);
    next_hint = env_str_cpy;

    do {
        hint = next_hint;
        deli = strchr(hint, ';');
        if (deli != NULL) {
            *deli = '\0'; /* add terminate char */
            next_hint = deli + 1;
        }
        else next_hint = "\0";
        if (hint_saved != NULL) free(hint_saved);

        /* skip all-blank hint */
        hint_saved = strdup(hint);
        if (strtok(hint, " \t") == NULL) continue;

        free(hint_saved);
        hint_saved = strdup(hint); /* save hint for error message */

        deli = strchr(hint, '=');
        if (deli == NULL) { /* ill-formed hint */
            printf("%s: '%s'\n", warn_str, hint_saved);
            continue;
        }
        *deli = '\0';

        /* hint key */
        key = strtok(hint, "= \t");
        if (key == NULL || NULL != strtok(NULL, "= \t")) {
            /* expect one token before = */
            printf("%s: '%s'\n", warn_str, hint_saved);
            continue;
        }

        /* hint value */
        val = strtok(deli+1, "= \t");
        if (NULL != strtok(NULL, "= \t")) { /* expect one token before = */
            printf("%s: '%s'\n", warn_str, hint_saved);
            continue;
        }

        if (!strcmp(key, "nc_header_align_size") && ncp->env_v_align == -1)
            ncp->env_v_align = atoll(val);
        else if (!strcmp(key, "nc_var_align_size"))
            ncp->env_v_align = atoll(val);
        else if (!strcmp(key, "nc_record_align_size"))
            ncp->env_r_align = atoll(val);

    } while (*next_hint != '\0');

    if (hint_saved != NULL) free(hint_saved);
    free(env_str_cpy);

    /* return no error as all hints are advisory */
}
#endif

/*----< ncmpio__enddef() >---------------------------------------------------*/
/* This is a collective subroutine.
 * h_minfree  Sets the pad at the end of the "header" section, i.e. at least
 *            this amount of free space includes at the end of header extent.
 * v_align    Controls the alignment of the beginning of the data section for
 *            fixed size variables. Its value is also the header extent.
 * v_minfree  Sets the pad at the end of the data section for fixed size
 *            variables, i.e. at least this amount of free space between the
 *            fixed-size variable section and record variable section.
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
    int i, mpireturn, err=NC_NOERR, status=NC_NOERR;
    char value[MPI_MAX_INFO_VAL];
    MPI_Offset saved_begin_var;
    NC *ncp = (NC*)ncdp;

    /* update the total number of record variables */
    ncp->vars.num_rec_vars = 0;
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.num_rec_vars += IS_RECVAR(ncp->vars.value[i]);

    /* h_minfree, v_align, v_minfree, r_align being -1 means this subroutine is
     * called from ncmpio_enddef().
     */

    /* Checking hints from environment variable PNETCDF_HINTS and MPI info has
     * been done at the dispatcher, combine_env_hints(), when create/open the
     * file. Calling read_hints(ncp); is no longer necessary here.
     */

    /* sanity check for NC_ENOTINDEFINE, NC_EINVAL, NC_EMULTIDEFINE_FNC_ARGS
     * has been done at dispatchers */
    ncp->h_minfree = (h_minfree < 0) ? NC_DEFAULT_H_MINFREE : h_minfree;
    ncp->v_minfree = (v_minfree < 0) ? NC_DEFAULT_V_MINFREE : v_minfree;

#ifdef ENABLE_SUBFILING
    if (ncp->num_subfiles > 1) {
        /* TODO: should return subfile-related msg when there's an error */
        err = ncmpio_subfile_partition(ncp);
        CHECK_ERROR(err)
    }
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
    saved_begin_var = ncp->begin_var;
    err = NC_begins(ncp, v_align, r_align);
    if (err != NC_NOERR) /* restore the original begin_var when failed */
        ncp->begin_var = saved_begin_var;
    CHECK_ERROR(err)

    if (ncp->safe_mode) {
        /* check whether variable begins are in an increasing order.
         * This check is for debugging purpose. */
        err = ncmpio_NC_check_voffs(ncp);
        CHECK_ERROR(err)
    }

#ifdef ENABLE_SUBFILING
    if (ncp->num_subfiles > 1) {
        /* get ncp info for the subfile */
        err = NC_begins(ncp->ncp_sf, v_align, r_align);
        CHECK_ERROR(err)

        if (ncp->safe_mode) {
            /* check whether variable begins are in an increasing order.
             * This check is for debugging purpose. */
            err = ncmpio_NC_check_voffs(ncp->ncp_sf);
            CHECK_ERROR(err)
        }
    }
#endif

    if (ncp->old != NULL && ncp->vars.ndefined > 0) {
        /* The current define mode was entered from ncmpi_redef, not from
         * ncmpi_create. We must individually check if the three sections of
         * header, fix-sized, and record variables have grown. This is only
         * required when there are variables defined.
         */
        int mov_done=0;
        MPI_Offset nbytes;

        assert(!NC_IsNew(ncp));
        assert(fIsSet(ncp->flags, NC_MODE_DEF));
        assert(ncp->begin_rec >= ncp->old->begin_rec);
        assert(ncp->begin_var >= ncp->old->begin_var);
        assert(ncp->vars.ndefined >= ncp->old->vars.ndefined);

        /* ncp->numrecs has already sync-ed in ncmpi_redef */

        /* Make sure all processes finish their I/O before any process starts
         * to read the data section.
         */
        if (ncp->nprocs > 1) MPI_Barrier(ncp->comm);

        mov_done = 0;

        /* move record variable section first */
        if (ncp->begin_rec > ncp->old->begin_rec ||
            ncp->vars.num_rec_vars > ncp->old->vars.num_rec_vars) {
            /* It is possible begin_rec remain the same after adding new record
             * variables, e.g. when both header extent and fix-sized variable
             * section did not grow.
             */
            if (ncp->vars.num_rec_vars == ncp->old->vars.num_rec_vars) {
                /* No new record variable is added. Move the entire record
                 * variable section as a single data chunk.
                 */
                nbytes = ncp->old->recsize * ncp->old->numrecs;

                err = move_file_block(ncp, ncp->begin_rec, ncp->old->begin_rec,
                                      nbytes);
                if (status == NC_NOERR) status = err;
            }
            else {
                /* Move one record variable at a time */
                err = move_record_vars(ncp, ncp->old);
                if (status == NC_NOERR) status = err;
            }
            mov_done = 1;
        }

        /* Move fix-sized variable section when starting offset grows and there
         * are fix-sized variables defined.
         */
        if (ncp->begin_var > ncp->old->begin_var &&
            ncp->old->vars.ndefined > ncp->old->vars.num_rec_vars) {

            nbytes = ncp->old->fix_end - ncp->old->begin_var;

            err = move_file_block(ncp, ncp->begin_var, ncp->old->begin_var,
                                  nbytes);
            if (status == NC_NOERR) status = err;
            mov_done = 1;
        }

        /* To prevent some ranks run faster than others and start to read after
         * exiting ncmpi_enddef(), while some processes are still moving the
         * data section
         */
        if (mov_done && ncp->nprocs > 1) MPI_Barrier(ncp->comm);

    } /* ... ncp->old != NULL */

    /* first sync header objects in memory across all processes, and then root
     * writes the header to file. Note safe_mode error check will be done in
     * write_NC().
     */
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

    if (ncp->mpiinfo != MPI_INFO_NULL) {
        /* reflect the hint changes to the MPI info object, so the user can
         * inquire what the true hint values are being used
         */
        sprintf(value, OFFFMT, ncp->v_align);
        MPI_Info_set(ncp->mpiinfo, "nc_var_align_size", value);
        sprintf(value, OFFFMT, ncp->r_align);
        MPI_Info_set(ncp->mpiinfo, "nc_record_align_size", value);

#ifdef ENABLE_SUBFILING
        sprintf(value, "%d", ncp->num_subfiles);
        MPI_Info_set(ncp->mpiinfo, "nc_num_subfiles", value);
#else
        MPI_Info_set(ncp->mpiinfo, "pnetcdf_subfiling", "disable");
        MPI_Info_set(ncp->mpiinfo, "nc_num_subfiles", "0");
#endif
    }

    return status;
}

/*----< ncmpio_enddef() >----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpio_enddef(void *ncdp)
{
    return ncmpio__enddef(ncdp, -1, -1, -1, -1);
}

