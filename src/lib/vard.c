/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include <mpi.h>

#include "nc.h"
#include "ncx.h"
#include "ncmpidtype.h"
#include "macro.h"
#ifdef ENABLE_SUBFILING
#include "subfile.h"
#endif

/* for write case, buf needs to swapped back if swapped previously */
#define FINAL_CLEAN_UP {                                                       \
    if (is_buf_swapped) /* byte-swap back to buf's original contents */        \
        ncmpii_in_swapn(buf, bnelems, ncmpix_len_nctype(varp->type));          \
                                                                               \
    if (cbuf != NULL && cbuf != buf) NCI_Free(cbuf);                           \
}


/*----< ncmpii_getput_vard() >------------------------------------------------*/
static int
ncmpii_getput_vard(NC               *ncp,
                   NC_var           *varp,
                   MPI_Datatype      filetype,  /* data type of the variable */
                   void             *buf,
                   MPI_Offset        bufcount,
                   MPI_Datatype      buftype,  /* data type of the bufer */
                   int               rw_flag,
                   int               io_method)
{
    void *cbuf=NULL;
    int i, isderived, el_size, mpireturn, err=NC_NOERR;
    int buftype_is_contig, filetype_is_contig, need_swap=0, is_buf_swapped=0;
    int filetype_size=0, buftype_size;
    MPI_Offset btnelems, bnelems=0, offset=0, orig_bufcount=bufcount;
    MPI_Status mpistatus;
    MPI_Datatype ptype, orig_buftype=buftype;
    MPI_File fh=MPI_FILE_NULL;

    if (filetype == MPI_DATATYPE_NULL) { /* this process does zero-length I/O */
        if (io_method == INDEP_IO) return NC_NOERR;
        bufcount = 0;
        goto err_check;
    }

    if (bufcount == 0 && buftype != MPI_DATATYPE_NULL) {
        /* if this process has nothing to read/write */
        if (io_method == INDEP_IO) return NC_NOERR;
        goto err_check;
    }

#ifdef ENABLE_SUBFILING
    /* call a separate routine if variable is stored in subfiles */
    if (varp->num_subfiles > 1) {
        printf("This feature for subfiling is yet to implement\n");
        return NC_ENOTSUPPORT;
    }
#endif

    /* PROBLEM: type size can only be a 4-byte integer */
    mpireturn = MPI_Type_size(filetype, &filetype_size);
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_handle_error(mpireturn, "MPI_Type_size");
        goto err_check;
    }

    if (filetype_size == 0) { /* zero-length request */
        if (io_method == INDEP_IO) return NC_NOERR;
        bufcount = 0;
        goto err_check;
    }

    MPI_Aint true_lb, true_extent;
    MPI_Type_get_true_extent(filetype, &true_lb, &true_extent);
    MPI_Aint lb, extent;
    MPI_Type_get_extent(filetype, &lb, &extent);

    if (!IS_RECVAR(varp)) {
        /* for fixed-size variable, extent should not be larger than the
         * variabe size */
        MPI_Offset var_size = varp->xsz;
        for (i=0; i<varp->ndims; i++)
            var_size *= varp->shape[i];

        if (extent > var_size) {
            err = NC_ETYPESIZE;
            goto err_check;
        }
    }

    cbuf = (void*) buf;

    /* find the element type of filetype */
    err = ncmpii_dtype_decode(filetype, &ptype, &el_size, &btnelems,
                              &isderived, &filetype_is_contig);
    if (err != NC_NOERR) goto err_check;

    /* element type of filetype must be the same as variable's type */
    if (ptype != ncmpii_nc2mpitype(varp->type)) {
        err = NC_ETYPE_MISMATCH;
        goto err_check;
    }

    if (buftype == MPI_DATATYPE_NULL) {
        /* In this case, bufcount is ignored and will be set to the size of
         * filetype. Note buf's data type must match the data type of variable
         * defined in the file - no data conversion will be done.
         */
        /* set buftype to the variable's data type */
        buftype = ncmpii_nc2mpitype(varp->type);
        MPI_Type_size(buftype, &buftype_size);
        bufcount = filetype_size / buftype_size;
        buftype_is_contig = 1;
        bnelems = bufcount;
    }
    else {
        /* find whether buftype is contiguous */
        err = ncmpii_dtype_decode(buftype, &ptype, &el_size, &btnelems,
                                  &isderived, &buftype_is_contig);
        if (err != NC_NOERR) goto err_check;

        bnelems = bufcount * btnelems;
        buftype_size = el_size * btnelems;
        if (filetype_size != bufcount * buftype_size)
            return NC_ETYPESIZE_MISMATCH;

        /* if buf is not contiguous, we need to pack it to one, cbuf */
        if (!buftype_is_contig && bnelems > 0) {
            int outsize = buftype_size * bufcount;
            cbuf = NCI_Malloc(outsize);

            if (rw_flag == WRITE_REQ) {
                /* pack buf into cbuf, a contiguous buffer */
                int position = 0;
                MPI_Pack(buf, bufcount, buftype, cbuf, outsize, &position,
                         MPI_COMM_SELF);
            }
            buftype = ptype;
            bufcount *= bnelems;
        }
    }

    /* Check if we need byte swap cbuf in-place or (into cbuf) */
    need_swap = ncmpii_need_swap(varp->type, ptype);
    if (need_swap) {
        if (rw_flag == WRITE_REQ) {
#ifdef DISABLE_IN_PLACE_SWAP
            if (cbuf == buf) {
                /* allocate cbuf and copy buf to cbuf, cbuf is to be freed */
                cbuf = NCI_Malloc(filetype_size);
                memcpy(cbuf, buf, filetype_size);
            }
#endif
            /* perform array in-place byte swap on cbuf */
            ncmpii_in_swapn(cbuf, bnelems, ncmpix_len_nctype(varp->type));
            is_buf_swapped = (cbuf == buf) ? 1 : 0;
            /* is_buf_swapped indicates if the contents of the original user
             * buffer, buf, have been changed, i.e. byte swapped. */
        }
    }
    /* no type conversion */

    /* set offset to the file offset of the variable */
    offset = varp->begin;

err_check:
    /* check API error from any proc before going into a collective call.
     * optimization: to avoid MPI_Allreduce to check parameters at
     * every call, we assume caller does the right thing most of the
     * time.  If caller passed in bad parameters, we'll still conduct a
     * zero-byte operation (everyone has to participate in the
     * collective I/O call) but return error */
    if (err != NC_NOERR || bufcount == 0 || filetype_size == 0) {
        if (io_method == INDEP_IO) {
            FINAL_CLEAN_UP  /* swap back put buffer and free temp buffers */
            return err;
        }
        /* else for COLL_IO, must participate successive collective calls */
    }

    if (io_method == COLL_IO)
        fh = ncp->nciop->collective_fh;
    else
        fh = ncp->nciop->independent_fh;

    /* MPI_File_set_view is a collective if (io_method == COLL_IO) */
    TRACE_IO(MPI_File_set_view)(fh, offset, MPI_BYTE, filetype,
                                  "native", MPI_INFO_NULL);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_handle_error(mpireturn, "MPI_File_set_view");

    if (rw_flag == WRITE_REQ) {
        if (io_method == COLL_IO) {
            TRACE_IO(MPI_File_write_all)(fh, cbuf, bufcount, buftype, &mpistatus);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_handle_error(mpireturn, "MPI_File_write_all");
        }
        else { /* io_method == INDEP_IO */
            TRACE_IO(MPI_File_write)(fh, cbuf, bufcount, buftype, &mpistatus);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_handle_error(mpireturn, "MPI_File_write");
        }
        int put_size;
        MPI_Get_count(&mpistatus, MPI_BYTE, &put_size);
        ncp->nciop->put_size += put_size;
    }
    else {  /* rw_flag == READ_REQ */
        if (io_method == COLL_IO) {
            TRACE_IO(MPI_File_read_all)(fh, cbuf, bufcount, buftype, &mpistatus);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_handle_error(mpireturn, "MPI_File_read_all");
        }
        else { /* io_method == INDEP_IO */
            TRACE_IO(MPI_File_read)(fh, cbuf, bufcount, buftype, &mpistatus);
            if (mpireturn != MPI_SUCCESS)
                return ncmpii_handle_error(mpireturn, "MPI_File_read");
        }
        int get_size;
        MPI_Get_count(&mpistatus, MPI_BYTE, &get_size);
        ncp->nciop->get_size += get_size;
    }

    /* reset the file view so the entire file is visible again */
    TRACE_IO(MPI_File_set_view)(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    if (rw_flag == READ_REQ) {
        if (need_swap)
            /* perform array in-place byte swap on cbuf */
            ncmpii_in_swapn(cbuf, bnelems, ncmpix_len_nctype(varp->type));

        if (!buftype_is_contig && bnelems > 0) {
            /* unpack cbuf, a contiguous buffer, to buf using buftype */
            int position = 0;
            int insize = bnelems * el_size;
            MPI_Unpack(cbuf, insize, &position, buf, orig_bufcount,
                       orig_buftype, MPI_COMM_SELF);
        }
    }
    else { /* WRITE_REQ */
        if (IS_RECVAR(varp)) {
            /* update header's number of records in memory */
            MPI_Offset new_numrecs;

            /* since filetype's LB is required to be == varp->begin for vard
             * API, we can simply use extent to calculate new_numrecs */
            new_numrecs = extent / ncp->recsize;
            if (extent % ncp->recsize) new_numrecs++;

            if (io_method == INDEP_IO) {
                /* For independent put, we delay the sync for numrecs until
                 * the next collective call, such as end_indep(), sync(),
                 * enddef(), or close(). This is because if we update numrecs
                 * to file now, race condition can happen. Note numrecs in
                 * memory may be inconsistent and obsolete till then.
                 */
                if (ncp->numrecs < new_numrecs) {
                    ncp->numrecs = new_numrecs;
                    set_NC_ndirty(ncp);
                }
            }
            else /* COLL_IO: sync numrecs in memory and file */
                err = ncmpii_sync_numrecs(ncp, new_numrecs);
        }

        if (NC_doFsync(ncp)) { /* NC_SHARE is set */
            TRACE_IO(MPI_File_sync)(fh);
            if (io_method == COLL_IO)
                TRACE_COMM(MPI_Barrier)(ncp->nciop->comm);
        }
    }

    FINAL_CLEAN_UP  /* swap back the put buffer and free temp buffers */

    return err;
}

/*----< ncmpi_get_vard() >---------------------------------------------------*/
int
ncmpi_get_vard(int           ncid,
               int           varid,
               MPI_Datatype  filetype,  /* data type of the variable in file */
               void         *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype)   /* data type of the bufer */
{
    int     status;
    NC     *ncp;
    NC_var *varp=NULL;

    SANITY_CHECK(ncid, ncp, varp, READ_REQ, INDEP_IO, status)

    return ncmpii_getput_vard(ncp, varp, filetype, buf, bufcount, buftype,
                              READ_REQ, INDEP_IO);
}
/*----< ncmpi_get_vard_all() >-----------------------------------------------*/
int
ncmpi_get_vard_all(int           ncid,
                   int           varid,
                   MPI_Datatype  filetype,  /* data type of the variable in file */
                   void         *buf,
                   MPI_Offset    bufcount,
                   MPI_Datatype  buftype)   /* data type of the bufer */
{
    int     status;
    NC     *ncp;
    NC_var *varp=NULL;

    SANITY_CHECK(ncid, ncp, varp, READ_REQ, COLL_IO, status)

    return ncmpii_getput_vard(ncp, varp, filetype, buf, bufcount, buftype,
                              READ_REQ, COLL_IO);
}

/*----< ncmpi_put_vard() >---------------------------------------------------*/
int
ncmpi_put_vard(int           ncid,
               int           varid,
               MPI_Datatype  filetype,  /* data type of the variable in file */
               const void   *buf,
               MPI_Offset    bufcount,
               MPI_Datatype  buftype)   /* data type of the bufer */
{
    int     status;
    NC     *ncp;
    NC_var *varp=NULL;

    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, INDEP_IO, status)

    return ncmpii_getput_vard(ncp, varp, filetype, (void*)buf, bufcount,
                              buftype, WRITE_REQ, INDEP_IO);
}
/*----< ncmpi_put_vard_all() >-----------------------------------------------*/
int
ncmpi_put_vard_all(int           ncid,
                   int           varid,
                   MPI_Datatype  filetype,  /* data type of the variable in file */
                   const void   *buf,
                   MPI_Offset    bufcount,
                   MPI_Datatype  buftype)   /* data type of the bufer */
{
    int     status;
    NC     *ncp;
    NC_var *varp=NULL;

    SANITY_CHECK(ncid, ncp, varp, WRITE_REQ, COLL_IO, status)

    return ncmpii_getput_vard(ncp, varp, filetype, (void*)buf, bufcount,
                              buftype, WRITE_REQ, COLL_IO);
}
