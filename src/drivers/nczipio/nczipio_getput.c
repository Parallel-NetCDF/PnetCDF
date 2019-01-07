/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>

static int get_block_idx(NC_var *varp, int* cord){
    int i, ret;
    
    ret = cord[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->stripesize[i - 1] + cord[i];
    }

    return ret;
}

int
nczipio_get_var(void             *ncdp,
              int               varid,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              void             *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               reqMode)
{
    int i, err;
    NC_zip *nczipp = (NC_zip*)ncdp;
    NC_var *varp;
    nc_type xtype;
    int *bstart, *bend, *bcord;
    int nb, bsize;
    int datavarid;
    int *tsize, *tssize, *tstart;
    int tpos;
    MPI_Datatype subarytype;
    char *rbuffer, *cbuffer;
    MPI_Offset cbsize;
    MPI_Offset **starts, **counts;

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    // Boundary of blocks involved
    *bstart = NCI_Malloc(sizeof(int) * varp->ndim);
    *bcord = NCI_Malloc(sizeof(int) * varp->ndim);
    *bend = NCI_Malloc(sizeof(int) * varp->ndim);
    for(i = 0; i < varp->ndim; i++){
        bstart[i] = start[i] / varp->stripesize;
        if (stride == NULL){
            bend[i] = (start[i] + count[i] - 1) / varp->stripesize;
        }
        else{
            bend[i] = (start[i] + (count[i] - 1) * stride[i]) / varp->stripesize + 1;
        }
    }
    
    // Number of blocks involved
    nb = 1;
    for(i = 0; i < varp->ndim; i++){
        nb *= bend[i] - bstart[i];
    }

    /* Use a varn call to read all compressed block involved 
     * Generate one request for each block
     */

    *bidx = NCI_Malloc(sizeof(int) * nb);
    **starts = NCI_Malloc(sizeof(MPI_Offset*) * nb);
    **counts = NCI_Malloc(sizeof(MPI_Offset*) * nb);
    // Iterate through all blocks involved
    i = 0;
    cbsize = 0;
    memcpy(bcord, bstart, sizeof(int) * varp->ndim);
    for(i = 0; i < nb; i++){
        j = get_block_idx(varp, bcord);   
        bidx[i] = j; // block idx
        cbsize += varp->lens[j];  // total buffer size of compressed data
        starts[i] = varp->offset + j;   // start of the block
        counts[i] = varp->lens + j; // count of the block

        // move on to next block
        bcord[varp->ndim - 1]++;
        for(j = varp->ndim - 1; j > 0; j--){
            if (bcord[j] >= bend[j]){
                bcord[j - 1]++;
                bcord[j] = bstart[j];
            }
        }
    }

    // Allocate buffers
    *cbuffer = NCI_Malloc(cbsize);  // Compressed data

    // Locate data var
    err = nczipp->driver->get_var(nczipp->ncp, varp->varid, NULL, NULL, NULL, NULL, &datavarid, 1, MPI_INT, reqMode); 
    if (err != NC_NOERR) return err;

    // read compressed data
    err = nczipp->driver->get_varn(nczipp->ncp, datavarid, nb, starts, counts, cbuffer, cbsize, MPI_BYTE, reqMode); 
    if (err != NC_NOERR) return err;

    // Decompression

    // Calculate block size
    // Original datatype
    err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_datatype", NC_INT, 1, &xtype, MPI_INT); 
    if (err != NC_NOERR) return err;

    // Calculate block size
    bsize = (int)NC_Type_size(xtype);
    for(i = 0; i < nblocks; i++){
        bsize *= varp->stripesize[i];
    }

    // Allocate buffers
    *rbuffer = NCI_Malloc(bsize * nb);  // Decompressed data

    // Decompress blocks
    cbsize = 0;
    for(i = 0; i < nb; i++){
        j = bidx[i];
        if (varp->lens[j] > 0){
            nczipp->zip->decompress(cbuffer + cbsize, varp->lens[j], rbuffer + bsize * i, NULL, carp->ndim, varp->dimsize, ncmpii_nc2mpitype(xtype));
        }
        else{
            memset(rbuffer + bsize * i, 0, bsize);
        }
        cbsize += varp->lens[j];  // move to next block location
    }

    // Copy data into user buffer

    // Create datatype of querying domain in the decompressed domain
    *tsize = NCI_Malloc(sizeof(int) * varp->ndim);
    *tssize = NCI_Malloc(sizeof(int) * varp->ndim);
    *tstart = NCI_Malloc(sizeof(int) * varp->ndim);
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = (bend[i] - bstart[i]) * varp->stripesize[i];
        tssize[i] = (int)count[i];
        tstart[i] = start[i] % varp->stripesize[i];
    }
    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, ncmpii_nc2mpitype(xtype), &subarytype);

    // Pack data into user buffer
    tpos = 0;
    MPI_Pack(rbuffer, bsize * nb, subarytype, buf, bsize * nb, &tpos, nczipp->comm);

    return NC_NOERR;
}

int
nczipio_get_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               void              *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{

    // Original datatype
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_datatype", NC_INT, 1, &xtype, MPI_INT); 
    if (err != NC_NOERR) return err;

    err = nczipp->driver->get_var(nczipp->ncp, varid, start, count, stride, imap,
                               buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipio_put_var(void             *ncdp,
              int               varid,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              const void       *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               reqMode)
{
    int err=NC_NOERR, status;
    void *cbuf=(void*)buf;
    NC_zip *nczipp = (NC_zip*)ncdp;

    if (imap != NULL || bufcount != -1) {
        /* pack buf to cbuf -------------------------------------------------*/
        /* If called from a true varm API or a flexible API, ncmpii_pack()
         * packs user buf into a contiguous cbuf (need to be freed later).
         * Otherwise, cbuf is simply set to buf. ncmpii_pack() also returns
         * etype (MPI primitive datatype in buftype), and nelems (number of
         * etypes in buftype * bufcount)
         */
        int ndims;
        MPI_Offset nelems;
        MPI_Datatype etype;

        err = nczipp->driver->inq_var(nczipp->ncp, varid, NULL, NULL, &ndims, NULL,
                                   NULL, NULL, NULL, NULL);
        if (err != NC_NOERR) goto err_check;

        err = ncmpii_pack(ndims, count, imap, (void*)buf, bufcount, buftype,
                          &nelems, &etype, &cbuf);
        if (err != NC_NOERR) goto err_check;

        imap     = NULL;
        bufcount = (nelems == 0) ? 0 : -1;  /* make it a high-level API */
        buftype  = etype;                   /* an MPI primitive type */
    }

err_check:
    if (err != NC_NOERR) {
        if (reqMode & NC_REQ_INDEP) return err;
        reqMode |= NC_REQ_ZERO; /* participate collective call */
    }

    status = nczipp->driver->put_var(nczipp->ncp, varid, start, count, stride, imap,
                                  cbuf, bufcount, buftype, reqMode);
    if (cbuf != buf) NCI_Free(cbuf);

    return (err == NC_NOERR) ? status : err; /* first error encountered */
}
