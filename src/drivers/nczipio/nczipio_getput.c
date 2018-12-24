/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_def_var()                  : dispatcher->def_var()
 * ncmpi_inq_varid()                : dispatcher->inq_varid()
 * ncmpi_inq_var()                  : dispatcher->inq_var()
 * ncmpi_rename_var()               : dispatcher->rename_var()
 *
 * ncmpi_get_var<kind>()            : dispatcher->get_var()
 * ncmpi_put_var<kind>()            : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>()     : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>()     : dispatcher->put_var()
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 *
 * ncmpi_iget_var<kind>()           : dispatcher->iget_var()
 * ncmpi_iput_var<kind>()           : dispatcher->iput_var()
 * ncmpi_iget_var<kind>_<type>()    : dispatcher->iget_var()
 * ncmpi_iput_var<kind>_<type>()    : dispatcher->iput_var()
 *
 * ncmpi_buffer_attach()            : dispatcher->buffer_attach()
 * ncmpi_buffer_detach()            : dispatcher->buffer_detach()
 * ncmpi_bput_var<kind>_<type>()    : dispatcher->bput_var()
 *
 * ncmpi_get_varn_<type>()          : dispatcher->get_varn()
 * ncmpi_put_varn_<type>()          : dispatcher->put_varn()
 *
 * ncmpi_iget_varn_<type>()         : dispatcher->iget_varn()
 * ncmpi_iput_varn_<type>()         : dispatcher->iput_varn()
 * ncmpi_bput_varn_<type>()         : dispatcher->bput_varn()
 *
 * ncmpi_get_vard()                 : dispatcher->get_vard()
 * ncmpi_put_vard()                 : dispatcher->put_vard()
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

static int get_block_id(NC_var *varp, int* cord){
    int i, ret;
    
    ret = 0;
    for(i = 0; i < varp->ndim; i++){
        ret = ret * 
    }
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
    int *bstart, *bend;
    MPI_Offset **starts,

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    *bstart = NCI_Malloc(sizeof(int) * varp->ndim);
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
