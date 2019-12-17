/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>
#include "nczipio_internal.h"
#include "../ncmpio/ncmpio_NC.h"

int
nczipioi_init(NC_zip *nczipp, int isnew){
    int err;

    nczipp->max_ndim = 0;
    nczipp->max_chunk_size = 0;
    nczipp->getsize = 0;
    nczipp->putsize = 0;
    nczipp->nmychunks = 0;
    nczipp->nwrite = 0;
    nczipp->cache_head = NULL;
    nczipp->cache_tail = NULL;
    nczipp->cache_used = 0;
    nczipp->cache_limit = 0;
    nczipp->cache_serial = 0;
    nczipp->ndim = 0;
    nczipp->chunkdim = NULL;

    err = nczipp->driver->inq(nczipp->ncp, NULL, NULL, NULL, &(nczipp->recdim));
    if (err != NC_NOERR) return err;

    if (isnew){
        nczipp->recsize = 0;
    }
    else{
        err = nczipp->driver->get_att(nczipp->ncp, NC_GLOBAL, "_recsize", &(nczipp->recsize), MPI_LONG_LONG); CHK_ERR // Mark this file as compressed
    }


    /* Initialize var list */
    err = nczipioi_var_list_init(&(nczipp->vars));
    if (err != NC_NOERR) return err;

    /* Initialize nonblocking list */
    err = nczipioi_req_list_init(&(nczipp->getlist));
    if (err != NC_NOERR) return err;
    err = nczipioi_req_list_init(&(nczipp->putlist));
    if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
    memset(&(nczipp->profile), 0, sizeof(NC_zip_timers));
    nczipp->sendsize = 0;
    nczipp->recvsize = 0;
    nczipp->nsend = 0;
    nczipp->nrecv = 0;
    nczipp->nremote = 0;
    nczipp->nreq = 0;
    nczipp->nlocal = 0;
#endif
}

int
nczipioi_parse_var_info(NC_zip *nczipp){
    int err;
    int vid;
    int i;
    int nvar;
    NC_zip_var var, *varp;
    
    int nread;
    int *lens;
    MPI_Aint *fdisps, *mdisps;
    MPI_Datatype ftype, mtype;
    MPI_Status status;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT_META)

    err = nczipp->driver->inq(nczipp->ncp, NULL, &nvar, NULL, &(nczipp->recdim));

    if (nvar > 0){
        lens = NCI_Malloc(sizeof(int) * nvar);
        fdisps = NCI_Malloc(sizeof(MPI_Aint) * nvar * 2);
        mdisps = fdisps + nvar;

        for(vid = 0; vid < nvar; vid++){
            memset(&var, 0, sizeof(var));
            
            err = nczipp->driver->get_att(nczipp->ncp, vid, "_varkind", &(var.varkind), MPI_INT);   // Comressed var?
            if (err != NC_NOERR || var.varkind == NC_ZIP_VAR_DATA){
                continue;
            }

            var.varid = vid;
            
            if (var.varkind == NC_ZIP_VAR_COMPRESSED){
                err = nczipp->driver->get_att(nczipp->ncp, var.varid, "_ndim", &(var.ndim), MPI_INT); // Original dimensions
                if (err != NC_NOERR) return err;

                var.dimids = (int*)NCI_Malloc(sizeof(int) * var.ndim);
                var.dimsize = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * var.ndim);

                err = nczipp->driver->get_att(nczipp->ncp, var.varid, "_dimids", var.dimids, MPI_INT);   // Dimensiona IDs
                if (err != NC_NOERR) return err;

                for(i = 0; i < var.ndim; i++){
                    nczipp->driver->inq_dim(nczipp->ncp, var.dimids[i], NULL, var.dimsize + i);
                }
                if (var.dimids[0] == nczipp->recdim){
                    var.isrec = 1;
                    if (var.dimsize[0] < nczipp->recsize){
                        var.dimsize[0] = nczipp->recsize;
                    }
                }
                else{
                    var.isrec = 0;
                }

                err = nczipp->driver->get_att(nczipp->ncp, var.varid, "_datatype", &(var.xtype), MPI_INT); // Original datatype
                if (err != NC_NOERR) return err;

                var.esize = NC_Type_size(var.xtype);
                var.etype = ncmpii_nc2mpitype(var.xtype);
                var.chunkdim = NULL;

                if (!(nczipp->delay_init)){
                    NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT_META)

                    nczipioi_var_init(nczipp, &var, 0, NULL, NULL);

                    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT_META)
                }
            }
        
            if (var.varkind == NC_ZIP_VAR_COMPRESSED || var.varkind == NC_ZIP_VAR_RAW){
                nczipioi_var_list_add(&(nczipp->vars), var);
            }
        }

        for(vid = 0; vid < nczipp->vars.cnt; vid++){
            varp = nczipp->vars.data + vid;
            err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_metaoffset", &(varp->metaoff), MPI_LONG_LONG);
            if (err == NC_NOERR){
                lens[nread] = sizeof(NC_zip_chunk_index_entry) * (varp->nchunk);
                fdisps[nread] = varp->metaoff;
                mdisps[nread++] = varp->chunk_index;
            }
            else{
                varp->metaoff = -1;;

                memset(varp->chunk_index, 0, sizeof(NC_zip_chunk_index_entry) * (varp->nchunk + 1));
            }
        }

        if (nread){
            nczipioi_sort_file_offset(nread, fdisps, mdisps, lens);

            MPI_Type_create_hindexed(nread, lens, fdisps, MPI_BYTE, &ftype);
            CHK_ERR_TYPE_COMMIT(&ftype);

            MPI_Type_create_hindexed(nread, lens, mdisps, MPI_BYTE, &mtype);
            CHK_ERR_TYPE_COMMIT(&mtype);

            // Set file view
            CHK_ERR_SET_VIEW(((NC*)(nczipp->ncp))->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
            
            // Read data
            CHK_ERR_READ_AT_ALL(((NC*)(nczipp->ncp))->collective_fh, 0, varp->chunk_index, 1, mtype, &status);
            
            // Restore file view
            CHK_ERR_SET_VIEW(ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

#ifndef WORDS_BIGENDIAN // Switch back to little endian
            //ncmpii_in_swapn(varp->chunk_index, varp->nchunk + 1, sizeof(long long));
            //ncmpii_in_swapn(varp->data_lens, varp->nchunk + 1, sizeof(int));
#endif

            MPI_Type_free(&ftype);
            MPI_Type_free(&mtype);
        }

        NCI_Free(lens);
        NCI_Free(fdisps);
    }

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT_META)

    return NC_NOERR;
}