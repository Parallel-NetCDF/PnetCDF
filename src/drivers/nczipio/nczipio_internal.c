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

int nczipioi_var_init(NC_zip *nczipp, NC_zip_var *varp) {
    int i, j, err;
    int valid;
    MPI_Offset len;
    NC_zip_var *var;

    if (varp->varkind == NC_ZIP_VAR_COMPRESSED){
        if (varp->chunkdim == NULL){    // This is a new uninitialized variable 
            // Determine its block size
            varp->chunkdim = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
            
            // First check attribute
            valid = 1;
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_chunkdim", NULL, &len);
            if (err == NC_NOERR && len == varp->ndim){
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_chunkdim", varp->chunkdim, MPI_UNSIGNED_LONG_LONG);
                if (err != NC_NOERR){
                    valid = 0;
                }
                //chunkdim must be at leasst 1
                for(j = 0; j < varp->ndim; j++){ 
                    if (varp->chunkdim[j] <= 0){
                        valid = 0;
                        printf("Warning: block size invalid, use default");
                        break;
                    }
                }
            }
            else{
                valid = 0;
            }
            
            // Default block size is same as dim size, only 1 blocks
            if (!valid){
                memcpy(varp->chunkdim, varp->dimsize, sizeof(MPI_Offset) * varp->ndim);
                err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_chunkdim", NC_INT, varp->ndim, varp->chunkdim, MPI_UNSIGNED_LONG_LONG);
                if (err != NC_NOERR){
                    return err;
                }
            }

            // Calculate block size
            varp->chunksize = 1;
            for(i = 0; i < varp->ndim; i++){
                varp->chunksize *= varp->chunkdim[i];
            }

            // Calculate number of blocks
            varp->nchunks = 1;
            len = 1;
            for(j = 0; j < varp->ndim; j++){ //chunkdim must be at leasst 1
                if (varp->chunkdim[j] % varp->chunkdim[j] == 0){
                    varp->nchunks *= varp->chunkdim[j] / varp->chunkdim[j];
                }
                else{
                    varp->nchunks *= varp->chunkdim[j] / varp->chunkdim[j] + 1;
                }
                len *= varp->chunkdim[j];   // Block size
            }

            // Determine block ownership
            varp->chunk_owner = (int*)NCI_Malloc(sizeof(int) * varp->nchunks);
            varp->chunk_cache = (char**)NCI_Malloc(sizeof(char*) * varp->nchunks);
            memset(varp->chunk_cache, 0, sizeof(char*) * varp->nchunks);
            varp->nmychunk = 0;
            if (nczipp->blockmapping == NC_ZIP_MAPPING_STATIC){
                for(j = 0; j < varp->nchunks; j++){ 
                    varp->chunk_owner[j] = j % nczipp->np;
                    if (varp->chunk_owner[j] == nczipp->rank){
                        varp->chunk_cache[j] = (void*)NCI_Malloc(varp->chunksize);  // Allocate buffer for blocks we own
                    }
                    varp->nmychunk++;
                }
            }

            // Build skip list of my own chunks
            varp->mychunks = (int*)NCI_Malloc(sizeof(int) * varp->nmychunk);
            varp->nmychunk = 0;
            for(j = 0; j < varp->nchunks; j++){ 
                if (varp->chunk_owner[j] == nczipp->rank){
                    varp->mychunks[varp->nmychunk++] = j;
                }
            }

            // Determine block offset
            varp->data_offs = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->nchunks);
            varp->data_lens = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->nchunks);
            // Try if there are offset recorded in attributes, it can happen after opening a file
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_offsets", NULL, &len);
            if (err == NC_NOERR && varp->nchunks == len - 1){
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_offsets", varp->data_offs, MPI_UNSIGNED_LONG_LONG);
            }
            // If not, 0 len no data avaiable
            if (err != NC_NOERR){
                memset(varp->data_offs, 0, sizeof(MPI_Offset) * varp->nchunks);
                memset(varp->data_lens, 0, sizeof(MPI_Offset) * varp->nchunks);
            }
        }   
    }

    return NC_NOERR;
}