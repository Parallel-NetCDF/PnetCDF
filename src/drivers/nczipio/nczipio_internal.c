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

int nczipioi_var_init(NC_zip *nczipp, NC_zip_var *varp) {
    int i;
    int valid;
    MPI_Offset len, bsize;
    NC_zip_var *var;

    if (varp->type == NC_ZIP_VAR_COMPRESSED){
        if (varp->stripesize == NULL){    // This is a new uninitialized variable 
            // Determine its block size
            varp->stripesize = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim);
            
            // First check attribute
            valid = 1;
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_stripesize", NULL, &len);
            if (err == NC_NOERR && len == varp->ndim){
                err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_stripesize", varp->stripesize, MPI_UNSIGNED_LONG_LONG);
                if (err != NC_NOERR){
                    valid = 0;
                }
                //stripesize must be at leasst 1
                for(j = 0; j < varp->ndim; j++){ 
                    if (varp->stripesize[j] <= 0){
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
                memcpy(varp->stripesize, varp->size, sizeof(MPI_Offset) * varp->ndim);
                err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_stripesize", varp->stripesize, MPI_UNSIGNED_LONG_LONG);
                if (err != NC_NOERR){
                    return err;
                }
            }

            // Calculate block size
            bsize = 1;
            for(i = 0; i < nblocks; i++){
                bsize *= varp->stripesize[i];
            }

            // Calculate number of blocks
            varp->nblocks = 1;
            len = 1;
            for(j = 0; j < varp->ndim; j++){ //stripesize must be at leasst 1
                if (varp->size[j] % varp->stripesize[j] == 0){
                    varp->nblocks *= varp->size[j] / varp->stripesize[j];
                }
                else{
                    varp->nblocks *= varp->size[j] / varp->stripesize[j] + 1;
                }
                len *= varp->stripesize[j];   // Block size
            }

            varp->owner = (int*)NCI_Malloc(sizeof(int) * varp->nblocks);
            // Determine block ownership
            if (nczipp->blockmapping == NC_ZIP_MAPPING_STATIC){
                for(j = 0; j < varp->nblocks; j++){ 
                    varp->owner[j] = j % nczipp->np;
                    if (varp->owner[j] == nczipp->rank){
                        varp->cache[j] = (void*)NCI_Malloc(bsize);  // Allocate buffer for blocks we own
                    }
                }
            }

            // Determine block offset
            varp->offset = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->nblocks + 1);
            // Try if there are offset recorded in attributes, it can happen after opening a file
            err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_offsets", NULL, &len);
            if (err == NC_NOERR && varp->nblocks == len - 1){
                err = nczipp->driver->inq_att(nczipp->ncp, varp->varid, "_offsets", varp->offset, &len);
            }
            // If not, -1 indicates no data avaiable
            if (err != NC_NOERR{
                memset(varp->offset, -1, sizeof(MPI_Offset) * varp->nblocks);
            }
        }   
    }

    return NC_NOERR;
}