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

int nczipioi_var_init(NC_zip *nczipp) {
    int i;
    int valid;
    MPI_Offset len;
    NC_zip_var *var;

    for(i = 0; i < list->cnt; i++){
        var = nczipp->vars.data + i;

        if (var->type == NC_ZIP_VAR_COMPRESSED){
            if (var->blocksize == NULL){    // This is a new variable after enddef

                // Determine its block size
                var->blocksize = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * var->ndim);
                
                // First check hint
                valid = 1;
                err = nczipp->driver->inq_att(nczipp->ncp, var->varid, "_blocksize", NULL, &len);
                if (err == NC_NOERR && len == var->ndim){
                    err = nczipp->driver->get_att(nczipp->ncp, var->varid, "_blocksize", var->blocksize, MPI_UNSIGNED_LONG_LONG);
                    if (err != NC_NOERR){
                        valid = 0;
                    }
                    //Blocksize must be at leasst 1
                    for(j = 0; j < var->ndim; j++){ 
                        if (var->blocksize[j] <= 0){
                            valid = 0;
                            break;
                        }
                    }
                }
                else{
                    valid = 0;
                }
                
                // Default block size is same as dim size, only 1 blocks
                if (!valid){
                    memcpy(var->blocksize, var->size, sizeof(MPI_Offset) * var->ndim);
                }

                // Calculate number of blocks
                var->nblocks = 1;
                for(j = 0; j < var->ndim; j++){ //Blocksize must be at leasst 1
                    if (var->size[j] % var->blocksize[j] == 0){
                        var->nblocks *= var->size[j] / var->blocksize[j];
                    }
                    else{
                        var->nblocks *= var->size[j] / var->blocksize[j] + 1;
                    }
                }

                var->owner = (int*)NCI_Malloc(sizeof(int) * var->nblocks);
                var->cache = (void**)NCI_Malloc(sizeof(void*) * var->nblocks);

                // Determine block ownership
                for(j = 0; j < var->nblocks; j++){ 
                    var->owner[j] = j % nczipp->np;
                    if (var->owner[j] == nczipp->rank){
                        
                    }
                }
            }   
        }
    }
    return NC_NOERR;
}