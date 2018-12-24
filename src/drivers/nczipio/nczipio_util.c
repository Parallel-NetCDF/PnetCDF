/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <ncbbio_driver.h>

/*
 * Extract mpi hints and set up the flags
 */
void nczipioi_extract_hint(NC_zip *nczipp, MPI_Info info){
    int flag;
    char value[MPI_MAX_INFO_VAL];

    // Compression alg
    MPI_Info_get(info, "nc_zip_method", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag) {
        if (strcmp(value, "dummy") == 0){
            nczipp->zipdriver = NC_ZIP_DRIVER_DUMMY;  
        }
        else{
            printf("Warning: Unknown zip method %s, using dummy\n", value);
        }
    }
    else {
        nczipp->zipdriver = NC_ZIP_DRIVER_DUMMY;    
    }

    // Block assignment
    MPI_Info_get(info, "nc_zip_block_mapping", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag) {
        if (strcmp(value, "static") == 0){
            nczipp->blockmapping = NC_ZIP_MAPPING_STATIC;  
        }
        else{
            printf("Warning: Unknown zip method %s, using dummy\n", value);
        }
    }
    else {
        nczipp->blockmapping = NC_ZIP_MAPPING_STATIC;    
    }
}

/*
 * Export hint based on flag
 * NOTE: We only set up the hint if it is not the default setting
 *       user hint maching the default behavior will be ignored
 */
void nczipioi_export_hint(NC_zip *nczipp, MPI_Info info){
    char value[MPI_MAX_INFO_VAL];

    MPI_Info_set(info, "nc_compression", "enable");
    
    if (nczipp->zipdriver == NC_ZIP_DRIVER_DUMMY) {
        MPI_Info_set(info, "nc_zip_method", "dummy");
    }
}
