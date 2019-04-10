/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
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
#include <math.h>
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>
#include "nczipio_internal.h"


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

int get_chunk_overlap(NC_zip_var *varp, MPI_Offset* cord, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *ostart, MPI_Offset *ocount){
    int i, ret;

    for(i = 0; i < varp->ndim; i++){
        ostart[i] = max(start[i], cord[i]);
        ocount[i] = min(start[i] + count[i], cord[i] + varp->chunkdim[i]) - ostart[i];
        if (ocount[i] < 0){
            ocount[i] = 0;
        }
    }

    return 0;
}

int get_chunk_id(NC_zip_var *varp, MPI_Offset *cord){
    int i, ret;
    
    ret = (int)(cord[0]) / varp->chunkdim[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->nchunks[i - 1] + (int)(cord[i]) / varp->chunkdim[i];
    }

    return ret;
}

int get_chunk_itr(NC_zip_var *varp, int idx, MPI_Offset* cord){
    int i;

    for(i = varp->ndim - 1; i >= 0; i--){
        cord[i] = (idx % varp->chunkdim[i]) * varp->chunkdim[i];
        idx /= varp->chunkdim[i];
    }

    return 0;
}

int nczipioi_chunk_itr_init(NC_zip_var *varp, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *citr, int *cid){
    int i;

    *cid = 0;
    for(i = 0; i < varp->ndim; i++){
        citr[i] = start[i] - (start[i] % varp->chunkdim[i]);
        *cid += citr[i] / varp->chunkdim[i] * varp->cidsteps[i];
    }

    return NC_NOERR;
}

int nczipioi_chunk_itr_next(NC_zip_var *varp, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *citr, int *cid){
    int i, j;
    int nchk = 1;

    i = varp->ndim - 1;
    citr[i] += varp->chunkdim[i];
    (*cid)++;
    for(; i > 0; i--){
        if (citr[i] >= start[i] + count[i]){
            citr[i - 1] += varp->chunkdim[i - 1];
            j = citr[i];
            citr[i] = start[i] - (start[i] % varp->chunkdim[i]);
            *cid += varp->cidsteps[i - 1] - varp->cidsteps[i] * (j - citr[i]) / varp->chunkdim[i];
        }
        else{
            break;
        }
    }

    if (citr[0] >= start[0] + count[0]){
        return 0;
    }

    return 1;
}