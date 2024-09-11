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
#include <ncchkio_driver.h>
#include "ncchkio_internal.h"


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

MPI_Offset get_chunk_overlap(NC_chk_var *varp, MPI_Offset* cord, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *ostart, MPI_Offset *ocount){
    int i;
    MPI_Offset ret = varp->esize;

    for(i = 0; i < varp->ndim; i++){
        ostart[i] = max(start[i], cord[i]);
        ocount[i] = min(start[i] + count[i], cord[i] + varp->chunkdim[i]) - ostart[i];
        if (ocount[i] <= 0){
            ocount[i] = 0;
        }
        ret *= ocount[i];
    }

    return ret;
}

int get_chunk_id(NC_chk_var *varp, MPI_Offset *cord){
    int i, ret;
    
    ret = (int)(cord[0]) / varp->chunkdim[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->nchunks[i] + (int)(cord[i]) / varp->chunkdim[i];
    }

    return ret;
}

int get_chunk_itr(NC_chk_var *varp, int idx, MPI_Offset* cord){
    int i;

    for(i = varp->ndim - 1; i >= 0; i--){
        cord[i] = (idx % varp->nchunks[i]) * varp->chunkdim[i];
        idx /= varp->nchunks[i];
    }

    return 0;
}

int ncchkioi_chunk_itr_init(NC_chk_var *varp, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *citr, int *cid){
    int i;

    *cid = 0;
    for(i = 0; i < varp->ndim; i++){
        citr[i] = start[i] - (start[i] % varp->chunkdim[i]);
        *cid += citr[i] / varp->chunkdim[i] * varp->cidsteps[i];
    }

    return NC_NOERR;
}

int ncchkioi_chunk_itr_next(NC_chk_var *varp, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *citr, int *cid){
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

int ncchkioi_chunk_itr_init_ex(NC_chk_var *varp, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *citr, int *cid, MPI_Offset *ostart, MPI_Offset *ocount){
    int i;

    *cid = 0;
    for(i = 0; i < varp->ndim; i++){
        citr[i] = start[i] - (start[i] % varp->chunkdim[i]);
        *cid += citr[i] / varp->chunkdim[i] * varp->cidsteps[i];
        ostart[i] = start[i];
        ocount[i] = min(count[i], citr[i] + varp->chunkdim[i] - ostart[i]);
    }

    return NC_NOERR;
}

int ncchkioi_chunk_itr_next_ex(NC_chk_var *varp, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *citr, int *cid, MPI_Offset *ostart, MPI_Offset *ocount){
    int i, j;
    int nchk = 1;

    i = varp->ndim - 1;
    citr[i] += varp->chunkdim[i];

    (*cid)++;
    for(; i > 0; i--){
        if (citr[i] >= start[i] + count[i]){
            citr[i - 1] += varp->chunkdim[i - 1];
            ostart[i - 1] += ocount[i - 1];
            ocount[i - 1] = min(varp->chunkdim[i - 1], start[i - 1] + count[i - 1] - ostart[i - 1]);
            j = citr[i];
            citr[i] = start[i] - (start[i] % varp->chunkdim[i]);
            ostart[i] = start[i];
            ocount[i] = min(count[i], citr[i] + varp->chunkdim[i] - ostart[i]);
            *cid += varp->cidsteps[i - 1] - varp->cidsteps[i] * (j - citr[i]) / varp->chunkdim[i];
        }
        else{
            break;
        }
    }

    if (citr[0] >= start[0] + count[0]){
        return 0;
    }

    if (i == varp->ndim - 1){
        ostart[i] += ocount[i];
        ocount[i] = min(varp->chunkdim[i], start[i] + count[i] - ostart[i]);
        for (i++; i < varp->ndim; i++) {
            ostart[i] = start[i];
            ocount[i] = min(count[i], citr[i] + varp->chunkdim[i] - ostart[i]);
        }
    }

    return 1;
}