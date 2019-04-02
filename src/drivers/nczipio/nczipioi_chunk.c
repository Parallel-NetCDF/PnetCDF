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

int get_chunk_idx(NC_zip_var *varp, int* cord){
    int i, ret;
    
    ret = cord[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->chunkdim[i - 1] + cord[i];
    }

    return ret;
}

int get_chunk_cord(NC_zip_var *varp, int idx, int* cord){
    int i;

    for(i = varp->ndim - 1; i >= 0; i--){
        cord[i] = idx % varp->chunkdim[i];
        idx /= varp->chunkdim[i];
    }

    return 0;
}

int get_chunk_overlap(NC_zip_var *varp, int* cord, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *ostart, MPI_Offset *ocount){
    int i, ret;

    for(i = 0; i < varp->ndim; i++){
        ostart[i] = max(start[i], cord[i] * varp->chunkdim[i]);
        ocount[i] = min(start[i] + count[i], (cord[i] + 1) * varp->chunkdim[i]) - ostart[i];
        if (ocount[i] < 0){
            ocount[i] = 0;
        }
    }

    return 0;
}

int get_chunk_overlap_str(NC_zip_var *varp, int* cord, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, int *ostart, int *ocount){
    int i, ret;

    for(i = 0; i < varp->ndim; i++){
        ostart[i] = max(start[i], cord[i] * varp->chunkdim[i]);
        ocount[i] = (min(start[i] + count[i] * stride[i], (cord[i] + 1) * varp->chunkdim[i]) - ostart[i]) / stride[i];
        if (ocount[i] < 0){
            ocount[i] = 0;
        }
    }

    return 0;
}

int nczipioi_chunk_itr_init(NC_zip_var *varp, MPI_Offset *start, MPI_Offset *count, int *cstart, int *cend, int *citr){
    int i;
    int nchk = 1;

    for(i = 0; i < varp->ndim; i++){
        cstart[i] = citr[i] = start[i] / varp->chunkdim[i];
        cend[i] = (start[i] + count[i] - 1) / varp->chunkdim[i] + 1;
        nchk *= (cend[i] - cstart[i] + 1);
    }

    return nchk;
}

int nczipioi_chunk_itr_init_str(NC_zip_var *varp, MPI_Offset *start, MPI_Offset *count, MPI_Offset *stride, int *cstart, int *cend, int *citr){
    int i;
    int nchk = 1;

    for(i = 0; i < varp->ndim; i++){
        cstart[i] = citr[i] = start[i] / varp->chunkdim[i];
        cend[i] = (start[i] + (count[i] - 1) * stride[i]) / varp->chunkdim[i] + 1;
        if (stride[i] > varp->chunkdim[i]){
            nchk *= count[i];
        }
        else{
            nchk *= (cend[i] - cstart[i] + 1);
        }
    }

    return nchk;
}

int nczipioi_chunk_itr_next(NC_zip_var *varp, MPI_Offset *start, MPI_Offset *count, int *cstart, int *cend, int *citr){
    int i;
    int nchk = 1;

    i = varp->ndim - 1;
    citr[i]++;
    for(; i > 0; i--){
        if (citr[i] >= cend[i]){
            citr[i - 1]++;
            citr[i] = cstart[i];
        }
        else{
            break;
        }
    }

    if (citr[0] >= cend[0]){
        return 0;
    }

    return 1;
}

int nczipioi_chunk_itr_next_str(NC_zip_var *varp, MPI_Offset *start, MPI_Offset *count, MPI_Offset *stride, int *cstart, int *cend, int *citr){
    int i;
    int nchk = 1;

    i = varp->ndim - 1;
    if (stride == NULL){
        citr[i]++;
    }
    else{
        if (stride[i] > varp->chunkdim[i]){
            citr[i] = (citr[i] * varp->chunkdim[i] - start[i]);
        }
        else{
            citr[i]++;
        }
    }

    for(i = varp->ndim - 1; i > 0; i--){
        if (citr[i] >= cend[i]){
            citr[i - 1]++;
            citr[i] = cstart[i];
        }
        else{
            break;
        }
    }

    if (citr[0] >= cend[0]){
        return 0;
    }

    return 1;
}

int nczipioi_chunk_itr_init_cord(NC_zip_var *varp, MPI_Offset *start, MPI_Offset *count, MPI_Offset *citr){
    int i;

    for(i = 0; i < varp->ndim; i++){
        citr[i] = start[i] - (start[i] % varp->chunkdim[i]);
    }

    return NC_NOERR;
}


int nczipioi_chunk_itr_next_cord(NC_zip_var *varp, MPI_Offset *start, MPI_Offset *count, MPI_Offset *citr){
    int i;
    int nchk = 1;

    i = varp->ndim - 1;
    citr[i] += varp->chunkdim[i];
    for(; i > 0; i--){
        if (citr[i] >= start[i] + count[i]){
            citr[i - 1] += varp->chunkdim[i - 1];
            citr[i] = start[i] - (start[i] % varp->chunkdim[i]);
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

int get_chunk_overlap_cord(NC_zip_var *varp, MPI_Offset* cord, const MPI_Offset *start, const MPI_Offset *count, MPI_Offset *ostart, MPI_Offset *ocount){
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

int get_chunk_idx_cord(NC_zip_var *varp, MPI_Offset *cord){
    int i, ret;
    
    ret = (int)(cord[0]) / varp->chunkdim[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->chunkdim[i - 1] + (int)(cord[i]) / varp->chunkdim[i];
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
