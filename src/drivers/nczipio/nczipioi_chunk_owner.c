/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

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

typedef struct MPI_Offset_int {
    int osize;
    int rank;
} MPI_Offset_int;

int nczipioi_calc_chunk_owner(NC_zip *nczipp, NC_zip_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts){
    int i, j, k;
    int cid;   // Chunk iterator
    int req;
    MPI_Offset overlapsize;
    MPI_Offset *ostart, *osize;
    MPI_Offset *citr; // Bounding box for chunks overlapping my own write region
    MPI_Offset_int *ocnt, *ocnt_all;

    
    ostart = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->ndim * 3);
    osize = ostart + varp->ndim;
    citr = osize + varp->ndim;

    ocnt = (MPI_Offset_int*)NCI_Malloc(sizeof(MPI_Offset_int) * varp->nchunk * 2);
    ocnt_all = ocnt + varp->nchunk;
    memset(ocnt, 0, sizeof(MPI_Offset_int) * varp->nchunk);

    // Count overlapsize of each request
    for(req = 0; req < nreq; req++){
        nczipioi_chunk_itr_init_ex(varp, starts[req], counts[req], citr, &cid, ostart, osize); // Initialize chunk iterator
        do{
            // Count overlap
            //get_chunk_overlap(varp, citr, starts[req], counts[req], ostart, osize);
            overlapsize = 1;
            for(i = 0; i < varp->ndim; i++){
                overlapsize *= osize[i];
            }
            ocnt[cid].osize += (int)overlapsize + (varp->ndim * 2 + 2) * sizeof(int);
        } while (nczipioi_chunk_itr_next_ex(varp, starts[req], counts[req], citr, &cid, ostart, osize));
    }
    for(i = 0; i < varp->nchunk; i++){
        ocnt[i].rank = nczipp->rank;
        if (ocnt[cid].osize > varp->chunksize){
            ocnt[cid].osize = varp->chunksize;
        }
        ocnt[i].osize -= nczipp->nmychunks;   // Penality for load ballance
    }
    // Noise to break tie
    for(i = (nczipp->rank + varp->varid) % nczipp->np; i < varp->nchunk; i += nczipp->np){
        ocnt[i].osize++;
    }

    CHK_ERR_ALLREDUCE(ocnt, ocnt_all, varp->nchunk, MPI_2INT, MPI_MAXLOC, nczipp->comm);

    for(i = 0; i < varp->nchunk; i++){
        varp->chunk_owner[i] = ocnt_all[i].rank;
    }

    NCI_Free(ostart);
    NCI_Free(ocnt);

    return NC_NOERR;
}