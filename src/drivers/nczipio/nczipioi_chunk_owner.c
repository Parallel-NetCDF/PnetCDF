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

int nczipioi_calc_chunk_owner(NC_zip *nczipp, NC_zip_var *varp, int nreq, MPI_Offset **starts, MPI_Offset **counts, int fixed){
    int err;
    int i, j, k;
    int cid;   // Chunk iterator
    int req;
    MPI_Offset overlapsize;
    MPI_Offset *ostart, *osize;
    MPI_Offset *citr; // Bounding box for chunks overlapping my own write region
    MPI_Offset_int *ocnt, *ocnt_all;

    NC_ZIP_TIMER_START(NC_ZIP_TIMER_INIT_COWN)
    
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
            overlapsize = 1;
            for(i = 0; i < varp->ndim; i++){
                overlapsize *= osize[i];
            }
            ocnt[cid].osize += (int)overlapsize;
            if (ocnt[cid].osize > varp->chunksize){
                ocnt[cid].osize = varp->chunksize;
            }
        } while (nczipioi_chunk_itr_next_ex(varp, starts[req], counts[req], citr, &cid, ostart, osize));
    }
    for(i = 0; i < varp->nchunk; i++){
        ocnt[i].rank = nczipp->rank;
        ocnt[i].osize -= nczipp->nmychunks;   // Penality for load ballance
    }
    // Noise to break tie
    for(i = (nczipp->rank + varp->varid) % nczipp->np; i < varp->nchunk; i += nczipp->np){
        ocnt[i].osize++;
    }

    CHK_ERR_ALLREDUCE(ocnt, ocnt_all, varp->nchunk, MPI_2INT, MPI_MAXLOC, nczipp->comm);

    for(i = fixed; i < varp->nchunk; i++){
        varp->chunk_owner[i] = ocnt_all[i].rank;
    }

#ifdef PNETCDF_PROFILING
    {
        if (nczipp->rank == 0){
            MPI_Status stat;
            FILE *pfile;
            char *pprefix = getenv("PNETCDF_OWNER_PREFIX");
            char fname[1024], ppath[1024];

            strcpy(fname, nczipp->path);
            for(i = strlen(fname); i > 0; i--){
                if (fname[i] == '.'){
                    fname[i] = '\0';
                }
                else if (fname[i] == '\\' || fname[i] == '/'){
                    i++;
                    break;
                }
            }
            sprintf(ppath, "%s%s_owner.csv", pprefix, fname + i);
            pfile = fopen(ppath, "a");

            fprintf(pfile, "Var:, %d\n", varp->varid);
            fprintf(pfile, "Rank\\Chunk, ");
            for(j = 0; j < varp->nchunk; j++){
                fprintf(pfile, "%d, ", j);
            }
            fprintf(pfile, "\n0, ");
            for(j = 0; j < varp->nchunk; j++){
                fprintf(pfile, "%d, ", ocnt[j].osize);
            }
            fprintf(pfile, "\n");
            for(i = 1; i < nczipp->np; i++){
                MPI_Recv(ocnt_all, varp->nchunk, MPI_2INT, i, 0, nczipp->comm, &stat);
                fprintf(pfile, "%d, ", i);
                for(j = 0; j < varp->nchunk; j++){
                    fprintf(pfile, "%d, ", ocnt_all[j].osize);
                }
                fprintf(pfile, "\n");
            }

            fclose(pfile);
        }
        else{
            MPI_Send(ocnt, varp->nchunk, MPI_2INT, 0, 0, nczipp->comm);
        }
    }
#endif

    NCI_Free(ostart);
    NCI_Free(ocnt);

    NC_ZIP_TIMER_STOP(NC_ZIP_TIMER_INIT_COWN)

    return NC_NOERR;
}