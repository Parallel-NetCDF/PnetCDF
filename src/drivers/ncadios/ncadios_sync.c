/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This file implements the helper function to synchronise header information
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strlen() */

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncadios_driver.h>
#include <ncadios_internal.h>

int ncadios_sync_header(NC_ad *ncadp) {
    int i, j;
    int bsize;
    int namelen;
    char *buf, *cur;

    if (ncadp->rank == 0){
        bsize = SIZEOF_INT * 2;   /* natt and nvar */
        for(i = 0; i < ncadp->dims.cnt; i++){
            bsize += strlen(ncadp->dims.data[i].name) + 1 + SIZEOF_INT * 2;
        }
        for(i = 0; i < ncadp->vars.cnt; i++){
            bsize += strlen(ncadp->vars.data[i].name) + 1 +
            SIZEOF_INT * 3 + ncadp->vars.data[i].ndim * SIZEOF_INT +
            ncadp->vars.data[i].atts.cnt * SIZEOF_INT + sizeof(nc_type);
        }
    }

    MPI_Bcast(&bsize, 1, MPI_INT, 0, ncadp->comm);

    buf = NCI_Malloc(bsize);
    cur = buf;

    if (ncadp->rank == 0){
        *((int*)cur) = ncadp->dims.cnt;
        cur += SIZEOF_INT;
        for(i = 0; i < ncadp->dims.cnt; i++){
            *((int*)cur) = ncadp->dims.data[i].len;
            cur += SIZEOF_INT;
            namelen = strlen(ncadp->dims.data[i].name);
            *((int*)cur) = namelen;
            cur += SIZEOF_INT;
            strcpy(cur, ncadp->dims.data[i].name);
            cur += namelen + 1;
        }
        *((int*)cur) = ncadp->vars.cnt;
        cur += SIZEOF_INT;
        for(i = 0; i < ncadp->vars.cnt; i++){
            *((nc_type*)cur) = ncadp->vars.data[i].type;
            cur += sizeof(nc_type);
            *((int*)cur) = ncadp->vars.data[i].ndim;
            cur += SIZEOF_INT;
            *((int*)cur) = ncadp->vars.data[i].atts.cnt;
            cur += SIZEOF_INT;
            namelen = strlen(ncadp->vars.data[i].name);
            *((int*)cur) = namelen;
            cur += SIZEOF_INT;
            memcpy(cur, ncadp->vars.data[i].dimids,
                    ncadp->vars.data[i].ndim * SIZEOF_INT);
            cur += ncadp->vars.data[i].ndim * 4;
            memcpy(cur, ncadp->vars.data[i].atts.data,
                    ncadp->vars.data[i].atts.cnt * SIZEOF_INT);
            cur += ncadp->vars.data[i].atts.cnt * 4;
            strcpy(cur, ncadp->vars.data[i].name);
            cur += namelen + 1;
        }
    }

    MPI_Bcast(buf, bsize, MPI_BYTE, 0, ncadp->comm);


    if (ncadp->rank != 0){
        int ndim, nvar, natt;
        int id, len;
        nc_type type;
        int *dimids, *attids;
        char *name;

        ndim = *((int*)cur);
        cur += 4;
        for(i = 0; i < ndim; i++){
            len = *((int*)cur);
            cur += SIZEOF_INT;
            namelen = *((int*)cur);
            cur += SIZEOF_INT;
            name = cur;
            cur += namelen + 1;
            ncadiosi_def_dim(ncadp, name, len, &id);
        }


        nvar = *((int*)cur);
        cur += 4;
        for(i = 0; i < nvar; i++){
            type = *((nc_type*)cur);
            cur += sizeof(nc_type);
            ndim = *((int*)cur);
            cur += SIZEOF_INT;
            natt = *((int*)cur);
            cur += SIZEOF_INT;
            namelen = *((int*)cur);
            cur += SIZEOF_INT;
            dimids = (int*)cur;
            cur += ndim * SIZEOF_INT;
            attids = (int*)cur;
            cur += natt * SIZEOF_INT;
            name = cur;
            cur += namelen + 1;
            ncadiosi_def_var(ncadp, name, type, ndim, dimids, &id);
            for(j = 0; j < natt; j++){
                ncadiosi_att_list_add(&(ncadp->vars.data[id].atts), attids[j]);
            }
        }

    }

    NCI_Free(buf);

    return NC_NOERR;
}
