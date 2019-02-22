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

int nczipioi_var_list_init(NC_zip_var_list *list) {
    list->cnt = 0;
    list->nalloc = 0;
    return NC_NOERR;
}

int nczipioi_var_list_free(NC_zip_var_list *list) {
    int i, j;
    if (list->nalloc > 0){
        for(i = 0; i < list->cnt; i++){
            if (list->data[i].dimsize != NULL){
                NCI_Free(list->data[i].dimsize);
            }
            if (list->data[i].chunkdim != NULL){
                NCI_Free(list->data[i].chunkdim);
            }
            if (list->data[i].dimids != NULL){
                NCI_Free(list->data[i].dimids);
            }
            if (list->data[i].nchunks >= 0){
                NCI_Free(list->data[i].data_offs);
                NCI_Free(list->data[i].chunk_cache);
                NCI_Free(list->data[i].chunk_owner);
                NCI_Free(list->data[i].data_lens);
                for(j = 0; j < list->data[i].nmychunk; j++){
                    NCI_Free(list->data[i].chunk_cache[list->data[i].mychunks[j]]);
                }
                NCI_Free(list->data[i].chunk_cache);
            }
        }
        NCI_Free(list->data);
    }
    return NC_NOERR;
}

int nczipioi_var_list_add(NC_zip_var_list *list, NC_zip_var data) {
    int id;

 //   return 0;

    id = list->cnt;

    if (list->nalloc == 0){
        list->nalloc = 16;
        list->data = NCI_Malloc(list->nalloc * sizeof(NC_zip_var));
    }
    else if (list->nalloc == id){
        list->nalloc *= 2;
        list->data = NCI_Realloc(list->data, list->nalloc * sizeof(NC_zip_var));
    }

    list->data[id] = data;
    list->cnt++;

    return id;
}
