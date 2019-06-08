/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncadios_driver.h>
#include <ncadios_internal.h>

int ncadiosi_var_list_init(NC_ad_var_list *list) {
    list->cnt = 0;
    list->nalloc = 0;
    return NC_NOERR;
}

int ncadiosi_dim_list_init(NC_ad_dim_list *list) {
    list->cnt = 0;
    list->nalloc = 0;
    return NC_NOERR;
}

int ncadiosi_att_list_init(NC_ad_att_list *list) {
    list->cnt = 0;
    list->nalloc = 0;
    return NC_NOERR;
}

int ncadiosi_var_list_free(NC_ad_var_list *list) {
    int i;
    if (list->nalloc > 0){
        for(i = 0; i < list->cnt; i++){
            NCI_Free(list->data[i].name);
            NCI_Free(list->data[i].dimids);
            ncadiosi_att_list_free(&(list->data[i].atts));
        }
        NCI_Free(list->data);
    }
    return NC_NOERR;
}

int ncadiosi_dim_list_free(NC_ad_dim_list *list) {
    int i;
    if (list->nalloc > 0){
        for(i = 0; i < list->cnt; i++){
            NCI_Free(list->data[i].name);
        }
        NCI_Free(list->data);
    }
    return NC_NOERR;
}

int ncadiosi_att_list_free(NC_ad_att_list *list) {
    if (list->nalloc > 0){
        NCI_Free(list->data);
    }
    return NC_NOERR;
}

int ncadiosi_var_list_add(NC_ad_var_list *list, NC_ad_var data) {
    int id;

    id = list->cnt;

    if (list->nalloc == 0){
        list->nalloc = 16;
        list->data = NCI_Malloc(list->nalloc * sizeof(NC_ad_var));
    }
    else if (list->nalloc == id){
        list->nalloc *= 2;
        list->data = NCI_Realloc(list->data, list->nalloc * sizeof(NC_ad_var));
    }

    list->data[id] = data;
    list->cnt++;

    return id;
}

int ncadiosi_dim_list_add(NC_ad_dim_list *list, NC_ad_dim data) {
    int id;

    id = list->cnt;

    if (list->nalloc == 0){
        list->nalloc = 16;
        list->data = NCI_Malloc(list->nalloc * sizeof(NC_ad_dim));
    }
    else if (list->nalloc == id){
        list->nalloc *= 2;
        list->data = NCI_Realloc(list->data, list->nalloc * sizeof(NC_ad_dim));
    }

    list->data[id] = data;
    list->cnt++;

    return id;
}

int ncadiosi_att_list_add(NC_ad_att_list *list, int data) {
    int id;

    id = list->cnt;

    if (list->nalloc == 0){
        list->nalloc = 16;
        list->data = NCI_Malloc(list->nalloc * SIZEOF_INT);
    }
    else if (list->nalloc == id){
        list->nalloc *= 2;
        list->data = NCI_Realloc(list->data, list->nalloc * SIZEOF_INT);
    }

    list->data[id] = data;
    list->cnt++;

    return id;
}

int ncadiosi_var_list_find(NC_ad_var_list *list, char *name) {
    int i;

    for(i = 0; i < list->cnt; i++){
        if (strcmp(name, list->data[i].name) == 0){
            return i;
        }
    }

    return -1;
}

int ncadiosi_dim_list_find(NC_ad_dim_list *list, char *name) {
    int i;

    for(i = 0; i < list->cnt; i++){
        if (strcmp(name, list->data[i].name) == 0){
            return i;
        }
    }

    return -1;
}


int ncadiosi_att_list_find(NC_ad_att_list *list, int data) {
    int i;

    for(i = 0; i < list->cnt; i++){
        if (list->data[i] == data){
            return i;
        }
    }

    return -1;
}
