/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>
#include <ctype.h>  /* isspace() */
#include <assert.h>
#include <mpi.h>

#include <pnetcdf.h>
#include <pnc_debug.h>
#include <common.h>

int hash_map_init(hash_map *map, int size, unsigned int (*hash)(const char* key)){
    map->hash = hash;
    map->table = (hash_map_node**)NCI_Calloc(size, sizeof(hash_map_node*));
    if (map->table == NULL){
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }
    map->size = size;

    return NC_NOERR;
}

int hash_map_free(hash_map *map) {
    int i;
    hash_map_node *nxt, *cur;

    /* Free the table */
    for(i = 0; i < map->size; i++){
        cur = map->table[i];
        while(cur != NULL){
            nxt = cur->next;
            NCI_Free(cur->key);
            NCI_Free(cur);
            cur = nxt;
        }
    }
    NCI_Free(map->table);

    return NC_NOERR;
}

int hash_map_add(hash_map *map, char *key, int val) {
    unsigned int idx;
    hash_map_node *pre = NULL, *cur;
    hash_map_node *new_node;

    /* Calculate has value */
    idx = map->hash(key) % map->size;

    /* Check if exist */
    cur = map->table[idx];
    while(cur != NULL){
        if (strcmp(key, cur->key) == 0){
            return NC_EEXIST;
        }
        pre = cur;
        cur = cur->next;
    }

    /* Create new node */
    new_node = (hash_map_node*)NCI_Malloc(sizeof(hash_map_node));
    if (new_node == NULL){
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }
    new_node->key = (char*)NCI_Malloc((strlen(key) + 1) * sizeof(char));
    if (new_node->key == NULL){
        NCI_Free(new_node);
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }
    strcpy(new_node->key, key);
    new_node->val = val;

    if (pre == NULL){
        map->table[idx] = new_node;
    }
    else{
        pre->next = new_node;
    }

    return NC_NOERR;
}


int hash_map_find(hash_map *map, char *key, int *val) {
    unsigned int idx;
    hash_map_node *cur;

    /* Calculate has value */
    idx = map->hash(key) % map->size;

    /* Check if exist */
    cur = map->table[idx];
    while(cur != NULL){
        if (strcmp(key, cur->key) == 0){
            *val = cur->val;
            return NC_NOERR;
        }
        cur = cur->next;
    }

    return NC_ENOTFOUND;
}
