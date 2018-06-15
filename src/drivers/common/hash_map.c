#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <common.h>
#include <errno.h>

int hash_map_init(hash_map *map, int size, unsigned int (*hash)(const char* key)){
    map->hash = hash;
    map->table = (hash_map_node**)NCI_Malloc(sizeof(hash_map_node*) * size);
    map->size = size;
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
}

int hash_map_add(hash_map *map, char *key, int val) {
    int i;
    unsigned int idx;
    hash_map_node *pre = NULL, *cur;
    hash_map_node *new_node;

    /* Calculate has value */
    idx = map->hash(key) % map->size;

    /* Check if exist */
    cur = map->table[idx];
    while(cur != NULL){
        if (strcmp(key, cur->key) == 0){
            return -1;
        }
        pre = cur;
        cur = cur->next;
    }

    /* Create new node */
    new_node = (hash_map_node*)NCI_Malloc(sizeof(hash_map_node));
    new_node->key = (char*)NCI_Malloc((strlen(key) + 1) * sizeof(char));
    strcpy(new_node->key, key);
    new_node->val = val;

    if (pre == NULL){
        map->table[idx] = new_node;
    }
    else{
        pre->next = new_node;
    }

    return 0;
}


int hash_map_find(hash_map *map, char *key, int *val) {
    int i;
    unsigned int idx;
    hash_map_node *pre = NULL, *cur;

    /* Calculate has value */
    idx = map->hash(key) % map->size;

    /* Check if exist */
    cur = map->table[idx];
    while(cur != NULL){
        if (strcmp(key, cur->key) == 0){
            *val = cur->val;
            return 0;
        }
        pre = cur;
        cur = cur->next;
    }

    return -1;
}