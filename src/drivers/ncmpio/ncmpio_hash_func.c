/*
 *  Copyright (C) 2016, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h> /* strlen() */
#include <assert.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< ncmpio_jenkins_one_at_a_time_hash() >--------------------------------*/
/* borrow Jenkins hash function:
 * https://en.wikipedia.org/wiki/Jenkins_hash_function
 */
int ncmpio_jenkins_one_at_a_time_hash(const char *str_name)
{
    unsigned int i, hash=0;
    for (i=0; i<strlen(str_name); ++i) {
        hash += (unsigned int)str_name[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);

#if 0
    ret = (int)hash; /* the return value will be used as an array index */
    return ((ret < 0) ? -ret : ret); /* make the value positive */
#endif
    /* this is to avoid expensive % operation, i.e. % HASH_TABLE_SIZE
    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
    */
    return (int)(hash & (HASH_TABLE_SIZE-1));
    /* return value will be used as an array index */
}

/*----< ncmpio_additive_hash() >---------------------------------------------*/
/* try different hash functions described in
 * http://www.burtleburtle.net/bob/hash/doobs.html
 */
int ncmpio_additive_hash(const char *str_name)
{
    size_t i, len = strlen(str_name);
    int hash = (int)len;
    for (i=0; i<len; ++i)
        hash += str_name[i]; /* additive hash */

    return (hash % 251); /* 251 is the largest prime <= 255 */
}

/*----< ncmpio_rotating_hash() >---------------------------------------------*/
int ncmpio_rotating_hash(const char *str_name)
{
    size_t i, len = strlen(str_name);
    unsigned int hash = (unsigned int)len;
    for (i=0; i<len; ++i)
        hash = (hash<<4)^(hash>>28)^(unsigned int)str_name[i];

    /* below is a clever way to replace (hash % prime) */
    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
}

/*----< ncmpio_Bernstein_hash() >--------------------------------------------*/
int ncmpio_Bernstein_hash(const char *str_name)
{
    size_t i, len = strlen(str_name);
    unsigned int hash = (unsigned int)len;
    for (i=0; i<len; ++i)
        /* hash = 65*hash+str_name[i]; */
        hash = hash+(hash<<6)+(unsigned int)str_name[i];

    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
}

/*----< ncmpio_Pearson_hash() >----------------------------------------------*/
int ncmpio_Pearson_hash(const char *str_name)
{
#if HASH_TABLE_SIZE == 256
    unsigned char T[256] = {
        251, 175, 119, 215, 81, 14, 79, 191, 103, 49, 181, 143, 186, 157,  0,
        232, 31, 32, 55, 60, 152, 58, 17, 237, 174, 70, 160, 144, 220, 90, 57,
        223, 59,  3, 18, 140, 111, 166, 203, 196, 134, 243, 124, 95, 222, 179,
        197, 65, 180, 48, 36, 15, 107, 46, 233, 130, 165, 30, 123, 161, 209, 23,
        97, 16, 40, 91, 219, 61, 100, 10, 210, 109, 250, 127, 22, 138, 29, 108,
        244, 67, 207,  9, 178, 204, 74, 98, 126, 249, 167, 116, 34, 77, 193,
        200, 121,  5, 20, 113, 71, 35, 128, 13, 182, 94, 25, 226, 227, 199, 75,
        27, 41, 245, 230, 224, 43, 225, 177, 26, 155, 150, 212, 142, 218, 115,
        241, 73, 88, 105, 39, 114, 62, 255, 192, 201, 145, 214, 168, 158, 221,
        148, 154, 122, 12, 84, 82, 163, 44, 139, 228, 236, 205, 242, 217, 11,
        187, 146, 159, 64, 86, 239, 195, 42, 106, 198, 118, 112, 184, 172, 87,
        2, 173, 117, 176, 229, 247, 253, 137, 185, 99, 164, 102, 147, 45, 66,
        231, 52, 141, 211, 194, 206, 246, 238, 56, 110, 78, 248, 63, 240, 189,
        93, 92, 51, 53, 183, 19, 171, 72, 50, 33, 104, 101, 69, 8, 252, 83, 120,
        76, 135, 85, 54, 202, 125, 188, 213, 96, 235, 136, 208, 162, 129, 190,
        132, 156, 38, 47, 1, 7, 254, 24, 4, 216, 131, 89, 21, 28, 133, 37, 153,
        149, 80, 170, 68, 6, 169, 234, 151
    };
    size_t i, len=strlen(str_name);
    unsigned char hash = (unsigned char)len;
    for (i=len; i>0; ) hash = T[hash ^ str_name[--i]];
    return (int)hash;
#else
    unsigned int i, hash=strlen(str_name);
    for (i=0; i<strlen(str_name); ++i)
        hash ^= str_name[i];

    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
#endif
}

/*----< ncmpio_update_name_lookup_table() >----------------------------------*/
/* remove the entry in lookup table for oldname and add a new entry for
 * newname
 */
int
ncmpio_update_name_lookup_table(NC_nametable *nameT,
                                const int     id,
                                const char   *oldname,  /*    normalized */
                                const char   *unewname) /* un-normalized */
{
    int i, key, err;
    char *name; /* normalized name string */

    /* remove the old name from the lookup table
     * hash the var name into a key for name lookup
     */
    key = HASH_FUNC(oldname);
    for (i=0; i<nameT[key].num; i++) {
        if (nameT[key].list[i] == id) break;
    }
    assert(i!=nameT[key].num);

    /* coalesce the id array */
    for (; i<nameT[key].num-1; i++)
        nameT[key].list[i] = nameT[key].list[i+1];

    /* decrease the number of IDs and free space if necessary */
    nameT[key].num--;
    if (nameT[key].num == 0) {
        NCI_Free(nameT[key].list);
        nameT[key].list = NULL;
    }

    /* normalized version of uname */
    err = ncmpii_utf8_normalize(unewname, &name);
    if (err != NC_NOERR) return err;

    /* hash the var name into a key for name lookup */
    key = HASH_FUNC(name);
    NCI_Free(name);

    /* add the new name to the lookup table
     * Note unewname must have already been checked for existence
     */
    if (nameT[key].num % NC_NAME_TABLE_CHUNK == 0)
        nameT[key].list = (int*) NCI_Realloc(nameT[key].list,
                          (size_t)(nameT[key].num+NC_NAME_TABLE_CHUNK) * SIZEOF_INT);
    nameT[key].list[nameT[key].num] = id;
    nameT[key].num++;

    return NC_NOERR;
}

/*----< ncmpio_hash_insert() >-----------------------------------------------*/
void
ncmpio_hash_insert(NC_nametable *nameT, /* var name lookup table */
                   const char   *name,
                   int           id)
{
    int key;

    /* hash the name into a key for name lookup */
    key = HASH_FUNC(name);

    /* allocate or expand the space for nameT[key].list */
    if (nameT[key].num % NC_NAME_TABLE_CHUNK == 0)
        nameT[key].list = (int*) NCI_Realloc(nameT[key].list,
                          (size_t)(nameT[key].num+NC_NAME_TABLE_CHUNK) *
                          SIZEOF_INT);

    /* add the ID to the name lookup table */
    nameT[key].list[nameT[key].num] = id;
    nameT[key].num++;
}

/*----< ncmpio_hash_delete() >-----------------------------------------------*/
/* only attributes can be deleted in NetCDF */
int
ncmpio_hash_delete(NC_nametable *nameT, /* var name lookup table */
                   const char   *name,
                   int           id)
{
    int i, j, key;

    /* hash the name into a key for name lookup */
    key = HASH_FUNC(name);

    /* find the entry from list[] that matches id */
    for (i=0; i<nameT[key].num; i++)
        if (nameT[key].list[i] == id)
            break;

    /* name is not found in nameT hash table */
    if (i == nameT[key].num) DEBUG_RETURN_ERROR(NC_ENOTATT)

    /* coalesce list[] */
    for (; i<nameT[key].num-1; i++)
        nameT[key].list[i] = nameT[key].list[i+1];

    nameT[key].num--;

    if (nameT[key].num == 0) {
        NCI_Free(nameT[key].list);
        nameT[key].list = NULL;
    }

    /* update all IDs that are > id */
    for (i=0; i<HASH_TABLE_SIZE; i++) {
        if (nameT[i].num == 0) continue;
        for (j=0; j<nameT[i].num; j++)
            if (nameT[i].list[j] > id)
                nameT[i].list[j]--;
    }

    return NC_NOERR;
}

/*----< ncmpio_hash_replace() >----------------------------------------------*/
/* remove old_name entry and add new_name entry */
int
ncmpio_hash_replace(NC_nametable *nameT, /* var name lookup table */
                    const char   *old_name,
                    const char   *new_name,
                    int           id)
{
    int i, key;

    /* hash the old name into a key for name lookup */
    key = HASH_FUNC(old_name);

    /* find the entry from list[] that matches id */
    for (i=0; i<nameT[key].num; i++)
        if (nameT[key].list[i] == id)
            break;

    /* name is not found in nameT hash table */
    if (i == nameT[key].num) DEBUG_RETURN_ERROR(NC_ENOTATT)

    /* coalesce list[] */
    for (; i<nameT[key].num-1; i++)
        nameT[key].list[i] = nameT[key].list[i+1];

    nameT[key].num--;

    if (nameT[key].num == 0) {
        NCI_Free(nameT[key].list);
        nameT[key].list = NULL;
    }

    /* hash the new name into a key for name lookup */
    key = HASH_FUNC(new_name);

    /* allocate or expand the space for nameT[key].list */
    if (nameT[key].num % NC_NAME_TABLE_CHUNK == 0)
        nameT[key].list = (int*) NCI_Realloc(nameT[key].list,
                          (size_t)(nameT[key].num+NC_NAME_TABLE_CHUNK) *
                          SIZEOF_INT);

    /* add the ID to the name lookup table */
    nameT[key].list[nameT[key].num] = id;
    nameT[key].num++;

    return NC_NOERR;
}

/*----< ncmpio_hash_table_copy() >-------------------------------------------*/
void
ncmpio_hash_table_copy(NC_nametable       *dest,
                       const NC_nametable *src)
{
    int i;

    for (i=0; i<HASH_TABLE_SIZE; i++) {
        dest[i].num = src[i].num;
        dest[i].list = NULL;
        if (dest[i].num > 0) {
            dest[i].list = NCI_Malloc((size_t)dest[i].num * SIZEOF_INT);
            memcpy(dest[i].list, src[i].list, (size_t)dest[i].num * SIZEOF_INT);
        }
    }
}

/*----< ncmpio_hash_table_free() >-------------------------------------------*/
/* free space allocated for name lookup table */
void
ncmpio_hash_table_free(NC_nametable *nameT)
{
    int i;
    for (i=0; i<HASH_TABLE_SIZE; i++) {
        if (nameT[i].num > 0)
            NCI_Free(nameT[i].list);
        nameT[i].num = 0;
    }
}

/*----< ncmpio_hash_table_populate_NC_dim() >--------------------------------*/
void
ncmpio_hash_table_populate_NC_dim(NC_dimarray *dimsp)
{
    int i;
    NC_nametable *nameT = dimsp->nameT;

    /* initialize dim name lookup table -------------------------------------*/
    memset(nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

    /* populate name lookup table */
    for (i=0; i<dimsp->ndefined; i++) {
        /* hash the dim name into a key for name lookup */
        int key = HASH_FUNC(dimsp->value[i]->name);
        nameT = &dimsp->nameT[key];
        if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
            nameT->list = (int*) NCI_Realloc(nameT->list,
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) *SIZEOF_INT);
        nameT->list[nameT->num] = i;
        nameT->num++;
    }
}

/*----< ncmpio_hash_table_populate_NC_var() >--------------------------------*/
void
ncmpio_hash_table_populate_NC_var(NC_vararray *varsp)
{
    int i;
    NC_nametable *nameT = varsp->nameT;

    /* initialize var name lookup table -------------------------------------*/
    memset(nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

    /* populate name lookup table */
    for (i=0; i<varsp->ndefined; i++) {
        /* hash the var name into a key for name lookup */
        int key = HASH_FUNC(varsp->value[i]->name);
        nameT = &varsp->nameT[key];
        if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
            nameT->list = (int*) NCI_Realloc(nameT->list,
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) *SIZEOF_INT);
        nameT->list[nameT->num] = i;
        nameT->num++;
    }
}

/*----< ncmpio_hash_table_populate_NC_attr() >-------------------------------*/
void
ncmpio_hash_table_populate_NC_attr(NC *ncp)
{
    int i, j;
    NC_nametable *nameT;

    /* populate name lookup table of global attributes */
    memset(ncp->attrs.nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

    for (i=0; i<ncp->attrs.ndefined; i++) {
        /* hash the var name into a key for name lookup */
        int key = HASH_FUNC(ncp->attrs.value[i]->name);
        nameT = &ncp->attrs.nameT[key];
        if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
            nameT->list = (int*) NCI_Realloc(nameT->list,
                          (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) *SIZEOF_INT);
        nameT->list[nameT->num] = i;
        nameT->num++;
    }

    /* populate name lookup table of each variable's attributes */
    for (j=0; j<ncp->vars.ndefined; j++) {
        NC_var *varp = ncp->vars.value[j];
        memset(varp->attrs.nameT, 0, sizeof(NC_nametable) * HASH_TABLE_SIZE);

        for (i=0; i<varp->attrs.ndefined; i++) {
            /* hash the var name into a key for name lookup */
            int key = HASH_FUNC(varp->attrs.value[i]->name);
            nameT = &varp->attrs.nameT[key];
            if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
                nameT->list = (int*) NCI_Realloc(nameT->list,
                              (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) *SIZEOF_INT);
            nameT->list[nameT->num] = i;
            nameT->num++;
        }
    }
}
