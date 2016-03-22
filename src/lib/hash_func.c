/*
 *  Copyright (C) 2016, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include "ncconfig.h"
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h> /* strlen() */

#include "nc.h"

/* borrow Jenkins hash function:
 * https://en.wikipedia.org/wiki/Jenkins_hash_function
 */
int ncmpii_jenkins_one_at_a_time_hash(const char *str_name)
{
int ncmpii_Pearson_hash(const char *);
return ncmpii_Pearson_hash(str_name);
    unsigned int i, hash=0;
    for (i=0; i<strlen(str_name); ++i) {
        hash += str_name[i];
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
    /* this is to avoid expensive % operation, i.e. % HASH_TABLE_SIZE */
    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
    /* return value will be used as an array index */
}

/* try different hash functions described in
 * http://www.burtleburtle.net/bob/hash/doobs.html
 */
int ncmpii_aditive_hash(const char *str_name)
{
    int i, hash=strlen(str_name);
    for (i=0; i<strlen(str_name); ++i)
        hash += str_name[i]; /* additive hash */

    return (hash % 251); /* 251 is the largest prime <= 255 */
}

int ncmpii_rotating_hash(const char *str_name)
{
    unsigned int i, hash=strlen(str_name);
    for (i=0; i<strlen(str_name); ++i)
        hash = (hash<<4)^(hash>>28)^str_name[i];

    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
}

int ncmpii_Bernstein_hash(const char *str_name)
{
    unsigned int i, hash=strlen(str_name);
    for (i=0; i<strlen(str_name); ++i)
        /* hash = 65*hash+str_name[i]; */
        hash = hash+(hash<<6)+str_name[i];

    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
}

int ncmpii_Pearson_hash(const char *str_name)
{
    unsigned int i, hash=strlen(str_name);
    for (i=0; i<strlen(str_name); ++i)
        hash ^= str_name[i];

    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
}

