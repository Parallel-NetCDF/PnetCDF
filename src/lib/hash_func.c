/*
 *  Copyright (C) 2016, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
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
    /* this is to avoid expensive % operation, i.e. % HASH_TABLE_SIZE
    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
    */
    return (int)(hash & (HASH_TABLE_SIZE-1));
    /* return value will be used as an array index */
}

/* try different hash functions described in
 * http://www.burtleburtle.net/bob/hash/doobs.html
 */
int ncmpii_additive_hash(const char *str_name)
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

    /* below is a clever way to replace (hash % prime) */
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
    unsigned char hash = len;
    for (i=len; i>0; ) hash = T[hash ^ str_name[--i]];
    return (int)hash;
#else
    unsigned int i, hash=strlen(str_name);
    for (i=0; i<strlen(str_name); ++i)
        hash ^= str_name[i];

    return (int)((hash ^ (hash>>10) ^ (hash>>20)) & (HASH_TABLE_SIZE-1));
#endif
}

