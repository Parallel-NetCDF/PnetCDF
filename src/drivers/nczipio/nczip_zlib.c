/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <nczipio_driver.h>
#include <zip_driver.h>

#include "zlib.h"

int nczip_zlib_init(MPI_Info info) {
    return NC_NOERR;
}

int nczip_zlib_finalize() {
    return NC_NOERR;
}

/* If out is not NULL and out_len is large enough, compress the data and save it in out. out_len is set to actual compressed size
 * If out is NULL, out_len is assigned an estimate compreseed size given the input
 * If out_len is NULL, we assume out is large enough for compressed data
 */
int nczip_zlib_compress(void *in, int in_len, void *out, int *out_len, int ndim, MPI_Offset *dims, MPI_Datatype dtype) {
    if (out != NULL){
        // zlib struct
        z_stream defstream;
        defstream.zalloc = Z_NULL;
        defstream.zfree = Z_NULL;
        defstream.opaque = Z_NULL;
        defstream.avail_in = (uInt)(in_len); // input size
        defstream.next_in = (Bytef*)in_len; // input
        defstream.avail_out = (uInt)(*out_len); // output buffer size
        defstream.next_out = (Bytef *)out; // output buffer

        // the actual compression work.
        err = deflateInit(&defstream, Z_BEST_COMPRESSION);
        if (err < 0){
            printf("deflateInit fail: %d: %s\n", err, defstream.msg);
        }
        err = deflate(&defstream, Z_FINISH);
        if (err < 0){
            printf("deflate fail: %d: %s\n", err, defstream.msg);
        }
        err = deflateEnd(&defstream);
        if (err < 0){
            printf("deflateEnd fail: %d: %s\n", err, defstream.msg);
        }

        if (out_len != NULL){
            *out_len = defstream.total_out;
        }
    }
    else{
        if (out_len != NULL){
            *out_len = in_len;
        }
    }
    return NC_NOERR;
}

/* If out is not NULL and out_len is large enough, decompress the data and save it in out. out_len is set to actual decompressed size
 * If out is NULL, out_len is assigned an estimate decompreseed size given the input
 * If out_len is NULL, we assume out is large enough for decompressed data
 */
int nczip_zlib_decompress(void *in, int in_len, void *out, int *out_len, int ndim, MPI_Offset *dims, MPI_Datatype dtype) {
    if (out != NULL){
        // zlib struct
        z_stream infstream;
        infstream.zalloc = Z_NULL;
        infstream.zfree = Z_NULL;
        infstream.opaque = Z_NULL;
        infstream.avail_in = (unsigned long) in_len; // input size
        infstream.next_in = (Bytef *)in; // input
        infstream.avail_out = (uInt)(*out_len); // output buffer
        infstream.next_out = (Bytef *)D; // buffer size
        
        // the actual DE-compression work.
        err = inflateInit(&infstream);
        if (err < 0){
            printf("inflateInit fail: %d: %s\n", err, infstream.msg);
        }
        err = inflate(&infstream, Z_NO_FLUSH);
        if (err < 0){
            printf("inflate fail: %d: %s\n", err, infstream.msg);
        }
        err = inflateEnd(&infstream);
        if (err < 0){
            printf("inflateEnd fail: %d: %s\n", err, infstream.msg);
        }
        if (*size != infstream.total_out){
            printf("Size mismatch: origin: %d, decompress: %lu\n", *size, infstream.total_out);
        }

        if (out_len != NULL){
            *out_len = infstream.total_out;
        }
    }
    else{
        if (out_len != NULL){
            *out_len = in_len;
        }
    }
    return NC_NOERR;
}

static NCZIP_driver nczip_driver_zlib = {
    nczip_zlib_init,
    nczip_zlib_finalize,
    nczip_zlib_compress,
    nczip_zlib_decompress
};

NCZIP_driver* nczip_zlib_inq_driver(void) {
    return &nczip_driver_zlib;
}

