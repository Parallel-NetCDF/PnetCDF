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
#include <common.h>
#include <nczipio_driver.h>
#include <zip_driver.h>

int nczip_dummy_init(MPI_Info info) {
    return NC_NOERR;
}

int nczip_dummy_finalize() {
    return NC_NOERR;
}

/* Return an estimated compressed data size
 * Actual compressed size should not exceed the estimation
 */
int nczip_dummy_inq_cpsize(void *in, int in_len, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    *out_len = in_len;
    return NC_NOERR;
}

/* If out_len is large enough, compress the data at in and save it to out. out_len is set to actual compressed data size
 * If out_len is NULL, we assume out is large enough for compressed data
 */
int nczip_dummy_compress(void *in, int in_len, void *out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    if (out_len != NULL){
        // Check output buffer size
        if ((*out_len) < in_len){
            DEBUG_RETURN_ERROR(NC_ENOMEM);
        }

        // Overwrite output buffer size with actual size
        *out_len = in_len;
    }
    
    // Copy data directly as dummy comrpession
    memcpy(out, in, in_len);

    return NC_NOERR;
}

/* Compress the data at in and save it to a newly allocated buffer at out. out_len is set to actual compressed data size
 * The caller is responsible to free the buffer
 * If out_len is not NULL, it will be set to buffer size allocated
 */
int nczip_dummy_compress_alloc(void *in, int in_len, void **out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    // Allocate output buffer
    *out = (void*)NCI_Malloc(in_len);

    // Buffer size
    if (out_len != NULL) {
        *out_len = in_len;
    }
    
    // Copy data directly as dummy comrpession
    memcpy(*out, in, in_len);

    return NC_NOERR;
}

/* Return an estimated decompressed data size
 * Actual decompressed size should not exceed the estimation
 */
int nczip_dummy_inq_dcsize(void *in, int in_len, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    *out_len = in_len;
    return NC_NOERR;
}

/* If out_len is large enough, decompress the data at in and save it to out. out_len is set to actual decompressed size
 * If out_len is NULL, we assume out is large enough for decompressed data
 */
int nczip_dummy_decompress(void *in, int in_len, void *out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    if (out_len != NULL){
        // Check output buffer size
        if ((*out_len) < in_len){
            DEBUG_RETURN_ERROR(NC_ENOMEM);
        }

        // Overwrite output buffer size with actual size
        *out_len = in_len;
    }
    
    // Copy data directly as dummy comrpession
    memcpy(out, in, in_len);

    return NC_NOERR;
}

/* Decompress the data at in and save it to a newly allocated buffer at out. out_len is set to actual decompressed data size
 * The caller is responsible to free the buffer
 * If out_len is not NULL, it will be set to buffer size allocated
 */
int nczip_dummy_decompress_alloc(void *in, int in_len, void **out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    // Allocate output buffer
    *out = (void*)NCI_Malloc(in_len);

    // Buffer size
    if (out_len != NULL) {
        *out_len = in_len;
    }
    
    // Copy data directly as dummy comrpession
    memcpy(*out, in, in_len);

    return NC_NOERR;
}

static NCZIP_driver nczipio_driver = {
    nczip_dummy_init,
    nczip_dummy_finalize,
    nczip_dummy_inq_cpsize,
    nczip_dummy_compress,
    nczip_dummy_compress_alloc,
    nczip_dummy_inq_dcsize,
    nczip_dummy_decompress,
    nczip_dummy_decompress_alloc
};

NCZIP_driver* nczip_dummy_inq_driver(void) {
    return &nczipio_driver;
}

