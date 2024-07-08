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
#include <ncchkio_driver.h>
#include <ncchk_filter_driver.h>

int ncchk_dummy_init(MPI_Info info) {
    return NC_NOERR;
}

int ncchk_dummy_finalize() {
    return NC_NOERR;
}

/* Return an estimated compressed data size
 * Actual compressed size should not exceed the estimation
 */
int ncchk_dummy_inq_cpsize(void *in, int in_len, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    *out_len = in_len;
    return NC_NOERR;
}

/* If out_len is large enough, compress the data at in and save it to out. out_len is set to actual compressed data size
 * If out_len is NULL, we assume out is large enough for compressed data
 */
int ncchk_dummy_compress(void *in, int in_len, void *out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
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
int ncchk_dummy_compress_alloc(void *in, int in_len, void **out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    // Allocate output buffer
    *out = (void*)malloc(in_len);

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
int ncchk_dummy_inq_dcsize(void *in, int in_len, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    *out_len = in_len;
    return NC_NOERR;
}

/* If out_len is large enough, decompress the data at in and save it to out. out_len is set to actual decompressed size
 * If out_len is NULL, we assume out is large enough for decompressed data
 */
int ncchk_dummy_decompress(void *in, int in_len, void *out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
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
int ncchk_dummy_decompress_alloc(void *in, int in_len, void **out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    // Allocate output buffer
    *out = (void*)malloc(in_len);

    // Buffer size
    if (out_len != NULL) {
        *out_len = in_len;
    }
    
    // Copy data directly as dummy comrpession
    memcpy(*out, in, in_len);

    return NC_NOERR;
}

static NCCHK_filter ncchkio_driver = {
    ncchk_dummy_init,
    ncchk_dummy_finalize,
    ncchk_dummy_inq_cpsize,
    ncchk_dummy_compress,
    ncchk_dummy_compress_alloc,
    ncchk_dummy_inq_dcsize,
    ncchk_dummy_decompress,
    ncchk_dummy_decompress_alloc
};

NCCHK_filter* ncchk_dummy_inq_driver(void) {
    return &ncchkio_driver;
}

