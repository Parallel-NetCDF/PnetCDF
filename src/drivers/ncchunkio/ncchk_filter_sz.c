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
#include <ncchkio_driver.h>
#include <ncchk_filter_driver.h>
#include <common.h>

#include <sz.h>

static int mpi_to_sz_type(MPI_Datatype dtype){
    if (dtype == MPI_FLOAT){
        return SZ_FLOAT;
    }
    else if (dtype == MPI_DOUBLE){
        return SZ_DOUBLE;
    }
    else if (dtype == MPI_BYTE){
        return SZ_UINT8;
    }
    else if (dtype == MPI_CHAR){
        return SZ_INT8;
    }
    else if (dtype == MPI_SHORT){
        return SZ_INT16;
    }
    else if (dtype == MPI_UNSIGNED_SHORT){
        return SZ_UINT16;
    }
    else if (dtype == MPI_INT){
        return SZ_INT32;
    }
    else if (dtype == MPI_UNSIGNED){
        return SZ_UINT32;
    }
    else if (dtype == MPI_LONG_LONG){
        return SZ_INT64;
    }
    else if (dtype == MPI_UNSIGNED_LONG_LONG){
        return SZ_UINT64;
    }

    return -1;
}

int ncchk_sz_init(MPI_Info info) {
    sz_params sz;

    memset(&sz, 0, sizeof(sz_params));
    sz.sol_ID = SZ;
    sz.sampleDistance = 50;
    sz.quantization_intervals = 0;
    sz.max_quant_intervals = 65536;
    sz.predThreshold = 0.98;
    sz.szMode = SZ_BEST_COMPRESSION;
    sz.losslessCompressor = ZSTD_COMPRESSOR;
    sz.gzipMode = 1;
    sz.errorBoundMode = ABS;
    sz.absErrBound = 1E-3;
    sz.relBoundRatio = 1E-5;
    SZ_Init_Params(&sz);

    return NC_NOERR;
}

int ncchk_sz_finalize() {
    SZ_Finalize();

    return NC_NOERR;
}

/* Return an estimated compressed data size
 * Actual compressed size should not exceed the estimation
 */
int ncchk_sz_inq_cpsize(void *in, int in_len, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    return NC_ENOTSUPPORT;  // sz has no size estimation
}

/* If out_len is large enough, compress the data at in and save it to out. out_len is set to actual compressed data size
 * If out_len is NULL, we assume out is large enough for compressed data
 */
int ncchk_sz_compress(void *in, int in_len, void *out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    int err;
    int i;
    int szdtype;
    size_t r[4];
    size_t outsize;
    void *buf = NULL;

    szdtype = mpi_to_sz_type(dtype);
    if (szdtype < 0){
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto out;
    }

    for(i = 0; i < 4; i++){
        if (i < ndim){
            r[i] = dims[i];
        }
        else{
            r[i] = 0;
        }
    }
    for(i = 4; i < ndim; i++){
        r[3] *= dims[i];
    }

    buf = SZ_compress(szdtype, in, &outsize, 0, r[3], r[2], r[1], r[0]);
    
    if (out_len != NULL){
        // If buffer not large enough
        if (*out_len < outsize){
            DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
            goto out;
        }

        // Size of comrpessed data
        *out_len = outsize;
    }

    memcpy(out, buf, outsize);

out:
    if (buf != NULL){
        free(buf);
    }

    return err;
}

/* Compress the data at in and save it to a newly allocated buffer at out. out_len is set to actual compressed data size
 * The caller is responsible to free the buffer
 * If out_len is not NULL, it will be set to buffer size allocated
 */
int ncchk_sz_compress_alloc(void *in, int in_len, void **out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    int err;
    int i;
    int szdtype;
    size_t r[4];
    size_t outsize;
    void *buf = NULL;

    szdtype = mpi_to_sz_type(dtype);
    if (szdtype < 0){
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto out;
    }

    for(i = 0; i < 4; i++){
        if (i < ndim){
            r[i] = dims[i];
        }
        else{
            r[i] = 0;
        }
    }
    for(i = 4; i < ndim; i++){
        r[3] *= dims[i];
    }

    *out = SZ_compress(szdtype, in, &outsize, 0, r[3], r[2], r[1], r[0]);
    
    if (out_len != NULL){
        // Size of comrpessed data
        *out_len = outsize;
    }

out:
    if (buf != NULL){
        free(buf);
    }

    return err;
}

/* Return an estimated decompressed data size
 * Actual decompressed size should not exceed the estimation
 */
int ncchk_sz_inq_dcsize(void *in, int in_len, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    return NC_ENOTSUPPORT;  // sz has no size estimation
}

/* If out_len is large enough, decompress the data at in and save it to out. out_len is set to actual decompressed size
 * If out_len is NULL, we assume out is large enough for decompressed data
 */
int ncchk_sz_decompress(void *in, int in_len, void *out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    int err;
    int i;
    size_t r[4];
    int szdtype;
    int outsize;
    void *buf = NULL;

    szdtype = mpi_to_sz_type(dtype);
    if (szdtype < 0){
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto out;
    }

    MPI_Type_size(dtype, &outsize);
    for(i = 0; i < 4; i++){
        if (i < ndim){
            r[i] = dims[i];
            outsize *= dims[i];
        }
        else{
            r[i] = 0;
        }
    }
    for(i = 4; i < ndim; i++){
        r[3] *= dims[i];
        outsize *= dims[i];
    }

    buf = SZ_decompress(szdtype, in, (size_t)in_len, 0, r[3], r[2], r[1], r[0]);
    
    if (out_len != NULL){
        // If buffer not large enough
        if (*out_len < outsize){
            DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
            goto out;
        }

        // Size of comrpessed data
        *out_len = outsize;
    }

    memcpy(out, buf, outsize);

out:
    if (buf != NULL){
        free(buf);
    }

    return err;
}

/* Decompress the data at in and save it to a newly allocated buffer at out. out_len is set to actual decompressed data size
 * The caller is responsible to free the buffer
 * If out_len is not NULL, it will be set to buffer size allocated
 */
int ncchk_sz_decompress_alloc(void *in, int in_len, void **out, int *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    int err;
    int i;
    size_t r[4];
    int szdtype;
    int outsize;
    void *buf = NULL;

    szdtype = mpi_to_sz_type(dtype);
    if (szdtype < 0){
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto out;
    }

    MPI_Type_size(dtype, &outsize);
    for(i = 0; i < 4; i++){
        if (i < ndim){
            r[i] = dims[i];
            outsize *= dims[i];
        }
        else{
            r[i] = 0;
        }
    }
    for(i = 4; i < ndim; i++){
        r[3] *= dims[i];
        outsize *= dims[i];
    }

    *out = SZ_decompress(szdtype, in, (size_t)in_len, 0, r[3], r[2], r[1], r[0]);
    
    if (out_len != NULL){
        // Size of comrpessed data
        *out_len = outsize;
    }

out:
    if (buf != NULL){
        free(buf);
    }

    return err;
}

static NCCHK_filter ncchk_driver_sz = {
    ncchk_sz_init,
    ncchk_sz_finalize,
    ncchk_sz_inq_cpsize,
    ncchk_sz_compress,
    ncchk_sz_compress_alloc,
    ncchk_sz_inq_dcsize,
    ncchk_sz_decompress,
    ncchk_sz_decompress_alloc
};

NCCHK_filter* ncchk_sz_inq_driver(void) {
    return &ncchk_driver_sz;
}

