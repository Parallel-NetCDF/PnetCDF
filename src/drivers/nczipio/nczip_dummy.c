/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <zip_driver.h>

int nczip_dummy_init(MPI_Info info) {
    return NC_NOERR;
}

int nczip_dummy_finalize() {
    return NC_NOERR;
}

/* If out is not NULL and out_len is large enough, compress the data and save it in out. out_len is set to actual compressed size
 * If out is NULL, out_len is assigned an estimate compreseed size given the input
 * If out_len is NULL, we assume out is large enough for compressed data
 */
int nczip_dummy_compress(void *in, MP_Offset in_len, void *out, MP_Offset *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    if (out_len != NULL){
        *out_len = in_len;
    }
    if (out != NULL){
        memcpy(out, in, in_len);
    }
    return NC_NOERR;
}

/* If out is not NULL and out_len is large enough, decompress the data and save it in out. out_len is set to actual decompressed size
 * If out is NULL, out_len is assigned an estimate decompreseed size given the input
 * If out_len is NULL, we assume out is large enough for decompressed data
 */
int nczip_dummy_decompress(void *in, MP_Offset in_len, void *out, MP_Offset *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    if (out_len != NULL){
        if (*out_len < in_len){
            DEBUG_RETURN_ERROR(NC_EINVAL)
        }
        *out_len = in_len;
    }
    if (out != NULL){
        memcpy(out, in, in_len);
    }
    return NC_NOERR;
}



static NCZIP_driver nczipio_driver = {
    nczip_dummy_init,
    nczip_dummy_finalize,
    nczip_dummy_compress,
    nczip_dummy_decompress
};

NCZIP_driver* nczip_dummy_inq_driver(void) {
    return &nczipio_driver;
}

