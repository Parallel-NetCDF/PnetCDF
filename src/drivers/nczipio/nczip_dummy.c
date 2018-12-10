/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <zip_driver.h>
#include <nczipio_driver.h>

int nczip_dummy_init(MPI_Info info) {
    return NC_NOERR;
}

int nczip_dummy_finalize() {
    return NC_NOERR;
}

int nczip_dummy_compress(void *in, MP_Offset in_len, void *out, MP_Offset *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    if (out_len != NULL){
        *out_len = in_len;
    }
    if (out != NULL){
        memcpy(out, in, in_len);
    }
    return NC_NOERR;
}

int nczip_dummy_decompress(void *in, MP_Offset in_len, void *out, MP_Offset *out_len, int ndim, int *dims, MPI_Datatype dtype) {
    if (out_len != NULL){
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

