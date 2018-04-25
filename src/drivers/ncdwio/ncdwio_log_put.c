/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <sys/types.h>
#include <dirent.h>
#include <assert.h>
#include "ncx.h"
#include <limits.h>
#include <fcntl.h>
#include <errno.h>
#include <stdint.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pnc_debug.h>
#include <common.h>
#include <pnetcdf.h>
#include <ncdwio_driver.h>

/*
 * Prepare a single log entry to be write to log
 * Used by ncmpii_getput_varm
 * IN    ncdwp:    log structure to log this entry
 * IN    varp:    NC_var structure associate to this entry
 * IN    start: start in put_var* call
 * IN    count: count in put_var* call
 * IN    stride: stride in put_var* call
 * IN    bur:    buffer of data to write
 * IN    buftype:    buftype from upper layer
 * IN    packedsize:    size of buf in byte
 */
int ncdwio_log_put_var(NC_dw *ncdwp, int varid, const MPI_Offset start[],
                       const MPI_Offset count[], const MPI_Offset stride[],
                       void *buf, MPI_Datatype buftype, MPI_Offset *putsize){
    int err, i, dim, elsize;
    int itype;    /* Type used in log file */
    int *dimids;
    char *buffer;
#ifdef PNETCDF_PROFILING
    double t1, t2, t3, t4, t5;
#endif
    MPI_Offset esize, dataoff, recsize;
    MPI_Offset *Start, *Count, *Stride;
    MPI_Offset size;
    NC_dw_metadataentry *entryp;
    NC_dw_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    t1 = MPI_Wtime();
#endif

    /* Calculate data size */
    /* Get ndims */
    err = ncdwp->ncmpio_driver->inq_var(ncdwp->ncp, varid, NULL, NULL, &dim,
                                        NULL, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR){
        return err;
    }

    /* Calcalate submatrix size */
    MPI_Type_size(buftype, &elsize);
    size = (MPI_Offset)elsize;
    for(i = 0; i < dim; i++){
        size *= count[i];
    }

    /* Return size */
    if (putsize != NULL){
        *putsize = size;
    }

    /* Record largest entry size */
    if (ncdwp->maxentrysize < size){
        ncdwp->maxentrysize = size;
    }

    /* Update record dimension size if is record variable */
    /* Get dimids */
    dimids = NCI_Malloc(SIZEOF_INT * dim);
    err = ncdwp->ncmpio_driver->inq_var(ncdwp->ncp, varid, NULL, NULL, NULL,
                                        dimids, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR){
        return err;
    }
    /* Update recdimsize if first dim is unlimited */
    if (dim > 0 && dimids[0] == ncdwp->recdimid) {
        /* Dim size after the put operation */
        if (stride == NULL) {
            recsize = start[0] + count[0];
        }
        else {
            recsize = start[0] + (count[0] - 1) * stride[0] + 1;
        }
        if (recsize > ncdwp->recdimsize) {
            ncdwp->recdimsize = recsize;
        }
    }
    NCI_Free(dimids);

    /* Convert to log types
     * Log spec has different enum of types than MPI
     */
    if (buftype == MPI_CHAR) {   /* put_*_text */
        itype = NC_LOG_TYPE_TEXT;
    }
    else if (buftype == MPI_SIGNED_CHAR) {    /* put_*_schar */
        itype = NC_LOG_TYPE_SCHAR;
    }
    else if (buftype == MPI_UNSIGNED_CHAR) {    /* put_*_uchar */
        itype = NC_LOG_TYPE_UCHAR;
    }
    else if (buftype == MPI_SHORT) { /* put_*_ushort */
        itype = NC_LOG_TYPE_SHORT;
    }
    else if (buftype == MPI_UNSIGNED_SHORT) { /* put_*_ushort */
        itype = NC_LOG_TYPE_USHORT;
    }
    else if (buftype == MPI_INT) { /* put_*_int */
        itype = NC_LOG_TYPE_INT;
    }
    else if (buftype == MPI_UNSIGNED) { /* put_*_uint */
        itype = NC_LOG_TYPE_UINT;
    }
    else if (buftype == MPI_FLOAT) { /* put_*_float */
        itype = NC_LOG_TYPE_FLOAT;
    }
    else if (buftype == MPI_DOUBLE) { /* put_*_double */
        itype = NC_LOG_TYPE_DOUBLE;
    }
    else if (buftype == MPI_LONG_LONG_INT) { /* put_*_longlong */
        itype = NC_LOG_TYPE_LONGLONG;
    }
    else if (buftype == MPI_UNSIGNED_LONG_LONG) { /* put_*_ulonglong */
        itype = NC_LOG_TYPE_ULONGLONG;
    }
    else { /* Unrecognized type */
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }

    /* Prepare metadata entry header */

    /* Find out the location of data in datalog
     * Which is current possition in data log
     * Datalog descriptor should always points to the end of file
     * Position must be recorded first before writing
     */
    dataoff = (MPI_Offset)ncdwp->datalogsize;

    /* Size of metadata entry
     * Include metadata entry header and variable size additional data
     * (start, count, stride)
     */
    esize = sizeof(NC_dw_metadataentry) + dim * 3 * SIZEOF_MPI_OFFSET;
    /* Allocate space for metadata entry header */
    buffer = (char*)ncdwio_log_buffer_alloc(&(ncdwp->metadata), esize);
    entryp = (NC_dw_metadataentry*)buffer;
    entryp->esize = esize; /* Entry size */
    entryp->itype = itype; /* Variable type */
    entryp->varid = varid;  /* Variable id */
    entryp->ndims = dim;  /* Number of dimensions of the variable*/
    /* The size of data in bytes. The size that will be write to data log */
	entryp->data_len = size;
    entryp->data_off = dataoff;

    /* Determine the api kind of original call
     * If stride is NULL, we log it as a vara call, otherwise, a vars call
     * Upper layer translates var1 and var to vara  and vars
     */
    if (stride == NULL){
        entryp->api_kind = NC_LOG_API_KIND_VARA;
    }
    else{
        entryp->api_kind = NC_LOG_API_KIND_VARS;
    }

    /* Calculate location of start, count, stride in metadata buffer */
    Start = (MPI_Offset*)(buffer + sizeof(NC_dw_metadataentry));
    Count = Start + dim;
    Stride = Count + dim;

    /* Fill up start, count, and stride */
    memcpy(Start, start, dim * SIZEOF_MPI_OFFSET);
    memcpy(Count, count, dim * SIZEOF_MPI_OFFSET);
    if(stride != NULL){
        memcpy(Stride, stride, dim * SIZEOF_MPI_OFFSET);
    }
    else{
        memset(Stride, 0, dim * SIZEOF_MPI_OFFSET);
    }

    /* Increase number of entry
     * This must be the final step of a log record
     * Increasing num_entries marks the completion of the record
     */

    /* Increase num_entries in the metadata buffer */
    headerp = (NC_dw_metadataheader*)ncdwp->metadata.buffer;
    headerp->num_entries++;

    //We only increase datalogsize by amount actually write
    ncdwp->datalogsize += size;

    /* Record data size */
    ncdwio_log_sizearray_append(&(ncdwp->entrydatasize), entryp->data_len);
    // Record in index
    // Entry address must be relative as metadata buffer can be reallocated
    ncdwio_metaidx_add(ncdwp, (NC_dw_metadataentry*)((char*)entryp -
                       (char*)(ncdwp->metadata.buffer)));

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
#endif

    /* Writing to data log
     * Note: Metadata record indicate completion, so data must go first
     */

    /*
     * Write data log
     */
    err = ncdwio_bufferedfile_write(ncdwp->datalog_fd, buf, size);
    if (err != NC_NOERR){
        return err;
    }

#ifdef PNETCDF_PROFILING
    t3 = MPI_Wtime();
#endif

    /* Seek to the head of metadata
     * Note: EOF may not be the place for next entry after a flush
     * Note: metadata size will be updated after allocating metadata buffer
     *       space, substract esize for original location
     */

    err = ncdwio_sharedfile_seek(ncdwp->metalog_fd, ncdwp->metadata.nused - esize,
                           SEEK_SET);
    if (err != NC_NOERR){
        return err;
    }

    /* Write meta data log */
    err = ncdwio_sharedfile_write(ncdwp->metalog_fd, buffer, esize);
    if (err != NC_NOERR){
        return err;
    }

#ifdef PNETCDF_PROFILING
    t4 = MPI_Wtime();
#endif

    /* Overwrite num_entries
     * This marks the completion of the record
     */
    err = ncdwio_sharedfile_pwrite(ncdwp->metalog_fd, &headerp->num_entries,
                            SIZEOF_MPI_OFFSET, 56);
    if (err != NC_NOERR){
        return err;
    }

#ifdef PNETCDF_PROFILING
    t5 = MPI_Wtime();
    ncdwp->put_data_wr_time += t3 - t2;
    ncdwp->put_meta_wr_time += t4 - t3;
    ncdwp->put_num_wr_time += t5 - t4;
    ncdwp->total_time += t5 - t1;
    ncdwp->put_time += t5 - t1;

    ncdwp->total_data += size;
    ncdwp->total_meta += esize;
#endif

    return NC_NOERR;
}
