/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <common.h>
#include <pnc_debug.h>
#include <ncbbio_driver.h>

/*
 * Prepare a single log entry to be write to log
 * Used by ncmpii_getput_varm
 * IN    ncbbp:    log structure to log this entry
 * IN    varp:    NC_var structure associate to this entry
 * IN    start: start in put_var* call
 * IN    count: count in put_var* call
 * IN    stride: stride in put_var* call
 * IN    bur:    buffer of data to write
 * IN    buftype:    buftype from upper layer
 * IN    packedsize:    size of buf in byte
 */
int ncbbio_log_put_var(NC_bb *ncbbp, int varid, const MPI_Offset start[],
                       const MPI_Offset count[], const MPI_Offset stride[],
                       void *buf, MPI_Datatype buftype, MPI_Offset *putsize){
    int err, i, dim, elsize;
    int itype;    /* Type used in log file */
    char *buffer;
    PNC *pncp;
#ifdef PNETCDF_PROFILING
    double t1, t2, t3, t4, t5;
#endif
    MPI_Offset esize, dataoff, recsize;
    MPI_Offset *Start, *Count, *Stride;
    MPI_Offset size;
    NC_bb_metadataentry *entryp;
    NC_bb_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    t1 = MPI_Wtime();
#endif

    /* Check parameters
     * Varid must be valid
     * Start, count must be valid
     * ECHAR must be detected
     */

    /* Get PNC */
    err = PNC_check_id(ncbbp->ncid, &pncp);
    if (err != NC_NOERR){
        return err;
    }

    /* Get ndims */
    dim = pncp->vars[varid].ndims;

    /* Calcalate submatrix size */
    if (buftype != MPI_DATATYPE_NULL){
        MPI_Type_size(buftype, &elsize);
    }
    else{
        elsize = 0;
    }
    size = (MPI_Offset)elsize;
    for(i = 0; i < dim; i++){
        size *= count[i];
    }

    /* Return size */
    if (putsize != NULL){
        *putsize = size;
    }

    /* Skip empty entries
     * Other arguments form upper layer may be invalid in case of 0 size request
     * skip to prevent unnecessary error
     */
    if (size == 0){
        return NC_NOERR;
    }

    /* Record largest entry size
     * This is used later to determine minimal buffer size required to flush the log
     */
    if (ncbbp->maxentrysize < size){
        ncbbp->maxentrysize = size;
    }

    /* Update record dimension size if is record variable */
    if (pncp->vars[varid].recdim >= 0) {
        /* Dim size after the put operation */
        if (stride == NULL) {
            recsize = start[0] + count[0];
        }
        else {
            recsize = start[0] + (count[0] - 1) * stride[0] + 1;
        }
        if (recsize > ncbbp->recdimsize) {
            ncbbp->recdimsize = recsize;
        }
    }

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
    else if (buftype == MPI_LONG) { /* put_*_int */
        itype = NC_LOG_TYPE_LONG;
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
        int name_len;
        char type_name[MPI_MAX_OBJECT_NAME];
        MPI_Type_get_name(buftype, type_name, &name_len);
        fprintf(stderr, "Rank: %d, Unrecognized type: %s\n", ncbbp->rank, type_name); fflush(stderr);
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }

    /* Prepare metadata entry header */

    /* Find out the location of data in datalog
     * Which is current possition in data log
     * Datalog descriptor should always points to the end of file
     * Position must be recorded first before writing
     */
    dataoff = (MPI_Offset)ncbbp->datalogsize;

    /* Size of metadata entry
     * Include metadata entry header and variable size additional data
     * (start, count, stride)
     */
    esize = sizeof(NC_bb_metadataentry) + dim * 3 * SIZEOF_MPI_OFFSET;
    /* Allocate space for metadata entry header */
    buffer = (char*)ncbbio_log_buffer_alloc(&(ncbbp->metadata), esize);
    entryp = (NC_bb_metadataentry*)buffer;
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
    Start = (MPI_Offset*)(buffer + sizeof(NC_bb_metadataentry));
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
    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;
    headerp->num_entries++;

    //We only increase datalogsize by amount actually write
    ncbbp->datalogsize += size;

    /* Record data size */
    ncbbio_log_sizearray_append(&(ncbbp->entrydatasize), entryp->data_len);
    
    // Record index
    // NOTE: Entry address must be relative as metadata buffer can be reallocated
    ncbbio_metaidx_add(ncbbp, (NC_bb_metadataentry*)((char*)entryp -
                       (char*)(ncbbp->metadata.buffer)));

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
#endif

    /* Writing to data log
     * Note: Metadata record indicate completion, so data must go first
     */

    /*
     * Write data log
     */
    err = ncbbio_sharedfile_write(ncbbp->datalog_fd, buf, size);
    if (err != NC_NOERR){
        return err;
    }

#ifdef PNETCDF_PROFILING
    t3 = MPI_Wtime();
#endif

    /* Write meta data log */
    err = ncbbio_sharedfile_write(ncbbp->metalog_fd, buffer, esize);
    if (err != NC_NOERR){
        return err;
    }

#ifdef PNETCDF_PROFILING
    t4 = MPI_Wtime();
#endif

    /* Overwrite num_entries
     * This marks the completion of the record
     */
    err = ncbbio_sharedfile_pwrite(ncbbp->metalog_fd, &headerp->num_entries,
                            SIZEOF_MPI_OFFSET, 56);
    if (err != NC_NOERR){
        return err;
    }

#ifdef PNETCDF_PROFILING
    t5 = MPI_Wtime();
    ncbbp->put_data_wr_time += t3 - t2;
    ncbbp->put_meta_wr_time += t4 - t3;
    ncbbp->put_num_wr_time += t5 - t4;
    ncbbp->total_time += t5 - t1;
    ncbbp->put_time += t5 - t1;

    ncbbp->total_data += size;
    ncbbp->total_meta += esize;
#endif

    return NC_NOERR;
}
