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

/* Convert MPI datatype to ncbb log types
 * Log spec has different enum of types than MPI
 */
static int mpitype2logtype(MPI_Datatype buftype, int *logtype){
    if (buftype == MPI_CHAR)                     /* put_*_text */
        *logtype = NC_LOG_TYPE_TEXT;
    else if (buftype == MPI_SIGNED_CHAR)         /* put_*_schar */
        *logtype = NC_LOG_TYPE_SCHAR;
    else if (buftype == MPI_UNSIGNED_CHAR)       /* put_*_uchar */
        *logtype = NC_LOG_TYPE_UCHAR;
    else if (buftype == MPI_SHORT)               /* put_*_ushort */
        *logtype = NC_LOG_TYPE_SHORT;
    else if (buftype == MPI_UNSIGNED_SHORT)      /* put_*_ushort */
        *logtype = NC_LOG_TYPE_USHORT;
    else if (buftype == MPI_INT)                 /* put_*_int */
        *logtype = NC_LOG_TYPE_INT;
    else if (buftype == MPI_UNSIGNED)            /* put_*_uint */
        *logtype = NC_LOG_TYPE_UINT;
    else if (buftype == MPI_LONG)                /* put_*_int */
        *logtype = NC_LOG_TYPE_LONG;
    else if (buftype == MPI_FLOAT)               /* put_*_float */
        *logtype = NC_LOG_TYPE_FLOAT;
    else if (buftype == MPI_DOUBLE)              /* put_*_double */
        *logtype = NC_LOG_TYPE_DOUBLE;
    else if (buftype == MPI_LONG_LONG_INT)       /* put_*_longlong */
        *logtype = NC_LOG_TYPE_LONGLONG;
    else if (buftype == MPI_UNSIGNED_LONG_LONG)  /* put_*_ulonglong */
        *logtype = NC_LOG_TYPE_ULONGLONG;
    else { /* Unrecognized type */
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }

    return NC_NOERR;
}

/*
 * Prepare a single log entry to be write to log
 * IN    ncbbp:  log structure to log this entry
 * IN    varp:   NC_var structure associate to this entry
 * IN    start:  start in put_var* call
 * IN    count:  count in put_var* call
 * IN    stride: stride in put_var* call
 * IN    bur:    buffer of data to write
 * IN    buftype: user buffer's internal data type, an MPI primitive type
 */
int ncbbio_log_put_var(NC_bb            *ncbbp,
                       int               varid,
                       const MPI_Offset  start[],  /* must not be NULL */
                       const MPI_Offset  count[],  /* may be NULL */
                       const MPI_Offset  stride[],
                       void             *buf,
                       MPI_Datatype      buftype)
{
    int err, i, ndims, elsize, itype;
    char *buffer;
    PNC *pncp;
    MPI_Offset esize, dataoff, recsize, put_size;
    MPI_Offset *Start, *Count, *Stride;
    NC_bb_metadataentry *entryp;
    NC_bb_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    double t1, t2, t3, t4, t5;
    t1 = MPI_Wtime();
#endif

    /* Get PNC */
    err = PNC_check_id(ncbbp->ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* Get ndims */
    ndims = pncp->vars[varid].ndims;

    /* Calcalate put request size */
    MPI_Type_size(buftype, &elsize); /* buftype is never MPI_DATATYPE_NULL */
    put_size = elsize;
    if (count != NULL) { /* if count == NULL, this is var1 request */
        for (i=0; i<ndims; i++)
            put_size *= count[i];
    }

    /* Skip zero-length request */
    if (put_size == 0) return NC_NOERR;

    /* Update the largest request size
     * This is used later to determine minimal buffer size for flushing
     */
    if (ncbbp->maxentrysize < put_size)
        ncbbp->maxentrysize = put_size;

    /* Update record dimension size if is record variable */
    if (pncp->vars[varid].recdim >= 0) {
        /* Dim size after this put request */
        if (stride == NULL)
            recsize = start[0] + ((count == NULL) ? 1 : count[0]);
        else if (count == NULL)
            recsize = start[0] + 1;
        else
            recsize = start[0] + (count[0] - 1) * stride[0] + 1;

        if (recsize > ncbbp->recdimsize)
            ncbbp->recdimsize = recsize;
    }

    err = mpitype2logtype(buftype, &itype);
    if (err != NC_NOERR){
#ifdef PNETCDF_DEBUG
        int name_len;
        char type_name[MPI_MAX_OBJECT_NAME];
        MPI_Type_get_name(buftype, type_name, &name_len);
        fprintf(stderr, "Rank: %d, Unrecognized buftype %s\n", ncbbp->rank,
                type_name);
        fflush(stderr);
#endif
        return err;
    }

    /* Prepare log metadata entry header */

    /* Find out the location of data in datalog
     * Which is the current position in data log
     * Datalog descriptor should always points to the end of file
     * Position must be recorded first before writing
     */
    dataoff = (MPI_Offset)ncbbp->datalogsize;

    /* Size of metadata entry
     * Include metadata entry header and variable size additional data
     * (start, count, stride)
     */
    esize = sizeof(NC_bb_metadataentry) + ndims * 3 * SIZEOF_MPI_OFFSET;
    /* Allocate space for metadata entry header */
    buffer = (char*)ncbbio_log_buffer_alloc(&(ncbbp->metadata), esize);
    entryp = (NC_bb_metadataentry*)buffer;
    entryp->esize = esize;  /* Entry size */
    entryp->itype = itype;  /* element data type in internal representation */
    entryp->varid = varid;  /* Variable id */
    entryp->ndims = ndims;  /* Number of dimensions of the variable*/
    /* The size of data in bytes. The size that will be write to data log */
    entryp->data_len = put_size;
    entryp->data_off = dataoff;

    /* Determine the api kind of original call
     * If stride is NULL, we log it as a vara call, otherwise, a vars call
     * Upper layer translates var1 and var to vara  and vars
     */
    if (stride == NULL)
        entryp->api_kind = NC_LOG_API_KIND_VARA;
    else
        entryp->api_kind = NC_LOG_API_KIND_VARS;

    /* Calculate location of start, count, stride in metadata buffer */
    Start = (MPI_Offset*)(buffer + sizeof(NC_bb_metadataentry));
    Count = Start + ndims;
    Stride = Count + ndims;

    /* Fill up start, count, and stride */
    for (i=0; i<ndims; i++) {
        Start[i]  = start[i];
        Count[i]  = (count  == NULL) ? 1 :  count[i];
        Stride[i] = (stride == NULL) ? 0 : stride[i];
    }

    /* Increment number of entry
     * This must be the final step of creating a log entry
     * Increasing num_entries marks the completion of the creation
     */

    /* Increment num_entries in the metadata buffer */
    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;
    headerp->num_entries++;

    /* We only increase datalogsize by amount actually write */
    ncbbp->datalogsize += put_size;

    /* Record data size */
    ncbbio_log_sizearray_append(&(ncbbp->entrydatasize), entryp->data_len);

    /* Record in index
     * Entry address must be relative as metadata buffer can be reallocated
     */
    ncbbio_metaidx_add(ncbbp, (NC_bb_metadataentry*)((char*)entryp -
                       (char*)(ncbbp->metadata.buffer)));

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
#endif

    /*
     * Write data log
     * Note: Metadata record indicate completion, so data must go first
     */
    err = ncbbio_sharedfile_write(ncbbp->datalog_fd, buf, put_size);
    if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
    t3 = MPI_Wtime();
#endif

    /* Seek to the head of metadata file
     * Note: EOF may not be the place for next entry after a flush
     * Note: metadata size will be updated after allocating metadata buffer
     *       space, subtract esize to get the starting offset
     */

    err = ncbbio_sharedfile_seek(ncbbp->metalog_fd,
                                 ncbbp->metadata.nused - esize, SEEK_SET);
    if (err != NC_NOERR) return err;

    /* Write the new entry to metadata log file */
    err = ncbbio_sharedfile_write(ncbbp->metalog_fd, buffer, esize);
    if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
    t4 = MPI_Wtime();
#endif

    /* Update num_entries
     * This marks the completion of the new entry creation
     */
    err = ncbbio_sharedfile_pwrite(ncbbp->metalog_fd, &headerp->num_entries,
                                   SIZEOF_MPI_OFFSET, 56);
    if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
    t5 = MPI_Wtime();
    ncbbp->put_data_wr_time += t3 - t2;
    ncbbp->put_meta_wr_time += t4 - t3;
    ncbbp->put_num_wr_time += t5 - t4;
    ncbbp->total_time += t5 - t1;
    ncbbp->put_time += t5 - t1;

    ncbbp->total_data += put_size;
    ncbbp->total_meta += esize;
#endif

    return NC_NOERR;
}

/*
 * Prepare a n log entry to be write to log
 * IN    ncbbp:  log structure to log this entry
 * IN    varp:   NC_var structure associate to this entry
 * IN    num:      Number of locations
 * IN    start:  start in put_var* call
 * IN    count:  count in put_var* call
 * IN    bur:    buffer of data to write
 * IN    buftype: user buffer's internal data type, an MPI primitive type
 */
int ncbbio_log_put_varn(NC_bb            *ncbbp,
                       int                varid,
                       int                num,
                       MPI_Offset* const *starts,  /* must not be NULL */
                       MPI_Offset* const *counts,  /* may be NULL */
                       void              *buf,
                       MPI_Datatype       buftype)
{
    int err, i, j, ndims, elsize, itype;
    char *buffer;
    PNC *pncp;
    MPI_Offset *start, *count;
    MPI_Offset esize, recsize, put_size, total_put_size;
    MPI_Offset *Start, *Count;
    NC_bb_metadataentry *entryp;
    NC_bb_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    double t1, t2, t3, t4, t5;
    t1 = MPI_Wtime();
#endif

    /* Get PNC */
    err = PNC_check_id(ncbbp->ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* Get ndims */
    ndims = pncp->vars[varid].ndims;

    /* Get element type size */
    MPI_Type_size(buftype, &elsize); /* buftype is never MPI_DATATYPE_NULL */

    /* Convert to ext type */
    err = mpitype2logtype(buftype, &itype);
    if (err != NC_NOERR){
#ifdef PNETCDF_DEBUG
        int name_len;
        char type_name[MPI_MAX_OBJECT_NAME];
        MPI_Type_get_name(buftype, type_name, &name_len);
        fprintf(stderr, "Rank: %d, Unrecognized buftype %s\n", ncbbp->rank,
                type_name);
        fflush(stderr);
#endif
        return err;
    }

    total_put_size = 0;
    for(j = 0; j < num; j++){
        start = (MPI_Offset*)starts[j];
        if (counts != NULL){
            count = (MPI_Offset*)counts[j];
        }
        else{
            count = NULL;
        }

        /* Calcalate put request size */
        put_size = elsize;
        if (count != NULL) { /* if count == NULL, this is var1 request */
            for (i=0; i<ndims; i++)
                put_size *= count[i];
        }
        total_put_size += put_size;

        /* Skip zero-length request */
        if (put_size == 0) continue;

        /* Update record dimension size if is record variable */
        if (pncp->vars[varid].recdim >= 0) {
            /* Dim size after this put request */
            recsize = start[0] + ((count == NULL) ? 1 : count[0]);
            if (recsize > ncbbp->recdimsize)
                ncbbp->recdimsize = recsize;
        }
    }

    /* Prepare log metadata entry header */

    /* Find out the location of data in datalog
    * Which is the current position in data log
    * Datalog descriptor should always points to the end of file
    * Position must be recorded first before writing
    */

    /* Size of metadata entry
    * Include metadata entry header and variable size additional data
    * (start, count, stride)
    */
    if (counts != NULL){
        esize = sizeof(NC_bb_metadataentry) + ndims * 2 * SIZEOF_MPI_OFFSET * num;
    }
    else{
        esize = sizeof(NC_bb_metadataentry) + ndims * SIZEOF_MPI_OFFSET * num;
    }
    /* Allocate space for metadata entry header */
    buffer = (char*)ncbbio_log_buffer_alloc(&(ncbbp->metadata), esize);
    entryp = (NC_bb_metadataentry*)buffer;
    entryp->esize = esize;  /* Entry size */
    entryp->itype = itype;  /* element data type in internal representation */
    entryp->varid = varid;  /* Variable id */
    entryp->ndims = ndims;  /* Number of dimensions of the variable*/
    entryp->api_kind = num; /* Positive number indicate varn */
    /* The size of data in bytes. The size that will be write to data log */
    entryp->data_len = total_put_size;
    /* Find out the location of data in datalog
    * Which is the current position in data log
    * Datalog descriptor should always points to the end of file
    * Position must be recorded first before writing
    */
    entryp->data_off = (MPI_Offset)ncbbp->datalogsize;


    /* Calculate location of start, count, stride in metadata buffer */
    Start = (MPI_Offset*)(buffer + sizeof(NC_bb_metadataentry));
    Count = Start + ndims * num;

    /* Fill up start, count, and stride */
    for(i = 0; i < num; i++) {
        memcpy(Start + i * ndims, starts[i], ndims * SIZEOF_MPI_OFFSET);
    }
    if (counts != NULL) {
        for (i=0; i<num; i++) {
            if (counts[i] != NULL){
                memcpy(Count + i * ndims, counts[i], ndims * SIZEOF_MPI_OFFSET);
            }
            else{
                for(j = 0; j < ndims; j++){
                    Count[i * ndims + j] = 1;
                }
            }
        }
    }

    /* Increment number of entry
    * This must be the final step of creating a log entry
    * Increasing num_entries marks the completion of the creation
    */

    /* Increment num_entries in the metadata buffer */
    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;
    headerp->num_entries++;

    /* We only increase datalogsize by amount actually write */
    ncbbp->datalogsize += total_put_size;

    /* Record data size */
    ncbbio_log_sizearray_append(&(ncbbp->entrydatasize), entryp->data_len);

    /* Record in index
    * Entry address must be relative as metadata buffer can be reallocated
    */
    ncbbio_metaidx_add(ncbbp, (NC_bb_metadataentry*)((char*)entryp -
                    (char*)(ncbbp->metadata.buffer)));

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
#endif

    /* Update the largest request size
    * This is used later to determine minimal buffer size for flushing
    */
    if (ncbbp->maxentrysize < total_put_size)
        ncbbp->maxentrysize = total_put_size;

    /*
     * Write data log
     * Note: Metadata record indicate completion, so data must go first
     */
    err = ncbbio_sharedfile_write(ncbbp->datalog_fd, buf, total_put_size);
    if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
    t3 = MPI_Wtime();
#endif

    /* Seek to the head of metadata file
     * Note: EOF may not be the place for next entry after a flush
     * Note: metadata size will be updated after allocating metadata buffer
     *       space, subtract esize to get the starting offset
     */

    err = ncbbio_sharedfile_seek(ncbbp->metalog_fd,
                                 ncbbp->metadata.nused - esize, SEEK_SET);
    if (err != NC_NOERR) return err;

    /* Write the new entry to metadata log file */
    err = ncbbio_sharedfile_write(ncbbp->metalog_fd, buffer, esize);
    if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
    t4 = MPI_Wtime();
#endif

    /* Update num_entries
     * This marks the completion of the new entry creation
     */
    err = ncbbio_sharedfile_pwrite(ncbbp->metalog_fd, &headerp->num_entries,
                                   SIZEOF_MPI_OFFSET, 56);
    if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
    t5 = MPI_Wtime();
    ncbbp->put_data_wr_time += t3 - t2;
    ncbbp->put_meta_wr_time += t4 - t3;
    ncbbp->put_num_wr_time += t5 - t4;
    ncbbp->total_time += t5 - t1;
    ncbbp->put_time += t5 - t1;

    ncbbp->total_data += total_put_size;
    ncbbp->total_meta += esize;
#endif

    return NC_NOERR;
}
