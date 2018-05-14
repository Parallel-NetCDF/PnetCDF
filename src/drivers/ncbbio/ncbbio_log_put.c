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
#include <ncbbio_driver.h>

/*----< check_EINVALCOORDS() >-----------------------------------------------*/
static
int check_EINVALCOORDS(MPI_Offset start,
                       MPI_Offset count,
                       MPI_Offset shape)
{
#ifdef RELAX_COORD_BOUND
    if (start < 0 || start > shape)
        DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
    if (start == shape && count > 0)
        DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
#else
    if (start < 0 || start >= shape)
        DEBUG_RETURN_ERROR(NC_EINVALCOORDS)
#endif
    return NC_NOERR;
}

typedef enum {
    API_GET,
    API_PUT,
    API_IGET,
    API_IPUT,
    API_BPUT
} IO_type;

/*----< check_EEDGE() >------------------------------------------------------*/
static
int check_EEDGE(const MPI_Offset *start,
                const MPI_Offset *count,
                const MPI_Offset *stride,
                const MPI_Offset *shape)
{
    if (*count > *shape || *start + *count > *shape)
        DEBUG_RETURN_ERROR(NC_EEDGE);
    if (stride == NULL) { /* vars APIs but stride is NULL */
        if (*count > *shape || *start + *count > *shape)
            DEBUG_RETURN_ERROR(NC_EEDGE)
    }
    else { /* for vars/varm APIs */
        if (*count > 0 && *start + (*count - 1) * (*stride) >= *shape)
            DEBUG_RETURN_ERROR(NC_EEDGE)
    }
    return NC_NOERR;
}

/*----< sanity_check() >-----------------------------------------------------*/
static
int sanity_check(PNC          *pncp,
                 int           varid,
                 IO_type       io,       /* get/put/iget/iput/bput */
                 MPI_Datatype  itype,    /* internal data type */
                 int           isColl)   /* collective or indepdnent API */
{
    /* check file write permission for put APIs */
    if (io == API_PUT || io == API_IPUT || io == API_BPUT)
        if (pncp->flag & NC_MODE_RDONLY) DEBUG_RETURN_ERROR(NC_EPERM)

    /* blocking get/put APIs must be called in data mode */
    if (io == API_PUT || io == API_GET)
        if (pncp->flag & NC_MODE_DEF) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* for blocking APIs, check if in collective or independent mode */
    if (io == API_PUT || io == API_GET) {
        if (isColl) { /* check if file is currently in collective data mode */
            if (pncp->flag & NC_MODE_INDEP) DEBUG_RETURN_ERROR(NC_EINDEP)
        }
        else { /* check if file is currently in independent data mode */
            if (!(pncp->flag & NC_MODE_INDEP)) DEBUG_RETURN_ERROR(NC_ENOTINDEP)
        }
    }

    /* variable NC_GLOBAL is illegal in get/put APIs */
    if (varid == NC_GLOBAL) DEBUG_RETURN_ERROR(NC_EGLOBAL)

    /* check whether variable ID is valid */
    if (varid < 0 || varid >= pncp->nvars) DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* MPI_DATATYPE_NULL in this case represent a flexible API */
    if (itype == MPI_DATATYPE_NULL) return NC_NOERR;

    /* check itype against xtype for NC_ECHAR */
    if (itype == MPI_CHAR) {
        if (pncp->vars[varid].xtype != NC_CHAR) DEBUG_RETURN_ERROR(NC_ECHAR)
    }
    else {
        if (pncp->vars[varid].xtype == NC_CHAR) DEBUG_RETURN_ERROR(NC_ECHAR)
    }
    return NC_NOERR;
}

/*----< check_start_count_stride() >-----------------------------------------*/
static
int check_start_count_stride(PNC              *pncp,
                             int               varid,
                             int               isRead,
                             NC_api            api_kind, /* var1/vara/vars */
                             const MPI_Offset *start,
                             const MPI_Offset *count,
                             const MPI_Offset *stride)
{
    /* only var1, vara, vars, and varm APIs will reach here */
    int i, err, ndims, firstDim;
    MPI_Offset *shape=NULL;

    shape = pncp->vars[varid].shape;
    /* if record variable, obtain the current size of record dimension */
    if (pncp->vars[varid].recdim >= 0) {
        err = pncp->driver->inq_dim(pncp->ncp, pncp->vars[varid].recdim, NULL,
                                    &shape[0]);
        if (err != NC_NOERR) return err;
    }

    /* Check NC_EINVALCOORDS error for argument start[]
     * for API var1/vara/vars/varm, start cannot be NULL, except for scalars
     * and negative start[] is illegal */
    if (start == NULL || start[0] < 0) DEBUG_RETURN_ERROR(NC_EINVALCOORDS)

    firstDim = 0;
    /* check NC_EINVALCOORDS for record dimension */
    if (pncp->vars[varid].recdim >= 0) {
        if (pncp->format < NC_FORMAT_CDF5 && start[0] > NC_MAX_UINT)
            DEBUG_RETURN_ERROR(NC_EINVALCOORDS) /* CDF-1 and 2 */

        /* for record variable, [0] is the NC_UNLIMITED dimension */
        /* read cannot go beyond current numrecs */
        if (isRead) {
            MPI_Offset len = (count == NULL) ? 1 : count[0];
            err = check_EINVALCOORDS(start[0], len, shape[0]);
            if (err != NC_NOERR) return err;
        }
        firstDim = 1; /* done for checking the record dimension */
    }

    /* continue to check NC_EINVALCOORDS for the rest dimensions */
    ndims = pncp->vars[varid].ndims;
    for (i=firstDim; i<ndims; i++) {
        MPI_Offset len = (count == NULL) ? 1 : count[i];
        err = check_EINVALCOORDS(start[i], len, shape[i]);
        if (err != NC_NOERR) return err;
    }

    /* check NC_EEDGE error for argument count[] */

    if (count == NULL) {
        if (api_kind == API_VARA || api_kind == API_VARS ||
            api_kind == API_VARM)
            /* vara/vars/varm, count cannot be NULL */
            DEBUG_RETURN_ERROR(NC_EEDGE)
    }
    else {
        firstDim = 0;
        /* check record dimension */
        if (pncp->vars[varid].recdim >= 0) {
            if (count[0] < 0)  /* no negative count[] */
                DEBUG_RETURN_ERROR(NC_ENEGATIVECNT)

            /* for record variable, [0] is the NC_UNLIMITED dimension */
            /* read cannot go beyond current numrecs */
            if (isRead) {
                err = check_EEDGE(start, count, stride, shape);
                if (err != NC_NOERR) return err;
            }
            firstDim = 1; /* skip checking the record dimension */
        }

        /* continue to check NC_EEDGE for the rest dimensions */
        for (i=firstDim; i<ndims; i++) {
            if (shape[i] < 0) DEBUG_RETURN_ERROR(NC_EEDGE)
            if (count[i] < 0) /* no negative count[] */
                DEBUG_RETURN_ERROR(NC_ENEGATIVECNT)
            if (stride == NULL)
                err = check_EEDGE(start+i, count+i, NULL, shape+i);
            else
                err = check_EEDGE(start+i, count+i, stride+i, shape+i);
            if (err != NC_NOERR) return err;
        }

        /* Check NC_ESTRIDE for non-positive values. We did not check
         * stride[i] >= shape[i], as it is caught as NC_EEDGE error above */
        if (stride != NULL) {
            for (i=0; i<ndims; i++) {
                if (stride[i] <= 0) DEBUG_RETURN_ERROR(NC_ESTRIDE)
            }
        }
    }
    return NC_NOERR;
}

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
    int *dimids;
    char *buffer;
    nc_type xtype;
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

    /* In case of collective I/O without safe mode, driver may still be called even argument are invalid
     * We must redo senity check to make sure everything in the log is valid
     */ 
    if (!ncbbp->isindep && !(pncp->flag & NC_MODE_SAFE)) {
        NC_api api_kind;

        if (stride == NULL){
            api_kind = API_VARA;
        }
        else{
            api_kind = API_VARS;
        }

        /* Sanity check */
        err = sanity_check(pncp, varid, API_PUT, buftype, 1);
        if (err != NC_NOERR){
            return err;
        }

        /* not-scalar variable checks start, count, stride */
        if (dim > 0){
            err = check_start_count_stride(pncp, varid, 0, api_kind,
                                        start, count, stride);
            if (err != NC_NOERR){
                return err;
            }
        }
    }

    /* Calculate data size */
    
    //xtype = pncp->vars[varid].xtype;
    /*
    err = ncbbp->ncmpio_driver->inq_var(ncbbp->ncp, varid, NULL, &xtype, &dim,
                                        NULL, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR){
        return err;
    }
    */

    /* Count should always be avaiable for non scalar var 
     * Upper layer may not fill up count if an error occurs
     */
    /*if (dim > 0){
    /    if (count == NULL){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
    }*/

    /* We need to check for ECHAR for scalar variables, they are not filtered by upper layer
        * Ncmpio driver will abort when it sees ECHAR, it should never be in the log
        */
    /*if (buftype == MPI_CHAR) {   
        if (xtype != NC_CHAR){
            DEBUG_RETURN_ERROR(NC_ECHAR);
        }
    }
    else{
        if (xtype == NC_CHAR){
            DEBUG_RETURN_ERROR(NC_ECHAR);
        }
    }*/


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

    /* Record largest entry size */
    if (ncbbp->maxentrysize < size){
        ncbbp->maxentrysize = size;
    }

    /* Update record dimension size if is record variable */
    /* Get dimids */
    dimids = NCI_Malloc(SIZEOF_INT * dim);
    err = ncbbp->ncmpio_driver->inq_var(ncbbp->ncp, varid, NULL, NULL, NULL,
                                        dimids, NULL, NULL, NULL, NULL);
    if (err != NC_NOERR){
        return err;
    }
    /* Update recdimsize if first dim is unlimited */
    //if (dim > 0 && dimids[0] == ncbbp->recdimid) {
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
        fprintf(stderr, "Rank: %d, Unrecognized type: %d\n", ncbbp->rank, buftype); fflush(stderr);
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
    // Record in index
    // Entry address must be relative as metadata buffer can be reallocated
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
    err = ncbbio_bufferedfile_write(ncbbp->datalog_fd, buf, size);
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

    err = ncbbio_sharedfile_seek(ncbbp->metalog_fd, ncbbp->metadata.nused - esize,
                           SEEK_SET);
    if (err != NC_NOERR){
        return err;
    }

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
