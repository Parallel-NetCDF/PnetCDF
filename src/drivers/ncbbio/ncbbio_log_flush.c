/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* "$Id$" */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <pnc_debug.h>
#include <common.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <ncbbio_driver.h>

/* Convert from log type to MPI type used by pnetcdf library
 * Log spec has different enum of types than MPI
 */
static int logtype2mpitype(int type, MPI_Datatype *buftype){
    /* Convert from log type to MPI type used by pnetcdf library
     * Log spec has different enum of types than MPI
     */
    if (type == NC_LOG_TYPE_TEXT) {
        *buftype = MPI_CHAR;
    }
    else if (type == NC_LOG_TYPE_SCHAR) {
        *buftype = MPI_SIGNED_CHAR;
    }
    else if (type == NC_LOG_TYPE_SHORT) {
        *buftype = MPI_SHORT;
    }
    else if (type == NC_LOG_TYPE_INT) {
        *buftype = MPI_INT;
    }
    else if (type == NC_LOG_TYPE_LONG) {
        *buftype = MPI_LONG;
    }
    else if (type == NC_LOG_TYPE_FLOAT) {
        *buftype = MPI_FLOAT;
    }
    else if (type == NC_LOG_TYPE_DOUBLE) {
        *buftype = MPI_DOUBLE;
    }
    else if (type == NC_LOG_TYPE_UCHAR) {
        *buftype = MPI_UNSIGNED_CHAR;
    }
    else if (type == NC_LOG_TYPE_USHORT) {
        *buftype = MPI_UNSIGNED_SHORT;
    }
    else if (type == NC_LOG_TYPE_UINT) {
        *buftype = MPI_UNSIGNED;
    }
    else if (type == NC_LOG_TYPE_LONGLONG) {
        *buftype = MPI_LONG_LONG_INT;
    }
    else if (type == NC_LOG_TYPE_ULONGLONG) {
        *buftype = MPI_UNSIGNED_LONG_LONG;
    }
    else {
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }

    return NC_NOERR;
}

/*
 * Commit log file into CDF file
 * Meta data is stored in memory, metalog is only used for restoration after abnormal shutdown
 * IN    ncbbp:    log structure
 */
int ncbbio_log_flush_core(NC_bb *ncbbp) {
    int i, j, lb, ub, err, status = NC_NOERR;
    int *reqids, *stats;
    //int ready, ready_all = 0;
    int nrounds, nrounds_all;
    int len = 0;
    size_t databufferused, databuffersize, dataread;
    NC_bb_metadataentry *entryp;
    MPI_Offset *start, *count, *stride;
    MPI_Offset **starts=NULL, **counts=NULL;
    MPI_Datatype buftype;
    char *databuffer, *databufferoff;
    NC_bb_metadataheader *headerp;
    NC_bb_metadataptr *ip;
#ifdef PNETCDF_PROFILING
    double t1, t2, t3, t4;

    t1 = MPI_Wtime();
#endif

    /* Read datalog in to memory */
    /*
     * Prepare data buffer
     * We determine the data buffer size according to:
     * hints, size of data log, the largest size of single record
     * 0 in hint means no limit
     * (Buffer size) = max((largest size of single record), min((size of data log), (size specified in hint)))
     */
    databuffersize = ncbbp->datalogsize;
    if (ncbbp->flushbuffersize > 0 && databuffersize > (size_t)(ncbbp->flushbuffersize)){
        databuffersize = ncbbp->flushbuffersize;
    }
    /* Without enabling large_req, we can not post requests larger than 2GiB */
#ifndef ENABLE_LARGE_SINGLE_REQ
    if (databuffersize > 2147483647){
        databuffersize = 2147483647;
    }
#endif
    /* We assume user will not issue single request larger than 2GiB wwithout enabling large_req */
    if (databuffersize < (size_t)(ncbbp->maxentrysize)){
        databuffersize = ncbbp->maxentrysize;
    }

#ifdef PNETCDF_PROFILING
    if ((size_t)(ncbbp->max_buffer) < databuffersize){
        ncbbp->max_buffer = databuffersize;
    }
#endif

    /* Sync rounds */
    databufferused = 0;
    dataread = 0;
    nrounds = 0;
    err = ncbbio_sharedfile_seek(ncbbp->datalog_fd, 8, SEEK_SET);
    if (err != NC_NOERR){
        return err;
    }
    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;
    entryp = (NC_bb_metadataentry*)(((char*)ncbbp->metadata.buffer) + headerp->entry_begin);
    for (i = 0; i < ncbbp->metaidx.nused; i++){
        if (ncbbp->metaidx.entries[i].valid){
            if (ncbbp->entrydatasize.values[i] + databufferused > databuffersize){
                nrounds++;
                databufferused = ncbbp->entrydatasize.values[i];
            }
            else{
                databufferused += ncbbp->entrydatasize.values[i];
            }
        }
    }
    nrounds++;
    if (!fIsSet(ncbbp->flag, NC_MODE_INDEP)){
        MPI_Allreduce(&nrounds, &nrounds_all, 1, MPI_INT, MPI_MAX, ncbbp->comm);
    }
    else{
        nrounds_all = nrounds;
    }

    /* Allocate buffer */
    databuffer = (char*)NCI_Malloc(databuffersize);
    if(databuffer == NULL){
        DEBUG_RETURN_ERROR(NC_ENOMEM);
    }

    /* Seek to the start position of first data record */
    err = ncbbio_sharedfile_seek(ncbbp->datalog_fd, 8, SEEK_SET);
    if (err != NC_NOERR){
        return err;
    }


    reqids = (int*)NCI_Malloc(ncbbp->entrydatasize.nused * SIZEOF_INT);
    stats = (int*)NCI_Malloc(ncbbp->entrydatasize.nused * SIZEOF_INT);

    /*
     * Iterate through meta log entries
     */
    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;
    entryp = (NC_bb_metadataentry*)(((char*)ncbbp->metadata.buffer) + headerp->entry_begin);
    for (lb = 0; lb < ncbbp->metaidx.nused;){
        /* Initialize buffer status */
        databufferused = 0;
        dataread = 0;

        for (ub = lb; ub < ncbbp->metaidx.nused; ub++) {
            if (ncbbp->metaidx.entries[ub].valid){
                if(ncbbp->entrydatasize.values[ub] + databufferused > databuffersize) {
                    break;  // Buffer full
                }
                else{
                    databufferused += ncbbp->entrydatasize.values[ub]; // Record size of entry
                }
            }
            else{
                // We encounter a canceled record
                // Read unread data into data buffer and jump through the gap
                /*
                 * Read data to buffer
                 * We read only what needed by pending requests
                 */
                if (dataread < databufferused){
#ifdef PNETCDF_PROFILING
                    t2 = MPI_Wtime();
#endif
                    err = ncbbio_sharedfile_read(ncbbp->datalog_fd, databuffer + dataread, databufferused - dataread);
                    if (err != NC_NOERR){
                        return err;
                    }
#ifdef PNETCDF_PROFILING
                    t3 = MPI_Wtime();
                    ncbbp->flush_data_rd_time += t3 - t2;
#endif
                    dataread = databufferused;
                }

                // Skip canceled entry
                err = ncbbio_sharedfile_seek(ncbbp->datalog_fd, ncbbp->entrydatasize.values[ub], SEEK_CUR);
                if (err != NC_NOERR){
                    return err;
                }
            }
        }

        /*
         * Read data to buffer
         * We read only what needed by pending requests
         */
        if (dataread < databufferused){
#ifdef PNETCDF_PROFILING
            t2 = MPI_Wtime();
#endif
            err = ncbbio_sharedfile_read(ncbbp->datalog_fd, databuffer + dataread, databufferused - dataread);
            if (err != NC_NOERR){
                return err;
            }
#ifdef PNETCDF_PROFILING
            t3 = MPI_Wtime();
            ncbbp->flush_data_rd_time += t3 - t2;
#endif
            dataread = databufferused;
        }

        // Pointer points to the data of current entry
        databufferoff = databuffer;

        j = 0;
        for(i = lb; i < ub; i++){
            ip = ncbbp->metaidx.entries + i;

            if (ip->valid) {
                if (entryp->api_kind < 0){ // VARA or VARS
                    /* start, count, stride */
                    start = (MPI_Offset*)(entryp + 1);
                    count = start + entryp->ndims;
                    stride = count + entryp->ndims;

                    // Convert from log type to MPI type
                    err = logtype2mpitype(entryp->itype, &buftype);
                    if (err != NC_NOERR){
                        return err;
                    }

                    /* Determine API_Kind */
                    if (entryp->api_kind == NC_LOG_API_KIND_VARA){
                        stride = NULL;
                    }

#ifdef PNETCDF_PROFILING
                    t2 = MPI_Wtime();
#endif

                    /* Replay event with non-blocking call */
                    err = ncbbp->ncmpio_driver->iput_var(ncbbp->ncp, entryp->varid, start, count, stride, NULL, (void*)(databufferoff), -1, buftype, reqids + j, NC_REQ_WR | NC_REQ_NBI | NC_REQ_HL);
                    if (status == NC_NOERR) {
                        status = err;
                    }

#ifdef PNETCDF_PROFILING
                    t3 = MPI_Wtime();
                    ncbbp->flush_put_time += t3 - t2;
#endif
                }
                else {  // VARN
                    int k;
                    int num;

                    num = entryp->api_kind;

                    /* Allocate starts and counts array
                     * Try to reuse allocated array if long enough
                     * Don't realloc because old data is not needed
                     */
                    if(len < num){
                        if (len > 0) NCI_Free(starts);
                        starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * num * 2);
                        counts = starts + num;
                        len = num;
                    }
                    start = (MPI_Offset*)(entryp + 1);
                    count = start + entryp->ndims * num;
                    for(k = 0; k < num; k++){
                        starts[k] = start + entryp->ndims * k;
                        counts[k] = count + entryp->ndims * k;
                    }

                    // Convert from log type to MPI type
                    err = logtype2mpitype(entryp->itype, &buftype);
                    if (err != NC_NOERR){
                        return err;
                    }

#ifdef PNETCDF_PROFILING
                    t2 = MPI_Wtime();
#endif

                    /* Replay event with non-blocking call */
                    if ((size_t)(entryp->esize) == sizeof(NC_bb_metadataentry) + entryp->ndims * SIZEOF_MPI_OFFSET * num){
                        err = ncbbp->ncmpio_driver->iput_varn(ncbbp->ncp, entryp->varid, num, starts, NULL, (void*)(databufferoff), -1, buftype, reqids + j, NC_REQ_WR | NC_REQ_NBI | NC_REQ_HL);
                    }
                    else{
                        err = ncbbp->ncmpio_driver->iput_varn(ncbbp->ncp, entryp->varid, num, starts, counts, (void*)(databufferoff), -1, buftype, reqids + j, NC_REQ_WR | NC_REQ_NBI | NC_REQ_HL);
                    }
                    if (status == NC_NOERR) {
                        status = err;
                    }

#ifdef PNETCDF_PROFILING
                    t3 = MPI_Wtime();
                    ncbbp->flush_put_time += t3 - t2;
#endif
                }

                // Move to next data location
                databufferoff += entryp->data_len;
                j++;
            }

            /* Move to next position */
            entryp = (NC_bb_metadataentry*)(((char*)entryp) + entryp->esize);
        }

#ifdef PNETCDF_PROFILING
        t2 = MPI_Wtime();
#endif
        /*
         * Wait must be called first or previous data will be corrupted
         */
        if (fIsSet(ncbbp->flag, NC_MODE_INDEP)){
            err = ncbbp->ncmpio_driver->wait(ncbbp->ncp, j, reqids, stats, NC_REQ_INDEP);
        }
        else{
            err = ncbbp->ncmpio_driver->wait(ncbbp->ncp, j, reqids, stats, NC_REQ_COLL);
        }
        if (status == NC_NOERR) {
            status = err;
        }

#ifdef PNETCDF_PROFILING
        t3 = MPI_Wtime();
        ncbbp->flush_wait_time += t3 - t2;
#endif

        // Fill up the status for nonblocking request
        for(i = lb; i < ub; i++){
            ip = ncbbp->metaidx.entries + i;
            j = 0;
            if (ip->valid) {
                if (ip->reqid >= 0){
                    ncbbp->putlist.reqs[ip->reqid].status = stats[j];
                    ncbbp->putlist.reqs[ip->reqid].ready = 1;
                }
                j++;
            }
        }

        /* Update batch status */
        databufferused = 0;

        // Mark as complete
        lb = ub;

        /*
         * In case of collective flush, we sync our status with other processes
         */
        nrounds_all--;
    }

    /*
     * In case of collective flush, we must continue to call wait until every process is ready
     */
    while(nrounds_all--){
        err = ncbbp->ncmpio_driver->wait(ncbbp->ncp, 0, NULL, NULL, NC_REQ_COLL);
        if (status == NC_NOERR) {
            status = err;
        }
    }


    /* Free the data buffer */
    NCI_Free(databuffer);
    NCI_Free(reqids);
    NCI_Free(stats);
    if (starts != NULL)
        NCI_Free(starts);

#ifdef PNETCDF_PROFILING
    t4 = MPI_Wtime();
    ncbbp->flush_replay_time += t4 - t1;
#endif

    return status;
}


