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
 * Create a new log structure
 * IN      info:    MPI info passed to ncmpi_create/ncmpi_open
 * INOUT   ncdwp:   NC_dw object holding the log structure 
 */
int ncdwio_log_create(NC_dw* ncdwp, MPI_Info info) {
    int i, rank, np, err, flag, masterrank;
    char logbase[NC_LOG_PATH_MAX], basename[NC_LOG_PATH_MAX];
    char *abspath, *fname;
    char *private_path = NULL, *stripe_path = NULL;
    char *logbasep = ".";
    int log_per_node = 0;
#ifdef PNETCDF_PROFILING
    double t1, t2;
#endif
    DIR *logdir;
    ssize_t headersize;
    NC_dw_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    t1 = MPI_Wtime();
#endif

    /* Get rank and number of processes */
    err = MPI_Comm_rank(ncdwp->comm, &rank);
    if (err != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(err, "MPI_Comm_rank");
        DEBUG_RETURN_ERROR(err);
    }
    err = MPI_Comm_size(ncdwp->comm, &np);
    if (err != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(err, "MPI_Comm_rank");
        DEBUG_RETURN_ERROR(err);
    }
    masterrank = rank;
       
    /* Initialize log structure */
    
    /* Determine log file name
     * Log file name is $(bufferdir)$(basename)_$(ncid)_$(rank).{meta/data}
     * filepath is absolute path to the cdf file
     */

    /* Read environment variable for burst buffer path */
    private_path = getenv("DW_JOB_PRIVATE");
    stripe_path = getenv("DW_JOB_STRIPED");
    
    /* Determine log base */
    if (ncdwp->logbase[0] != '\0'){
        logbasep = ncdwp->logbase;
    }
    else if (private_path != NULL){
        logbasep = private_path;
    }
    else if (stripe_path != NULL){
        logbasep = stripe_path;
    }

    /* 
     * Make sure bufferdir exists 
     * NOTE: Assume directory along netcdf file path exists
     */
    logdir = opendir(logbasep);
    if (logdir == NULL) {
        /* Log base does not exist or not accessible */
        DEBUG_RETURN_ERROR(NC_EBAD_FILE);
    }
    closedir(logdir);

    /* Resolve absolute path */    
    abspath = realpath(ncdwp->path, basename);
    if (abspath == NULL){
        /* Can not resolve absolute path */
        DEBUG_RETURN_ERROR(NC_EBAD_FILE);
    }
    abspath = realpath(logbasep, logbase);
    if (abspath == NULL){
        /* Can not resolve absolute path */
        DEBUG_RETURN_ERROR(NC_EBAD_FILE);
    }

    /* Warn if log base not set by user */
    if (rank == 0){
        if (ncdwp->logbase[0] == '\0'){
            printf("Warning: Log directory not set. Using %s.\n", logbase);
            fflush(stdout);
        }
    }

    /* Determine log to process mapping */
    if (rank == 0){
        int j;
        char abs_private_path[NC_LOG_PATH_MAX], abs_stripe_path[NC_LOG_PATH_MAX]; 

        /* Resolve DW_JOB_PRIVATE and DW_JOB_STRIPED into absolute path*/
        memset(abs_private_path, 0, sizeof(abs_private_path));
        memset(abs_stripe_path, 0, sizeof(abs_stripe_path));
        if (private_path != NULL){
            abspath = realpath(private_path, abs_private_path);
            if (abspath == NULL){
                /* Can not resolve absolute path */
                memset(abs_private_path, 0, sizeof(abs_private_path));
            }
        }
        if (stripe_path != NULL){
            abspath = realpath(stripe_path, abs_stripe_path);
            if (abspath == NULL){
                /* Can not resolve absolute path */
                memset(abs_stripe_path, 0, sizeof(abs_stripe_path));
            }
        }

        /* Match against logbase */
        for(i = 0; i < NC_LOG_PATH_MAX; i++){
            if (logbase[i] == '\0' || abs_private_path[i] == '\0'){
                break;
            }
            if (logbase[i] != abs_private_path[i]){
                break;
            }
        }
        for(j = 0; j < NC_LOG_PATH_MAX; j++){
            if (logbase[j] == '\0' || abs_stripe_path[j] == '\0'){
                break;
            }
            if (logbase[j] != abs_stripe_path[j]){
                break;
            }
        }

        /* Whichever has longer matched prefix is considered a match 
         * Use log per node only when striped mode wins
         */
        if (j > i) {
            log_per_node = 1;
        }
        else {
            log_per_node = 0;
        }

        /* Hints can overwrite the default action */
        if (ncdwp->hints & NC_LOG_HINT_LOG_SHARE){
            log_per_node = 1;
        }
    }
    MPI_Bcast(&log_per_node, 1, MPI_INT, 0, ncdwp->comm);
    
    if (log_per_node){
        MPI_Comm_split_type(ncdwp->comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &(ncdwp->logcomm));
        MPI_Bcast(&masterrank, 1, MPI_INT, 0, ncdwp->logcomm); 
    }
    else{
        ncdwp->logcomm = MPI_COMM_SELF;
        masterrank = rank;
    }

    /* Extract file anme 
     * Search for first / charactor from the tail
     * Absolute path should always contains one '/', return error otherwise
     * We return the string including '/' for convenience
     */
    fname = strrchr(basename, '/');
    if (fname == NULL){
        DEBUG_RETURN_ERROR(NC_EBAD_FILE);  
    }
    
    /* Log file path may also contain non-existing directory
     * We need to create them before we can search for usable id
     * As log file name hasn't been determined, we need to use a dummy one here
     */
    sprintf(ncdwp->metalogpath, "%s%s_%d_%d.meta", logbase, fname,
            ncdwp->ncid, masterrank);
    sprintf(ncdwp->datalogpath, "%s%s_%d_%d.data", logbase, fname,
            ncdwp->ncid, masterrank);
     
    /* Initialize metadata buffer */
    err = ncdwio_log_buffer_init(&(ncdwp->metadata));
    if (err != NC_NOERR){
        return err;
    }
       
    /* Initialize metadata entry array */
    err = ncdwio_log_sizearray_init(&(ncdwp->entrydatasize));
    if (err != NC_NOERR){
        return err;
    }

    /* Set log file descriptor to NULL */

#ifdef PNETCDF_PROFILING
    /* Performance counters */
    ncdwp->total_time = 0;
    ncdwp->create_time = 0;
    ncdwp->enddef_time = 0;
    ncdwp->put_time = 0;
    ncdwp->flush_time = 0;
    ncdwp->close_time = 0;
    ncdwp->flush_replay_time = 0;
    ncdwp->flush_data_rd_time = 0;
    ncdwp->flush_put_time = 0;
    ncdwp->flush_wait_time = 0;
    ncdwp->put_data_wr_time = 0;
    ncdwp->put_meta_wr_time = 0;
    ncdwp->put_num_wr_time = 0;
    ncdwp->max_buffer = 0;
#endif

    /* Misc */
    ncdwp->rank = rank;
    ncdwp->np = np;
    ncdwp->maxentrysize = 0;

    /* Initialize metadata header */
    
    /*
     * Allocate space for metadata header
     * Header consists of a fixed size info and variable size basename
     */
    headersize = sizeof(NC_dw_metadataheader) + strlen(basename);
    if (headersize % 4 != 0){
        headersize += 4 - (headersize % 4);
    }
    headerp = (NC_dw_metadataheader*)ncdwio_log_buffer_alloc(&(ncdwp->metadata),
                                                             headersize);
    memset(headerp, 0, headersize);

    /* Fill up the metadata log header */
    memcpy(headerp->magic, NC_LOG_MAGIC, sizeof(headerp->magic));
    memcpy(headerp->format, NC_LOG_FORMAT_CDF_MAGIC, sizeof(headerp->format));
    strncpy((char*)headerp->basename, basename,
            headersize - sizeof(NC_dw_metadataheader) + 1);
    headerp->rank_id = rank;   /* Rank */
    headerp->num_ranks = np;   /* Number of processes */
    /* Without convertion before logging, data in native representation */
    headerp->is_external = 0;    
    /* Determine endianess */
#ifdef WORDS_BIGENDIAN 
    headerp->big_endian = NC_LOG_TRUE;
#else 
    headerp->big_endian = NC_LOG_FALSE;
#endif
    headerp->num_entries = 0;    /* The log is empty */
    /* Highest dimension among all variables */
    headerp->max_ndims = 0;    
    /* Location of the first entry */
    headerp->entry_begin = ncdwp->metadata.nused;  
    headerp->basenamelen = strlen(basename);

    /* Create log files */
    flag = O_RDWR | O_CREAT;
    if (!(ncdwp->hints & NC_LOG_HINT_LOG_OVERWRITE)) {
        flag |= O_EXCL;
    }
    //ncdwp->datalog_fd = ncdwp->metalog_fd = -1;
    err = ncdwio_sharedfile_open(ncdwp->logcomm, ncdwp->metalogpath, flag,
                           MPI_INFO_NULL, &(ncdwp->metalog_fd));
    if (err != NC_NOERR) {
        return err;
    }
    err = ncdwio_bufferedfile_open(ncdwp->logcomm, ncdwp->datalogpath, flag,
                           MPI_INFO_NULL, &(ncdwp->datalog_fd));
    if (err != NC_NOERR) {
        return err;
    }
   
    /* Write metadata header to file
     * Write from the memory buffer to file
     */
    err = ncdwio_sharedfile_write(ncdwp->metalog_fd, headerp, headersize);
    if (err != NC_NOERR){
        return err;
    }

    /* Write data header to file 
     * Data header consists of a fixed sized string PnetCDF0
     */
    err = ncdwio_bufferedfile_write(ncdwp->datalog_fd, "PnetCDF0", 8);
    if (err != NC_NOERR){
        return err;
    }

    ncdwp->datalogsize = 8;

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
    ncdwp->total_time += t2 - t1;
    ncdwp->create_time += t2 - t1;

    ncdwp->total_meta += headersize;
    ncdwp->total_data += 8;
#endif

    return NC_NOERR;
}

/*
 * Update information used by bb layer on enddef
 * IN    ncdwp:    NC_dw structure
 */
int ncdwio_log_enddef(NC_dw *ncdwp){   
    int err;
#ifdef PNETCDF_PROFILING
    double t1, t2; 
#endif
    NC_dw_metadataheader *headerp;
    
#ifdef PNETCDF_PROFILING
    t1 = MPI_Wtime();
#endif

    headerp = (NC_dw_metadataheader*)ncdwp->metadata.buffer;
    
    /* 
     * Update the header if max ndim increased 
     * For now, max_ndims is the only field that can change
     */
    if (ncdwp->max_ndims > headerp->max_ndims){
        headerp->max_ndims = ncdwp->max_ndims;

        /* Overwrite maxndims
         * This marks the completion of the record
         */
        err = ncdwio_sharedfile_pwrite(ncdwp->metalog_fd, &headerp->max_ndims,
                                SIZEOF_MPI_OFFSET, sizeof(NC_dw_metadataheader) - 
                               sizeof(headerp->basename) - 
                               sizeof(headerp->basenamelen) - 
                               sizeof(headerp->num_entries) - 
                               sizeof(headerp->max_ndims));
        if (err != NC_NOERR){
            return err;
        }
    }

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
    ncdwp->total_time += t2 - t1;
    ncdwp->enddef_time += t2 - t1;
#endif

    return NC_NOERR;
}

/*
 * Flush the log to CDF file and clean up the log structure 
 * Used by ncmpi_close()
 * IN    ncdwp:    log structure
 */
int ncdwio_log_close(NC_dw *ncdwp) {
    int err;
#ifdef PNETCDF_PROFILING
    double t1, t2; 
    unsigned long long total_data;
    unsigned long long total_meta;
    unsigned long long buffer_size;
    double total_time;
    double create_time;
    double enddef_time;
    double put_time;
    double flush_time;
    double close_time;
    double flush_replay_time;
    double flush_data_rd_time;
    double flush_put_time;
    double flush_wait_time;
    double put_data_wr_time;
    double put_meta_wr_time;
    double put_num_wr_time;
#endif
    NC_dw_metadataheader* headerp;

#ifdef PNETCDF_PROFILING
    t1 = MPI_Wtime();
#endif

    headerp = (NC_dw_metadataheader*)ncdwp->metadata.buffer;

    /* If log file is created, flush the log */
    if (ncdwp->metalog_fd >= 0){
        /* Commit to CDF file */
        if (headerp->num_entries > 0){
            log_flush(ncdwp);
        }

        /* Close log file */
        err = ncdwio_sharedfile_close(ncdwp->metalog_fd);
        if (err != NC_NOERR){
            return err;
        }
        err = ncdwio_bufferedfile_close(ncdwp->datalog_fd);
        if (err != NC_NOERR){
            return err;
        }

        /* Delete log files if delete flag is set */
        if (ncdwp->hints & NC_LOG_HINT_DEL_ON_CLOSE){
            unlink(ncdwp->datalogpath);
            unlink(ncdwp->metalogpath);
        }
    }

    /* Free meta data buffer and metadata offset list*/
    ncdwio_log_buffer_free(&(ncdwp->metadata));
    ncdwio_log_sizearray_free(&(ncdwp->entrydatasize));

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
    ncdwp->total_time += t2 - t1;
    ncdwp->close_time += t2 - t1;

#ifdef PNETCDF_DEBUG
    /* Print accounting info in debug build */
    MPI_Reduce(&(ncdwp->total_time), &total_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncdwp->comm);
    MPI_Reduce(&(ncdwp->create_time), &create_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncdwp->comm);
    MPI_Reduce(&(ncdwp->enddef_time), &enddef_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncdwp->comm);
    MPI_Reduce(&(ncdwp->put_time), &put_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncdwp->comm);
    MPI_Reduce(&(ncdwp->flush_time), &flush_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncdwp->comm);
    MPI_Reduce(&(ncdwp->close_time), &close_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncdwp->comm);
    MPI_Reduce(&(ncdwp->flush_replay_time), &flush_replay_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->flush_data_rd_time), &flush_data_rd_time, 1,
                MPI_DOUBLE, MPI_MAX, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->flush_put_time), &flush_put_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->flush_wait_time), &flush_wait_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->put_data_wr_time), &put_data_wr_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->put_meta_wr_time), &put_meta_wr_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->put_num_wr_time), &put_num_wr_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->total_meta), &total_meta, 1, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->total_data), &total_data, 1, MPI_UNSIGNED_LONG_LONG, 
                MPI_SUM, 0, ncdwp->comm);
    MPI_Reduce(&(ncdwp->flushbuffersize), &buffer_size, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, ncdwp->comm);
    
    if (ncdwp->rank == 0){ 
        printf("==========================================================\n");
        printf("File: %s\n", ncdwp->path);
        printf("Data writen to variable: %llu\n", total_data);
        printf("Metadata generated: %llu\n", total_meta);
        printf("Flush buffer size: %llu\n", buffer_size);
        printf("Time in log: %lf\n", total_time);
        printf("\tTime in log_create: %lf\n", create_time);
        printf("\tTime in log_enddef: %lf\n", enddef_time);
        printf("\tTime in log_put: %lf\n", put_time);
        printf("\t\tTime writing data log: %lf\n", put_data_wr_time);
        printf("\t\tTime writing metadata log: %lf\n", put_meta_wr_time);
        printf("\t\tTime updating numrecs: %lf\n", put_num_wr_time);
        printf("\tTime in log_flush: %lf\n", flush_time);
        printf("\tTime in log_close: %lf\n", close_time);
        printf("\tTime replaying the log: %lf\n", flush_replay_time);
        printf("\t\tTime reading data log: %lf\n", flush_data_rd_time);
        printf("\t\tTime calling iput: %lf\n", flush_put_time);
        printf("\t\tTime calling wait: %lf\n", flush_wait_time);
        printf("==========================================================\n");
    }
#endif
#endif

    return NC_NOERR;
}

/*
 * Commit the log into cdf file and delete the log
 * User can call this to force a commit without closing
 * It work by flush and re-initialize the log structure
 * IN    ncdwp:    log structure
 */
int ncdwio_log_flush(NC_dw* ncdwp) {
    int err, status = NC_NOERR;
#ifdef PNETCDF_PROFILING
    double t1, t2; 
#endif
    //NC_req *putlist;
    NC_dw_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    t1 = MPI_Wtime();
#endif

    headerp = (NC_dw_metadataheader*)ncdwp->metadata.buffer;

    /* Nothing to replay if nothing have been written */
    if (headerp->num_entries == 0){
        return NC_NOERR;
    }

    /* Replay log file */
    err = log_flush(ncdwp);
    if (err != NC_NOERR) {
        if (status == NC_NOERR){
            DEBUG_ASSIGN_ERROR(status, err);
        }
    }
 
    /* Reset log status */
    
    /* Set num_entries to 0 */
    headerp->num_entries = 0;

    /* Overwrite num_entries
     * This marks the completion of flush
     */
    err = ncdwio_sharedfile_pwrite(ncdwp->metalog_fd, &headerp->num_entries,
                            SIZEOF_MPI_OFFSET, 56);
    if (err != NC_NOERR){
        return err;
    }

    /* Reset metadata buffer and entry array status */
    ncdwp->metadata.nused = headerp->entry_begin;
    ncdwp->entrydatasize.nused = 0;
    ncdwp->metaidx.nused = 0;

    /* Rewind data log file descriptors and reset the size */
    err = ncdwio_bufferedfile_seek(ncdwp->datalog_fd, 8, SEEK_SET);
    if (err != NC_NOERR){
        return err;
    }

    ncdwp->datalogsize = 8;

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
    ncdwp->total_time += t2 - t1;
    ncdwp->flush_time += t2 - t1;
#endif

    return status;
}
