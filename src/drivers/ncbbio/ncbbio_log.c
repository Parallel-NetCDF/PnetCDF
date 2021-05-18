/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <limits.h> /* realpath() */
#include <stdlib.h> /* realpath() */
#include <string.h>
#include <sys/types.h> /* opendir() closedir() */
#include <dirent.h>    /* opendir() closedir() */
#include <unistd.h>
#include <fcntl.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncbbio_driver.h>

/*
 * Create a new log structure
 * IN      info:    MPI info passed to ncmpi_create/ncmpi_open
 * INOUT   ncbbp:   NC_bb object holding the log structure
 */
int ncbbio_log_create(NC_bb* ncbbp,
                      __attribute__((unused)) MPI_Info info)
{
    int i, rank, np, err, flag, masterrank, procname_len;
    char logbase[NC_LOG_MAX_PATH], basename[NC_LOG_MAX_PATH];
    char procname[MPI_MAX_PROCESSOR_NAME];
    char *abspath, *fname, *path, *fdir = NULL;
    char *logbasep = ".";
#ifdef PNETCDF_PROFILING
    double t1, t2;
#endif
    DIR *logdir;
    ssize_t headersize;
    NC_bb_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    t1 = MPI_Wtime();
#endif

    /* Get rank and number of processes */
    err = MPI_Comm_rank(ncbbp->comm, &rank);
    if (err != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(err, "MPI_Comm_rank");
        DEBUG_RETURN_ERROR(err);
    }
    err = MPI_Comm_size(ncbbp->comm, &np);
    if (err != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(err, "MPI_Comm_rank");
        DEBUG_RETURN_ERROR(err);
    }
    masterrank = rank;

    err = MPI_Get_processor_name(procname, &procname_len);
    if (err != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(err, "MPI_Get_processor_name");
        DEBUG_RETURN_ERROR(err);
    }

    /* Initialize log structure */

    /* Determine log file name
     * Log file name is $(bufferdir)$(basename)_$(ncid)_$(rank).{meta/data}
     * filepath is absolute path to the cdf file
     * If buffer directory is not set, we use the same directory as the NetCDF file
     */

    /* Read environment variable for burst buffer path */

    /* Remove romio driver specifier form the beginning of path */
    path = strchr(ncbbp->path, ':');
    if (path == NULL) {
        /* No driver specifier, use full path */
        path = ncbbp->path;
    }
    else {
        /* Skip until after the first ':' */
        path += 1;
    }

    /* Determine log base */
    if (ncbbp->logbase[0] != '\0') {
        /* We don't need driver specifier in logbase as well */
        logbasep = strchr(ncbbp->logbase, ':');
        if (logbasep == NULL) {
            /* No driver specifier, use full path */
            logbasep = ncbbp->logbase;
        }
        else {
            /* Skip until after the first ':' */
            logbasep += 1;
        }
    }
    else {
        i = strlen(path);
        fdir = (char*)NCI_Malloc((i + 1) * sizeof(char));
        strncpy(fdir, path, i + 1);
        /* Search for first '\' from the back */
        for (i--; i > -1; i--) {
            if (fdir[i] == '/') {
                fdir[i + 1] = '\0';
                break;
            }
        }

        /* If directory is fund, use it as logbase */
        if (i >= 0) {
            logbasep = fdir;
        }

        /* Warn if log base not set by user */
        if (rank == 0) {
            printf("Warning: Log directory not set. Using %s\n", logbasep);
            fflush(stdout);
        }
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
    abspath = realpath(path, basename);
    if (abspath == NULL) {
        /* Can not resolve absolute path */
        DEBUG_RETURN_ERROR(NC_EBAD_FILE);
    }
    abspath = realpath(logbasep, logbase);
    if (abspath == NULL) {
        /* Can not resolve absolute path */
        DEBUG_RETURN_ERROR(NC_EBAD_FILE);
    }

    if (fdir != NULL) {
        NCI_Free(fdir);
    }

    /* Communicator for processes sharing log files */
    if (ncbbp->hints & NC_LOG_HINT_LOG_SHARE) {
#if MPI_VERSION >= 3
        err = MPI_Comm_split_type(ncbbp->comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &(ncbbp->logcomm));
        if (err != NC_NOERR) {
            DEBUG_RETURN_ERROR(NC_EMPI);
        }
#else
        err = ncbbio_get_node_comm(ncbbp->comm, &ncbbp->logcomm);
        if (err != NC_NOERR) {
            return err;
        }
#endif
        err = MPI_Bcast(&masterrank, 1, MPI_INT, 0, ncbbp->logcomm);
        if (err != NC_NOERR) {
            DEBUG_RETURN_ERROR(NC_EMPI);
        }
    }
    else
    {
        ncbbp->logcomm = MPI_COMM_SELF;
        masterrank = rank;
    }

    /* Extract file anme
     * Search for last / character
     * Absolute path should always contains one '/', return error otherwise
     * We return the string including '/' for convenience
     */
    fname = strrchr(basename, '/');
    if (fname == NULL) {
        DEBUG_RETURN_ERROR(NC_EBAD_FILE);
    }

    /* Log file path may also contain non-existing directory
     * We need to create them before we can search for usable id
     * As log file name hasn't been determined, we need to use a dummy one here
     */
    sprintf(ncbbp->metalogpath, "%s%s_%d_%d.meta", logbase, fname,
            ncbbp->ncid, masterrank);
    sprintf(ncbbp->datalogpath, "%s%s_%d_%d.data", logbase, fname,
            ncbbp->ncid, masterrank);

    /* Initialize metadata buffer */
    err = ncbbio_log_buffer_init(&(ncbbp->metadata));
    if (err != NC_NOERR) {
        return err;
    }

    /* Initialize metadata entry array */
    err = ncbbio_log_sizearray_init(&(ncbbp->entrydatasize));
    if (err != NC_NOERR) {
        return err;
    }

    /* Set log file descriptor to NULL */

#ifdef PNETCDF_PROFILING
    /* Performance counters */
    ncbbp->total_data = 0;
    ncbbp->total_meta = 0;
    ncbbp->max_buffer = 0;
    ncbbp->total_time = 0;
    ncbbp->create_time = 0;
    ncbbp->enddef_time = 0;
    ncbbp->put_time = 0;
    ncbbp->flush_time = 0;
    ncbbp->close_time = 0;
    ncbbp->flush_replay_time = 0;
    ncbbp->flush_data_rd_time = 0;
    ncbbp->flush_put_time = 0;
    ncbbp->flush_wait_time = 0;
    ncbbp->put_data_wr_time = 0;
    ncbbp->put_meta_wr_time = 0;
    ncbbp->put_num_wr_time = 0;
    ncbbp->max_buffer = 0;
#endif

    /* Misc */
    ncbbp->rank = rank;
    ncbbp->np = np;
    ncbbp->maxentrysize = 0;

    /* Initialize metadata header */

    /*
     * Allocate space for metadata header
     * Header consists of a fixed-size info and variable size basename
     */
    headersize = sizeof(NC_bb_metadataheader) + strlen(basename) + 1 + SIZEOF_INT + procname_len + 1;
    if (headersize % 16 != 0) {
        headersize += 16 - (headersize % 16);
    }
    headerp = (NC_bb_metadataheader*)ncbbio_log_buffer_alloc(&(ncbbp->metadata),
                                                             headersize);
    memset(headerp, 0, headersize);

    /* Fill up the metadata log header */
    memcpy(headerp->magic, NC_LOG_MAGIC, sizeof(headerp->magic));
    memcpy(headerp->format, NC_LOG_FORMAT_CDF_MAGIC, sizeof(headerp->format));
    strncpy((char*)headerp->basename, basename,
            headersize - SIZEOF_INT - sizeof(NC_bb_metadataheader) - 1);
    headerp->rank_id = rank;   /* Rank */
    headerp->num_ranks = np;   /* Number of processes */
    /* Without conversion before logging, data in native representation */
    headerp->is_external = 0;
    /* Determine Endianness */
#ifdef WORDS_BIGENDIAN
    headerp->big_endian = NC_LOG_TRUE;
#else
    headerp->big_endian = NC_LOG_FALSE;
#endif
    headerp->num_entries = 0;    /* The log is empty */
    /* Highest dimension among all variables */
    headerp->max_ndims = 0;
    /* Location of the first entry */
    headerp->entry_begin = ncbbp->metadata.nused;
    headerp->basenamelen = strlen(basename);

    /* Process anme */
    *(int*)((char*)headerp->basename + headerp->basenamelen + 1) = procname_len;
    strncpy((char*)headerp->basename + headerp->basenamelen + 5,
            procname, procname_len + 1);

    /* Create log files */
    flag = O_RDWR | O_CREAT;
    if (!(ncbbp->hints & NC_LOG_HINT_LOG_OVERWRITE)) {
        flag |= O_EXCL;
    }
    err = ncbbio_sharedfile_open(ncbbp->logcomm, ncbbp->metalogpath, flag,
                           MPI_INFO_NULL, &(ncbbp->metalog_fd));
    if (err != NC_NOERR) {
        return err;
    }
    err = ncbbio_sharedfile_open(ncbbp->logcomm, ncbbp->datalogpath, flag,
                           MPI_INFO_NULL, &(ncbbp->datalog_fd));
    if (err != NC_NOERR) {
        return err;
    }

    /* Write metadata header to file
     * Write from the memory buffer to file
     */
    err = ncbbio_sharedfile_write(ncbbp->metalog_fd, headerp, headersize);
    if (err != NC_NOERR) {
        return err;
    }

    /* Write data header to file
     * Data header consists of a fixed-length string "PnetCDF0"
     */
    err = ncbbio_sharedfile_write(ncbbp->datalog_fd, "PnetCDF0", 8);
    if (err != NC_NOERR) {
        return err;
    }

    ncbbp->datalogsize = 8;

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
    ncbbp->total_time += t2 - t1;
    ncbbp->create_time += t2 - t1;

    ncbbp->total_meta += headersize;
    ncbbp->total_data += 8;
#endif

    return NC_NOERR;
}

/*
 * Update information used by bb layer on enddef
 * IN    ncbbp:    NC_bb structure
 */
int ncbbio_log_enddef(NC_bb *ncbbp)
{
    int err;
    NC_bb_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    double t1, t2;
    t1 = MPI_Wtime();
#endif

    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;

    /*
     * Update the header if max ndim increased
     * For now, max_ndims is the only field that can change
     */
    if (ncbbp->max_ndims > headerp->max_ndims) {
        headerp->max_ndims = ncbbp->max_ndims;

        /* Overwrite maxndims
         * This marks the completion of the record
         */
        err = ncbbio_sharedfile_pwrite(ncbbp->metalog_fd, &headerp->max_ndims,
                                SIZEOF_MPI_OFFSET, sizeof(NC_bb_metadataheader) -
                               sizeof(headerp->basename) -
                               sizeof(headerp->basenamelen) -
                               sizeof(headerp->num_entries) -
                               sizeof(headerp->max_ndims));
        if (err != NC_NOERR) {
            return err;
        }
    }

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
    ncbbp->total_time += t2 - t1;
    ncbbp->enddef_time += t2 - t1;
#endif

    return NC_NOERR;
}

/*
 * Flush the log to CDF file and clean up the log structure
 * Used by ncmpi_close()
 * IN    ncbbp:    log structure
 */
int ncbbio_log_close(NC_bb *ncbbp,
                     int    replay)
{
    int err;
    NC_bb_metadataheader* headerp;

#ifdef PNETCDF_PROFILING
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
    double t1, t2;
    t1 = MPI_Wtime();
#endif

    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;

    /* If log file is created, flush the log */
    if (ncbbp->metalog_fd != NULL) {
        /* Commit to CDF file */
        if (replay && (headerp->num_entries > 0 || !(fIsSet(ncbbp->flag, NC_MODE_INDEP)))) {
            ncbbio_log_flush_core(ncbbp);
        }

        /* Close log file */
        err = ncbbio_sharedfile_close(ncbbp->metalog_fd);
        if (err != NC_NOERR) {
            return err;
        }
        err = ncbbio_sharedfile_close(ncbbp->datalog_fd);
        if (err != NC_NOERR) {
            return err;
        }

        /* Delete log files if delete flag is set */
        if (ncbbp->hints & NC_LOG_HINT_DEL_ON_CLOSE) {
            unlink(ncbbp->datalogpath);
            unlink(ncbbp->metalogpath);
        }
    }

    /* Free meta data buffer and metadata offset list*/
    ncbbio_log_buffer_free(&(ncbbp->metadata));
    ncbbio_log_sizearray_free(&(ncbbp->entrydatasize));

    /* Close shared log communicator */
    if (ncbbp->logcomm != MPI_COMM_SELF) {
        MPI_Comm_free(&(ncbbp->logcomm));
    }

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
    ncbbp->total_time += t2 - t1;
    ncbbp->close_time += t2 - t1;

    /* Print accounting info in debug build */
    MPI_Reduce(&(ncbbp->total_time), &total_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncbbp->comm);
    MPI_Reduce(&(ncbbp->create_time), &create_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncbbp->comm);
    MPI_Reduce(&(ncbbp->enddef_time), &enddef_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncbbp->comm);
    MPI_Reduce(&(ncbbp->put_time), &put_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncbbp->comm);
    MPI_Reduce(&(ncbbp->flush_time), &flush_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncbbp->comm);
    MPI_Reduce(&(ncbbp->close_time), &close_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                ncbbp->comm);
    MPI_Reduce(&(ncbbp->flush_replay_time), &flush_replay_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->flush_data_rd_time), &flush_data_rd_time, 1,
                MPI_DOUBLE, MPI_MAX, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->flush_put_time), &flush_put_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->flush_wait_time), &flush_wait_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->put_data_wr_time), &put_data_wr_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->put_meta_wr_time), &put_meta_wr_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->put_num_wr_time), &put_num_wr_time, 1, MPI_DOUBLE,
                MPI_MAX, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->total_meta), &total_meta, 1, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->total_data), &total_data, 1, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, 0, ncbbp->comm);
    MPI_Reduce(&(ncbbp->flushbuffersize), &buffer_size, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, ncbbp->comm);
#if 0
    if (ncbbp->rank == 0) {
        printf("==========================================================\n");
        printf("File: %s\n", ncbbp->path);
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
        printf("\tTime in ncbbio_log_flush_core: %lf\n", flush_time);
        printf("\tTime in ncbbio_log_close: %lf\n", close_time);
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
 * IN    ncbbp:    log structure
 */
int ncbbio_log_flush(NC_bb* ncbbp)
{
    int err, status = NC_NOERR;
    NC_bb_metadataheader *headerp;

#ifdef PNETCDF_PROFILING
    double t1, t2;
    t1 = MPI_Wtime();
#endif

    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;

    /* Nothing to replay if nothing have been written
     * We still need to participate the flush in collective mode
     * We assume some processes will have things to flush to save communication cost
     */
    if (headerp->num_entries == 0 && fIsSet(ncbbp->flag, NC_MODE_INDEP)) {
        return NC_NOERR;
    }

    /* Replay log file */
    err = ncbbio_log_flush_core(ncbbp);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) {
            DEBUG_ASSIGN_ERROR(status, err);
        }
    }

    /* Reset log status */

    /* Set num_entries to 0 */
    headerp->num_entries = 0;

    /* Overwrite num_entries
     * This marks the completion of flush
     */
    err = ncbbio_sharedfile_pwrite(ncbbp->metalog_fd, &headerp->num_entries,
                            SIZEOF_MPI_OFFSET, 56);
    if (err != NC_NOERR) {
        return err;
    }

    /* Reset metadata buffer and entry array status */
    ncbbp->metadata.nused = headerp->entry_begin;
    ncbbp->entrydatasize.nused = 0;
    ncbbp->metaidx.nused = 0;

    /* Rewind data log file descriptors and reset the size */
    err = ncbbio_sharedfile_seek(ncbbp->datalog_fd, 8, SEEK_SET);
    if (err != NC_NOERR) {
        return err;
    }

    ncbbp->datalogsize = 8;

#ifdef PNETCDF_PROFILING
    t2 = MPI_Wtime();
    ncbbp->total_time += t2 - t1;
    ncbbp->flush_time += t2 - t1;
#endif

    return status;
}
