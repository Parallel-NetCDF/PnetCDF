/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <stdio.h>
#include <math.h>
#include <ncbbio_driver.h>

/*
 * Extract hints related to the burst buffering feature
 */
void ncbbio_extract_hint(NC_bb *ncbbp, MPI_Info info) {
    int flag;
    char value[MPI_MAX_INFO_VAL];

    ncbbp->hints = NC_LOG_HINT_DEL_ON_CLOSE | NC_LOG_HINT_FLUSH_ON_READ |
                   NC_LOG_HINT_FLUSH_ON_SYNC;

    /* Directory to store log files */
    MPI_Info_get(info, "nc_burst_buf_dirname", MPI_MAX_INFO_VAL - 1,
                 ncbbp->logbase, &flag);
    if (!flag)
        ncbbp->logbase[0] = '\0';

    /* Overwrite the log file if already exists (disable) */
    MPI_Info_get(info, "nc_burst_buf_overwrite", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag && strcasecmp(value, "enable") == 0){
        ncbbp->hints |= NC_LOG_HINT_LOG_OVERWRITE;
    }

    /* Use shared logfiles among processes on the same compute node (default is
     * disabled). This feature depends on the availability of MPI constant
     * MPI_COMM_TYPE_SHARED, which is first defined in MPI standard version 3.0
     */
    MPI_Info_get(info, "nc_burst_buf_shared_logs", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag && strcasecmp(value, "enable") == 0){
        ncbbp->hints |= NC_LOG_HINT_LOG_SHARE;
    }

    /* Delete the log file after file closing (enable) */
    MPI_Info_get(info, "nc_burst_buf_del_on_close", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag && strcasecmp(value, "disable") == 0){
        ncbbp->hints ^= NC_LOG_HINT_DEL_ON_CLOSE;
    }

    /* Buffer size used to flush the log (0 (unlimited)) */
    MPI_Info_get(info, "nc_burst_buf_flush_buffer_size", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag){
        long int bsize = strtol(value, NULL, 0);
        if (bsize < 0) {
            bsize = 0;
        }
        ncbbp->flushbuffersize = (size_t)bsize; /* Unit: byte */
    }
    else{
        ncbbp->flushbuffersize = 0; /* 0 means unlimited} */
    }
}

/*
 * Export I/O hints to user info object.
 * NOTE: Here, we only export hints related to burst buffer feature.
 */
void ncbbio_export_hint(NC_bb *ncbbp, MPI_Info *info) {

    MPI_Info_set(*info, "nc_burst_buf", "enable");
    if (ncbbp->hints & NC_LOG_HINT_LOG_OVERWRITE)
        MPI_Info_set(*info, "nc_burst_buf_overwrite", "enable");

    if (ncbbp->hints & NC_LOG_HINT_LOG_SHARE)
        MPI_Info_set(*info, "nc_burst_buf_shared_logs", "enable");

    if (!(ncbbp->hints & NC_LOG_HINT_DEL_ON_CLOSE))
        MPI_Info_set(*info, "nc_burst_buf_del_on_close", "disable");

    if (ncbbp->logbase[0] != '\0')
        MPI_Info_set(*info, "nc_burst_buf_dirname", ncbbp->logbase);

    if (ncbbp->flushbuffersize > 0) {
        char value[MPI_MAX_INFO_VAL];
        sprintf(value, "%llu", ncbbp->flushbuffersize);
        MPI_Info_set(*info, "nc_burst_buf_flush_buffer_size", value);
    }
}

