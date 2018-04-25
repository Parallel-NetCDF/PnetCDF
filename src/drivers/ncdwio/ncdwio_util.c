/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <sys/types.h>
#include "ncx.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <pnc_debug.h>
#include <common.h>
#include <pnetcdf.h>
#include <ncdwio_driver.h>

/*
 * Extract mpi hints and set up the flags
 */
void ncdwio_extract_hint(NC_dw *ncdwp, MPI_Info info){
    int flag;
    char value[MPI_MAX_INFO_VAL];

    ncdwp->hints = NC_LOG_HINT_DEL_ON_CLOSE | NC_LOG_HINT_FLUSH_ON_READ |
                    NC_LOG_HINT_FLUSH_ON_SYNC;
    // Directory to place log files
    MPI_Info_get(info, "nc_dw_dirname", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag) {
        strncpy(ncdwp->logbase, value, PATH_MAX);
    }
    else {
        memset(ncdwp->logbase, 0, sizeof(ncdwp->logbase));
    }
    // Overwrite the logfile is already exists (disable)
    MPI_Info_get(info, "nc_dw_overwrite", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag && strcasecmp(value, "enable") == 0){
        ncdwp->hints |= NC_LOG_HINT_LOG_OVERWRITE;
    }
    // Use shared logfile (disable)
    MPI_Info_get(info, "nc_dw_shared_logs", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag && strcasecmp(value, "enable") == 0){
        ncdwp->hints |= NC_LOG_HINT_LOG_SHARE;
    }
    // Delete the logfile after file closing (enable)
    MPI_Info_get(info, "nc_dw_del_on_close", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag && strcasecmp(value, "disable") == 0){
        ncdwp->hints ^= NC_LOG_HINT_DEL_ON_CLOSE;
    }
    // Buffer size used to flush the log (0 (unlimited))
    MPI_Info_get(info, "nc_dw_flush_buffer_size", MPI_MAX_INFO_VAL - 1,
                 value, &flag);
    if (flag){
        long int bsize = strtol(value, NULL, 0);
        if (bsize < 0) {
            bsize = 0;
        }
        ncdwp->flushbuffersize = (size_t)bsize; // Unit: byte
    }
    else{
        ncdwp->flushbuffersize = 0; // 0 means unlimited}
    }
}

/*
 * Export hint based on flag
 * NOTE: We only set up the hint if it is not the default setting
 *       user hint maching the default behavior will be ignored
 */
void ncdwio_export_hint(NC_dw *ncdwp, MPI_Info info){
    char value[MPI_MAX_INFO_VAL];

    MPI_Info_set(info, "nc_dw", "enable");
    if (ncdwp->hints & NC_LOG_HINT_LOG_OVERWRITE) {
        MPI_Info_set(info, "nc_dw_overwrite", "enable");
    }
    if (ncdwp->hints & NC_LOG_HINT_LOG_SHARE) {
        MPI_Info_set(info, "nc_dw_shared_logs", "enable");
    }
    if (!(ncdwp->hints & NC_LOG_HINT_DEL_ON_CLOSE)) {
        MPI_Info_set(info, "nc_dw_del_on_close", "disable");
    }
    if (ncdwp->logbase[0] != '\0') {
        MPI_Info_set(info, "nc_dw_dirname", ncdwp->logbase);
    }
    if (ncdwp->flushbuffersize > 0) {
        sprintf(value, "%llu", ncdwp->flushbuffersize);
        MPI_Info_set(info, "nc_dw_flush_buffer_size", value);
    }
}
