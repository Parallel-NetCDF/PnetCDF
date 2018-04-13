/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <sys/types.h>
#include <dirent.h>
#include <assert.h>
#include "nc.h"
#include "ncx.h"
#include <limits.h>
#include <fcntl.h>
#include <errno.h>
#include <stdint.h>
#include <pnetcdf.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <pnc_debug.h>
#include <common.h>
#include <pnetcdf.h>
#include <log.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int ncbbio_log_check_header(NC* ncp, int num_entries){
    int i, err, np, rank, max_ndims;
    int data_fd, meta_fd;
    ssize_t ioret;
    size_t metasize, headersize;
    char *buffer, *abspath, *fname;
    char basename[NC_LOG_PATH_MAX];
    NC_bb_metadataheader *headerp;
    NC_bb *ncbbp = ncp->ncbbp;

    if (ncbbp == NULL){
        return NC_NOERR;
    }

    /* Get rank and number of processes */
    err = MPI_Comm_rank(ncp->comm, &rank);
    if (err != MPI_SUCCESS) {
        err = ncmpii_handle_error(err, "MPI_Comm_rank");
        DEBUG_RETURN_ERROR(err);
    }
    err = MPI_Comm_size(ncp->comm, &np);
    if (err != MPI_SUCCESS) {
        err = ncmpii_handle_error(err, "MPI_Comm_rank");
        DEBUG_RETURN_ERROR(err);
    }

    /* Open log file */
    meta_fd = open(ncbbp->metalogpath, O_RDONLY);
    if (meta_fd < 0){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    data_fd = open(ncbbp->datalogpath, O_RDONLY);
    if (data_fd < 0){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }

    /* Calculate metadata log header size */
    headerp = (NC_bb_metadataheader*)ncbbp->metadata.buffer;
    metasize = sizeof(NC_bb_metadataheader) + headerp->basenamelen;

    /* Allocate buffer */
    buffer = (char*)malloc(metasize);

    /* Read metadata log header */
    ioret = read(meta_fd, buffer, metasize);
    if (ioret < 0) {
        err = ncmpii_handle_io_error("read");
        if (err == NC_EFILE){
            err = NC_EREAD;
        }
        DEBUG_RETURN_ERROR(err);
    }
    if (ioret != metasize){
        DEBUG_RETURN_ERROR(NC_EBADLOG);
    }

    /* Resolve absolute path */
    abspath = realpath(ncp->path, basename);
    if (abspath == NULL){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }

    /* Check metadata log header */
    headerp = (NC_bb_metadataheader*)buffer;
    if (strncmp(headerp->magic, NC_LOG_MAGIC, NC_LOG_MAGIC_SIZE) != 0) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (strncmp(headerp->format, NC_LOG_FORMAT_CDF_MAGIC, NC_LOG_FORMAT_SIZE) != 0) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
#ifdef WORDS_BIGENDIAN
    if (!headerp->big_endian){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
#else
    if (headerp->big_endian){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
#endif
    if (headerp->is_external){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (headerp->num_ranks != np) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (headerp->rank_id != rank) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    headersize = sizeof(NC_bb_metadataheader) + strlen(basename);
    if (headersize % 4 != 0){
        headersize += 4 - (headersize % 4);
    }
    if (headerp->entry_begin != headersize){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    max_ndims = 0;
    for(i = 0; i < ncp->vars.ndefined; i++){
        if (ncp->vars.value[i]->ndims > max_ndims){
            max_ndims = ncp->vars.value[i]->ndims;
        }
    }
    if (headerp->max_ndims != max_ndims){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (num_entries >= 0) {
        if (headerp->num_entries != num_entries){
            DEBUG_RETURN_ERROR(NC_ELOGCHECK);
        }
    }

    if (headerp->basenamelen != strlen(basename)) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (strncmp(headerp->basename, basename, strlen(basename)) != 0) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }

    /* Read data log header */
    ioret = read(data_fd, buffer, 8);
    if (ioret < 0) {
        err = ncmpii_handle_io_error("read");
        if (err == NC_EFILE){
            err = NC_EREAD;
        }
        DEBUG_RETURN_ERROR(err);
    }
    if (ioret != 8){
        DEBUG_RETURN_ERROR(NC_EBADLOG);
    }
    /* Check data log header */
    if (strncmp(buffer, "PnetCDF0", 8) != 0) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }

    close(meta_fd);
    close(data_fd);

    return NC_NOERR;
}

int ncbbio_log_check_put(NC* ncp, int varid, int api_kind, int itype, int packedsize, MPI_Offset *start, MPI_Offset *count, MPI_Offset *stride, int num_entries){
    int i, err;
    int meta_fd, data_fd;
    ssize_t ioret;
    size_t datasize, metasize, esize;
    char *buffer, *abspath, *fname;
    char basename[NC_LOG_PATH_MAX];
    struct stat logstat;
    MPI_Offset *_start, *_count, *_stride;
    NC_bb_metadataentry *entryp;
    NC_bb_metadataheader *headerp;
    NC_bb *ncbbp = ncp->ncbbp;
    NC_var *varp;

    if (ncbbp == NULL){
        return NC_NOERR;
    }

    /* Open log file */
    meta_fd = open(ncbbp->metalogpath, O_RDONLY);
    if (meta_fd < 0){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    data_fd = open(ncbbp->datalogpath, O_RDONLY);
    if (meta_fd < 0){
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }

    /* Calculate size of metadata log and data log */
    metasize = ncbbp->metadata.nused;
    ioret = lseek(ncbbp->datalog_fd, 0, SEEK_CUR);
    if (ioret < 0){
        DEBUG_RETURN_ERROR(ncmpii_handle_io_error("lseek"));
    }
    datasize = ioret;

    /* Get variable */
    err = ncmpii_NC_lookupvar(ncp, varid, &varp);
    if (err != NC_NOERR) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }

    /* Calculate log entry size and allocate buffer */
    esize = sizeof(NC_bb_metadataentry) + varp->ndims * SIZEOF_MPI_OFFSET * 3;
    if (esize > sizeof(NC_bb_metadataheader)){
        buffer = (char*)malloc(esize);
    }
    else{
        buffer = (char*)malloc(sizeof(NC_bb_metadataheader));
    }

    /* Check num_entries in header */
    if (num_entries >= 0) {
        ioret = read(meta_fd, buffer, sizeof(NC_bb_metadataheader));
        if (ioret < 0) {
            err = ncmpii_handle_io_error("read");
            if (err == NC_EFILE){
                err = NC_EREAD;
            }
            DEBUG_RETURN_ERROR(err);
        }
        if (ioret != sizeof(NC_bb_metadataheader)){
            DEBUG_RETURN_ERROR(NC_EBADLOG);
        }
        headerp = (NC_bb_metadataheader*)buffer;

        if (headerp->num_entries != num_entries){
            DEBUG_RETURN_ERROR(NC_ELOGCHECK);
        }
    }

    /* Read last entry */
    ioret = lseek(meta_fd, metasize - esize, SEEK_SET);
    if (ioret < 0){
        DEBUG_RETURN_ERROR(ncmpii_handle_io_error("lseek"));
    }
    ioret = read(meta_fd, buffer, esize);
    if (ioret < 0) {
        err = ncmpii_handle_io_error("read");
        if (err == NC_EFILE){
            err = NC_EREAD;
        }
        DEBUG_RETURN_ERROR(err);
    }
    if (ioret != esize){
        DEBUG_RETURN_ERROR(NC_EBADLOG);
    }

    /* Check entry */
    entryp = (NC_bb_metadataentry*)buffer;
    _start = (MPI_Offset*)(entryp + 1);
    _count = _start + entryp->ndims;
    _stride = _count + entryp->ndims;

    if (entryp->ndims != varp->ndims) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (entryp->api_kind != api_kind) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (entryp->itype != itype) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (entryp->varid != varid) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (entryp->data_len != packedsize) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }
    if (entryp->data_off != datasize - entryp->data_len) {
        DEBUG_RETURN_ERROR(NC_ELOGCHECK);
    }

    for(i = 0; i < varp->ndims; i++){
        if(start[i] != _start[i]){
            DEBUG_RETURN_ERROR(NC_ELOGCHECK);
        }
        if(count[i] != _count[i]){
            DEBUG_RETURN_ERROR(NC_ELOGCHECK);
        }
    }

    if (api_kind == NC_LOG_API_KIND_VARS) {
        for(i = 0; i < varp->ndims; i++){
            if(stride[i] != _stride[i]){
                DEBUG_RETURN_ERROR(NC_ELOGCHECK);
            }
        }
    }

    close(meta_fd);
    close(data_fd);

    return NC_NOERR;
}



