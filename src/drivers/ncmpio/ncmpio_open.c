/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_open() : dispatcher->open()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy() */
#ifdef HAVE_ACCESS
#include <unistd.h>  /* access() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"
#ifdef ENABLE_SUBFILING
#include "ncmpio_subfile.h"
#endif

/*----< ncmpio_open() >------------------------------------------------------*/
int
ncmpio_open(MPI_Comm     comm,
            const char  *path,
            int          omode,
            int          ncid,
            int          env_mode,
            MPI_Info     user_info, /* user's and env info combined */
            void       **ncpp)
{
    char *filename, value[MPI_MAX_INFO_VAL + 1], *mpi_name;
    int i, rank, nprocs, mpiomode, err, status=NC_NOERR, mpireturn, flag;
    int striping_unit;
    MPI_File fh=MPI_FILE_NULL;
    NC *ncp=NULL;

    *ncpp = NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* Note path's validity and omode consistency have been checked in
     * ncmpi_open() in src/dispatchers/file.c and path consistency will be done
     * in MPI_File_open.
     */

    /* First, check whether omode is valid or supported ---------------------*/

    /* NC_DISKLESS is not supported yet */
    if (omode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* NC_MMAP is not supported yet */
    if (omode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* allocate buffer for header object NC and initialize its contents */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    *ncpp = (void*)ncp;

    ncp->ncid     = ncid;
    ncp->comm     = comm;     /* reuse comm duplicated in dispatch layer */
    ncp->rank     = rank;
    ncp->nprocs   = nprocs;
    ncp->mpiinfo  = MPI_INFO_NULL;

    /* Extract hints from user_info. Two hints must be extracted now in order
     * to continue:
     *     nc_pncio: whether to user MPI-IO or PnetCDF's PNCIO driver.
     *     nc_num_aggrs_per_node: number of processes per node to be INA
     *     aggregators.
     *
     * ncp->fstype will be initialized in ncmpio_hint_extract() and set in
     * PNCIO_FileSysType().
     */
    ncmpio_hint_extract(ncp, user_info);

    if (ncp->fstype == PNCIO_FSTYPE_CHECK) {
        if (rank == 0)
            /* Check file system type. If the given file does not exist, check
             * its parent folder. Currently PnetCDF's PNCIO drivers support
             * Lustre (PNCIO_LUSTRE) and Unix File System (PNCIO_UFS).
             */
            ncp->fstype = PNCIO_FileSysType(path);

        MPI_Bcast(&ncp->fstype, 1, MPI_INT, 0, ncp->comm);
    }

#ifdef WKL_DEBUG
if (rank == 0) printf("%s at %d fstype=%s\n", __func__,__LINE__,(ncp->fstype == PNCIO_FSTYPE_MPIIO)? "PNCIO_FSTYPE_MPIIO" : (ncp->fstype == PNCIO_LUSTRE) ? "PNCIO_LUSTRE" : "PNCIO_UFS");
#endif

    /* Remove the file system type prefix name if there is any. For example,
     * when path = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * filename pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * access() below
     */
    filename = ncmpii_remove_file_system_type_prefix(path);

    ncp->path     = path;  /* reuse path duplicated in dispatch layer */
    ncp->pncio_fh = NULL;
    ncp->iomode   = omode;

    ncp->collective_fh  = MPI_FILE_NULL;
    ncp->independent_fh = MPI_FILE_NULL;

    /* Setting file open mode in mpiomode which may later be needed in
     * ncmpi_begin_indep_data() to open file for independent data mode.
     */
    mpiomode = fIsSet(omode, NC_WRITE) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;
    ncp->mpiomode = mpiomode;

    /* PnetCDF default fill mode is no fill */
    fClr(ncp->flags, NC_MODE_FILL);

    /* set read-only mode */
    if (!fIsSet(omode, NC_WRITE)) fSet(ncp->flags, NC_MODE_RDONLY);

    fSet(ncp->flags, env_mode);

    /* Construct a list of unique IDs of compute nodes allocated to this job
     * and save it in ncp->node_ids[nprocs], which contains node IDs of each
     * rank. The node IDs are used either when intra-node aggregation is
     * enabled or when using PnetCDF's PNCIO driver.
     *
     * When intra-node aggregation is enabled, node IDs are used to create a
     * new MPI communicator consisting of the intra-node aggregators only. The
     * communicator will be used to call file open in MPI-IO or PnetCDF's PNCIO
     * driver. This means only intra-node aggregators will perform file I/O in
     * PnetCDF collective put and get operations.
     */
    ncp->node_ids = NULL;
    if (ncp->fstype != PNCIO_FSTYPE_MPIIO || ncp->num_aggrs_per_node != 0) {
        err = ncmpii_construct_node_list(comm, &ncp->num_nodes, &ncp->node_ids);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);

        /* When the total number of aggregators >= number of processes, disable
         * intra-node aggregation.
         */
        if (ncp->num_aggrs_per_node * ncp->num_nodes >= ncp->nprocs)
            ncp->num_aggrs_per_node = 0;
    }

    /* ncp->num_aggrs_per_node = 0, or > 0 indicates whether this feature
     * is disabled or enabled globally for all processes.
     */
    ncp->my_aggr = -1;
    ncp->ina_comm = MPI_COMM_NULL;
    ncp->ina_nprocs = 0;
    ncp->ina_rank = -1;
    ncp->ina_node_list = NULL;
    if (ncp->num_aggrs_per_node > 0) {
        /* Divide all ranks into groups. Each group is assigned with one
         * intra-node aggregator. The following metadata related to intra-node
         * aggregation will be set up.
         * ncp->my_aggr is the aggregator's rank ID of this group. When ==
         *     ncp->rank, this rank is an aggregator.
         * ncp->num_nonaggrs is the number of non-aggregators assigned to this
         *     rank (an aggregator)
         * ncp->ina_comm will be created consisting of only intra-node
         *     aggregators, which will be used when calling MPI_File_open().
         *     For non-aggregator, ncp->ina_comm == MPI_COMM_NULL.
         * ncp->node_ids[] will be modified to contain the nodes IDs of
         *     intra-node aggregators only, which will be passed to pncio_fh.
         */
        err = ncmpio_ina_init(ncp);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);

        /* As non-aggregators will not perform any file I/O, we now can replace
         * comm with ina_comm. Same for nprocs.
         */
        comm = ncp->ina_comm;
        nprocs = ncp->ina_nprocs;

        /* For non-aggregators, comm is MPI_COMM_NULL. As the remaining task of
         * this subroutine is to open the file and obtain the file handler,
         * non-aggregators can skip.
         */
        if (comm == MPI_COMM_NULL) {
            if (user_info != MPI_INFO_NULL)
                MPI_Info_dup(user_info, &ncp->mpiinfo);
            goto fn_exit;
        }
    }

    /* open file collectively ---------------------------------------------- */
    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        TRACE_IO(MPI_File_open, (comm, path, mpiomode, user_info, &fh));
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            DEBUG_FOPEN_ERROR(err);
        }

        /* Now the file has been successfully opened */
        ncp->collective_fh  = fh;
        ncp->independent_fh = (nprocs > 1) ? MPI_FILE_NULL : fh;

        /* get the I/O hints used/modified by MPI-IO */
        TRACE_IO(MPI_File_get_info, (fh, &ncp->mpiinfo));
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            DEBUG_FOPEN_ERROR(err);
        }
    }
    else {
        /* When ncp->fstype != PNCIO_FSTYPE_MPIIO, use PnetCDF's PNCIO driver */
        ncp->pncio_fh = (PNCIO_File*) NCI_Calloc(1,sizeof(PNCIO_File));
        ncp->pncio_fh->file_system = ncp->fstype;
        ncp->pncio_fh->num_nodes   = ncp->num_nodes;
        ncp->pncio_fh->node_ids    = ncp->node_ids;

        err = PNCIO_File_open(comm, filename, mpiomode, user_info,
                              ncp->pncio_fh);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);

        /* Now the file has been successfully opened, obtain the I/O hints
         * used/modified by PNCIO driver.
         */
        err = PNCIO_File_get_info(ncp->pncio_fh, &ncp->mpiinfo);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);
    }

fn_exit:
    striping_unit = -1;

    if (ncp->num_aggrs_per_node > 0) {
        /* When intra-node aggregation is enabled, it is necessary to make sure
         * non-aggregators obtain consistent values of file striping hints.
         *
         * non-aggregator do not have hints returned from MPI_File_get_info()
         */
        int striping_info[2];
        if (ncp->rank == 0) {
            MPI_Info_get(ncp->mpiinfo, "striping_unit", MPI_MAX_INFO_VAL-1,
                         value, &flag);
            striping_info[0] = 0;
            if (flag) {
                errno = 0;  /* errno must set to zero before calling strtoll */
                striping_info[0] = (int)strtol(value,NULL,10);
                if (errno != 0) striping_info[0] = 0;
            }

            MPI_Info_get(ncp->mpiinfo, "striping_factor", MPI_MAX_INFO_VAL-1,
                         value, &flag);
            striping_info[1] = 0;
            if (flag) {
                errno = 0;  /* errno must set to zero before calling strtoll */
                striping_info[1] = (int)strtol(value,NULL,10);
                if (errno != 0) striping_info[1] = 0;
            }
        }

        MPI_Bcast(striping_info, 2, MPI_INT, 0, ncp->comm);

        if (ncp->my_aggr != ncp->rank) {
            sprintf(value, "%d", striping_info[0]);
            MPI_Info_set(ncp->mpiinfo, "striping_unit", value);
            sprintf(value, "%d", striping_info[1]);
            MPI_Info_set(ncp->mpiinfo, "striping_factor", value);
        }

        striping_unit = striping_info[0];
    }
    else {
        MPI_Info_get(ncp->mpiinfo, "striping_unit", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) {
            errno = 0;  /* errno must set to zero before calling strtoll */
            striping_unit = (int)strtol(value,NULL,10);
            if (errno != 0) striping_unit = -1;
        }
    }

    if (ncp->data_chunk == -1)
        /* if not set by user hint, nc_data_move_chunk_size */
        ncp->data_chunk = (striping_unit > 0) ? striping_unit
                                              : PNC_DATA_MOVE_CHUNK_SIZE;

    /* Copy MPI-IO hints into ncp->mpiinfo */
    ncmpio_hint_set(ncp, ncp->mpiinfo);

    /* ina_node_list is no longer needed */
    if (ncp->ina_node_list != NULL) {
        NCI_Free(ncp->ina_node_list);
        ncp->ina_node_list = NULL;
    }
    /* node_ids is no longer needed */
    if (ncp->node_ids != NULL) {
        NCI_Free(ncp->node_ids);
        ncp->node_ids = NULL;
    }
    if (ncp->pncio_fh != NULL)
        ncp->pncio_fh->node_ids = NULL;

    /* read header from file into NC object pointed by ncp -------------------*/
    err = ncmpio_hdr_get_NC(ncp);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) { /* fatal error */
        ncmpio_file_close(ncp);
        if (ncp->ina_comm != MPI_COMM_NULL) MPI_Comm_free(&ncp->ina_comm);
        ncmpio_free_NC(ncp);
        DEBUG_RETURN_ERROR(err);
    }

#ifdef ENABLE_SUBFILING
    if (ncp->subfile_mode) {
        /* check subfiling attribute */
        err = ncmpio_get_att(ncp, NC_GLOBAL, "_PnetCDF_SubFiling.num_subfiles",
                             &ncp->num_subfiles, MPI_INT);
        if (err == NC_NOERR && ncp->num_subfiles > 1) {
            /* ignore error NC_ENOTATT if this attribute is not defined */
            for (i=0; i<ncp->vars.ndefined; i++) {
                /* variables may have different numbers of subfiles */
                err = ncmpio_get_att(ncp, i, "_PnetCDF_SubFiling.num_subfiles",
                             &ncp->vars.value[i]->num_subfiles,MPI_INT);
                if (err == NC_ENOTATT) continue;
                if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);
                if (ncp->vars.value[i]->num_subfiles > 1) {
                    /* find the orginal ndims of variable i */
                    err = ncmpio_get_att(ncp,i,"_PnetCDF_SubFiling.ndims_org",
                                 &ncp->vars.value[i]->ndims_org,MPI_INT);
                    if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);
                    ncp->vars.value[i]->dimids_org = (int*) NCI_Malloc(
                              ncp->vars.value[i]->ndims_org * SIZEOF_INT);
                    err = ncmpio_get_att(ncp,i,"_PnetCDF_SubFiling.dimids_org",
                              ncp->vars.value[i]->dimids_org, MPI_INT);
                    if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);
                }
            }
            /* open subfile */
            err = ncmpio_subfile_open(ncp);
            if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);
        }
        else ncp->num_subfiles = 0;
    }
    else
        ncp->num_subfiles = 0;
#endif

#ifndef SEARCH_NAME_LINEARLY
    /* initialize and populate name lookup tables ---------------------------*/
    ncmpio_hash_table_populate_NC_dim(&ncp->dims, ncp->dims.hash_size);
    ncmpio_hash_table_populate_NC_var(&ncp->vars, ncp->vars.hash_size);
    ncmpio_hash_table_populate_NC_attr(ncp);
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.value[i]->attrs.hash_size = ncp->hash_size_attr;
#endif

    return status;
}

