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
ncmpio_open(MPI_Comm         comm,
            const char      *path,
            int              omode,
            int              ncid,
            int              env_mode,
            MPI_Info         user_info, /* user's and env info combined */
            PNC_comm_attr    comm_attr, /* node IDs and INA metadata */
            void           **ncpp)
{
    char *filename, value[MPI_MAX_INFO_VAL + 1], *mpi_name;
    int i, rank, nprocs, mpi_amode, err, status=NC_NOERR, mpireturn, flag;
    int striping_info[2];
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
     *     nc_driver: whether to user MPI-IO or PnetCDF's PNCIO driver.
     *     nc_num_aggrs_per_node: number of processes per node to be INA
     *     aggregators.
     *
     * ncp->driver is initialized in ncmpio_hint_extract().
     * ncp->fstype is set in PNCIO_FileSysType().
     */
    ncmpio_hint_extract(ncp, user_info);

    if (rank == 0)
        /* Check file system type. If the given file does not exist, check
         * its parent folder. Currently PnetCDF's PNCIO drivers support
         * Lustre (PNCIO_FS_LUSTRE) and Unix File System (PNCIO_FS_UFS).
         */
        ncp->fstype = PNCIO_FileSysType(path);

    MPI_Bcast(&ncp->fstype, 1, MPI_INT, 0, ncp->comm);

    /* Remove the file system type prefix name if there is any. For example,
     * when path = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * filename pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * access() below
     */
    filename = ncmpii_remove_file_system_type_prefix(path);

    ncp->path     = path;  /* reuse path duplicated in dispatch layer */
    ncp->nc_amode = omode;

    if (ncp->driver == PNC_DRIVER_PNCIO) {
        /* Initialize pncio_fh, PNCIO file handler, with common metadata shared
         * among all processes, including non-INA aggregators when INA is
         * enabled. This is necessary for non-INA aggregators to perform
         * independent I/O.
         */
        ncp->pncio_fh = (PNCIO_File*) NCI_Calloc(1, sizeof(PNCIO_File));
        ncp->pncio_fh->comm           = comm;
        ncp->pncio_fh->fstype         = ncp->fstype;
        ncp->pncio_fh->comm_attr      = comm_attr;
        ncp->pncio_fh->file_view.size = -1;
        ncp->pncio_fh->filename       = filename;
        ncp->pncio_fh->info           = MPI_INFO_NULL;
        ncp->pncio_fh->amode = fIsSet(omode, NC_WRITE) ? O_RDWR : O_RDONLY;
    }
    else
        ncp->pncio_fh  = NULL; /* used only when using PNCIO driver */

    ncp->collective_fh  = MPI_FILE_NULL;
    ncp->independent_fh = MPI_FILE_NULL;

    /* Setting file open mode in mpi_amode which may later be needed in
     * ncmpi_begin_indep_data() to open file for independent data mode.
     */
    mpi_amode = fIsSet(omode, NC_WRITE) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;
    ncp->mpi_amode = mpi_amode;

    /* PnetCDF default fill mode is no fill */
    fClr(ncp->flags, NC_MODE_FILL);

    /* set read-only mode */
    if (!fIsSet(omode, NC_WRITE)) fSet(ncp->flags, NC_MODE_RDONLY);

    fSet(ncp->flags, env_mode);

    /* comm_attr.ids[] stores a list of unique IDs of compute nodes of all MPI
     * ranks in the MPI communicator passed from the user application. It is a
     * keyval attribute cached in the communicator. See src/dispatchers/file.c
     * for details.
     *
     * commm_attr also stores the INA metadata. The INA communicator has been
     * created at the dispatcher.
     *
     * When intra-node aggregation (INA) is enabled, node IDs are used to
     * create a new MPI communicator consisting of the intra-node aggregators
     * only. The communicator will be used to call file open in MPI-IO or
     * PnetCDF's PNCIO driver. This means only intra-node aggregators will
     * perform file I/O in PnetCDF collective put and get operations.
     *
     * comm_attr.ids[] will be used to calculate cb_nodes, the number of
     * MPI-IO/PNCIO aggregators (not INA aggregators).
     */
    ncp->comm_attr = comm_attr;

    /* When the total number of aggregators >= number of processes, disable
     * intra-node aggregation.
     */
    if (ncp->num_aggrs_per_node * comm_attr.num_nodes >= ncp->nprocs)
        ncp->num_aggrs_per_node = 0;

    /* ncp->num_aggrs_per_node = 0, or > 0 indicates whether this feature
     * is disabled or enabled globally for all processes.
     */
    ncp->ina_nprocs = 0;
    ncp->ina_rank = -1;

    if (ncp->num_aggrs_per_node > 0) {
        if (ncp->rank == comm_attr.my_aggr) {
            int *ids;

            MPI_Comm_size(comm_attr.ina_inter_comm, &ncp->ina_nprocs);
            MPI_Comm_rank(comm_attr.ina_inter_comm, &ncp->ina_rank);

            /* overwrite comm_attr.ids[] to make it to contain the node IDs of
             * processes in the INA communicator.
             */
            ids = (int*) NCI_Malloc(sizeof(int) * ncp->ina_nprocs);

            for (i=0; i<ncp->ina_nprocs; i++) {
                int j;
                /* j is the process rank in comm passed into ncmpi_create() */
                j = comm_attr.ina_ranks[i];
                ids[i] = comm_attr.ids[j];
                /* Now ids[] store the node IDs of the processes in the INA
                 * communicator. ncp->ids[] will be used by PnetCDF's PNCIO
                 * driver only.
                 */
            }
            ncp->comm_attr.ids = ids;
        }

        /* As non-aggregators will not perform any file I/O, we now can replace
         * comm with ina_inter_comm. Same for nprocs.
         */
        comm = comm_attr.ina_inter_comm;
        nprocs = ncp->ina_nprocs;

        /* For non-INA aggregators, we keep comm, a local variable in this
         * subroutine, to be comm_attr.ina_inter_comm, a MPI_COMM_NULL for
         * non-INA aggregators. Because the remaining lines of codes below till
         * label 'fn_exit' is for INA aggregators to open the file and obtain a
         * file handler, in which non-INA aggregators do not participate and
         * thus do not make use of comm.
         */
        if (comm == MPI_COMM_NULL) {
            if (user_info != MPI_INFO_NULL)
                MPI_Info_dup(user_info, &ncp->mpiinfo);
            goto fn_exit;
        }
    }

    /* open file collectively ---------------------------------------------- */
    if (ncp->driver == PNC_DRIVER_MPIIO) {
#ifdef MPICH_VERSION
        /* MPICH recognizes file system type acronym prefixed to file names */
        TRACE_IO(MPI_File_open, (comm, path, mpi_amode, user_info, &fh));
#else
        TRACE_IO(MPI_File_open, (comm, filename, mpi_amode, user_info, &fh));
#endif
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
        /* When ncp->driver == PNC_DRIVER_PNCIO, use PnetCDF's PNCIO driver */
        int amode;

        ncp->pncio_fh->fstype = ncp->fstype;
        ncp->pncio_fh->comm_attr = ncp->comm_attr;

        amode = fIsSet(omode, NC_WRITE) ? O_RDWR : O_RDONLY;

        err = PNCIO_File_open(comm, filename, amode, user_info, ncp->pncio_fh);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);

        /* Now the file has been successfully opened, obtain the I/O hints
         * used/modified by PNCIO driver.
         */
        err = PNCIO_File_get_info(ncp->pncio_fh, &ncp->mpiinfo);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err);
    }

fn_exit:
    ncp->striping_unit = 0;
    ncp->striping_factor = 0;

    /* All processes must obtain striping_unit and striping_factor consistent
     * across all processes. They are used to fulfill ncmpi_inq_striping() and
     * striping_unit is also used to set ncp->data_chunk if hint
     * nc_data_move_chunk_size is not set by the user.
     */
    if (ncp->num_aggrs_per_node == 0 ||
        (ncp->num_aggrs_per_node > 0 && ncp->rank == 0)) {
        /* When INA is enabled, only INA aggregators have valid striping_unit
         * and striping_factor in their ncp->mpiinfo. This is because non-INA
         * aggregators do not participate the call to MPI_File_open() or
         * PNCIO_File_open(). We need to have root to broadcast these 2 info.
         */
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

    if (ncp->num_aggrs_per_node > 0) {
        MPI_Bcast(striping_info, 2, MPI_INT, 0, ncp->comm);

        if (ncp->ina_rank < 0) {
            /* Only non-INA aggregators need to add these to ncp->mpiinfo */
            sprintf(value, "%d", striping_info[0]);
            MPI_Info_set(ncp->mpiinfo, "striping_unit", value);
            sprintf(value, "%d", striping_info[1]);
            MPI_Info_set(ncp->mpiinfo, "striping_factor", value);

            if (ncp->driver == PNC_DRIVER_PNCIO) {
                /* Must initialize hints for non-INA aggregators, as they will
                 * be used when non-INA aggregators perform independent
                 * reads/writes. Currently, the following hints are used.
                 *      romio_ds_write, ind_wr_buffer_size
                 *      romio_ds_read,  ind_rd_buffer_size
                 */
                ncp->pncio_fh->comm = MPI_COMM_SELF;
                ncp->pncio_fh->hints = (PNCIO_Hints*)
                                        NCI_Calloc(1, sizeof(PNCIO_Hints));
                PNCIO_File_set_info(ncp->pncio_fh, ncp->mpiinfo);
            }
        }
    }

    ncp->striping_unit = striping_info[0];
    ncp->striping_factor = striping_info[1];

    if (ncp->data_chunk == -1)
        /* if hint nc_data_move_chunk_size is not set by the user */
        ncp->data_chunk = (ncp->striping_unit > 0) ? ncp->striping_unit
                                                   : PNC_DATA_MOVE_CHUNK_SIZE;

    /* Copy MPI-IO hints into ncp->mpiinfo */
    ncmpio_hint_set(ncp, ncp->mpiinfo);

    if (ncp->num_aggrs_per_node > 0 && ncp->rank == comm_attr.my_aggr) {
        /* comm_attr.ids[] is no longer needed. Note it has been duplicated
         * above from the MPI communicator's cached keyval attribute when
         * ncp->num_aggrs_per_node > 0.
         */
        NCI_Free(ncp->comm_attr.ids);
        ncp->comm_attr.ids = NULL;
    }
    if (ncp->pncio_fh != NULL)
        ncp->pncio_fh->comm_attr.ids = NULL;

    /* read header from file into NC object pointed by ncp ------------------*/
    err = ncmpio_hdr_get_NC(ncp);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) { /* fatal error */
        ncmpio_file_close(ncp);
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

