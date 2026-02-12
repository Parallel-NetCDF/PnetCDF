/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_create() : dispatcher->create()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>    /* strcpy(), strchr() */

#if defined(HAVE_LSTAT) || defined(HAVE_ACCESS) || defined(HAVE_OPEN) || defined(HAVE_UNLINK) || defined(HAVE_CLOSE)
#include <sys/types.h> /* lstat(), open() */
#include <sys/stat.h>  /* lstat(), open() */
#include <unistd.h>    /* lstat(), access(), unlink(), open(), close() */
#include <fcntl.h>     /* open() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< ncmpio_create() >----------------------------------------------------*/
int
ncmpio_create(MPI_Comm         comm,
              const char      *path,
              int              cmode,
              int              ncid,
              int              env_mode,
              MPI_Info         user_info, /* user's and env info combined */
              PNCIO_node_ids   node_ids,  /* node IDs of all processes */
              void           **ncpp)      /* OUT */
{
    char *filename, value[MPI_MAX_INFO_VAL + 1], *mpi_name;
    int rank, nprocs, mpiomode, err, mpireturn, default_format, file_exist=1;
    int use_trunc=1, flag, striping_unit;
    MPI_File fh=MPI_FILE_NULL;
    NC *ncp=NULL;

    *ncpp = NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* Note path's validity and cmode consistency have been checked in
     * ncmpi_create() in src/dispatchers/file.c and path consistency will be
     * done in MPI_File_open.
     */

    /* First, check whether cmode is valid or supported ---------------------*/

    /* NC_DISKLESS is not supported yet */
    if (cmode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* NC_MMAP is not supported yet */
    if (cmode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* Check cmode for other illegal flags already done in dispatcher layer */

    /* Get default format, in case cmode does not include either
     * NC_64BIT_OFFSET or NC_64BIT_DATA.
     */
    ncmpi_inq_default_format(&default_format);

    /* allocate buffer for header object NC and initialize its contents */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    *ncpp = (void*)ncp;

    ncp->ncid     = ncid;
    ncp->comm     = comm;     /* reuse comm duplicated in dispatch layer */
    ncp->rank     = rank;
    ncp->nprocs   = nprocs;

    /* Extract hints from user_info. Two hints must be extracted now in order
     * to continue:
     *     nc_pncio: whether to user MPI-IO or PnetCDF's PNCIO driver.
     *     nc_num_aggrs_per_node: number of processes per node to be the INA
     *     aggregators.
     *
     * ncp->fstype will be initialized in ncmpio_hint_extract(), and set in
     * PNCIO_FileSysType().
     */
    ncmpio_hint_extract(ncp, user_info);

    if (ncp->fstype == PNCIO_FSTYPE_CHECK && rank == 0)
        /* Check file system type. If the given file does not exist, check its
         * folder. Currently PnetCDF's PNCIO drivers support Lustre
         * (PNCIO_LUSTRE) and Unix File System (PNCIO_UFS).
         */
        ncp->fstype = PNCIO_FileSysType(path);

#ifdef WKL_DEBUG
if (rank == 0) printf("%s at %d fstype=%s\n", __func__,__LINE__,(ncp->fstype == PNCIO_FSTYPE_MPIIO)? "PNCIO_FSTYPE_MPIIO" : (ncp->fstype == PNCIO_LUSTRE) ? "PNCIO_LUSTRE" : "PNCIO_UFS");
#endif

    /* Setting file open mode in mpiomode which may later be needed in
     * ncmpi_begin_indep_data() to open file for independent data mode.
     */
    mpiomode = MPI_MODE_RDWR | MPI_MODE_CREATE;

    /* Remove the file system type prefix name if there is any. For example,
     * when path = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * filename pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * access() below.
     */
    filename = ncmpii_remove_file_system_type_prefix(path);

    /* In case of clobber mode, first check if the file already exists, through
     * a call to lstat() or access() if they are is available. If not, we
     * assume the file exists and will add some MPI flag to open mode argument
     * of MPI_File_open to either delete or truncate the file first.
     */
#ifdef HAVE_LSTAT
    /* Call lstat() to check the file if exists and if is a symbolic link */
    if (rank == 0) {
        struct stat st_buf;
        st_buf.st_mode = 0;

        if (lstat(filename, &st_buf) == -1) file_exist = 0;
        errno = 0; /* reset errno */

        /* If the file is a regular file, not a symbolic link, then we delete
         * the file first and later create it when calling MPI_File_open() with
         * MPI_MODE_CREATE. If the file is a regular file, not a symbolic link,
         * it is faster to delete it and then re-create the file, as truncating
         * it to zero size is more expensive.
         *
         * If the file is a symbolic link, then we cannot delete the file, as
         * the link will be gone. If the file is deleted and there are other
         * files symbolically linked to this file, then their links will become
         * invalid.
         */
        if (S_ISREG(st_buf.st_mode)) use_trunc = 0;
    }
#elif defined HAVE_ACCESS
    /* If access() is available, use it to check whether file already exists,
     * by having rank 0 to call access() and broadcast file_exist.
     */
    if (rank == 0) {
        if (access(filename, F_OK) == -1) file_exist = 0;
        errno = 0; /* reset errno */
    }
#endif

    if (fIsSet(cmode, NC_NOCLOBBER)) {
        /* Error NC_EEXIST will be returned, if the file already exists and
         * NC_NOCLOBBER mode is set in ncmpi_create.
         */
#ifdef HAVE_ACCESS
        if (nprocs > 1) {
            int msg[2] = {file_exist, ncp->fstype};
            TRACE_COMM(MPI_Bcast)(msg, 2, MPI_INT, 0, comm);
            file_exist  = msg[0];
            ncp->fstype = msg[1];
        }
        if (file_exist) {
            NCI_Free(ncp);
            DEBUG_RETURN_ERROR(NC_EEXIST)
        }
#else
        if (nprocs > 1)
            TRACE_COMM(MPI_Bcast)(&ncp->fstype, 1, MPI_INT, 0, comm);

        /* Add MPI_MODE_EXCL mode for MPI_File_open, so it can error out, if
         * the file exists.
         */
        fSet(mpiomode, MPI_MODE_EXCL);
        errno = 0; /* reset errno, as MPI_File_open may change it */
#endif
    }
    else {
        /* NC_CLOBBER is the default mode in ncmpi_create(). Below, rank 0
         * truncates or deletes the file and ignores error code.  Note in some
         * implementation of MPI-IO, calling MPI_File_set_size is expensive as
         * it may call truncate() by all ranks.
         */
        err = NC_NOERR;
        if (rank == 0 && file_exist) {
            if (!use_trunc) { /* delete the file */
#ifdef HAVE_UNLINK
                /* unlink() is likely faster then truncate(). However, unlink()
                 * can be expensive when the file size is large. For example,
                 * it taook 1.1061 seconds to delete a file of size 27.72 GiB
                 * on Perlmutter at NERSC.
                 */
                err = unlink(filename);
                if (err < 0 && errno != ENOENT)
                    /* ignore ENOENT: file not exist */
                    DEBUG_ASSIGN_ERROR(err, NC_EFILE) /* report other error */
                else
                    err = NC_NOERR;
#else
                err = NC_NOERR;
                if (ncp->fstype != PNCIO_FSTYPE_MPIIO)
                    err = PNCIO_File_delete(filename);
                else {
                    TRACE_IO(MPI_File_delete, (path, MPI_INFO_NULL));
                    if (mpireturn != MPI_SUCCESS) {
                        int errorclass;
                        MPI_Error_class(mpireturn, &errorclass);
                        if (errorclass != MPI_ERR_NO_SUCH_FILE)
                            /* ignore file not exist */
                            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    }
                }
#endif
            }
            else {
                /* If file is not a regular file (e.g. a symbolic link), we
                 * cannot delete it and must truncate it to zero size. In this
                 * case, file open mode needs to remove MPI_MODE_CREATE.
                 */
                mpiomode = MPI_MODE_RDWR;

#ifdef HAVE_TRUNCATE
                err = truncate(filename, 0); /* This may be expensive */
                if (err < 0 && errno != ENOENT)
                    /* ignore ENOENT: file not exist */
                    DEBUG_ASSIGN_ERROR(err, NC_EFILE) /* report other error */
                else
                    err = NC_NOERR;
#elif defined HAVE_OPEN
                int fd = open(filename, O_TRUNC | O_WRONLY);
                if (fd < 0)
                    DEBUG_ASSIGN_ERROR(err, NC_EFILE)
                else {
                    err = close(fd);
                    if (err < 0)
                        DEBUG_ASSIGN_ERROR(err, NC_EFILE)
                }
#else
                /* When all POSIX system calls are not available, the last
                 * resort is to call MPI_File_set_size() to truncate the file.
                 * Note for some ROMIO versions that have all processes call
                 * truncate(), this option can be expensive.
                 */
                err = NC_NOERR;
                if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
                    PNCIO_File pncio_fh;
                    pncio_fh = (PNCIO_File*) NCI_Calloc(1,sizeof(PNCIO_File));
                    err = PNCIO_File_open(MPI_COMM_SELF, filename,
                                          MPI_MODE_RDWR, MPI_INFO_NULL,
                                          pncio_fh);
                    if (err == NC_NOERR)
                        PNCIO_File_set_size(pncio_fh, 0); /* can be expensive */
                    else
                        PNCIO_File_close(&pncio_fh);
                    NCI_Free(pncio_fh);
                }
                else {
                    TRACE_IO(MPI_File_open, (MPI_COMM_SELF, path, MPI_MODE_RDWR, MPI_INFO_NULL, &fh));
                    if (mpireturn != MPI_SUCCESS) {
                        int errorclass;
                        MPI_Error_class(mpireturn, &errorclass);
                        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    }
                    else {
                        TRACE_IO(MPI_File_set_size, (fh, 0)); /* can be expensive */
                        if (mpireturn != MPI_SUCCESS) {
                            int errorclass;
                            MPI_Error_class(mpireturn, &errorclass);
                            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                        }
                        else {
                            TRACE_IO(MPI_File_close, (&fh));
                            if (mpireturn != MPI_SUCCESS) {
                                int errorclass;
                                MPI_Error_class(mpireturn, &errorclass);
                                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                            }
                        }
                    }
                }
#endif
            }
            if (errno == ENOENT) errno = 0; /* reset errno */
        }
        /* All processes must wait here until clobbering file by root process
         * is completed. Note mpiomode may be changed to remove MPI_MODE_CREATE
         * when the file to be clobbered is a symbolic link.
         */
        if (nprocs > 1) {
            int msg[3] = {err, mpiomode, ncp->fstype};
            TRACE_COMM(MPI_Bcast)(&msg, 3, MPI_INT, 0, comm);
            err         = msg[0];
            mpiomode    = msg[1];
            ncp->fstype = msg[2];
        }
        if (err != NC_NOERR) return err;
    }
    /* Now file has been clobbered, i.e. deleted if it is not a symbolic link.
     * If it is a symbolic link, it now has been truncated to zero size.
     */

    ncp->path     = path;     /* reuse path duplicated in dispatch layer */
    ncp->pncio_fh = NULL;     /* non-aggregators have NULL pncio_fh */
    ncp->mpiomode = mpiomode;
    ncp->mpiinfo  = MPI_INFO_NULL;

    /* For file create, ignore NC_NOWRITE if set in cmode argument. */
    ncp->iomode   = cmode | NC_WRITE;

    ncp->collective_fh  = MPI_FILE_NULL;
    ncp->independent_fh = MPI_FILE_NULL;

    /* set the file format version based on the create mode, cmode */
         if (fIsSet(cmode, NC_64BIT_DATA))   ncp->format = 5;
    else if (fIsSet(cmode, NC_64BIT_OFFSET)) ncp->format = 2;
    else {
             if (default_format == NC_FORMAT_CDF5) ncp->format = 5;
        else if (default_format == NC_FORMAT_CDF2) ncp->format = 2;
        else                                       ncp->format = 1;
    }

    /* indicate this is from ncmpi_create */
    fSet(ncp->flags, NC_MODE_CREATE);
    /* create automatically enter write mode */
    fClr(ncp->flags, NC_MODE_RDONLY);
    /* create automatically enter define mode */
    fSet(ncp->flags, NC_MODE_DEF);
    /* PnetCDF default mode is no fill */
    fClr(ncp->flags, NC_MODE_FILL);

    /* incorporate modes set in environment variables */
    fSet(ncp->flags, env_mode);

    /* initialize unlimited_id as no unlimited dimension yet defined */
    ncp->dims.unlimited_id = -1;

    /* node_ids stores a list of unique IDs of compute nodes of all MPI ranks
     * in the MPI communicator passed from the user application. It is a keyval
     * attribute cached in the communicator. See src/dispatchers/file.c for
     * details. The node IDs will be used when the intra-node aggregation (INA)
     * is enabled and when PnetCDF's PNCIO driver is used.
     *
     * When intra-node aggregation (INA) is enabled, node IDs are used to
     * create a new MPI communicator consisting of the intra-node aggregators
     * only. The communicator will be used to call file open in MPI-IO or
     * PnetCDF's PNCIO driver. This means only intra-node aggregators will
     * perform file I/O in PnetCDF collective put and get operations.
     *
     * node_ids will be used to calculate cb_nodes, the number of MPI-IO/PNCIO
     * aggregators (not INA aggregators).
     */
    ncp->node_ids = node_ids;

    /* When the total number of aggregators >= number of processes, disable
     * intra-node aggregation.
     */
    if (ncp->num_aggrs_per_node * node_ids.num_nodes >= ncp->nprocs)
        ncp->num_aggrs_per_node = 0;

    /* ncp->num_aggrs_per_node = 0, or > 0 is an indicator of whether the INA
     * feature is disabled or enabled globally for all processes.
     */
    ncp->my_aggr = -1;
    ncp->ina_comm = MPI_COMM_NULL;
    ncp->ina_nprocs = 0;
    ncp->ina_rank = -1;
    ncp->ina_node_list = NULL;
    if (ncp->num_aggrs_per_node > 0) {
        /* Must duplicate node_ids, as node_ids.ids[] will be modified by
         * ncmpio_ina_init().
         */
        ncp->node_ids.ids = (int*) NCI_Malloc(sizeof(int) * ncp->nprocs);
        memcpy(ncp->node_ids.ids, node_ids.ids, sizeof(int) * ncp->nprocs);

        /* Divide all ranks into groups. Each group is assigned one intra-node
         * aggregator. The following metadata related to intra-node aggregation
         * will be set up in ncmpio_ina_init().
         * ncp->my_aggr is the aggregator's rank ID (related to ncp->comm) of
         *     this group. When == ncp->rank, this rank is an aggregator.
         * ncp->num_nonaggrs is the number of non-aggregators assigned to this
         *     rank (an aggregator)
         * ncp->ina_comm is an MPI communicator consisting of only intra-node
         *     aggregators across all nodes, which will be used when calling
         *     MPI_File_open(). For non-aggregator, it == MPI_COMM_NULL.
         * ncp->node_ids.ids[] will be modified to contain the nodes IDs of all
         *     intra-node aggregators, and will be passed to pncio_fh.
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

    /* create file collectively -------------------------------------------- */
    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        /* If hint nc_striping is set to "auto" and hint striping_factor is not
         * set by the user, then set hint striping_factor to ncp->num_nodes.
         */
        if (ncp->nc_striping == PNCIO_STRIPING_AUTO) {
            int striping_factor=0;
            if (user_info != MPI_INFO_NULL) {
                MPI_Info_get(user_info, "striping_factor", MPI_MAX_INFO_VAL-1,
                            value, &flag);
                if (flag)
                    striping_factor = atoi(value);
            }
            if (striping_factor == 0) {
                sprintf(value, "%d", ncp->node_ids.num_nodes);
                MPI_Info_set(user_info, "striping_factor", value);
            }
        }

        TRACE_IO(MPI_File_open, (comm, path, mpiomode, user_info, &fh));
        if (mpireturn != MPI_SUCCESS) {
#ifndef HAVE_ACCESS
            if (fIsSet(cmode, NC_NOCLOBBER)) {
                /* This is the case when NC_NOCLOBBER is used in file creation
                 * and function access() is not available. MPI_MODE_EXCL is set
                 * in open mode. When MPI_MODE_EXCL is used and the file
                 * already exists, MPI-IO should return error class
                 * MPI_ERR_FILE_EXISTS. But, some MPI-IO implementations (older
                 * ROMIO) do not correctly return this error class. In this
                 * case, we can do the followings: check errno to see if it set
                 * to EEXIST. Note usually rank 0 makes the file open call and
                 * can be the only one having errno set.
                 */
                if (nprocs > 1)
                    TRACE_COMM(MPI_Bcast)(&errno, 1, MPI_INT, 0, comm);
                if (errno == EEXIST) {
                    NCI_Free(ncp);
                    DEBUG_FOPEN_ERROR(NC_EEXIST)
                }
            }
#endif
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_File_open");
            DEBUG_FOPEN_ERROR(err);
            /* for NC_NOCLOBBER, MPI_MODE_EXCL was added to mpiomode. If the
             * file already exists, MPI-IO should return error class
             * MPI_ERR_FILE_EXISTS which PnetCDF will return error code
             * NC_EEXIST. This is checked inside of ncmpii_error_mpi2nc()
             */
        }
        else
            /* reset errno, as MPI_File_open may change it, even if it returns
             * MPI_SUCCESS
             */
            errno = 0;

        /* Now the file has been successfully created */
        ncp->collective_fh  = fh;
        ncp->independent_fh = (nprocs == 1) ? fh : MPI_FILE_NULL;

        /* get the I/O hints used/modified by MPI-IO */
        TRACE_IO(MPI_File_get_info, (fh, &ncp->mpiinfo));
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            DEBUG_FOPEN_ERROR(err);
        }
    }
    else {
        /* When ncp->fstype != PNCIO_FSTYPE_MPIIO, use PnetCDF's PNCIO driver */
        ncp->pncio_fh = (PNCIO_File*) NCI_Calloc(1, sizeof(PNCIO_File));
        ncp->pncio_fh->file_system = ncp->fstype;
        ncp->pncio_fh->node_ids    = ncp->node_ids;

        err = PNCIO_File_open(comm, filename, mpiomode, user_info,
                              ncp->pncio_fh);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err)

        /* Now the file has been successfully created, obtain the I/O hints
         * used/modified by PNCIO driver.
         */
        err = PNCIO_File_get_info(ncp->pncio_fh, &ncp->mpiinfo);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err)
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

/*
if (ncp->rank == 0) {
    int  i, nkeys;
    MPI_Info_get_nkeys(ncp->mpiinfo, &nkeys);
    printf("%s line %d: MPI File Info: nkeys = %d\n",__func__,__LINE__,nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY];
        int  valuelen;
        MPI_Info_get_nthkey(ncp->mpiinfo, i, key);
        MPI_Info_get_valuelen(ncp->mpiinfo, key, &valuelen, &flag);
        MPI_Info_get(ncp->mpiinfo, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %25s, value = %s\n",i,key,value);
    }
}
*/

    /* ina_node_list is no longer needed */
    if (ncp->ina_node_list != NULL) {
        NCI_Free(ncp->ina_node_list);
        ncp->ina_node_list = NULL;
    }
    if (ncp->num_aggrs_per_node > 0) {
        /* node_ids is no longer needed. Note node_ids is duplicated above from
         * the MPI communicator's cached keyval attribute when
         * ncp->num_aggrs_per_node > 0.
         */
        NCI_Free(ncp->node_ids.ids);
        ncp->node_ids.ids = NULL;
    }
    if (ncp->pncio_fh != NULL)
        ncp->pncio_fh->node_ids.ids = NULL;

    return NC_NOERR;
}

