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
#include <fcntl.h>     /* open(), O_CREAT */
#endif
#include <libgen.h>    /* dirname() */
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"


#ifdef HAVE_LUSTRE
/* /usr/include/lustre/lustreapi.h
 * /usr/include/linux/lustre/lustre_user.h
 */
#include <lustre/lustreapi.h>

// #define PNETCDF_LUSTRE_DEBUG

#define PATTERN_STR(pattern, int_str) ( \
    (pattern == LLAPI_LAYOUT_DEFAULT)      ? "LLAPI_LAYOUT_DEFAULT" : \
    (pattern == LLAPI_LAYOUT_RAID0)        ? "LLAPI_LAYOUT_RAID0" : \
    (pattern == LLAPI_LAYOUT_WIDE)         ? "LLAPI_LAYOUT_WIDE" : \
    (pattern == LLAPI_LAYOUT_MDT)          ? "LLAPI_LAYOUT_MDT" : \
    (pattern == LLAPI_LAYOUT_OVERSTRIPING) ? "LLAPI_LAYOUT_OVERSTRIPING" : \
    (pattern == LLAPI_LAYOUT_SPECIFIC)     ? "LLAPI_LAYOUT_SPECIFIC" : \
    int_str)

/*----< lustre_get_striping() >----------------------------------------------*/
static
void lustre_get_striping(const char *path,
                         uint64_t   *stripe_count,
                         uint64_t   *stripe_size)
{
    char *dirc, *dname;
    int err, fd;
    struct llapi_layout *layout;

    /* Retrieve the file striping settings from the parent folder. */
    dirc = NCI_Strdup(path);
    dname = dirname(dirc); /* folder name */

    fd = open(dname, O_RDONLY, PNCIO_PERM);

    layout = llapi_layout_get_by_fd(fd, LLAPI_LAYOUT_GET_COPY);
    if (layout == NULL) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_layout_get_by_fd() fails\n",
                __FILE__, __LINE__);
#endif
        goto err_out;
    }

    if (stripe_count != NULL && *stripe_count == 0) {
        /* obtain file striping count */
        err = llapi_layout_stripe_count_get(layout, stripe_count);
        if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
            char int_str[32];
            snprintf(int_str, 32, "%lu", *stripe_count);
            printf("Error at %s (%d) llapi_layout_stripe_count_get() fails to get stripe count %s\n",
                __FILE__, __LINE__, PATTERN_STR(*stripe_count, int_str));
#endif
            goto err_out;
        }
    }

    if (stripe_size != NULL && *stripe_size == 0) {
        /* obtain file striping unit size */
        err = llapi_layout_stripe_size_get(layout, stripe_size);
        if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
            char int_str[32];
            snprintf(int_str, 32, "%lu", *stripe_size);
            printf("Error at %s (%d) llapi_layout_stripe_size_get() fails to get stripe size %s\n",
                __FILE__,__LINE__, PATTERN_STR(*stripe_size, int_str));
#endif
            goto err_out;
        }
    }

err_out:
    if (layout != NULL) llapi_layout_free(layout);

    close(fd);
}
#endif

/*----< ncmpio_create() >----------------------------------------------------*/
int
ncmpio_create(MPI_Comm         comm,
              const char      *path,
              int              cmode,
              int              ncid,
              int              env_mode,
              MPI_Info         user_info, /* user's and env info combined */
              PNC_comm_attr    comm_attr, /* node IDs and INA metadata */
              void           **ncpp)      /* OUT */
{
    char *filename, value[MPI_MAX_INFO_VAL + 1], *mpi_name;
    int rank, nprocs, mpi_amode, err, mpireturn, default_format, file_exist=1;
    int use_trunc=1, flag, striping_info[2];
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

    ncp->ncid   = ncid;
    ncp->comm   = comm;     /* reuse comm duplicated in dispatch layer */
    ncp->rank   = rank;
    ncp->nprocs = nprocs;

    /* Extract hints from user_info.
     *
     *   Three hints must be extracted before calling PNCIO_File_open() or
     *   MPI_File_open().
     *      + nc_driver: whether to use MPI-IO or PnetCDF's PNCIO driver.
     *      + nc_num_aggrs_per_node: number of processes per node to be the INA
     *          aggregators, which determines the MPI communicator to be passed
     *          to PNCIO_File_open() and MPI_File_open().
     *      + nc_striping: whether to inherit parent folder's striping.
     *
     *   Four hints must be extracted before ncmpio_enddef(). They are used
     *   when creating dimensions, attributes, and variables.
     *      + nc_hash_size_dim: hash table size for dimensions
     *      + nc_hash_size_var: hash table size for variables
     *      + nc_hash_size_gattr: hash table size for global attributes
     *      + nc_hash_size_vattr: hash table size for non-global attributes
     *
     *   The remaining PnetCDF hints are used at ncmpio_enddef() and after.
     *      + nc_var_align_size             ncp->info_v_align
     *      + nc_header_align_size          ncp->info_v_align
     *      + nc_record_align_size          ncp->info_r_align
     *      + nc_header_read_chunk_size     ncp->hdr_chunk
     *      + nc_in_place_swap              fSet(ncp->flags, NC_MODE_SWAP_ON)
     *      + nc_ibuf_size                  ncp->ibuf_size
     *      + pnetcdf_subfiling             ncp->subfile_mode/num_subfiles
     *      + nc_data_move_chunk_size       ncp->data_chunk
     *      + romio_no_indep_rw             fSet(ncp->flags, NC_HCOLL)
     *
     *   Note after a call to MPI_File_open() returns, any hints that are not
     *   used by MPI_File_open() will be discarded, which includes all PnetCDF
     *   hints. Therefore, we must at first call MPI_File_get_info() to obtain
     *   an info object containing all the hints used by MPI-IO, and then add
     *   the PnetCDF hints into the info object (ncp->mpiinfo). So ncp->mpiinfo
     *   can be used in ncmpi_inq_file_info() to return all hints to users.
     *
     *   Note MPI_File_open() may also add new hints, such as those related to
     *   file striping and I/O aggregators.
     *      + striping_unit
     *      + striping_factor
     *      + start_iodevice
     *      + cb_nodes
     *
     *   PNCIO_File_open() may also add new hints, such as those related to
     *   file striping and I/O aggregators.
     *      + lustre_overstriping_ratio
     *      + lustre_num_osts
     *      + cb_nodes
     *      + cb_node_list
     */
    ncmpio_hint_extract(ncp, user_info);

    if (rank == 0)
        /* Check file system type. If the given file does not exist, check its
         * parent folder. Currently PnetCDF's PNCIO drivers support Lustre
         * (PNCIO_FS_LUSTRE) and Unix File System (PNCIO_FS_UFS). This info
         * will also be used when MPI-IO driver is used to configure file
         * striping settings.
         */
        ncp->fstype = PNCIO_FileSysType(path);

    /* Setting file open mode in mpi_amode which may later be needed in
     * ncmpi_begin_indep_data() to open file in the independent data mode.
     */
    mpi_amode = MPI_MODE_RDWR | MPI_MODE_CREATE;

    /* Remove the file system type prefix name if there is any. For example,
     * when path = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * filename pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * access() below.
     */
    filename = ncmpii_remove_file_system_type_prefix(path);

    /* In case of file clobber mode, we first check if the file already exists,
     * through a call to lstat() or access() if they are is available. If not,
     * we call MPI_File_open() to try opening the file. If the file exists, we
     * will cal MPI_File_set_size() to truncate the file.
     */
#ifdef HAVE_LSTAT
    /* Call lstat() to check the file if exists and if is a symbolic link */
    if (rank == 0) {
        struct stat st_buf;
        st_buf.st_mode = 0;

        if (lstat(filename, &st_buf) == -1) file_exist = 0;
        errno = 0; /* reset errno */

        /* If the file is a regular file, not a symbolic link, then we delete
         * the file first and later create it (by calling MPI_File_open() with
         * MPI_MODE_CREATE). In this case, it is faster to delete and re-create
         * the file, as truncating a file to zero size is more expensive.
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
        /* Error NC_EEXIST will be returned, if the file already exists. */
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
        fSet(mpi_amode, MPI_MODE_EXCL);
        errno = 0; /* reset errno, as MPI_File_open may change it */
#endif
    }
    else {
        /* NC_CLOBBER is the default mode in ncmpi_create(). Below, rank 0
         * truncates or deletes the file and ignores error code.  Note in some
         * implementation of MPI-IO, calling MPI_File_set_size() can be very
         * expensive as it may call truncate() by all ranks.
         */
        err = NC_NOERR;
        if (rank == 0 && file_exist) {
            if (!use_trunc) { /* delete the file */
#ifdef HAVE_UNLINK
                /* unlink() is likely faster then truncate(). However, unlink()
                 * can be expensive when the file size is large. For example,
                 * it took 1.1061 seconds to delete a file of size 27.72 GiB
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
                if (ncp->driver == PNC_DRIVER_PNCIO)
                    err = PNCIO_File_delete(filename);
                else {
#ifdef MPICH_VERSION
                    /* MPICH recognizes file system type acronym prefixed to
                     * the file name. Other MPI implementations may not.
                     */
                    TRACE_IO(MPI_File_delete, (path, MPI_INFO_NULL));
#else
                    TRACE_IO(MPI_File_delete, (filename, MPI_INFO_NULL));
#endif
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
                mpi_amode = MPI_MODE_RDWR;

#ifdef HAVE_TRUNCATE
                err = truncate(filename, 0); /* truncate() may be expensive */
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
                 * Note in some ROMIO versions that make all processes to call
                 * truncate(), this option may be expensive.
                 */
                err = NC_NOERR;
                if (ncp->driver == PNC_DRIVER_PNCIO) {
                    PNCIO_File pncio_fh;
                    pncio_fh = (PNCIO_File*) NCI_Calloc(1,sizeof(PNCIO_File));
                    err = PNCIO_File_open(MPI_COMM_SELF, filename, O_RDWR,
                                          MPI_INFO_NULL, pncio_fh);
                    if (err == NC_NOERR)
                        PNCIO_File_set_size(pncio_fh, 0); /* may be expensive */
                    else
                        PNCIO_File_close(&pncio_fh);
                    NCI_Free(pncio_fh);
                }
                else {
#ifdef MPICH_VERSION
                    /* MPICH recognizes file system type acronym prefixed to
                     * the file name. Other MPI implementations may not.
                     */
                    TRACE_IO(MPI_File_open, (MPI_COMM_SELF, path,
                             MPI_MODE_RDWR, MPI_INFO_NULL, &fh));
#else
                    TRACE_IO(MPI_File_open, (MPI_COMM_SELF, filename,
                             MPI_MODE_RDWR, MPI_INFO_NULL, &fh));
#endif
                    if (mpireturn != MPI_SUCCESS) {
                        int errorclass;
                        MPI_Error_class(mpireturn, &errorclass);
                        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    }
                    else {
                        /* MPI_File_set_size() can be expensive */
                        TRACE_IO(MPI_File_set_size, (fh, 0));
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
        /* All processes must wait here until file clobbering by root process
         * is completed. Note mpi_amode may be modified by removing
         * MPI_MODE_CREATE when the file to be clobbered is a symbolic link.
         */
        if (nprocs > 1) {
            int msg[3] = {err, mpi_amode, ncp->fstype};
            TRACE_COMM(MPI_Bcast)(&msg, 3, MPI_INT, 0, comm);
            err         = msg[0];
            mpi_amode   = msg[1];
            ncp->fstype = msg[2];
        }
        if (err != NC_NOERR) return err;
    }
    /* Now file has been clobbered, i.e. deleted if it is not a symbolic link.
     * If it is a symbolic link, it now has been truncated to zero size.
     */

    ncp->path      = path;     /* reuse path duplicated in dispatch layer */
    ncp->mpi_amode = mpi_amode;
    ncp->mpiinfo   = MPI_INFO_NULL;

    if (ncp->driver == PNC_DRIVER_PNCIO) {
        /* Initialize pncio_fh, PNCIO file handler, with common metadata shared
         * among all processes, including non-INA aggregators when INA is
         * enabled. This is necessary for non-INA aggregators to perform
         * independent I/O.
         */
        ncp->pncio_fh = (PNCIO_File*) NCI_Calloc(1, sizeof(PNCIO_File));
        ncp->pncio_fh->comm           = comm;
        ncp->pncio_fh->fstype         = ncp->fstype;
        ncp->pncio_fh->comm_attr      = ncp->comm_attr;
        ncp->pncio_fh->comm_attr      = ncp->comm_attr;
        ncp->pncio_fh->file_view.size = -1;
        ncp->pncio_fh->filename       = filename;
        ncp->pncio_fh->info           = MPI_INFO_NULL;
        ncp->pncio_fh->amode          = O_CREAT|O_RDWR;
    }
    else
        ncp->pncio_fh  = NULL; /* used only when using PNCIO driver */

    /* For file create, ignore NC_NOWRITE if set in cmode argument. */
    ncp->nc_amode = cmode | NC_WRITE;

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

    /* comm_attr has been constructed at the dispatchers.
     *
     * comm_attr.ids[] stores a list of unique IDs of compute nodes of all MPI
     * ranks of the MPI communicator passed from the user application. It is a
     * keyval attribute cached in the communicator. See src/dispatchers/file.c
     * for details.
     *
     * commm_attr also stores the intra-node aggregation (INA) metadata,
     * including the INA communicator.
     *
     * When INA is enabled, comm_attr.ids[] is used to create a new MPI
     * communicator consisting of the intra-node aggregators only. The
     * communicator will be used to call PNCIO_File_open() or MPI_File_open()
     * This means only intra-node aggregators will perform file I/O in PnetCDF
     * collective put and get operations.
     *
     * In addition, comm_attr.ids[] will be used to calculate cb_nodes, the
     * number of file I/O aggregators (not INA aggregators).
     */
    ncp->comm_attr = comm_attr;

    /* When the total number of aggregators >= number of processes, disable
     * intra-node aggregation.
     */
    if (ncp->num_aggrs_per_node * comm_attr.num_nodes >= ncp->nprocs)
        ncp->num_aggrs_per_node = 0;

    /* ncp->num_aggrs_per_node = 0, or > 0 is an indicator of whether the INA
     * feature is disabled or enabled globally for all processes.
     */
    ncp->ina_nprocs = 0;
    ncp->ina_rank = -1;

    if (ncp->num_aggrs_per_node > 0) {

        if (ncp->rank == comm_attr.my_aggr) {
            /* this rank is an INA aggregator */
            int i, j, *ids;

            MPI_Comm_size(comm_attr.ina_comm, &ncp->ina_nprocs);
            MPI_Comm_rank(comm_attr.ina_comm, &ncp->ina_rank);

            /* overwrite comm_attr.ids[] by condensing it to contain only the
             * node IDs of processes in the INA communicator.
             */
            ids = (int*) NCI_Malloc(sizeof(int) * ncp->ina_nprocs);

            for (i=0; i<ncp->ina_nprocs; i++) {
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
         * comm with ina_comm. Same for nprocs.
         */
        comm = comm_attr.ina_comm;
        nprocs = ncp->ina_nprocs;

        /* For non-INA aggregators, we keep comm, a local variable in this
         * subroutine, to be comm_attr.ina_comm, a MPI_COMM_NULL for non-INA
         * aggregators. Because the remaining lines of codes below till label
         * 'fn_exit' is for INA aggregators to open the file and obtain a file
         * handler, in which non-INA aggregators do not participate and thus do
         * not make use of comm.
         */
        if (comm == MPI_COMM_NULL) {
            if (user_info != MPI_INFO_NULL)
                MPI_Info_dup(user_info, &ncp->mpiinfo);
            goto fn_exit;
        }
    }

    /* create file collectively -------------------------------------------- */
    if (ncp->driver == PNC_DRIVER_MPIIO) {
        /* If hint file_striping is set to "auto" and hint striping_factor is
         * not set by the user, then set hint striping_factor to
         * ncp->comm_attr.num_nodes.
         */
        if (ncp->file_striping == PNCIO_STRIPING_AUTO) {
            int striping_factor=0;
            if (user_info != MPI_INFO_NULL) {
                MPI_Info_get(user_info, "striping_factor", MPI_MAX_INFO_VAL-1,
                            value, &flag);
                if (flag)
                    striping_factor = atoi(value);
            }
            if (striping_factor == 0) {
                sprintf(value, "%d", ncp->comm_attr.num_nodes);
                MPI_Info_set(user_info, "striping_factor", value);
            }
        }
        else { /* ncp->file_striping == PNCIO_STRIPING_INHERIT */
            if (user_info != MPI_INFO_NULL) {
                striping_info[0] = striping_info[1] = 0;

                /* check if hint striping_factor is set by the user */
                MPI_Info_get(user_info, "striping_factor", MPI_MAX_INFO_VAL-1,
                            value, &flag);
                if (flag)
                    striping_info[0] = atoi(value);

                /* check if hint striping_unit is set by the user */
                MPI_Info_get(user_info, "striping_unit", MPI_MAX_INFO_VAL-1,
                            value, &flag);
                if (flag)
                    striping_info[1] = atoi(value);

#ifdef HAVE_LUSTRE
                uint64_t striping_factor, striping_unit;
                striping_factor = striping_info[0];
                striping_unit   = striping_info[1];
                /* When either striping_factor or striping_unit is not set, but
                 * not both, retrieve folder's striping factor or unit in order
                 * to inherit the missing one.
                 */
                if (ncp->rank == 0 &&
                    striping_factor * striping_unit == 0 &&
                    striping_factor + striping_unit > 0) {
                    /* rank 0 retrieves folder's striping settings */
                    lustre_get_striping(filename, &striping_factor,
                                        &striping_unit);
                    /* error is ignored, if there is any */
                    striping_info[0] = striping_factor;
                    striping_info[1] = striping_unit;
                }
                MPI_Bcast(striping_info, 2, MPI_INT, 0, comm);
#endif

                if (striping_info[0] > 0) {
                    sprintf(value, "%d", striping_info[0]);
                    MPI_Info_set(user_info, "striping_factor", value);
                }

                if (striping_info[1] > 0) {
                    sprintf(value, "%d", striping_info[1]);
                    MPI_Info_set(user_info, "striping_info", value);
                }
            }
        }

#ifdef MPICH_VERSION
        /* MPICH recognizes file system type acronym prefixed to file names */
        TRACE_IO(MPI_File_open, (comm, path, mpi_amode, user_info, &fh));
#else
        TRACE_IO(MPI_File_open, (comm, filename, mpi_amode, user_info, &fh));
#endif
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
            /* for NC_NOCLOBBER, MPI_MODE_EXCL was added to mpi_amode. If the
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

        /* Now the file has been successfully created, obtain the I/O hints
         * used/modified by MPI-IO.
         */
        TRACE_IO(MPI_File_get_info, (fh, &ncp->mpiinfo));
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            DEBUG_FOPEN_ERROR(err);
        }
    }
    else {
        /* When ncp->driver == PNC_DRIVER_PNCIO, use PnetCDF's PNCIO driver.
         * When INA is enabled, only the INA aggregators can reach here.
         * Non-INA aggregators have gone to fn_exit from above, with their comm
         * remain MPI_COMM_NULL (ncp->comm is always assigned to comm).
         */
        ncp->pncio_fh->fstype = ncp->fstype;
        ncp->pncio_fh->comm_attr = ncp->comm_attr;

        /* use_trunc also indicates whether the file has already existed as a
         * symbolic link, For a symbolic link file, we cannot add O_CREAT.
         */
        int amode = (mpi_amode & MPI_MODE_CREATE) ? O_CREAT|O_RDWR : O_RDWR;
        err = PNCIO_File_open(comm, filename, amode, user_info, ncp->pncio_fh);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err)

        /* Now the file has been successfully created, obtain the I/O hints
         * used/modified by PNCIO driver.
         */
        err = PNCIO_File_get_info(ncp->pncio_fh, &ncp->mpiinfo);
        if (err != NC_NOERR) DEBUG_FOPEN_ERROR(err)
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
                /* Must initialize hints for non-INA aggregators, as hints will
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

    return NC_NOERR;
}

