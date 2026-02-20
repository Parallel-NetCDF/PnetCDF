/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>      /* open(), O_CREAT */
#include <sys/types.h>  /* open(), umask(), fstat() */

#if defined(HAVE_SYS_STAT_H) && HAVE_SYS_STAT_H == 1
#include <sys/stat.h>   /* fstat() */
#endif
#include <unistd.h>     /* fstat() */

#include <assert.h>
#include <sys/errno.h>

#include <mpi.h>

#include "pncio.h"

/*----< GEN_set_cb_node_list() >---------------------------------------------*/
/* Construct the list of I/O aggregators. It sets the followings.
 *   fd->hints->ranklist[].
 *   fd->hints->cb_nodes and set file info for hint cb_nodes.
 *   fd->is_agg: indicating whether this rank is an I/O aggregator
 *   fd->my_cb_nodes_index: index into fd->hints->ranklist[]. -1 if N/A
 */
static
int GEN_set_cb_node_list(PNCIO_File *fd)
{
    int i, j, k, nprocs, rank, *nprocs_per_node, **ranks_per_node;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &rank);

    if (fd->hints->cb_nodes == 0)
        /* If hint cb_nodes is not set by user, select one rank per node to be
         * an I/O aggregator
         */
        fd->hints->cb_nodes = fd->node_ids.num_nodes;
    else if (fd->hints->cb_nodes > nprocs)
        /* cb_nodes must be <= nprocs */
        fd->hints->cb_nodes = nprocs;

    fd->hints->ranklist = (int *) NCI_Malloc(sizeof(int) * fd->hints->cb_nodes);
    if (fd->hints->ranklist == NULL)
        return NC_ENOMEM;

    /* number of MPI processes running on each node */
    nprocs_per_node = (int *) NCI_Calloc(fd->node_ids.num_nodes, sizeof(int));

    for (i=0; i<nprocs; i++) nprocs_per_node[fd->node_ids.ids[i]]++;

    /* construct rank IDs of MPI processes running on each node */
    ranks_per_node = (int **) NCI_Malloc(sizeof(int*) * fd->node_ids.num_nodes);
    ranks_per_node[0] = (int *) NCI_Malloc(sizeof(int) * nprocs);
    for (i=1; i<fd->node_ids.num_nodes; i++)
        ranks_per_node[i] = ranks_per_node[i - 1] + nprocs_per_node[i - 1];

    for (i=0; i<fd->node_ids.num_nodes; i++) nprocs_per_node[i] = 0;

    /* Populate ranks_per_node[], list of MPI ranks running on each node.
     * Populate nprocs_per_node[], number of MPI processes on each node.
     */
    for (i=0; i<nprocs; i++) {
        k = fd->node_ids.ids[i];
        ranks_per_node[k][nprocs_per_node[k]] = i;
        nprocs_per_node[k]++;
    }

    /* select process ranks from nodes in a round-robin fashion to be I/O
     * aggregators
     */
    k = j = 0;
    for (i=0; i<fd->hints->cb_nodes; i++) {
        if (j >= nprocs_per_node[k]) { /* if run out of ranks in this node k */
            k++;
            if (k == fd->node_ids.num_nodes) { /* round-robin to first node */
                k = 0;
                j++;
            }
        }
        /* select jth rank of node k as an I/O aggregator */
        fd->hints->ranklist[i] = ranks_per_node[k++][j];
        if (rank == fd->hints->ranklist[i]) {
            fd->is_agg = 1;
            fd->my_cb_nodes_index = i;
        }
        if (k == fd->node_ids.num_nodes) { /* round-robin to first node */
            k = 0;
            j++;
        }
    }
    NCI_Free(ranks_per_node[0]);
    NCI_Free(ranks_per_node);
    NCI_Free(nprocs_per_node);

    return 0;
}

/*----< GEN_create() >-------------------------------------------------------*/
/*   1. root creates the file
 *   2. root sets and obtains striping info
 *   3. root broadcasts striping info
 *   4. non-root processes receive striping info from root
 *   5. non-root processes opens the file
 */
static int
GEN_create(PNCIO_File *fd,
           int         mpi_io_mode)
{
    int err=NC_NOERR, rank, amode, perm, old_mask;
    int stripin_info[4] = {-1, -1, -1, -1};

    MPI_Comm_rank(fd->comm, &rank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
if (world_rank == 0) { printf("\nxxxx %s at %d: ---- %s\n",__func__,__LINE__,fd->filename); fflush(stdout);}
#endif

    amode = O_CREAT;
    if (mpi_io_mode & MPI_MODE_RDWR) amode |= O_RDWR;

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ PNCIO_PERM;

    /* root process creates the file first, followed by all processes open the
     * file.
     */
    if (rank > 0) goto err_out;

    fd->fd_sys = open(fd->filename, amode, perm);
    if (fd->fd_sys == -1) {
        fprintf(stderr,"%s line %d: rank %d fails to create file %s (%s)\n",
                __func__,__LINE__, rank, fd->filename, strerror(errno));
        err = ncmpii_error_posix2nc("open");
        goto err_out;
    }
    fd->is_open = 1;
    stripin_info[0] = 1048576; /* default to 1 MiB */

    /* Only root obtains the striping information and bcast to all other
     * processes. For UFS, file striping is the file system block size.
     */
#if defined(HAVE_SYS_STAT_H) && HAVE_SYS_STAT_H == 1
    struct stat statbuf;
    err = fstat(fd->fd_sys, &statbuf);
    if (err >= 0)
        /* file system block size usually < MAX_INT */
        stripin_info[0] = (int)statbuf.st_blksize;
#endif

err_out:
    MPI_Bcast(stripin_info, 4, MPI_INT, 0, fd->comm);

    fd->hints->striping_unit   = stripin_info[0];
    fd->hints->striping_factor = stripin_info[1];
    fd->hints->start_iodevice  = stripin_info[2];

    if (rank > 0) { /* non-root processes */
        fd->fd_sys = open(fd->filename, O_RDWR, perm);
        if (fd->fd_sys == -1) {
            fprintf(stderr,"%s line %d: rank %d failure to open file %s (%s)\n",
                    __func__,__LINE__, rank, fd->filename, strerror(errno));
            return ncmpii_error_posix2nc("ioctl");
        }
        fd->is_open = 1;
    }

    /* construct cb_nodes rank list */
    GEN_set_cb_node_list(fd);
    MPI_Info_set(fd->info, "romio_filesystem_type", "UFS:");

    return err;
}

/*----< GEN_open() >---------------------------------------------------------*/
/*   1. all processes open the file.
 *   2. root obtains striping info and broadcasts to all others
 */
static int
GEN_open(PNCIO_File *fd)
{
    int err=NC_NOERR, rank, perm, old_mask, omode;
    int stripin_info[4] = {1048576, -1, -1, -1};

    MPI_Comm_rank(fd->comm, &rank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
if (world_rank == 0) { printf("\nxxxx %s at %d: ---- %s\n",__func__,__LINE__,fd->filename); fflush(stdout);}
#endif

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ PNCIO_PERM;

    if (fIsSet(fd->access_mode, MPI_MODE_RDWR))
        omode = O_RDWR;
    else
        omode = O_RDONLY;

    /* All processes open the file. */
    fd->fd_sys = open(fd->filename, omode, perm);
    if (fd->fd_sys == -1) {
        fprintf(stderr, "%s line %d: rank %d failure to open file %s (%s)\n",
                __func__,__LINE__, rank, fd->filename, strerror(errno));
        err = ncmpii_error_posix2nc("open");
        goto err_out;
    }
    fd->is_open = 1;
    stripin_info[0] = 1048576; /* default to 1 MiB */

    /* Only root obtains the striping information and bcast to all other
     * processes. For UFS, file striping is the file system block size.
     */
#if defined(HAVE_SYS_STAT_H) && HAVE_SYS_STAT_H == 1
    if (rank == 0) {
        /* Get the underlying file system block size as file striping_unit */
        struct stat statbuf;
        err = fstat(fd->fd_sys, &statbuf);
        if (err >= 0)
            /* file system block size usually < MAX_INT */
            stripin_info[0] = (int)statbuf.st_blksize;
    }
#endif

err_out:
    MPI_Bcast(stripin_info, 4, MPI_INT, 0, fd->comm);
    fd->hints->striping_unit   = stripin_info[0];
    fd->hints->striping_factor = stripin_info[1];
    fd->hints->start_iodevice  = stripin_info[2];

    /* construct cb_nodes rank list */
    GEN_set_cb_node_list(fd);
    MPI_Info_set(fd->info, "romio_filesystem_type", "UFS:");

    return err;
}

/*----< PNCIO_File_open() >---------------------------------------------------*/
int PNCIO_File_open(MPI_Comm    comm,
                    const char *filename,
                    int         amode,
                    MPI_Info    info,
                    PNCIO_File *fd)
{
    /* Before reaching to this subroutine, PNCIO_FileSysType() should have been
     * called to check the file system type.
     */
    char value[MPI_MAX_INFO_VAL + 1], int_str[16];
    int i, err, min_err, status=NC_NOERR;

    fd->comm        = comm;
    fd->filename    = filename;  /* without file system type name prefix */
    fd->atomicity   = 0;
    fd->filetype    = MPI_BYTE;
    fd->is_open     = 0;
    fd->access_mode = amode;
    fd->io_buf      = NULL; /* collective buffer used by aggregators only */

    fd->flat_file.count = 0; /* flattened fileview in offset-length pairs */
    fd->flat_file.size = -1;
    fd->flat_file.is_contig = 1;
    fd->flat_file.off = NULL;
    fd->flat_file.len = NULL;

    /* create and initialize info object */
    fd->hints = (PNCIO_Hints*) NCI_Calloc(1, sizeof(PNCIO_Hints));
    status = PNCIO_File_SetInfo(fd, info);
    if (status != NC_NOERR && status != NC_EMULTIDEFINE_HINTS) {
        /* Inconsistent I/O hints is not a fatal error.
         * In PNCIO_File_SetInfo(), root's hints overwrite local's.
         */
        goto err_out;
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    for (i=0; i<NMEASURES; i++) {
        fd->write_timing[i]  = fd->read_timing[i]  = 0;
        fd->write_counter[i] = fd->read_counter[i] = 0;
    }
#endif

    assert(fd->file_system != PNCIO_FSTYPE_MPIIO);

    /* TODO: When hint romio_no_indep_rw hint is set to true, only aggregators open
     * the file.
     * Note because fd->is_agg is set at the end of create/open call.
     */
    if (fd->file_system == PNCIO_LUSTRE) {
        if (amode & MPI_MODE_CREATE)
            err = PNCIO_Lustre_create(fd, amode);
        else
            err = PNCIO_Lustre_open(fd);
    }
    else {
        if (amode & MPI_MODE_CREATE)
            err = GEN_create(fd, amode);
        else
            err = GEN_open(fd);
    }
    if (err != NC_NOERR) { /* fatal error */
        status = err;
        goto err_out;
    }

    /* set file striping hints */
    snprintf(int_str, 16, "%d", fd->hints->striping_unit);
    MPI_Info_set(fd->info, "striping_unit", int_str);

    snprintf(int_str, 16, "%d", fd->hints->striping_factor);
    MPI_Info_set(fd->info, "striping_factor", int_str);

    snprintf(int_str, 16, "%d", fd->hints->start_iodevice);
    MPI_Info_set(fd->info, "start_iodevice", int_str);

    /* set file striping hints */
    snprintf(int_str, 16, "%d", fd->hints->cb_nodes);
    MPI_Info_set(fd->info, "cb_nodes", int_str);

    /* add hint "cb_node_list", list of aggregators' rank IDs */
    snprintf(value, 16, "%d", fd->hints->ranklist[0]);
    for (i=1; i<fd->hints->cb_nodes; i++) {
        snprintf(int_str, 16, " %d", fd->hints->ranklist[i]);
        if (strlen(value) + strlen(int_str) >= MPI_MAX_INFO_VAL-5) {
            strcat(value, " ...");
            break;
        }
        strcat(value, int_str);
    }
    MPI_Info_set(fd->info, "cb_node_list", value);

    /* collective buffer size must be at least file striping size */
    if (fd->hints->cb_buffer_size < fd->hints->striping_unit) {
        fd->hints->cb_buffer_size = fd->hints->striping_unit;
        snprintf(int_str, 16, " %d", fd->hints->cb_buffer_size);
        MPI_Info_set(fd->info, "cb_buffer_size", int_str);
    }

    /* collective buffer is used only by I/O aggregators only */
    if (fd->is_agg) {
        fd->io_buf = NCI_Calloc(1, fd->hints->cb_buffer_size);
        if (fd->io_buf == NULL) /* fatal error */
            status = NC_ENOMEM;
    }

err_out:
    MPI_Allreduce(&status, &min_err, 1, MPI_INT, MPI_MIN, comm);
    /* All NC errors are < 0 */

    if (min_err != NC_NOERR) {
        if (status == NC_NOERR && fd->is_open)
            /* close file if opened successfully */
            close(fd->fd_sys);
        NCI_Free(fd->hints);
        if (fd->info != MPI_INFO_NULL)
            MPI_Info_free(&(fd->info));
        if (fd->io_buf != NULL)
            NCI_Free(fd->io_buf);
    }
    return status;
}

