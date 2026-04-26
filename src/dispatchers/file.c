/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>     /* getenv() */
#include <string.h>     /* strtok(), strtok_r(), strchr(), strcpy() */
#include <strings.h>    /* strcasecmp() */
#include <fcntl.h>      /* open() */
#include <sys/types.h>  /* lseek() */
#include <unistd.h>     /* read(), close(), lseek() */
#include <assert.h>     /* assert() */
#include <errno.h>      /* errno */

#ifdef ENABLE_THREAD_SAFE
#include<pthread.h>
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef ENABLE_NETCDF4
#include <netcdf.h>
#endif

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

#ifndef MAX_INT_LEN
#define MAX_INT_LEN 24
#endif

#ifdef ENABLE_ADIOS
#include "adios_read.h"
#include <arpa/inet.h>
#define BP_MINIFOOTER_SIZE 28
#define BUFREAD64(buf,var) memcpy(&var, buf, 8); if (diff_endian) swap_64(&var);
#endif

/* Note accessing the following 5 global variables must be protected by a
 * mutex, otherwise it will not be thread safe.
 */

/* static variables are initialized to NULLs */
static PNC *pnc_filelist[NC_MAX_NFILES];
static int  pnc_numfiles;

/* This is the default create format for ncmpi_create and nc__create.
 * The use of this file scope variable is not thread-safe.
 */
static int ncmpi_default_create_format = NC_FORMAT_CLASSIC;

/* attribute to be cached in all communicators */
static int ncmpi_comm_keyval = MPI_KEYVAL_INVALID;

/* attribute to be cached in MPI_COMM_SELF */
static int ncmpi_init_keyval = MPI_KEYVAL_INVALID;

#define NCMPII_HANDLE_ERROR(func)                                  \
    if (mpireturn != MPI_SUCCESS) {                                \
        int errorStringLen;                                        \
        char errorString[MPI_MAX_ERROR_STRING];                    \
        MPI_Error_string(mpireturn, errorString, &errorStringLen); \
        printf("Error: file %s line %d calling func %s: (%s)\n",   \
               __FILE__, __LINE__, func, errorString);             \
    }

#define CHECK_ERRNO(err, func) {                                   \
    if (err != 0) {                                                \
        printf("Error: file %s line %d calling func %s: (%s)\n",   \
               __FILE__, __LINE__, func, strerror(err));           \
        goto err_out;                                              \
    }                                                              \
}

/* struct PNC_comm_attr is defined in dispatch.h */

/*----< ina_init() >---------------------------------------------------------*/
/* When the intra-node write aggregation (INA) hint is enabled, this subroutine
 * initializes the metadata to be used in intra- and inter-node communication,
 * including an intra-node communicator for communication between INA
 * aggregators and non-INA aggregators, and an inter-node communicator
 * consisting of all the INA aggregators to be used when calling GIO/MPI-IO
 * file open and their collective and independent I/O calls.
 *
 * Processes on the same NUMA node will first be divided into disjoined groups.
 * The passed in communicator will be split into sub-communicators, referred to
 * as intra-node communicators, one for each INA group. Within a group, the
 * process with the lowest rank ID is selected as the group's INA aggregator.
 *
 * A new MPI communicator consisting of all INA aggregators across all nodes
 * will be created, which is referred to as the inter-node communicator. The
 * inter-node communicator will be used to open the file, i.e. in the call to
 * GIO_open()/MPI_File_open() later. Only the INA aggregators make calls to
 * the GIO/MPI-IO APIs to access the file. Thus, this subroutine must be
 * called before opening the file and should be called only once in
 * ncmpi_create() or ncmpi_open().
 *
 * This subroutine performs the following tasks.
 * 1. Makes use of the affinity of each MPI process to its NUMA node to
 *    calculate the number of INA groups within a node and identify the
 *    membership of every process to its INA group.
 *    + comm_attr->num_NUMAs is the number of NUMA nodes.
 *    + comm_attr->NUMA_IDs[nprocs] contains NUMA node IDs of all processes.
 *    Note comm_attr should have already been established during a call to
 *    ncmpii_construct_node_list() at the beginning of ncmpi_create() and
 *    ncmpi_open().
 * 2. Based on hint num_aggrs_per_node and the number of processes per NUMA
 *    node, calculates the number of INA groups per node and divides processes
 *    into groups. Select the process with the lowest rank in a group as the
 *    INA aggregator.
 *    + comm_attr->is_ina_aggr indicates whether this rank is an INA aggregator
 * 3. Create a new MPI communicator by splitting 'comm' into sub-communicators,
 *    each consisting of processes belonging to the same INA group.
 *    + comm_attr->ina_intra_comm is the INA intra-node communicator.
 *    + MPI_Comm_size(comm_attr->ina_intra_comm, &size) is the number of
 *      processes within an INA group.
 *    + Rank 0 in comm_attr->ina_intra_comm is the INA aggregator.
 * 4. Create a new MPI communicator consisting of all the INA aggregators.
 *    + comm_attr->ina_inter_comm is the INA inter-node communicator.
 *    + MPI_Comm_size(comm_attr->ina_inter_comm, &size) is the total number of
 *      INA aggregators.
 *    + comm_attr->ina_inter_comm will be used when calling GIO_open()/
 *      MPI_File_open().
 */
static
int ina_init(MPI_Comm        comm,
             int             num_aggrs_per_node,
             PNC_comm_attr  *comm_attr)
{
    int i, j, k, mpireturn, nprocs, rank, my_node_naggrs, aggr_rank;
    int my_node_nprocs, my_node_rank;
    int *ina_flags, grp_nprocs, rem;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif

    if (num_aggrs_per_node == 0) return NC_NOERR;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

#ifdef PNETCDF_DEBUG
    /* Note that ill value of num_aggrs_per_node has been checked before
     * entering this subroutine. Thus num_aggrs_per_node must be > 0.
     */
    assert(num_aggrs_per_node > 0);
    assert(comm_attr->numa_comm != MPI_COMM_NULL);
#endif

    /* comm_attr->NUMA_IDs[] has been set in ncmpii_construct_node_list()
     * called earlier in ncmpio_create() or ncmpio_open() before entering this
     * subroutine.
     */

    /* my_node_nprocs is the number of processes in my NUMA compute node. */
    MPI_Comm_size(comm_attr->numa_comm, &my_node_nprocs);

    /* my_node_rank is the rank ID of this process in my NUMA compute node. */
    MPI_Comm_rank(comm_attr->numa_comm, &my_node_rank);

    /* Make sure the actual number of INA aggregators per node, initially set
     * in hint num_aggrs_per_node, is <= my_node_nprocs. In some cases, the
     * number of processes allocated to the last few compute nodes can be less
     * than others.
     */
    my_node_naggrs = MIN(num_aggrs_per_node, my_node_nprocs);

    /* Divide processes in a NUMA node into INA groups and calculate the number
     * of processes in each INA group, grp_nprocs. Select the INA aggregator,
     * as the process with loweset rank in an INA group, whose local rank in
     * comm_attr->numa_comm is 'aggr_rank'.
     */
    grp_nprocs = my_node_nprocs / my_node_naggrs; /* no. processes per group */
    rem        = my_node_nprocs % my_node_naggrs;
    if (rem > 0) { /* non-divisible case */
        grp_nprocs++;
        if (my_node_rank < grp_nprocs * rem)
            /* Select the first rank of my INA group as INA aggregator. */
            aggr_rank = my_node_rank - my_node_rank % grp_nprocs;
       else {
            aggr_rank = grp_nprocs * rem;
            grp_nprocs--;
            aggr_rank = my_node_rank
                       - (my_node_rank - aggr_rank) % grp_nprocs;
        }
    }
    else /* divisible case */
        aggr_rank = my_node_rank - my_node_rank % grp_nprocs;

    /* whether this rank is an INA aggregator */
    comm_attr->is_ina_aggr = (my_node_rank == aggr_rank);

#ifdef PNETCDF_DEBUG
    /* Make sure the number of processes in my INA group does not go beyond
     * my_node_nprocs.
     */
    assert(grp_nprocs <= my_node_nprocs - aggr_rank);
#endif

    if (comm_attr->ina_intra_comm != MPI_COMM_NULL) {
        /* free ina_intra_comm if previous created */
        TRACE_COMM(MPI_Comm_free)(&comm_attr->ina_intra_comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_free");
    }
    comm_attr->ina_intra_comm = MPI_COMM_NULL;

    /* Split NUMA comm into INA intra-node comm based on the assigned INA
     * aggregator's rank ID, i.e. processes sharing the same INA aggregator
     * form an ina_intra_comm. This process's local rank on the NUMA node
     * is my_node_rank and its INA aggregator's rank is aggr_rank.
     *
     * Special note on when the number of processes in this INA group is 1. In
     * this case, there is only one process in this group and thus the
     * intra-node aggregation of this group will not perform. However, the
     * intra-node communicator of this INA group must still be created. This is
     * because MPI_Comm_split is a collective call with the processes running
     * on the same NUMA node, i.e. on comm_attr->numa_comm.
     *
     * The case of grp_nprocs == 1, does not mean intra-node aggregation is
     * disabled globally. It just means this group will not perform INA
     * aggregation. The indicator of whether intra-node aggregation is globally
     * enabled or disabled is 'num_aggrs_per_node', whose value must be kept
     * consistent across all processes. It is possible for some groups
     * containing only one process, in which case the aggregation is not
     * necessarily to perform within those groups.
     */
    TRACE_COMM(MPI_Comm_split)(comm_attr->numa_comm, aggr_rank, my_node_rank,
                               &comm_attr->ina_intra_comm);

    /* Next step is to construct an inter-node MPI communicator consisting of
     * all INA aggregators. It will later be used to call MPI_File_open(), and
     * successively only INA aggregators call MPI-IO functions to access the
     * file.
     */

    /* construct an array containing ranks of aggregators */
    ina_flags = (int*) malloc(sizeof(int) * nprocs);
    TRACE_COMM(MPI_Allgather)(&comm_attr->is_ina_aggr, 1, MPI_INT, ina_flags,
                              1, MPI_INT, comm);

    /* Given a comm_attr->NUMA_IDs[], the rank IDs of INA aggregators will
     * dependent on the layout of MPI process allocation to the compute nodes.
     * The common layouts can be two kinds:
     *   + cyclic - MPI ranks are assigned to nodes round-robin-ly,
     *   + block - MPI ranks are assigned to a node and then move on to next.
     *
     * Below uses an example of nodes=3, nprocs=10, * num_aggrs_per_node=2.
     * comm_attr->NUMA_IDs[] should be
     *     block  process allocation: 0,0,0,0,1,1,1,2,2,2
     *     cyclic process allocation: 0,1,2,0,1,2,0,1,2,0
     * Accordingly, rank IDs of INA aggregators can be two kinds
     *     block  process allocation: 1,0,1,0,1,0,1,1,0,1
     *     cyclic process allocation: 1,1,1,0,0,0,1,1,1,0
     */

    /* calculate actual number of INA aggregators */
    comm_attr->num_ina_aggrs = 0;
    for (j=0; j<nprocs; j++)
        comm_attr->num_ina_aggrs += ina_flags[j];

    /* Collect aggregators' rank IDs and store them in an increasing order of
     * node IDs. Note rank IDs in ina_ranks[] are relative to comm (not
     * inter-node comm or intra-node comm).
     */
    comm_attr->ina_ranks = (int*)malloc(sizeof(int) * comm_attr->num_ina_aggrs);
    for (k=0, i=0; i<comm_attr->num_NUMAs; i++) {
        for (j=0; j<nprocs; j++) {
            if (comm_attr->NUMA_IDs[j] == i && ina_flags[j] > 0)
                comm_attr->ina_ranks[k++] = j;
        }
    }
    free(ina_flags);

    if (comm_attr->ina_inter_comm != MPI_COMM_NULL) {
        /* free comm_attr->ina_inter_comm if created previously */
        TRACE_COMM(MPI_Comm_free)(&comm_attr->ina_inter_comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_free");
    }
    comm_attr->ina_inter_comm = MPI_COMM_NULL;

    /* create an inter-node communicator consisting of all INA aggregators */
    MPI_Group origin_group, ina_group;
    TRACE_COMM(MPI_Comm_group)(comm, &origin_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_group");
    TRACE_COMM(MPI_Group_incl)(origin_group, comm_attr->num_ina_aggrs,
                               comm_attr->ina_ranks, &ina_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_incl");
    TRACE_COMM(MPI_Comm_create)(comm, ina_group, &comm_attr->ina_inter_comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_create");
    TRACE_COMM(MPI_Group_free)(&ina_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_free");
    TRACE_COMM(MPI_Group_free)(&origin_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_free");

    /* TODO: automatically determine whether or not to enable intra-node
     * aggregation.
     *
     * The ideal case is it can be determined right before each collective
     * write call, because only at that time, the communication pattern is
     * known. If the pattern can cause contention, then enable it. Otherwise,
     * disable it.
     *
     * Such mechanism may depends on the followings.
     *   1. MPI-IO hint cb_noddes, and striping_unit
     *   2. calculate aggregate access region
     *   3. If the number of senders to each cb_nodes is very large, then
     *      intra-node aggregation should be enabled.
     *   4. Average of nprocs_per_node across all processes may be a factor for
     *      determining whether to enable intra-node aggregation. It indicates
     *      whether the high number of processes are allocated on the same
     *      node.
     */

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    pnc_ina_init = MPI_Wtime() - timing;
#endif

    return NC_NOERR;
}

/*----< comm_attr_copy() >---------------------------------------------------*/
/* A function to be invoked when a communicator is duplicated, which adds a
 * reference to the already allocated memory space storing node ID array.
 */
static
int comm_attr_copy(MPI_Comm  comm,
                   int       keyval,
                   void     *extra,
                   void     *attr_inP,
                   void     *attr_outP,
                   int      *flag)
{
    PNC_comm_attr *attr_in   = (PNC_comm_attr*) attr_inP;
    PNC_comm_attr **attr_out = (PNC_comm_attr**)attr_outP;

    if (attr_in == NULL)
        return MPI_ERR_KEYVAL;
    else
        attr_in->ref_count++;

    *attr_out = attr_in;

    *flag = 1;  /* make a copy in the new communicator */

    return MPI_SUCCESS;
}

/*----< comm_attr_delete() >-------------------------------------------------*/
/* Callback function to be called when a communicator is freed, which frees the
 * allocated memory space of node ID array.
 */
static
int comm_attr_delete(MPI_Comm  comm,
                     int       keyval,
                     void     *attr_val,
                     void     *extra)
{
    PNC_comm_attr *attr = (PNC_comm_attr*) attr_val;

    if (attr == NULL)
        return MPI_ERR_KEYVAL;
    else
        attr->ref_count--;

    if (attr->ref_count <= 0) {
        /* free the allocated array */
        if (attr->NUMA_IDs != NULL)
            free(attr->NUMA_IDs);
        if (attr->ina_ranks != NULL)
            free(attr->ina_ranks);
        if (attr->numa_comm != MPI_COMM_NULL)
            MPI_Comm_free(&attr->numa_comm);
        if (attr->ina_inter_comm != MPI_COMM_NULL)
            MPI_Comm_free(&attr->ina_inter_comm);
        if (attr->ina_intra_comm != MPI_COMM_NULL)
            MPI_Comm_free(&attr->ina_intra_comm);
        free(attr);
    }
    return MPI_SUCCESS;
}

/*----< end_call() >---------------------------------------------------------*/
/* Callback function to be called at MPI_Finalize(), which frees all cached
 * attributes.
 */
static
int end_call(MPI_Comm  comm,
             int       keyval,
             void     *attribute_val,
             void     *extra_state)
{
    /* Free all keyvals used by PnetCDF */

    MPI_Comm_free_keyval(&keyval); /* free ncmpi_init_keyval */

    if (ncmpi_comm_keyval != MPI_KEYVAL_INVALID)
        MPI_Comm_free_keyval(&ncmpi_comm_keyval);

    return MPI_SUCCESS;
}

/*----< set_get_comm_attr() >------------------------------------------------*/
/* Create/set/get attributes into/from the MPI communicators passed in from
 * the user application.
 */
static
int set_get_comm_attr(MPI_Comm        comm,
                      int             num_aggrs_per_node,
                      PNC_comm_attr  *attrP)
{
    int err, nprocs;
    PNC_comm_attr *attr;

    MPI_Comm_size(comm, &nprocs);

    if (ncmpi_init_keyval == MPI_KEYVAL_INVALID) {
        /* This is the first call ever to PnetCDF API. Creating key
         * ncmpi_init_keyval is necessary for MPI_Finalize() to free key
         * ncmpi_comm_keyval.
         */
        err = MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, end_call,
                                     &ncmpi_init_keyval, (void*)0);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Comm_create_keyval");

        err = MPI_Comm_set_attr(MPI_COMM_SELF, ncmpi_init_keyval, (void*)0);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Comm_set_attr");
    }

    if (ncmpi_comm_keyval == MPI_KEYVAL_INVALID) {
        err = MPI_Comm_create_keyval(comm_attr_copy,
                                     comm_attr_delete,
                                     &ncmpi_comm_keyval, NULL);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Comm_create_keyval");
    }

    if (ncmpi_comm_keyval != MPI_KEYVAL_INVALID) {
        int found;

        err = MPI_Comm_get_attr(comm, ncmpi_comm_keyval, &attr, &found);
        if (err != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(err, "MPI_Comm_get_attr");

        if (!found) {
            /* Construct an array storing node IDs of all processes. Note the
             * memory allocated for setting the attribute will be freed in
             * comm_attr_delete(), a callback function invoked when the
             * MPI communicator is freed.
             */
            attr = (PNC_comm_attr*) malloc(sizeof(PNC_comm_attr));
            attr->ref_count = 1;

            /* initialize INA metadata of intra-node aggregation */
            attr->num_aggrs_per_node = 0;
            attr->num_ina_aggrs = 0;
            attr->is_ina_aggr = 0;
            attr->ina_ranks = NULL;
            attr->numa_comm      = MPI_COMM_NULL;
            attr->ina_inter_comm = MPI_COMM_NULL;
            attr->ina_intra_comm = MPI_COMM_NULL;

            if (nprocs == 1) {
                attr->num_NUMAs = 1;
                attr->NUMA_IDs = (int*) malloc(sizeof(int));
                attr->NUMA_IDs[0] = 0;
            }
            else {
                /* Constructing NUMA compute node IDs requires communication
                 * calls to MPI_Comm_split_type(), MPI_Bcast(), and
                 * MPI_Allgather().
                 */
                err = ncmpii_construct_node_list(comm, &attr->num_NUMAs,
                                                 &attr->NUMA_IDs,
                                                 &attr->numa_comm);
                if (err != NC_NOERR)
                    return err;

                /* If INA is enabled, construct INA metadata */
                err = ina_init(comm, num_aggrs_per_node, attr);
                if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
                attr->num_aggrs_per_node = num_aggrs_per_node;
            }

            /* FYI. The same key ncmpi_comm_keyval can be added to different
             * MPI communicators with same or different values.
             */
            err = MPI_Comm_set_attr(comm, ncmpi_comm_keyval, attr);
            if (err != MPI_SUCCESS)
                return ncmpii_error_mpi2nc(err, "MPI_Comm_set_attr");
        }
        else {
            if (num_aggrs_per_node == 0) {
                /* When INA is disabled, there is no need to retrieved cached
                 * INA metadata.
                 */
                *attrP = *attr; /* retrieve num_NUMAs and ids */
                attrP->num_aggrs_per_node = 0;
                attrP->num_ina_aggrs = 0;
                attrP->is_ina_aggr = 01;
                attrP->ina_ranks = NULL;
                attrP->ina_inter_comm = MPI_COMM_NULL;
                attrP->ina_intra_comm = MPI_COMM_NULL;
                return NC_NOERR;
            }

            /* When num_aggrs_per_node has changed, the cached attribute must
             * be reset. MPI standard says calling MPI_Comm_set_attr() will
             * trigger MPI_Comm_delete_attr() first, which will free all
             * allocated space, e.g. NUMA_IDs[] and numa_comm. Thus, we must
             * preserve attr->num_NUMAs, attr->NUMA_IDs[], and attr->numa_comm.
             */
            if (num_aggrs_per_node != attr->num_aggrs_per_node) {
                PNC_comm_attr *new_attr;

                /* updating an attribute must delete it first and set again */
                new_attr = (PNC_comm_attr*) malloc(sizeof(PNC_comm_attr));

                /* copy out attributes that remains with the communicator */
                new_attr->ref_count = attr->ref_count;
                new_attr->num_NUMAs = attr->num_NUMAs;
                new_attr->NUMA_IDs = (int*) malloc(sizeof(int) * nprocs);
                memcpy(new_attr->NUMA_IDs, attr->NUMA_IDs, sizeof(int)*nprocs);

                /* Must duplicate numa_comm that has been established, because
                 * it will be freed at the call to MPI_Comm_set_attr() below.
                 * Once created, numa_comm remains fixed with 'comm'. Unlike
                 * ina_inter_comm and ina_intra_comm, numa_comm does not change
                 * when num_aggrs_per_node changes.
                 */
                MPI_Comm_dup(attr->numa_comm, &new_attr->numa_comm);

                /* update with new INA metadata */
                new_attr->num_aggrs_per_node = num_aggrs_per_node;
                new_attr->num_ina_aggrs = 0;
                new_attr->is_ina_aggr = 0;
                new_attr->ina_ranks = NULL;
                new_attr->ina_inter_comm = MPI_COMM_NULL;
                new_attr->ina_intra_comm = MPI_COMM_NULL;

                /* re-construct INA metadata */
                err = ina_init(comm, num_aggrs_per_node, new_attr);
                if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

                err = MPI_Comm_set_attr(comm, ncmpi_comm_keyval, new_attr);
                if (err != MPI_SUCCESS)
                    return ncmpii_error_mpi2nc(err, "MPI_Comm_set_attr");
                attr = new_attr;
            }
        }

        /* copy contents */
        *attrP = *attr;
    }

    return NC_NOERR;
}

/*----< new_id_PNCList() >---------------------------------------------------*/
/* Return a new ID (array index) from the PNC list, pnc_filelist[] that is
 * not used. Note the used elements in pnc_filelist[] may not be contiguous.
 * For example, some files created/opened later may be closed earlier than
 * others, leaving those array elements NULL in the middle.
 */
static int
new_id_PNCList(int *new_id, PNC *pncp)
{
    int i, err=NC_NOERR, perr=0;

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_lock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_lock")
#endif
    *new_id = -1;
    if (pnc_numfiles == NC_MAX_NFILES) { /* Too many files open */
        DEBUG_ASSIGN_ERROR(err, NC_ENFILE)
    }
    else {
        err = NC_NOERR;
        for (i=0; i<NC_MAX_NFILES; i++) { /* find the first unused element */
            if (pnc_filelist[i] == NULL) {
                *new_id = i;
                pnc_filelist[i] = pncp;
                pnc_numfiles++; /* increment number of files opened */
                break;
            }
        }
    }
#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_unlock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_unlock")

err_out:
#endif
    return (err != NC_NOERR) ? err : perr;
}

/*----< del_from_PNCList() >-------------------------------------------------*/
static int
del_from_PNCList(int ncid)
{
    int perr=0;

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_lock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_lock")
#endif

    /* validity of ncid should have been checked already */
    pnc_filelist[ncid] = NULL;
    pnc_numfiles--;

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_unlock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_unlock")

err_out:
#endif
    return perr;
}

/*----< PNC_check_id() >-----------------------------------------------------*/
int
PNC_check_id(int ncid, PNC **pncp)
{
    int err=NC_NOERR, perr=0;

    assert(pncp != NULL);

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_lock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_lock")
#endif

    if (pnc_numfiles == 0 || ncid < 0 || ncid >= NC_MAX_NFILES)
        DEBUG_ASSIGN_ERROR(err, NC_EBADID)
    else
        *pncp = pnc_filelist[ncid];

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_unlock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_unlock")

err_out:
#endif
    return (err != NC_NOERR) ? err : perr;
}

/*----< construct_info() >---------------------------------------------------*/
/* This subroutine reads I/O hints from the environment variable PNETCDF_HINTS,
 * if set at the run time. The value of PNETCDF_HINTS is a character string
 * consisting of one or more hints separated by ";" and each hint is in the
 * form of "hint=value". E.g. "cb_nodes=16;cb_config_list=*:6".
 *
 * Hints set in PNETCDF_HINTS environment variable takes the highest precedence
 * over hints set in the MPI info object passed from the application programs.
 *
 * When use_info is MPI_INFO_NULL and PNETCDF_HINTS is empty, MPI_INFO_NULL
 * will be returned as new_info.
 */
static void
combine_env_hints(MPI_Info  user_info,  /* IN */
                  MPI_Info *new_info)   /* OUT: may be MPI_INFO_NULL */
{
    char *warn_str="Warning: skip ill-formed hint set in PNETCDF_HINTS";
    char *env_str;
    char *hdr_align_val=NULL, *var_align_val=NULL;

    if (user_info != MPI_INFO_NULL)
        MPI_Info_dup(user_info, new_info); /* ignore error */
    else
        *new_info = MPI_INFO_NULL;

    /* get environment variable PNETCDF_HINTS */
    if ((env_str = getenv("PNETCDF_HINTS")) != NULL) {
#ifdef USE_STRTOK_R
        char *env_str_cpy, *env_str_saved, *hint, *key;
        env_str_cpy = NCI_Strdup(env_str);
        env_str_saved = env_str_cpy;
        hint = strtok_r(env_str_cpy, ";", &env_str_saved);
        while (hint != NULL) {
            char *hint_saved = NCI_Strdup(hint);
            char *val = strchr(hint, '=');
            if (val == NULL) { /* ill-formed hint */
                if (NULL != strtok(hint, " \t"))
                    printf("%s: '%s'\n", warn_str, hint_saved);
                /* else case: ignore white-spaced hints */
                NCI_Free(hint_saved);
                hint = strtok_r(NULL, ";", &env_str_saved); /* get next hint */
                continue;
            }
            key = strtok(hint, "= \t");
            val = strtok(NULL, "= \t");
            if (NULL != strtok(NULL, "= \t")) /* expect no more token */
                printf("%s: '%s'\n", warn_str, hint_saved);
            else {
                if (*new_info == MPI_INFO_NULL)
                    MPI_Info_create(new_info); /* ignore error */

                if (!strcmp(key, "nc_header_align_size"))
                    hdr_align_val = NCI_Strdup(val);
                else if (!strcmp(key, "nc_var_align_size"))
                    var_align_val = NCI_Strdup(val);
                else
                    MPI_Info_set(*new_info, key, val); /* override or add */
            }
            /* printf("env hint: key=%s val=%s\n",key,val); */
            hint = strtok_r(NULL, ";", &env_str_saved);
            NCI_Free(hint_saved);
        }
        NCI_Free(env_str_cpy);
#else
        char *env_str_cpy, *hint, *next_hint, *key, *val, *deli;
        char *hint_saved=NULL;

        env_str_cpy = NCI_Strdup(env_str);
        next_hint = env_str_cpy;

        do {
            hint = next_hint;
            deli = strchr(hint, ';');
            if (deli != NULL) {
                *deli = '\0'; /* add terminate char */
                next_hint = deli + 1;
            }
            else next_hint = "\0";
            if (hint_saved != NULL) NCI_Free(hint_saved);

            /* skip all-blank hint */
            hint_saved = NCI_Strdup(hint);
            if (strtok(hint, " \t") == NULL) continue;

            NCI_Free(hint_saved);
            hint_saved = NCI_Strdup(hint); /* save hint for error message */

            deli = strchr(hint, '=');
            if (deli == NULL) { /* ill-formed hint */
                printf("%s: '%s'\n", warn_str, hint_saved);
                continue;
            }
            *deli = '\0';

            /* hint key */
            key = strtok(hint, "= \t");
            if (key == NULL || NULL != strtok(NULL, "= \t")) {
                /* expect one token before = */
                printf("%s: '%s'\n", warn_str, hint_saved);
                continue;
            }

            /* hint value */
            val = strtok(deli+1, "= \t");
            if (NULL != strtok(NULL, "= \t")) { /* expect one token before = */
                printf("%s: '%s'\n", warn_str, hint_saved);
                continue;
            }
            if (*new_info == MPI_INFO_NULL)
                MPI_Info_create(new_info); /* ignore error */

            if (!strcmp(key, "nc_header_align_size"))
                hdr_align_val = NCI_Strdup(val);
            else if (!strcmp(key, "nc_var_align_size"))
                var_align_val = NCI_Strdup(val);
            else
                MPI_Info_set(*new_info, key, val); /* override or add */

        } while (*next_hint != '\0');

        if (hint_saved != NULL) NCI_Free(hint_saved);
        NCI_Free(env_str_cpy);
#endif

        /* nc_var_align_size supersedes nc_header_align_size */
        if (var_align_val != NULL) {
            MPI_Info_set(*new_info, "nc_var_align_size", var_align_val);
            MPI_Info_set(*new_info, "nc_header_align_size", var_align_val);
        }
        else if (hdr_align_val != NULL) {
            MPI_Info_set(*new_info, "nc_var_align_size", hdr_align_val);
            MPI_Info_set(*new_info, "nc_header_align_size", hdr_align_val);
        }
    }
    /* return no error as all hints are advisory */

    if (hdr_align_val != NULL) NCI_Free(hdr_align_val);
    if (var_align_val != NULL) NCI_Free(var_align_val);
}

/*----< set_env_mode() >-----------------------------------------------------*/
static
void set_env_mode(int *env_mode)
{
    char *env_str;

#ifdef PNETCDF_DEBUG
    fSet(*env_mode, NC_MODE_SAFE);
    /* When debug mode is enabled at the configure time, safe mode is by
     * default enabled. This can be overwritten by the run-time environment
     * variable PNETCDF_SAFE_MODE.
     */
#endif
    /* get environment variable PNETCDF_SAFE_MODE
     * if it is set to 1, then we perform a strict parameter consistent test
     */
    if ((env_str = getenv("PNETCDF_SAFE_MODE")) != NULL) {
        if (*env_str == '0') fClr(*env_mode, NC_MODE_SAFE);
        else                 fSet(*env_mode, NC_MODE_SAFE);
        /* if PNETCDF_SAFE_MODE is set but without a value, *env_str can
         * be '\0' (null character). In this case, safe mode is enabled */
    }

    /* get environment variable PNETCDF_RELAX_COORD_BOUND
     * if it is set to 0, then we perform a strict start bound check
     */
#ifdef RELAX_COORD_BOUND
    fSet(*env_mode, NC_MODE_STRICT_COORD_BOUND);
#endif
    if ((env_str = getenv("PNETCDF_RELAX_COORD_BOUND")) != NULL) {
        if (*env_str == '0') fClr(*env_mode, NC_MODE_STRICT_COORD_BOUND);
        else                 fSet(*env_mode, NC_MODE_STRICT_COORD_BOUND);
        /* if PNETCDF_RELAX_COORD_BOUND is set but without a value, *env_str
         * can be '\0' (null character). This is equivalent to setting
         * PNETCDF_RELAX_COORD_BOUND to 1 */
    }
}

/*----< ncmpi_create() >-----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_create(MPI_Comm    comm,
             const char *path,
             int         cmode,
             MPI_Info    info,
             int        *ncidp)
{
    char value[MPI_MAX_INFO_VAL];
    int rank, nprocs, status=NC_NOERR, err, flag;
    int env_mode=0, mpireturn, format, num_aggrs_per_node;
    MPI_Info combined_info=MPI_INFO_NULL;
    void *ncp;
    PNC *pncp;
    PNC_driver *driver;
    PNC_comm_attr comm_attr;
#ifdef BUILD_DRIVER_FOO
    int enable_foo_driver=0;
#endif
#ifdef ENABLE_BURST_BUFFER
    int enable_bb_driver=0;
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    {
        int i;
        pnc_num_aggrs_per_node = 0;
        pnc_ina_init = 0;
        pnc_ina_flatten = 0;
        pnc_ina_npairs_put = 0;
        pnc_ina_npairs_get = 0;

        for (i=0; i<NTIMERS; i++) {
            pnc_ina_put[i]     = pnc_ina_get[i] = 0;
            pnc_ina_mem_put[i] = pnc_ina_mem_get[i] = 0;
        }
    }
#endif

     *ncidp = -1;

    /* allocate a new PNC object */
    pncp = (PNC*) NCI_Malloc(sizeof(PNC));
    if (pncp == NULL)
        DEBUG_RETURN_ERROR(NC_ENOMEM) /* fatal error */

    pncp->path = NULL;

    /* The first thing is to duplicate the MPI communicator (even if comm is
     * MPI_COMM_WORLD or MPI_COMM_SELF) before any communication can be made
     * within PnetCDF. This is because users may use 'comm' to do other
     * point-to-point communication while PnetCDF also makes some MPI
     * communication calls. When this happened, user's communication can mess
     * up with the PnetCDF's communication, particularly when in independent
     * data mode. Once comm is duplicated, we pass pncp->comm to PnetCDF
     * drivers, so there is no need for a driver to duplicate it again.
     */
    pncp->comm = MPI_COMM_NULL;
    mpireturn = MPI_Comm_dup(comm, &pncp->comm);
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_dup");
        goto err_out;
    }

    /* Rudimentary check path's validity.
     *
     * Note MPI standard's requirement for filename argument is: " ... all
     * processes must provide filenames that reference the same file. ... The
     * user is responsible for ensuring that a single file is referenced by the
     * filename argument, as it may be impossible for an implementation to
     * detect this type of namespace error."
     */
    if (path == NULL || *path == '\0') {
        if (status == NC_NOERR) status = NC_EBAD_FILE;
        goto err_out;
    }

    MPI_Comm_rank(pncp->comm, &rank);
    MPI_Comm_size(pncp->comm, &nprocs);

    if (rank == 0)
        set_env_mode(&env_mode);

    /* duplicate file path */
    pncp->path = (char*) NCI_Strdup(path);
    if (pncp->path == NULL) {
        DEBUG_ASSIGN_ERROR(status, NC_ENOMEM)
        goto err_out;
    }

    if (nprocs > 1) { /* Check cmode consistency */
        int modes[2] = {cmode, env_mode}; /* only root's matters */

        TRACE_COMM(MPI_Bcast)(&modes, 2, MPI_INT, 0, pncp->comm);
        NCMPII_HANDLE_ERROR("MPI_Bcast")

        /* Overwrite cmode with root's cmode */
        if (modes[0] != cmode) {
            cmode = modes[0];
            DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE_CMODE)
        }

        env_mode = modes[1];
        if (fIsSet(env_mode, NC_MODE_SAFE)) {
            /* sync status among all processes */
            err = status;
            TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,
                                      pncp->comm);
            NCMPII_HANDLE_ERROR("MPI_Allreduce")
        }
        /* continue to use root's cmode to create the file, but will report
         * cmode inconsistency error, if there is any */
    }

    /* combine user's MPI info and PNETCDF_HINTS env variable */
    combine_env_hints(info, &combined_info);

    num_aggrs_per_node = 0;
    if (nprocs > 1 && combined_info != MPI_INFO_NULL) {
        /* check if INA hint is enabled */
        MPI_Info_get(combined_info, "nc_num_aggrs_per_node",
                     MPI_MAX_INFO_VAL-1, value, &flag);
        if (flag) {
            int ival;
            errno = 0;  /* errno must set to zero before calling atoi */
            ival = atoi(value);
            if (errno == 0 && ival >= 0)
                num_aggrs_per_node = ival;
        }
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    pnc_num_aggrs_per_node = num_aggrs_per_node;
#endif

#ifdef ENABLE_THREAD_SAFE
    int perr;
    perr = pthread_mutex_lock(&lock);
    if (perr != 0)
        printf("Warning in file %s line %d: pthread_mutex_lock() failed (%s)\n",
               __FILE__, __LINE__, strerror(perr));
#endif

    /* creating communicator attributes must be protected by a mutex */
    set_get_comm_attr(pncp->comm, num_aggrs_per_node, &comm_attr);
    /* ignore error, as it is not a critical error */

#if PNETCDF_DEBUG_MODE == 1
    if (num_aggrs_per_node == 0)
        assert(comm_attr.ina_intra_comm == MPI_COMM_NULL);
    else
        assert(comm_attr.ina_intra_comm != MPI_COMM_NULL);
#endif

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_unlock(&lock);
    if (perr != 0)
        printf("Warning in file %s line %d: pthread_mutex_unlock() failed (%s)\n",
               __FILE__, __LINE__, strerror(perr));
#endif

    if (combined_info == MPI_INFO_NULL)
        MPI_Info_create(&combined_info);

    /* add this rank's NUMA node ID */
    snprintf(value, MAX_INT_LEN, "%d", comm_attr.NUMA_IDs[rank]);
    MPI_Info_set(combined_info, "NUMA_ID", value);

#ifdef BUILD_DRIVER_FOO
    if (combined_info != MPI_INFO_NULL) {
        /* check if nc_foo_driver is enabled */
        MPI_Info_get(combined_info, "nc_foo_driver", MPI_MAX_INFO_VAL-1,
                    value, &flag);
        if (flag && strcasecmp(value, "enable") == 0)
            enable_foo_driver = 1;
    }
#endif
#ifdef ENABLE_BURST_BUFFER
    if (combined_info != MPI_INFO_NULL) {
        /* check if nc_burst_buf is enabled */
        MPI_Info_get(combined_info, "nc_burst_buf", MPI_MAX_INFO_VAL-1,
                    value, &flag);
        if (flag && strcasecmp(value, "enable") == 0)
            enable_bb_driver = 1;
    }
#endif

    /* Use environment variable and cmode to tell the file format
     * which is later used to select the right driver.
     */

#ifdef ENABLE_NETCDF4
    /* It is illegal to have NC_64BIT_OFFSET & NC_64BIT_DATA & NC_NETCDF4 */
    if ((cmode & (NC_64BIT_OFFSET|NC_NETCDF4)) ==
                 (NC_64BIT_OFFSET|NC_NETCDF4) ||
        (cmode & (NC_64BIT_DATA|NC_NETCDF4)) ==
                 (NC_64BIT_DATA|NC_NETCDF4)) {
        DEBUG_ASSIGN_ERROR(status, NC_EINVAL_CMODE)
        goto err_out;
    }
#else
    if (cmode & NC_NETCDF4) {
        DEBUG_ASSIGN_ERROR(status, NC_ENOTBUILT)
        goto err_out;
    }
#endif

    /* It is illegal to have both NC_64BIT_OFFSET & NC_64BIT_DATA */
    if ((cmode & (NC_64BIT_OFFSET|NC_64BIT_DATA)) ==
                 (NC_64BIT_OFFSET|NC_64BIT_DATA)) {
        DEBUG_ASSIGN_ERROR(status, NC_EINVAL_CMODE)
        goto err_out;
    }

    /* Check if cmode contains format specific flag */
    if (fIsSet(cmode, NC_64BIT_DATA))
        format = NC_FORMAT_CDF5;
    else if (fIsSet(cmode, NC_64BIT_OFFSET))
        format = NC_FORMAT_CDF2;
    else if (fIsSet(cmode, NC_NETCDF4)) {
        if (fIsSet(cmode, NC_CLASSIC_MODEL))
            format = NC_FORMAT_NETCDF4_CLASSIC;
        else
            format = NC_FORMAT_NETCDF4;
    }
    else if (fIsSet(cmode, NC_CLASSIC_MODEL))
        format = NC_FORMAT_CLASSIC;
    else {
        /* if no file format flag is set in cmode, use default */
        ncmpi_inq_default_format(&format);
        if (format == NC_FORMAT_CDF5)
            cmode |= NC_64BIT_DATA;
        else if (format == NC_FORMAT_CDF2)
            cmode |= NC_64BIT_OFFSET;
        else if (format == NC_FORMAT_NETCDF4)
            cmode |= NC_NETCDF4;
        else if (format == NC_FORMAT_NETCDF4_CLASSIC)
            cmode |= NC_NETCDF4 | NC_CLASSIC_MODEL;
    }

#ifdef ENABLE_NETCDF4
    if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC) {
        driver = nc4io_inq_driver();
#ifdef ENABLE_BURST_BUFFER
        /* Burst buffering does not support NetCDF-4 files yet.
         * If hint nc_burst_buf is enabled in combined_info, disable it.
         */
        if (enable_bb_driver == 1 && combined_info != MPI_INFO_NULL)
            MPI_Info_set(combined_info, "nc_burst_buf", "disable");
        enable_bb_driver = 0;
#endif
    }
    else
#endif
#ifdef BUILD_DRIVER_FOO
    if (enable_foo_driver)
        driver = ncfoo_inq_driver();
    else
#endif
#ifdef ENABLE_BURST_BUFFER
    if (enable_bb_driver)
        driver = ncbbio_inq_driver();
    else
#endif
        /* default is the driver built on top of MPI-IO */
        driver = ncmpio_inq_driver();

    pncp->flag = NC_MODE_DEF | NC_MODE_CREATE;
    fSet(pncp->flag, env_mode);

    /* generate a new nc file ID from NCPList */
    err = new_id_PNCList(ncidp, pncp);
    if (err != NC_NOERR) {
        if (combined_info != MPI_INFO_NULL)
            MPI_Info_free(&combined_info);
        DEBUG_ASSIGN_ERROR(status, err)
        goto err_out;
    }

    /* calling the driver's create subroutine */
    err = driver->create(pncp->comm, pncp->path, cmode, *ncidp, env_mode,
                         combined_info, comm_attr, &ncp);
    if (status == NC_NOERR) status = err;
    if (status != NC_NOERR && status != NC_EMULTIDEFINE_CMODE) {
        del_from_PNCList(*ncidp);
        goto err_out;
    }

    pncp->mode       = cmode;
    pncp->driver     = driver;
    pncp->ndims      = 0;
    pncp->unlimdimid = -1;
    pncp->nvars      = 0;
    pncp->nrec_vars  = 0;
    pncp->vars       = NULL;
    pncp->ncp        = ncp;
    pncp->format     = format;

    if (fIsSet(env_mode, NC_MODE_SAFE))
        pncp->flag |= NC_MODE_SAFE;

    if (fIsSet(env_mode, NC_MODE_STRICT_COORD_BOUND))
        pncp->flag |= NC_MODE_STRICT_COORD_BOUND;

err_out:
    if (combined_info != MPI_INFO_NULL)
        MPI_Info_free(&combined_info);

    if (status != NC_NOERR && status != NC_EMULTIDEFINE_CMODE) {
        if (pncp->comm != MPI_COMM_NULL)
            MPI_Comm_free(&pncp->comm); /* a collective call */
        if (pncp->path != NULL)
            NCI_Free(pncp->path);
        NCI_Free(pncp);
        *ncidp = -1;
    }

    return status;
}

#define NDIMS_ 16

/*----< ncmpi_open() >-------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_open(MPI_Comm    comm,
           const char *path,
           int         omode,
           MPI_Info    info,
           int        *ncidp)  /* OUT */
{
    char value[MPI_MAX_INFO_VAL];
    int i, j, nalloc, rank, nprocs, format, status=NC_NOERR, err, flag;
    int env_mode=0, mpireturn, DIMIDS[NDIMS_], *dimids, num_aggrs_per_node;
    MPI_Info combined_info=MPI_INFO_NULL;
    void *ncp;
    PNC *pncp;
    PNC_driver *driver;
    PNC_comm_attr comm_attr;
#ifdef BUILD_DRIVER_FOO
    int enable_foo_driver=0;
#endif
#ifdef ENABLE_BURST_BUFFER
    int enable_bb_driver=0;
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    pnc_num_aggrs_per_node = 0;
    pnc_ina_init = 0;
    pnc_ina_flatten = 0;
    pnc_ina_npairs_put = 0;
    pnc_ina_npairs_get = 0;

    for (i=0; i<NTIMERS; i++) {
        pnc_ina_put[i]     = pnc_ina_get[i] = 0;
        pnc_ina_mem_put[i] = pnc_ina_mem_get[i] = 0;
    }
#endif

    *ncidp = -1;

    /* allocate a PNC object */
    pncp = (PNC*) NCI_Malloc(sizeof(PNC));
    if (pncp == NULL)
        DEBUG_RETURN_ERROR(NC_ENOMEM) /* fatal error */

    pncp->path = NULL;

    /* The first thing is to duplicate the MPI communicator (even if comm is
     * MPI_COMM_WORLD or MPI_COMM_SELF) before any communication can be made
     * within PnetCDF. This is because users may use 'comm' to do other
     * point-to-point communication while PnetCDF also makes some MPI
     * communication calls. When this happened, user's communication can mess
     * up with the PnetCDF's communication, particularly when in independent
     * data mode. Once comm is duplicated, we pass pncp->comm to PnetCDF
     * drivers, so there is no need for a driver to duplicate it again.
     */
    pncp->comm = MPI_COMM_NULL;
    mpireturn = MPI_Comm_dup(comm, &pncp->comm);
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_dup");
        if (status == NC_NOERR) status = err;
        goto err_out;
    }

    /* Rudimentary check path's validity.
     *
     * Note MPI standard's requirement for filename argument is: " ... all
     * processes must provide filenames that reference the same file. ... The
     * user is responsible for ensuring that a single file is referenced by the
     * filename argument, as it may be impossible for an implementation to
     * detect this type of namespace error."
     */
    if (path == NULL || *path == '\0') {
        if (status == NC_NOERR) status = NC_EBAD_FILE;
        goto err_out;
    }

    MPI_Comm_rank(pncp->comm, &rank);
    MPI_Comm_size(pncp->comm, &nprocs);

    if (rank == 0)
        set_env_mode(&env_mode);

    /* Check the file signature to tell the file format which is later used to
     * select the right driver.
     */
    format = NC_FORMAT_UNKNOWN;
    if (rank == 0) {
        err = ncmpi_inq_file_format(path, &format);
        if (err != NC_NOERR) {
            if (nprocs == 1) {
                if (status == NC_NOERR) status = err;
                goto err_out;
            }
            format = err;
        }
        else if (format == NC_FORMAT_UNKNOWN) {
            if (nprocs == 1) {
                if (status == NC_NOERR) status = NC_ENOTNC;
                goto err_out;
            }
            format = NC_ENOTNC;
        }
#ifndef ENABLE_NETCDF4
        else if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC) {
            if (nprocs == 1) {
                if (status == NC_NOERR) status = NC_ENOTBUILT;
                goto err_out;
            }
            format = NC_ENOTBUILT;
        }
#endif
    }

    if (nprocs > 1) { /* root broadcasts format and omode */
        int modes[3] = {format, omode, env_mode};

        /* Check consistency:
         * Note only root's values matter, format, omode, env_mode.
         * Note if omode contains NC_NOWRITE, it is equivalent to NC_CLOBBER.
         * In pnetcdf.h, they both are defined the same value, 0.
         */

        TRACE_COMM(MPI_Bcast)(&modes, 3, MPI_INT, 0, pncp->comm);
        NCMPII_HANDLE_ERROR("MPI_Bcast")

        /* check format error (a fatal error, must return now) */
        format = modes[0];
        if (format < 0) { /* all netCDF errors are negative */
            if (status == NC_NOERR) status = format;
            goto err_out;
        }

        /* check omode consistency */
        if (modes[1] != omode) {
            omode = modes[1];
            if (status == NC_NOERR) status = NC_EMULTIDEFINE_OMODE;
        }

        env_mode = modes[2];
        if (fIsSet(env_mode, NC_MODE_SAFE)) {
            /* sync status among all processes */
            err = status;
            TRACE_COMM(MPI_Allreduce)(&err, &status, 1, MPI_INT, MPI_MIN,
                                      pncp->comm);
            NCMPII_HANDLE_ERROR("MPI_Allreduce")
        }
        /* continue to use root's omode to open the file, but will report omode
         * inconsistency error, if there is any
         */
    }

    /* combine user's MPI info and PNETCDF_HINTS env variable */
    combine_env_hints(info, &combined_info);

    num_aggrs_per_node = 0;
    if (nprocs > 1 && combined_info != MPI_INFO_NULL) {
        /* check if INA hint is enabled */
        MPI_Info_get(combined_info, "nc_num_aggrs_per_node",
                     MPI_MAX_INFO_VAL-1, value, &flag);
        if (flag) {
            int ival;
            errno = 0;  /* errno must set to zero before calling atoi */
            ival = atoi(value);
            if (errno == 0 && ival >= 0)
                num_aggrs_per_node = ival;
        }
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    pnc_num_aggrs_per_node = num_aggrs_per_node;
#endif

#ifdef ENABLE_THREAD_SAFE
    int perr;
    perr = pthread_mutex_lock(&lock);
    if (perr != 0)
        printf("Warning in file %s line %d: pthread_mutex_lock() failed (%s)\n",
               __FILE__, __LINE__, strerror(perr));
#endif

    /* creating communicator attributes must be protected by a mutex */
    set_get_comm_attr(pncp->comm, num_aggrs_per_node, &comm_attr);
    /* ignore error, as it is not a critical error */

#if PNETCDF_DEBUG_MODE == 1
    if (num_aggrs_per_node == 0)
        assert(comm_attr.ina_intra_comm == MPI_COMM_NULL);
    else
        assert(comm_attr.ina_intra_comm != MPI_COMM_NULL);
#endif

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_unlock(&lock);
    if (perr != 0)
        printf("Warning in file %s line %d: pthread_mutex_unlock() failed (%s)\n",
               __FILE__, __LINE__, strerror(perr));
#endif

    if (combined_info == MPI_INFO_NULL)
        MPI_Info_create(&combined_info);

    /* add this rank's NUMA node ID */
    snprintf(value, MAX_INT_LEN, "%d", comm_attr.NUMA_IDs[rank]);
    MPI_Info_set(combined_info, "NUMA_ID", value);

#ifdef BUILD_DRIVER_FOO
    if (combined_info != MPI_INFO_NULL) {
        /* check if nc_foo_driver is enabled */
        MPI_Info_get(combined_info, "nc_foo_driver", MPI_MAX_INFO_VAL-1,
                    value, &flag);
        if (flag && strcasecmp(value, "enable") == 0)
            enable_foo_driver = 1;
    }
#endif
#ifdef ENABLE_BURST_BUFFER
    if (combined_info != MPI_INFO_NULL) {
        /* check if nc_burst_buf is enabled */
        MPI_Info_get(combined_info, "nc_burst_buf", MPI_MAX_INFO_VAL-1,
                    value, &flag);
        if (flag && strcasecmp(value, "enable") == 0)
            enable_bb_driver = 1;
    }
#endif

#ifdef ENABLE_NETCDF4
    if (format == NC_FORMAT_NETCDF4_CLASSIC || format == NC_FORMAT_NETCDF4) {
        driver = nc4io_inq_driver();
#ifdef ENABLE_BURST_BUFFER
        /* Burst buffering does not support NetCDF-4 files yet.
         * If hint nc_burst_buf is enabled in combined_info, disable it.
         */
        if (enable_bb_driver == 1) {
            if (combined_info != MPI_INFO_NULL)
                MPI_Info_set(combined_info, "nc_burst_buf", "disable");
        }
        enable_bb_driver = 0;
#endif
    }
    else
#else
    if (format == NC_FORMAT_NETCDF4_CLASSIC || format == NC_FORMAT_NETCDF4) {
        if (status == NC_NOERR) status = NC_ENOTBUILT;
        goto err_out;
    }
    else
#endif
#ifdef BUILD_DRIVER_FOO
    if (enable_foo_driver)
        driver = ncfoo_inq_driver();
    else
#endif
#ifdef ENABLE_BURST_BUFFER
    if (enable_bb_driver)
        driver = ncbbio_inq_driver();
    else
#endif
    {
        /* ncmpio driver */
        if (format == NC_FORMAT_CLASSIC ||
            format == NC_FORMAT_CDF2 ||
            format == NC_FORMAT_CDF5) {
            driver = ncmpio_inq_driver();
        }
#ifdef ENABLE_ADIOS
        else if (format == NC_FORMAT_BP) {
            driver = ncadios_inq_driver();
        }
#endif
        else { /* unrecognized file format */
            if (status == NC_NOERR) status = NC_ENOTNC;
            goto err_out;
        }
    }

    pncp->flag = 0;
    fSet(pncp->flag, env_mode);

    /* duplicate file path */
    pncp->path = (char*) NCI_Strdup(path);
    if (pncp->path == NULL) {
        if (status == NC_NOERR) status = NC_ENOMEM;
        goto err_out;
    }

    /* generate a new nc file ID from NCPList */
    err = new_id_PNCList(ncidp, pncp);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        goto err_out;
    }

    /* calling the driver's open subroutine */
    err = driver->open(pncp->comm, pncp->path, omode, *ncidp, env_mode,
                       combined_info, comm_attr, &ncp);
    if (status == NC_NOERR) status = err;
    if (status != NC_NOERR && status != NC_EMULTIDEFINE_OMODE &&
        status != NC_ENULLPAD) {
        /* NC_EMULTIDEFINE_OMODE and NC_ENULLPAD are not fatal error. We
         * continue the rest open procedure */
        del_from_PNCList(*ncidp);
        goto err_out;
    }

    /* fill in pncp members */
    pncp->mode       = omode;
    pncp->driver     = driver;
    pncp->ndims      = 0;
    pncp->unlimdimid = -1;
    pncp->nvars      = 0;
    pncp->nrec_vars  = 0;
    pncp->vars       = NULL;
    pncp->ncp        = ncp;
    pncp->format     = format;

    if (!fIsSet(omode, NC_WRITE))
        pncp->flag |= NC_MODE_RDONLY;

    if (fIsSet(env_mode, NC_MODE_SAFE))
        pncp->flag |= NC_MODE_SAFE;

    if (fIsSet(env_mode, NC_MODE_STRICT_COORD_BOUND))
        pncp->flag |= NC_MODE_STRICT_COORD_BOUND;

    /* inquire number of dimensions, variables defined and rec dim ID */
    err = driver->inq(pncp->ncp, &pncp->ndims, &pncp->nvars, NULL,
                      &pncp->unlimdimid);
    if (err != NC_NOERR) {
        driver->close(ncp); /* close file and ignore error */
        del_from_PNCList(*ncidp);
        if (status == NC_NOERR) status = err;
        goto err_out;
    }

    if (pncp->nvars == 0) /* no variable defined in the file */
        goto err_out;

    /* make a copy of variable metadata at the dispatcher layer, because sanity
     * check is done at the dispatcher layer
     */

    /* allocate chunk size for pncp->vars[] */
    nalloc = PNETCDF_RNDUP(pncp->nvars, PNC_VARS_CHUNK);
    pncp->vars = NCI_Malloc(sizeof(PNC_var) * nalloc);
    if (pncp->vars == NULL) {
        driver->close(ncp); /* close file and ignore error */
        del_from_PNCList(*ncidp);
        if (status == NC_NOERR) status = NC_ENOMEM;
        goto err_out;
    }

    dimids = DIMIDS;

    /* construct array of PNC_var for all variables */
    for (i=0; i<pncp->nvars; i++) {
        int ndims, max_ndims=NDIMS_;
        pncp->vars[i].shape  = NULL;
        pncp->vars[i].recdim = -1;   /* if fixed-size variable */
        err = driver->inq_var(pncp->ncp, i, NULL, &pncp->vars[i].xtype, &ndims,
                              NULL, NULL, NULL, NULL, NULL);
        if (err != NC_NOERR) break; /* loop i */
        pncp->vars[i].ndims = ndims;

        if (ndims > 0) {
            pncp->vars[i].shape = (MPI_Offset*)
                                  NCI_Malloc(sizeof(MPI_Offset) * ndims);
            if (ndims > max_ndims) { /* avoid repeated malloc */
                if (dimids == DIMIDS) dimids = NULL;
                dimids = (int*) NCI_Realloc(dimids, sizeof(int) * ndims);
                max_ndims = ndims;
            }
            err = driver->inq_var(pncp->ncp, i, NULL, NULL, NULL,
                                  dimids, NULL, NULL, NULL, NULL);
            if (err != NC_NOERR) break; /* loop i */
            if (dimids[0] == pncp->unlimdimid)
                pncp->vars[i].recdim = pncp->unlimdimid;
            for (j=0; j<ndims; j++) {
                /* obtain size of dimension j */
                err = driver->inq_dim(pncp->ncp, dimids[j], NULL,
                                      pncp->vars[i].shape+j);
                if (err != NC_NOERR) break; /* loop i */
            }
        }
        if (pncp->vars[i].recdim >= 0) pncp->nrec_vars++;
    }

    if (err != NC_NOERR) { /* error happens in loop i */
        assert(i < pncp->nvars);
        for (j=0; j<=i; j++) {
            if (pncp->vars[j].shape != NULL)
                NCI_Free(pncp->vars[j].shape);
        }
        NCI_Free(pncp->vars);
        driver->close(ncp); /* close file and ignore error */
        del_from_PNCList(*ncidp);
        if (status == NC_NOERR) status = err;
    }
    if (dimids != DIMIDS) NCI_Free(dimids);

err_out:
    if (combined_info != MPI_INFO_NULL)
        MPI_Info_free(&combined_info);

    if (status != NC_NOERR && status != NC_EMULTIDEFINE_OMODE &&
        status != NC_ENULLPAD) {
        if (pncp->comm != MPI_COMM_NULL)
            MPI_Comm_free(&pncp->comm); /* a collective call */
        if (pncp->path != NULL)
            NCI_Free(pncp->path);
        NCI_Free(pncp);
        *ncidp = -1;
    }

    return status;
}

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)

int       pnc_num_aggrs_per_node;
double    pnc_ina_init;
double    pnc_ina_flatten;
double    pnc_ina_put[NTIMERS];
double    pnc_ina_get[NTIMERS];
MPI_Count pnc_ina_npairs_put;
MPI_Count pnc_ina_npairs_get;
MPI_Count pnc_ina_mem_put[NTIMERS];
MPI_Count pnc_ina_mem_get[NTIMERS];

static
void print_profiled(MPI_Comm comm)
{
    int i, rank;
    double max_t[NTIMERS];
    MPI_Count max_c[NTIMERS];

    MPI_Comm_rank(comm, &rank);

    /* print intra-node aggregation timing breakdown */
    if (pnc_num_aggrs_per_node > 0) {
        double timing, max_MiB[NTIMERS], wr_total, rd_total;
        MPI_Count count;

        wr_total = pnc_ina_init + pnc_ina_flatten;
        for (i=0; i<NTIMERS; i++) wr_total += pnc_ina_put[i];
        wr_total -= pnc_ina_put[5]; /* exclude file write time */

        rd_total = pnc_ina_init + pnc_ina_flatten;
        for (i=0; i<NTIMERS; i++) rd_total += pnc_ina_get[i];
        rd_total -= pnc_ina_get[2]; /* exclude file read time */

        MPI_Reduce(&pnc_ina_init, &timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        pnc_ina_init = timing;

        MPI_Reduce(&pnc_ina_flatten, &timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        pnc_ina_flatten = timing;

        MPI_Reduce(&wr_total, &timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        wr_total = timing;

        MPI_Reduce(&rd_total, &timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        rd_total = timing;

        MPI_Reduce(&pnc_ina_npairs_put, &count, 1, MPI_COUNT, MPI_MAX, 0, comm);
        pnc_ina_npairs_put = count;

        MPI_Reduce(pnc_ina_put, max_t, NTIMERS, MPI_DOUBLE, MPI_MAX, 0, comm);
        for (i=0; i<NTIMERS; i++) pnc_ina_put[i] = max_t[i];

        MPI_Reduce(pnc_ina_mem_put, max_c, NTIMERS, MPI_COUNT, MPI_MAX, 0, comm);
        for (i=0; i<NTIMERS; i++) max_MiB[i] = (float)max_c[i] / 1048576.0;

        if (rank == 0 && pnc_ina_npairs_put > 0) {
            printf("PNC INA put npairs=%lld mem=%.1f %.1f %.1f %.1f %.1f %.1f (MiB)\n",
                   pnc_ina_npairs_put,
                   max_MiB[0],max_MiB[1],max_MiB[2],max_MiB[3],max_MiB[4],max_MiB[5]);
            printf("PNC INA put time: init %.2f flat %.2f MD %.2f sort %.2f post %.2f wait %.2f setview %.2f total %.2f (write %.2f)\n",
                   pnc_ina_init,pnc_ina_flatten,
                   pnc_ina_put[0],pnc_ina_put[1],pnc_ina_put[2],pnc_ina_put[3],pnc_ina_put[4],
                   wr_total,pnc_ina_put[5]);
        }

        MPI_Reduce(&pnc_ina_npairs_get, &count, 1, MPI_COUNT, MPI_MAX, 0, comm);
        pnc_ina_npairs_get = count;

        MPI_Reduce(pnc_ina_get, max_t, NTIMERS, MPI_DOUBLE, MPI_MAX, 0, comm);
        for (i=0; i<NTIMERS; i++) pnc_ina_get[i] = max_t[i];

        MPI_Reduce(pnc_ina_mem_get, max_c, NTIMERS, MPI_COUNT, MPI_MAX, 0, comm);
        for (i=0; i<NTIMERS; i++) max_MiB[i] = (float)max_c[i] / 1048576.0;

        if (rank == 0 && pnc_ina_npairs_get > 0) {
            printf("PNC INA get npairs=%lld mem=%.1f %.1f %.1f %.1f %.1f %.1f (MiB)\n",
                   pnc_ina_npairs_get,
                   max_MiB[0],max_MiB[1],max_MiB[2],max_MiB[3],max_MiB[4],max_MiB[5]);
            printf("PNC INA get time: init %.2f flat %.2f MD %.2f sort %.2f post %.2f wait %.2f total %.2f (read %.2f)\n",
                   pnc_ina_init,pnc_ina_flatten,
                   pnc_ina_get[0],pnc_ina_get[1],pnc_ina_get[3],pnc_ina_get[4],
                   rd_total,pnc_ina_get[2]);
        }
    }
}
#endif

/*----< ncmpi_close() >------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_close(int ncid)
{
    int i, err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_close() */
    err = pncp->driver->close(pncp->ncp);

    /* Remove from the PNCList, even if err != NC_NOERR */
    del_from_PNCList(ncid);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    print_profiled(pncp->comm);
#endif

    /* free the PNC object */
    MPI_Comm_free(&pncp->comm); /* a collective call */

    NCI_Free(pncp->path);
    for (i=0; i<pncp->nvars; i++)
        if (pncp->vars[i].shape != NULL)
            NCI_Free(pncp->vars[i].shape);
    if (pncp->vars != NULL)
        NCI_Free(pncp->vars);
    NCI_Free(pncp);

    return err;
}

/*----< ncmpi_enddef() >-----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_enddef(int ncid) {
    int err=NC_NOERR;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (!(pncp->flag & NC_MODE_DEF)) DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)

    if (pncp->flag & NC_MODE_SAFE) { /* safe mode */
        int minE, mpireturn;
        /* check the error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }
    else if (err != NC_NOERR) return err; /* fatal error */

    /* calling the subroutine that implements ncmpi_enddef() */
    err = pncp->driver->enddef(pncp->ncp);
    if (err != NC_NOERR) return err;

    fClr(pncp->flag, NC_MODE_INDEP); /* default enters collective data mode */
    fClr(pncp->flag, NC_MODE_DEF);
    return NC_NOERR;
}

/*----< ncmpi__enddef() >----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi__enddef(int        ncid,
              MPI_Offset h_minfree,
              MPI_Offset v_align,
              MPI_Offset v_minfree,
              MPI_Offset r_align)
{
    int err=NC_NOERR;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (!(pncp->flag & NC_MODE_DEF)) {
        DEBUG_ASSIGN_ERROR(err, NC_ENOTINDEFINE)
        goto err_check;
    }

    if (h_minfree < 0 || v_align < 0 || v_minfree < 0 || r_align < 0) {
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
        goto err_check;
    }

err_check:
    if (pncp->flag & NC_MODE_SAFE) { /* safe mode */
        int minE, mpireturn;
        MPI_Offset root_args[4];

        /* first check the error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;

        /* check if h_minfree, v_align, v_minfree, and r_align are consistent
         * among all processes */
        root_args[0] = h_minfree;
        root_args[1] = v_align;
        root_args[2] = v_minfree;
        root_args[3] = r_align;
        TRACE_COMM(MPI_Bcast)(&root_args, 4, MPI_OFFSET, 0, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");

        if (root_args[0] != h_minfree ||
            root_args[1] != v_align   ||
            root_args[2] != v_minfree ||
            root_args[3] != r_align)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

        /* find min error code across processes */
        TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, pncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
        if (minE != NC_NOERR) return minE;
    }
    else if (err != NC_NOERR) return err; /* fatal error */

    /* calling the subroutine that implements ncmpi__enddef() */
    err = pncp->driver->_enddef(pncp->ncp, h_minfree, v_align,
                                           v_minfree, r_align);
    if (err != NC_NOERR) return err;

    fClr(pncp->flag, NC_MODE_INDEP); /* default enters collective data mode */
    fClr(pncp->flag, NC_MODE_DEF);
    return NC_NOERR;
}

/*----< ncmpi_redef() >------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_redef(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (fIsSet(pncp->flag, NC_MODE_RDONLY)) /* read-only */
        DEBUG_RETURN_ERROR(NC_EPERM)
    /* if open mode is inconsistent, then this return might cause parallel
     * program to hang */

    /* cannot be in define mode, must enter from data mode */
    if (fIsSet(pncp->flag, NC_MODE_DEF)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* calling the subroutine that implements ncmpi_redef() */
    err = pncp->driver->redef(pncp->ncp);
    if (err != NC_NOERR) return err;

    fSet(pncp->flag, NC_MODE_DEF);
    return NC_NOERR;
}

/*----< ncmpi_sync() >-------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpi_sync(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_sync() */
    return pncp->driver->sync(pncp->ncp);
}

/*----< ncmpi_flush() >-------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpi_flush(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_flush() */
    return pncp->driver->flush(pncp->ncp);
}

/*----< ncmpi_abort() >------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_abort(int ncid)
{
    int i, err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_abort() */
    err = pncp->driver->abort(pncp->ncp);

    /* Remove from the PNCList, even if err != NC_NOERR */
    del_from_PNCList(ncid);

    /* free the PNC object */
    MPI_Comm_free(&pncp->comm); /* a collective call */

    NCI_Free(pncp->path);
    for (i=0; i<pncp->nvars; i++)
        if (pncp->vars[i].shape != NULL)
            NCI_Free(pncp->vars[i].shape);
    if (pncp->vars != NULL)
        NCI_Free(pncp->vars);
    NCI_Free(pncp);

    return err;
}

/*----< ncmpi_set_fill() >---------------------------------------------------*/
/* This is a collective subroutine.
 * This subroutine serves both purposes of setting and inquiring the fill mode.
 */
int
ncmpi_set_fill(int  ncid,
               int  fill_mode,     /* mode to be changed by user */
               int *old_fill_mode) /* current fill mode */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (fIsSet(pncp->flag, NC_MODE_RDONLY)) /* read-only */
        DEBUG_RETURN_ERROR(NC_EPERM)

    /* not allowed to call in data mode for classic formats */
    if ((pncp->format != NC_FORMAT_NETCDF4) && !(pncp->flag & NC_MODE_DEF))
        DEBUG_RETURN_ERROR(NC_ENOTINDEFINE)

    /* calling the subroutine that implements ncmpi_set_fill() */
    err = pncp->driver->set_fill(pncp->ncp, fill_mode, old_fill_mode);
    if (err != NC_NOERR) return err;

    if (fill_mode == NC_FILL)
        fSet(pncp->flag, NC_MODE_FILL);
    else /* NC_NOFILL */
        fClr(pncp->flag, NC_MODE_FILL);

    return NC_NOERR;
}

/*----< ncmpi_inq_format() >-------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_format(int  ncid,
                 int *formatp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (formatp != NULL) *formatp = pncp->format;

    return NC_NOERR;
}

#ifdef ENABLE_ADIOS
static void swap_64(void *data)
{
    uint64_t *dest = (uint64_t*) data;
    uint64_t tmp;
    memcpy(&tmp, dest, 8);
    *dest = ((tmp & 0x00000000000000FFULL) << 56) |
            ((tmp & 0x000000000000FF00ULL) << 40) |
            ((tmp & 0x0000000000FF0000ULL) << 24) |
            ((tmp & 0x00000000FF000000ULL) <<  8) |
            ((tmp & 0x000000FF00000000ULL) >>  8) |
            ((tmp & 0x0000FF0000000000ULL) >> 24) |
            ((tmp & 0x00FF000000000000ULL) >> 40) |
            ((tmp & 0xFF00000000000000ULL) >> 56);
}

static int adios_parse_endian(char *footer, int *diff_endianness) {
    unsigned int version;
    unsigned int test = 1; /* If high bit set, big endian */

    version = ntohl (*(uint32_t *) (footer + BP_MINIFOOTER_SIZE - 4));
    char *v = (char *) (&version);
    if ((*v && !*(char *) &test) /* Both writer and reader are big endian */
        || (!*(v+3) && *(char *) &test)){ /* Both are little endian */
        *diff_endianness = 0; /* No need to change endianness */
    }
    else{
        *diff_endianness = 1;
    }

    return 0;
}
#endif

/*----< ncmpi_inq_file_format() >--------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_file_format(const char *filename,
                      int        *formatp) /* out */
{
    const char *cdf_signature="CDF";
    const char *hdf5_signature="\211HDF\r\n\032\n";
    const char *path;
    char signature[8];
    int fd;
    ssize_t rlen;

    if (formatp == NULL) return NC_NOERR;

    *formatp = NC_FORMAT_UNKNOWN;

    /* remove the file system type prefix name if there is any.  For example,
     * when filename = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * path pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * open() below
     */
    path = ncmpii_remove_file_system_type_prefix(filename);

    /* must include config.h on 32-bit machines, as AC_SYS_LARGEFILE is called
     * at the configure time and it defines _FILE_OFFSET_BITS to 64 if large
     * file feature is supported.
     */
    if ((fd = open(path, O_RDONLY, 00400)) == -1) { /* open for read */
             if (errno == ENOENT)       DEBUG_RETURN_ERROR(NC_ENOENT)
        else if (errno == EACCES)       DEBUG_RETURN_ERROR(NC_EACCESS)
        else if (errno == ENAMETOOLONG) DEBUG_RETURN_ERROR(NC_EBAD_FILE)
        else {
            fprintf(stderr,"Error on opening file %s (%s)\n",
                    filename,strerror(errno));
            DEBUG_RETURN_ERROR(NC_EFILE)
        }
    }
    /* get first 8 bytes of file, which contains the file signature */
    errno = 0;
    rlen = read(fd, signature, 8);
    if (rlen != 8) {
        close(fd); /* ignore error */
        if (rlen == 0 && errno == 0)
            fprintf(stderr, "Error in %s at %d: empty file %s\n",
                    __func__,__LINE__,filename);
        else
            fprintf(stderr, "Error in %s at %d: fail to read signature of file %s\n",
                    __func__,__LINE__,filename);
        DEBUG_RETURN_ERROR(NC_EFILE)
    }
    if (close(fd) == -1)
        DEBUG_RETURN_ERROR(NC_EFILE)

    if (memcmp(signature, cdf_signature, 3) == 0) {
             if (signature[3] == 5)  *formatp = NC_FORMAT_CDF5;
        else if (signature[3] == 2)  *formatp = NC_FORMAT_CDF2;
        else if (signature[3] == 1)  *formatp = NC_FORMAT_CLASSIC;
    }

    /* check if the file is an HDF5. */
    if (*formatp == NC_FORMAT_UNKNOWN) {
        /* The HDF5 superblock is located by searching for the HDF5 format
         * signature at byte offset 0, byte offset 512, and at successive
         * locations in the file, each a multiple of two of the previous
         * location; in other words, at these byte offsets: 0, 512, 1024, 2048,
         * and so on. The space before the HDF5 superblock is referred as to
         * "user block".
         */
        off_t offset=0;

        fd = open(path, O_RDONLY, 00400); /* error check already done */
        /* get first 8 bytes of file */
        rlen = read(fd, signature, 8); /* error check already done */

        while (rlen == 8 && memcmp(signature, hdf5_signature, 8)) {
            offset = (offset == 0) ? 512 : offset * 2;
            lseek(fd, offset, SEEK_SET);
            rlen = read(fd, signature, 8);
        }
        close(fd); /* ignore error */

        if (rlen == 8) { /* HDF5 signature found */
            /* TODO: whether the file is NC_FORMAT_NETCDF4_CLASSIC is
             * determined by HDF5 attribute "_nc3_strict" which requires a call
             * to H5Aget_name(). For now, we do not distinguish
             * NC_CLASSIC_MODEL, but simply return NETCDF4 format.
             */
#ifdef ENABLE_NETCDF4
            int err, ncid;
            err = nc_open(path, NC_NOWRITE, &ncid);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
            err = nc_inq_format(ncid, formatp);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
            err = nc_close(ncid);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
#else
            *formatp = NC_FORMAT_NETCDF4;
#endif
        }
    }

#ifdef ENABLE_ADIOS
    /* check if the file is a BP. */
    if (*formatp == NC_FORMAT_UNKNOWN) {
        off_t fsize;
        int diff_endian;
        char footer[BP_MINIFOOTER_SIZE];
        off_t h1, h2, h3;

        /* test if the file footer follows BP specification */
        if ((fd = open(path, O_RDONLY, 00400)) == -1) {
                 if (errno == ENOENT)       DEBUG_RETURN_ERROR(NC_ENOENT)
            else if (errno == EACCES)       DEBUG_RETURN_ERROR(NC_EACCESS)
            else if (errno == ENAMETOOLONG) DEBUG_RETURN_ERROR(NC_EBAD_FILE)
            else {
                fprintf(stderr,"Error on opening file %s (%s)\n",
                        filename,strerror(errno));
                DEBUG_RETURN_ERROR(NC_EFILE)
            }
        }

        /* Seek to end of file */
        fsize = lseek(fd, (off_t)(-(BP_MINIFOOTER_SIZE)), SEEK_END);

        /* read footer */
        rlen = read(fd, footer, BP_MINIFOOTER_SIZE);
        if (rlen != BP_MINIFOOTER_SIZE) {
            close(fd);
            DEBUG_RETURN_ERROR(NC_EFILE)
        }
        if (close(fd) == -1) {
            DEBUG_RETURN_ERROR(NC_EFILE)
        }

        /* check endianness of file and this running system */
        adios_parse_endian(footer, &diff_endian);

        BUFREAD64(footer,      h1) /* file offset of process group index table */
        BUFREAD64(footer + 8,  h2) /* file offset of variable index table */
        BUFREAD64(footer + 16, h3) /* file offset of attribute index table */

        /* All index tables must fall within the range of file size.
         * Process group index table must comes before variable index table.
         * Variable index table must comes before attribute index table.
         */
        if (0 < h1 && h1 < fsize &&
            0 < h2 && h2 < fsize &&
            0 < h3 && h3 < fsize &&
            h1 < h2 && h2 < h3){
            /* basic footer check is passed, now we try to open the file with
             * ADIOS library to make sure it is indeed a BP file
             */
            ADIOS_FILE *fp;
            fp = adios_read_open_file(path, ADIOS_READ_METHOD_BP,
                                        MPI_COMM_SELF);
            if (fp != NULL) {
                *formatp = NC_FORMAT_BP;
                adios_read_close(fp);
            }
        }
    }
#endif

    return NC_NOERR;
}

/*----< ncmpi_inq_version() >------------------------------------------------*/
int
ncmpi_inq_version(int ncid, int *nc_mode)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (nc_mode == NULL) return NC_NOERR;

    if (pncp->format == NC_FORMAT_CDF5)
        *nc_mode = NC_64BIT_DATA;
    else if (pncp->format == NC_FORMAT_CDF2)
        *nc_mode = NC_64BIT_OFFSET;
    else if (pncp->format == NC_FORMAT_CLASSIC)
        *nc_mode = NC_CLASSIC_MODEL;

#ifdef ENABLE_NETCDF4
    else if (pncp->format == NC_FORMAT_NETCDF4)
        *nc_mode = NC_NETCDF4;
    else if (pncp->format == NC_FORMAT_NETCDF4_CLASSIC)
        *nc_mode = NC_NETCDF4 | NC_CLASSIC_MODEL;
#endif

#ifdef ENABLE_ADIOS
    else if (pncp->format == NC_FORMAT_BP)
        *nc_mode = NC_BP;
#endif

    return NC_NOERR;
}

/*----< ncmpi_inq() >--------------------------------------------------------*/
int
ncmpi_inq(int  ncid,
          int *ndimsp,
          int *nvarsp,
          int *nattsp,
          int *xtendimp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq() */
    return pncp->driver->inq(pncp->ncp, ndimsp, nvarsp, nattsp, xtendimp);
}

/*----< ncmpi_inq_ndims() >--------------------------------------------------*/
int
ncmpi_inq_ndims(int  ncid,
                int *ndimsp)
{
    return ncmpi_inq(ncid, ndimsp, NULL, NULL, NULL);
}

/*----< ncmpi_inq_nvars() >--------------------------------------------------*/
int
ncmpi_inq_nvars(int  ncid,
                int *nvarsp)
{
    return ncmpi_inq(ncid, NULL, nvarsp, NULL, NULL);
}

/*----< ncmpi_inq_natts() >--------------------------------------------------*/
int
ncmpi_inq_natts(int  ncid,
                int *nattsp)
{
    return ncmpi_inq(ncid, NULL, NULL, nattsp, NULL);
}

/*----< ncmpi_inq_unlimdim() >-----------------------------------------------*/
int
ncmpi_inq_unlimdim(int  ncid,
                   int *unlimdimidp)
{
    return ncmpi_inq(ncid, NULL, NULL, NULL, unlimdimidp);
}

/*----< ncmpi_inq_path() >---------------------------------------------------*/
/* Get the file pathname which was used to open/create the ncid's file.
 * pathlen and path must already be allocated. Ignored if NULL.
 * This is an independent subroutine.
 */
int
ncmpi_inq_path(int   ncid,
               int  *pathlen,/* Ignored if NULL */
               char *path)   /* must have already been allocated. Ignored if NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (pathlen != NULL) {
        if (pncp->path == NULL) *pathlen = 0;
        else                    *pathlen = (int)strlen(pncp->path);
    }
    if (path != NULL) {
        if (pncp->path == NULL) *path = '\0';
        else                    strcpy(path, pncp->path);
    }
    return NC_NOERR;
}

/*----< ncmpi_inq_num_fix_vars() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_num_fix_vars(int ncid, int *num_fix_varsp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (num_fix_varsp == NULL) return NC_NOERR;

#ifdef ENABLE_NETCDF4
    if (pncp->format == NC_FORMAT_NETCDF4 ||
        pncp->format == NC_FORMAT_NETCDF4_CLASSIC) {
        /* calling the subroutine that implements ncmpi_inq_num_fix_vars() */
        return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, num_fix_varsp,
                                      NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                      NULL, NULL, NULL, NULL, NULL);
    }
#endif

    *num_fix_varsp = pncp->nvars - pncp->nrec_vars;

    /* number of fixed-size variables can also be calculated below.
    int i;
    *num_fix_varsp = 0;
    for (i=0; i<pncp->nvars; i++) {
        if (pncp->vars[i].recdim < 0)
            (*num_fix_varsp)++;
    }
    */

    return NC_NOERR;
}

/*----< ncmpi_inq_num_rec_vars() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_num_rec_vars(int ncid, int *num_rec_varsp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (num_rec_varsp == NULL) return NC_NOERR;

#ifdef ENABLE_NETCDF4
    if (pncp->format == NC_FORMAT_NETCDF4 ||
        pncp->format == NC_FORMAT_NETCDF4_CLASSIC) {
        /* calling the subroutine that implements ncmpi_inq_num_rec_vars() */
        return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL,
                                      num_rec_varsp, NULL, NULL, NULL, NULL,
                                      NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        }
#endif

    *num_rec_varsp = pncp->nrec_vars;

    /* number of record variables can also be calculated below.
    int i;
    *num_rec_varsp = 0;
    for (i=0; i<pncp->nvars; i++) {
        if (pncp->vars[i].recdim >= 0)
            (*num_rec_varsp)++;
    }
    */

    return NC_NOERR;
}

/*----< ncmpi_inq_striping() >-----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_striping(int ncid, int *striping_size, int *striping_count)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_striping() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  striping_size, striping_count, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_header_size() >--------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_header_size(int ncid, MPI_Offset *header_size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (header_size == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_header_size() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, header_size, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_header_extent() >------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_header_extent(int ncid, MPI_Offset *header_extent)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (header_extent == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_header_extent() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, header_extent, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_recsize() >------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_recsize(int ncid, MPI_Offset *recsize)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (recsize == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_recsize() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, recsize, NULL,
                                  NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_put_size() >-----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_put_size(int ncid, MPI_Offset *put_size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (put_size == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_put_size() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, put_size,
                                  NULL, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_get_size() >-----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_get_size(int ncid, MPI_Offset *get_size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (get_size == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_get_size() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  get_size, NULL, NULL, NULL, NULL);
}

/*----< ncmpi_inq_file_info() >----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_file_info(int ncid, MPI_Info *info)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (info == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_file_info() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, info, NULL, NULL, NULL);
}

/* ncmpi_get_file_info() is now deprecated, replaced by ncmpi_inq_file_info() */
int
ncmpi_get_file_info(int ncid, MPI_Info *info)
{
    return ncmpi_inq_file_info(ncid, info);
}

/*----< ncmpi_begin_indep_data() >-------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_begin_indep_data(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_begin_indep_data() */
    err = pncp->driver->begin_indep_data(pncp->ncp);
    if (err != NC_NOERR) return err;

    fSet(pncp->flag, NC_MODE_INDEP);
    return NC_NOERR;
}

/*----< ncmpi_end_indep_data() >---------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_end_indep_data(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_end_indep_data() */
    err = pncp->driver->end_indep_data(pncp->ncp);
    if (err != NC_NOERR) return err;

    fClr(pncp->flag, NC_MODE_INDEP);
    return NC_NOERR;
}

/*----< ncmpi_sync_numrecs() >-----------------------------------------------*/
/* this API is collective, but can be called in independent data mode.
 * Note numrecs (number of records) is always sync-ed in memory and file in
 * collective data mode.
 */
int
ncmpi_sync_numrecs(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_sync_numrecs() */
    return pncp->driver->sync_numrecs(pncp->ncp);
}

/*----< ncmpi_set_default_format() >-----------------------------------------*/
/* This function sets a default create file format.
 * Valid formats are NC_FORMAT_CLASSIC, NC_FORMAT_CDF2, and NC_FORMAT_CDF5
 * This API is NOT collective, as there is no way to check against an MPI
 * communicator. It should be called by all MPI processes that intend to
 * create a file later. Consistency check will have to be done in other APIs.
 */
int
ncmpi_set_default_format(int format, int *old_formatp)
{
    int err=NC_NOERR, perr=0;

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_lock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_lock")
#endif

    /* Return existing format if desired. */
    if (old_formatp != NULL)
        *old_formatp = ncmpi_default_create_format;

    /* Make sure only valid format is set. */
    if (format != NC_FORMAT_CLASSIC &&
        format != NC_FORMAT_CDF2 &&
        format != NC_FORMAT_NETCDF4 &&
        format != NC_FORMAT_NETCDF4_CLASSIC &&
        format != NC_FORMAT_CDF5) {
        DEBUG_ASSIGN_ERROR(err, NC_EINVAL)
    }
    else {
        ncmpi_default_create_format = format;
        err = NC_NOERR;
    }

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_unlock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_unlock")

err_out:
#endif
    return (err != NC_NOERR) ? err : perr;
}

/*----< ncmpi_inq_default_format() >-----------------------------------------*/
/* returns a value suitable for a create flag.  Will return one or more of the
 * following values OR-ed together: NC_64BIT_OFFSET, NC_CLOBBER */
int
ncmpi_inq_default_format(int *formatp)
{
    int perr=0;

    if (formatp == NULL) DEBUG_RETURN_ERROR(NC_EINVAL)

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_lock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_lock")
#endif

    *formatp = ncmpi_default_create_format;

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_unlock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_unlock")

err_out:
#endif
    return perr;
}

/*----< ncmpi_inq_files_opened() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_files_opened(int *num,    /* cannot be NULL */
                       int *ncids)  /* can be NULL */
{
    int i, perr=0;

    if (num == NULL) DEBUG_RETURN_ERROR(NC_EINVAL)

#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_lock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_lock")
#endif

    *num = pnc_numfiles;

    if (ncids != NULL) { /* ncids can be NULL */
        *num = 0;
        for (i=0; i<NC_MAX_NFILES; i++) {
            if (pnc_filelist[i] != NULL) {
                ncids[*num] = i;
                (*num)++;
            }
        }
    }
#ifdef ENABLE_THREAD_SAFE
    perr = pthread_mutex_unlock(&lock);
    CHECK_ERRNO(perr, "pthread_mutex_unlock")

err_out:
#endif
    return perr;
}

/*----< ncmpi_inq_nreqs() >--------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_nreqs(int  ncid,
                int *nreqs) /* number of pending nonblocking requests */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (nreqs == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_nreqs() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, NULL, nreqs, NULL, NULL);
}

/*----< ncmpi_inq_buffer_usage() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_buffer_usage(int         ncid,
                       MPI_Offset *usage) /* amount of space used so far */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (usage == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_buffer_usage() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, usage, NULL);
}

/*----< ncmpi_inq_buffer_size() >--------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_buffer_size(int         ncid,
                      MPI_Offset *buf_size) /* amount of space attached */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (buf_size == NULL) return NC_NOERR;

    /* calling the subroutine that implements ncmpi_inq_buffer_size() */
    return pncp->driver->inq_misc(pncp->ncp, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, buf_size);
}

/*----< ncmpi_buffer_attach() >----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_buffer_attach(int        ncid,
                    MPI_Offset bufsize) /* amount of memory space allowed for
                                           PnetCDF library to buffer the
                                           nonblocking requests */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_buffer_attach() */
    return pncp->driver->buffer_attach(pncp->ncp, bufsize);
}

/*----< ncmpi_buffer_detach() >----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_buffer_detach(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_buffer_detach() */
    return pncp->driver->buffer_detach(pncp->ncp);
}

/*----< ncmpi_delete() >-----------------------------------------------------*/
/*
 * filename: the name of the file we will remove.
 * info: MPI info object, in case underlying file system needs hints.
 *
 * This API is implemented in src/driver/ncmpio/ncmpio_file.c
 *
 */

/*----< ncmpi_wait() >-------------------------------------------------------*/
/* This API is an independent subroutine. */
int
ncmpi_wait(int  ncid,
           int  num_reqs, /* number of requests */
           int *req_ids,  /* [num_reqs]: IN/OUT */
           int *statuses) /* [num_reqs], can be NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_wait() */
    return pncp->driver->wait(pncp->ncp, num_reqs, req_ids, statuses,
                              NC_REQ_INDEP);
}

/*----< ncmpi_wait_all() >---------------------------------------------------*/
/* This API is a collective subroutine. */
int
ncmpi_wait_all(int  ncid,
               int  num_reqs, /* number of requests */
               int *req_ids,  /* [num_reqs]: IN/OUT */
               int *statuses) /* [num_reqs], can be NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_wait_all() */
    return pncp->driver->wait(pncp->ncp, num_reqs, req_ids, statuses,
                              NC_REQ_COLL);
}

/*----< ncmpi_cancel() >-----------------------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_cancel(int  ncid,
             int  num_reqs, /* number of requests */
             int *req_ids,  /* [num_reqs]: IN/OUT */
             int *statuses) /* [num_reqs], can be NULL */
{
    int err;
    PNC *pncp;

    /* check if ncid is valid.
     * For invalid ncid, we must return error now, as there is no way to
     * continue with invalid ncp. However, collective APIs might hang if this
     * error occurs only on a subset of processes
     */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_cancel() */
    return pncp->driver->cancel(pncp->ncp, num_reqs, req_ids, statuses);
}

