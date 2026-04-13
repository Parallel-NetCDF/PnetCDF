/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strlen(), strcpy(), strchr(), strncmp() */
#include <assert.h>

#include <mpi.h>

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

/*----< ncmpii_nc2mpitype() >------------------------------------------------*/
MPI_Datatype
ncmpii_nc2mpitype(nc_type xtype)
{
    switch(xtype){
        case NC_CHAR :   return MPI_CHAR;
        case NC_BYTE :   return MPI_SIGNED_CHAR;
        case NC_SHORT :  return MPI_SHORT;
        case NC_INT :    return MPI_INT;
        case NC_FLOAT :  return MPI_FLOAT;
        case NC_DOUBLE : return MPI_DOUBLE;
        case NC_UBYTE :  return MPI_UNSIGNED_CHAR;
        case NC_USHORT : return MPI_UNSIGNED_SHORT;
        case NC_UINT :   return MPI_UNSIGNED;
        case NC_INT64 :  return MPI_LONG_LONG_INT;
        case NC_UINT64 : return MPI_UNSIGNED_LONG_LONG;
        default:         return MPI_DATATYPE_NULL;
    }
}

/*----< ncmpii_xlen_nc_type() >----------------------------------------------*/
/* return the length of external NC data type */
int
ncmpii_xlen_nc_type(nc_type xtype, int *size)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  *size = 1; return NC_NOERR;
        case NC_SHORT:
        case NC_USHORT: *size = 2; return NC_NOERR;
        case NC_INT:
        case NC_UINT:
        case NC_FLOAT:  *size = 4; return NC_NOERR;
        case NC_DOUBLE:
        case NC_INT64:
        case NC_UINT64: *size = 8; return NC_NOERR;
        default: DEBUG_RETURN_ERROR(NC_EBADTYPE);
    }
}

/* File system types recognized by ROMIO in MPICH 4.0.0, and by PnetCDF */
static const char* fstypes[] = {"ufs", "nfs", "xfs", "pvfs2", "gpfs", "panfs", "lustre", "daos", "testfs", "ime", "quobyte", NULL};

/* Return a pointer to filename by removing the file system type prefix name if
 * there is any.  For example, when filename = "lustre:/home/foo/testfile.nc",
 * remove "lustre:" to return a pointer to "/home/foo/testfile.nc", so the name
 * can be used in POSIX open() calls.
 */
char* ncmpii_remove_file_system_type_prefix(const char *filename)
{
    char *ret_filename = (char*)filename;

    if (filename == NULL) return NULL;

    if (strchr(filename, ':') != NULL) { /* there is a prefix end with ':' */
        /* check if prefix is one of recognized file system types */
        int i=0;
        while (fstypes[i] != NULL) {
            size_t prefix_len = strlen(fstypes[i]);
            if (!strncmp(filename, fstypes[i], prefix_len)) { /* found */
                ret_filename += prefix_len + 1;
                break;
            }
            i++;
        }
    }

    return ret_filename;
}

/*----< ncmpii_construct_node_list() >---------------------------------------*/
/* This subroutine is a collective call. It finds the affinity of MPI processes
 * to their shared-memory compute nodes (NUMA) and returns the followings:
 *   num_NUMAs: Number of NUMA nodes
 *   numa_ids[nprocs]: node IDs of each rank, must be freed by the caller.
 *   hwcomm: sub-communicator split based on NUMA hardware. It is the caller's
 *           responsibility to free it.
 */
int
ncmpii_construct_node_list(MPI_Comm   comm,
                           int       *num_NUMAs, /* OUT: */
                           int      **numa_ids,  /* OUT: [nprocs] */
                           MPI_Comm  *hwcomm)    /* OUT: NUMA communicator */
{
    char *err_msg="No error";
    int i, err, rank, nprocs, numa_id, *ids;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    *num_NUMAs = 0;
    *numa_ids = NULL;

#if 1
    /* split comm based on NUMA nodes (processes sharing memory) */
    err = MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                              hwcomm);
#else
    /* Below code fragment is from MPI standard 4.0's example 7.3:
     * Splitting MPI_COMM_WORLD into NUMANode subcommunicators.
     */
    MPI_Info_set(info, "mpi_hw_resource_type" , "NUMANode");
    err = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_HW_GUIDED,
                              rank, info, &hwcomm);

    /* Below code fragment is from MPI standard 5.0's example 7.3:
     * Splitting MPI_COMM_WORLD into NUMANode subcommunicators.
     */
    MPI_Info_set(info, "mpi_hw_resource_type" , "hwloc://NUMANode");
    err = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_HW_GUIDED,
                              rank, info, &hwcomm);

    /* Below code fragment is from MPI standard 5.0's example 7.4:
     * Splitting MPI_COMM_WORLD into NUMANode subcommunicators.
     */
    MPI_Info_set(info, "mpi_hw_resource_type", "hwloc://NUMANode");
    err = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_RESOURCE_GUIDED,
                              rank, info, &hwcomm);
#endif

    if (err != MPI_SUCCESS) {
        err_msg = "MPI_Comm_split_type()";
        goto err_out;
    }

    if (*hwcomm == MPI_COMM_NULL) {
        err_msg = "MPI_Comm_split_type() hwcomm NULL";
        goto err_out;
    }

    /* Use hwcomm's root's rank as this process's NUMA node ID */
    numa_id = rank;
    MPI_Bcast(&numa_id, 1, MPI_INT, 0, *hwcomm);

    /* Gather all NUMA node IDs */
    *numa_ids = (int*) malloc(sizeof(int) * nprocs);
    MPI_Allgather(&numa_id, 1, MPI_INT, *numa_ids, 1, MPI_INT, comm);

    /* Count number of unique IDs and reassign NUMA ID */
    ids = (int*) calloc(nprocs, sizeof(int));
    *num_NUMAs = 0;
    for (i=0; i<nprocs; i++) {
        if (ids[(*numa_ids)[i]] == 0) {
            (*num_NUMAs)++; /* unique count */
            ids[(*numa_ids)[i]] = (*num_NUMAs); /* New ID, starting from 0 */
        }
        (*numa_ids)[i] = ids[(*numa_ids)[i]] - 1;
    }
    free(ids);

err_out:
    if (err != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(err, err_msg);

    return NC_NOERR;
}

