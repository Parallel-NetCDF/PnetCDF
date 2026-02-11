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
/* This subroutine is a collective call. It finds the affinity of each MPI
 * process to the compute node and returns the followings:
 *   num_nodes_ptr  Number of unique nodes (host names)
 *   node_ids_ptr   [nprocs] node IDs of each rank, must be freed by caller.
 */
int
ncmpii_construct_node_list(MPI_Comm   comm,
                           int       *num_nodes_ptr, /* OUT: */
                           int      **node_ids_ptr)  /* OUT: [nprocs] */
{
    char my_procname[MPI_MAX_PROCESSOR_NAME], **all_procnames=NULL;
    int i, j, k, rank, nprocs, num_nodes, my_procname_len, root=0;
    int *node_ids=NULL, *all_procname_lens=NULL;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    /* Collect host name of alocated compute nodes. Note my_procname is null
     * character terminated, but my_procname_len does not include the null
     * character.
     */
    MPI_Get_processor_name(my_procname, &my_procname_len);
#if 0
#ifdef MIMIC_LUSTRE
#define MIMIC_NUM_NODES 1
    /* mimic number of compute nodes = MIMIC_NUM_NODES */
    int node_id, np_per_node = nprocs / MIMIC_NUM_NODES;
    if (nprocs % MIMIC_NUM_NODES > 0) np_per_node++;
    if (rank < np_per_node * (nprocs % MIMIC_NUM_NODES))
        node_id = rank / np_per_node;
    else
        node_id = (rank - np_per_node * (nprocs % MIMIC_NUM_NODES)) / (nprocs / MIMIC_NUM_NODES) + (nprocs % MIMIC_NUM_NODES);

    sprintf(my_procname,"compute.node.%d", node_id);
    my_procname_len = (int)strlen(my_procname);
#endif
#endif

    my_procname_len++; /* to include terminate null character */

    if (rank == root) {
        /* root collects all procnames */
        all_procnames = (char **) NCI_Malloc(sizeof(char*) * nprocs);
        if (all_procnames == NULL)
            DEBUG_RETURN_ERROR(NC_ENOMEM)

        all_procname_lens = (int *) NCI_Malloc(sizeof(int) * nprocs);
        if (all_procname_lens == NULL) {
            NCI_Free(all_procnames);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }
    }
    /* gather process name lengths from all processes first */
    MPI_Gather(&my_procname_len, 1, MPI_INT, all_procname_lens, 1, MPI_INT,
               root, comm);

    if (rank == root) {
        int *disp;
        size_t alloc_size = 0;

        for (i=0; i<nprocs; i++)
            alloc_size += all_procname_lens[i];

        all_procnames[0] = (char *) NCI_Malloc(alloc_size);
        if (all_procnames[0] == NULL) {
            NCI_Free(all_procname_lens);
            NCI_Free(all_procnames);
            DEBUG_RETURN_ERROR(NC_ENOMEM)
        }

        /* Construct displacement array for the MPI_Gatherv, as each process
         * may have a different length for its process name.
         */
        disp = (int *) NCI_Malloc(sizeof(int) * nprocs);
        disp[0] = 0;
        for (i=1; i<nprocs; i++) {
            all_procnames[i] = all_procnames[i - 1] + all_procname_lens[i - 1];
            disp[i] = disp[i - 1] + all_procname_lens[i - 1];
        }

        /* gather all process names */
        MPI_Gatherv(my_procname, my_procname_len, MPI_CHAR,
                    all_procnames[0], all_procname_lens, disp, MPI_CHAR,
                    root, comm);

        NCI_Free(disp);
        NCI_Free(all_procname_lens);
    } else
        /* send process name to root */
        MPI_Gatherv(my_procname, my_procname_len, MPI_CHAR,
                    NULL, NULL, NULL, MPI_CHAR, root, comm);

    /* node_ids is an array storing the compute node IDs of all MPI processes
     * in the MPI communicator supplied by the application program. Here, we
     * use malloc() instead of NCI_Malloc, because node_ids will be freed when
     * the communicator is freed. When communicator is MPI_COMM_WORLD or
     * MPI_COMM_SELF, it is freed at MPI_Finalize() whose calls to free()
     * cannot be tracked by PnetCDF.
     */
    node_ids = (int *) malloc(sizeof(int) * (nprocs + 1));

    if (rank == root) {
        /* all_procnames[] can tell us the number of nodes and number of
         * processes per node.
         */
        char **node_names;
        int last;

        /* array of pointers pointing to unique host names (compute nodes) */
        node_names = (char **) NCI_Malloc(sizeof(char*) * nprocs);

        /* calculate node_ids[] */
        last = 0;
        num_nodes = 0; /* number of unique compute nodes */
        for (i=0; i<nprocs; i++) {
            k = last;
            for (j=0; j<num_nodes; j++) {
                /* check if [i] has already appeared in [] */
                if (!strcmp(all_procnames[i], node_names[k])) { /* found */
                    node_ids[i] = k;
                    break;
                }
                k = (k == num_nodes - 1) ? 0 : k + 1;
            }
            if (j < num_nodes)  /* found, next iteration, start with node n */
                last = k;
            else {      /* not found, j == num_nodes, add a new node */
                node_names[j] = strdup(all_procnames[i]);
                node_ids[i] = j;
                last = j;
                num_nodes++;
            }
        }
        /* num_nodes is now the number of compute nodes (unique node names) */

        for (i=0; i<num_nodes; i++)
            free(node_names[i]); /* allocated by strdup() */
        NCI_Free(node_names);
        NCI_Free(all_procnames[0]);
        NCI_Free(all_procnames);

        /* piggyback num_nodes to MPI_Bcast */
        node_ids[nprocs] = num_nodes;
    }

    /* broadcast compute node IDs of each MPI process */
    MPI_Bcast(node_ids, nprocs+1, MPI_INT, root, comm);

    *node_ids_ptr = node_ids;
    *num_nodes_ptr = node_ids[nprocs];

    return NC_NOERR;
}

