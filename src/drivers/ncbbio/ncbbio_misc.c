/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <pnc_debug.h>
#include <stdlib.h>
#include <string.h>
#include <common.h>

static unsigned int hash(const char* key){
    int i;
    unsigned int val = 0;

    for(i = 0; key[i] != '\0'; i++){
        val = (val * 256 + (unsigned int)key[i]) % 65535;
    }

    return val;
}

/*----< Node_local_comm() >-----------------------------------------------------*/
int ncbbio_get_node_comm(MPI_Comm global_comm, MPI_Comm *node_comm)
{
    int err;
    int i;
    int rank, np;
    int namelen, color, nnode = 0;
    char *myname, *buf, *name;
    hash_map map;

    MPI_Comm_rank(global_comm, &rank);
    MPI_Comm_size(global_comm, &np);

    hash_map_init(&map, 65535, hash);

    myname = (char*)NCI_Malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));
    buf = (char*)NCI_Malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char) * np);

    err = MPI_Get_processor_name(myname, &namelen);
    if (err != MPI_SUCCESS){
        /* Set name to empty string and participate collective call */
        memset(myname, 0, MPI_MAX_PROCESSOR_NAME * sizeof(char));
    }
    myname[namelen] = '\0';

    err = MPI_Allgather(myname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, buf, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, global_comm);
    if (err != MPI_SUCCESS){
        DEBUG_RETURN_ERROR(NC_EMPI);
    }

    for(i = 0; i < np; i++){
        name = buf + MPI_MAX_PROCESSOR_NAME * i;
        err = hash_map_add(&map, name, nnode);
        if (err == NC_NOERR){
            /* New host name, increase node number */
            nnode += 1;
        }
    }

    err = hash_map_find(&map, myname, &color);
    if (err != NC_NOERR){
        /* Our own name is not in the map, not possible */
        MPI_Abort(global_comm, -1);
    }

    err = MPI_Comm_split(global_comm, color, 0, node_comm);
    if (err != MPI_SUCCESS){
        DEBUG_RETURN_ERROR(NC_EMPI);
    }

    return NC_NOERR;
}
