/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */
#ifndef _SUBFILE_H
#define _SUBFILE_H

#include "ncmpio_NC.h"

/* structure for storing access info of this process's request
   to the subfiles on all other processes, and vice-versa. used
   as array of structures indexed by subfile index. */
typedef struct {
    MPI_Offset *start;
    MPI_Offset *count;
    MPI_Offset *start_org;
} NC_subfile_access;

#define CEIL(x) ( (x - (int)x)==0 ? (int)x : (int)x+1 )
#define FLOOR(x) ( (x - (int)x)==0 ? (int)x : (int)x-1 )
#define ROUND(x) ( x >= 0 ? (int)(x+0.5) : (int)(x-0.5) )
#define ABS(a) (((a) < 0) ? -(a) : (a))

extern int ncmpio_subfile_open(NC *ncp);

extern int ncmpio_subfile_close(NC *ncp);

extern int ncmpio_subfile_partition(NC *ncp);

extern int ncmpio_subfile_getput_vars(NC *ncp, NC_var *varp,
           const MPI_Offset start[], const MPI_Offset count[],
           const MPI_Offset  stride[], void *buf, MPI_Offset bufcount,
           MPI_Datatype buftype, int reqMode);

#endif /* _SUBFILE_H */
