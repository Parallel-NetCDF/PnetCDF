/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _NCIO_H
#define _NCIO_H

#include <stddef.h>       /* size_t */
#include <sys/types.h>    /* off_t */
#include <mpi.h>
#include "pnetcdf.h"

#define MPI_COLLECTIVE_FH 2
#define MPI_INDEPENDENT_FH 8

#define NC_collectiveFhOpened(nciop) \
    (fIsSet(nciop->mpioflags, MPI_COLLECTIVE_FH))

#define NC_independentFhOpened(nciop) \
    (fIsSet(nciop->mpioflags, MPI_INDEPENDENT_FH))

#define set_NC_collectiveFh(nciop) \
    fSet((nciop)->mpioflags, MPI_COLLECTIVE_FH)

#define set_NC_independentFh(nciop) \
    fSet((nciop)->mpioflags, MPI_INDEPENDENT_FH)


/*
 * netcdf i/o abstraction
 */
typedef struct {
    /*
     * A copy of the ioflags argument passed in to ncmpiio_open()
     * or ncmpiio_create().
     */
    int ioflags;
    int mpioflags; /* MPI_COLLECTIVE_FH or MPI_INDEPENDENT_FH */

    int mpiomode;  /* mode used in MPI_File_open, passed from collective open
                      to independent open */

    /*
     * The MPI File handle and the communicator
     */
    MPI_File collective_fh;
    MPI_File independent_fh;
    MPI_Comm comm;
    MPI_Info mpiinfo;

    int striping_unit;   /* file stripe size of the file */

    MPI_Offset put_size;  /* amount of writes committed so far in bytes */
    MPI_Offset get_size;  /* amount of reads  committed so far in bytes */

    /*
     * A copy of the 'path' argument passed in to ncmpiio_open()
     * or ncmpiio_create(). Used by ncabort() to remove (unlink)
     * the file and by error messages.
     */
    char *path;
} ncio;


#endif /* _NCIO_H */
