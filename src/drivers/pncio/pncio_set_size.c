/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* strdup() */
#include <assert.h>
#include <sys/errno.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h> /* ftruncate(), lseek() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "pncio.h"

/*----< PNCIO_File_set_size() >-----------------------------------------------*/
int PNCIO_File_set_size(PNCIO_File *fd,
                        MPI_Offset  size)
{
    int err = NC_NOERR, rank;

    MPI_Comm_rank(fd->comm, &rank);

    if (rank == 0) {
        err = ftruncate(fd->fd_sys, (off_t) size);
        if (err != 0)
            err = ncmpii_error_posix2nc("ftruncate");
    }

    MPI_Bcast(&err, 1, MPI_INT, 0, fd->comm);

    return err;
}

/*----< PNCIO_File_get_size() >-----------------------------------------------*/
int PNCIO_File_get_size(PNCIO_File *fd,
                        MPI_Offset *size)
{
    int err = NC_NOERR, rank;
    MPI_Offset msg[2];

    MPI_Comm_rank(fd->comm, &rank);

    if (rank == 0) {
        *size = lseek(fd->fd_sys, 0, SEEK_END);
        if (*size == -1)
            err = ncmpii_error_posix2nc("lseek");
        msg[0] = err;
        msg[1] = *size;
    }

    MPI_Bcast(msg, 2, MPI_OFFSET, 0, fd->comm);
    err = (int)msg[0];
    *size = msg[1];

    return err;
}

