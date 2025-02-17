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
#include <unistd.h> /* fsync(), unlink(), ftruncate(), lseek() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "pncio.h"

/*----< PNCIO_File_sync() >---------------------------------------------------*/
int PNCIO_File_sync(PNCIO_File *fd)
{
    int err = NC_NOERR;

    if (fd->is_open > 0) {
        err = fsync(fd->fd_sys);
        if (err != 0)
            err = ncmpii_error_posix2nc("fsync");
    }

    return err;
}

