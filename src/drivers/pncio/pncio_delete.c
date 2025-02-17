/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h> /* unlink() */
#endif

#include <common.h>
#include "pncio.h"

/*----< PNCIO_File_delete() >-------------------------------------------------*/
int PNCIO_File_delete(const char *filename)
{
    int err = NC_NOERR;
    char *path = ncmpii_remove_file_system_type_prefix(filename);

    err = unlink(path);
    if (err != 0)
        err = ncmpii_error_posix2nc("unlink");

    return err;
}

