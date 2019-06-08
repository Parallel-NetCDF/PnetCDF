/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <common.h>
#include <errno.h>

/*----< ncmpii_error_posix2nc() ------------------------------------------------*/
/* translate posix io error codes to PnetCDF/netCDF error codes */
int ncmpii_error_posix2nc(char *err_msg)       /* extra error message */
{
#if defined(HAVE_STRERROR) && (HAVE_STRERROR == 1)
    char *errorString= strerror(errno);
#else
    char *errorString="Other I/O error";
#endif

    /* check for specific error codes understood by PnetCDF */
    switch (errno){
        case ENOSPC :
            return NC_ENO_SPACE;
        case ENAMETOOLONG:
        case ENOTDIR :
        case EISDIR:
            return NC_EBAD_FILE;
        case EDQUOT:
            return NC_EQUOTA;
        case ENOENT:
            return NC_ENOENT;
        case EEXIST:
            return NC_EEXIST;
    }

    /* other errors that currently have no corresponding PnetCDF error codes */
    if (err_msg == NULL) err_msg = "";
#ifdef PNETCDF_DEBUG
    /* report the world rank */
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("rank %d: IO error (%s) : %s\n", rank, err_msg, errorString);
#else
    printf("IO error (%s) : %s\n", err_msg, errorString);
#endif

    return NC_EFILE; /* other unknown file I/O error */
}

