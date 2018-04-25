#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <common.h>
#include <errno.h>

/*----< ncmpii_error_posix2nc() ------------------------------------------------*/
/* translate posix io error codes to PnetCDF/netCDF error codes */
int ncmpii_error_posix2nc(char *err_msg)       /* extra error message */
{
    int io_errorcode;
    char errorString[MPI_MAX_ERROR_STRING];

    /* check for specific error codes understood by PnetCDF */
    io_errorcode = errno;
    switch (io_errorcode){
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
    strerror_r(io_errorcode, errorString, MPI_MAX_ERROR_STRING);
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

