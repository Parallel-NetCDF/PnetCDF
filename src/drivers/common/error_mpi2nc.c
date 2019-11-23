/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <mpi.h>

#include <pnetcdf.h>
#include <common.h>

/*----< ncmpii_error_mpi2nc() -----------------------------------------------*/
/* translate MPI error codes to PnetCDF/netCDF error codes */
int ncmpii_error_mpi2nc(int   mpi_errorcode, /* returned value from MPI call */
                        char *err_msg)       /* extra error message */
{
    int errorclass, errorStringLen;
    char errorString[MPI_MAX_ERROR_STRING];

    /* check for specific error codes understood by PnetCDF */

    /* When NC_NOCLOBBER is used in ioflags(cmode) for open to create,
     * netCDF requires NC_EEXIST returned if the file already exists.
     * In MPI standard 2.1, if MPI_File_open uses MPI_MODE_EXCL and the file has
     * already existed, the error class MPI_ERR_FILE_EXISTS should be returned.
     * For opening an existing file but the file does not exist, MPI 2.1
     * will return MPI_ERR_NO_SUCH_FILE
     * Note for MPI 2.1 and prior, we return MPI_ERR_IO, as these error classes
     * have not been defined.
     */
    MPI_Error_class(mpi_errorcode, &errorclass);

    if (errorclass == MPI_ERR_FILE_EXISTS)  return NC_EEXIST;
    if (errorclass == MPI_ERR_NO_SUCH_FILE) return NC_ENOENT;
    /* MPI-IO should return MPI_ERR_NOT_SAME when one or more arguments of a
     * collective MPI call are different. However, MPI-IO may not report this
     * error code correctly. For instance, some MPI-IO returns MPI_ERR_AMODE
     * instead when amode is found inconsistent. MPI_ERR_NOT_SAME can also
     * report inconsistent file name. */
    if (errorclass == MPI_ERR_NOT_SAME) return NC_EMULTIDEFINE_FNC_ARGS;
    /* MPI-IO may or may not report MPI_ERR_AMODE if inconsistent amode is
     * detected. MPI_ERR_AMODE can also indicate other conflict amode used
     * on each process. But in PnetCDF, MPI_ERR_AMODE can only be caused by
     * inconsistent file open/create mode. So, if MPI-IO returns this error
     * we are sure it is because of the inconsistent mode */
    if (errorclass == MPI_ERR_AMODE)     return NC_EMULTIDEFINE_OMODE;
    if (errorclass == MPI_ERR_READ_ONLY) return NC_EPERM;
    if (errorclass == MPI_ERR_ACCESS)    return NC_EACCESS;
    if (errorclass == MPI_ERR_BAD_FILE)  return NC_EBAD_FILE;
    if (errorclass == MPI_ERR_NO_SPACE)  return NC_ENO_SPACE;
    if (errorclass == MPI_ERR_QUOTA)     return NC_EQUOTA;

    /* other errors that currently have no corresponding PnetCDF error codes,
     * or the error class is MPI_ERR_IO (Other I/O error). For example,
     * MPI_ERR_INFO_VALUE (MPI info Value longer than MPI_MAX_INFO_VAL).
     */

    MPI_Error_string(mpi_errorcode, errorString, &errorStringLen);
    if (err_msg == NULL) err_msg = "";
#ifdef PNETCDF_DEBUG
    /* report the world rank */
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("rank %d: MPI error (%s) : %s\n", rank, err_msg, errorString);
#else
    printf("MPI error (%s) : %s\n", err_msg, errorString);
#endif

    return NC_EFILE; /* other unknown file I/O error */
}


