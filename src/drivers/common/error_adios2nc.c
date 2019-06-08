/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
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
#include <adios_read.h>
#include <adios_error.h>

/*----< ncmpii_error_adios2nc() ------------------------------------------------*/
/* translate posix io error codes to PnetCDF/netCDF error codes */
int ncmpii_error_adios2nc(int adios_err, char *err_msg)       /* extra error message */
{
    const char *errstr;

    /* check for specific error codes understood by PnetCDF */
    switch (adios_err){
        case err_file_not_found:
        case err_invalid_file_pointer:
            return NC_EBAD_FILE;
        case err_no_memory:
            return NC_ENOMEM;
        case err_invalid_varid:
        case err_invalid_varname:
            return NC_ENOTVAR;
        case err_invalid_attrid:
        case err_invalid_attrname:
            return NC_ENOTATT;
        case err_invalid_attribute_reference:
        case err_invalid_timestep:
        case err_invalid_read_method:
        case err_invalid_group:
        case err_invalid_group_struct:
            return NC_EINVAL;
        case err_out_of_bound:
            return NC_EINVALCOORDS;
        case err_operation_not_supported:
            return NC_ENOTSUPPORT;
        case err_file_read_error:
            return NC_EREAD;
        case err_corrupted_variable:
        case err_corrupted_attribute:
            return NC_ETRUNC;
        default:
            return NC_EADIOS;
    }

    /* other errors that currently have no corresponding PnetCDF error codes */
    errstr = adios_errmsg();
    if (err_msg == NULL) err_msg = "";

#ifdef PNETCDF_DEBUG
    /* report the world rank */
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("rank %d: IO error (%s) : %s\n", rank, err_msg, errstr);
#else
    printf("IO error (%s) : %s\n", err_msg, errstr);
#endif

    return NC_EFILE; /* other unknown file I/O error */
}

