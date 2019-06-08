/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This file implements helper functions used by the ADIOS driver.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <ncadios_driver.h>
#include <ncadios_internal.h>

nc_type ncadios_to_nc_type(enum ADIOS_DATATYPES atype){
    switch (atype) {
        case adios_unsigned_byte:
            return NC_BYTE;
        case adios_byte:
            return NC_BYTE;
        case adios_short:
            return NC_SHORT;
        case adios_unsigned_short:
            return NC_USHORT;
        case adios_integer:
            return NC_INT;
        case adios_unsigned_integer:
            return NC_UINT;
        case adios_long:
            return NC_INT64;
        case adios_unsigned_long:
            return NC_UINT64;
        case adios_real:
            return NC_FLOAT;
        case adios_double:
            return NC_DOUBLE;
        case adios_long_double:
            return NC_DOUBLE;
        case adios_string:
            return NC_CHAR;
        case adios_complex:
#ifdef PNETCDF_DEBUG
            printf("Warning: unsupported adios type: adios_string_array\n");
            fflush(stdout);
#endif
            return NC_BYTE;
        case adios_double_complex:
#ifdef PNETCDF_DEBUG
            printf("Warning: unsupported adios type: adios_string_array\n");
            fflush(stdout);
#endif
            return NC_BYTE;
        case adios_string_array:
#ifdef PNETCDF_DEBUG
            printf("Warning: unsupported adios type: adios_string_array\n");
            fflush(stdout);
#endif
            return NC_BYTE;
        case adios_unknown:
#ifdef PNETCDF_DEBUG
            printf("Warning: unsupported adios type: adios_unknown\n");
            fflush(stdout);
#endif
            return NC_BYTE;
    }

    return NC_NAT;
}

MPI_Datatype ncadios_to_mpi_type(enum ADIOS_DATATYPES atype){
    switch (atype) {
        case adios_unsigned_byte:
            return MPI_BYTE;
        case adios_byte:
            return MPI_BYTE;
        case adios_short:
            return MPI_SHORT;
        case adios_unsigned_short:
            return MPI_UNSIGNED_SHORT;
        case adios_integer:
            return MPI_INT;
        case adios_unsigned_integer:
            return MPI_UNSIGNED;
        case adios_long:
            return MPI_LONG_LONG ;
        case adios_unsigned_long:
            return MPI_UNSIGNED_LONG_LONG ;
        case adios_real:
            return MPI_FLOAT;
        case adios_double:
            return MPI_DOUBLE;
        case adios_long_double:
            return MPI_DOUBLE;
        case adios_string:
            return MPI_CHAR;
        case adios_string_array:
            return MPI_CHAR;
        case adios_complex:
#ifdef PNETCDF_DEBUG
            printf("Warning: unsupported adios type: adios_string_array\n");
            fflush(stdout);
#endif
            return NC_BYTE;
        case adios_double_complex:
#ifdef PNETCDF_DEBUG
            printf("Warning: unsupported adios type: adios_string_array\n");
            fflush(stdout);
#endif
            return NC_BYTE;
        case adios_unknown:
#ifdef PNETCDF_DEBUG
            printf("Warning: unsupported adios type: adios_unknown\n");
            fflush(stdout);
#endif
            return NC_BYTE;
    }

    return NC_NAT;
}

MPI_Datatype ncadios_nc_to_mpi_type(nc_type atype){
    switch (atype) {
        case NC_BYTE:
            return MPI_BYTE;
        case NC_CHAR:
            return MPI_CHAR;
        case NC_SHORT:
            return MPI_SHORT;
        case NC_INT:
            return MPI_INT;
        case NC_FLOAT:
            return MPI_FLOAT;
        case NC_DOUBLE:
            return MPI_DOUBLE;
    }

    return NC_NAT;
}
