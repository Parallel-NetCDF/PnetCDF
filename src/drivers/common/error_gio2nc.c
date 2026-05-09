/*
 *  Copyright (C) 2026, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <gio.h>

#include <pnetcdf.h>
#include <common.h>

/*----< ncmpii_error_gio2nc() -----------------------------------------------*/
/* translate GIO error codes to PnetCDF/netCDF error codes */
int ncmpii_error_gio2nc(int         gio_err, /* returned value from GIO call */
                        const char *err_msg) /* extra error message */
{
    const char *dump_str = (err_msg == NULL) ? "" : err_msg;

    switch (gio_err) {
        case GIO_NOERR:         return NC_NOERR;
        case GIO_EEXIST:        return NC_EEXIST;
        case GIO_EINVAL:        return NC_EINVAL;
        case GIO_EPERM:         return NC_EPERM;
        case GIO_ENOMEM:        return NC_ENOMEM;
        case GIO_EACCESS:       return NC_EACCESS;
        case GIO_ENEGATIVECNT:  return NC_ENEGATIVECNT;
        case GIO_EFILE:         return NC_EFILE;
        case GIO_ENOENT:        return NC_ENOENT;
        case GIO_EBAD_FILE:     return NC_EBAD_FILE;
        case GIO_ENO_SPACE:     return NC_ENO_SPACE;
        case GIO_EQUOTA:        return NC_EQUOTA;
        case GIO_EFSTYPE:       return NC_EFSTYPE;
        case GIO_EMULTIDEFINE_AMODE:    return NC_EMULTIDEFINE_OMODE;
        case GIO_EMULTIDEFINE_FNC_ARGS: return NC_EMULTIDEFINE_FNC_ARGS;
        case GIO_EMULTIDEFINE_HINTS:    return NC_EMULTIDEFINE_HINTS;
        case GIO_EFILEVIEW:     return NC_EFILEVIEW;
    }

#if PNETCDF_DEBUG_MODE == 1
    /* report the world rank */
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fprintf(stderr,"rank %d: GIO error (%s) : %s\n", rank, dump_str, GIO_strerror(gio_err));
#else
    fprintf(stderr,"GIO error (%s) : %s\n", dump_str, GIO_strerror(gio_err));
#endif

    return NC_EFILE; /* other unknown file I/O error */
}


