/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pnetcdf.h>

/* The const string below is for the RCS ident(1) command to find a string like
 * "\044Id: \100(#) PnetCDF library version 1.4.0 of 16 Nov 2013 $"
 * in the library file (libpnetcdf.a).
 *
 * This string must be made a global variable. Otherwise, it won't work when
 * compiled with optimization options, e.g. -O2
 *
 * Contents of string pnetcdf_libvers is slightly different from the one to be
 * returned from ncmpi_inq_libvers(). pnetcdf_libvers is to be used for command
 * "ident" to identify the RCS keyword strings. Note command "ident' looks for
 * a specific keyword pattern and print it. See man page of ident.  One can run
 * command "ident libfoo.a" to obtain the version string of library named foo
 * (or an executable compiled from that library). In the PnetCDF case, the
 * command "ident libpnetcdf.so" will print the contents of pnetcdf_libvers.
 */
char const pnetcdf_libvers[128] =
        "\044Id: \100(#) PnetCDF library version "PNETCDF_VERSION" of "PNETCDF_RELEASE_DATE" $";


/* Below pnetcdf_lib_vers is a cleaner version of pnetcdf_libvers. It is for
 * running command "strings", e.g.
 *    % strings libpnetcdf.a | grep "^PnetCDF library version"
 * or
 *    % strings a.out | grep "^PnetCDF library version"
 */
char const pnetcdf_lib_vers[128] =
        "PnetCDF library version "PNETCDF_VERSION" of "PNETCDF_RELEASE_DATE;

/* The API ncmpi_inq_libvers() below returns the contents of pnetcdf_lib_vers,
 * instead of pnetcdf_libvers, which is used by the PnetCDF utility tools like
 * ncmpidump and ncmpidiff. Check the last line from the output of commands
 * "ncmpidump -h" and "ncmpidiff -h".
 */

/*----< ncmpi_inq_libvers() >------------------------------------------------*/
const char*
ncmpi_inq_libvers(void) {

    /* match the style used by netCDF API nc_inq_libvers()
     * for example, "4.3.0 of Jun 16 2013 12:11:30 $"
     * we need some silly operation so the compiler will emit the otherwise
     * unused pnetcdf_libvers
     */
    if ((void *)pnetcdf_libvers != (void *)ncmpi_inq_libvers) {
	; /* do nothing */
    }
    return PNETCDF_VERSION " of " PNETCDF_RELEASE_DATE;
}

