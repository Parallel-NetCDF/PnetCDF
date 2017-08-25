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
 * This string must be made a global variable. Otherwise, it won't work
 * when compiled with optimization options, e.g. -O2
 */
char const pnetcdf_libvers[] =
        "\044Id: \100(#) PnetCDF library version "PNETCDF_VERSION" of "PNETCDF_RELEASE_DATE" $";

/* a cleaner version for running command "strings", e.g.
 * % strings libpnetcdf.a | grep "PnetCDF library version"
 * or
 * % strings a.out | grep "PnetCDF library version"
 */
char const pnetcdf_lib_vers[] = "PnetCDF library version "PNETCDF_VERSION" of "PNETCDF_RELEASE_DATE;

/* pnetcdf_libvers is slightly different from the one returned from
 * ncmpi_inq_libvers(). The string pnetcdf_libvers is for command "ident" to
 * use. People can run command ident libpnetcdf.a to obtain the version of a
 * library (or an executable built from that library). In PnetCDF case, the
 * command will print the string of pnetcdf_libvers. Command "ident' looks for
 * a specific keyword pattern and print it. See man page of ident.
 *
 * The API ncmpi_inq_libvers() below on the other hand returns a string to be
 * used by the utility tools like ncmpidump, ncmpigen, etc. Check the last line
 * of output from command "ncmpidump -v".
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

