/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header$
 *********************************************************************/

#include <netcdf.h>
#include "tests.h"

  /*
   * Test driver for netCDF implementation.  This program performs
   * tests against the netCDF specification for all user-level
   * functions in an implementation of the netCDF library.  Must be
   * invoked from a directory in which the invoker has write
   * permission.
   */

int
main()
{
    extern int ncopts;		/* netCDF error options */

    static char testfile[] = "testfile.nc";

    ncopts &= ~NC_FATAL;	/* make errors nonfatal */
    ncopts &= ~NC_VERBOSE;	/* turn off error messages */

    test_nccreate(testfile);

    test_ncopen(testfile);

    test_ncredef(testfile);

    test_ncendef(testfile);

    test_ncclose(testfile);

    test_ncinquire(testfile);

    test_ncsync(testfile);

    test_ncabort(testfile);

    test_ncdimdef(testfile);

    test_ncdimid(testfile);

    test_ncdiminq(testfile);

    test_ncdimrename(testfile);

    test_ncvardef(testfile);

    test_ncvarid(testfile);

    test_ncvarinq(testfile);

    test_ncvarput1(testfile);

    test_ncvarget1(testfile);

    test_ncvarput(testfile);

    test_ncvarget(testfile);

    test_ncvarputg(testfile);

    test_ncvargetg(testfile);

    test_ncrecinq(testfile);

    test_ncrecput(testfile);

    test_ncrecget(testfile);

    test_ncvarrename(testfile);

    test_ncattput(testfile);

    test_ncattinq(testfile);

    test_ncattget(testfile);

    test_ncattcopy(testfile, "test2.nc");

    test_ncattname(testfile);

    test_ncattrename(testfile);

    test_ncattdel(testfile);

    test_nctypelen();

#ifdef vms
#define EXIT_SUCCESS 1
#else
#define EXIT_SUCCESS 0
#endif
    return EXIT_SUCCESS;
}
