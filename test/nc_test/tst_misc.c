/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* This program is based on the test program tst_misc.c of the netCDF package */

/*
  Copyright 2007, UCAR/Unidata
  See COPYRIGHT file for copying and redistribution conditions.

  This is part of netCDF.

  This program runs some extra tests.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

int
main(int argc, char **argv)
{
    char *cmd_str, *path, filename[256];
    int rank, nprocs, err, nerrs=0, ncid;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (rank > 0) goto fn_exit;
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    cmd_str = (char*)malloc(strlen(argv[0]) + 256);
    sprintf(cmd_str, "*** TESTING C   %s for emulating netCDF t_misc ", basename(argv[0]));
    if (rank == 0) printf("%-66s ------ ", cmd_str);
    free(cmd_str);

    /* remove the file system type prefix name if there is any.  For example,
     * when filename = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * path pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * fopen() below
     */
    path = remove_file_system_type_prefix(filename);

/*
   printf("\n*** Testing some extra stuff.\n");
   printf("*** Trying to open non-netCDF files of tiny length...");
*/
    {
#define DATA_LEN 32
      char dummy_data[DATA_LEN];
      int i, openstat;
      FILE *file;

      /* Appease valgrind by initializing our data. */
      for (i = 0; i < DATA_LEN; i++)
	 dummy_data[i] = i;

      for (i = DATA_LEN; i >= 0; i--)
      {
	 /* Create a small file which is not a netCDF file. */
	 if (!(file = fopen(path, "w+"))) nerrs++;
         else {
	     if (fwrite(dummy_data, 1, i, file) != i) nerrs++;
	     if (fclose(file)) nerrs++;
         }

	 /* Make sure that netCDF rejects this file politely. */
	 openstat = ncmpi_open(MPI_COMM_SELF, path, NC_NOWRITE, MPI_INFO_NULL, &ncid);
	 /* Some platforms (OSX, buddy) return stat = 2 (file not found)
	    for index i == 2.  Not sure why, but this is a work around. */
	 if(openstat != NC_ENOTNC && openstat != NC_ENOENT && openstat != NC_EFILE) {
            /* older version of OpenMPI and MPICH may return MPI_ERR_IO instead of MPI_ERR_NO_SUCH_FILE */
            printf("Error %s line %d: Expecting error code %d or %d but got %d\n",
                   __FILE__,__LINE__,NC_ENOTNC,NC_ENOENT,openstat);
            nerrs++;
         }
      }
   }

   err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL,
                      &ncid); CHECK_ERR
   err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}
