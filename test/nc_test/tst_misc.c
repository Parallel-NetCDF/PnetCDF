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
#include <pnetcdf.h>

#define FAIL_COLOR "\x1b[31mfail\x1b[0m\n"
#define PASS_COLOR "\x1b[32mpass\x1b[0m\n"

int
main(int argc, char **argv) 
{
    char filename[128];
    int rank, nprocs, err, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    strcpy(filename, "testfile.nc");
    if (argc == 2) strcpy(filename, argv[1]);
    if (rank > 0) goto fn_exit;

    char cmd_str[256];
    sprintf(cmd_str, "*** TESTING C   %s for emulating netCDF t_misc ", argv[0]);
    if (rank == 0) printf("%-66s ------ ", cmd_str);
/*
   printf("\n*** Testing some extra stuff.\n");
   printf("*** Trying to open non-netCDF files of tiny length...");
*/
   {
#define DATA_LEN 32    
     int ncid,openstat;
      char dummy_data[DATA_LEN];
      FILE *file;
      int i, nerrs=0;

      /* Appease valgrind by initializing our data. */
      for (i = 0; i < DATA_LEN; i++)
	 dummy_data[i] = i;

      for (i = DATA_LEN; i >= 0; i--)
      {
	 /* Create a small file which is not a netCDF file. */
	 if (!(file = fopen(filename, "w+"))) nerrs++;
	 if (fwrite(dummy_data, 1, i, file) != i) nerrs++;
	 if (fclose(file)) nerrs++;
	 
	 /* Make sure that netCDF rejects this file politely. */
	 openstat = ncmpi_open(MPI_COMM_SELF, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
	 /* Some platforms (OSX, buddy) return stat = 2 (file not found)
	    for index i == 2.  Not sure why, but this is a work around. */
	 if(openstat != NC_ENOTNC && openstat != NC_ENOENT) {
            printf("Expecting error code %d or %d but got %d\n",NC_ENOTNC,NC_ENOENT,openstat);
            nerrs++;
         }
      }
   }

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    if (rank == 0) {
        if (nerrs) printf(FAIL_COLOR);
        else       printf(PASS_COLOR);
    }

    MPI_Finalize();
    return 0;
}
