/*
  Copyright 2004-2006, UCAR/Unidata
  See COPYRIGHT file for copying and redistribution conditions.

  This is part of netCDF.
   
  This program also takes a long time to run - it writes some data in
  a very large file, and then reads it all back to be sure it's
  correct.

  This program is an add-on test to check very large 64-bit offset
  files (8 GB, so make sure you have the disk space!).

  $Id: large_files.c 531 2007-02-17 00:02:41Z robl $
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>

#define FILE_NAME "./large_files.nc"

void
check_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
	   (void) fprintf(stderr, "line %d of %s: %s\n", line, file, ncmpi_strerror(stat));
        exit(1);
    }
}

int
main(int argc, char **argv) {

   int  stat;			/* return status */
   int  ncid;			/* netCDF id */
   int rec, i, j, k;
   signed char x[] = {42, 21};

   /* dimension ids */
   int rec_dim;
   int i_dim;
   int j_dim;
   int k_dim;
   int n_dim;
 
#define NUMRECS 1
#define I_LEN 4104
#define J_LEN 1023
#define K_LEN 1023
#define N_LEN 2

   /* dimension lengths */
   MPI_Offset rec_len = NC_UNLIMITED;
   MPI_Offset i_len = I_LEN;
   MPI_Offset j_len = J_LEN;
   MPI_Offset k_len = K_LEN;
   MPI_Offset n_len = N_LEN;

   /* variable ids */
   int var1_id;
   int x_id;

   /* rank (number of dimensions) for each variable */
#  define RANK_var1 4
#  define RANK_x 2

   /* variable shapes */
   int var1_dims[RANK_var1];
   int x_dims[RANK_x];

    printf("\n*** Testing large files, slowly.\n");

    printf("*** Creating large file %s...", FILE_NAME);

    MPI_Init(&argc, &argv);

   /* enter define mode */
   stat = ncmpi_create(MPI_COMM_WORLD, FILE_NAME, NC_CLOBBER|NC_64BIT_OFFSET, 
		   MPI_INFO_NULL, &ncid);
   check_err(stat,__LINE__,__FILE__);
 
   /* define dimensions */
   stat = ncmpi_def_dim(ncid, "rec", rec_len, &rec_dim);
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_def_dim(ncid, "i", i_len, &i_dim);
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_def_dim(ncid, "j", j_len, &j_dim);
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_def_dim(ncid, "k", k_len, &k_dim);
   check_err(stat,__LINE__,__FILE__);
   stat = ncmpi_def_dim(ncid, "n", n_len, &n_dim);
   check_err(stat,__LINE__,__FILE__);

   /* define variables */

   var1_dims[0] = rec_dim;
   var1_dims[1] = i_dim;
   var1_dims[2] = j_dim;
   var1_dims[3] = k_dim;
   stat = ncmpi_def_var(ncid, "var1", NC_BYTE, RANK_var1, var1_dims, &var1_id);
   check_err(stat,__LINE__,__FILE__);

   x_dims[0] = rec_dim;
   x_dims[1] = n_dim;
   stat = ncmpi_def_var(ncid, "x", NC_BYTE, RANK_x, x_dims, &x_id);
   check_err(stat,__LINE__,__FILE__);
   /* don't initialize variables with fill values */
   stat = ncmpi_set_fill(ncid, NC_NOFILL, 0);
   check_err(stat,__LINE__,__FILE__);

   /* leave define mode */
   stat = ncmpi_enddef (ncid);
   check_err(stat,__LINE__,__FILE__);

   {			/* store var1 */
       int n = 0;
       static signed char var1[J_LEN][K_LEN];
       static MPI_Offset var1_start[RANK_var1] = {0, 0, 0, 0};
       static MPI_Offset var1_count[RANK_var1] = {1, 1, J_LEN, K_LEN};
       static MPI_Offset x_start[RANK_x] = {0, 0};
       static MPI_Offset x_count[RANK_x] = {1, N_LEN};
       for(rec=0; rec<NUMRECS; rec++) {
	   var1_start[0] = rec;
	   x_start[0] = rec;
	   for(i=0; i<I_LEN; i++) {
	       for(j=0; j<J_LEN; j++) {
		   for (k=0; k<K_LEN; k++) {
		       var1[j][k] = n++;
		   }
	       }
	       var1_start[1] = i;
	       stat = ncmpi_put_vara_schar_all(ncid, var1_id, var1_start, var1_count, &var1[0][0]);
	       check_err(stat,__LINE__,__FILE__);
	   }
       }
       stat = ncmpi_put_vara_schar_all(ncid, x_id, x_start, x_count, x);
       check_err(stat,__LINE__,__FILE__);
   }

   stat = ncmpi_close(ncid);
   check_err(stat,__LINE__,__FILE__);

   printf("ok\n");
   printf("*** Reading large file %s...", FILE_NAME);

   stat = ncmpi_open(MPI_COMM_WORLD, FILE_NAME, NC_NOWRITE, 
		   MPI_INFO_NULL, &ncid);
   check_err(stat,__LINE__,__FILE__);

   {			/* read var1 */
       int n = 0;
       static signed char var1[J_LEN][K_LEN];
       static MPI_Offset var1_start[RANK_var1] = {0, 0, 0, 0};
       static MPI_Offset var1_count[RANK_var1] = {1, 1, J_LEN, K_LEN};
       static MPI_Offset x_start[RANK_x] = {0, 0};
       static MPI_Offset x_count[RANK_x] = {1, N_LEN};
       for(rec=0; rec<NUMRECS; rec++) {
	   var1_start[0] = rec;
	   x_start[0] = rec;
	   for(i=0; i<I_LEN; i++) {
	       var1_start[1] = i;
	       stat = ncmpi_get_vara_schar_all(ncid, var1_id, var1_start, var1_count, &var1[0][0]);
	       check_err(stat,__LINE__,__FILE__);
	       for(j=0; j<J_LEN; j++) {
		   for (k=0; k<K_LEN; k++) {
		       if (var1[j][k] != (signed char) n) {
			   printf("Error on read, var1[%d, %d, %d, %d] = %d wrong, "
				  "should be %d !\n", rec, i, j, k, var1[j][k], (signed char) n); 
			   return 1;
		       }
		       n++;
		   }
	       }
	   }
	   ncmpi_get_vara_schar_all(ncid, x_id, x_start, x_count, x);
	   if(x[0] != 42 || x[1] != 21) {
	       printf("Error on read, x[] = %d, %d\n", x[0], x[1]);
	       return 1;
	   }
       }
   }
   stat = ncmpi_close(ncid);
   check_err(stat,__LINE__,__FILE__);

   printf("ok\n");
   printf("*** Tests successful!\n");

   /* Delete the file. */
   (void) remove(FILE_NAME);
   MPI_Finalize();
   return 0;
}
