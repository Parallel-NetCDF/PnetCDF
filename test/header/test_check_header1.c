/***********************************************************
 *
 * This test program writes a netCDF file using the parallel
 * netCDF library using MPI-IO. 
 *
 * The output file is: "testwrite.nc"
 *
 *  The CDL notation for the test dataset is shown below:
 *
 *    netcdf test {
 *
 *       dimensions:
 *
 *            x = 100, y = 100, z = 100, time = NC_UNLIMITED;
 *
 *
 *       variables:  // variable types, names, shapes, attributes
 *
 *            int   square(x, y);
 *                     squre: description = "2-D integer array";
 *
 *            int   cube(x,y,z);
 *
 *            int   time(time);  // coordinate & record variable
 *
 *            int   xytime(time, x, y);  // record variable
 *
 *
 *      // global attributes
 *
 *           :title = "example netCDF dataset";
 *
 *
 *      data:  // data written for variables
 *          square  = 0, 1, 2, 3,  ... , 9999;
 *          cube    = 0, 1, 2, 3,  ... , 999999;
 *	    time    = 0, 1, 2, 3,  ... , 99;    // 100 records
 *          xytime  = 0, 1, 2, 3,  ... , 9999;  // 100 records
 *   }
 *
 *
 *
 * This test uses collective APIs to write variable data and only
 * deals with integer variables. 
 *
 * This test assume # of processors = 4
 *
 **********************************************************/



#include <mpi.h>
#include <stdio.h>
#include <pnetcdf.h>
#include <string.h>
#include "testutils.h"

int main(int argc, char **argv) {
  MPI_Offset i, j, k;
  int err, nerrs=0;
  int ncid;
  int dimid1, dimid2, dimid3, udimid;
  int square_dim[2], cube_dim[3], xytime_dim[3], time_dim[1];
  MPI_Offset square_start[2], cube_start[3] = {0, 0, 0};
  MPI_Offset square_count[2] = {50, 50}, cube_count[3] = {100, 50, 50};
  MPI_Offset xytime_start[3] = {0, 0, 0};
  MPI_Offset xytime_count[3] = {100, 50, 50};
  MPI_Offset time_start[1], time_count[1] = {25};
  int square_id, cube_id, xytime_id, time_id;
  static char title[] = "example netCDF dataset";
  static char description[] = "2-D integer array";
  int data[100][50][50], buffer[100];
  int rank;
  int nprocs;
  MPI_Comm comm = MPI_COMM_WORLD;
  double TotalWriteTime;
  params opts;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) 
	  fprintf(stderr, "Testing write ... \n");
  parse_write_args(argc, argv, rank, &opts);

  MPI_Barrier(MPI_COMM_WORLD);
  TotalWriteTime = MPI_Wtime();

  /**********  START OF NETCDF ACCESS **************/

  /**
   * Create the dataset
   *   File name: "testwrite.nc"
   *   Dataset API: Collective
   */

  err = ncmpi_create(comm, opts.outfname, NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncid);
  CHECK_ERR


  /**
   * Create a global attribute:
   *    :title = "example netCDF dataset";
   */
  sprintf(title, "%s:%d of %d", title, rank, nprocs);
  printf("title:%s\n", title);
  err = ncmpi_put_att_text (ncid, NC_GLOBAL, "title", strlen(title), title);
  CHECK_ERR
  
   
  /**
   * Add 4 pre-defined dimensions:
   *   x = 100, y = 100, z = 100, time = NC_UNLIMITED
   */
  err = ncmpi_def_dim(ncid, "x", 100L, &dimid1);
  CHECK_ERR
  err = ncmpi_def_dim(ncid, "y", 100L, &dimid2);
  CHECK_ERR
  err = ncmpi_def_dim(ncid, "z", 100L, &dimid3);
  CHECK_ERR
  err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &udimid);
  CHECK_ERR

  /**
   * Define the dimensionality and then add 4 variables:
   *    square(x, y), cube(x,y,z), time(time), xytime(time, x, y)  
   */

  square_dim[0] = cube_dim[0] = xytime_dim[1] = dimid1;
  square_dim[1] = cube_dim[1] = xytime_dim[2] = dimid2;
  cube_dim[2] = dimid3;
  xytime_dim[0] = udimid;
  time_dim[0] = udimid;
  err = ncmpi_def_var (ncid, "square", NC_INT, 2, square_dim, &square_id);
  CHECK_ERR
  err = ncmpi_def_var (ncid, "cube", NC_INT, 3, cube_dim, &cube_id);
  CHECK_ERR
 // err = ncmpi_def_var (ncid, "time", NC_INT, 1, time_dim, &time_id);
  err = ncmpi_def_var (ncid, "time", NC_INT, 1, time_dim, &time_id);
  CHECK_ERR
  err = ncmpi_def_var (ncid, "xytime", NC_INT, 3, xytime_dim, &xytime_id);
  CHECK_ERR

  /**
   * Add an attribute for variable: 
   *    square: description = "2-D integer array"
   */

  err = ncmpi_put_att_text (ncid, square_id, "description",
                          strlen(description), description);
  CHECK_ERR

  /**
   * End Define Mode (switch to data mode)
   *   Dataset API: Collective
   */
  
  err = ncmpi_enddef(ncid);
  CHECK_ERR

  /**
   * Data Partition (Assume 4 processors):
   *   square: 2-D, (Block, Block), 50*50 from 100*100 
   *   cube:   3-D, (*, Block, Block), 100*50*50 from 100*100*100
   *   xytime: 3-D, (*, Block, Block), 100*50*50 from 100*100*100
   *   time:   1-D, Block-wise, 25 from 100
   */
  else {
	  square_start[0] = cube_start[1] = xytime_start[1] = (rank/2) * 50;
	  square_start[1] = cube_start[2] = xytime_start[2] = (rank%2) * 50;
	  time_start[0] = (rank%4) * 25;


  /**
   * Packing data in the buffer 
   */

  /* Data for variable: time */
	  for ( i = time_start[0]; i < time_start[0] + time_count[0]; i++ )
	    buffer[i - time_start[0]] = i;   

  /* Data for variable: square, cube and xytime */
	  for ( i = 0; i < 100; i++ )
	    for ( j = square_start[0]; j < square_start[0]+square_count[0]; j++ )
	      for ( k = square_start[1]; k < square_start[1]+square_count[1]; k++ )
	        data[i][j-square_start[0]][k-square_start[1]] = i*100*100 + j*100 + k;

  /**
   * Write data into variables: square, cube, time and xytime  
   *   Access Method: subarray
   *   Data Mode API: collective
   */ 
  
	  err = ncmpi_put_vara_int_all(ncid, square_id,
                    square_start, square_count,
                    &data[0][0][0]);
	  CHECK_ERR
	  err = ncmpi_put_vara_int_all(ncid, cube_id,
                    cube_start, cube_count,
                    &data[0][0][0]);
	  CHECK_ERR
	  err = ncmpi_put_vara_int_all(ncid, time_id,
                    time_start, time_count,
                    (void *)buffer);
	  CHECK_ERR
	  err = ncmpi_put_vara_int_all(ncid, xytime_id,
                    xytime_start, xytime_count,
                    &data[0][0][0]);
	  CHECK_ERR

/*
err = ncmpi_sync(ncid);
CHECK_ERR
err = ncmpi_redef(ncid);
CHECK_ERR
err = ncmpi_del_att(ncid, square_id, "description");
CHECK_ERR
err = ncmpi_enddef(ncid);
CHECK_ERR
*/

  /**
   * Close the dataset
   *   Dataset API:  collective
   */

	  err = ncmpi_close(ncid);
	  CHECK_ERR
  /*******************  END OF NETCDF ACCESS  ****************/

	MPI_Barrier(MPI_COMM_WORLD); 
	TotalWriteTime = MPI_Wtime() - TotalWriteTime;

	if (rank == 0) {
	  fprintf(stderr, "OK\nFile written to: %s!\n", opts.outfname);
	  fprintf(stderr, "Total Write Time = %10.8f\n", TotalWriteTime);
	}
  }
  MPI_Finalize();
  return (nerrs > 0);
}

