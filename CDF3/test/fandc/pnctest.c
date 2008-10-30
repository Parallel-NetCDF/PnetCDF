#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <pnetcdf.h>

/* Test program thanks to From: John Tannahill <tannahill1@llnl.gov> */


int main (int argc, char *argv[]) {

  const int TOTSIZ_3D[3] = { 10, 20, 30 };


  int reorder = 0;

  int isperiodic[3] = { 0, 0, 0 };

  MPI_Comm comm_cart;
  int ierr;
  int lat_id, lev_id, lon_id;
  int ncid;
  int totpes;
  int tt_id;

  int dim_id[3];

  int numpes[3] = { 0, 1, 1 };  /* number of PEs along axes;
                                   determined by MPI where a
                                   zero is specified */


  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &totpes);

  MPI_Dims_create (totpes, 3, numpes);
  MPI_Cart_create (MPI_COMM_WORLD, 3, numpes, isperiodic, reorder, &comm_cart);


  ierr = ncmpi_create (comm_cart, "pnc_test.nc", NC_CLOBBER, MPI_INFO_NULL,
                       &ncid);

  ierr = ncmpi_def_dim (ncid, "level",     (MPI_Offset) TOTSIZ_3D[0], &lev_id);
  ierr = ncmpi_def_dim (ncid, "latitude",  (MPI_Offset) TOTSIZ_3D[1], &lat_id);
  ierr = ncmpi_def_dim (ncid, "longitude", (MPI_Offset) TOTSIZ_3D[2], &lon_id);

  dim_id[0] = lev_id;
  dim_id[1] = lat_id;
  dim_id[2] = lon_id;

  ierr = ncmpi_def_var (ncid, "tt", NC_FLOAT, 3, dim_id, &tt_id);

  ierr = ncmpi_enddef (ncid);

  ierr = ncmpi_close (ncid);


  MPI_Comm_free (&comm_cart);
  MPI_Finalize ( );
  return 0;
}
