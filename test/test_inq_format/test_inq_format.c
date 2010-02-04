#include <stdio.h>
#include <mpi.h>
#include <pnetcdf.h>
#include <stdlib.h>

/* Prototype for functions used only in this file */
int main(int argc, char **argv) {

  int status;
  int ncid;
  int format;
    

  MPI_Init(&argc, &argv);

  status = ncmpi_open(MPI_COMM_WORLD, "../data/test_int.nc", 0, MPI_INFO_NULL, &ncid);
 
  status = ncmpi_inq_format(ncid, &format);
  printf("../data/test.nc format:%d\n", format);

  status = ncmpi_close(ncid);
  
  status = ncmpi_inq_file_format("../data/test_int_cdf5.nc", &format);
  printf("../data/test_int_cdf5.nc format:%d\n", format);

  MPI_Finalize();
  return 0;
}
