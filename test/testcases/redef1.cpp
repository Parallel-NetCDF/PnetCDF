#include <pnetcdf.h>
#include <vector>
#include <iostream>
#include <stdlib.h>

void PNCDF_Error(const int status, const char *const msg)  {  
  if (status != NC_NOERR) {
    std::cerr << msg << std::endl;
    exit(status);
  }  
}

char dset_name[] = "redef1.nc";

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);


  int ncid;
  int cmode = NC_64BIT_OFFSET;
  MPI_Comm comm = MPI_COMM_WORLD;

  int commsize;
  MPI_Comm_size(comm, &commsize);

  if ( commsize  > 1 ) {
    std::cerr << "run it with one process" << std::endl;
    MPI_Abort(comm, -1);
  }
  

  int status;
  status = ncmpi_create(comm, dset_name, cmode, MPI_INFO_NULL, &ncid);
  PNCDF_Error(status, "open");
  
  int dim0id;
  size_t len0 = 10;
  status = ncmpi_def_dim(ncid, "dim0", len0, &dim0id);
  PNCDF_Error(status, "def_dim");

  int dim1id;
  size_t len1 = 3;
  status = ncmpi_def_dim(ncid, "dim1", len1, &dim1id);
  PNCDF_Error(status, "def_dim");


  int dim5id;
  size_t len5 = 5;
  status = ncmpi_def_dim(ncid, "dim5", len5, &dim5id);
  PNCDF_Error(status, "def_dim5");


  int dim9id;
  size_t len9 = 9;
  status = ncmpi_def_dim(ncid, "dim9", len9, &dim9id);
  PNCDF_Error(status, "def_dim9");
  
  int varid;
  int dimsid[2] = { dim0id, dim1id };
  status = ncmpi_def_var(ncid, "xyz", NC_INT, 2, dimsid, &varid);
  PNCDF_Error(status, "def_var");

 
  int var3id;
  {
    int dimsid[2] = { dim0id, dim5id };
    status = ncmpi_def_var(ncid, "connect", NC_INT, 2, dimsid, &var3id);
    PNCDF_Error(status, "def_var3");
  }

  int var4id;
  {
    int dimsid[2] = { dim0id, dim9id };
    status = ncmpi_def_var(ncid, "connect_exterior", NC_INT, 2, dimsid, &var4id);
    PNCDF_Error(status, "def_var4");
  }



  status = ncmpi_enddef(ncid);
  PNCDF_Error(status, "enddef");  

  //put data
  MPI_Offset start[2] = {0, 0};
  MPI_Offset count[2] = {len0, len1};
  std::vector<int> data(len0*len1);
  int k=0;
  for (size_t i=0; i<len0; i++)
    for (size_t j=0; j<len1; j++)
      data[i*len1+j] = k++;
  status = ncmpi_put_vara_int_all(ncid, varid, start, count, &data[0]);
  PNCDF_Error(status, "put1");  
    
  {
    MPI_Offset count[2] = {len0, len5};
    std::vector<int> data(len0*len5);
    int k=0;
    for (size_t i=0; i<len0; i++)
      for (size_t j=0; j<len5; j++)
        data[i*len5+j] = k++;
    status = ncmpi_put_vara_int_all(ncid, var3id, start, count, &data[0]);
    PNCDF_Error(status, "put3");  
  }


  {
    MPI_Offset count[2] = {len0, len9};
    std::vector<int> data(len0*len9);
    int k=0;
    for (size_t i=0; i<len0; i++)
      for (size_t j=0; j<len9; j++)
        data[i*len9+j] = k++;
    status = ncmpi_put_vara_int_all(ncid, var4id, start, count, &data[0]);
    PNCDF_Error(status, "put4");  
  }
#if 1
  status = ncmpi_close(ncid);
  PNCDF_Error(status, "close");  

  status = ncmpi_open(comm, dset_name, NC_WRITE|NC_64BIT_OFFSET, 
		  MPI_INFO_NULL, &ncid);
#endif

  status = ncmpi_redef(ncid);
  PNCDF_Error(status, "redef");  

  int dim2id;
  size_t len2 = 10;
  status = ncmpi_def_dim(ncid, "dim2", len2, &dim2id);
  PNCDF_Error(status, "def_dim");
  
  int var2id;
  int dims2id[2] = { dim0id, dim2id };
  status = ncmpi_def_var(ncid, "xyz_r", NC_DOUBLE, 2, dims2id, &var2id);
  PNCDF_Error(status, "def_var");

  status = ncmpi_enddef(ncid);
  PNCDF_Error(status, "enddef");  


  {
    MPI_Offset start[2] = {0, 0};
    MPI_Offset count[2] = {len0, len2};
    int k=0;
    std::vector<double> data(len0*len2);
    for (size_t i=0; i<len0; i++)
      for (size_t j=0; j<len2; j++) {
        data[i*len2+j] = (k*k);
        k++;
      }
    status = ncmpi_put_vara_double_all(ncid, var2id, start, count, &data[0]);
    PNCDF_Error(status, "put2");  
  }


  status = ncmpi_close(ncid);
  PNCDF_Error(status, "close");  

  MPI_Finalize();

  return 0;
}
