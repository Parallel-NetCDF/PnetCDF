#include <mpi.h>
#include <stdio.h>
#include <pnetcdf.h>
#include <string.h>

#include "test_dtype.h"

#define TEST_DEFAULT_FILE "test_darray.nc"
#define TEST_DIMS 3 /* number of dimensions */
#define TEST_N 16 /* size for the shortest dimension of the memory array */
#define TEST_ARRAY_ORDER MPI_ORDER_C

void partition_array(	/* input parameters : */
		     int nprocs, int myrank, int ndims, int *array_of_sizes,
			/* output parameters : */
		     int *array_of_distribs, int *array_of_dargs,
		     MPI_Offset *local_starts, 
		     MPI_Offset *local_edges, 
		     MPI_Offset *local_strides,
		     int *array_of_psizes)
{
  int i;
  int cycle, remain_cycle;
  int pcoord;

  for (i=0; i<ndims; i++) {
   /* make it a strided-subarray access pattern */
    array_of_distribs[i] = MPI_DISTRIBUTE_CYCLIC;
    array_of_dargs[i] = MPI_DISTRIBUTE_DFLT_DARG;

   /* all dimensions will be partitioned */
    array_of_psizes[i] = 0;
  }
  
  MPI_Dims_create(nprocs, ndims, array_of_psizes);

  for (i=ndims-1; i>=0; i--) {
   /* compute the netCDF strided-subarray access parameters */
    pcoord = myrank % array_of_psizes[i];
    myrank /= array_of_psizes[i];
    local_starts[i] = (MPI_Offset)( (pcoord>=array_of_sizes[i])?0:pcoord );
    local_edges[i] = (MPI_Offset)( array_of_sizes[i]/array_of_psizes[i] );
    if ( local_edges[i]*array_of_psizes[i]+pcoord < array_of_sizes[i] )
      local_edges[i]++;
    local_strides[i] = (MPI_Offset)array_of_psizes[i];
  }
  
  return;
}

int main(int argc, char** argv) {

  int i;
  int power2;
  int ndims;
  int *array_of_sizes;
  int *array_of_distribs, *array_of_dargs;
  int *array_of_psizes;
  int ncid, *dimids, varid_1, varid_2;
  MPI_Offset *local_starts, *local_edges, *local_strides;
  char dimname[20], *filename;
  int order;
  ncmpi_type nc_etype;
  MPI_Datatype mpi_etype, mpi_darray;
  TEST_NATIVE_ETYPE *buf1, *buf2, *tmpbuf;
  void *packbuf;
  int packsz;
  int packpos;
  int total_sz, local_sz;
  int nprocs, rank;
  int status;
  int success, successall;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    printf("testing memory darray layout ...\n");
  }

 /* test initializations: nc_file, nc_variables, dimensions, arrays */

  if (argc < 2)
    filename = TEST_DEFAULT_FILE;
  else
    filename = argv[1];

  TEST_SET_NCMPI_ETYPE(nc_etype, mpi_etype);
#ifdef TEST_NCTYPE						
  nc_etype = TEST_NCTYPE;					
#endif								

  status = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, 
			MPI_INFO_NULL, &ncid);
  TEST_HANDLE_ERR(status);

  ndims = TEST_DIMS;
  order = TEST_ARRAY_ORDER;
  array_of_sizes = (int *)
		   malloc(sizeof(int)*ndims*5 + sizeof(MPI_Offset)*ndims*3);
  array_of_distribs = array_of_sizes + ndims;
  array_of_dargs = array_of_distribs + ndims;
  array_of_psizes = array_of_dargs + ndims;
  dimids = array_of_psizes + ndims;
  local_starts = (MPI_Offset *)(dimids + ndims);
  local_edges = local_starts + ndims;
  local_strides = local_edges + ndims;

  total_sz = 1;
  power2 = 1;
  for (i=0; i<ndims; i++, power2*=2) {
    total_sz *= (array_of_sizes[i] = TEST_N*power2);
    sprintf(dimname, "dim_%d", i);
    status = ncmpi_def_dim(ncid, dimname,
			   (MPI_Offset)array_of_sizes[i], dimids+i);
    TEST_HANDLE_ERR(status);
  }

  status = ncmpi_def_var(ncid, "var_1", nc_etype, ndims, dimids, &varid_1);
  TEST_HANDLE_ERR(status);

  status = ncmpi_enddef(ncid);
  TEST_HANDLE_ERR(status);

  if (rank == 0) {
    printf("\t Filesize = %2.3fMB, MAX_Memory_needed = %2.3fMB\n\n", 
	   1*total_sz*TEST_NCTYPE_LEN(nc_etype)/1024.0/1024.0, 
	   ( (2*total_sz + 4*total_sz/nprocs)*sizeof(TEST_NATIVE_ETYPE)
	   + total_sz*TEST_NCTYPE_LEN(nc_etype) )/1024.0/1024.0 );
  }

  buf1 = (TEST_NATIVE_ETYPE *)malloc(total_sz*sizeof(TEST_NATIVE_ETYPE)*2);
  buf2 = buf1 + total_sz;
  for (i=0; i<total_sz; i++)
    /* just make sure any type can represent the number */
    /* and make it irregular (cycle != power of 2) for random test */
    buf1[i] = buf2[i] = (TEST_NATIVE_ETYPE)(i%127); 

 /* PARTITION and calculate the local target region */

  partition_array(nprocs, rank, ndims, array_of_sizes, 
		  array_of_distribs, array_of_dargs,
		  local_starts, local_edges, local_strides,
		  array_of_psizes);

  local_sz = 1;
  for (i=0; i<ndims; i++) {
    local_sz *= (int)local_edges[i];
  }

 /* CREATE local darray memory view */

  MPI_Type_create_darray(nprocs, rank, ndims, array_of_sizes, 
			 array_of_distribs, array_of_dargs,
			 array_of_psizes, 
			 order, 
			 mpi_etype, 
			 &mpi_darray);
  MPI_Type_commit(&mpi_darray);

 /* PRINT stats */
  if (rank == 0) {
    printf("Initialization:  NDIMS = %d, NATIVE_ETYPE = %s, NC_TYPE = %s\n\n",
	   ndims, TEST_NATIVE_ETYPE_STR, TEST_GET_NCTYPE_STR(nc_etype));

    printf("\t NC Var_1 Shape:\t [");
    TEST_PRINT_LIST(array_of_sizes, 0, ndims-1, 1);
    printf("]\n");

    printf("\t Memory Array Shape:\t [");
    TEST_PRINT_LIST(array_of_sizes, 0, ndims-1, 1);
    printf("]\n");
    printf("\t Memory Array Copys: buf1 for write, buf2 for read back (and compare)\n");

    printf("\n");

    printf("Logical Array Partition:\t CYCLIC darray partition\n\n");
    printf("\t Process Grid Cartesian:\t [");
    TEST_PRINT_LIST(array_of_psizes, 0, ndims-1, 1);
    printf("]\n\n");

    printf("Local Memory Darray Regions:  NPROCS = %d\n\n", nprocs);

  }

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  for (i=0; i<nprocs; i++) {
    if (rank == i) {
      printf("\t Proc %2d of %2d:  \t starts = [", rank, nprocs);
      TEST_PRINT_LIST(local_starts, 0, ndims-1, 1);
      printf("]\n");
      printf("\t \t\t\t counts = [");
      TEST_PRINT_LIST(local_edges, 0, ndims-1, 1);
      printf("]\n");
      printf("\t \t\t\t strides = [");
      TEST_PRINT_LIST(local_strides, 0, ndims-1, 1);
      printf("]\n\n");
 
      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
    } else {
      /* Synchronizer: processes should print out their stuffs in turn :) */
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  fflush(stdout); sleep(2);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    printf("\n");

 /* RESET the target region of buf2 */

  MPI_Pack_size(local_sz, mpi_etype, MPI_COMM_SELF, &packsz);
  tmpbuf = (TEST_NATIVE_ETYPE *)
	   malloc(local_sz*sizeof(TEST_NATIVE_ETYPE) + packsz);
  packbuf = (void *)(tmpbuf + local_sz);
  for (i=0; i<local_sz; i++)
    tmpbuf[i] = (TEST_NATIVE_ETYPE)(-1);
  packpos = 0;
  MPI_Pack((void *)tmpbuf, local_sz, mpi_etype, 
	   packbuf, packsz, &packpos, MPI_COMM_SELF);
  packpos = 0;
  MPI_Unpack(packbuf, packsz, &packpos, buf2, 1, mpi_darray, MPI_COMM_SELF);

/* Begin of TEST1: test local write-n-readback */

  if (rank == 0) {
    printf("TEST1: \n");
  }

 /* WRITE target region from buf1 */

  if (rank == 0) {
    printf("\t all procs collectively writing their darrays into Var_1 ...\n");
  }

  status = ncmpi_put_vars_all(ncid, varid_1, 
			      local_starts, local_edges, local_strides,
			      (const void *)buf1, 1, mpi_darray);
  TEST_HANDLE_ERR(status);

 /* READ target region back into buf2 */

  if (rank == 0) {
    printf("\t all procs collectively reading their darrays from Var_1 ...\n");
  }

  status = ncmpi_get_vars_all(ncid, varid_1, 
			      local_starts, local_edges, local_strides,
			      (void *)buf2, 1, mpi_darray);
  TEST_HANDLE_ERR(status);

 /* COMPARE buf1 and buf2 for equality (test1) */

  if ( memcmp((void *)buf1, (void *)buf2, total_sz*sizeof(TEST_NATIVE_ETYPE)) )
    success = 0;
  else
    success = 1;

  MPI_Reduce(&success, &successall, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    if (successall)
      printf("\t PASS: memory darray layout passes test1! \n\n");
    else
      printf("\t ERROR: memory darray layout fails test1! \n\n");
  }

/* End of TEST1 */

/* Begin of TEST2: varify global/whole file data */

  if (rank == 0) {
    printf("TEST2: \n");
  }

 /* READ whole array back into resetted buf2 */

  if (rank == 0) {
    printf("\t reading whole array from Var_1 ...\n");
  }

  for (i=0; i<total_sz; i++) {
    buf2[i] = (TEST_NATIVE_ETYPE)(-1);
  }
  status = ncmpi_get_var_all(ncid, varid_1, (void *)buf2, total_sz, mpi_etype);
  TEST_HANDLE_ERR(status);

 /* COMPARE buf1 and buf2 for equality (test2) */

  if ( memcmp((void *)buf1, (void *)buf2, total_sz*sizeof(TEST_NATIVE_ETYPE)) )
    success = 0;
  else
    success = 1;

  MPI_Reduce(&success, &successall, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    if (successall)
      printf("\t PASS: memory darray layout passes test2! \n\n");
    else
      printf("\t ERROR: memory darray layout fails test2! \n\n");
  }

/* End of TEST2 */

 /* test finalization */

  ncmpi_close(ncid);

  MPI_Type_free(&mpi_darray);
  free(tmpbuf);
  free(buf1);
  free(array_of_sizes);

  MPI_Finalize();

  return 0;
}

