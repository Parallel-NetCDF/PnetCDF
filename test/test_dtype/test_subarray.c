#include <mpi.h>
#include <stdio.h>
#include <pnetcdf.h>
#include <string.h>

#include "test_dtype.h"

#define TEST_DEFAULT_FILE "test_subarray.nc"
#define TEST_DIMS 3 /* number of dimensions */
#define TEST_N 16 /* size for the shortest dimension of the memory array */
#define TEST_ARRAY_ORDER MPI_ORDER_C

void partition_array(int ndims, 
		     int *total_sizes, 
		     int *local_subsizes, 
		     int *local_starts, 
		     int nprocs, 
		     int myrank)
{
  int k;

  for (k=0; k<ndims; k++) {
    local_subsizes[k] = total_sizes[k];
    local_starts[k] = 0;
  }

  /* starting from the last dimension */
  k = ndims - 1; 

  while (nprocs > 1) {

    /* cut this dimension into two halves:
	- [0 : nprocs/2) get first half
	- [nprocs/2 : nprocs) get second half 
    */

    local_subsizes[k] /= 2; 
    nprocs /= 2;
    if (myrank >= nprocs) 
      local_starts[k] += local_subsizes[k];
    myrank %= nprocs;

    /* move to the next dimension */
    if (--k < 0)
      k = ndims-1;
  }

}

int main(int argc, char** argv) {

  int i;
  int power2;
  int ndims;
  int *array_of_sizes, *array_of_subsizes, *array_of_starts;
  int ncid, *dimids, varid_1, varid_2;
  MPI_Offset *local_starts, *local_edges, *stride, *imap;
  char dimname[20], *filename;
  int order;
  ncmpi_type nc_etype;
  MPI_Datatype mpi_etype, mpi_subarray;
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
    printf("testing memory subarray layout ...\n");
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
		   malloc(sizeof(int)*ndims*4 + sizeof(MPI_Offset)*ndims*4);
  array_of_subsizes = array_of_sizes + ndims;
  array_of_starts = array_of_subsizes + ndims;
  dimids = array_of_starts + ndims;
  local_starts = (MPI_Offset *)(dimids + ndims);
  local_edges = local_starts + ndims;
  stride = local_edges + ndims;
  imap = stride + ndims;

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

  TEST_REVERSE(dimids, ndims, int);
  status = ncmpi_def_var(ncid, "var_2", nc_etype, ndims, dimids, &varid_2);
  TEST_HANDLE_ERR(status);

  status = ncmpi_enddef(ncid);
  TEST_HANDLE_ERR(status);

  if (rank == 0) {
    printf("\t Filesize = %2.3fMB, MAX_Memory_needed = %2.3fMB\n\n", 
	   2*total_sz*TEST_NCTYPE_LEN(nc_etype)/1024.0/1024.0, 
	   ( (2*total_sz + 4*total_sz/nprocs)*sizeof(TEST_NATIVE_ETYPE)
	   + total_sz*TEST_NCTYPE_LEN(nc_etype) )/1024.0/1024.0);
  }

  buf1 = (TEST_NATIVE_ETYPE *)malloc(total_sz*sizeof(TEST_NATIVE_ETYPE)*2);
  buf2 = buf1 + total_sz;
  for (i=0; i<total_sz; i++)
    /* just make sure any type can represent the number */
    /* and make it irregular (cycle != power of 2) for random test */
    buf1[i] = buf2[i] = (TEST_NATIVE_ETYPE)(i%127); 

 /* PARTITION and calculate the local target region */

  partition_array(ndims, 
		  array_of_sizes, array_of_subsizes, array_of_starts,
                  nprocs, rank);

  local_sz = array_of_subsizes[0];
  local_edges[0] = (MPI_Offset)array_of_subsizes[0];
  local_starts[0] = (MPI_Offset)array_of_starts[0];
  for (i=1; i<ndims; i++) {
    local_sz *= array_of_subsizes[i];
    local_edges[i] = (MPI_Offset)array_of_subsizes[i];
    local_starts[i] = (MPI_Offset)array_of_starts[i];
  }

 /* CREATE local subarray memory view */

  MPI_Type_create_subarray(ndims,
                           array_of_sizes,
                           array_of_subsizes,
                           array_of_starts,
                           order,
                           mpi_etype,
                           &mpi_subarray);
  MPI_Type_commit(&mpi_subarray);

 /* PRINT stats */
  if (rank == 0) {
    printf("Initialization:  NDIMS = %d, NATIVE_ETYPE = %s, NC_TYPE = %s\n\n",
	   ndims, TEST_NATIVE_ETYPE_STR, TEST_GET_NCTYPE_STR(nc_etype));

    printf("\t NC Var_1 Shape:\t [");
    TEST_PRINT_LIST(array_of_sizes, 0, ndims-1, 1);
    printf("]\n");

    printf("\t NC Var_2 Shape:\t [");
    TEST_PRINT_LIST(array_of_sizes, ndims-1, 0, -1);
    printf("]\n");

    printf("\t Memory Array Shape:\t [");
    TEST_PRINT_LIST(array_of_sizes, 0, ndims-1, 1);
    printf("]\n");
    printf("\t Memory Array Copys: buf1 for write, buf2 for read back (and compare)\n");

    printf("\n");
    printf("Local Memory Subarray Regions:  NPROCS = %d\n\n", nprocs);

  }

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  for (i=0; i<nprocs; i++) {
    if (rank == i) {
      printf("\t Proc %2d of %2d:  starts = [", rank, nprocs);
      TEST_PRINT_LIST(array_of_starts, 0, ndims-1, 1);
      printf("], ");
      printf("counts = [");
      TEST_PRINT_LIST(array_of_subsizes, 0, ndims-1, 1);
      printf("] \n");
 
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
  MPI_Unpack(packbuf, packsz, &packpos, buf2, 1, mpi_subarray, MPI_COMM_SELF);

/* Begin of TEST1: test local write-n-readback */

  if (rank == 0) {
    printf("TEST1: \n");
  }

 /* WRITE target region from buf1 */

  if (rank == 0) {
    printf("\t all procs collectively writing their subarrays into Var_1 ...\n");
  }

  status = ncmpi_put_vara_all(ncid, varid_1, local_starts, local_edges, 
			      (const void *)buf1, 1, mpi_subarray);
  TEST_HANDLE_ERR(status);

 /* READ target region back into buf2 */

  if (rank == 0) {
    printf("\t all procs collectively reading their subarrays from Var_1 ...\n");
  }

  status = ncmpi_get_vara_all(ncid, varid_1, local_starts, local_edges,
			      (void *)buf2, 1, mpi_subarray);
  TEST_HANDLE_ERR(status);

 /* COMPARE buf1 and buf2 for equality (test1) */

  if ( memcmp((void *)buf1, (void *)buf2, total_sz*sizeof(TEST_NATIVE_ETYPE)) )
    success = 0;
  else
    success = 1;

  MPI_Reduce(&success, &successall, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    if (successall)
      printf("\t PASS: memory subarray layout passes test1! \n\n");
    else
      printf("\t ERROR: memory subarray layout fails test1! \n\n");
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
      printf("\t PASS: memory subarray layout passes test2! \n\n");
    else
      printf("\t ERROR: memory subarray layout fails test2! \n\n");
  }

/* End of TEST2 */

/* Begin of TEST3: read/write local subarray from/into transposed file array */

  if (rank == 0) {
    printf("TEST3: \n");
    printf("\t testing the transposed subarray access through put/get_varm ...\n");
  }

 /* TRANSPOSE the subarray fileview */

  TEST_REVERSE(local_starts, ndims, MPI_Offset);
  TEST_REVERSE(local_edges, ndims, MPI_Offset);

 /* SET the transpose imap (and stride) */

  stride[0] = 1;
  imap[0] = 1;
  for (i=1; i<ndims; i++) {
    stride[i] = 1;
    imap[i] = imap[i-1]*(int)local_edges[i-1]; /* a transpose imap */
  }

 /* WRITE local subarray (through transpose imap) into file */

  if (rank == 0) {
    printf("\t all procs collectively writing the transposed array into Var_2 ...\n");
  }

  status = ncmpi_put_varm_all(ncid, varid_2, 
			      local_starts, local_edges, stride, imap,
			      buf1, 1, mpi_subarray);
  TEST_HANDLE_ERR(status);

 /* RESET the target region of buf2 */

  packpos = 0;
  MPI_Unpack(packbuf, packsz, &packpos, buf2, 1, mpi_subarray, MPI_COMM_SELF);

 /* READ back the local subarray (through transpose imap) from file */

  if (rank == 0) {
    printf("\t each proc reading subarray from transposed array of Var_2 ...\n");
  }

  status = ncmpi_get_varm_all(ncid, varid_2, 
			      local_starts, local_edges, stride, imap,
			      buf2, 1, mpi_subarray);
  TEST_HANDLE_ERR(status);

 /* COMPARE buf1 and buf2 for equality (test3) */

  if ( memcmp((void *)buf1, (void *)buf2, total_sz*sizeof(TEST_NATIVE_ETYPE)) )
    success = 0;
  else
    success = 1;

  MPI_Reduce(&success, &successall, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    if (successall)
      printf("\t PASS: memory subarray layout + varm passes test3! \n\n");
    else
      printf("\t ERROR: memory subarray layout + varm fails test3! \n\n");
  }

/* End of TEST3 */

 /* test finalization */

  ncmpi_close(ncid);

  MPI_Type_free(&mpi_subarray);
  free(tmpbuf);
  free(buf1);
  free(array_of_sizes);

  MPI_Finalize();

  return 0;
}

