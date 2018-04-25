#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <string.h>
#include <unistd.h>

#include "test_dtype.h"

/* Configurable test parameters: edit it as you wish */
#define TEST_DEFAULT_FILE "test_subarray.nc" /* file to be written/read */
#define TEST_DIMS 3 /* number of dimensions */
#define TEST_N 17 /* size for the shortest dimension of the array */
#define TEST_M 3 /* size of each dimension is increasingly multiplied by M */
#define TEST_ARRAY_ORDER MPI_ORDER_C

int ndims = TEST_DIMS;
char *filename = TEST_DEFAULT_FILE;
int test_n = TEST_N;
double test_m = TEST_M;
int order = TEST_ARRAY_ORDER;

static
void parse_args(int argc, char **argv, int rank) {
  extern char * optarg;
  int c;

  if (argc == 1)
    return;

  while ((c=getopt(argc,argv,"h:f:d:n:m:o:"))!=-1) {
    switch(c) {
      case 'f':
        filename = optarg;
        break;

      case 'd':
        ndims = (int)strtol(optarg,NULL,10);
        break;

      case 'n':
        test_n = (int)strtol(optarg,NULL,10);
        break;

      case 'm':
        test_m = atof(optarg);
        break;

      case 'o':
        if ( !strncmp(optarg, "MPI_ORDER_C", 11) || !strcmp(optarg, "0") )
          order = MPI_ORDER_C;
        else if ( !strncmp(optarg, "MPI_ORDER_FORTRAN", 17) ||
                  !strcmp(optarg, "1") )
          order = MPI_ORDER_FORTRAN;
        else {
          if (rank == 0) {
            fprintf(stderr,
                  "Invalid order_of_array specified in command line args! \n");
            fprintf(stderr,
                  "[-o <order_of_array>] 0/1: MPI_ORDER_C/MPI_ORDER_FORTRAN\n");
          }
          TEST_EXIT(-1);
        }
        break;

      case 'h':
      default:
        if (rank == 0) {
          /* print help */
          printf("Usage: %s [options]\n", argv[0]);
          printf("or     <mpirun -np 8> %s [options]\n", argv[0]);
          printf("   [-h] [-f <scratch_file>] [-d <ndims_of_array>]\n");
          printf("   [-n <first_dimsz>] [-m <multiplier_for_each_dimsz>]\n");
          printf("   [-o <order_of_array>]\n\n");
          printf("\t[-help] print this help.\n");
          printf("\t[-f] filename to be used.\n");
          printf("\t     NOTE: This test will create/overwrite the file!\n");
          printf("\t[-d] ndims of array to be tested. (int)\n");
          printf("\t[-n] size of the first array dimension. (int)\n");
          printf("\t[-m] multiplier for the size of each later dimension\n");
          printf("\t     over the size of previous dimension. (int/float)\n");
          printf("\t[-o] memory layout order of the array\n");
          printf("\t     0/1 for MPI_ORDER_C/MPI_ORDER_FORTRAN.\n\n");
	}

        TEST_EXIT( c != 'h' );
    }
  }

  return;
}

static
void partition_array(int ndims,
		     int *total_sizes,
		     int *local_subsizes,
		     int *local_starts,
		     int nprocs,
		     int myrank)
{
  int k;
  int odd;

  for (k=0; k<ndims; k++) {
    local_subsizes[k] = total_sizes[k];
    local_starts[k] = 0;
  }

  /* starting from the last dimension */
  k = ndims - 1;

  while (nprocs > 1) {

    /* cut this dimension into two halves:
	- [0 : ceiling{nprocs/2}) get first half
	- [ceiling{nprocs/2} : nprocs) get second half
    */

    odd = local_subsizes[k] % 2;
    local_subsizes[k] /= 2;
    if (myrank >= (nprocs+1)/2) {
      local_starts[k] += local_subsizes[k];
      local_subsizes[k] += odd;
      myrank -= (nprocs+1)/2; /* rank_in_myhalfgroup */
      nprocs /= 2; /* nprocs_in_myhalfgroup */
    } else {
      nprocs = (nprocs+1)/2; /* nprocs_in_myhalfgroup */
    }

    /* move to the next dimension */
    if (--k < 0)
      k = ndims-1;
  }

}

int main(int argc, char** argv) {

  int i;
  double power_M;
  int *array_of_sizes, *array_of_subsizes, *array_of_starts;
  int ncid, *dimids, varid_1, varid_2;
  MPI_Offset *local_starts, *local_edges;
  char dimname[20];
  nc_type nc_etype;
  MPI_Datatype mpi_etype, mpi_subarray;
  TEST_NATIVE_ETYPE *buf1, *buf2, *tmpbuf;
  void *packbuf;
  int packsz;
  int packpos;
  int total_sz, local_sz;
  int nprocs, rank;
  int status;
  int success, successall;
  int malloc_failure, malloc_failure_any;
  int request;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 /* test initializations: nc_file, nc_variables, dimensions, arrays */

  parse_args(argc, argv, rank);
  TEST_SET_NCMPI_ETYPE(nc_etype, mpi_etype);
#ifdef TEST_NCTYPE
  nc_etype = TEST_NCTYPE;
#endif
  if (rank == 0) {
    printf("testing memory subarray layout ...\n");
  }

  status = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER,
			MPI_INFO_NULL, &ncid);
  TEST_HANDLE_ERR(status);

  array_of_sizes = (int *)
		   malloc(sizeof(int)*ndims*4 + sizeof(MPI_Offset)*ndims*4);
  array_of_subsizes = array_of_sizes + ndims;
  array_of_starts = array_of_subsizes + ndims;
  dimids = array_of_starts + ndims;
  local_starts = (MPI_Offset *)(dimids + ndims);
  local_edges = local_starts + ndims;

  total_sz = 1;
  power_M = 1;
  for (i=0; i<ndims; i++, power_M*=test_m) {
    array_of_sizes[i] = (int)(test_n*power_M);
    if (array_of_sizes[i] < 1) {
      /* lower bound check */
      array_of_sizes[i] = 1;
    } else if ( (double)total_sz*array_of_sizes[i] > (double)TEST_MAX_INT ){
      /* upper bound check */
      if (rank == 0) {
        fprintf(stderr, "Total size of array is too big to be represented\n");
        fprintf(stderr, "Current size = %f, Max size = %d\n",
                (double)total_sz*array_of_sizes[i], TEST_MAX_INT);
      }
      TEST_EXIT(-1);
    }
    total_sz *= array_of_sizes[i];
    sprintf(dimname, "dim_%d", i);
    status = ncmpi_def_dim(ncid, dimname,
			   (MPI_Offset)array_of_sizes[i], dimids+i);
    TEST_HANDLE_ERR(status);
  }

  if (order == MPI_ORDER_FORTRAN) {
    /* reverse the filearray dimension, since NC always use C ORDER */
    TEST_REVERSE(dimids, ndims, int);
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
  malloc_failure = (buf1 == NULL ||
		    (float)total_sz*sizeof(TEST_NATIVE_ETYPE)*2 > TEST_MAX_INT);
  MPI_Allreduce(&malloc_failure, &malloc_failure_any, 1, MPI_INT,
                MPI_LOR, MPI_COMM_WORLD);
  if (malloc_failure_any) {
    if (rank == 0) {
      fprintf(stderr, "malloc(%2.3fMB) failed!\n",
              (float)total_sz*sizeof(TEST_NATIVE_ETYPE)*2/1024/1024);
      fprintf(stderr, "The whole array may be too big for malloc to handle!\n");
      fprintf(stderr, "Please choose smaller array size.\n");
    }
    TEST_EXIT(-1);
  }

  buf2 = buf1 + total_sz;
  for (i=0; i<total_sz; i++)
    /* just make sure any type can represent the number */
    /* and make it irregular (cycle != power of 2) for random test */
    buf1[i] = buf2[i] = (TEST_NATIVE_ETYPE)(i%127);

 /* PARTITION and calculate the local target region */

  partition_array(ndims,
		  array_of_sizes, array_of_subsizes, array_of_starts,
                  nprocs, rank);

  local_sz = 1;
  for (i=0; i<ndims; i++) {
    local_sz *= array_of_subsizes[i];
    local_edges[i] = (MPI_Offset)array_of_subsizes[i];
    local_starts[i] = (MPI_Offset)array_of_starts[i];
  }

  if (order == MPI_ORDER_FORTRAN) {
    /* reverse the filearray dimension, since NC always use C ORDER */
    TEST_REVERSE(local_edges, ndims, MPI_Offset);
    TEST_REVERSE(local_starts, ndims, MPI_Offset);
  }

 /* CREATE local subarray memory view */

  if (local_sz == 0)
    MPI_Type_contiguous(0, mpi_etype, &mpi_subarray);
  else
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
    if (order == MPI_ORDER_C) {
      TEST_PRINT_LIST(array_of_sizes, 0, ndims-1, 1);
    } else {
      TEST_PRINT_LIST(array_of_sizes, ndims-1, 0, -1);
    }
    printf("] Always ORDER_C\n");

    printf("\t NC Var_2 Shape:\t [");
    if (order == MPI_ORDER_C) {
      TEST_PRINT_LIST(array_of_sizes, ndims-1, 0, -1);
    } else {
      TEST_PRINT_LIST(array_of_sizes, 0, ndims-1, 1);
    }
    printf("] Always ORDER_C\n");

    printf("\t Memory Array Shape:\t [");
    TEST_PRINT_LIST(array_of_sizes, 0, ndims-1, 1);
    printf("] %s\n", ((order==MPI_ORDER_C)?"MPI_ORDER_C":"MPI_ORDER_FORTRAN"));
    printf("\t Memory Array Copys: buf1 for write, buf2 for read back (and compare)\n");

    printf("\n");

    printf("Logical Array Partition:\t BLOCK partition along all dimensions\n\n");

    printf("Access Pattern (subarray):  NPROCS = %d\n\n", nprocs);

  }

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  for (i=0; i<nprocs; i++) {
    if (rank == i) {
      printf("\t Proc %2d of %2d:  starts = [", rank, nprocs);
      TEST_PRINT_LIST(local_starts, 0, ndims-1, 1);
      printf("], ");
      printf("counts = [");
      TEST_PRINT_LIST(local_edges, 0, ndims-1, 1);
      printf("] \n");

    }
    fflush(stdout);
    /* Synchronizer: processes should print out their stuffs in turn :) */
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (rank == 0) {
    printf("\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);

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

  fflush(stdout);

  if (rank == 0) {
    printf("TEST1: \n");
  }

 /* WRITE target region from buf1 */

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf("\t [nonblocking] all procs writing their subarrays into Var_1 ...\n");
  }

  status = ncmpi_begin_indep_data(ncid);
  TEST_HANDLE_ERR(status);

  status = ncmpi_iput_vara(ncid, varid_1, local_starts, local_edges,
			   (const void *)buf1, 1, mpi_subarray, &request);
  TEST_HANDLE_ERR(status);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf("\t nonblocking I/O returns ...\n");
    fflush(stdout);
  }

  ncmpi_wait(ncid, 1, &request, &status);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf("\t nonblocking I/O finishes ...\n");
  }

  status = ncmpi_end_indep_data(ncid);
  TEST_HANDLE_ERR(status);

 /* READ target region back into buf2 */

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf("\t [nonblocking] all procs reading their subarrays from Var_1 ...\n");
  }

  status = ncmpi_begin_indep_data(ncid);
  TEST_HANDLE_ERR(status);

  status = ncmpi_iget_vara(ncid, varid_1, local_starts, local_edges,
			   (void *)buf2, 1, mpi_subarray, &request);
  TEST_HANDLE_ERR(status);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf("\t nonblocking I/O returns ...\n");
    fflush(stdout);
  }

  ncmpi_wait(ncid, 1, &request, &status);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf("\t nonblocking I/O finishes ...\n");
  }

  status = ncmpi_end_indep_data(ncid);
  TEST_HANDLE_ERR(status);

 /* COMPARE buf1 and buf2 for equality */

  if (memcmp((void *)buf1, (void *)buf2, total_sz*sizeof(TEST_NATIVE_ETYPE)))
    success = 0;
  else
    success = 1;

  MPI_Allreduce(&success, &successall, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

  if (rank == 0) {
    if (successall)
      printf("\t PASS: memory subarray layout passes test1! \n\n");
    else
      printf("\t ERROR: memory subarray layout fails test1! \n\n");
  }

/* End of TEST1 */

 /* test finalization */

  ncmpi_close(ncid);

  MPI_Type_free(&mpi_subarray);
  free(tmpbuf);
  free(buf1);
  free(array_of_sizes);

  MPI_Finalize();

  return !successall;
}

