#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#ifndef MPI_OFFSET
#define MPI_OFFSET MPI_LONG_LONG_INT
#endif

#define NDIMS 3

#define ERR { if (err != NC_NOERR) { printf("Error at line = %d: %s\n", __LINE__, ncmpi_strerror(err)); nerrs++; } }

int main(int argc, char** argv)
{
  int i, j, rank, nprocs, err, nerrs = 0;
  int ncid, varid, num_reqs;
  double *buffer;
  float *fbuffer;
  MPI_Offset r_len, **starts = NULL, **counts = NULL;
  MPI_Offset st[3], ct[3];
  char       filename[256];
  int dimids[3];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
  }
  if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
  else           strcpy(filename, "testfile.nc");
  MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

  if (nprocs != 4) {
    if (rank == 0)
      printf("Error: this program is intended to run on 4 processes\n");
    MPI_Finalize();
    return 1;
  }



/*    original test case read from a file like this:
 *    % ncdump -h lnfm.nc
 *        netcdf lnfm {
 *        dimensions:
 *          time = UNLIMITED ; // (3 currently)
 *          lat = 94 ;
 *           lon = 192 ;
 *         variables:
 *          float lnfm(time, lat, lon) ;
 *        }
 */
  err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_WRITE|NC_FORMAT_CDF5, MPI_INFO_NULL, &ncid); CHECK_ERR
  err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimids[0]); CHECK_ERR
  err = ncmpi_def_dim(ncid, "lat", 94, &dimids[1]); CHECK_ERR
  err = ncmpi_def_dim(ncid, "lon", 192, &dimids[2]); CHECK_ERR
  err = ncmpi_def_var(ncid, "lnfm", NC_FLOAT, 3, dimids, &varid); CHECK_ERR
  err = ncmpi_enddef(ncid); CHECK_ERR

  st[0] = rank*2;
  st[1] = 0;
  st[2] = 0;

  ct[0] = 2;
  ct[1] = 1;
  ct[2] = 1;
  float scramble[1024];

  err = ncmpi_put_vara_float_all(ncid, varid, st, ct, scramble); CHECK_ERR
  err = ncmpi_close(ncid); CHECK_ERR

  /* now we can finally exercise the read path of this record varable */

  err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
  ERR

  /* pick 2 requests for 4 processes */
  /* num_reqs = 1; => works fine*/
  num_reqs = 2;

  starts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
  counts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
  starts[0] = (MPI_Offset*)  calloc(num_reqs*NDIMS, sizeof(MPI_Offset));
  counts[0] = (MPI_Offset*)  calloc(num_reqs*NDIMS, sizeof(MPI_Offset));
  for (i = 1; i < num_reqs; i++) {
    starts[i] = starts[i - 1] + NDIMS;
    counts[i] = counts[i - 1] + NDIMS;
  }

  /* assign specific starts and counts */
  if (num_reqs > 0){
    if (rank == 1) {
      starts[0][0] = 0; starts[0][1] = 0; starts[0][2] = 1;
      counts[0][0] = 1; counts[0][1] = 93; counts[0][2] = 1;
    }
 
  }

  if(num_reqs > 1){
    if (rank == 1) {
      starts[1][0] = 1; starts[1][1] = 0; starts[1][2] = 1;
      counts[1][0] = 1; counts[1][1] = 93; counts[1][2] = 2;
    }
  
  }

  r_len = 0; /* total read length for this process */
  for (i = 0; i < num_reqs; i++) {
    MPI_Offset r_req_len = 1;
    for (j = 0; j < NDIMS; j++)
      r_req_len *= counts[i][j];
    r_len += r_req_len;
  }

  /* allocate I/O buffer */
  buffer = (double*) calloc(r_len, sizeof(double));
  fbuffer = (float*) calloc(r_len, sizeof(float));

  /* set the buffer pointers to different offsets to the I/O buffer */
  varid = 0; /* only one variable in lnfm.nc */
  err = ncmpi_get_varn_double_all(ncid, varid, num_reqs, starts, counts, buffer);
  /* err = ncmpi_get_varn_float_all(ncid, varid, num_reqs, starts, counts, fbuffer); */
  ERR

  err = ncmpi_close(ncid);
  ERR

  if (rank == 3) {
    printf("Dumping some double type data read by rank 3 (count = 10) ...\n");
    for (i = 0; i < (r_len > 10)?10:r_len; i++)
      printf("%lf, ", buffer[i]);
    printf("\n");
  }

  free(buffer);
  free(fbuffer);

  free(starts[0]);
  free(counts[0]);
  free(starts);
  free(counts);

   MPI_Finalize();

  return nerrs;
}
