#include <mpi.h>
#include <pnetcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* This test program check if a collective API can be nicely aborted without
 * causing program hang. It uses a case that one process deliberately cause
 * an erro, while the other does not.
 */

#define CHECK_ERROR(fn) { \
   if (rank == 0 && err != NC_NOERR) \
       printf("PE %d: %s error is %s\n",rank,fn,ncmpi_strerror(err)); \
   if (rank == 1 && err != NC_EINVALCOORDS) \
       printf("PE %d: %s error code should be %d, but got %d",rank,fn,NC_EINVALCOORDS,err); \
}

int main(int argc, char *argv[])
{
   int rank, nproc, ncid, err, dim_len, varid, dimids[1];
   int req, status;
   MPI_Offset start[1], count[1], memCountScalar;
   double buf[2];

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);

   if (nproc != 2) {
       printf("Program requires 2 processors\n");
       MPI_Finalize();
       return 1;
   }

   /* Create a 2 element vector of doubles */
   err = ncmpi_create(MPI_COMM_WORLD, "testfile.nc", NC_CLOBBER, MPI_INFO_NULL, &ncid);
   assert(err == NC_NOERR);

   err = ncmpi_def_dim(ncid, "dim", 2, &dimids[0]);
   assert(err == NC_NOERR);

   err = ncmpi_def_var(ncid, "var", NC_DOUBLE, 1, dimids, &varid);
   assert(err == NC_NOERR);

   err = ncmpi_enddef(ncid);
   assert(err == NC_NOERR);

   if (rank == 0) {
       start[0] = 0;
       count[0] = 2;
   } else if (rank == 1) {
       start[0] = 2; /* illegal for a start > defined shape */
       count[0] = 0;
   }

   err = ncmpi_put_vara_all(ncid, varid, start, count,
			    buf, count[0], MPI_DOUBLE);
   CHECK_ERROR("ncmpi_put_vara_all")

   err = ncmpi_put_vara_double_all(ncid, varid, start, count, buf);
   CHECK_ERROR("ncmpi_put_vara_double_all")

   err = ncmpi_iput_vara_double(ncid, varid, start, count, buf, &req);
   CHECK_ERROR("ncmpi_iput_vara_double")

   err = ncmpi_wait_all(ncid, 1, &req, &status);
   if (err != NC_NOERR)
       printf("PE %d: ncmpi_wait_all error is %s\n",rank,ncmpi_strerror(err));

   err = ncmpi_get_vara_all(ncid, varid, start, count,
			    buf, count[0], MPI_DOUBLE);
   CHECK_ERROR("ncmpi_get_vara_all")

   err = ncmpi_get_vara_double_all(ncid, varid, start, count, buf);
   CHECK_ERROR("ncmpi_get_vara_double_all")

   err = ncmpi_iget_vara_double(ncid, varid, start, count, buf, &req);
   CHECK_ERROR("ncmpi_iget_vara_double")

   err = ncmpi_wait_all(ncid, 1, &req, &status);
   if (err != NC_NOERR)
       printf("PE %d: ncmpi_wait_all error is %s\n",rank,ncmpi_strerror(err));

   err = ncmpi_close(ncid);
   assert(err == NC_NOERR);

   MPI_Finalize();
   return 0;
}
