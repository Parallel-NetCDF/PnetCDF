/* simple demonstration of pnetcdf 
 * text attribute on dataset
 * write out rank into 1-d array collectively.  The most basic way to do
 * parallel i/o with pnetcdf */

#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>
#include <stdio.h>

static void handle_error(int status)
{
	fprintf(stderr, "%s\n", ncmpi_strerror(status));
	exit(-1);
}

#define DSET_NAME_LEN 1024

int main(int argc, char **argv) {

	int ret, ncfile, nprocs, rank, dimid, varid1, varid2, ndims=1;
	// MPI_Offset count=1;
	char buf[13] = "Hello World\n";
	int data;
	char filename[DSET_NAME_LEN];

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	/* Many applications find "one file per process" easy, but there are
	 * several deficincies with that approach:
	 * - here we need to construct a unique file name for each processor */
	ret = snprintf(filename, DSET_NAME_LEN, "%s.%d-%d.nc", argv[1], rank, nprocs);
	if (ret >= DSET_NAME_LEN) {
		fprintf(stderr, "name too long \n");
		exit(-1);
	}

	/* note that the communicator is still needed but now it is
	 * MPI_COMM_SELF: since each processor opens its own file we cannot use
	 * MPI_COMM_WORLD */
	ret = ncmpi_create(MPI_COMM_SELF, filename,
			    NC_WRITE|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
	if (ret != NC_NOERR) handle_error(ret);

	/* each processor writes its data to a file, so instead of an "nprocs"
	 * sized array, we just have an array big enough to hold one
	 * processor's data */
	ret = ncmpi_def_dim(ncfile, "d1", 1, &dimid);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_def_var(ncfile, "v1", NC_INT, ndims, &dimid, &varid1);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_def_var(ncfile, "v2", NC_INT, ndims, &dimid, &varid2);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_put_att_text(ncfile, NC_GLOBAL, "string", 13, buf);
	if (ret != NC_NOERR) handle_error(ret);

	/* ncmpi_enddef writes the header out as in other examples, but because
	 * each processor opened the file independently, there can be no "write
	 * and broadcast" optimization.  Instead, every processor does header
	 * i/o.  */	
	ret = ncmpi_enddef(ncfile); if (ret != NC_NOERR) handle_error(ret);

	/* the one advantage to this approach: data decomposistion is easy the
	 * application does not need to worry aobut the shape and location of
	 * the data (the 'start' and 'count' parmeters in the 'vara' family of
	 * functions) and can instead just write the enitre (small) variable */

	data=rank;

	/* in this simple example every process writes its rank to two 1d
	 * variables */
	/* When each processor writes to its own file, a whole host of
	 * optimizations cannot take place.   */
	
	ret = ncmpi_put_var_int_all(ncfile, varid1, &data);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_put_var_int_all(ncfile, varid2, &data);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_close(ncfile);
	if (ret != NC_NOERR) handle_error(ret);

	MPI_Finalize();

	return 0;
}
