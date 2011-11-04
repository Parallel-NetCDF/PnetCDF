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


int main(int argc, char **argv) {

	int ret, ncfile, nprocs, rank, dimid, varid1, varid2, ndims=1;
	MPI_Offset start, count=1;
	char buf[13] = "Hello World\n";
	int data;
	int requests[2], statuses[2];

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	ret = ncmpi_create(MPI_COMM_WORLD, argv[1],
			    NC_WRITE|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_def_dim(ncfile, "d1", nprocs, &dimid);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_def_var(ncfile, "v1", NC_INT, ndims, &dimid, &varid1);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_def_var(ncfile, "v2", NC_INT, ndims, &dimid, &varid2);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_put_att_text(ncfile, NC_GLOBAL, "string", 13, buf);
	if (ret != NC_NOERR) handle_error(ret);

	/* all processors defined the dimensions, attributes, and variables,
	 * but here in ncmpi_enddef is the one place where metadata I/O
	 * happens.  Behind the scenes, rank 0 takes the information and writes
	 * the netcdf header.  All processes communicate to ensure they have
	 * the same (cached) view of the dataset */

	ret = ncmpi_enddef(ncfile); if (ret != NC_NOERR) handle_error(ret);

	start=rank, count=1, data=rank;

	/* in this simple example every process writes its rank to two 1d variables */

	/* we used a basic MPI_INT type to this flexible mode call, but could
	 * have used any derived MPI datatype that describes application data
	 * structures */

	/* furthermore, we use the non-blocking interface to essentially
	 * schedule the two write operations.  No i/o actually happens here,
	 * which is why these routines do not need to be collective.   */
	ret = ncmpi_iput_vara(ncfile, varid1, &start, &count, &data, count,
			MPI_INT, &(requests[0]) );
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_iput_vara(ncfile, varid2, &start, &count, &data, count, 
			MPI_INT, &(requests[1]));
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_wait_all(ncfile, 2, requests, statuses);
	if (ret != NC_NOERR) handle_error(ret);



	ret = ncmpi_close(ncfile);
	if (ret != NC_NOERR) handle_error(ret);

	MPI_Finalize();

	return 0;
}
