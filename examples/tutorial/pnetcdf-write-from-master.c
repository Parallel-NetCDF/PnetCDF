/* simple demonstration of pnetcdf 
 * text attribute on dataset
 * write out rank into 1-d array after sending to rank 0.  This is a dumb way
 * to do parallel I/O but folks do this sometimes... */

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
	int *data;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (rank == 0) {
	    ret = ncmpi_create(MPI_COMM_SELF, argv[1],
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
	
	    ret = ncmpi_enddef(ncfile);
	    if (ret != NC_NOERR) handle_error(ret);

	    /* first reason this approach is not scalable:  need to allocate
	     * enough memory to hold data from all processors */
	    data = calloc(nprocs, sizeof(int));
	}

	/* second reason this approch is not scalable: sending to rank 0
	 * introduces a serialization point, even if using an optimized
	 * collective routine */
	MPI_Gather(&rank, 1, MPI_INT, 
			data, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
	    /* and lastly, the third reason this approach is not scalable: I/O
	     * happens from a single processor.  This approach can be ok if the
	     * amount of data is quite small, but almost always the underlying
	     * MPI-IO library can do a better job */
	    start=0, count=nprocs;
	    ret = ncmpi_put_vara_int_all(ncfile, varid1, &start, &count, data);
	    if (ret != NC_NOERR) handle_error(ret);

	    ret = ncmpi_put_vara_int_all(ncfile, varid2, &start, &count, data);
	    if (ret != NC_NOERR) handle_error(ret);

	    ret = ncmpi_close(ncfile);
	    if (ret != NC_NOERR) handle_error(ret);
	}

	MPI_Finalize();

	return 0;
}
