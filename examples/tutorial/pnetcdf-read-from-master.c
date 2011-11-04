/* simple demonstration of pnetcdf 
 * text attribute on dataset
 * rank 0 reads into 1-d array, broadcasts to all.  This is a dumb way
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

    int rank, nprocs;
    int ret, ncfile, ndims, nvars, ngatts, unlimited;
    int var_ndims, var_natts;;
    MPI_Offset *dim_sizes, var_size;

    char varname[NC_MAX_NAME+1];
    int dimids[NC_MAX_VAR_DIMS];
    nc_type type;

    int i, j;

    int *data;


    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (rank == 0) {
	ret = ncmpi_open(MPI_COMM_SELF, argv[1],
		NC_NOWRITE, MPI_INFO_NULL, &ncfile);
	if (ret != NC_NOERR) handle_error(ret);

	/* reader knows nothing about dataset, but we can interrogate with
	 * query routines: ncmpi_inq tells us how many of each kind of
	 * "thing" (dimension, variable, attribute) we will find in the
	 * file  */

	ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
	if (ret != NC_NOERR) handle_error(ret);

	/* we do not really need the name of the dimension or the variable
	 * for reading in this example.  we could, in a different example,
	 * take the name of a variable on the command line and read just
	 * that one */

	dim_sizes = calloc(ndims, sizeof(MPI_Offset));
	/* netcdf dimension identifiers are allocated sequentially starting
	 * at zero; same for variable identifiers */
	for(i=0; i<ndims; i++)  {
	    ret = ncmpi_inq_dimlen(ncfile, i, &(dim_sizes[i]) );
	    if (ret != NC_NOERR) handle_error(ret);
	}
    }

    /* need to inform other MPI processors how many variables we will send */
    MPI_Bcast(&nvars, 1, MPI_INT, 0, MPI_COMM_WORLD);

    for(i=0; i<nvars; i++) { 
	/* rank 0 will find out the size of each variable, read it, and
	 * broadcast it to the rest of the processors */
	if (rank == 0) {
	    ret = ncmpi_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
		    &var_natts);
	    if (ret != NC_NOERR) handle_error(ret);

	    for (j=0, var_size=1; j<var_ndims; j++)  {
		var_size *= dim_sizes[dimids[j]];
	    }
	}
	/* oddity: there's no predefined MPI_Offset type */
	MPI_Bcast(&var_size, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

	data = calloc(var_size, sizeof(int));

	if (rank == 0) {
	    switch(type) {
		case NC_INT:
		    /* now we have the variable identifiers and we know how big
		     * they are */

		    /* this approach is not scalable: i/o happens from a single
		     * processor.  This approach can be ok if the amount of
		     * data is quite small, but almost always the underlying
		     * MPI-IO library can do a better job */

		    ret = ncmpi_get_var_int_all(ncfile, j, data);
		    if (ret != NC_NOERR) handle_error(ret);
		    break;
		default:
		    /* we can do this for all the known netcdf types but this
		     * example is already getting too long  */
		    fprintf(stderr, "unsupported NetCDF type \n");
	    }
	}

	/*and finally all processors have the data */
	MPI_Bcast(data, var_size, MPI_INT, 0, MPI_COMM_WORLD);

    }

    if (rank == 0) {
	    ret = ncmpi_close(ncfile);
	    if (ret != NC_NOERR) handle_error(ret);
    }

    MPI_Finalize();

    return 0;
}

/*
 *vim: ts=8 sts=4 sw=4 noexpandtab */
