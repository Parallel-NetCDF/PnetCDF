#include <mpi.h>
#include <pnetcdf.h>
#include <stdio.h>
#include <stdlib.h>

static void handle_error(int status)
{
	fprintf(stderr, "%s\n", ncmpi_strerror(status));
	exit(-1);
}
#define NDIMS 1
int main(int argc, char **argv)
{
	int ret, ncfile, dimid, varid, ndims=NDIMS;
	int i, nprocs, rank;
	MPI_Offset start[NDIMS] = {0};
	MPI_Offset count[NDIMS] = {0};
	int buf[512];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	ret = ncmpi_create(MPI_COMM_WORLD, argv[1], 0, MPI_INFO_NULL, &ncfile);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_def_dim(ncfile, "d1", nprocs, &dimid);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_def_var(ncfile, "v1", NC_INT, ndims, &dimid, &varid);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_enddef(ncfile);
	if (ret != NC_NOERR) handle_error(ret);

	start[0] = rank;
	count[0] = 1;
	for (i=0; i<512; i++) buf[i] = rank;
	ret = ncmpi_put_vars_int_all(ncfile, varid, start, count, NULL, buf);
	if (ret != NC_NOERR) handle_error(ret);

	ret = ncmpi_close(ncfile);
	if (ret != NC_NOERR) handle_error(ret);

	MPI_Finalize();
	return 0;
}
