#include <mpi.h>
#include <pnetcdf.h>
#include <stdlib.h>

#include <string.h>
#include <stdio.h>

#define NC_CHECK(fn) {int ncstat; ncstat = (fn); if (ncstat != NC_NOERR) handle_error_nc(ncstat, NULL); }

static void handle_error_nc(int ncerr, char *str)
{
	        fprintf(stderr, "%s: %s\n", str, ncmpi_strerror(ncerr));
		        MPI_Abort(MPI_COMM_WORLD, 1);
}

#define VECCOUNT 4
#define BLOCKLEN  3
#define STRIDE   5
int main(int argc, char ** argv)
{
	int ncid, dimid, varid;
	MPI_Init(&argc, &argv);
	MPI_Datatype vtype, rtype, usertype;
	MPI_Aint lb, extent;
	int userbufsz, *userbuf, *cmpbuf, i, errs=0;
	int count = 25;
	double pi = 3.14159;
	MPI_Offset start, acount;

	ncmpi_create(MPI_COMM_WORLD, "vectors.nc", NC_CLOBBER, MPI_INFO_NULL,
			&ncid);
	ncmpi_def_dim(ncid, "50k", 1024*50, &dimid);
	ncmpi_def_var(ncid, "vector", NC_DOUBLE, 1, &dimid, &varid);

	ncmpi_enddef(ncid);


	MPI_Type_vector(VECCOUNT, BLOCKLEN, STRIDE, MPI_INT, &vtype);
	MPI_Type_create_resized(vtype, 0, STRIDE*VECCOUNT*sizeof(int), &rtype);
	MPI_Type_contiguous(count, rtype, &usertype);
	MPI_Type_commit(&usertype);

	MPI_Type_free(&vtype);
	MPI_Type_free(&rtype);

	MPI_Type_get_extent(usertype, &lb, &extent);
	userbufsz = extent;
	userbuf = malloc(userbufsz);
	cmpbuf = calloc(userbufsz, 1);
	for (i=0; i< userbufsz/sizeof(int); i++) {
		userbuf[i] = pi*i;
	}


	start = 10; acount = count*12;
	ncmpi_begin_indep_data(ncid);
	ncmpi_put_vara(ncid, varid, &start, &acount, 
			userbuf, 1, usertype);

	ncmpi_close(ncid);

	NC_CHECK(ncmpi_open(MPI_COMM_WORLD, "vectors.nc", NC_NOWRITE,
				MPI_INFO_NULL, &ncid));
	ncmpi_begin_indep_data(ncid);
	NC_CHECK(ncmpi_inq_varid(ncid, "vector", &varid));
	NC_CHECK(ncmpi_get_vara(ncid, varid, &start, &acount,
			cmpbuf, 1, usertype));
	ncmpi_close(ncid);

	for (i=0; errs < 10 &&  i < acount; i++) {
		/* vector of 4,3,5, so skip 4th and 5th items of every block */
		if (i%STRIDE >= BLOCKLEN) continue;
		if (userbuf[i] != cmpbuf[i]) {
			errs++;
			fprintf(stderr, "%d: expected 0x%x got 0x%x\n", 
					i, userbuf[i], cmpbuf[i]);
		}
	}
	free(userbuf);
	free(cmpbuf);
	MPI_Type_free(&usertype);
	MPI_Finalize();
	return 0;
}
