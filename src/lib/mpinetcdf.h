#include <mpi.h>
#include "nc.h"

int
ncmpi_create(MPI_Comm comm, const char *path, int cmode, MPI_Info info, int *ncidp);

int
ncmpi_open(MPI_Comm comm, const char *path, int omode, MPI_Info info, int *ncidp);


int
ncmpi_enddef(int ncid);



int
ncmpi_close(int ncid);

/* Begin _dim */


/* End _dim */
/* Begin {put,get}_att */


/* End {put,get}_att */
/* Begin _var */


/* End _var */
/* Begin {put,get}_vara */

int
ncmpi_put_vara_int_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const int *op);

int
ncmpi_put_vara_int(int ncid, int varid,
                const size_t start[], const size_t count[],
                const int *op);

int
ncmpi_put_vara_float_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const float *op);

int
ncmpi_put_vara_float(int ncid, int varid,
                const size_t start[], const size_t count[],
                const float *op);

int
ncmpi_put_vara_double_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const double *op);

int
ncmpi_put_vara_double(int ncid, int varid,
                const size_t start[], const size_t count[],
                const double *op);

int
ncmpi_get_vara_int_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    int *ip);

int
ncmpi_get_vara_int(int ncid, int varid,
                const size_t start[], const size_t count[],
                int *ip);

int
ncmpi_get_vara_float_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    float *ip);

int
ncmpi_get_vara_float(int ncid, int varid,
                const size_t start[], const size_t count[],
                float *ip);

int
ncmpi_get_vara_double_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    double *ip);

int
ncmpi_get_vara_double(int ncid, int varid,
                const size_t start[], const size_t count[],
                double *ip);


/* End {put,get}_vara */
