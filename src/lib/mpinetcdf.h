#include <mpi.h>
#include "nc.h"

/* Begin Prototypes */
int ncmpi_create(MPI_Comm comm, const char *path, int cmode, MPI_Info info, int *ncidp);

int ncmpi_open(MPI_Comm comm, const char *path, int omode, MPI_Info info, int *ncidp);

int ncmpi_enddef(int ncid);

int ncmpi_close(int ncid);

/* Begin _dim */


/* End _dim */
/* Begin {put,get}_att */


/* End {put,get}_att */
/* Begin _var */


/* End _var */
/* Begin {put,get}_vara */

int ncmpi_put_vara_int_all(int ncid, int varid,
                    const size_t [], const size_t [],
                    const int *op);

int ncmpi_put_vara_int(int ncid, int varid,
                const size_t [], const size_t [],
                const int *op);

int ncmpi_put_vara_float_all(int ncid, int varid,
                    const size_t [], const size_t [],
                    const float *op);

int ncmpi_put_vara_float(int ncid, int varid,
                const size_t [], const size_t [],
                const float *op);

int ncmpi_put_vara_double_all(int ncid, int varid,
                    const size_t [], const size_t [],
                    const double *op);

int ncmpi_put_vara_double(int ncid, int varid,
                const size_t [], const size_t [],
                const double *op);

int ncmpi_get_vara_int_all(int ncid, int varid,
                    const size_t [], const size_t [],
                    int *ip);

int ncmpi_get_vara_int(int ncid, int varid,
                const size_t [], const size_t [],
                int *ip);

int ncmpi_get_vara_float_all(int ncid, int varid,
                    const size_t [], const size_t [],
                    float *ip);

int ncmpi_get_vara_float(int ncid, int varid,
                const size_t [], const size_t [],
                float *ip);

int ncmpi_get_vara_double_all(int ncid, int varid,
                    const size_t [], const size_t [],
                    double *ip);

int ncmpi_get_vara_double(int ncid, int varid,
                const size_t [], const size_t [],
                double *ip);

/* End Prototypes */

/* End {put,get}_vara */
