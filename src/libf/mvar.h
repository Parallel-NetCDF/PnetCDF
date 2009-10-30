#include <stdio.h>
#include <malloc.h>
#include "pnetcdf.h"

typedef struct nc_req   
{    NCMPI_Request req;
     MPI_Aint reqid;
     struct nc_req *next;
} lnc_req;

