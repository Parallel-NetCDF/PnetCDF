/*********************************************************************
 *
 *  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* simple demonstration of pnetcdf
 * knowing nothing about the file, read in the variables.
 *
 * This example demonstrates the flexible interface */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <pnetcdf.h>

static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char **argv) {

    int i, j, rank, nprocs, ret;
    int ncfile, ndims, nvars, ngatts, unlimited;
    int var_ndims, var_natts;;
    MPI_Offset *dim_sizes, var_size;
    MPI_Offset *start, *count;
    char filename[256], varname[NC_MAX_NAME+1];
    int *dimids=NULL;
    nc_type type;
    int *data=NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (rank == 0) printf("Usage: %s filename\n", argv[0]);
        MPI_Finalize();
        exit(-1);
    }
    if (argc > 1) snprintf(filename, 256, "%s", argv[1]);
    else          strcpy(filename, "testfile.nc");

    ret = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    /* reader knows nothing about dataset, but we can interrogate with query
     * routines: ncmpi_inq tells us how many of each kind of "thing"
     * (dimension, variable, attribute) we will find in the file  */

    /* no communication needed after ncmpi_open: all processors have a cached
     * view of the metadata once ncmpi_open returns */

    ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    /* we do not really need the name of the dimension or the variable for
     * reading in this example.  we could, in a different example, take the
     * name of a variable on the command line and read just that one */

    dim_sizes = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));
    /* netcdf dimension identifiers are allocated sequentially starting
     * at zero; same for variable identifiers */
    for(i=0; i<ndims; i++)  {
        ret = ncmpi_inq_dimlen(ncfile, i, &(dim_sizes[i]) );
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
    }

    for(i=0; i<nvars; i++) {
        /* obtain the number of dimensions of variable i, so we can allocate
         * the dimids array */
        ret = ncmpi_inq_varndims(ncfile, i, &var_ndims);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        dimids = (int*) malloc(var_ndims * sizeof(int));

        /* much less coordination in this case compared to rank 0 doing all
         * the i/o: everyone already has the necessary information */
        ret = ncmpi_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                &var_natts);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        start = (MPI_Offset*) calloc(var_ndims, sizeof(MPI_Offset));
        count = (MPI_Offset*) calloc(var_ndims, sizeof(MPI_Offset));

        /* we will simply decompose along one dimension.  Generally the
         * application has some algorithm for domain decomposition.  Note
         * that data decomposition can have an impact on i/o performance.
         * Often it's best just to do what is natural for the application,
         * but something to consider if performance is not what was
         * expected/desired */

        start[0] = (dim_sizes[dimids[0]]/nprocs)*rank;
        count[0] = (dim_sizes[dimids[0]]/nprocs);
        var_size = count[0];

        for (j=1; j<var_ndims; j++) {
            start[j] = 0;
            count[j] = dim_sizes[dimids[j]];
            var_size *= count[j];
        }

        switch(type) {
            case NC_INT:
                data = (int*) calloc(var_size, sizeof(int));
                ret = ncmpi_get_vara_all(ncfile, i, start, count, data,
                        var_size, MPI_INT);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                break;
            default:
                /* we can do this for all the known netcdf types but this
                 * example is already getting too long  */
                fprintf(stderr, "unsupported NetCDF type \n");
        }

        free(start);
        free(count);
        free(dimids);
        if (data != NULL) free(data);
    }

    ret = ncmpi_close(ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);
    free(dim_sizes);

    MPI_Finalize();
    return 0;
}

/*
 *vim: ts=8 sts=4 sw=4 noexpandtab */
