/*********************************************************************
 *
 *  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* simple demonstration of pnetcdf
 * text attribute on dataset
 * rank 0 reads into 1-d array, broadcasts to all.  This is a dumb way
 * to do parallel I/O, but folks do this sometimes... */

/* This program reads a file created by pnetcdf-write-from-master.c, say file
   named output.nc with the following contents, shown by running ncmpidump command .

    % mpiexec -n 4 pnetcdf-read-from-master /orangefs/wkliao/output.nc

    % ncmpidump /orangefs/wkliao/output.nc
    netcdf output {
    // file format: CDF-2 (large file)
    dimensions:
            d1 = 4 ;
    variables:
            int v1(d1) ;
            int v2(d1) ;

    // global attributes:
                :string = "Hello World\n",
        "" ;
    data:

         v1 = 0, 1, 2, 3 ;

         v2 = 0, 1, 2, 3 ;
    }
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <pnetcdf.h>

#define MPI_ERR(err) { \
    if (err != MPI_SUCCESS) { \
        char err_string[MPI_MAX_ERROR_STRING+1]; \
        int  err_len; \
        MPI_Error_string(err, err_string, &err_len); \
        fprintf(stderr, "Rank %d: Error at line %d (%s)\n",rank,__LINE__,err_string); \
    } \
}

static void handle_error(int err, int lineno)
{
    fprintf(stderr, "Error at line %d of %s: %s\n", lineno, __FILE__, ncmpi_strerror(err));
    MPI_Abort(MPI_COMM_WORLD, 1);
}


int main(int argc, char **argv) {

    int i, j=0, rank, nprocs, err;
    int ncfile, ndims, nvars, ngatts, unlimited, var_ndims, var_natts;;
    int *dimids=NULL;
    char filename[256], varname[NC_MAX_NAME+1];
    MPI_Offset *dim_sizes=NULL, var_size;
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

    if (rank == 0) {
        err = ncmpi_open(MPI_COMM_SELF, filename,
                         NC_NOWRITE, MPI_INFO_NULL, &ncfile);
        if (err != NC_NOERR) handle_error(err, __LINE__);

        /* reader knows nothing about dataset, but we can interrogate with
         * query routines: ncmpi_inq tells us how many of each kind of
         * "thing" (dimension, variable, attribute) we will find in the
         * file  */

        err = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
        if (err != NC_NOERR) handle_error(err, __LINE__);

        /* we do not really need the name of the dimension or the variable
         * for reading in this example.  we could, in a different example,
         * take the name of a variable on the command line and read just
         * that one */

        dim_sizes = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));
        /* netcdf dimension identifiers are allocated sequentially starting
         * at zero; same for variable identifiers */
        for(i=0; i<ndims; i++)  {
            err = ncmpi_inq_dimlen(ncfile, i, &(dim_sizes[i]) );
            if (err != NC_NOERR) handle_error(err, __LINE__);
        }
    }

    /* need to inform other MPI processors how many variables we will send */
    err = MPI_Bcast(&nvars, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_ERR(err)

    for(i=0; i<nvars; i++) {
        /* rank 0 will find out the size of each variable, read it, and
         * broadcast it to the rest of the processors */
        if (rank == 0) {
            /* obtain the number of dimensions of variable i, so we can
             * allocate the dimids array */
            err = ncmpi_inq_varndims(ncfile, i, &var_ndims);
            if (err != NC_NOERR) handle_error(err, __LINE__);
            dimids = (int*) malloc(var_ndims * sizeof(int));

            err = ncmpi_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                    &var_natts);
            if (err != NC_NOERR) handle_error(err, __LINE__);

            for (j=0, var_size=1; j<var_ndims; j++)  {
                var_size *= dim_sizes[dimids[j]];
            }
            free(dimids);
        }
        /* oddity: there's no predefined MPI_Offset type */
        err = MPI_Bcast(&var_size, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
        MPI_ERR(err)

        data = (int*) calloc(var_size, sizeof(int));

        if (rank == 0) {
            switch(type) {
                case NC_INT:
                    /* now we have the variable identifiers and we know how big
                     * they are */

                    /* this approach is not scalable: i/o happens from a single
                     * processor.  This approach can be ok if the amount of
                     * data is quite small, but almost always the underlying
                     * MPI-IO library can do a better job */

                    err = ncmpi_get_var_int_all(ncfile, j, data);
                    if (err != NC_NOERR) handle_error(err, __LINE__);
                    break;
                default:
                    /* we can do this for all the known netcdf types but this
                     * example is already getting too long  */
                    fprintf(stderr, "unsupported NetCDF type \n");
            }
        }

        /*and finally all processors have the data */
        err = MPI_Bcast(data, var_size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_ERR(err)

        /* Here, every process can do computation on the local buffer, data,
           or copy the contents to somewhere else */
        free(data);
    }

    if (rank == 0) {
        err = ncmpi_close(ncfile);
        if (err != NC_NOERR) handle_error(err, __LINE__);
        free(dim_sizes);
    }

    MPI_Finalize();
    return 0;
}

/*
 *vim: ts=8 sts=4 sw=4 noexpandtab */
