/*********************************************************************
 *
 *  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* simple demonstration of pnetcdf
 * text attribute on dataset
 * write out rank into 1-d array after sending to rank 0.  This is a dumb way
 * to do parallel I/O, but folks do this sometimes... */

/* This program creates a file, say named output.nc, with the following
   contents, shown by running ncmpidump command .

    % mpiexec -n 4 pnetcdf-write-from-master /orangefs/wkliao/output.nc

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

static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d of %s: %s\n", lineno, __FILE__, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char **argv) {
    int ret, ncid=0, nprocs, rank, dimid, varid1=0, varid2=0, ndims=1;
    char filename[256], buf[13] = "Hello World\n";
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
        ret = ncmpi_create(MPI_COMM_SELF, filename,
                           NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncid);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_def_dim(ncid, "d1", nprocs, &dimid);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_def_var(ncid, "v1", NC_INT, ndims, &dimid, &varid1);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_def_var(ncid, "v2", NC_INT, ndims, &dimid, &varid2);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_put_att_text(ncid, NC_GLOBAL, "string", 13, buf);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_enddef(ncid);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        /* first reason this approach is not scalable:  need to allocate
        * enough memory to hold data from all processors */
        data = (int*)calloc(nprocs, sizeof(int));
    }

    /* second reason this approach is not scalable: sending to rank 0
     * introduces a serialization point, even if using an optimized
     * collective routine */
    MPI_Gather(&rank, 1, MPI_INT, data, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        /* and lastly, the third reason this approach is not scalable: I/O
         * happens from a single processor.  This approach can be ok if the
         * amount of data is quite small, but almost always the underlying
         * MPI-IO library can do a better job */
        MPI_Offset start[1], count[1];
        start[0]=0, count[0]=nprocs;
        ret = ncmpi_put_vara_int_all(ncid, varid1, start, count, data);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_put_vara_int_all(ncid, varid2, start, count, data);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_close(ncid);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        free(data);
    }

    MPI_Finalize();
    return 0;
}
