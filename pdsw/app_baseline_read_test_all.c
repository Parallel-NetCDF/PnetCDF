/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>


static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}
const char *source_name = NULL;



/*----< read_metadata_test() >-------------------------------------------------------*/
static int
read_metadata_test(MPI_Comm comm, const char *filename, int cmode)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, varid, blkid, dimid[2];
    char str_att[128];
    float float_att[100];

    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    
    double total_read_time = 0;
    double read_start_time = MPI_Wtime();
    /* open the newly created file for read only -----------------------------*/
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    ERR

    char var_name[20];    
    int nvars;
    int ndims, v_ndims;
    int* v_dimids;
    MPI_Offset* dim_sizes;
    MPI_Offset var_size = 1;



    err = ncmpi_inq(ncid, NULL, &ndims, &nvars, NULL);
    ERR;
    for (varid = 0; varid < nvars; varid++){
        err = ncmpi_inq_var(ncid, varid, var_name, NULL, &v_ndims, NULL, NULL);
        v_dimids = (int*) malloc(v_ndims * sizeof(int));
        dim_sizes = (MPI_Offset*) malloc(v_ndims * sizeof(MPI_Offset));
        err = ncmpi_inq_var(ncid, varid, NULL, NULL, NULL, v_dimids, NULL);
        for (int j = 0; j < v_ndims; j++){
            err = ncmpi_inq_dimlen(ncid, v_dimids[j], &dim_sizes[j]);
        }
        free(dim_sizes);
        free(v_dimids);
    }



    err = ncmpi_close(ncid);
    ERR
    total_read_time += MPI_Wtime() - read_start_time;
    double max_time, min_time;
    MPI_Reduce(&total_read_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_read_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Max read time: %f seconds\n",  max_time);
        printf("Min read time: %f seconds\n", min_time);
    }

    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    int i, rank, kind=0, cmode=0, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;

    if (argc < 2) {
            if (rank == 0)
                fprintf(stderr, "Usage: %s <source_file> \n", argv[0]);
            MPI_Finalize();
            return 1;
        }
        source_name = argv[1];

#ifndef PNETCDF_DRIVER_NETCDF4
    /* netcdf4 driver is not enabled, skip */
    if (kind == 3 || kind == 4) {
        MPI_Finalize();
        return 0;
    }
#endif
    nerrs += read_metadata_test(MPI_COMM_WORLD, source_name, cmode);

    MPI_Finalize();
    return (nerrs > 0);
}