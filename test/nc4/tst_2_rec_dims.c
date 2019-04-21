/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests reading files with more than 2 record dimensions with
 * NetCDF4 driver
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 8
#define NX 5

int main(int argc, char** argv) {
    char filename[512];
    int i, rank, nprocs, err, nerrs=0, buf[16];
    int ncid, mode, dimids[2], varid, num_rec_vars;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Offset recsize;
    size_t start[2], count[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for reading file with 2 rec dims ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* create a NetCDF-4 file. Note NC_MPIIO is used in NetCDF 4.6.1 and
     * earlier, but ignored in 4.6.2 and after. */
    mode = NC_CLOBBER | NC_NETCDF4 | NC_MPIIO;
    err = nc_create_par(filename, mode, comm, info, &ncid); CHECK_ERR

    err = nc_def_dim(ncid, "y", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = nc_def_dim(ncid, "x", NC_UNLIMITED, &dimids[1]); CHECK_ERR
    err = nc_def_var(ncid, "data", NC_INT, 2, dimids, &varid); CHECK_ERR

    err = nc_var_par_access(ncid, varid, NC_COLLECTIVE); CHECK_ERR

    for (i=0; i<16; i++) buf[i] = rank + i;

    /* No need to explicitly end define mode for netCDF-4 files */
    start[0] = 0;
    start[1] = 0;
    count[0] = 2;
    count[1] = 1;
    err = nc_put_vara_int(ncid, varid, start, count, buf); CHECK_ERR

    /* Close the file. */
    err = nc_close(ncid); CHECK_ERR

    /* reopen the file using PnetCDF */
    err = ncmpi_open(comm, filename, NC_NOWRITE, info, &ncid); CHECK_ERR

    err = ncmpi_inq_num_rec_vars(ncid, &num_rec_vars); CHECK_ERR
    if (num_rec_vars != 1){
        printf("Error at line %d of %s: expect num_rec_vars %d but got %d\n",
                __LINE__,__FILE__,2,num_rec_vars);
        nerrs++;
    }
    
    err = ncmpi_inq_num_fix_vars(ncid, &num_rec_vars); CHECK_ERR
    if (num_rec_vars != 0){
        printf("Error at line %d of %s: expect num_fix_vars %d but got %d\n",
                __LINE__,__FILE__,0,num_rec_vars);
        nerrs++;
    }       

    err = ncmpi_inq_recsize(ncid, &recsize);  CHECK_ERR
    if (recsize != 4){
        printf("Error at line %d of %s: expect recsize %lld but got %lld\n",
                __LINE__,__FILE__, (MPI_Offset)4, recsize);
        nerrs++;
    }  

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, comm);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

