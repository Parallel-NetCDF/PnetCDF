/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
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
#include <pnetcdf.h>

#include <testutils.h>

#define NY 8
#define NX 5

#define FILENAME "rec2.nc"

int main(int argc, char** argv) {
    int rank, nprocs, err, nerrs=0;
    int ncid, num_rec_vars;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Offset recsize;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for def_var_fill ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* reopen the file and read data back */
    err = ncmpi_open(comm, FILENAME, NC_NOWRITE | NC_NETCDF4, info, &ncid); CHECK_ERR

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

