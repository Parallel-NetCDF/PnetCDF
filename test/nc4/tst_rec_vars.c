/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests writing record variables with NetCDF4 driver
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

int main(int argc, char** argv) {
    char filename[256];
    int rank, nprocs, err, nerrs=0;
    int ncid, cmode, varid[2], dimid[2], buf[1];
    int v1;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Offset start[1], count[1];

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
        sprintf(cmd_str, "*** TESTING C   %s for record variables to NetCDF4 file ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* create a new file for writing ------------------------------------*/
    cmode = NC_CLOBBER | NC_NETCDF4;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR

    /* define dimension */
    err = ncmpi_def_dim(ncid, "X", NC_UNLIMITED, &dimid[0]); CHECK_ERR

    /* define 2 variables of size NY x NX */
    err = ncmpi_def_var(ncid, "U", NC_INT, 1, dimid, varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "V", NC_INT, 1, dimid, varid + 1); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* write a subarray to each variable */
    start[0] = rank;
    count[0] = 1;
    buf[0] = rank + 1;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR

    start[0] = nprocs - rank - 1;
    err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, dimid[0], start); CHECK_ERR
    if (start[0] != (MPI_Offset)nprocs){
        printf("Error at line %d of %s: expect NC_UNLIMITED dimension X of len %lld but got %lld\n",
                __LINE__,__FILE__,(MPI_Offset)nprocs,start[0]);
        nerrs++;
    }

    err = ncmpi_inq_num_rec_vars(ncid, &v1); CHECK_ERR
    if (v1 != 2){
        printf("Error at line %d of %s: expect num_rec_vars %d but got %d\n",
                __LINE__,__FILE__,2,v1);
        nerrs++;
    }
    
    err = ncmpi_inq_num_fix_vars(ncid, &v1); CHECK_ERR
    if (v1 != 0){
        printf("Error at line %d of %s: expect num_fix_vars %d but got %d\n",
                __LINE__,__FILE__,0,v1);
        nerrs++;
    }       

    err = ncmpi_inq_recsize(ncid, start);  CHECK_ERR
    if (start[0] != 8){
        printf("Error at line %d of %s: expect recsize %lld but got %lld\n",
                __LINE__,__FILE__, (MPI_Offset)8, start[0]);
        nerrs++;
    }  

    err = ncmpi_close(ncid); CHECK_ERR

    /* reopen the file and read data back */
    err = ncmpi_open(comm, filename, NC_WRITE, info, &ncid); CHECK_ERR

    /* Read a subarray from each variable */
    start[0] = rank;
    count[0] = 1;
    buf[0] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[0], start, count, buf); CHECK_ERR
    if (buf[0] != rank + 1){
        printf("Error at line %d of %s: expect buf = %d but got %d\n",
                __LINE__,__FILE__, rank + 1, buf[0]);
        nerrs++;
    }  

    start[0] = nprocs - rank - 1;
    err = ncmpi_get_vara_int_all(ncid, varid[1], start, count, buf); CHECK_ERR
    if (buf[0] != rank + 1){
        printf("Error at line %d of %s: expect buf = %d but got %d\n",
                __LINE__,__FILE__, rank + 1, buf[0]);
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

