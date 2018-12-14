/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This example tests if ncmpi_end_indep_data() works properly when nonblocking
 * APIs are called in the independent data mode, but the wait call is made later
 * in collective data mode.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 4
#define NX 10
#define NDIMS 2

int main(int argc, char** argv)
{
    char filename[256];
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, varid, dimid[2], req, st;
    MPI_Offset start[2], count[2], stride[2];
    unsigned char buffer[NY][NX];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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
        sprintf(cmd_str, "*** TESTING C   %s for ncmpi_end_indep_data ", basename(argv[0]));
        printf("%-66s ------ ",cmd_str);
        free(cmd_str);
    }

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_DATA,
                       MPI_INFO_NULL, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs,    &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_UBYTE, NDIMS, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    for (i=0; i<NY; i++) for (j=0; j<NX; j++) buffer[i][j] = rank+10;

     start[0] = 0;   start[1] = rank;
     count[0] = NY;  count[1] = NX;
    stride[0] = 1;  stride[1] = nprocs;
    err = ncmpi_buffer_attach(ncid, NY*NX); CHECK_ERR

    err = ncmpi_begin_indep_data(ncid); CHECK_ERR
    err = ncmpi_bput_vars_uchar(ncid, varid, start, count, stride,
                                &buffer[0][0], &req);
    CHECK_ERR

    /* check if write buffer contents have been altered */
    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++) {
            if (buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: put buffer[%d][%d]=%hhu altered, should be %d\n",
                       __LINE__,__FILE__,i,j,buffer[i][j],rank+10);
                nerrs++;
            }
        }

    err = ncmpi_end_indep_data(ncid); CHECK_ERR

    /* calling wait API after exiting independent data mode on purpose */
    err = ncmpi_wait_all(ncid, 1, &req, &st); CHECK_ERR
    err = st; CHECK_ERR

    /* check if write buffer contents have been altered */
    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++) {
            if (buffer[i][j] != rank+10) {
                printf("Error at line %d in %s: put buffer[%d][%d]=%hhu altered, should be %d\n",
                       __LINE__,__FILE__,i,j,buffer[i][j],rank+10);
                nerrs++;
            }
        }

    err = ncmpi_buffer_detach(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0) {
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
            ncmpi_inq_malloc_list();
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

