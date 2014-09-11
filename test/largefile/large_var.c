/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests using ncmpi_put_vara_int_all() and ncmpi_iput_vara_int()
 * to write multiple subarray of a large 3D 4-byte integer array. It first
 * defines a netCDF variable of size 4 x 10 x 4294967296.
 * 1st write: subarray of 1 x 2 x 10 at the start 1 x 8 x (2G + rank * 10)
 * 2nd write: subarray of 1 x 2 x 5  at the start 3 x 8 x (2G + rank * 10)
 * 3rd write: subarray of 1 x 1 x 5  at the start 3 x 8 x (2G + rank * 10 + 5)
 * 4th write: subarray of 1 x 1 x 5  at the start 3 x 9 x (2G + rank * 10 + 5)
 *
 * The written contents are read back by a different process to check
 * correctness.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <mpi.h>
#include <pnetcdf.h>

#define FOUR_G 4294967296
#define TWO_G  2147483648

#define NZ 4
#define NY 10
#define NX FOUR_G

#define ERR if (err!=NC_NOERR) {printf("Error at line %d: %s\n", __LINE__,ncmpi_strerror(err)); exit(-1);}

int main(int argc, char** argv)
{
    char filename[128];
    int i, rank, nprocs, err, pass, bufsize;
    int ncid, cmode, varid, dimid[3], req[3], st[3], *buf;
    MPI_Offset start[3], count[3];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    strcpy(filename, "testfile.nc");
    if (argc == 2) strcpy(filename, argv[1]);
    MPI_Bcast(filename, 128, MPI_CHAR, 0, MPI_COMM_WORLD);

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* define dimensions Z, Y, and X */
    err = ncmpi_def_dim(ncid, "Z", NZ, &dimid[0]);
    ERR
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[1]);
    ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[2]);
    ERR

    /* define a 3D variable of integer type */
    err = ncmpi_def_var(ncid, "var", NC_INT, 3, dimid, &varid);
    ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid);
    ERR

    /* now we are in data mode */
    start[0] = 1;
    start[1] = NY - 2;
    start[2] = TWO_G + 10 * rank;
    count[0] = 1;
    count[1] = 2;
    count[2] = 10;

    bufsize = count[0]*count[1]*count[2];
    buf = (int*) malloc(bufsize * sizeof(int));
    for (i=0; i<bufsize; i++) buf[i] = rank*100 + i;

    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf);
    ERR

    /* now test nonblocking put */
    /* rearrange buffer contents */
    for (i=5;  i<10; i++) buf[i] = rank*100 + i + 5;
    for (i=10; i<15; i++) buf[i] = rank*100 + i - 5;

    start[0] = NZ-1;
    count[2] = 5;
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf, &req[0]);
    ERR

    start[2] += 5;
    count[1]  = 1;
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf+10, &req[1]);
    ERR

    start[1]++;
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf+15, &req[2]);
    ERR

    err = ncmpi_wait_all(ncid, 3, req, st);
    ERR

    err = ncmpi_close(ncid);
    ERR

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); ERR

    err = ncmpi_inq_varid(ncid, "var", &varid); ERR

    /* initialize the contents of the array to a different value */
    for (i=0; i<bufsize; i++) buf[i] = -1;

    /* read back subarray written by the process (rank+1)%nprocs */
    start[0] = 1;
    start[1] = NY - 2;
    start[2] = TWO_G + ((rank == nprocs - 1) ? 0 : 10 * (rank + 1));
    count[0] = 1;
    count[1] = 2;
    count[2] = 10;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); ERR

    /* check if the contents of buf are expected */
    pass = 1;
    for (i=0; i<bufsize; i++) {
        int expected = (rank == nprocs - 1) ? i : (rank+1)*100 + i;
        if (buf[i] != expected) {
            printf("%d: Unexpected read buf[%d]=%d, should be %d\n",
                   rank, i, buf[i], expected);
            pass = 0;
        }
    }

    for (i=0; i<bufsize; i++) buf[i] = -1;

    start[0] = NZ-1;
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        int expected = (rank == nprocs - 1) ? i : (rank+1)*100 + i;
        if (buf[i] != expected) {
            printf("%d: 2 Unexpected read buf[%d]=%d, should be %d\n",
                   rank, i, buf[i], expected);
            pass = 0;
        }
    }

    err = ncmpi_close(ncid); ERR

    MPI_Allreduce(MPI_IN_PLACE, &pass, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    free(buf);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    char cmd_str[80];
    sprintf(cmd_str, "*** TESTING C   %s for writing to a large variable ", argv[0]);
    if (rank == 0) {
        if (pass) printf("%-66s ------ pass\n", cmd_str);
        else      printf("%-66s ------ failed\n", cmd_str);
    }

    MPI_Finalize();
    return 0;
}

