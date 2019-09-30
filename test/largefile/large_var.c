/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program is to test if MPI filetypes are defined correctly for accessing
 * arrays with more than 2G elements.
 *
 * This program calls ncmpi_put_vara_int_all() and ncmpi_iput_vara_int()
 * to write multiple subarray of a large 3D 4-byte integer array. It first
 * defines a netCDF variable of size 4 x 10 x 4294967296.
 * 1st write: subarray of 1 x 2 x 10 at the start 1 x 8 x (2G + rank * 10)
 * 2nd write: subarray of 1 x 2 x 5  at the start 3 x 8 x (2G + rank * 10)
 * 3rd write: subarray of 1 x 1 x 5  at the start 3 x 8 x (2G + rank * 10 + 5)
 * 4th write: subarray of 1 x 1 x 5  at the start 3 x 9 x (2G + rank * 10 + 5)
 *
 * The written contents are read back by a different process to check
 * correctness. Two reads are performed: PnetCDF and MPI-IO.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

#define FOUR_G 4294967296LL
#define TWO_G  2147483648LL
#define ONE_G  1073741824LL

#define NZ 4
#define NY 10
#define NX FOUR_G

#ifndef WORDS_BIGENDIAN
/* Endianness byte swap: done in-place */
#define SWAP(x,y) {tmp = (x); (x) = (y); (y) = tmp;}
static void
swapn(void       *buf,
           MPI_Offset  nelems)
{
    int  i;
    unsigned char tmp, *op = (unsigned char*)buf;

    while (nelems-- > 0) {
        for (i=0; i<2; i++)
            SWAP(op[i], op[3-i])
        op += 4;
    }
}
#endif

int main(int argc, char** argv)
{
    char filename[256];
    int i, j, rank, nprocs, err, nerrs=0, bufsize, expected;
    int ncid, cmode, varid, dimid[3], req[3], st[3], *buf, *buf_ptr;
    MPI_Offset offset, var_offset, start[3], count[3];
    MPI_File fh;
    MPI_Status status;
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
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
        sprintf(cmd_str, "*** TESTING C   %s for writing to a large variable ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_ds_write", "disable");
    MPI_Info_set(info, "romio_ds_read", "disable");

#ifdef ENABLE_NETCDF4
    /* Test NetCDF-4 feature */
    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_NETCDF4;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    /* define dimensions Z, Y, and X */
    err = ncmpi_def_dim(ncid, "Z", NZ, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[2]); CHECK_ERR

    /* define a big 1D variable of ubyte type, so there are 2 big variables */
    err = ncmpi_def_var(ncid, "big_1D", NC_UBYTE, 1, &dimid[2], &varid);
    CHECK_ERR

    /* define a big 3D variable of integer type */
    err = ncmpi_def_var(ncid, "var", NC_INT, 3, dimid, &varid); CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* now we are in data mode */
    start[0] = 1;
    start[1] = 8;
    start[2] = TWO_G + 10 * rank;
    count[0] = 1;
    count[1] = 2;
    count[2] = 10;

    buf = (int*) malloc(40 * sizeof(int));
    for (i=0; i<40; i++) buf[i] = rank*100 + i;

    /* write: subarray of 1 x 2 x 10 at start 1 x 8 x (2G + rank * 10) */
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* now test nonblocking put */
    /* write: subarray of 1 x 2 x 5 at start 3 x 8 x (2G + rank * 10) */
    start[0] = 3;
    count[2] = 5;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf+20);
    CHECK_ERR

    /* write: subarray of 1 x 1 x 5 at start 3 x 8 x (2G + rank * 10 + 5) */
    start[2] += 5;
    count[1]  = 1;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf+30);
    CHECK_ERR

    /* write: subarray of 1 x 1 x 5 at start 3 x 9 x (2G + rank * 10 + 5) */
    start[1] = 9;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf+35);
    CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the same file and read back for validation */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    err = ncmpi_inq_varid(ncid, "var", &varid); CHECK_ERR

    /* read back subarray written by the process (rank+1)%nprocs */
    expected = (rank == nprocs - 1) ? 0 : (rank+1)*100;

    /* read back ((rank+1)%nprocs)'s 1st write:
     *                      subarray of 1x2x10 at start 1x8x(2G+rank*10) */
    start[0] = 1;
    start[1] = 8;
    start[2] = TWO_G + ((rank == nprocs - 1) ? 0 : 10 * (rank + 1));
    count[0] = 1;
    count[1] = 2;
    count[2] = 10;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    /* read back ((rank+1)%nprocs)'s 2nd write:
     *                      subarray of 1x2x5 at start 3x8x(2G+rank*10) */
    start[0] = 3;
    count[2] = 5;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    /* read back ((rank+1)%nprocs)'s 3rd write:
     *                      subarray of 1x1x5 at start 3x8x(2G+rank*10+5) */
    start[2] += 5;
    count[1] = 1;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    /* read back ((rank+1)%nprocs)'s 4th write:
     *                      subarray of 1x1x5 at start 3x9x(2G+rank*10+5) */
    start[1] = 9;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);
#endif

    /* Test classic format */
    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    /* define dimensions Z, Y, and X */
    err = ncmpi_def_dim(ncid, "Z", NZ, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[2]); CHECK_ERR

    /* define a big 1D variable of ubyte type, so there are 2 big variables */
    err = ncmpi_def_var(ncid, "big_1D", NC_UBYTE, 1, &dimid[2], &varid);
    CHECK_ERR

    /* define a big 3D variable of integer type */
    err = ncmpi_def_var(ncid, "var", NC_INT, 3, dimid, &varid); CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* get the beginning of file offset for the variable */
    err = ncmpi_inq_varoffset(ncid, varid, &var_offset); CHECK_ERR

    /* now we are in data mode */
    start[0] = 1;
    start[1] = 8;
    start[2] = TWO_G + 10 * rank;
    count[0] = 1;
    count[1] = 2;
    count[2] = 10;

    buf = (int*) malloc(40 * sizeof(int));
    for (i=0; i<40; i++) buf[i] = rank*100 + i;

    /* write: subarray of 1 x 2 x 10 at start 1 x 8 x (2G + rank * 10) */
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* now test nonblocking put */
    /* write: subarray of 1 x 2 x 5 at start 3 x 8 x (2G + rank * 10) */
    start[0] = 3;
    count[2] = 5;
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf+20, &req[0]);
    CHECK_ERR

    /* write: subarray of 1 x 1 x 5 at start 3 x 8 x (2G + rank * 10 + 5) */
    start[2] += 5;
    count[1]  = 1;
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf+30, &req[1]);
    CHECK_ERR

    /* write: subarray of 1 x 1 x 5 at start 3 x 9 x (2G + rank * 10 + 5) */
    start[1] = 9;
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf+35, &req[2]);
    CHECK_ERR

    err = ncmpi_wait_all(ncid, 3, req, st); CHECK_ERR
    for (i=0; i<3; i++) {
        err = st[i];
        CHECK_ERR
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the same file and read back for validation */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    err = ncmpi_inq_varid(ncid, "var", &varid); CHECK_ERR

    /* read back subarray written by the process (rank+1)%nprocs */
    expected = (rank == nprocs - 1) ? 0 : (rank+1)*100;

    /* read back ((rank+1)%nprocs)'s 1st write:
     *                      subarray of 1x2x10 at start 1x8x(2G+rank*10) */
    start[0] = 1;
    start[1] = 8;
    start[2] = TWO_G + ((rank == nprocs - 1) ? 0 : 10 * (rank + 1));
    count[0] = 1;
    count[1] = 2;
    count[2] = 10;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    /* read back ((rank+1)%nprocs)'s 2nd write:
     *                      subarray of 1x2x5 at start 3x8x(2G+rank*10) */
    start[0] = 3;
    count[2] = 5;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    /* read back ((rank+1)%nprocs)'s 3rd write:
     *                      subarray of 1x1x5 at start 3x8x(2G+rank*10+5) */
    start[2] += 5;
    count[1] = 1;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    /* read back ((rank+1)%nprocs)'s 4th write:
     *                      subarray of 1x1x5 at start 3x9x(2G+rank*10+5) */
    start[1] = 9;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* MPI file open the same file and read back for validation */
    err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (err != MPI_SUCCESS) {
        int errorStringLen;
        char errorString[MPI_MAX_ERROR_STRING];
        MPI_Error_string(err, errorString, &errorStringLen);
        printf("MPI error MPI_File_open : %s\n", errorString);
    }

    /* read back subarray written by the process (rank+1)%nprocs */
    start[0] = 1;
    start[1] = 8;
    start[2] = TWO_G + ((rank == nprocs - 1) ? 0 : 10 * (rank + 1));
    count[0] = 1;
    count[1] = 2;
    count[2] = 10;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    buf_ptr = buf;
    for (i=0; i<count[0]; i++) {
        for (j=0; j<count[1]; j++) {
            offset = var_offset + ((start[0] + i) * NY * NX + (start[1] + j) * NX + start[2]) * sizeof(int);
            MPI_File_read_at(fh, offset, buf_ptr, count[2], MPI_INT, &status);
#ifndef WORDS_BIGENDIAN
            swapn(buf_ptr, count[2]);
#endif
            buf_ptr += count[2];
        }
    }

    expected = (rank == nprocs - 1) ? 0 : (rank+1)*100;
    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    start[0] = 3;
    count[2] = 5;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    buf_ptr = buf;
    for (i=0; i<count[0]; i++) {
        for (j=0; j<count[1]; j++) {
            offset = var_offset + ((start[0] + i) * NY * NX + (start[1] + j) * NX + start[2]) * sizeof(int);
            MPI_File_read_at(fh, offset, buf_ptr, count[2], MPI_INT, &status);
#ifndef WORDS_BIGENDIAN
            swapn(buf_ptr, count[2]);
#endif
            buf_ptr += count[2];
        }
    }

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    start[2] += 5;
    count[1]  = 1;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    buf_ptr = buf;
    for (i=0; i<count[0]; i++) {
        for (j=0; j<count[1]; j++) {
            offset = var_offset + ((start[0] + i) * NY * NX + (start[1] + j) * NX + start[2]) * sizeof(int);
            MPI_File_read_at(fh, offset, buf_ptr, count[2], MPI_INT, &status);
#ifndef WORDS_BIGENDIAN
            swapn(buf_ptr, count[2]);
#endif
            buf_ptr += count[2];
        }
    }

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    start[1] = 9;

    /* initialize the contents of the array to a different value */
    bufsize = count[0]*count[1]*count[2];
    for (i=0; i<bufsize; i++) buf[i] = -1;

    buf_ptr = buf;
    for (i=0; i<count[0]; i++) {
        for (j=0; j<count[1]; j++) {
            offset = var_offset + ((start[0] + i) * NY * NX + (start[1] + j) * NX + start[2]) * sizeof(int);
            MPI_File_read_at(fh, offset, buf_ptr, count[2], MPI_INT, &status);
#ifndef WORDS_BIGENDIAN
            swapn(buf_ptr, count[2]);
#endif
            buf_ptr += count[2];
        }
    }

    /* check if the contents of buf are expected */
    for (i=0; i<bufsize; i++) {
        if (buf[i] != expected) {
            printf("Error file %s line %d: expect buf[%d]=%d, but got %d\n",
                   __FILE__, __LINE__, i, expected, buf[i]);
            nerrs++;
        }
        expected++;
    }

    MPI_File_close(&fh);

    free(buf);
    MPI_Info_free(&info);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

