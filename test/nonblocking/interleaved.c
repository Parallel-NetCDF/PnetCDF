/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests nonblocking APIs for handling interleaved file types.
 * In the first case,
 *    It first defines a netCDF variable of size 10 x 18.
 *    3 calls to ncmpi_iput_vara_int(), where the first one interleaves with
 *    the second and third writes and its access region covers both of the
 *    two writes. The 2nd and 3rd writes do not interleave with each other.
 *
 * In the second case,
 *    It first defines a netCDF variable of size 10 x 18.
 *    4 calls to ncmpi_iput_vara_int(),
 *    First  write: a subarray of size 8 x 2 at the start offsets of 0 x 3
 *    Second write: a subarray of size 3 x 5 at the start offsets of 6 x 7
 *    Third  write: a subarray of size 2 x 5 at the start offsets of 6 x 12
 *    Fourth write: a subarray of size 1 x 5 at the start offsets of 8 x 12
 *    Note the first write interleaves with the second and third writes and
 *    the fourth write does not interleave with any other writes.
 *
 * % mpiexec -n 1 ./interleaved
 * % ncmpidump testfile.nc
 * netcdf testfile {
 * // file format: CDF-5 (big variables)
 * dimensions:
 * 	Y = 10 ;
 * 	X = 18 ;
 * variables:
 * 	int var0(Y, X) ;
 * 	int var1(Y, X) ;
 *
 * data:
 *
 *  var0 =
 *    -1, -1, -1, 50, 51, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 52, 53, -1, -1, -1, 60, 61, 62, 63, 64, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 54, 55, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 56, 57, -1, -1, 70, 71, 72, 73, 74, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 58, 59, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 60, 61, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 62, 63, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 64, 65, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ;
 * }
 *  var1 =
 *    -1, -1, -1, 50, 51, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 52, 53, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 54, 55, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 56, 57, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 58, 59, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 60, 61, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
 *    -1, -1, -1, 62, 63, -1, -1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, -1,
 *    -1, -1, -1, 64, 65, -1, -1, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, -1,
 *    -1, -1, -1, -1, -1, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, -1,
 *    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ;
 * }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 10
#define NX 18

#define CHECK_CONTENTS(buf, exp) { \
    if (buf[i] != (exp)) { \
        printf("Error at line %d in %s: put buffer[%d] altered to %d, expect %d\n", \
        __LINE__,__FILE__,i,buf[i],(exp)); \
        nerrs++; \
    } \
}
int main(int argc, char** argv)
{
    char filename[256];
    int i, j, rank, nprocs, err, nerrs=0, expected;
    int ncid, cmode, varid[2], dimid[2], req[4], st[4], *buf;
    int *buf0, *buf1, *buf2;
    size_t len;
    MPI_Offset start[2], count[2];
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* this program is intended to run on one process */
    if (rank) goto fn_exit;

    /* get command-line arguments */
    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for writing interleaved fileviews ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_cb_write", "disable");
    MPI_Info_set(info, "ind_wr_buffer_size", "8");
    /* these 2 hints are required to cause a core dump if r1758 fix is not
     * presented */

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_SELF, filename, cmode, info, &ncid); CHECK_ERR

    MPI_Info_free(&info);

    /* define dimensions Y and X */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]); CHECK_ERR

    /* define 2D variables of integer type */
    err = ncmpi_def_var(ncid, "var0", NC_INT, 2, dimid, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 2, dimid, &varid[1]); CHECK_ERR

    /* enable fill mode */
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* now we are in data mode */
    buf = (int*) malloc(NY*NX * sizeof(int));

    /* fill the entire variable var0 with -1s */
    for (i=0; i<NY*NX; i++) buf[i] = -1;
    err = ncmpi_put_var_int_all(ncid, varid[0], buf); CHECK_ERR

    /* When using burst buffering, flush the log to prevent new value being
     * skipped due to overlaping domain
     */
    err = ncmpi_flush(ncid); CHECK_ERR

    /* write 8 x 2 elements so this only interleaves the next two
     * iput requests */
    start[0] = 0; start[1] = 3;
    count[0] = 8; count[1] = 2;
    len = (size_t)(count[0] * count[1]);
    buf0 = (int*) malloc(len * sizeof(int));
    for (i=0; i<len; i++) buf0[i] = 50+i;
    err = ncmpi_iput_vara_int(ncid, varid[0], start, count, buf0, &req[0]);
    CHECK_ERR

    /* write 1 x 3 elements */
    start[0] = 1; start[1] = 8;
    count[0] = 1; count[1] = 5;
    len = (size_t)(count[0] * count[1]);
    buf1 = (int*) malloc(len * sizeof(int));
    for (i=0; i<len; i++) buf1[i] = 60+i;
    err = ncmpi_iput_vara_int(ncid, varid[0], start, count, buf1, &req[1]);
    CHECK_ERR

    /* write 1 x 3 elements */
    start[0] = 3; start[1] = 7;
    count[0] = 1; count[1] = 5;
    len = (size_t)(count[0] * count[1]);
    buf2 = (int*) malloc(len * sizeof(int));
    for (i=0; i<len; i++) buf2[i] = 70+i;
    err = ncmpi_iput_vara_int(ncid, varid[0], start, count, buf2, &req[2]);
    CHECK_ERR

    err = ncmpi_wait_all(ncid, 3, req, st); CHECK_ERR
    free(buf0); free(buf1); free(buf2);

    /* fill the entire variable var1 with -1s */
    for (i=0; i<NY*NX; i++) buf[i] = -1;
    err = ncmpi_put_var_int_all(ncid, varid[1], buf); CHECK_ERR

    /* When using burst buffering, flush the log to prevent new value being
     * skipped due to overlaping domain
     */
    err = ncmpi_flush(ncid); CHECK_ERR

    /* write 8 x 2 elements so this only interleaves the next two iput
     * requests */
    start[0] = 0; start[1] = 3;
    count[0] = 8; count[1] = 2;
    len = (size_t)(count[0] * count[1]);
    buf0 = (int*) malloc(len * sizeof(int));
    for (i=0; i<count[0]*count[1]; i++) buf0[i] = 50+i;
    err = ncmpi_iput_vara_int(ncid, varid[1], start, count, buf0, &req[0]);
    CHECK_ERR

    /* rearrange buffer contents, as buf is 2D */
    for (i=0;  i<5;  i++) buf[i] = 10 + i;
    for (i=5;  i<10; i++) buf[i] = 10 + i +  5;
    for (i=10; i<15; i++) buf[i] = 10 + i + 10;
    start[0] = 6; start[1] = 7;
    count[0] = 3; count[1] = 5;
    err = ncmpi_iput_vara_int(ncid, varid[1], start, count, buf, &req[1]);
    CHECK_ERR

    for (i=15; i<20; i++) buf[i] = 10 + i - 10;
    for (i=20; i<25; i++) buf[i] = 10 + i -  5;
    start[0] = 6; start[1] = 12;
    count[0] = 2; count[1] = 5;
    err = ncmpi_iput_vara_int(ncid, varid[1], start, count, buf+15, &req[2]);
    CHECK_ERR

    for (i=25; i<30; i++) buf[i] = 10 + i;
    start[0] = 8; start[1] = 12;
    count[0] = 1; count[1] = 5;
    err = ncmpi_iput_vara_int(ncid, varid[1], start, count, buf+25, &req[3]);
    CHECK_ERR

    err = ncmpi_wait_all(ncid, 4, req, st); CHECK_ERR

    /* check if write buffer contents have been altered */
    for (i=0;  i<16; i++) CHECK_CONTENTS(buf0, 50 + i)
    for (i=0;  i<5;  i++) CHECK_CONTENTS(buf, 10 + i)
    for (i=5;  i<10; i++) CHECK_CONTENTS(buf, 10 + i +  5)
    for (i=10; i<15; i++) CHECK_CONTENTS(buf, 10 + i + 10)
    for (i=15; i<20; i++) CHECK_CONTENTS(buf, 10 + i - 10)
    for (i=20; i<25; i++) CHECK_CONTENTS(buf, 10 + i -  5)
    for (i=25; i<30; i++) CHECK_CONTENTS(buf, 10 + i)

    err = ncmpi_close(ncid); CHECK_ERR
    free(buf0);

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_SELF, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "var0", &varid[0]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "var1", &varid[1]); CHECK_ERR

    /* read the entire array */
    for (i=0; i<NY*NX; i++) buf[i] = -1;
    err = ncmpi_get_var_int_all(ncid, varid[0], buf); CHECK_ERR

    /* check if the contents of buf are expected */
    expected = 50;
    for (j=0; j<8; j++) {
        for (i=3; i<5; i++) {
            if (buf[j*NX+i] != expected) {
                printf("%d: Unexpected read buf[%d][%d]=%d, should be %d\n",
                       rank, j, i, buf[j*NX+i], expected);
                nerrs++;
            }
            expected++;
        }
    }
    expected = 60;
    j = 1;
    for (i=8; i<13; i++) {
        if (buf[j*NX+i] != expected) {
            printf("%d: Unexpected read buf[%d][%d]=%d, should be %d\n",
                   rank, j, i, buf[j*NX+i], expected);
            nerrs++;
        }
        expected++;
    }
    expected = 70;
    j = 3;
    for (i=7; i<12; i++) {
        if (buf[j*NX+i] != expected) {
            printf("%d: Unexpected read buf[%d][%d]=%d, should be %d\n",
                   rank, j, i, buf[j*NX+i], expected);
            nerrs++;
        }
        expected++;
    }

    /* initialize the contents of the array to a different value */
    for (i=0; i<NY*NX; i++) buf[i] = -1;

    /* read the entire array */
    err = ncmpi_get_var_int_all(ncid, varid[1], buf); CHECK_ERR

    /* check if the contents of buf are expected */
    expected = 10;
    for (j=6; j<9; j++) {
        for (i=7; i<17; i++) {
            if (buf[j*NX+i] != expected) {
                printf("%d: Unexpected read buf[%d]=%d, should be %d\n",
                       rank, i, buf[j*NX+i], expected);
                nerrs++;
            }
            expected++;
        }
    }
    expected = 50;
    for (j=0; j<8; j++) {
        for (i=3; i<5; i++) {
            if (buf[j*NX+i] != expected) {
                printf("%d: Unexpected read buf[%d][%d]=%d, should be %d\n",
                       rank, j, i, buf[j*NX+i], expected);
                nerrs++;
            }
            expected++;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    free(buf);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR && malloc_size > 0) {
        printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n", malloc_size);
        ncmpi_inq_malloc_list();
    }

fn_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

