/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

static int
check_read_contents(float *rh)
{
    int i, j;
    float k;

    /* check the contents of read */
    k = 0.0;
    for (i=0; i<6; i++) {
        for (j=0; j<4; j++) {
            if (rh[j*6+i] != k) {
#ifdef PRINT_ERR_ON_SCREEN
                printf("Error at %s:%d : expect rh[%d][%d]=%f but got %f\n",
                __FILE__,__LINE__,j,i,k,rh[j*6+i]);
#endif
                return 1;
            }
            k += 1.0;
        }
    }
#ifdef PRINT_ON_SCREEN
    /* print the contents of read */
    for (j=0; j<4; j++) {
        printf("[%2d]: ",j);
        for (i=0; i<6; i++) {
            printf("%5.1f",rh[j][i]);
        }
        printf("\n");
    }
#endif
    /* the stdout should be:
           [ 0]:   0.0  4.0  8.0 12.0 16.0 20.0
           [ 1]:   1.0  5.0  9.0 13.0 17.0 21.0
           [ 2]:   2.0  6.0 10.0 14.0 18.0 22.0
           [ 3]:   3.0  7.0 11.0 15.0 19.0 23.0
     */
    return 0;
}

static int
check_write_contents(signed char *varT)
{
    int i, j;
    /* the output from command "ncmpidump -v var test.nc" should be:
           var =
            50, 56, 62, 68,
            51, 57, 63, 69,
            52, 58, 64, 70,
            53, 59, 65, 71,
            54, 60, 66, 72,
            55, 61, 67, 73 ;
     */

    /* check if the contents of write buffer have been altered */
    for (j=0; j<4; j++) {
        for (i=0; i<6; i++) {
            if (varT[j*6+i] != j*6+i + 50) {
#ifdef PRINT_ERR_ON_SCREEN
                /* this error is a pnetcdf internal error, if occurs */
                printf("Error at line %d in %s: expecting varT[%d][%d]=%d but got %d\n",
                __LINE__,__FILE__,j,i,j*6+i + 50,varT[j*6+i]);
#endif
                return 1;
            }
        }
    }
    return 0;
}

static int
tst_fmt(char *filename, int cmode)
{
    int i, j, err, nerrs=0, rank, nprocs;
    int ncid, dimid[2], varid, req, status;
    MPI_Offset start[2], count[2], stride[2], imap[2];
    int          var[6][4];
    float         rh[4][6];
    signed char varT[4][6];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* define a variable of a 6 x 4 integer array in the nc file */
    err = ncmpi_def_dim(ncid, "Y", 6, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 4, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* create a 6 x 4 integer variable in the file with contents:
           0,  1,  2,  3,
           4,  5,  6,  7,
           8,  9, 10, 11,
          12, 13, 14, 15,
          16, 17, 18, 19,
          20, 21, 22, 23
     */
    for (j=0; j<6; j++) for (i=0; i<4; i++) var[j][i] = j*4+i;

    start[0] = 0; start[1] = 0;
    count[0] = 6; count[1] = 4;
    if (rank > 0) count[0] = count[1] = 0;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, &var[0][0]); CHECK_ERR

    if (nprocs > 1) MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "var", &varid); CHECK_ERR

    /* read the variable back in the matrix transposed way, rh is 4 x 6 */
     count[0] = 6;  count[1] = 4;
    stride[0] = 1; stride[1] = 1;
      imap[0] = 1;   imap[1] = 6;   /* would be {4, 1} if not transposing */

    if (cmode & NC_NETCDF4) {
        for (i=0; i<6; i++) for (j=0; j<4; j++) rh[j][i] = -1.0;
        err = ncmpi_get_varm_float_all(ncid, varid, start, count, stride, imap,
                                       &rh[0][0]); CHECK_ERR
        nerrs += check_read_contents(&rh[0][0]);

        /* test when stride == NULL and imap != NULL */
        for (i=0; i<6; i++) for (j=0; j<4; j++) rh[j][i] = -1.0;
        err = ncmpi_get_varm_float_all(ncid, varid, start, count, NULL, imap,
                                       &rh[0][0]); CHECK_ERR
        nerrs += check_read_contents(&rh[0][0]);
    }
    else {
        /* test nonblocking API */
        for (i=0; i<6; i++) for (j=0; j<4; j++) rh[j][i] = -1.0;
        err = ncmpi_iget_varm_float(ncid, varid, start, count, stride, imap,
                                    &rh[0][0], &req); CHECK_ERR
        err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERR
        err = status; CHECK_ERR
        nerrs += check_read_contents(&rh[0][0]);

        /* test when stride == NULL and imap != NULL */
        for (i=0; i<6; i++) for (j=0; j<4; j++) rh[j][i] = -1.0;
        err = ncmpi_iget_varm_float(ncid, varid, start, count, NULL, imap,
                                    &rh[0][0], &req); CHECK_ERR
        err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERR
        err = status; CHECK_ERR
        nerrs += check_read_contents(&rh[0][0]);

        /* test blocking API */
        for (i=0; i<6; i++) for (j=0; j<4; j++) rh[j][i] = -1.0;
        err = ncmpi_get_varm_float_all(ncid, varid, start, count, stride, imap,
                                       &rh[0][0]); CHECK_ERR
        nerrs += check_read_contents(&rh[0][0]);

        /* test when stride == NULL and imap != NULL */
        for (i=0; i<6; i++) for (j=0; j<4; j++) rh[j][i] = -1.0;
        err = ncmpi_get_varm_float_all(ncid, varid, start, count, NULL, imap,
                                       &rh[0][0]); CHECK_ERR
        nerrs += check_read_contents(&rh[0][0]);
    }

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_WRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "var", &varid); CHECK_ERR

    /* testing get_varm(), first zero-out the variable in the file */
    memset(&var[0][0], 0, 6*4*sizeof(int));
    start[0] = 0; start[1] = 0;
    count[0] = 6; count[1] = 4;
    if (rank > 0) count[0] = count[1] = 0;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, &var[0][0]); CHECK_ERR

    /* set the contents of the write buffer varT, a 4 x 6 char array
          50, 51, 52, 53, 54, 55,
          56, 57, 58, 59, 60, 61,
          62, 63, 64, 65, 66, 67,
          68, 69, 70, 71, 72, 73
     */
    for (j=0; j<4; j++) for (i=0; i<6; i++) varT[j][i] = j*6+i + 50;

    /* write varT to the NC variable in the matrix transposed way */
    start[0]  = 0; start[1]  = 0;
    count[0]  = 6; count[1]  = 4;
    stride[0] = 1; stride[1] = 1;
    imap[0]   = 1; imap[1]   = 6;   /* would be {4, 1} if not transposing */
    if (rank > 0) count[0] = count[1] = 0;

    if (cmode & NC_NETCDF4) {
        err = ncmpi_put_varm_schar_all(ncid, varid, start, count, stride,
                                       imap, &varT[0][0]); CHECK_ERR
        nerrs += check_write_contents(&varT[0][0]);
    }
    else {
        /* test nonblocking API */
        err = ncmpi_iput_varm_schar(ncid, varid, start, count, stride,
                                    imap, &varT[0][0], &req); CHECK_ERR
        err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERR
        err = status; CHECK_ERR
        nerrs += check_write_contents(&varT[0][0]);

        /* test when stride == NULL and imap != NULL */
        err = ncmpi_iput_varm_schar(ncid, varid, start, count, NULL,
                                    imap, &varT[0][0], &req); CHECK_ERR
        err = ncmpi_wait_all(ncid, 1, &req, &status); CHECK_ERR
        err = status; CHECK_ERR
        nerrs += check_write_contents(&varT[0][0]);

        /* test blocking API */
        err = ncmpi_put_varm_schar_all(ncid, varid, start, count, stride, imap,
                                       &varT[0][0]); CHECK_ERR
        nerrs += check_write_contents(&varT[0][0]);

        /* test when stride == NULL and imap != NULL */
        err = ncmpi_put_varm_schar_all(ncid, varid, start, count, NULL, imap,
                                       &varT[0][0]); CHECK_ERR
        nerrs += check_write_contents(&varT[0][0]);
    }

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char filename[256], *hint_value;
    int err, nerrs=0, rank, nprocs, bb_enabled=0;

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

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for get/put varm ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

#ifdef DEBUG
    if (nprocs > 1 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n", argv[0]);
#endif

    /* check whether burst buffering is enabled */
    if (inq_env_hint("nc_burst_buf", &hint_value)) {
        if (strcasecmp(hint_value, "enable") == 0) bb_enabled = 1;
        free(hint_value);
    }

    nerrs += tst_fmt(filename, 0);
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += tst_fmt(filename, NC_NETCDF4);
        nerrs += tst_fmt(filename, NC_NETCDF4 | NC_CLASSIC_MODEL);
#endif
    }
    nerrs += tst_fmt(filename, NC_64BIT_DATA);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

