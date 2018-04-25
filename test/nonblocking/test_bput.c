/*********************************************************************
 *
 *  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/*  $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <assert.h>
#include <pnetcdf.h>

#include <testutils.h>

#define FILE_NAME "testfile.nc"

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    int i, j, ncid, dimid[2], varid, err, nerrs=0, rank, nprocs;
    int req[2], status[2];
    float  var[4][6];
    char filename[256];
    MPI_Offset bufsize,  start[2], count[2], stride[2], imap[2];
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#ifdef DEBUG
    if (nprocs > 1 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n", argv[0]);
#endif

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for bput API ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    MPI_Info_create(&info);
    /* MPI_Info_set(info, "romio_pvfs2_posix_write","enable"); */

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER | NC_64BIT_DATA, info, &ncid); CHECK_ERR
    MPI_Info_free(&info);

    /* define a variable of a 6 x 4 integer array in the nc file */
    err = ncmpi_def_dim(ncid, "Y", 6, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 4, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT64, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* set the contents of the write buffer var, a 4 x 6 float array
          50, 51, 52, 53, 54, 55,
          56, 57, 58, 59, 60, 61,
          62, 63, 64, 65, 66, 67,
          68, 69, 70, 71, 72, 73
     */
    for (j=0; j<4; j++) for (i=0; i<6; i++) var[j][i] = 50.5 + j*6+i;

    /* bufsize must be max of data type converted before and after */
    bufsize = 4*6*sizeof(long long);
    err = ncmpi_buffer_attach(ncid, bufsize); CHECK_ERR

    /* write var to the NC variable in the matrix transposed way */
    count[0]  = 6; count[1]  = 2;
    stride[0] = 1; stride[1] = 1;
    imap[0]   = 1; imap[1]   = 6;   /* would be {4, 1} if not transposing */

    if (rank > 0) /* non-root processes just participate the call */
        count[0] = count[1] = 0;

    /* write the first two columns of the NC variable in the matrix transposed way */
    start[0]  = 0; start[1]  = 0;
    err = ncmpi_bput_varm_float(ncid, varid, start, count, stride, imap, &var[0][0], &req[0]); CHECK_ERR

    /* check if write buffer contents have been altered */
    for (j=0; j<4; j++)
        for (i=0; i<6; i++) {
            if (var[j][i] != 50.5 + j*6+i) {
                printf("Error at line %d in %s: put buffer[%d][%d]=%f altered, should be %f\n",
                       __LINE__,__FILE__,j,i,var[j][i],50.5+j*6+i);
                nerrs++;
            }
        }

    /* write the second two columns of the NC variable in the matrix transposed way */
    start[0]  = 0; start[1]  = 2;
    err = ncmpi_bput_varm_float(ncid, varid, start, count, stride, imap, &var[2][0], &req[1]); CHECK_ERR

    /* check if write buffer contents have been altered */
    for (j=0; j<4; j++)
        for (i=0; i<6; i++) {
            if (var[j][i] != 50.5 + j*6+i) {
                printf("Error at line %d in %s: put buffer[%d][%d]=%f altered, should be %f\n",
                       __LINE__,__FILE__,j,i,var[j][i],50.5+j*6+i);
                nerrs++;
            }
        }

    err = ncmpi_wait_all(ncid, 2, req, status); CHECK_ERR

    /* check each bput status */
    for (i=0; i<2; i++) {
        err = status[i];
        CHECK_ERR
    }

    err = ncmpi_buffer_detach(ncid); CHECK_ERR

    /* the output from command "ncmpidump -v var test.nc" should be:
           var =
            50, 56, 62, 68,
            51, 57, 63, 69,
            52, 58, 64, 70,
            53, 59, 65, 71,
            54, 60, 66, 72,
            55, 61, 67, 73 ;
     */

    /* check if the contents of write buffer have been altered (should not be) */
    for (j=0; j<4; j++) {
        for (i=0; i<6; i++) {
            if (var[j][i] != 50.5+j*6+i) {
                /* this error is a pnetcdf internal error, if occurs */
                printf("Error at line %d in %s: put buffer[%d][%d]=%f altered, should be %f\n",
                       __LINE__,__FILE__,j,i,var[j][i],50.5+j*6+i);
                nerrs++;
                break;
            }
        }
    }
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

