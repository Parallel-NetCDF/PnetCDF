/*
 *  Copyright (C) 2012, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pnetcdf.h>

#define FILE_NAME "test.nc"

#define ERRCODE 2
#define ERR(e) {printf("Error at line %d: err=%d %s\n", __LINE__, e, ncmpi_strerror(e)); exit(ERRCODE);}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    int i, j, ncid, dimid[2], varid, retval, err=0, rank, nprocs;
    int req[2], status[2]; 
    MPI_Offset start[2];
    MPI_Offset count[2];
    MPI_Offset stride[2];
    MPI_Offset imap[2];
    float  var[4][6];
    MPI_Offset bufsize;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs > 1)
        printf("This test is designed to run on one process\n");
    if (rank > 0) {
        MPI_Finalize();
        return 0;
    }

    if (NC_NOERR != (retval = ncmpi_create(MPI_COMM_WORLD, FILE_NAME,
        NC_CLOBBER | NC_64BIT_DATA, MPI_INFO_NULL, &ncid)))
       ERR(retval);

    /* define a variable of a 6 x 4 integer array in the nc file */
    if (NC_NOERR != (retval = ncmpi_def_dim(ncid, "Y", 6, &dimid[0]))) ERR(retval);
    if (NC_NOERR != (retval = ncmpi_def_dim(ncid, "X", 4, &dimid[1]))) ERR(retval);
    if (NC_NOERR != (retval = ncmpi_def_var(ncid, "var", NC_INT64, 2, dimid, &varid)))
        ERR(retval);
    if (NC_NOERR != (retval = ncmpi_enddef(ncid))) ERR(retval);

    /* set the contents of the write buffer var, a 4 x 6 float array
          50, 51, 52, 53, 54, 55,
          56, 57, 58, 59, 60, 61,
          62, 63, 64, 65, 66, 67,
          68, 69, 70, 71, 72, 73
     */
    for (j=0; j<4; j++) for (i=0; i<6; i++) var[j][i] = j*6+i + 50;

    /* bufsize must be max of data type converted before and after */
    bufsize = 4*6*sizeof(long long);
    if (NC_NOERR != (retval = ncmpi_buffer_attach(ncid, bufsize)))
        ERR(retval);

    /* write var to the NC variable in the matrix transposed way */
    count[0]  = 6; count[1]  = 2;
    stride[0] = 1; stride[1] = 1;
    imap[0]   = 1; imap[1]   = 6;   /* would be {4, 1} if not transposing */

    /* write the first two columns of the NC variable in the matrix transposed way */
    start[0]  = 0; start[1]  = 0;
    if (NC_NOERR != (retval = ncmpi_bput_varm_float(ncid, varid, start, count, stride, imap, &var[0][0], &req[0])))
        ERR(retval);

    /* write the second two columns of the NC variable in the matrix transposed way */
    start[0]  = 0; start[1]  = 2;
    if (NC_NOERR != (retval = ncmpi_bput_varm_float(ncid, varid, start, count, stride, imap, &var[2][0], &req[1])))
        ERR(retval);

    if (NC_NOERR != (retval = ncmpi_wait_all(ncid, 2, req, status)))
        ERR(retval);

    /* check each bput status */
    for (i=0; i<2; i++)
        if (status[i] != NC_NOERR) ERR(status[i]);

    if (NC_NOERR != (retval = ncmpi_buffer_detach(ncid)))
        ERR(retval);

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
            if (var[j][i] != j*6+i + 50) {
#ifdef PRINT_ERR_ON_SCREEN
                /* this error is a pntecdf internal error, if occurs */
                printf("Error: bput_varm write buffer has been altered at j=%d i=%d\n",j,i);
#endif
                err++;
                break;
            }
        }
    }
    if (NC_NOERR != (retval = ncmpi_close(ncid))) ERR(retval);

    MPI_Finalize();

    char cmd_str[80];
    sprintf(cmd_str, "*** TESTING %s for bput API ", argv[0]);

    if (err)
        printf("%-66s ------ failed\n", cmd_str);
    else
        printf("%-66s ------ pass\n", cmd_str);
    return err;
}

