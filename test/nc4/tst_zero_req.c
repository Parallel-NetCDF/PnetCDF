/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* This program tests a HDF5 1.10.2 bug that causes test program to hang when
 * writing to and reading back a 2D record variable in collective mode with
 * some of the processes making zero-length requests.
 */

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define Y_LEN 2
#define X_LEN 4

static int
tst_fmt(char *filename, int cmode)
{
    int i, j, err, nerrs=0, nprocs, rank, ncid, dimid[2], varid, *buf;
    MPI_Offset start[2], count[2], num_records=0;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        printf("Error at %s:%d ncmpi_create %s\n",__FILE__,__LINE__,
               ncmpi_strerrno(err));
        return 1;
    }

    buf = (int*)malloc(Y_LEN*X_LEN*sizeof(int));
    for (i=0; i<Y_LEN*X_LEN; i++) buf[i] = rank*100+i;

    /* To test HDF5 1.10.2 bug, it must use a record variable, fixed-size
     * variable will not work. In addition, the varialble must be 2D. 1D
     * won't reveal HDF5 hanging problem.
     */
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", X_LEN*nprocs, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid); CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* all processes write data of size Y_LEN*X_LEN
     * Y_LEN must be at least 2. 1 won't reveal HDF5 hanging problem */
    start[0] = 0;
    start[1] = X_LEN * rank;
    count[0] = Y_LEN;
    count[1] = X_LEN;

    /* all zero-lenghth requests will fail HDF5 1.10.2 */
    count[0] = 0;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    count[0] = Y_LEN;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, dimid[0], &num_records); CHECK_ERR
    if (num_records != Y_LEN) {
        printf("Error at %s:%d: expect num_records=%d but got %lld\n",
               __FILE__,__LINE__,Y_LEN,num_records);
        nerrs++;
    }

    /* Must call ncmpi_sync() to cause HDF5 hanging problem */
    err = ncmpi_sync(ncid); CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, dimid[0], &num_records); CHECK_ERR
    if (num_records != Y_LEN) {
        printf("Error at %s:%d: expect num_records=%d but got %lld\n",
               __FILE__,__LINE__,Y_LEN,num_records);
        nerrs++;
    }

    /* test zero-length read from all processes */
    count[0] = 0;
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    /* only one process reads back what it wrote */
    for (i=0; i<Y_LEN*X_LEN; i++) buf[i] = -1;
    if (rank == 0) count[0] = Y_LEN;
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    if (rank == 0) {
        for (i=0; i<Y_LEN; i++) for (j=0; j<X_LEN; j++) {
            int k= i*X_LEN+j;
            int expect=k+(rank * 100);
            if (buf[k] != expect) {
                printf("Error at %s:%d: expect buf[%d]=%d but got %d\n",
                       __FILE__,__LINE__,k,expect,buf[k]);
                nerrs++;
                i = Y_LEN;
                break;
            }
        }
    }

    /* process i reads data written by process i+1 */
    for (i=0; i<Y_LEN*X_LEN; i++) buf[i] = -1;
    start[1] = X_LEN * ((rank+1) % nprocs);
    count[0] = Y_LEN;
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf); CHECK_ERR

    for (i=0; i<Y_LEN; i++) for (j=0; j<X_LEN; j++) {
        int k=i*X_LEN+j;
        int expect=k+((rank+1)%nprocs * 100);
        if (buf[k] != expect) {
            printf("Error at %s:%d: expect buf[%d]=%d but got %d\n",
                   __FILE__,__LINE__,k,expect,buf[k]);
            nerrs++;
            i = Y_LEN;
            break;
        }
    }

    /* expect program to hang at ncmpi_close() if HDF5 1.10.2 is used */
    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);

    return nerrs;
}

int main(int argc, char **argv) {
    char filename[256];
    int rank=0, nerrs=0;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
        sprintf(cmd_str, "*** TESTING C   %s for zero-length request ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* test all file formats separately */
    nerrs += tst_fmt(filename, 0);
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
/* #if defined(HDF5_VERSION) && (H5_VERS_MAJOR > 1 || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 10) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 4)) */
#if defined(HDF5_VER_GE_1_10_4) && HDF5_VER_GE_1_10_4 == 1
    /* this test requires a bug fix in HDF5 1.11 */
    /* Note NC_MPIIO is used in NetCDF 4.6.1 and earlier, but ignored in 4.6.2
     * and after. */
    nerrs += tst_fmt(filename, NC_MPIIO | NC_NETCDF4);
    nerrs += tst_fmt(filename, NC_MPIIO | NC_NETCDF4 | NC_CLASSIC_MODEL);
#endif
    nerrs += tst_fmt(filename, NC_64BIT_DATA);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    int err = ncmpi_inq_malloc_size(&malloc_size);
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
