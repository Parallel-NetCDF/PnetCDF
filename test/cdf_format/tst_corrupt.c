/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program reads CDF files with corrupted header and check if the
 * expected error code can be returned.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_corrupt tst_corrupt.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_corrupt dir_name
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */

/* This program can also be used to test NetCDF.
 * Add #define TEST_NETCDF and compile with command:
 * gcc -I/netcdf/path/include last_large_var.c -o last_large_var -L/netcdf/path/lib -lnetcdf
 */
#ifdef TEST_NETCDF
#include <netcdf.h>
#include <netcdf_meta.h>
#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,nc_strerror(err)); \
    } \
}
#define EXP_ERR(exp) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expecting %d but got %d\n", \
        __LINE__,__FILE__,exp, err); \
    } \
}
#define FileOpen(a,b,c,d,e) nc_open(b,c,e)
#define MPI_Init(a,b)
#define MPI_Comm_rank(a,b)
#define MPI_Bcast(a,b,c,d,e)
#define MPI_Finalize()
#else
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>
#define FileOpen ncmpi_open
#endif

int main(int argc, char** argv) {
    char *bad_xtype[3] ={"bad_xtype.nc1",  "bad_xtype.nc2",  "bad_xtype.nc5"};
    char *bad_ndims[3] ={"bad_ndims.nc1",  "bad_ndims.nc2",  "bad_ndims.nc5"};
    char *bad_dimid[3] ={"bad_dimid.nc1",  "bad_dimid.nc2",  "bad_dimid.nc5"};
    char *bad_nattrs[3]={"bad_nattrs.nc1", "bad_nattrs.nc2", "bad_nattrs.nc5"};

    char filename[1024], dirname[512];
    int i, rank, err, ncid, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 2) {
        if (!rank) printf("Usage: %s [directory name]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    snprintf(dirname, 512, "%s", argv[1]);
    MPI_Bcast(dirname, 512, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str,
        "*** TESTING C   %s for checking corrupted file header ",
        basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* The 3 CDF files in bad_xtype[3] have been created on purpose to contain
     * a corrupted header entry "nc_type" with an invalid value larger than
     * the maximum allowed atomic data type. For CDF-1 and CDF-2, the maximum
     * atomic data type is NC_DOUBLE and for CDF-5, the maximum is NC_UINT64.
     * The corrupted entry is indicated below.
     *
     * CDF format specification:
     *     var = name nelems [dimid ...] vatt_list nc_type vsize begin
     *                                             ^^^^^^^
     *                                             corrupted
     */
    for (i=0; i<3; i++) {
        sprintf(filename, "%s/%s", dirname, bad_xtype[i]);
        err = FileOpen(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                         &ncid);
        EXP_ERR(NC_EBADTYPE)
    }

    /* The 3 CDF files in bad_ndims[3] have been created on purpose to contain
     * a corrupted header entry "number of dimensions" with a value larger than
     * NC_MAX_DIMS. The corrupted entry is indicated below.
     *
     * CDF format specification:
     *     dim_list = ABSENT | NC_DIMENSION  nelems  [dim ...]
     *                                       ^^^^^^
     *                                       corrupted
     */
    for (i=0; i<3; i++) {
        sprintf(filename, "%s/%s", dirname, bad_ndims[i]);
        err = FileOpen(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                         &ncid);
        EXP_ERR(NC_EMAXDIMS)
    }

    /* The 3 CDF files in bad_dimid[3] have been created on purpose to contain
     * a corrupted header entry "dimension ID" of a variable that is beyond
     * the number of dimensions defined in the file.
     *
     * CDF format specification:
     *     var = name nelems [dimid ...] vatt_list nc_type vsize begin
     *                        ^^^^^
     *                        corrupted
     */
    for (i=0; i<3; i++) {
        sprintf(filename, "%s/%s", dirname, bad_dimid[i]);
        err = FileOpen(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                         &ncid);
        EXP_ERR(NC_EBADDIM)
    }

    /* The 3 CDF files in bad_nattrs[3] have been created on purpose to contain
     * a corrupted header entry "number of attributes" with a value larger than
     * NC_MAX_ATTRS. The corrupted entry is indicated below.
     *
     * CDF format specification:
     *     att_list = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     *                                       ^^^^^^
     *                                       corrupted
     */
    for (i=0; i<3; i++) {
        sprintf(filename, "%s/%s", dirname, bad_nattrs[i]);
        err = FileOpen(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                         &ncid);
        EXP_ERR(NC_EMAXATTS)
    }

#ifdef TEST_NETCDF
    if (nerrs) printf("fail with %d mismatches\n",nerrs);
    else       printf("pass\n");
#else
    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0)
            ncmpi_inq_malloc_list();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }
#endif

    MPI_Finalize();
    return (nerrs > 0);
}

