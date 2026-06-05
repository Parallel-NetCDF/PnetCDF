/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests the use of nonblocking API.
 * The write buffer is a 2D array of size NY x NX
 * It writes the 2nd row of the memory buffer to the 1st row of the variable
 * array in file. Then it writes the 1st row of the memory buffer to the
 * 2nd row of the variable array in file.
 *
 * The expected reults from the output file contents are:
 * (when running on 1 MPI process)
 *
 *  % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 *         Y = UNLIMITED ; // (2 currently)
 *         X = 5 ;
 *    variables:
 *         int VAR(Y, X) ;
 *    data:
 *
 *    var =
 *      1, 1, 1, 1, 1,
 *      0, 0, 0, 0, 0 ;
 *    }
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 4
#define NX 5

static
int tst_mode(const char *filename,
             int         mode,
             MPI_Info    info)
{
    int i, j, err, ncid, varid, dimids[2], req[2], st[2], nerrs=0;
    int rank, nprocs, buf[NY+1][NX];
    MPI_Offset start[2], count[2];

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);
    CHECK_FATAL_ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX,    &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (mode == MODE_INDEP) {
        err = ncmpi_sync(ncid); CHECK_ERR
    }

    /* initialize the contents of the array */
    for (j=0; j<NY+1; j++) for (i=0; i<NX; i++) buf[j][i] = j;

    start[0] = 2*rank; start[1] = 0;
    count[0] = 1;      count[1] = NX;

    /* call nonblocking API */
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf[1], &req[0]);
    CHECK_ERR

    start[0] += 1;
    err = ncmpi_iput_vara_int(ncid, varid, start, count, buf[0], &req[1]);
    CHECK_ERR

    st[0] = st[1] = NC_NOERR;

    if (mode == MODE_COLL) {
        err = ncmpi_wait_all(ncid, 2, req, st);
        CHECK_ERR
    }
    else {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
        err = ncmpi_wait(ncid, 2, req, st);
        CHECK_ERR
    }
    err = st[0]; CHECK_ERR
    err = st[1]; CHECK_ERR

    /* check if the contents of buf are altered */
    for (j=0; j<NY; j++)
        for (i=0; i<NX; i++)
            if (buf[j][i] != j) {
                printf("Error at line %d in %s: buf[%d][%d]=%d != %d\n",
                __LINE__,__FILE__,j,i,buf[j][i],j);
                nerrs++;
                goto fn_exit;
            }

    /* check if root process can write to file header in data mode */
    err = ncmpi_rename_var(ncid, varid, "VAR"); CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_FATAL_ERR

    err = ncmpi_inq_varid(ncid, "VAR", &varid); CHECK_ERR

    if (mode == MODE_INDEP) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR_ALL
    }

    /* initialize the contents of the array to a different value */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = -1;

    /* read back variable */
    start[0] = 2*rank; start[1] = 0;
    count[0] = 2;      count[1] = NX;
    if (mode == MODE_COLL)
        err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf[0]);
    else
        err = ncmpi_get_vara_int(ncid, varid, start, count, buf[0]);
    CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if the contents of buf are expected */
    for (j=0; j<2; j++) {
        int val = (j == 0) ? 1 : 0;
        for (i=0; i<NX; i++)
            if (buf[j][i] != val) {
                printf("Error: unexpected read buf[%d][%d]=%d, should be %d\n",
                       j,i,buf[j][i],val);
                nerrs++;
                goto fn_exit;
            }
    }

fn_exit:
    return nerrs;
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {

    char filename[256];
    int err, nerrs=0, rank;
    MPI_Info info;

    MPI_Init(&argc, &argv);
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
        sprintf(cmd_str, "*** TESTING C   %s for using ncmpi_iput_vara_int() ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    MPI_Info_create(&info);
    /* When using PVFS2, unexpected buffer value error message might occur.
     * This is due to  a possible bug in ADIOI_PVFS2_OldWriteStrided() when
     * filetype is contiguous and buftype is non-contiguous.
     * Fix: Add ROMIO hint to force MPI-IO to use POSIX I/O driver */
    /* MPI_Info_set(info, "romio_pvfs2_posix_write", "enable"); */

    /* disable internal buffering for small non-blocking APIs */
    MPI_Info_set(info, "nc_ibuf_size", "0");

    nerrs = tst_mode(filename, MODE_COLL,  MPI_INFO_NULL);
    if (nerrs > 0) goto err_out;

    nerrs = tst_mode(filename, MODE_INDEP, MPI_INFO_NULL);
    if (nerrs > 0) goto err_out;

    nerrs = tst_mode(filename, MODE_COLL, info);
    if (nerrs > 0) goto err_out;

    nerrs = tst_mode(filename, MODE_INDEP, info);
    if (nerrs > 0) goto err_out;

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has "OFFFMT" bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

err_out:
    MPI_Info_free(&info);

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}
