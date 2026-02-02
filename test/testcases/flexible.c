/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests the use of flexible API.
 * The write buffer is a 2D array of size NY x NX
 * The MPI data type for the buffer is defined by swapping the 1st and 2nd
 * rows of the array. It uses MPI_Type_create_hindex(). After the write, this
 * test reads back the array using regular and flexible get APIs (blocking and
 * nonblocking) and check the contents.
 *
 * The expected results from the output file contents are:
 * (when running on 1 MPI process)
 *
 *  % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 * 	   Y = UNLIMITED ; // (2 currently)
 * 	   X = 5 ;
 *    variables:
 * 	   int VAR(Y, X) ;
 *    data:
 *
 *    var =
 *      11, 11, 11, 11, 11,
 *      10, 10, 10, 10, 10 ;
 *    }
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 2
#define NX 70

#define INDEP_MODE 0
#define COLL_MODE 1

static int debug;

/*----< tst_io() >-----------------------------------------------------------*/
static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int          i, j, err, ncid, varid1, varid2, varid3, dimids[2];
    int          rank, nprocs, blocklengths[2], buf[NY][NX], *bufptr;
    int         *ncbuf, req, st, nerrs=0;
    int          array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    MPI_Offset   start[2], count[2];
    MPI_Aint     a0, a1, disps[2];
    MPI_Datatype buftype;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs,    &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_int", NC_INT, 2, dimids, &varid1); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_dbl", NC_DOUBLE, 2, dimids, &varid2); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_byte", NC_BYTE, 2, dimids, &varid3); CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR /* enable fill mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

#ifdef STRONGER_CONSISTENCY
    ncmpi_sync(ncid);
    MPI_Barrier(MPI_COMM_WORLD);
    ncmpi_sync(ncid);
#endif

    /* fill 2 records with default fill values */
    err = ncmpi_fill_var_rec(ncid, varid1, 0); CHECK_ERR
    err = ncmpi_fill_var_rec(ncid, varid1, 1); CHECK_ERR
    err = ncmpi_fill_var_rec(ncid, varid2, 0); CHECK_ERR
    err = ncmpi_fill_var_rec(ncid, varid2, 1); CHECK_ERR
    err = ncmpi_fill_var_rec(ncid, varid3, 0); CHECK_ERR
    err = ncmpi_fill_var_rec(ncid, varid3, 1); CHECK_ERR

#ifdef STRONGER_CONSISTENCY
    ncmpi_sync(ncid);
    MPI_Barrier(MPI_COMM_WORLD);
    ncmpi_sync(ncid);
#endif

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* initialize the contents of the array */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = j+rank+10;

    /* construct an MPI derived data type for swapping 1st row with 2nd row */
    blocklengths[0] = blocklengths[1] = NX;
    MPI_Get_address(buf[1], &a0);
    MPI_Get_address(buf[0], &a1);
    disps[0] = 0;
    disps[1] = a1 - a0;
    bufptr = buf[1];
    err = MPI_Type_create_hindexed(2, blocklengths, disps, MPI_INT, &buftype);
    if (err != MPI_SUCCESS) printf("MPI error MPI_Type_create_hindexed\n");
    MPI_Type_commit(&buftype);

    start[0] = 0; start[1] = NX*rank;
    count[0] = 2; count[1] = NX;
    if (debug) printf("put start="OFFFMT" "OFFFMT" count="OFFFMT" "OFFFMT"\n",start[0],start[1],count[0],count[1]);

    /* call flexible API */
    if (coll_io)
        err = ncmpi_put_vara_all(ncid, varid1, start, count, bufptr, 1, buftype);
    else
        err = ncmpi_put_vara(ncid, varid1, start, count, bufptr, 1, buftype);
    CHECK_ERR

    /* check if the contents of buf are altered */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        if (buf[j][i] != j+rank+10) {
            printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                   __FILE__,__LINE__,j,i,j+rank+10,buf[j][i]);
            return 1;
        }

    /* call flexible API that do type conversion */
    if (coll_io)
        err = ncmpi_put_vara_all(ncid, varid2, start, count, bufptr, 1, buftype);
    else
        err = ncmpi_put_vara(ncid, varid2, start, count, bufptr, 1, buftype);
    CHECK_ERR

    /* check if the contents of buf are altered */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        if (buf[j][i] != j+rank+10) {
            printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                   __FILE__,__LINE__,j,i,j+rank+10,buf[j][i]);
            return 1;
        }

    MPI_Type_free(&buftype);

    /* call flexible API that do no type conversion */
    signed char *schar_buf = (signed char *) malloc(NY*NX);
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) schar_buf[j*NX+i] = j+rank+10;
    disps[0] = 0;
    disps[1] = -NX;
    err = MPI_Type_create_hindexed(2, blocklengths, disps, MPI_SIGNED_CHAR, &buftype);
    if (err != MPI_SUCCESS) printf("MPI error MPI_Type_create_hindexed\n");
    MPI_Type_commit(&buftype);

    if (coll_io)
        err = ncmpi_put_vara_all(ncid, varid3, start, count, schar_buf+NX, 1, buftype);
    else
        err = ncmpi_put_vara(ncid, varid3, start, count, schar_buf+NX, 1, buftype);
    CHECK_ERR

    /* check if the contents of buf are altered */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        if (schar_buf[j*NX+i] != j+rank+10) {
            printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                   __FILE__,__LINE__,j,i,j+rank+10,schar_buf[j*NX+i]);
            return 1;
        }

    for (j=0; j<NY; j++) for (i=0; i<NX; i++) schar_buf[j*NX+i] = -1;

    if (coll_io)
        err = ncmpi_get_vara_all(ncid, varid3, start, count, schar_buf+NX, 1, buftype);
    else
        err = ncmpi_get_vara(ncid, varid3, start, count, schar_buf+NX, 1, buftype);
    CHECK_ERR

    /* check read contents */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        if (schar_buf[j*NX+i] != j+rank+10) {
            printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                   __FILE__,__LINE__,j,i,j+rank+10,schar_buf[j*NX+i]);
            return 1;
        }

    free(schar_buf);
    MPI_Type_free(&buftype);

    /* check if root process can write to file header in data mode */
    err = ncmpi_rename_var(ncid, varid1, "VAR"); CHECK_ERR

    /* Ensure all write data committed to file system before any rank open
     * the file for read.
     */
    ncmpi_sync(ncid);
    MPI_Barrier(MPI_COMM_WORLD);
    ncmpi_sync(ncid);

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_varid(ncid, "VAR", &varid1); CHECK_ERR

    /* initialize the contents of the array to a different value */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = -1;

    /* read back variable */
    start[0] = 0; start[1] = NX*rank;
    count[0] = 2; count[1] = NX;
    if (debug) printf("get start="OFFFMT" "OFFFMT" count="OFFFMT" "OFFFMT"\n",start[0],start[1],count[0],count[1]);

    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid1, start, count, buf[0]);
    else
        err = ncmpi_get_vara_int(ncid, varid1, start, count, buf[0]);
    CHECK_ERR

    /* check if the contents of buf are expected */
    for (j=0; j<2; j++) {
        int expect = (j == 0) ? 1 : 0;
        for (i=0; i<NX; i++)
            if (buf[j][i] != expect+rank+10) {
                printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                       __FILE__,__LINE__,j,i,expect+rank+10,buf[j][i]);
                return 1;
            }
    }

    /* initialize the contents of the array to a different value */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = -1;

    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid2, start, count, buf[0]);
    else
        err = ncmpi_get_vara_int(ncid, varid2, start, count, buf[0]);
    CHECK_ERR

    /* check if the contents of buf are expected */
    for (j=0; j<2; j++) {
        int expect = (j == 0) ? 1 : 0;
        for (i=0; i<NX; i++)
            if (buf[j][i] != expect+rank+10) {
                printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                       __FILE__,__LINE__,j,i,expect+rank+10,buf[j][i]);
                return 1;
            }
    }

    /* create a buftype with ghost cells on each side */
    ncbuf = (int *) malloc(sizeof(int) * (count[0]+4) * (count[1]+4));
    array_of_sizes[0] = (int)count[0]+4;
    array_of_sizes[1] = (int)count[1]+4;
    array_of_subsizes[0] = (int)count[0];
    array_of_subsizes[1] = (int)count[1];
    array_of_starts[0] = 2;
    array_of_starts[1] = 2;
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C,
                             MPI_INT, &buftype);
    MPI_Type_commit(&buftype);
    if (coll_io)
        err = ncmpi_get_vara_all(ncid, varid1, start, count, ncbuf, 1, buftype);
    else
        err = ncmpi_get_vara(ncid, varid1, start, count, ncbuf, 1, buftype);
    CHECK_ERR

    for (j=0; j<count[0]; j++) for (i=0; i<count[1]; i++) {
        int getValue = ncbuf[(j+2)*(count[1]+4)+(i+2)];
        if (buf[j][i] != getValue) {
            printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                   __FILE__,__LINE__,j,i,buf[j][i],getValue);
            return 1;
        }
    }

    for (i=0; i<(count[0]+4)*(count[1]+4); i++) ncbuf[i] = -1;
    if (coll_io)
        err = ncmpi_get_vara_all(ncid, varid2, start, count, ncbuf, 1, buftype);
    else
        err = ncmpi_get_vara(ncid, varid2, start, count, ncbuf, 1, buftype);
    CHECK_ERR

    for (j=0; j<count[0]; j++) for (i=0; i<count[1]; i++) {
        int getValue = ncbuf[(j+2)*(count[1]+4)+(i+2)];
        if (buf[j][i] != getValue) {
            printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                   __FILE__,__LINE__,j,i,buf[j][i],getValue);
            return 1;
        }
    }

    for (i=0; i<(count[0]+4)*(count[1]+4); i++) ncbuf[i] = -1;

    err = ncmpi_iget_vara(ncid, varid1, start, count, ncbuf, 1, buftype, &req); CHECK_ERR
    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, &st);
    else
        err = ncmpi_wait(ncid, 1, &req, &st);
    CHECK_ERR
    err = st; CHECK_ERR

    for (j=0; j<count[0]; j++) for (i=0; i<count[1]; i++) {
        int getValue = ncbuf[(j+2)*(count[1]+4)+(i+2)];
        if (buf[j][i] != getValue) {
            printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                   __FILE__,__LINE__,j,i,buf[j][i],getValue);
            return 1;
        }
    }

    for (i=0; i<(count[0]+4)*(count[1]+4); i++) ncbuf[i] = -1;

    err = ncmpi_iget_vara(ncid, varid2, start, count, ncbuf, 1, buftype, &req); CHECK_ERR
    if (coll_io)
        err = ncmpi_wait_all(ncid, 1, &req, &st);
    else
        err = ncmpi_wait(ncid, 1, &req, &st);
    CHECK_ERR
    err = st; CHECK_ERR

    for (j=0; j<count[0]; j++) for (i=0; i<count[1]; i++) {
        int getValue = ncbuf[(j+2)*(count[1]+4)+(i+2)];
        if (buf[j][i] != getValue) {
            printf("Error at %s line %d: expect buf[%d][%d]=%d but got %d\n",
                   __FILE__,__LINE__,j,i,buf[j][i],getValue);
            return 1;
        }
    }

    MPI_Type_free(&buftype);
    free(ncbuf);

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};

    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "flexible APIs", opt, test_io);

    MPI_Finalize();

    return err;
}
