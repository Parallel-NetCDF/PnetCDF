/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests the vard API.
 * The write buffer is a 2D array of size NY x NX
 * The MPI data type for the buffer is defined by swapping the 1st and 2nd
 * rows of the array using a buftype constructed by MPI_Type_create_hindex().
 * It also writes a fixed-size variable using a buftype constructed by
 * MPI_Type_create_subarray(). Both record and foxed-size variables are read
 * back using various filetypes and buftypes and check the contents.
 *
 * The expected results from the output file contents are:
 * (when running on 1 MPI process)
 *
 *  % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 *           REC_DIM = UNLIMITED ; // (2 currently)
 *           X = 5 ;
 *           FIX_DIM = 2 ;
 *    variables:
 *           int rec_var(REC_DIM, X) ;
 *           int dummy_rec(REC_DIM, X) ;
 *           int fix_var(FIX_DIM, X) ;
 *    data:
 *
 *    rec_var =
 *      10, 11, 12, 13, 14,
 *      0, 1, 2, 3, 4 ;
 *
 *    dummy_rec =
 *      0, 0, 0, 0, 0,
 *      0, 0, 0, 0, 0 ;
 *
 *    fix_var =
 *      10, 11, 12, 13, 14,
 *      0, 1, 2, 3, 4 ;
 * }
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 2
#define NX 5

#define CHECK_VALUE_PERMUTED { \
    for (j=0; j<count[0]; j++) { \
        int val = (j == 0) ? 1 : 0; \
        val = rank * 100 + val * 10; \
        for (i=0; i<count[1]; i++) \
            if (buf[j][i] != val+i) { \
                printf("line %d: expecting buf[%d][%d]=%d but got %d\n",__LINE__,j,i,val+i,buf[j][i]); \
                nerrs++; \
                goto fn_exit; \
            } \
    } \
} \

#define CHECK_VALUE(buf) { \
    for (j=0; j<count[0]; j++) { \
        for (i=0; i<count[1]; i++) \
            if (buf[j][i] != rank*100+j*10+i) { \
                printf("line %d: expecting buf[%d][%d]=%d but got %d\n",__LINE__,j,i,rank*100+j*10+i,(int)buf[j][i]); \
                nerrs++; \
                goto fn_exit; \
            } \
    } \
}

static
int get_var_and_verify(int ncid,
                       int coll_io,
                       int varid,
                       MPI_Offset *start,
                       MPI_Offset *count,
                       int **buf,
                       MPI_Datatype buftype,
                       MPI_Datatype ghost_buftype,
                       MPI_Datatype filetype)
{
    int i, j, rank, err, *ncbuf, nerrs=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ncbuf = (int *) malloc(sizeof(int) * (count[0]+4) * (count[1]+4));

    /* clear the contents of the read buffer */
    for (j=0; j<count[0]; j++) for (i=0; i<count[1]; i++) buf[j][i] = -1;

    /* read back using regular vara API */
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid, start, count, buf[0]);
    else
        err = ncmpi_get_vara_int(ncid, varid, start, count, buf[0]);
    CHECK_ERR

    /* check if the contents of buf are expected */
    CHECK_VALUE_PERMUTED

    /* clear the contents of the read buffer */
    for (j=0; j<count[0]; j++) for (i=0; i<count[1]; i++) buf[j][i] = -1;

    /* read back using flexible vara API */
    if (coll_io)
        err = ncmpi_get_vara_all(ncid, varid, start, count, buf[1], 1, buftype);
    else
        err = ncmpi_get_vara(ncid, varid, start, count, buf[1], 1, buftype);
    CHECK_ERR

    /* check if the contents of buf are expected */
    CHECK_VALUE(buf)

    /* clear the contents of the read buffer */
    for (j=0; j<count[0]; j++) for (i=0; i<count[1]; i++) buf[j][i] = -1;

    /* read back using vard API and permuted buftype */
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid, filetype, buf[1], 1, buftype);
    else
        err = ncmpi_get_vard(ncid, varid, filetype, buf[1], 1, buftype);
    CHECK_ERR

    /* check if the contents of buf are expected */
    CHECK_VALUE(buf)

    /* clear the contents of the read buffer */
    for (j=0; j<count[0]; j++) for (i=0; i<count[1]; i++) buf[j][i] = -1;

    /* read back using vard API and no buftype */
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid, filetype, buf[0], 0, MPI_DATATYPE_NULL);
    else
        err = ncmpi_get_vard(ncid, varid, filetype, buf[0], 0, MPI_DATATYPE_NULL);
    CHECK_ERR

    /* check if the contents of buf are expected */
    CHECK_VALUE_PERMUTED

    /* clear the contents of the read buffer */
    for (i=0; i<(count[0]+4)*(count[1]+4); i++) ncbuf[i] = -1;

    /* read back using ghost buftype */
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid, filetype, ncbuf, 1, ghost_buftype);
    else
        err = ncmpi_get_vard(ncid, varid, filetype, ncbuf, 1, ghost_buftype);
    CHECK_ERR

    for (j=0; j<count[0]; j++) {
        for (i=0; i<count[1]; i++)
            if (buf[j][i] != ncbuf[(j+2)*(count[1]+4)+(i+2)]) {
                printf("Error at line %d in %s: expecting ncbuf[%d][%d]=%d but got %d\n",
                       __LINE__,__FILE__,j,i,buf[j][i],ncbuf[(j+2)*(count[1]+4)+(i+2)]);
                nerrs++;
            }
    }

fn_exit:
    free(ncbuf);

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int          i, j, err, ncid, varid0, varid1, varid2, dimids[2], nerrs=0;
    int          rank, nprocs, blocklengths[2], **buf, *bufptr;
    int          array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    int          buftype_size, expected_put_size, fmt;
    float        **flt_buf=NULL, *flt_bufptr;
    double       **dbl_buf=NULL, *dbl_bufptr;
    MPI_Offset   start[2], count[2], header_size, put_size, new_put_size;
    MPI_Aint     a0, a1, disps[2];
    MPI_Datatype buftype, ghost_buftype, rec_filetype, fix_filetype;
    MPI_Datatype flt_buftype, dbl_buftype;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* construct various MPI derived data types */

    buf = (int**)malloc(sizeof(int*) * NY);
    buf[0] = (int*)malloc(sizeof(int) * NY * NX);
    for (i=1; i<NY; i++) buf[i] = buf[i-1] + NX;

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

    /* create a file type for the fixed-size variable */
    array_of_sizes[0] = 2;
    array_of_sizes[1] = NX*nprocs;
    array_of_subsizes[0] = (int)count[0];
    array_of_subsizes[1] = (int)count[1];
    array_of_starts[0] = (int)start[0];
    array_of_starts[1] = (int)start[1];
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C,
                             MPI_INT, &fix_filetype);
    MPI_Type_commit(&fix_filetype);

    /* create a buftype with ghost cells on each side */
    array_of_sizes[0] = (int)count[0]+4;
    array_of_sizes[1] = (int)count[1]+4;
    array_of_subsizes[0] = (int)count[0];
    array_of_subsizes[1] = (int)count[1];
    array_of_starts[0] = 2;
    array_of_starts[1] = 2;
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C,
                             MPI_INT, &ghost_buftype);
    MPI_Type_commit(&ghost_buftype);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for write */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X",       NX*nprocs,    &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "rec_var", NC_INT, 2, dimids, &varid0); CHECK_ERR
    err = ncmpi_def_var(ncid, "dummy_rec", NC_BYTE, 2, dimids, &varid2); CHECK_ERR
    err = ncmpi_def_dim(ncid, "FIX_DIM", 2, &dimids[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "fix_var", NC_INT, 2, dimids, &varid1); CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR /* enable fill mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* fill 2 records with default fill values */
    err = ncmpi_fill_var_rec(ncid, varid0, 0); CHECK_ERR
    err = ncmpi_fill_var_rec(ncid, varid0, 1); CHECK_ERR
    err = ncmpi_fill_var_rec(ncid, varid2, 0); CHECK_ERR
    err = ncmpi_fill_var_rec(ncid, varid2, 1); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* create a file type for the record variable */
    int *array_of_blocklengths=(int*) malloc(sizeof(int) * count[0]);
    MPI_Aint *array_of_displacements=(MPI_Aint*) malloc(sizeof(MPI_Aint) * count[0]);
    MPI_Offset recsize;
    err = ncmpi_inq_recsize(ncid, &recsize);
    for (i=0; i<count[0]; i++) {
        array_of_blocklengths[i] = (int)count[1];
        array_of_displacements[i] = (int)start[1]*sizeof(int) + recsize * i;
    }
    MPI_Type_create_hindexed(2, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &rec_filetype);
    MPI_Type_commit(&rec_filetype);

    /* initialize the contents of the array */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = rank*100 + j*10 + i;

    /* NULL argument test */
    err = ncmpi_inq_put_size(ncid, NULL); CHECK_ERR
    err = ncmpi_inq_get_size(ncid, NULL); CHECK_ERR

    /* get header size and put size by far */
    err = ncmpi_inq_header_size(ncid, &header_size); CHECK_ERR
    err = ncmpi_inq_put_size(ncid, &put_size); CHECK_ERR

    /* write the record variable */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid0, rec_filetype, bufptr, 1, buftype);
    else
        err = ncmpi_put_vard(ncid, varid0, rec_filetype, bufptr, 1, buftype);
    CHECK_ERR

    /* check if put_size is correctly reported */
    err = ncmpi_inq_put_size(ncid, &new_put_size); CHECK_ERR
    MPI_Type_size(buftype, &buftype_size);
    err = ncmpi_inq_format(ncid, &fmt); CHECK_ERR
    expected_put_size = buftype_size;

    /* for writing a record variable, root process will update numrec to the
     * file header, However, because the first 2 records have been filled
     * above, root process need not write to file header.
    if (rank == 0) expected_put_size += (fmt == NC_FORMAT_CDF5) ? 8 : 4;
     */
    if (expected_put_size != new_put_size - put_size) {
        printf("Error at line %d in %s: unexpected put size ("OFFFMT") reported, expecting %d\n",
               __LINE__,__FILE__,new_put_size-put_size, expected_put_size);
        nerrs++;
    }

    /* check if the contents of buf are altered */
    CHECK_VALUE(buf)

    /* check if root process can write to file header in data mode */
    err = ncmpi_rename_var(ncid, varid0, "rec_VAR"); CHECK_ERR

    err = ncmpi_inq_put_size(ncid, &put_size); CHECK_ERR

    /* write the fixed-size variable */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid1, fix_filetype, bufptr, 1, buftype);
    else
        err = ncmpi_put_vard(ncid, varid1, fix_filetype, bufptr, 1, buftype);
    CHECK_ERR

    /* check if put_size is correctly reported */
    err = ncmpi_inq_put_size(ncid, &new_put_size); CHECK_ERR
    expected_put_size = buftype_size;
    if (expected_put_size != new_put_size - put_size) {
        printf("Error at line %d in %s: unexpected put size ("OFFFMT") reported, expecting %d\n",
               __LINE__,__FILE__,new_put_size-put_size, expected_put_size);
        nerrs++;
    }

    /* check if the contents of buf are altered */
    CHECK_VALUE(buf)

    /* check if root process can write to file header in data mode */
    err = ncmpi_rename_var(ncid, varid0, "rec_var"); CHECK_ERR

    /* test the same routines in independent data mode */
    err = ncmpi_begin_indep_data(ncid); CHECK_ERR
    err = ncmpi_put_vard(ncid, varid0, rec_filetype, bufptr, 1, buftype); CHECK_ERR
    CHECK_VALUE(buf)
    err = ncmpi_rename_var(ncid, varid0, "rec_VAR"); CHECK_ERR
    err = ncmpi_put_vard(ncid, varid1, fix_filetype, bufptr, 1, buftype); CHECK_ERR
    CHECK_VALUE(buf)
    err = ncmpi_rename_var(ncid, varid0, "rec_var"); CHECK_ERR
    err = ncmpi_end_indep_data(ncid); CHECK_ERR

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_WRITE, info, &ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_varid(ncid, "rec_var", &varid0); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "fix_var", &varid1); CHECK_ERR

    nerrs += get_var_and_verify(ncid, coll_io, varid0, start, count, buf, buftype, ghost_buftype, rec_filetype);
    nerrs += get_var_and_verify(ncid, coll_io, varid1, start, count, buf, buftype, ghost_buftype, fix_filetype);

    /* test type conversion from float to int */
    flt_buf = (float**)malloc(sizeof(float*) * NY);
    flt_buf[0] = (float*)malloc(sizeof(float) * NY * NX);
    for (i=1; i<NY; i++) flt_buf[i] = flt_buf[i-1] + NX;
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        flt_buf[j][i] = rank*100.0 + j*10.0 + i;

    /* construct an MPI derived data type for swapping 1st row with 2nd row */
    blocklengths[0] = blocklengths[1] = NX;
    MPI_Get_address(flt_buf[1], &a0);
    MPI_Get_address(flt_buf[0], &a1);
    disps[0] = 0;
    disps[1] = a1 - a0;
    flt_bufptr = flt_buf[1];
    err = MPI_Type_create_hindexed(2, blocklengths, disps, MPI_FLOAT, &flt_buftype);
    if (err != MPI_SUCCESS) printf("MPI error MPI_Type_create_hindexed\n");
    MPI_Type_commit(&flt_buftype);

    /* write the record variable with type conversion from float to int */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid0, rec_filetype, flt_bufptr, 1, flt_buftype);
    else
        err = ncmpi_put_vard(ncid, varid0, rec_filetype, flt_bufptr, 1, flt_buftype);
    CHECK_ERR
    CHECK_VALUE(flt_buf)
    nerrs += get_var_and_verify(ncid, coll_io, varid0, start, count, buf, buftype, ghost_buftype, rec_filetype);

    /* read the record variable with type conversion from int to float */
    for (j=0; j<NY*NX; j++) flt_buf[0][j] = -1;
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid0, rec_filetype, flt_bufptr, 1, flt_buftype);
    else
        err = ncmpi_get_vard(ncid, varid0, rec_filetype, flt_bufptr, 1, flt_buftype);
    CHECK_ERR
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) {
        float expected = rank*100.0 + j*10.0 + i;
        if (flt_buf[j][i] != expected) {
            printf("Error at line %d in %s: expecting flt_buf[%d][%d]=%.1f but got %.1f\n",
            __LINE__,__FILE__,j,i,expected,flt_buf[j][i]);
            nerrs++;
        }
    }

    /* write the fixed-size variable with type conversion from float to int */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid1, fix_filetype, flt_bufptr, 1, flt_buftype);
    else
        err = ncmpi_put_vard(ncid, varid1, fix_filetype, flt_bufptr, 1, flt_buftype);
    CHECK_ERR
    CHECK_VALUE(flt_buf)
    nerrs += get_var_and_verify(ncid, coll_io, varid1, start, count, buf, buftype, ghost_buftype, fix_filetype);

    /* read the fixed-size variable with type conversion from int to float */
    for (j=0; j<NY*NX; j++) flt_buf[0][j] = -1;
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid1, fix_filetype, flt_bufptr, 1, flt_buftype);
    else
        err = ncmpi_get_vard(ncid, varid1, fix_filetype, flt_bufptr, 1, flt_buftype);
    CHECK_ERR

    for (j=0; j<NY; j++) for (i=0; i<NX; i++) {
        float expected = rank*100.0 + j*10.0 + i;
        if (flt_buf[j][i] != expected) {
            printf("Error at line %d in %s: expecting flt_buf[%d][%d]=%.1f but got %.1f\n",
            __LINE__,__FILE__,j,i,expected,flt_buf[j][i]);
            nerrs++;
        }
    }

    /* test type conversion from double to int */
    dbl_buf = (double**)malloc(sizeof(double*) * NY);
    dbl_buf[0] = (double*)malloc(sizeof(double) * NY * NX);
    for (i=1; i<NY; i++) dbl_buf[i] = dbl_buf[i-1] + NX;
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        dbl_buf[j][i] = rank*100.0 + j*10.0 + i;

    /* construct an MPI derived data type for swapping 1st row with 2nd row */
    blocklengths[0] = blocklengths[1] = NX;
    MPI_Get_address(dbl_buf[1], &a0);
    MPI_Get_address(dbl_buf[0], &a1);
    disps[0] = 0;
    disps[1] = a1 - a0;
    dbl_bufptr = dbl_buf[1];
    err = MPI_Type_create_hindexed(2, blocklengths, disps, MPI_DOUBLE, &dbl_buftype);
    if (err != MPI_SUCCESS) printf("MPI error MPI_Type_create_hindexed\n");
    MPI_Type_commit(&dbl_buftype);

    /* write the record variable with type conversion from double to int */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid0, rec_filetype, dbl_bufptr, 1, dbl_buftype);
    else
        err = ncmpi_put_vard(ncid, varid0, rec_filetype, dbl_bufptr, 1, dbl_buftype);
    CHECK_ERR
    CHECK_VALUE(dbl_buf)
    nerrs += get_var_and_verify(ncid, coll_io, varid0, start, count, buf, buftype, ghost_buftype, rec_filetype);

    /* read the record variable with type conversion from int to double */
    for (j=0; j<NY*NX; j++) dbl_buf[0][j] = -1;
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid0, rec_filetype, dbl_bufptr, 1, dbl_buftype);
    else
        err = ncmpi_get_vard(ncid, varid0, rec_filetype, dbl_bufptr, 1, dbl_buftype);
    CHECK_ERR
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) {
        double expected = rank*100.0 + j*10.0 + i;
        if (dbl_buf[j][i] != expected) {
            printf("Error at line %d in %s: expecting dbl_buf[%d][%d]=%.1f but got %.1f\n",
            __LINE__,__FILE__,j,i,expected,dbl_buf[j][i]);
            nerrs++;
        }
    }

    /* write the fixed-size variable with type conversion from double to int */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid1, fix_filetype, dbl_bufptr, 1, dbl_buftype);
    else
        err = ncmpi_put_vard(ncid, varid1, fix_filetype, dbl_bufptr, 1, dbl_buftype);
    CHECK_ERR
    CHECK_VALUE(dbl_buf)
    nerrs += get_var_and_verify(ncid, coll_io, varid1, start, count, buf, buftype, ghost_buftype, fix_filetype);

    /* read the fixed-size variable with type conversion from int to double */
    for (j=0; j<NY*NX; j++) dbl_buf[0][j] = -1;
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid1, fix_filetype, dbl_bufptr, 1, dbl_buftype);
    else
        err = ncmpi_get_vard(ncid, varid1, fix_filetype, dbl_bufptr, 1, dbl_buftype);
    CHECK_ERR
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) {
        double expected = rank*100.0 + j*10.0 + i;
        if (dbl_buf[j][i] != expected) {
            printf("Error at line %d in %s: expecting dbl_buf[%d][%d]=%.1f but got %.1f\n",
            __LINE__,__FILE__,j,i,expected,dbl_buf[j][i]);
            nerrs++;
        }
    }

    /* test data type that need no conversion */
    MPI_Type_free(&ghost_buftype);
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C,
                             MPI_SIGNED_CHAR, &ghost_buftype);
    MPI_Type_commit(&ghost_buftype);

    MPI_Type_free(&rec_filetype);

    /* construct rec_filetype for dummy_rec of type NC_BYTE */
    for (i=0; i<count[0]; i++) {
        array_of_blocklengths[i] = (int)count[1];
        array_of_displacements[i] = start[1]*sizeof(signed char) + recsize * i;
    }
    MPI_Type_create_hindexed(2, array_of_blocklengths, array_of_displacements,
                             MPI_SIGNED_CHAR, &rec_filetype);
    MPI_Type_commit(&rec_filetype);

    signed char *schar_buf;
    schar_buf = (signed char*) malloc(array_of_sizes[0]*array_of_sizes[1]);
    for (i=0; i<array_of_sizes[0]*array_of_sizes[1]; i++) schar_buf[i] = i + rank*10;

    /* write to dummy_rec */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid2, rec_filetype, schar_buf, 1, ghost_buftype);
    else
        err = ncmpi_put_vard(ncid, varid2, rec_filetype, schar_buf, 1, ghost_buftype);
    CHECK_ERR

    /* read from dummy_rec */
    for (i=0; i<array_of_sizes[0]*array_of_sizes[1]; i++) schar_buf[i] = -1;
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid2, rec_filetype, schar_buf, 1, ghost_buftype);
    else
        err = ncmpi_get_vard(ncid, varid2, rec_filetype, schar_buf, 1, ghost_buftype);
    CHECK_ERR
    for (i=0; i<array_of_sizes[0]; i++)
    for (j=0; j<array_of_sizes[1]; j++) {
        signed char expected = i*array_of_sizes[1] + j + rank*10;
        if (i<2 || i >= 4)         expected = -1;
        else if (j<2 || j >= NX+2) expected = -1;
        if (schar_buf[i*array_of_sizes[1]+j] != expected) {
            printf("Error at line %d in %s: expecting schar_buf[%d][%d]=%d but got %d\n",
            __LINE__,__FILE__,i,j,expected,schar_buf[i*array_of_sizes[1]+j]);
            nerrs++;
        }
    }
    free(schar_buf);

fn_exit:
    MPI_Type_free(&rec_filetype);
    MPI_Type_free(&fix_filetype);
    MPI_Type_free(&buftype);
    MPI_Type_free(&flt_buftype);
    MPI_Type_free(&dbl_buftype);
    MPI_Type_free(&ghost_buftype);
    free(array_of_blocklengths);
    free(array_of_displacements);
    free(buf[0]); free(buf);
    if (flt_buf != NULL) {
        free(flt_buf[0]);
        free(flt_buf);
    }
    if (dbl_buf != NULL) {
        free(dbl_buf[0]);
        free(dbl_buf);
    }

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
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 0; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "vard APIs", opt, test_io);

    MPI_Finalize();

    return err;
}
