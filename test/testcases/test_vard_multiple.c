/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program tests the vard API for writing/reading 2 consecutive variables
 * in a single put/get call.
 * Two datatypes are constructed for individual variables first using
 * MPI_Type_create_subarray and MPI_Type_create_hindexed. The two datatypes
 * are then concatenated into a single filetype.
 *
 * The expected results from the output file contents are:
 * (when running on 1 MPI process)
 *
 *  % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 *           REC_DIM = UNLIMITED ; // (2 currently)
 *           Y = 2 ;
 *           X = 5 ;
 *    variables:
 *           int fix_var0(Y, X) ;
 *           int fix_var1(Y, X) ;
 *           int rec_var2(REC_DIM, Y, X) ;
 *           int rec_var3(REC_DIM, Y, X) ;
 *    data:
 *
 *    fix_var0 =
 *      0, 1, 2, 3, 4,
 *      10, 11, 12, 13, 14 ;
 *
 *    fix_var1 =
 *      1000, 1001, 1002, 1003, 1004,
 *      1010, 1011, 1012, 1013, 1014 ;
 *
 *    rec_var2 =
 *      0, 0, 0, 0, 0,
 *      0, 0, 0, 0, 0,
 *      2000, 2001, 2002, 2003, 2004,
 *      2010, 2011, 2012, 2013, 2014 ;
 *
 *    rec_var3 =
 *      0, 0, 0, 0, 0,
 *      0, 0, 0, 0, 0,
 *      3000, 3001, 3002, 3003, 3004,
 *      3010, 3011, 3012, 3013, 3014 ;
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

#define CHECK_VALUE(buf,base) { \
    for (j=0; j<NY; j++) { \
        for (i=0; i<NX; i++) \
            if ((buf)[j*NX+i] != (base)+rank*100+j*10+i) { \
                printf("line %d: expecting buf[%d*NX+%d]=%d but got %d\n",\
                       __LINE__,j,i,(base)+rank*100+j*10+i,(buf)[j*NX+i]); \
                nerrs++; \
            } \
    } \
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {

    char         filename[256];
    int          i, j, err, ncid, varid[4], dimids[3], nerrs=0, unlimit_dimid;
    int          rank, nprocs, *buf[2];
    int          array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    int          array_of_blocklengths[NY];
    MPI_Offset   len, recsize, start[2], count[2], offset[2];
    MPI_Aint     a0, a1, array_of_displacements[NY];
    MPI_Datatype buftype, vtype[2], filetype;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
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
        sprintf(cmd_str, "*** TESTING C   %s for vard to 2 variables ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    buf[0] = (int*)malloc(NY * NX * sizeof(int));
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        buf[0][j*NX+i] = rank*100 + j*10 + i;

    buf[1] = (int*)malloc(NY * NX * sizeof(int));
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        buf[1][j*NX+i] = 1000 + rank*100 + j*10 + i;

    /* create a new NetCDF file */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL,
                       &ncid); CHECK_ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y",       NY,           &dimids[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X",       NX*nprocs,    &dimids[2]); CHECK_ERR

    /* define 2D fixed-size variables */
    err = ncmpi_def_var(ncid, "fix_var0", NC_INT, 2, dimids+1, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "fix_var1", NC_INT, 2, dimids+1, &varid[1]); CHECK_ERR
    /* define 3D record variables */
    err = ncmpi_def_var(ncid, "rec_var2", NC_INT, 3, dimids, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "rec_var3", NC_INT, 3, dimids, &varid[3]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* MPI_DATATYPE_NULL means this is a zero-length request */
    err = ncmpi_put_vard_all(ncid, varid[0], MPI_DATATYPE_NULL, NULL, 0,
                             MPI_INT); CHECK_ERR

    /* when filetype is MPI_DATATYPE_NULL, buftype is ignored */
    err = ncmpi_put_vard_all(ncid, varid[0], MPI_DATATYPE_NULL, NULL, 0,
                             MPI_DATATYPE_NULL); CHECK_ERR

    /* bufcount is ignored when buftype is MPI_DATATYPE_NULL */
    err = ncmpi_put_vard_all(ncid, varid[0], MPI_INT, buf[0], -1,
                             MPI_DATATYPE_NULL); CHECK_ERR

    /* filetype and buftype are matched */
    err = ncmpi_put_vard_all(ncid, varid[0], MPI_INT, buf[0], 1,
                             MPI_INT); CHECK_ERR

    /* create a datatype contains 2 different etypes */
    array_of_blocklengths[0]  = 1;
    array_of_blocklengths[1]  = 1;
    array_of_displacements[0] = 0;
    array_of_displacements[1] = sizeof(int);
    vtype[0] = MPI_INT;
    vtype[1] = MPI_FLOAT;
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &filetype);
    MPI_Type_commit(&filetype);
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &buftype);
    MPI_Type_commit(&buftype);

    /* test NC_EMULTITYPES */
    err = ncmpi_put_vard_all(ncid, varid[0], filetype, buf[0], 2, MPI_INT); EXP_ERR(NC_EMULTITYPES)
    err = ncmpi_put_vard_all(ncid, varid[0], MPI_INT,  buf[0], 2, buftype); EXP_ERR(NC_EMULTITYPES)
    /* test NC_ETYPESIZE_MISMATCH */
    err = ncmpi_put_vard_all(ncid, varid[0], MPI_INT, buf[0], 2, MPI_INT);   EXP_ERR(NC_EIOMISMATCH)
    err = ncmpi_put_vard_all(ncid, varid[0], MPI_INT, buf[0], 2, MPI_SHORT); EXP_ERR(NC_EIOMISMATCH)
    MPI_Type_free(&filetype);
    MPI_Type_free(&buftype);

    /* construct buftype: concatenate two separated buffers */
    array_of_blocklengths[0] = NY*NX;
    array_of_blocklengths[1] = NY*NX;
    array_of_displacements[0] = 0;
    MPI_Get_address(buf[0], &a0);
    MPI_Get_address(buf[1], &a1);
    array_of_displacements[1] = a1 - a0;
    vtype[0] = vtype[1] = MPI_INT;
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &buftype);
    MPI_Type_commit(&buftype);

    /* construct two MPI derived data types */
    start[0] = 0;  start[1] = NX*rank;
    count[0] = NY; count[1] = NX;

    /* create the first datatype using subarray */
    array_of_sizes[0]    = 2;
    array_of_sizes[1]    = NX*nprocs;
    array_of_subsizes[0] = count[0];
    array_of_subsizes[1] = count[1];
    array_of_starts[0]   = start[0];
    array_of_starts[1]   = start[1];
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C, MPI_INT, &vtype[0]);
    MPI_Type_commit(&vtype[0]);

    /* create the second datatype using hindexed */
    for (i=0; i<NY; i++)
        array_of_blocklengths[i] = NX;

    array_of_displacements[0] = NX*rank*sizeof(int);
    for (i=1; i<NY; i++)
        array_of_displacements[i] = array_of_displacements[i-1]
                                  + NX*nprocs*sizeof(int);

    MPI_Type_create_hindexed(NY, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &vtype[1]);
    MPI_Type_commit(&vtype[1]);

    /* concatenate two datatypes into filetype */
    array_of_blocklengths[0] = 1;
    array_of_blocklengths[1] = 1;
    array_of_displacements[0] = 0;

    /* calculate distance between 2 fixed-size variables */
    err = ncmpi_inq_varoffset(ncid, varid[0], &offset[0]); CHECK_ERR
    err = ncmpi_inq_varoffset(ncid, varid[1], &offset[1]); CHECK_ERR
    array_of_displacements[1] = offset[1] - offset[0];

    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &filetype);
    MPI_Type_commit(&filetype);

    /* write 2 consecutive fixed-size variables */
    err = ncmpi_put_vard_all(ncid, varid[0], filetype, buf[0], 1, buftype); CHECK_ERR

    /* check if the contents of write buf are altered */
    CHECK_VALUE(buf[0], 0)
    CHECK_VALUE(buf[1], 1000)

    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[0][j*NX+i] = -1;
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[1][j*NX+i] = -1;

    /* read back fixed-size variables */
    err = ncmpi_get_vard_all(ncid, varid[0], filetype, buf[0], 1, buftype); CHECK_ERR

    /* check the contents of read buf */
    CHECK_VALUE(buf[0], 0)
    CHECK_VALUE(buf[1], 1000)

    MPI_Type_free(&filetype);

    /* obtain record size */
    err = ncmpi_inq_recsize(ncid, &recsize); CHECK_ERR

    /* calculate distance between 2 record variables */
    err = ncmpi_inq_varoffset(ncid, varid[2], &offset[0]); CHECK_ERR
    err = ncmpi_inq_varoffset(ncid, varid[3], &offset[1]); CHECK_ERR
    array_of_displacements[1] = offset[1] - offset[0];

    /* write to 2nd record of two consecutive variables */
    array_of_displacements[0] += recsize;
    array_of_displacements[1] += recsize;

    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &filetype);
    MPI_Type_commit(&filetype);

    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        buf[0][j*NX+i] = 2000 + rank*100 + j*10 + i;

    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        buf[1][j*NX+i] = 3000 + rank*100 + j*10 + i;

    /* write 2 consecutive record variables */
    err = ncmpi_put_vard_all(ncid, varid[2], filetype, buf[0], 1, buftype); CHECK_ERR

    /* check if the contents of write buf are altered */
    CHECK_VALUE(buf[0], 2000)
    CHECK_VALUE(buf[1], 3000)

    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[0][j*NX+i] = -1;
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[1][j*NX+i] = -1;

    /* read back record variables */
    err = ncmpi_get_vard_all(ncid, varid[2], filetype, buf[0], 1, buftype); CHECK_ERR

    /* check the contents of read buf */
    CHECK_VALUE(buf[0], 2000)
    CHECK_VALUE(buf[1], 3000)

    /* write to 1st record of two consecutive variables */
    MPI_Type_free(&filetype);
    array_of_displacements[0] -= recsize;
    array_of_displacements[1] -= recsize;
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &filetype);
    MPI_Type_commit(&filetype);

    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[0][j*NX+i] = 0;
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[1][j*NX+i] = 0;

    /* write 2 consecutive record variables */
    err = ncmpi_put_vard_all(ncid, varid[2], filetype, buf[0], 1, buftype); CHECK_ERR

    /* check if the number of records is still 2 */
    err = ncmpi_inq_unlimdim(ncid, &unlimit_dimid); CHECK_ERR
    err = ncmpi_inq_dimlen(ncid, unlimit_dimid, &len); CHECK_ERR
    if (len != 2)
        printf("Error at line %d in %s: number of records should be 2 but got %lld\n",
        __LINE__,__FILE__,len);

    MPI_Type_free(&vtype[0]);
    MPI_Type_free(&vtype[1]);
    MPI_Type_free(&filetype);
    MPI_Type_free(&buftype);
    free(buf[0]); free(buf[1]);

    err = ncmpi_close(ncid); CHECK_ERR

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
