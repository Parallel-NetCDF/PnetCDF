/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests ncmpi_put_vard APIs for writing a record variable with
 * one record at a time. Check the number of records after each write. This is
 * to test the fix to bug reported by Jim Edwards in r3675.
 *
 *  % mpiexec -n 4 test_vard_rec
 *  % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 *           REC_DIM = UNLIMITED ; // (2 currently)
 *           X = 12 ;
 *    variables:
 *           int rec_var(REC_DIM, X) ;
 *    data:
 *
 *     rec_var =
 *       0, 1, 2, 100, 101, 102, 200, 201, 202, 300, 301, 302,
 *       10, 11, 12, 110, 111, 112, 210, 211, 212, 310, 311, 312 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h> /* strcpy(), strncpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 2
#define NX 3

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {

    char         filename[256];
    int          i, j, err, nerrs=0, ncid, varid, dimids[2], unlimit_dimid;
    int          rank, nprocs, array_of_blocklengths[2], buf[NY][NX];
    MPI_Offset   recsize, len;
    MPI_Aint     array_of_displacements[2];
    MPI_Datatype rec_filetype;

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
        sprintf(cmd_str, "*** TESTING C   %s for vard put on record var ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* create a new file for write */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL,
                       &ncid); CHECK_ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X",       NX*nprocs,    &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "rec_var", NC_INT, 2, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* initialize the contents of the array */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = rank*100 + j*10 + i;

    /* obtain record size: sum of single records of all record variables */
    err = ncmpi_inq_recsize(ncid, &recsize); CHECK_ERR
    err = ncmpi_inq_unlimdim(ncid, &unlimit_dimid); CHECK_ERR

    /* create a file type for writing 1st record */
    array_of_blocklengths[0] = NX;
    array_of_displacements[0] = NX*rank*sizeof(int);
    MPI_Type_create_hindexed(1, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &rec_filetype);
    MPI_Type_commit(&rec_filetype);

    /* write the record variable */
    err = ncmpi_put_vard_all(ncid, varid, rec_filetype, buf, NX, MPI_INT);
    CHECK_ERR

    MPI_Type_free(&rec_filetype);

    /* check if the number of records changed to 1 */
    err = ncmpi_inq_dimlen(ncid, unlimit_dimid, &len); CHECK_ERR
    if (len != 1)
        printf("Error at line %d in %s: number of records should be 1 but got %lld\n",
        __LINE__,__FILE__,len);

    /* create a file type for writing 2nd record */
    array_of_blocklengths[0] = NX;
    array_of_displacements[0] = NX*rank*sizeof(int) + recsize;
    MPI_Type_create_hindexed(1, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &rec_filetype);
    MPI_Type_commit(&rec_filetype);

    /* write the record variable */
    err = ncmpi_put_vard_all(ncid, varid, rec_filetype, buf[1], NX, MPI_INT);
    CHECK_ERR

    MPI_Type_free(&rec_filetype);

    /* check if the number of records changed to 2 */
    err = ncmpi_inq_dimlen(ncid, unlimit_dimid, &len); CHECK_ERR
    if (len != 2)
        printf("Error at line %d in %s: number of records should be 2 but got %lld\n",
        __LINE__,__FILE__,len);

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "rec_var", &varid); CHECK_ERR

    /* create a file type for writing 2nd record */
    array_of_blocklengths[0] = NX;
    array_of_blocklengths[1] = NX;
    array_of_displacements[0] = NX*rank*sizeof(int);
    array_of_displacements[1] = NX*rank*sizeof(int) + recsize;
    MPI_Type_create_hindexed(2, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &rec_filetype);
    MPI_Type_commit(&rec_filetype);

    /* reset contents of buf before read */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = -1;

    /* read back record variable */
    err = ncmpi_get_vard_all(ncid, varid, rec_filetype, buf, NX*2, MPI_INT);
    CHECK_ERR

    MPI_Type_free(&rec_filetype);

    /* check read contents */
    for (j=0; j<NY; j++)
    for (i=0; i<NX; i++) {
        if (buf[j][i] != rank*100 + j*10 + i) {
            printf("Error: read buf[%d][%d] expect %d but got %d\n",j,i,
                   rank*100 + j*10 + i, buf[j][i]);
            nerrs++;
            break;
        }
    }

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
