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
 *
 * When setting NX to 3, below shows the expected file contents.
 *
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
#define NX 100

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int          i, j, err, nerrs=0, ncid, varid, dimids[2], unlimit_dimid;
    int          rank, nprocs, verbose, array_of_blocklengths[2], buf[NY][NX];
    MPI_Offset   recsize, len;
    MPI_Aint     array_of_displacements[2];
    MPI_Datatype rec_filetype;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 0;

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for write */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X",       NX*nprocs,    &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "rec_var", NC_INT, 2, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

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

    if (verbose)
        printf("%d: %s line %d: displacements=%ld blocklengths=%d\n",rank,
               __func__,__LINE__, array_of_displacements[0], array_of_blocklengths[0]);

    /* write the record variable */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid, rec_filetype, buf, NX, MPI_INT);
    else
        err = ncmpi_put_vard(ncid, varid, rec_filetype, buf, NX, MPI_INT);
    CHECK_ERR

    MPI_Type_free(&rec_filetype);

    /* check if the number of records changed to 1 */
    err = ncmpi_inq_dimlen(ncid, unlimit_dimid, &len); CHECK_ERR
    if (len != 1)
        printf("Error at line %d in %s: number of records should be 1 but got "OFFFMT"\n",
        __LINE__,__FILE__,len);

    /* create a file type for writing 2nd record */
    array_of_blocklengths[0] = NX;
    array_of_displacements[0] = NX*rank*sizeof(int) + recsize;
    MPI_Type_create_hindexed(1, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &rec_filetype);
    MPI_Type_commit(&rec_filetype);

    if (verbose)
        printf("%d: %s line %d: displacements=%ld blocklengths=%d\n",rank,
               __func__,__LINE__, array_of_displacements[0], array_of_blocklengths[0]);

    /* write the record variable */
    if (coll_io)
        err = ncmpi_put_vard_all(ncid, varid, rec_filetype, buf[1], NX, MPI_INT);
    else
        err = ncmpi_put_vard(ncid, varid, rec_filetype, buf[1], NX, MPI_INT);
    CHECK_ERR

    MPI_Type_free(&rec_filetype);

    /* check if the number of records changed to 2 */
    err = ncmpi_inq_dimlen(ncid, unlimit_dimid, &len); CHECK_ERR
    if (len != 2)
        printf("Error at line %d in %s: number of records should be 2 but got "OFFFMT"\n",
        __LINE__,__FILE__,len);

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the same file and read back for validate */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_varid(ncid, "rec_var", &varid); CHECK_ERR

    /* create a file type for writing 2nd record */
    array_of_blocklengths[0] = NX;
    array_of_blocklengths[1] = NX;
    array_of_displacements[0] = NX*rank*sizeof(int);
    array_of_displacements[1] = NX*rank*sizeof(int) + recsize;
    MPI_Type_create_hindexed(2, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &rec_filetype);
    MPI_Type_commit(&rec_filetype);

    if (verbose)
        printf("%d: %s line %d: displacements=%ld %ld blocklengths=%d %d\n",rank,
               __func__,__LINE__,
               array_of_displacements[0], array_of_displacements[1],
               array_of_blocklengths[0], array_of_blocklengths[1]);

    /* reset contents of buf before read */
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[j][i] = -1;

    /* read back record variable */
    if (coll_io)
        err = ncmpi_get_vard_all(ncid, varid, rec_filetype, buf, NX*2, MPI_INT);
    else
        err = ncmpi_get_vard(ncid, varid, rec_filetype, buf, NX*2, MPI_INT);
    CHECK_ERR

    MPI_Type_free(&rec_filetype);

    /* check read contents */
    for (j=0; j<NY; j++)
    for (i=0; i<NX; i++) {
        int exp = rank*100 + j*10 + i;
        if (buf[j][i] != exp) {
            printf("Error at %d: read buf[%d][%d] expect %d but got %d\n",
                  __LINE__,j,i, exp, buf[j][i]);
            nerrs++;
            break;
        }
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

    err = tst_main(argc, argv, "vard put on record var", opt, test_io);

    MPI_Finalize();

    return err;
}
