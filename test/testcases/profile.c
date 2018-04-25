/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/*
 * This program is to be used by PnetCDF developers to profile the sequence of
 * internal subroutine calls (mainly MPI calls).
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> /* getopt() */
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 2
#define NX 5
#define TRC(x) if(verbose) printf("%d: ---- before %s() ----\n",rank,#x);err=x

static int verbose;

static int test_vara(int ncid, int *varid)
{
    int rank, err, nerrs=0, buf[2][NX];
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    start[0] = 0; start[1] = NX*rank;
    count[0] = 2; count[1] = NX;

    TRC(ncmpi_put_vara_int_all)(ncid, varid[1], start, count, &buf[0][0]); CHECK_ERR
    TRC(ncmpi_rename_att)(ncid, varid[0], "att_name", "att_NAME"); CHECK_ERR
    TRC(ncmpi_put_vara_int_all)(ncid, varid[0], start, count, &buf[0][0]); CHECK_ERR
    TRC(ncmpi_get_vara_int_all)(ncid, varid[1], start, count, &buf[0][0]); CHECK_ERR
    TRC(ncmpi_get_vara_int_all)(ncid, varid[0], start, count, &buf[0][0]); CHECK_ERR

    TRC(ncmpi_begin_indep_data)(ncid); CHECK_ERR
    TRC(ncmpi_put_vara_int)(ncid, varid[1], start, count, &buf[0][0]); CHECK_ERR
    TRC(ncmpi_put_vara_int)(ncid, varid[0], start, count, &buf[0][0]); CHECK_ERR
    TRC(ncmpi_rename_att)(ncid, varid[0], "att_NAME", "att_name"); CHECK_ERR
    TRC(ncmpi_get_vara_int)(ncid, varid[1], start, count, &buf[0][0]); CHECK_ERR
    TRC(ncmpi_get_vara_int)(ncid, varid[0], start, count, &buf[0][0]); CHECK_ERR
    TRC(ncmpi_end_indep_data)(ncid); CHECK_ERR

    return nerrs;
}

static int test_ivara(int ncid, int *varid)
{
    int rank, err, nerrs=0, buf1[2][NX], buf2[2][NX], req[2], st[2];
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    start[0] = 0; start[1] = NX*rank;
    count[0] = 2; count[1] = NX;

    TRC(ncmpi_iput_vara_int)(ncid, varid[1], start, count, &buf1[0][0], &req[0]); CHECK_ERR
    TRC(ncmpi_wait_all)(ncid, 1, req, st); CHECK_ERR
    TRC(ncmpi_iput_vara_int)(ncid, varid[0], start, count, &buf2[0][0], &req[1]); CHECK_ERR
    TRC(ncmpi_wait_all)(ncid, 1, req, st); CHECK_ERR
    TRC(ncmpi_iget_vara_int)(ncid, varid[1], start, count, &buf1[0][0], &req[0]); CHECK_ERR
    TRC(ncmpi_iget_vara_int)(ncid, varid[0], start, count, &buf2[0][0], &req[1]); CHECK_ERR
    TRC(ncmpi_wait_all)(ncid, 2, req, st); CHECK_ERR

    TRC(ncmpi_begin_indep_data)(ncid); CHECK_ERR
    TRC(ncmpi_iput_vara_int)(ncid, varid[1], start, count, &buf1[0][0], &req[0]); CHECK_ERR
    TRC(ncmpi_iput_vara_int)(ncid, varid[0], start, count, &buf2[0][0], &req[1]); CHECK_ERR
    TRC(ncmpi_wait)(ncid, 2, req, st); CHECK_ERR
    TRC(ncmpi_iget_vara_int)(ncid, varid[1], start, count, &buf1[0][0], &req[0]); CHECK_ERR
    TRC(ncmpi_iget_vara_int)(ncid, varid[0], start, count, &buf2[0][0], &req[1]); CHECK_ERR
    TRC(ncmpi_wait)(ncid, 2, req, st); CHECK_ERR
    TRC(ncmpi_end_indep_data)(ncid); CHECK_ERR

    return nerrs;
}

static int test_vard(int ncid, int *varid)
{
    int          rank, nprocs, err, nerrs=0, i, buf[NY+4][NX+4];
    int          array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    MPI_Offset   start[2], count[2];
    MPI_Datatype buftype, rec_filetype, fix_filetype;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    start[0] = 0; start[1] = NX*rank;
    count[0] = 2; count[1] = NX;

    /* create a buftype with ghost cells on each side */
    array_of_sizes[0]    = count[0]+4;
    array_of_sizes[1]    = count[1]+4;
    array_of_subsizes[0] = count[0];
    array_of_subsizes[1] = count[1];
    array_of_starts[0]   = 2;
    array_of_starts[1]   = 2;
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C,
                             MPI_INT, &buftype);
    MPI_Type_commit(&buftype);

    /* create a file type for the fixed-size variable */
    array_of_sizes[0]    = 2;
    array_of_sizes[1]    = NX*nprocs;
    array_of_subsizes[0] = count[0];
    array_of_subsizes[1] = count[1];
    array_of_starts[0]   = start[0];
    array_of_starts[1]   = start[1];
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C,
                             MPI_INT, &fix_filetype);
    MPI_Type_commit(&fix_filetype);

    /* create a file type for the record variable */
    int *array_of_blocklengths=(int*) malloc(count[0]*sizeof(int));
    MPI_Aint *array_of_displacements=(MPI_Aint*) malloc(count[0]*sizeof(MPI_Aint));
    MPI_Offset recsize;
    err = ncmpi_inq_recsize(ncid, &recsize);
    for (i=0; i<count[0]; i++) {
        array_of_blocklengths[i] = count[1];
        array_of_displacements[i] = start[1]*sizeof(int) + recsize * i;
    }
    MPI_Type_create_hindexed(2, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &rec_filetype);
    MPI_Type_commit(&rec_filetype);
    free(array_of_blocklengths);
    free(array_of_displacements);

    TRC(ncmpi_put_vard_all)(ncid, varid[0], rec_filetype, &buf[0][0], 1, buftype); CHECK_ERR
    TRC(ncmpi_rename_var)(ncid, varid[0], "rec_VAR"); CHECK_ERR
    TRC(ncmpi_put_vard_all)(ncid, varid[1], fix_filetype, &buf[0][0], 1, buftype); CHECK_ERR
    TRC(ncmpi_rename_var)(ncid, varid[0], "rec_var"); CHECK_ERR

    TRC(ncmpi_begin_indep_data)(ncid); CHECK_ERR
    TRC(ncmpi_put_vard)(ncid, varid[0], rec_filetype, &buf[0][0], 1, buftype); CHECK_ERR
    TRC(ncmpi_rename_var)(ncid, varid[0], "rec_VAR"); CHECK_ERR
    TRC(ncmpi_put_vard)(ncid, varid[1], fix_filetype, &buf[0][0], 1, buftype); CHECK_ERR
    TRC(ncmpi_rename_var)(ncid, varid[0], "rec_var"); CHECK_ERR
    TRC(ncmpi_end_indep_data)(ncid); CHECK_ERR

    MPI_Type_free(&rec_filetype);
    MPI_Type_free(&fix_filetype);
    MPI_Type_free(&buftype);
    return nerrs;
}

static int test_varn(int ncid)
{
    int i, rank, nprocs, err, nerrs=0, num_reqs, dimids[2], varid, *buffer;
    MPI_Offset **starts, **counts;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* pick arbitrary numbers of requests for 4 processes */
    num_reqs = 0;
    if (rank == 0)      num_reqs = 4;
    else if (rank == 1) num_reqs = 6;
    else if (rank == 2) num_reqs = 5;
    else if (rank == 3) num_reqs = 4;

    starts    = (MPI_Offset**) malloc(num_reqs *    sizeof(MPI_Offset*));
    counts    = (MPI_Offset**) malloc(num_reqs *    sizeof(MPI_Offset*));
    starts[0] = (MPI_Offset*)  calloc(num_reqs * 2, sizeof(MPI_Offset));
    counts[0] = (MPI_Offset*)  calloc(num_reqs * 2, sizeof(MPI_Offset));
    for (i=1; i<num_reqs; i++) {
        starts[i] = starts[i-1] + 2;
        counts[i] = counts[i-1] + 2;
    }

    /* assign arbitrary starts and counts */
    const int y=0, x=1;
    if (rank == 0) {
        starts[0][y] = 0; starts[0][x] = 5; counts[0][y] = 1; counts[0][x] = 2;
        starts[1][y] = 1; starts[1][x] = 0; counts[1][y] = 1; counts[1][x] = 1;
        starts[2][y] = 2; starts[2][x] = 6; counts[2][y] = 1; counts[2][x] = 2;
        starts[3][y] = 3; starts[3][x] = 0; counts[3][y] = 1; counts[3][x] = 3;
        /* rank 0 is writing the followings: ("-" means skip)
                  -  -  -  -  -  0  0  -  -  -
                  0  -  -  -  -  -  -  -  -  -
                  -  -  -  -  -  -  0  0  -  -
                  0  0  0  -  -  -  -  -  -  -
         */
    } else if (rank ==1) {
        starts[0][y] = 0; starts[0][x] = 3; counts[0][y] = 1; counts[0][x] = 2;
        starts[1][y] = 0; starts[1][x] = 8; counts[1][y] = 1; counts[1][x] = 2;
        starts[2][y] = 1; starts[2][x] = 5; counts[2][y] = 1; counts[2][x] = 2;
        starts[3][y] = 2; starts[3][x] = 0; counts[3][y] = 1; counts[3][x] = 2;
        starts[4][y] = 2; starts[4][x] = 8; counts[4][y] = 1; counts[4][x] = 2;
        starts[5][y] = 3; starts[5][x] = 4; counts[5][y] = 1; counts[5][x] = 3;
        /* rank 1 is writing the followings: ("-" means skip)
                  -  -  -  1  1  -  -  -  1  1
                  -  -  -  -  -  1  1  -  -  -
                  1  1  -  -  -  -  -  -  1  1
                  -  -  -  -  1  1  1  -  -  -
         */
    } else if (rank ==2) {
        starts[0][y] = 0; starts[0][x] = 7; counts[0][y] = 1; counts[0][x] = 1;
        starts[1][y] = 1; starts[1][x] = 1; counts[1][y] = 1; counts[1][x] = 3;
        starts[2][y] = 1; starts[2][x] = 7; counts[2][y] = 1; counts[2][x] = 3;
        starts[3][y] = 2; starts[3][x] = 2; counts[3][y] = 1; counts[3][x] = 1;
        starts[4][y] = 3; starts[4][x] = 3; counts[4][y] = 1; counts[4][x] = 1;
        /* rank 2 is writing the followings: ("-" means skip)
                  -  -  -  -  -  -  -  2  -  -
                  -  2  2  2  -  -  -  2  2  2
                  -  -  2  -  -  -  -  -  -  -
                  -  -  -  2  -  -  -  -  -  -
         */
    } else if (rank ==3) {
        starts[0][y] = 0; starts[0][x] = 0; counts[0][y] = 1; counts[0][x] = 3;
        starts[1][y] = 1; starts[1][x] = 4; counts[1][y] = 1; counts[1][x] = 1;
        starts[2][y] = 2; starts[2][x] = 3; counts[2][y] = 1; counts[2][x] = 3;
        starts[3][y] = 3; starts[3][x] = 7; counts[3][y] = 1; counts[3][x] = 3;
        /* rank 3 is writing the followings: ("-" means skip)
                  3  3  3  -  -  -  -  -  -  -
                  -  -  -  -  3  -  -  -  -  -
                  -  -  -  3  3  3  -  -  -  -
                  -  -  -  -  -  -  -  3  3  3
         */
    }
    buffer = (int*) malloc(4*10 * sizeof(int));

    TRC(ncmpi_redef)(ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "M",  4, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "N", 10, &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimids, &varid); CHECK_ERR
    TRC(ncmpi_enddef)(ncid); CHECK_ERR

    TRC(ncmpi_begin_indep_data)(ncid); CHECK_ERR
    TRC(ncmpi_put_varn_int)(ncid, varid, num_reqs, starts, counts, buffer); CHECK_ERR
    TRC(ncmpi_get_varn_int)(ncid, varid, num_reqs, starts, counts, buffer); CHECK_ERR
    TRC(ncmpi_end_indep_data)(ncid); CHECK_ERR

    TRC(ncmpi_put_varn_int_all)(ncid, varid, num_reqs, starts, counts, buffer); CHECK_ERR
    TRC(ncmpi_get_varn_int_all)(ncid, varid, num_reqs, starts, counts, buffer); CHECK_ERR

    free(buffer);
    free(starts[0]);
    free(counts[0]);
    free(starts);
    free(counts);
    return nerrs;
}

static int test_ivarn(int ncid)
{
    int i, j, rank, nprocs, err, nerrs=0, num_reqs[4], dimids[2], varid[4], *buffer[4];
    int req[5], st[5];
    MPI_Offset **starts[4], **counts[4];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* pick arbitrary numbers of requests for 4 processes */
    num_reqs[0] = 4;
    num_reqs[1] = 6;
    num_reqs[2] = 5;
    num_reqs[3] = 4;

    for (i=0; i<4; i++) {
        starts[i]    = (MPI_Offset**) malloc(num_reqs[i] *    sizeof(MPI_Offset*));
        counts[i]    = (MPI_Offset**) malloc(num_reqs[i] *    sizeof(MPI_Offset*));
        starts[i][0] = (MPI_Offset*)  calloc(num_reqs[i] * 2, sizeof(MPI_Offset));
        counts[i][0] = (MPI_Offset*)  calloc(num_reqs[i] * 2, sizeof(MPI_Offset));
        for (j=1; j<num_reqs[i]; j++) {
            starts[i][j] = starts[i][j-1] + 2;
            counts[i][j] = counts[i][j-1] + 2;
        }
    }

    /* assign arbitrary starts and counts */
    const int y=0, x=1;
    starts[0][0][y] = 0; starts[0][0][x] = 5; counts[0][0][y] = 1; counts[0][0][x] = 2;
    starts[0][1][y] = 1; starts[0][1][x] = 0; counts[0][1][y] = 1; counts[0][1][x] = 1;
    starts[0][2][y] = 2; starts[0][2][x] = 6; counts[0][2][y] = 1; counts[0][2][x] = 2;
    starts[0][3][y] = 3; starts[0][3][x] = 0; counts[0][3][y] = 1; counts[0][3][x] = 3;

    starts[1][0][y] = 0; starts[1][0][x] = 3; counts[1][0][y] = 1; counts[1][0][x] = 2;
    starts[1][1][y] = 0; starts[1][1][x] = 8; counts[1][1][y] = 1; counts[1][1][x] = 2;
    starts[1][2][y] = 1; starts[1][2][x] = 5; counts[1][2][y] = 1; counts[1][2][x] = 2;
    starts[1][3][y] = 2; starts[1][3][x] = 0; counts[1][3][y] = 1; counts[1][3][x] = 2;
    starts[1][4][y] = 2; starts[1][4][x] = 8; counts[1][4][y] = 1; counts[1][4][x] = 2;
    starts[1][5][y] = 3; starts[1][5][x] = 4; counts[1][5][y] = 1; counts[1][5][x] = 3;

    starts[2][0][y] = 0; starts[2][0][x] = 7; counts[2][0][y] = 1; counts[2][0][x] = 1;
    starts[2][1][y] = 1; starts[2][1][x] = 1; counts[2][1][y] = 1; counts[2][1][x] = 3;
    starts[2][2][y] = 1; starts[2][2][x] = 7; counts[2][2][y] = 1; counts[2][2][x] = 3;
    starts[2][3][y] = 2; starts[2][3][x] = 2; counts[2][3][y] = 1; counts[2][3][x] = 1;
    starts[2][4][y] = 3; starts[2][4][x] = 3; counts[2][4][y] = 1; counts[2][4][x] = 1;

    starts[3][0][y] = 0; starts[3][0][x] = 0; counts[3][0][y] = 1; counts[3][0][x] = 3;
    starts[3][1][y] = 1; starts[3][1][x] = 4; counts[3][1][y] = 1; counts[3][1][x] = 1;
    starts[3][2][y] = 2; starts[3][2][x] = 3; counts[3][2][y] = 1; counts[3][2][x] = 3;
    starts[3][3][y] = 3; starts[3][3][x] = 7; counts[3][3][y] = 1; counts[3][3][x] = 3;

    for (i=0; i<4; i++) {
        buffer[i] = (int*) malloc(4*10 * sizeof(int));
        for (j=0; j<4*10; j++) buffer[i][j] = rank+100;
    }

    TRC(ncmpi_redef)(ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "M",  4, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "N", 10, &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 2, dimids, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 2, dimids, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var2", NC_INT, 2, dimids, &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var3", NC_INT, 2, dimids, &varid[3]); CHECK_ERR
    TRC(ncmpi_enddef)(ncid); CHECK_ERR

    for (i=0; i<3; i++) {
        j = (nprocs > 1) ? (i + rank) % nprocs : i;
        TRC(ncmpi_iput_varn_int)(ncid, varid[j], num_reqs[j], starts[j], counts[j], buffer[j], &req[i]); CHECK_ERR
    }
    TRC(ncmpi_wait_all)(ncid, 3, req, st); CHECK_ERR

    j = (nprocs > 1) ? (3 + rank) % nprocs : 3;
    TRC(ncmpi_iput_varn_int)(ncid, varid[j], num_reqs[j], starts[j], counts[j], buffer[j], &req[3]); CHECK_ERR
    for (i=0; i<3; i++) {
        j = (nprocs > 1) ? (i + rank) % nprocs : i;
        TRC(ncmpi_iget_varn_int)(ncid, varid[j], num_reqs[j], starts[j], counts[j], buffer[j], &req[i]); CHECK_ERR
    }
    TRC(ncmpi_wait_all)(ncid, 4, req, st); CHECK_ERR
    if (err != NC_NOERR) {
        for (i=0; i<4; i++) {
            if (st[i] != NC_NOERR) {
                printf("Error at line %d in %s: st[%d] %s\n",
                __FILE__,__LINE__,i,ncmpi_strerror(st[i]));
            }
        }
    }

    for (i=0; i<4; i++) {
        free(buffer[i]);
        free(starts[i][0]);
        free(counts[i][0]);
        free(starts[i]);
        free(counts[i]);
    }
    return nerrs;
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {

    extern int optind;
    char   filename[256];
    int    i, nerrs=0, err, rank, nprocs, ncid, varid[2], dimids[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for profiling ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = verbose = 0;
    if (rank == 0) {
        /* get command-line arguments */
        while ((i = getopt(argc, argv, "vh")) != EOF) {
            switch(i) {
                case 'v': verbose = 1;
                          break;
                case 'h':
                default:  printf("Usage: %s [-v] [filename]\n",argv[0]);
                          err = 1;
            }
        }
        if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
        else                      snprintf(filename, 256, "%s", argv[optind]);
    }
    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (err == 1) {
        MPI_Finalize();
        return 1;
    }
    MPI_Bcast(&verbose, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    /* create a new file for write */
    TRC(ncmpi_create)(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR
    if (verbose) printf("%d: ---- after ncmpi_create\n",rank);

    /* define a 2D array */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X",       NX*nprocs,    &dimids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "rec_var", NC_INT, 2, dimids, &varid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y",       2,            &dimids[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "fix_var", NC_INT, 2, dimids, &varid[1]); CHECK_ERR

    /* add attributes to the variable */
    err = ncmpi_put_att_text(ncid, varid[0], "att_name", 14, "attribute text"); CHECK_ERR
    TRC(ncmpi_enddef)(ncid); CHECK_ERR

    // nerrs += test_vara(ncid, varid);
    // nerrs += test_ivara(ncid, varid);
    // nerrs += test_vard(ncid, varid);
    // nerrs += test_varn(ncid);
    nerrs += test_ivarn(ncid);

    TRC(ncmpi_close)(ncid); CHECK_ERR

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

