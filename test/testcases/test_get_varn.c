/*
 *  Copyright (C) 2022, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  Test ncmpi_get_varn_double_all() using E3SM-IO pattern.
 *  See Pull Request #90
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h> /* basename() */

#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#ifndef MPI_OFFSET
#define MPI_OFFSET MPI_LONG_LONG_INT
#endif

#define NDIMS 3

int main(int argc, char** argv)
{
    int i, j, rank, nprocs, err, nerrs = 0;
    int ncid, varid, num_reqs;
    double *buffer;
    float *fbuffer;
    MPI_Offset r_len, **starts = NULL, **counts = NULL;
    MPI_Offset st[3], ct[3];
    char filename[256];
    int dimids[3];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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
        sprintf(cmd_str, "*** TESTING C   %s for get_varn ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

#ifdef DEBUG
    if (nprocs != 4 && rank == 0)
        printf("Warning: %s is designed to run on 4 process\n",argv[0]);
#endif

/*    original test case read from a file like this:
 *    % ncdump -h lnfm.nc
 *        netcdf lnfm {
 *        dimensions:
 *          time = UNLIMITED ; // (3 currently)
 *          lat = 94 ;
 *          lon = 192 ;
 *        variables:
 *          float lnfm(time, lat, lon) ;
 *        }
 */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_DATA, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "lat", 94, &dimids[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "lon", 192, &dimids[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "lnfm", NC_FLOAT, 3, dimids, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    st[0] = rank*2;
    st[1] = 0;
    st[2] = 0;

    ct[0] = 2;
    ct[1] = 94;
    ct[2] = 192;
    float *scramble = (float*) calloc(ct[0]*ct[1]*ct[2], sizeof(float));

    err = ncmpi_put_vara_float_all(ncid, varid, st, ct, scramble); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    free(scramble);

    /* now we can finally exercise the read path of this record varable */

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* pick 2 requests for 4 processes */
    /* num_reqs = 1; => works fine*/
    num_reqs = 2;

    starts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
    counts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
    starts[0] = (MPI_Offset*)  calloc(num_reqs*NDIMS, sizeof(MPI_Offset));
    counts[0] = (MPI_Offset*)  calloc(num_reqs*NDIMS, sizeof(MPI_Offset));
    for (i = 1; i < num_reqs; i++) {
        starts[i] = starts[i - 1] + NDIMS;
        counts[i] = counts[i - 1] + NDIMS;
    }

    /* assign specific starts and counts */
    if (num_reqs > 0){
        if (rank == 1) {
            starts[0][0] = 0; starts[0][1] = 0;  starts[0][2] = 1;
            counts[0][0] = 1; counts[0][1] = 93; counts[0][2] = 1;
        }
    }

    if (num_reqs > 1){
        if (rank == 1) {
            starts[1][0] = 1; starts[1][1] = 0;  starts[1][2] = 1;
            counts[1][0] = 1; counts[1][1] = 93; counts[1][2] = 2;
        }
    }

    r_len = 0; /* total read length for this process */
    for (i = 0; i < num_reqs; i++) {
        MPI_Offset r_req_len = 1;
        for (j = 0; j < NDIMS; j++)
            r_req_len *= counts[i][j];
        r_len += r_req_len;
    }

    /* allocate I/O buffer */
    buffer = (double*) calloc(r_len, sizeof(double));
    fbuffer = (float*) calloc(r_len, sizeof(float));

    /* set the buffer pointers to different offsets to the I/O buffer */
    varid = 0; /* only one variable in lnfm.nc */
    err = ncmpi_get_varn_double_all(ncid, varid, num_reqs, starts, counts, buffer);
    /* err = ncmpi_get_varn_float_all(ncid, varid, num_reqs, starts, counts, fbuffer); */
    CHECK_ERR

    err = ncmpi_close(ncid);
    CHECK_ERR

#ifdef DEBUG
    if (rank == 3) {
        printf("Dumping some double type data read by rank 3 (count = 10) ...\n");
        for (i = 0; i < ((r_len > 10)?10:r_len); i++)
            printf("%lf, ", buffer[i]);
        printf("\n");
    }
#endif

    free(buffer);
    free(fbuffer);

    free(starts[0]);
    free(counts[0]);
    free(starts);
    free(counts);

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

