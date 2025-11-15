/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests entering define modes multiple times.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_multi_redefine tst_multi_redefine.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./tst_multi_redefine testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h>  /* basename() */
#include <unistd.h>  /* getopt() */

#include <pnetcdf.h>

#include <testutils.h>

#define NROUNDS 2
#define NY 10
#define NX 10
#define NVARS 2

static int verbose;

static int
tst_fmt(char *filename,
        int   cmode,
        int   nRounds,
        int   len_y,
        int   len_x)
{
    int i, j, k, rank, nprocs, ncid, err, nerrs=0, dimids[3];
    MPI_Offset old_hdr_size, new_hdr_size;
    MPI_Offset old_hdr_ext, new_hdr_ext;

    MPI_Info info=MPI_INFO_NULL;

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (verbose && rank == 0) {
        char *fmt;
        if (cmode == 0)
            fmt = "NC_CLASSIC";
        else if (cmode == NC_64BIT_OFFSET)
            fmt = "NC_64BIT_OFFSET";
        else if (cmode == NC_64BIT_DATA)
            fmt = "NC_64BIT_DATA";
        printf("\n---- Testing file format %s ----\n", fmt);
    }

    /* create a new file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERROUT

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimids[0]); CHECK_ERROUT
    err = ncmpi_def_dim(ncid, "Y",    len_y,        &dimids[1]); CHECK_ERROUT
    err = ncmpi_def_dim(ncid, "X",    len_x*nprocs, &dimids[2]); CHECK_ERROUT

    /* Default v_align of 512 is pretty bad choice, as almost every iteration
     * grows the file extension.
     */
    // err = ncmpi__enddef(ncid, 0, 4096, 0, 0);
    err = ncmpi_enddef(ncid);
    CHECK_ERROUT

    err = ncmpi_inq_header_size(ncid, &old_hdr_size); CHECK_ERROUT
    err = ncmpi_inq_header_extent(ncid, &old_hdr_ext); CHECK_ERROUT
    if (verbose && rank == 0)
        printf("Newly created file with header size %lld, extension %lld\n",
               old_hdr_size, old_hdr_ext);

    k = 0;
    for (i=0; i<nRounds; i++) {
        char str[64];
        int *fix_varids, *rec_varids;

        fix_varids = (int*) malloc(sizeof(int) * NVARS);
        rec_varids = (int*) malloc(sizeof(int) * NVARS);

        err = ncmpi_redef(ncid); CHECK_ERROUT

        /* add new fix-sized variables */
        for (j=0; j<NVARS; j++) {
            sprintf(str, "fix_var_%d", k);
            err = ncmpi_def_var(ncid, str, NC_DOUBLE, 2, dimids+1, &fix_varids[j]);
            CHECK_ERROUT
            sprintf(str, "attribute of variable %d", k);
            err = ncmpi_put_att_text(ncid, fix_varids[j], "attr", strlen(str), str); CHECK_ERROUT
            k++;
        }

        /* add new record variables */
        for (j=0; j<NVARS; j++) {
            sprintf(str, "rec_var_%d", k);
            err = ncmpi_def_var(ncid, str, NC_DOUBLE, 3, dimids, &rec_varids[j]);
            CHECK_ERROUT
            sprintf(str, "attribute of variable %d", k);
            err = ncmpi_put_att_text(ncid, rec_varids[j], "attr", strlen(str), str); CHECK_ERROUT
            k++;
        }

        err = ncmpi_enddef(ncid); CHECK_ERROUT
        err = ncmpi_inq_header_size(ncid, &new_hdr_size); CHECK_ERROUT
        err = ncmpi_inq_header_extent(ncid, &new_hdr_ext); CHECK_ERROUT
        if (verbose && rank == 0) {
            printf("Iternation %d: file header size grows from %lld to %lld\n",
                   i, old_hdr_size, new_hdr_size);
            if (new_hdr_ext > old_hdr_ext)
                printf("Iternation %d: file header extension grows from %lld to %lld\n",
                    i, old_hdr_ext, new_hdr_ext);
        }
        old_hdr_size = new_hdr_size;
        old_hdr_ext  = new_hdr_ext;

        /* write to new variables */
        MPI_Offset start[3], count[3];
        start[0] = 0;
        start[1] = 0;
        start[2] = rank * len_x;
        count[0] = 1;
        count[1] = len_y;
        count[2] = len_x;

        /* write to fix-sized variables */
        int *int_buf = (int*) malloc(sizeof(int) * count[0] * count[1]);
        for (j=0; j<count[0] * count[1]; j++) int_buf[j] = rank + j;

        for (j=0; j<NVARS; j++) {
            err = ncmpi_put_vara_int_all(ncid, fix_varids[j], start+1, count+1, int_buf);
            CHECK_ERROUT
        }
        free(int_buf);

        /* write to record variables */
        double *dbl_buf = (double*) malloc(sizeof(double) * count[0] * count[1]);
        for (j=0; j<count[0] * count[1]; j++) dbl_buf[j] = rank + j;

        for (j=0; j<NVARS; j++) {
            err = ncmpi_put_vara_double_all(ncid, rec_varids[j], start, count, dbl_buf);
            CHECK_ERROUT
        }
        free(dbl_buf);

        free(fix_varids);
        free(rec_varids);
    }

    err = ncmpi_close(ncid); CHECK_ERROUT

err_out:
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    return nerrs;
}

#define FILE_NAME "testfile.nc"

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTIONS]...[filename]\n"
    "       [-h] Print help\n"
    "       [-v]: enable verbose output mode (default: no)\n"
    "       [-n num]: number of rounds entering define mode (default: %d)\n"
    "       [-y num]: Y dimension size of global array (default: %d)\n"
    "       [-x num]: X dimension size of local array (default: %d)\n"
    "       [filename]: output netCDF file name (default: %s)\n";
    fprintf(stderr, help, argv0, FILE_NAME, NROUNDS, NY, NX);
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, fmt, rank, err, nerrs=0, cmode[3], nRounds, len_y, len_x;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 0;
    nRounds = NROUNDS;
    len_y = NY;
    len_x = NX;
    while ((i = getopt(argc, argv, "hvn:y:x:")) != EOF)
        switch(i) {
            case 'v': verbose = 1;
                      break;
            case 'n': nRounds = atoi(optarg);
                      break;
            case 'y': len_y = atoi(optarg);
                      break;
            case 'x': len_x = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, FILE_NAME);
    else                      snprintf(filename, 256, "%s", argv[optind]);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for entering define mode ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }
    cmode[0] = 0;
    cmode[1] = NC_64BIT_OFFSET;
    cmode[2] = NC_64BIT_DATA;

    for (fmt=0; fmt<3; fmt++) {
        nerrs += tst_fmt(filename, cmode[fmt], nRounds, len_y, len_x);
        if (nerrs > 0) goto err_out;
    }

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
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

