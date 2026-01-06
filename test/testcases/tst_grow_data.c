/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This program tests adding new fix-sized and record variables by re-entering
 * the define mode, without growing file header.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

static int debug;

#define LON 100
#define LAT 100
#define NVARS 10

#define PRINT_VAR_OFF \
    for (i=0; i<nvars; i++) { \
        err = ncmpi_inq_varndims(ncid, varid[i], &ndims); \
        CHECK_ERR \
        kind = (ndims == 2) ? "fix" : "rec"; \
        err = ncmpi_inq_varoffset(ncid, varid[i], &new_var_off[i]); \
        CHECK_ERR \
        if (debug && rank == 0) \
            printf("nvars=%d: %s var %d old offset %6lld new offset %6lld\n", \
                   nvars, kind, i, old_var_off[i], new_var_off[i]); \
        old_var_off[i] = new_var_off[i]; \
    }

static
int read_back_check(int         ncid,
                    MPI_Offset *start,
                    MPI_Offset *count)
{
    char name[64];
    int i, j, k, err, nerrs=0, rank, nvars, ndims, tdim, *int_buf;
    float *flt_buf;
    MPI_Offset nelems, nrecords;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    nelems = count[1] * count[2];
    int_buf = (int*) malloc(sizeof(int) * nelems);
    flt_buf = (float*) malloc(sizeof(float) * nelems);

    err = ncmpi_inq_nvars(ncid, &nvars);
    CHECK_ERR

    err = ncmpi_inq_dimid(ncid, "time", &tdim);
    CHECK_ERR

    err = ncmpi_inq_dimlen(ncid, tdim, &nrecords);
    CHECK_ERR

    if (debug && rank == 0) printf("number of records=%lld\n",nrecords);

    for (i=0; i<nvars; i++) {
        err = ncmpi_inq_varname(ncid, i, name);
        CHECK_ERR

        err = ncmpi_inq_varndims(ncid, i, &ndims);
        CHECK_ERR

        if (ndims == 2) { /* fix-sized variable of type int */
            err = ncmpi_get_vara_int_all(ncid, i, start+1, count+1, int_buf);
            CHECK_ERR

            for (j=0; j<nelems; j++) {
                int exp = rank + j + i;
                if (int_buf[j] != exp) {
                    printf("Error at rank %d: fix var %s buf[%d] expects %d but got %d\n",
                           rank, name, j, exp, int_buf[j]);
                    return 1;
                }
            }
        }
        else if (ndims == 3) { /* record variable of type float */
            for (k=0; k<nrecords; k++) {
                start[0] = k;
                err = ncmpi_get_vara_float_all(ncid, i, start, count, flt_buf);
                CHECK_ERR

                for (j=0; j<nelems; j++) {
                    float exp = rank + j + i;
                    if (flt_buf[j] != exp) {
                        printf("Error at rank %d: rec %d var %s buf[%d] expects %.f but got %.f\n",
                               rank, k, name, j, exp, flt_buf[j]);
                        return 1;
                    }
                }
            }
        }
        else {
            printf("Error at rank %d: var %s unexpected ndims %d\n",
                   rank, name, ndims);
            return 1;
        }
    }
    return nerrs;
}

static int
tst_fmt(char *filename, int cmode)
{
    char *kind, *fname = basename(__FILE__);
    int i, err, nerrs=0, nprocs, rank, psize[2], rank_y, rank_x;
    int ncid, ndims, dimids[3], nvars, varid[NVARS], *int_buf=NULL;
    float *flt_buf=NULL;
    MPI_Offset start[3], count[3], nelems;
    MPI_Offset old_var_off[NVARS], new_var_off[NVARS];
    MPI_Offset h_size, h_extent;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (i=0; i<NVARS; i++) old_var_off[i] = new_var_off[i] = -1;

    /* calculate number of processes along each dimension */
    psize[0] = psize[1] = 0;
    MPI_Dims_create(nprocs, 2, psize);
    rank_y = rank / psize[1];
    rank_x = rank % psize[1];
    if (debug) printf("rank %d: rank_y=%d rank_x=%d\n",rank,rank_y,rank_x);

    start[0] = 0;
    start[1] = LON * rank_y;
    start[2] = LAT * rank_x;
    count[0] = 1;
    count[1] = LON;
    count[2] = LAT;
    nelems = count[1] * count[2];

    /* allocate buffers */
    int_buf = (int*) malloc(sizeof(int) * nelems);
    flt_buf = (float*) malloc(sizeof(float) * nelems);

    /* create a file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* define 3 dimensions */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimids[0]); CHECK_ERR
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "lon", LON*psize[0], &dimids[1]); CHECK_ERR
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "lat", LAT*psize[1], &dimids[2]); CHECK_ERR
    CHECK_ERR

    /* define 4 variables: 2 fix-sized and 2 record */
    nvars = 0;

    err = ncmpi_def_var(ncid, "lon", NC_INT, 2, dimids+1, &varid[0]);
    CHECK_ERR
    nvars++;

    err = ncmpi_def_var(ncid, "lat", NC_INT, 2, dimids+1, &varid[1]);
    CHECK_ERR
    nvars++;

    err = ncmpi_def_var(ncid, "wind", NC_FLOAT, 3, dimids, &varid[2]);
    CHECK_ERR
    nvars++;

    err = ncmpi_def_var(ncid, "temperature", NC_FLOAT, 3, dimids, &varid[3]);
    CHECK_ERR
    nvars++;

    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld\n", fname,__LINE__,
               h_size, h_extent);

    for (i=0; i<nvars; i++) {
        err = ncmpi_inq_varndims(ncid, varid[i], &ndims);
        CHECK_ERR
        kind = (ndims == 2) ? "fix" : "rec";
        err = ncmpi_inq_varoffset(ncid, varid[i], &old_var_off[i]);
        CHECK_ERR
        if (debug && rank == 0)
            printf("nvars=%d: %s var %d offset %6lld\n", nvars, kind, i,
                   old_var_off[i]);
    }

    /* write to file */
    for (i=0; i<nelems; i++) int_buf[i] = rank + i;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start+1, count+1, int_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) int_buf[i] = rank + i + 1;
    err = ncmpi_put_vara_int_all(ncid, varid[1], start+1, count+1, int_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 2;
    err = ncmpi_put_vara_float_all(ncid, varid[2], start, count, flt_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 3;
    err = ncmpi_put_vara_float_all(ncid, varid[3], start, count, flt_buf);
    CHECK_ERR

    /* write 2nd record */
    start[0] = 1;

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 2;
    err = ncmpi_put_vara_float_all(ncid, varid[2], start, count, flt_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 3;
    err = ncmpi_put_vara_float_all(ncid, varid[3], start, count, flt_buf);
    CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the file */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_WRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld\n", fname,__LINE__,
               h_size, h_extent);

    err = read_back_check(ncid, start, count);
    if (err > 0) {
        printf("Error at rank %d line %d: read_back_check() failed\n",rank, __LINE__);
        nerrs++;
        goto err_out;
    }

    for (i=0; i<nvars; i++) {
        err = ncmpi_inq_varndims(ncid, varid[i], &ndims);
        CHECK_ERR
        kind = (ndims == 2) ? "fix" : "rec";
        err = ncmpi_inq_varoffset(ncid, varid[i], &new_var_off[i]);
        CHECK_ERR
        if (new_var_off[i] != old_var_off[i]) {
            printf("Error: %s var %d offset expect %6lld but got %6lld\n",
                   kind, i, old_var_off[i], new_var_off[i]);
            nerrs++;
            goto err_out;
        }
        new_var_off[i] = old_var_off[i];
    }

    /* enter define mode and add a new fix-sized variable */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* add a new fix-sized variable */
    err = ncmpi_def_var(ncid, "lev", NC_INT, 2, dimids+1, &varid[4]);
    CHECK_ERR
    nvars++;

    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    if (debug && rank == 0) {
        printf("%s at %d: after adding a new fix=sized variable\n",
               fname,__LINE__);
        printf("%s at %d: header size=%lld extent=%lld\n", fname,__LINE__,
               h_size, h_extent);
    }

    for (i=0; i<nelems; i++) int_buf[i] = rank + i + 4;
    err = ncmpi_put_vara_int_all(ncid, varid[4], start+1, count+1, int_buf);
    CHECK_ERR

    PRINT_VAR_OFF

    err = read_back_check(ncid, start, count);
    if (err > 0) {
        printf("Error at rank %d line %d: read_back_check() failed\n",rank, __LINE__);
        nerrs++;
        goto err_out;
    }

    /* enter define mode and add a new record variable */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* add a new record variable */
    err = ncmpi_def_var(ncid, "snow", NC_FLOAT, 3, dimids, &varid[5]);
    CHECK_ERR
    nvars++;

    err = ncmpi_enddef(ncid); CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld\n", fname,__LINE__,
               h_size, h_extent);

    start[0] = 0;
    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 5;
    err = ncmpi_put_vara_float_all(ncid, varid[5], start, count, flt_buf);
    CHECK_ERR

    start[0] = 1;
    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 5;
    err = ncmpi_put_vara_float_all(ncid, varid[5], start, count, flt_buf);
    CHECK_ERR

    PRINT_VAR_OFF

    err = read_back_check(ncid, start, count);
    if (err > 0) {
        printf("Error at rank %d line %d: read_back_check() failed\n",rank, __LINE__);
        nerrs++;
        goto err_out;
    }

err_out:
    err = ncmpi_close(ncid); CHECK_ERR

    if (int_buf != NULL) free(int_buf);
    if (flt_buf != NULL) free(flt_buf);

    return nerrs;
}

int main(int argc, char **argv) {
    char filename[256];
    int  err, nerrs=0, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    debug = 0;

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
        sprintf(cmd_str, "*** TESTING C   %s for growing data section ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    nerrs += tst_fmt(filename, 0);
    if (nerrs > 0) goto err_out;
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
    if (nerrs > 0) goto err_out;
    nerrs += tst_fmt(filename, NC_64BIT_DATA);
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
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}
