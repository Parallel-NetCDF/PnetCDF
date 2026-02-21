/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This program checks whether the data movement subroutines are called in the
 * following scenario.
 * 1. There is a sufficiently large free space in the header extent and fix-sized
 *    variable section, such that all successive adding new variables will not
 *    change the starting offsets of both fix-sized and record variable sections.
 * 2. Add new fix-sized variables should not call the data movement subroutines.
 * 3. When there is only one record, i.e. unlimited dimension size == 1, adding
 *    a new record variable should not call the data movement subroutines.
 * 4. When there are two records, i.e. unlimited dimension size == 2, adding a
 *    new a new record variable SHOULD call the data movement subroutines, but
 *    only to move the 2nd record to a higher offset, and by moving the entire
 *    record, not individual record variables separately.
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
        if (old_var_off[i] > 0 && old_var_off[i] != new_var_off[i]) { \
            printf("Error: %s var %d old offset %6lld != new offset %6lld\n", \
                   kind, i, old_var_off[i], new_var_off[i]); \
            nerrs++; \
            goto err_out; \
        } \
        if (debug && rank == 0) \
            printf("nvars=%d: %s var %d old offset %6lld new offset %6lld\n", \
                   nvars, kind, i, old_var_off[i], new_var_off[i]); \
        old_var_off[i] = new_var_off[i]; \
    }


static
int read_back_check(int         ncid,
                    int         coll_io,
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
            if (coll_io)
                err = ncmpi_get_vara_int_all(ncid, i, start+1, count+1, int_buf);
            else
                err = ncmpi_get_vara_int(ncid, i, start+1, count+1, int_buf);
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
                if (coll_io)
                    err = ncmpi_get_vara_float_all(ncid, i, start, count, flt_buf);
                else
                    err = ncmpi_get_vara_float(ncid, i, start, count, flt_buf);
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
    free(int_buf);
    free(flt_buf);

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char *kind, *fname = basename(__FILE__);
    int i, err, nerrs=0, nprocs, rank, psize[2], rank_y, rank_x;
    int ncid, ndims, dimids[3], nvars, varid[NVARS], *int_buf=NULL;
    int varid_1st_fix, varid_1st_rec, fix_nvars=0, rec_nvars=0;
    float *flt_buf=NULL;
    MPI_Offset start[3], count[3], nelems;
    MPI_Offset old_var_off[NVARS], new_var_off[NVARS];
    MPI_Offset h_size, h_extent, var_size, fix_vars_size, rec_size;
    MPI_Offset h_minfree, v_align, v_minfree, r_align;

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
    var_size = sizeof(int) * nelems * nprocs;

    /* allocate buffers */
    int_buf = (int*) malloc(sizeof(int) * nelems);
    flt_buf = (float*) malloc(sizeof(float) * nelems);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a file */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define 3 dimensions */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimids[0]);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "lon", LON*psize[0], &dimids[1]);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "lat", LAT*psize[1], &dimids[2]);
    CHECK_ERR

    /* define 4 variables: 2 fix-sized and 2 record */
    nvars = 0;

    err = ncmpi_def_var(ncid, "lon", NC_INT, 2, dimids+1, &varid[0]);
    CHECK_ERR
    varid_1st_fix = varid[0];
    nvars++;
    fix_nvars++;

    err = ncmpi_def_var(ncid, "lat", NC_INT, 2, dimids+1, &varid[1]);
    CHECK_ERR
    nvars++;
    fix_nvars++;

    err = ncmpi_def_var(ncid, "wind", NC_FLOAT, 3, dimids, &varid[2]);
    CHECK_ERR
    varid_1st_rec = varid[2];
    nvars++;
    rec_nvars++;

    err = ncmpi_def_var(ncid, "temperature", NC_FLOAT, 3, dimids, &varid[3]);
    CHECK_ERR
    nvars++;
    rec_nvars++;

    err = ncmpi_inq_header_size(ncid, &h_size);
    CHECK_ERR

    /* allocate a sufficient large space for header extent and fix-sized
     * variable section
     */
    h_minfree = h_size;   /* a free space enough for adding 4 more variables */
    v_align = 1024;
    v_minfree = var_size; /* a free space enough for a new fix-size variable */
    r_align = 1024;

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    if (debug && rank == 0) {
        printf("\n%s at %d: Create file\n", fname,__LINE__);
        printf("\t\t\th_minfree=%lld v_align=%lld v_minfree=%lld r_align=%lld h_size=%lld\n",
               h_minfree,v_align,v_minfree,r_align,h_size);
    }

    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    err = ncmpi_inq_recsize(ncid, &rec_size); CHECK_ERR
    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld var_size=%lld rec_size=%lld\n",
               fname,__LINE__, h_size, h_extent, var_size, rec_size);

    fix_vars_size = 0;
    for (i=0; i<nvars; i++) {
        err = ncmpi_inq_varndims(ncid, varid[i], &ndims);
        CHECK_ERR
        kind = (ndims == 2) ? "fix" : "rec";
        err = ncmpi_inq_varoffset(ncid, varid[i], &old_var_off[i]);
        CHECK_ERR
        if (ndims == 2) fix_vars_size += var_size;
        if (debug && rank == 0)
            printf("nvars=%d: %s var %d offset %6lld\n", nvars, kind, i,
                   old_var_off[i]);
    }
    if (debug && rank == 0) {
        MPI_Offset fix_off, rec_off;
        err = ncmpi_inq_varoffset(ncid, varid_1st_fix, &fix_off);
        CHECK_ERR
        printf("Line %d: fix_nvars=%d rec_nvars=%d \n", __LINE__,fix_nvars,rec_nvars);
        printf("Line %d: first fix-sized var offset %6lld\n", __LINE__,fix_off);
        err = ncmpi_inq_varoffset(ncid, varid_1st_rec, &rec_off);
        CHECK_ERR
        printf("Line %d: sum of all fix-sized vars  %6lld\n", __LINE__,fix_vars_size);
        printf("Line %d: fix-sized var free space   %6lld\n", __LINE__,
                rec_off - (fix_off + fix_vars_size));
        printf("Line %d: 1st record offset          %6lld\n", __LINE__,rec_off);
    }

    /* write to file */
    for (i=0; i<nelems; i++) int_buf[i] = rank + i;
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid[0], start+1, count+1, int_buf);
    else
        err = ncmpi_put_vara_int(ncid, varid[0], start+1, count+1, int_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) int_buf[i] = rank + i + 1;
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid[1], start+1, count+1, int_buf);
    else
        err = ncmpi_put_vara_int(ncid, varid[1], start+1, count+1, int_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 2;
    if (coll_io)
        err = ncmpi_put_vara_float_all(ncid, varid[2], start, count, flt_buf);
    else
        err = ncmpi_put_vara_float(ncid, varid[2], start, count, flt_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 3;
    if (coll_io)
        err = ncmpi_put_vara_float_all(ncid, varid[3], start, count, flt_buf);
    else
        err = ncmpi_put_vara_float(ncid, varid[3], start, count, flt_buf);
    CHECK_ERR

    /* file sync before reading */
    if (!coll_io) {
        err = ncmpi_sync(ncid);
        CHECK_ERR
    }
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* open the file */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_WRITE, info, &ncid);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    if (debug && rank == 0) {
        printf("\n%s at %d: Re-open the file\n", fname,__LINE__);
        printf("\t\t\theader size=%lld extent=%lld\n", h_size, h_extent);
    }

    err = read_back_check(ncid, coll_io, start, count);
    if (err > 0) {
        printf("Error at rank %d line %d: read_back_check() failed\n",rank, __LINE__);
        nerrs++;
        goto err_out;
    }

    fix_vars_size = 0;
    for (i=0; i<nvars; i++) {
        err = ncmpi_inq_varndims(ncid, varid[i], &ndims);
        CHECK_ERR
        kind = (ndims == 2) ? "fix" : "rec";
        if (ndims == 2) fix_vars_size += var_size;
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

    if (debug && rank == 0) {
        MPI_Offset fix_off, rec_off;
        err = ncmpi_inq_varoffset(ncid, varid_1st_fix, &fix_off);
        CHECK_ERR
        printf("Line %d: fix_nvars=%d rec_nvars=%d \n", __LINE__,fix_nvars,rec_nvars);
        printf("Line %d: first fix-sized var offset %6lld\n", __LINE__,fix_off);
        err = ncmpi_inq_varoffset(ncid, varid_1st_rec, &rec_off);
        CHECK_ERR
        printf("Line %d: sum of all fix-sized vars  %6lld\n", __LINE__,fix_vars_size);
        printf("Line %d: fix-sized var free space   %6lld\n", __LINE__,
                rec_off - (fix_off + fix_vars_size));
        printf("Line %d: 1st record offset          %6lld\n", __LINE__,rec_off);
    }

    /* enter define mode and add a new fix-sized variable */
    err = ncmpi_redef(ncid); CHECK_ERR

    if (debug && rank == 0) {
        printf("\n===================================================\n");
        printf("%s at %d: fix_nvars=%d adding a new fix-sized variable\n",
               fname,__LINE__, fix_nvars);
    }

    /* add a new fix-sized variable */
    err = ncmpi_def_var(ncid, "lev", NC_INT, 2, dimids+1, &varid[4]);
    CHECK_ERR
    nvars++;
    fix_nvars++;

    fix_vars_size += var_size;

    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    if (debug && rank == 0) {
        printf("%s at %d: after adding a new fix-sized variable\n",
               fname,__LINE__);
        printf("\t\t\theader size=%lld extent=%lld\n", h_size, h_extent);
    }

    for (i=0; i<nelems; i++) int_buf[i] = rank + i + 4;
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid[4], start+1, count+1, int_buf);
    else
        err = ncmpi_put_vara_int(ncid, varid[4], start+1, count+1, int_buf);
    CHECK_ERR

    PRINT_VAR_OFF

    err = read_back_check(ncid, coll_io, start, count);
    if (err > 0) {
        printf("Error at rank %d line %d: read_back_check() failed\n",rank, __LINE__);
        nerrs++;
        goto err_out;
    }

    if (debug && rank == 0) {
        MPI_Offset fix_off, rec_off;
        err = ncmpi_inq_varoffset(ncid, varid_1st_fix, &fix_off);
        CHECK_ERR
        printf("Line %d: fix_nvars=%d rec_nvars=%d \n", __LINE__,fix_nvars,rec_nvars);
        printf("Line %d: first fix-sized var offset %6lld\n", __LINE__,fix_off);
        err = ncmpi_inq_varoffset(ncid, varid_1st_rec, &rec_off);
        CHECK_ERR
        printf("Line %d: sum of all fix-sized vars  %6lld\n", __LINE__,fix_vars_size);
        printf("Line %d: fix-sized var free space   %6lld\n", __LINE__,
                rec_off - (fix_off + fix_vars_size));
        printf("Line %d: 1st record offset          %6lld\n", __LINE__,rec_off);
    }

    /* enter define mode and add a new record variable */
    err = ncmpi_redef(ncid); CHECK_ERR

    if (debug && rank == 0) {
        printf("\n===================================================\n");
        printf("%s at %d: rec_nvars=%d adding a new record variable\n",
                fname,__LINE__,rec_nvars);
    }

    /* add a new record variable */
    err = ncmpi_def_var(ncid, "snow", NC_FLOAT, 3, dimids, &varid[5]);
    CHECK_ERR
    nvars++;
    rec_nvars++;

    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    err = ncmpi_inq_recsize(ncid, &rec_size); CHECK_ERR
    if (debug && rank == 0) {
        printf("%s at %d: after adding a new record variable\n",
               fname,__LINE__);
        printf("\t\t\theader size=%lld extent=%lld rec_size=%lld\n",
               h_size, h_extent,rec_size);
    }

    start[0] = 0;
    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 5;
    if (coll_io)
        err = ncmpi_put_vara_float_all(ncid, varid[5], start, count, flt_buf);
    else
        err = ncmpi_put_vara_float(ncid, varid[5], start, count, flt_buf);
    CHECK_ERR

    PRINT_VAR_OFF

    err = read_back_check(ncid, coll_io, start, count);
    if (err > 0) {
        printf("Error at rank %d line %d: read_back_check() failed\n",rank, __LINE__);
        nerrs++;
        goto err_out;
    }

    if (debug && rank == 0) {
        MPI_Offset fix_off, rec_off;
        err = ncmpi_inq_varoffset(ncid, varid_1st_fix, &fix_off);
        CHECK_ERR
        printf("Line %d: fix_nvars=%d rec_nvars=%d \n", __LINE__,fix_nvars,rec_nvars);
        printf("Line %d: first fix-sized var offset %6lld\n", __LINE__,fix_off);
        err = ncmpi_inq_varoffset(ncid, varid_1st_rec, &rec_off);
        CHECK_ERR
        printf("Line %d: sum of all fix-sized vars  %6lld\n", __LINE__,fix_vars_size);
        printf("Line %d: fix-sized var free space   %6lld\n", __LINE__,
                rec_off - (fix_off + fix_vars_size));
        printf("Line %d: 1st record offset          %6lld\n", __LINE__,rec_off);
        printf("Line %d: 2nd record offset          %6lld\n", __LINE__,rec_off+rec_size);
    }

    /* write 2nd record */
    if (debug && rank == 0)
        printf("\n%s at %d: increment to 2 records\n", fname,__LINE__);

    start[0] = 1;

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 2;
    if (coll_io)
        err = ncmpi_put_vara_float_all(ncid, varid[2], start, count, flt_buf);
    else
        err = ncmpi_put_vara_float(ncid, varid[2], start, count, flt_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 3;
    if (coll_io)
        err = ncmpi_put_vara_float_all(ncid, varid[3], start, count, flt_buf);
    else
        err = ncmpi_put_vara_float(ncid, varid[3], start, count, flt_buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 5;
    if (coll_io)
        err = ncmpi_put_vara_float_all(ncid, varid[5], start, count, flt_buf);
    else
        err = ncmpi_put_vara_float(ncid, varid[5], start, count, flt_buf);
    CHECK_ERR

    /* enter define mode and add a new record variable */
    err = ncmpi_redef(ncid); CHECK_ERR

    if (debug && rank == 0) {
        printf("\n===================================================\n");
        printf("%s at %d: rec_nvars=%d adding a new record variable\n",
               fname,__LINE__,rec_nvars);
    }

    /* add a new record variable */
    err = ncmpi_def_var(ncid, "ice", NC_FLOAT, 3, dimids, &varid[6]);
    CHECK_ERR
    nvars++;
    rec_nvars++;

    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_header_size(ncid, &h_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &h_extent); CHECK_ERR
    err = ncmpi_inq_recsize(ncid, &rec_size); CHECK_ERR
    if (debug && rank == 0)
        printf("%s at %d: header size=%lld extent=%lld rec_size=%lld\n",
               fname,__LINE__, h_size, h_extent,rec_size);

    /* write 1st record */
    start[0] = 0;
    for (i=0; i<nelems; i++) flt_buf[i] = rank + i + 6;
    if (coll_io)
        err = ncmpi_put_vara_float_all(ncid, varid[6], start, count, flt_buf);
    else
        err = ncmpi_put_vara_float(ncid, varid[6], start, count, flt_buf);
    CHECK_ERR
    /* write 2nd record */
    start[0] = 1;
    if (coll_io)
        err = ncmpi_put_vara_float_all(ncid, varid[6], start, count, flt_buf);
    else
        err = ncmpi_put_vara_float(ncid, varid[6], start, count, flt_buf);
    CHECK_ERR

    PRINT_VAR_OFF

    err = read_back_check(ncid, coll_io, start, count);
    if (err > 0) {
        printf("Error at rank %d line %d: read_back_check() failed\n",rank, __LINE__);
        nerrs++;
        goto err_out;
    }

    if (debug && rank == 0) {
        MPI_Offset fix_off, rec_off;
        err = ncmpi_inq_varoffset(ncid, varid_1st_fix, &fix_off);
        CHECK_ERR
        printf("Line %d: fix_nvars=%d rec_nvars=%d \n", __LINE__,fix_nvars,rec_nvars);
        printf("Line %d: first fix-sized var offset %6lld\n", __LINE__,fix_off);
        err = ncmpi_inq_varoffset(ncid, varid_1st_rec, &rec_off);
        CHECK_ERR
        printf("Line %d: sum of all fix-sized vars  %6lld\n", __LINE__,fix_vars_size);
        printf("Line %d: fix-sized var free space   %6lld\n", __LINE__,
                rec_off - (fix_off + fix_vars_size));
        printf("Line %d: 1st record offset          %6lld\n", __LINE__,rec_off);
        printf("Line %d: 2nd record offset          %6lld\n", __LINE__,rec_off+rec_size);
    }

err_out:
    err = ncmpi_close(ncid); CHECK_ERR

    if (int_buf != NULL) free(int_buf);
    if (flt_buf != NULL) free(flt_buf);
    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1;    /* test intra-node aggregation */
    opt.drv      = 1;    /* test PNCIO driver */
    opt.ind      = 1;    /* test hint romio_no_indep_rw */
    opt.chk      = 5432; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1;    /* test burst-buffering feature */
    opt.mod      = 1;    /* test independent data mode */
    opt.hdr_diff = 1;    /* run ncmpidiff for file header only */
    opt.var_diff = 1;    /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "growing data section", opt, test_io);

    MPI_Finalize();

    return err;
}
