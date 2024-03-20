/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests header extent alignment if entering the redefine mode.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_redefine tst_redefine.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./tst_redefine testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

static int verbose;

#define NUM_RNDUP(x, unit) ((((x) + (unit) - 1) / (unit)) * (unit))

#define LEN 16

#define CHECK_VAL(ncid, varid, ii, val, expect) {                     \
    if (val != expect) {                                              \
        char name[16];                                                \
        err = ncmpi_inq_varname(ncid, varid, name);                   \
        CHECK_ERROUT                                                  \
        printf("%s line %d: var %s i=%d expecting %d but got %d\n",   \
               __func__,__LINE__,name,ii,expect,val);                 \
        nerrs++;                                                      \
        goto err_out;                                                 \
    }                                                                 \
}

static int
check_fix_vars(MPI_Comm comm, int ncid, int *varid)
{
    int i, nerrs=0, err, rank, *buf;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);

    start[0] = 0; start[1] = rank * LEN;
    count[0] = 1; count[1] = LEN;

    buf = (int*) malloc(sizeof(int) * count[0] * count[1]);

    for (i=0; i<count[0]*count[1]; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[0], start+1, count+1, buf); CHECK_ERROUT
    for (i=0; i<LEN; i++)
        CHECK_VAL(ncid, varid[0], i, buf[i], rank+i)

    for (i=0; i<count[0]*count[1]; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[1], start+1, count+1, buf); CHECK_ERROUT
    for (i=0; i<LEN; i++)
        CHECK_VAL(ncid, varid[1], i, buf[i], rank+i+100)

err_out:
    free(buf);
    return nerrs;
}

static int
check_rec_vars(MPI_Comm comm, int ncid, int *varid)
{
    int i, nerrs=0, err, rank, *buf;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);

    start[0] = 0; start[1] = rank * LEN;
    count[0] = 2; count[1] = LEN;

    buf = (int*) malloc(sizeof(int) * count[0] * count[1]);

    for (i=0; i<count[0]*count[1]; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[0], start,   count,   buf); CHECK_ERROUT
    for (i=0; i<2*LEN; i++)
        CHECK_VAL(ncid, varid[0], i, buf[i], rank+i+1000)

    for (i=0; i<count[0]*count[1]; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[1], start,   count,   buf); CHECK_ERROUT
    for (i=0; i<2*LEN; i++)
        CHECK_VAL(ncid, varid[1], i, buf[i], rank+i+10000)

err_out:
    free(buf);
    return nerrs;
}

static int
tst_fmt(char *filename, int cmode)
{
    int i, rank, nprocs, ncid, err, nerrs=0;
    int *buf, dimid[3], varid[4], minfree=100;
    MPI_Offset start[2], count[2];
    MPI_Offset old_hsize, hsize;
    MPI_Offset old_extent, extent;
    MPI_Offset old_var_off, var_off;
    MPI_Offset v_align, r_align, exp_var_off;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    unsetenv("PNETCDF_HINTS");

    /* create a new file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim", LEN*nprocs, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "fa", NC_INT, 1, dimid+1, &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "fb", NC_INT, 1, dimid+1, &varid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "ta", NC_INT, 2, dimid,   &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "tb", NC_INT, 2, dimid,   &varid[3]); CHECK_ERR

    /* disable header alignment */
    err = ncmpi__enddef(ncid, 0, 4, 0, 4); CHECK_ERR

    start[0] = 0; start[1] = rank * LEN;
    count[0] = 2; count[1] = LEN;

    buf = (int*) malloc(sizeof(int) * count[0] * count[1]);

    for (i=0; i<LEN; i++) buf[i] = rank + i;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start+1, count+1, buf); CHECK_ERR
    for (i=0; i<LEN; i++) buf[i] = rank + i + 100;
    err = ncmpi_put_vara_int_all(ncid, varid[1], start+1, count+1, buf); CHECK_ERR
    for (i=0; i<2*LEN; i++) buf[i] = rank + i + 1000;
    err = ncmpi_put_vara_int_all(ncid, varid[2], start,   count,   buf); CHECK_ERR
    for (i=0; i<2*LEN; i++) buf[i] = rank + i + 10000;
    err = ncmpi_put_vara_int_all(ncid, varid[3], start,   count,   buf); CHECK_ERR

    nerrs += check_fix_vars(comm, ncid, varid);
    if (nerrs > 0) goto err_out;
    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

    err = ncmpi_close(ncid); CHECK_ERR

    /* reopen the file and check file header size and extent */
    err = ncmpi_open(comm, filename, NC_WRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "fa", &varid[0]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "fb", &varid[1]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "ta", &varid[2]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "tb", &varid[3]); CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR
    if (verbose)
        printf("File header size = %lld extent = %lld\n", hsize, extent);

    /* make sure header size == extent */
    if (hsize != extent) {
        nerrs++;
        printf("Error at line %d in %s: File header size %lld != extent %lld\n",
               __LINE__,__FILE__, hsize, extent);
    }

    /* enter redefine mode and add a dimension to grow the header */
    err = ncmpi_redef(ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "dim_x", 25, &dimid[2]); CHECK_ERR

    /* disable header alignment */
    err = ncmpi__enddef(ncid, 0, 4, 0, 4); CHECK_ERR

    old_hsize = hsize;
    old_extent = extent;
    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR
    if (verbose) {
        printf("File old header size = %6lld old extent = %6lld\n", old_hsize, old_extent);
        printf("File new header size = %6lld new extent = %6lld\n", hsize, extent);
    }

    /* make sure header size grows */
    if (hsize <= old_hsize) {
        nerrs++;
        printf("Error at line %d in %s: File header size %lld fails to grow from %lld\n",
               __LINE__,__FILE__, hsize, old_hsize);
    }
    /* make sure header extent grows */
    if (extent <= old_extent) {
        nerrs++;
        printf("Error at line %d in %s: File header extent %lld fails to grow from %lld\n",
               __LINE__,__FILE__, extent, old_extent);
    }

    nerrs += check_fix_vars(comm, ncid, varid);
    if (nerrs > 0) goto err_out;
    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

    /* enter redefine mode and add nothing */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* add free space into header extent */
    err = ncmpi__enddef(ncid, minfree, 0, 0, 0); CHECK_ERR

    old_hsize = hsize;
    old_extent = extent;
    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR
    if (verbose) {
        printf("File old header size = %6lld old extent = %6lld\n", old_hsize, old_extent);
        printf("File new header size = %6lld new extent = %6lld\n", hsize, extent);
    }

    /* make sure header size does not change */
    if (hsize != old_hsize) {
        nerrs++;
        printf("Error at line %d in %s: File header size %lld changed from %lld\n",
               __LINE__,__FILE__, hsize, old_hsize);
    }
    /* make sure header extent grows minfree bytes */
    if (extent != old_extent + minfree) {
        nerrs++;
        printf("Error at line %d in %s: File header extent %lld fails to grow into %lld\n",
               __LINE__,__FILE__, extent, old_extent + minfree);
    }

    nerrs += check_fix_vars(comm, ncid, varid);
    if (nerrs > 0) goto err_out;
    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

    /* enter redefine mode and add nothing */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* align header extent to 512 bytes */
    err = ncmpi__enddef(ncid, 0, 512, 0, 0); CHECK_ERR

    old_hsize = hsize;
    old_extent = extent;
    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR
    if (verbose) {
        printf("File old header size = %6lld old extent = %6lld\n", old_hsize, old_extent);
        printf("File new header size = %6lld new extent = %6lld\n", hsize, extent);
    }

    /* make sure header size does not change */
    if (hsize != old_hsize) {
        nerrs++;
        printf("Error at line %d in %s: File header size %lld changed from %lld\n",
               __LINE__,__FILE__, hsize, old_hsize);
    }
    /* make sure header extent grows to 512 bytes */
    if (extent != 512) {
        nerrs++;
        printf("Error at line %d in %s: File header extent %lld fails to grow into 512\n",
               __LINE__,__FILE__, extent);
    }

    nerrs += check_fix_vars(comm, ncid, varid);
    if (nerrs > 0) goto err_out;
    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

    err = ncmpi_close(ncid); CHECK_ERR


    setenv("PNETCDF_HINTS", "nc_header_align_size=700", 1);

    /* reopen the file and check file header size and extent */
    err = ncmpi_open(comm, filename, NC_WRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "fa", &varid[0]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "fb", &varid[1]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "ta", &varid[2]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "tb", &varid[3]); CHECK_ERR

    /* enter redefine mode and add nothing */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* align header extent to 800 bytes, this should take no effect, as
     * run-time hints precedes function calls. */
    err = ncmpi__enddef(ncid, 0, 800, 0, 0); CHECK_ERR

    old_hsize = hsize;
    old_extent = extent;
    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR
    if (verbose) {
        printf("File old header size = %6lld old extent = %6lld\n", old_hsize, old_extent);
        printf("File new header size = %6lld new extent = %6lld\n", hsize, extent);
    }

    /* make sure header size does not change */
    if (hsize != old_hsize) {
        nerrs++;
        printf("Error at line %d in %s: File header size %lld changed from %lld\n",
               __LINE__,__FILE__, hsize, old_hsize);
    }
    /* make sure header extent grows to 700 bytes */
    if (extent != 700) {
        nerrs++;
        printf("Error at line %d in %s: File header extent %lld fails to grow into 700\n",
               __LINE__,__FILE__, extent);
    }

    unsetenv("PNETCDF_HINTS");

    nerrs += check_fix_vars(comm, ncid, varid);
    if (nerrs > 0) goto err_out;
    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

    /* enter redefine mode and add nothing */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* set minfree to the free space, this should not grow the extent */
    err = ncmpi__enddef(ncid, extent-hsize, 0, 0, 0); CHECK_ERR

    old_hsize = hsize;
    old_extent = extent;
    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR
    if (verbose) {
        printf("File old header size = %6lld old extent = %6lld\n", old_hsize, old_extent);
        printf("File new header size = %6lld new extent = %6lld\n", hsize, extent);
    }

    /* make sure header size does not change */
    if (hsize != old_hsize) {
        nerrs++;
        printf("Error at line %d in %s: File header size %lld changed from %lld\n",
               __LINE__,__FILE__, hsize, old_hsize);
    }
    /* make sure header extent remain the same */
    if (extent != old_extent) {
        nerrs++;
        printf("Error at line %d in %s: File header extent %lld changed from %lld\n",
               __LINE__,__FILE__, extent, old_extent);
    }

    /* obtain 1st record variable's file offset */
    err = ncmpi_inq_varoffset(ncid, varid[2], &old_var_off); CHECK_ERR

    nerrs += check_fix_vars(comm, ncid, varid);
    if (nerrs > 0) goto err_out;
    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

    /* enter redefine mode and add nothing */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* set v_minfree of size 400 */
    err = ncmpi__enddef(ncid, 0, 0, 400, 0); CHECK_ERR

    /* obtain 1st record variable's file offset */
    err = ncmpi_inq_varoffset(ncid, varid[2], &var_off); CHECK_ERR

    /* var_off should grows 400 bytes */
    if (var_off != old_var_off + 400) {
        nerrs++;
        printf("Error at line %d in %s: 1st record variable offset %lld (expecting %lld)\n",
               __LINE__,__FILE__, var_off, old_var_off+400);
    }

    nerrs += check_fix_vars(comm, ncid, varid);
    if (nerrs > 0) goto err_out;
    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

#if 0
    err = ncmpi_close(ncid); CHECK_ERR

    /* reopen the file and set r_align */
    err = ncmpi_open(comm, filename, NC_WRITE, MPI_INFO_NULL, &ncid); CHECK_ERR
#endif

    /* obtained the old offset of 1st record variable */
    err = ncmpi_inq_varoffset(ncid, varid[2], &old_var_off); CHECK_ERR

    /* enter redefine mode and add nothing */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* set r_align to 1500 */
    r_align = 1500;
    err = ncmpi__enddef(ncid, 0, 0, 0, r_align); CHECK_ERR

    /* obtain 1st record variable's file offset */
    err = ncmpi_inq_varoffset(ncid, varid[2], &var_off); CHECK_ERR

    /* round up to r_align */
    exp_var_off = NUM_RNDUP(old_var_off, r_align);

    /* var_off should grows into 1500 bytes */
    if (var_off != exp_var_off) {
        nerrs++;
        printf("Error at line %d in %s: 1st record variable offset %lld (expecting %lld)\n",
               __LINE__,__FILE__, var_off, exp_var_off);
    }

    nerrs += check_fix_vars(comm, ncid, varid);
    if (nerrs > 0) goto err_out;
    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

    err = ncmpi_close(ncid); CHECK_ERR

    unsetenv("PNETCDF_HINTS");

    /* re-create a new file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* define only record variables */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim", LEN*nprocs, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "ta", NC_INT, 2, dimid,   &varid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "tb", NC_INT, 2, dimid,   &varid[3]); CHECK_ERR

    /* set v_align and r_align. Without record variables, r_align will take
     * effect over v_align
     */
    v_align = 100;
    r_align = 512;
    err = ncmpi__enddef(ncid, 0, v_align, 0, r_align); CHECK_ERR

    start[0] = 0; start[1] = rank * LEN;
    count[0] = 2; count[1] = LEN;

    for (i=0; i<2*LEN; i++) buf[i] = rank + i + 1000;
    err = ncmpi_put_vara_int_all(ncid, varid[2], start,   count,   buf); CHECK_ERR
    for (i=0; i<2*LEN; i++) buf[i] = rank + i + 10000;
    err = ncmpi_put_vara_int_all(ncid, varid[3], start,   count,   buf); CHECK_ERR

    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR
    if (verbose)
        printf("File header size = %lld extent = %lld\n", hsize, extent);

    /* obtain 1st record variable's file offset */
    err = ncmpi_inq_varoffset(ncid, varid[2], &var_off); CHECK_ERR
    if (var_off != extent) {
        nerrs++;
        printf("Error at line %d in %s: 1st record variable offset %lld (expecting %lld)\n",
               __LINE__,__FILE__, var_off, extent);
    }

    /* extent should be r_align */
    if (extent != r_align) {
        nerrs++;
        printf("Error at line %d in %s: file extent %lld (expecting %lld)\n",
               __LINE__,__FILE__, extent, r_align);
    }

    nerrs += check_rec_vars(comm, ncid, varid+2);
    if (nerrs > 0) goto err_out;

err_out:
    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);

    return nerrs;
}

int main(int argc, char** argv)
{
    char filename[256];
    int commsize, rank, err, nerrs=0;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &commsize);
    MPI_Comm_rank(comm, &rank);

    verbose = 0;

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "tst_redefine.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for header alignment ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    nerrs += tst_fmt(filename, 0);
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
    nerrs += tst_fmt(filename, NC_64BIT_DATA);

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
