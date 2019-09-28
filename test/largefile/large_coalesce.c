/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Given two large nonblocking requests and they are actually contiguous in
 * fileview and buffer type, this program tests whether filetype and buftype
 * coalescing is skipped when calling MPI_Type_create_hindexed internally.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

#define FOUR_G 4294967296LL
#define TWO_G  2147483648LL
#define ONE_G  1073741824LL

int main(int argc, char** argv)
{
    char filename[256];
    unsigned char *buf;
    int rank, nprocs, err, nerrs=0;
    int ncid, cmode, varid, dimid[2], req[3], st[3];
    MPI_Offset start[2], count[2];
    MPI_Info info;
    size_t i;
    int bb_enabled=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
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
        sprintf(cmd_str, "*** TESTING C   %s for skip filetype buftype coalesce ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    buf = (unsigned char*) calloc(TWO_G+1024,1);
    if (buf == NULL) {
        printf("malloc failed for size %lld\n", TWO_G+1024);
        MPI_Finalize();
        return 1;
    }

    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_cb_write", "enable");
    MPI_Info_set(info, "romio_ds_read", "disable"); /* run slow without it */

#if defined(ENABLE_LARGE_SINGLE_REQ) || defined(ENABLE_BURST_BUFFER)
#else
    /* silence iternal debug messages */
    setenv("PNETCDF_SAFE_MODE", "0", 1);
#endif

#ifdef ENABLE_NETCDF4
    /* Test for NetCDF 4 first as ncvalidator checks only read classic files */

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_NETCDF4;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "NPROCS", nprocs, &dimid[0]);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "X", TWO_G+1024, &dimid[1]);
    CHECK_ERR

    /* define a big 1D variable of ubyte type */
    err = ncmpi_def_var(ncid, "big_var", NC_UBYTE, 2, dimid, &varid);
    CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    /* now we are in data mode */
    for (i=0; i<20; i++) buf[ONE_G-10+i] = 'a'+i;
    for (i=0; i<20; i++) buf[TWO_G-10+i] = 'A'+i;

    start[0] = rank;
    count[0] = 1;

    start[1] = 0;
    count[1] = 10;
    err = ncmpi_put_vara_uchar_all(ncid, varid, start, count, buf);
    CHECK_ERR

    /* 2nd request is not contiguous from the first */
    start[1] = 1024;
    count[1] = ONE_G-1024;
    err = ncmpi_put_vara_uchar_all(ncid, varid, start, count, buf+1024);
    CHECK_ERR

    /* make file access and write buffer of 3rd request contiguous from the 2nd
     * request to check whether the internal fileview and buftype coalescing
     * are skipped */
    start[1] = ONE_G;
    count[1] = ONE_G+1024;
    err = ncmpi_put_vara_uchar_all(ncid, varid, start, count, buf+ONE_G);
    CHECK_ERR

    start[1] = 0;
    count[1] = 10;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf);
    CHECK_ERR

    start[1] = 1024;
    count[1] = ONE_G-1024;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf+1024);
    CHECK_ERR

    start[1] = ONE_G;
    count[1] = ONE_G+1024;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf+ONE_G);
    CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if open to read header fine */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
#endif
    /* Test classic format */

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR
    MPI_Info_free(&info);

    {
        int flag;
        char hint[MPI_MAX_INFO_VAL];
        MPI_Info infoused;

        ncmpi_inq_file_info(ncid, &infoused);
        MPI_Info_get(infoused, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, hint, &flag);
        if (flag && strcasecmp(hint, "enable") == 0)
            bb_enabled = 1;
        MPI_Info_free(&infoused);
    }

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "NPROCS", nprocs, &dimid[0]);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "X", TWO_G+1024, &dimid[1]);
    CHECK_ERR

    /* define a big 1D variable of ubyte type */
    err = ncmpi_def_var(ncid, "big_var", NC_UBYTE, 2, dimid, &varid);
    CHECK_ERR

    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    /* now we are in data mode */
#if defined(ENABLE_LARGE_SINGLE_REQ) || defined(ENABLE_BURST_BUFFER)
#ifndef ENABLE_LARGE_SINGLE_REQ
    if (bb_enabled) {
#endif
    for (i=0; i<20; i++) buf[ONE_G-10+i] = 'a'+i;
    for (i=0; i<20; i++) buf[TWO_G-10+i] = 'A'+i;
#ifndef ENABLE_LARGE_SINGLE_REQ
    }
#endif
#endif

    start[0] = rank;
    count[0] = 1;

    start[1] = 0;
    count[1] = 10;
    err = ncmpi_iput_vara_uchar(ncid, varid, start, count, buf, &req[0]);
    CHECK_ERR

    /* 2nd request is not contiguous from the first */
    start[1] = 1024;
    count[1] = ONE_G-1024;
    err = ncmpi_iput_vara_uchar(ncid, varid, start, count, buf+1024, &req[1]);
    CHECK_ERR

    /* make file access and write buffer of 3rd request contiguous from the 2nd
     * request to check whether the internal fileview and buftype coalescing
     * are skipped */
    start[1] = ONE_G;
    count[1] = ONE_G+1024;
    err = ncmpi_iput_vara_uchar(ncid, varid, start, count, buf+ONE_G, &req[2]);
    CHECK_ERR

    err = ncmpi_wait_all(ncid, 3, req, st);
#if defined(ENABLE_LARGE_SINGLE_REQ) || defined(ENABLE_BURST_BUFFER)
#ifndef ENABLE_LARGE_SINGLE_REQ
    if (bb_enabled) {
#endif
    CHECK_ERR

    /* read back to check contents */
    start[1] = ONE_G-10;
    count[1] = 20;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf);
    CHECK_ERR
    for (i=0; i<20; i++) {
        if (buf[i] != 'a'+i) {
            printf("%d (at line %d): expect buf[%lld]=%zd but got %d\n",
                   rank, __LINE__, ONE_G-10+i, i+'a', buf[i]);
            nerrs++;
        }
    }

    /* read back to check contents */
    start[1] = TWO_G-10;
    count[1] = 20;
    err = ncmpi_get_vara_uchar_all(ncid, varid, start, count, buf);
    CHECK_ERR
    for (i=0; i<20; i++) {
        if (buf[i] != 'A'+i) {
            printf("%d (at line %d): expect buf[%lld]=%zd but got %d\n",
                   rank, __LINE__, TWO_G-10+i, i+'A', buf[i]);
            nerrs++;
        }
    }

    /* test the same pattern but for iget */
    for (i=0; i<TWO_G+1024; i++) buf[i] = 0;
#ifndef ENABLE_LARGE_SINGLE_REQ
    }
    else
        EXP_ERR(NC_EMAX_REQ)
#endif
#else
    EXP_ERR(NC_EMAX_REQ)
#endif

    start[1] = 0;
    count[1] = 10;
    err = ncmpi_iget_vara_uchar(ncid, varid, start, count, buf, &req[0]);
    CHECK_ERR

    start[1] = 1024;
    count[1] = ONE_G-1024;
    err = ncmpi_iget_vara_uchar(ncid, varid, start, count, buf+1024, &req[1]);
    CHECK_ERR

    start[1] = ONE_G;
    count[1] = ONE_G+1024;
    err = ncmpi_iget_vara_uchar(ncid, varid, start, count, buf+ONE_G, &req[2]);
    CHECK_ERR

    err = ncmpi_wait_all(ncid, 3, req, st);
#ifndef ENABLE_LARGE_SINGLE_REQ
    EXP_ERR(NC_EMAX_REQ)
#else
    CHECK_ERR

    for (i=0; i<20; i++) {
        if (buf[ONE_G-10+i] != 'a'+i) {
            printf("%d (at line %d): expect buf[%lld]=%zd but got %d\n",
                   rank, __LINE__, ONE_G-10+i, i+'a', buf[i]);
            nerrs++;
        }
    }

    for (i=0; i<20; i++) {
        if (buf[TWO_G-10+i] != 'A'+i) {
            printf("%d (at line %d): expect buf[%lld]=%zd but got %d\n",
                   rank, __LINE__, TWO_G-10+i, i+'A', buf[i]);
            nerrs++;
        }
    }
#endif

    err = ncmpi_close(ncid); CHECK_ERR

    /* check if open for reading header */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    free(buf);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

