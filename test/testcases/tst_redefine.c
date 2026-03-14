/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests all alignment features available from PnetCDF:
 * 1. h_minfree: free space in the header section, i.e.
 *      (header extent) - (header size) >= h_minfree
 * 2. v_align: alignment of the beginning of the fix-size variable section, i.e.
 *       (header extent) % v_align == 0
 *       If no fixed-size variable is defined, v_align is ignored.
 *       Default value of v_align is 512.
 * 3. v_minfree: free space between the end of last fix-sized variable and the
 *       record variable section.
 *       If no fixed-size variable is defined, v_minfree is ignored.
 * 4. r_align: alignment of the beginning of the record variable section.
 *       If no fixed-size variable is defined, default value of r_align is 512.
 *       Otherwise, default value of r_align is 4.
 *
 * Tests are done by reentering the define mode multiple times.
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
#include <assert.h>

#include <pnetcdf.h>

#include <testutils.h>

#define LEN 101

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif
#define RNDUP(x, unit) ((((x) + (unit) - 1) / (unit)) * (unit))

static int verbose;

#define CHECK_VAL(ncid, varid, ii, val, expect) {                           \
    if (val != expect) {                                                    \
        char name[16];                                                      \
        int ndims;                                                          \
        err = ncmpi_inq_varndims(ncid, varid, &ndims);                      \
        CHECK_ERROUT                                                        \
        err = ncmpi_inq_varname(ncid, varid, name);                         \
        CHECK_ERROUT                                                        \
        if (ndims == 1)                                                     \
            printf("%s line %d: var %s[%d] expecting %d but got %d\n",      \
                   __func__,__LINE__,name,ii,expect,val);                   \
        else /* record variable */                                          \
            printf("%s line %d: var %s[%d][%d] expecting %d but got %d\n",  \
                   __func__,__LINE__,name,ii/LEN,ii%LEN,expect,val);        \
        nerrs++;                                                            \
        goto err_out;                                                       \
    }                                                                       \
}

/*----< check_vars() >-------------------------------------------------------*/
/* read back variables from file and check their contents */
static int
check_vars(MPI_Comm comm, int ncid, int *varid, int coll_io)
{
    int i, j, nerrs=0, err, rank, *buf[4], nvars, off_val;
    MPI_Offset start[2], count[2], bufLen[4];

    MPI_Comm_rank(comm, &rank);

    err = ncmpi_inq_nvars(ncid, &nvars); CHECK_ERROUT

    /* check fix-sized variables */
    start[0] = 0;
    start[1] = rank * LEN;
    count[1] = LEN;

    buf[0] = (int*) malloc(sizeof(int) * 2 * LEN * 4);
    for (i=0; i<2*LEN*4; i++) buf[0][i] = -1;

    /* check record variables */
    count[0] = 2;

    bufLen[0] = count[0]*count[1];

#ifdef USE_BLOCKING_APIS
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid[0], start, count, buf[0]);
    else
        err = ncmpi_get_vara_int(ncid, varid[0], start, count, buf[0]);
#else
    err = ncmpi_iget_vara_int(ncid, varid[0], start, count, buf[0], NULL);
#endif
    CHECK_ERROUT

    buf[1] = buf[0] + bufLen[0];
    bufLen[1] = count[0]*count[1];

#ifdef USE_BLOCKING_APIS
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid[1], start, count, buf[1]);
    else
        err = ncmpi_get_vara_int(ncid, varid[1], start, count, buf[1]);
#else
    err = ncmpi_iget_vara_int(ncid, varid[1], start, count, buf[1], NULL);
#endif
    CHECK_ERROUT

    if (nvars == 2) goto cmp_val;

    /* check fix-sized variables */
    count[0] = 1;

    buf[2] = buf[1] + bufLen[1];
    bufLen[2] = count[0]*count[1];

#ifdef USE_BLOCKING_APIS
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid[2], start+1, count+1, buf[2]);
    else
        err = ncmpi_get_vara_int(ncid, varid[2], start+1, count+1, buf[2]);
#else
    err = ncmpi_iget_vara_int(ncid, varid[2], start+1, count+1, buf[2], NULL);
#endif
    CHECK_ERROUT

    buf[3] = buf[2] + bufLen[2];
    bufLen[3] = count[0]*count[1];

#ifdef USE_BLOCKING_APIS
    if (coll_io)
        err = ncmpi_get_vara_int_all(ncid, varid[3], start+1, count+1, buf[3]);
    else
        err = ncmpi_get_vara_int(ncid, varid[3], start+1, count+1, buf[3]);
#else
    err = ncmpi_iget_vara_int(ncid, varid[3], start+1, count+1, buf[3], NULL);
#endif
    CHECK_ERROUT

cmp_val:
#ifndef USE_BLOCKING_APIS
    if (coll_io)
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    else
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERROUT
#endif

    off_val = 1;
    for (j=0; j<nvars; j++) {
        for (i=0; i<bufLen[j]; i++)
            CHECK_VAL(ncid, varid[j], i, buf[j][i], rank+i+off_val)
        off_val *= 10;
    }
    free(buf[0]);

err_out:
    return nerrs;
}

#define GET_HEADER_SIZE { \
    old_hsize   = hsize; \
    old_extent  = extent; \
    old_r_begin = r_begin; \
    old_h_free  = h_free; \
    old_v_free  = v_free; \
    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR \
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR \
    h_free = extent - hsize; \
    err = ncmpi_inq_varoffset(ncid, varid[0], &r_begin); CHECK_ERR \
    v_free = r_begin - (extent + fix_v_size); \
    if (verbose && rank == 0) { \
        printf("Line %d: header size   old = "OFFFMT" new = "OFFFMT"\n", \
               __LINE__,old_hsize, hsize); \
        printf("Line %d: header extent old = "OFFFMT" new = "OFFFMT"\n", \
               __LINE__,old_extent, extent); \
        printf("Line %d: header free   old = "OFFFMT" new = "OFFFMT"\n", \
               __LINE__,old_h_free, h_free); \
        printf("Line %d: record begin  old = "OFFFMT" new = "OFFFMT"\n", \
               __LINE__,old_r_begin, r_begin); \
        printf("Line %d: var free      old = "OFFFMT" new = "OFFFMT"\n", \
               __LINE__,old_v_free, v_free); \
    } \
}

#define CHECK_HEADER_SIZE { \
    if (hsize != exp_hsize) { \
        nerrs++; \
        printf("Error at line %d in %s: header size expecting "OFFFMT" but got "OFFFMT"\n", \
               __LINE__,__FILE__, exp_hsize, hsize); \
    } \
    if (extent != exp_extent) { \
        nerrs++; \
        printf("Error at line %d in %s: header extent expecting "OFFFMT" but got "OFFFMT"\n", \
               __LINE__,__FILE__, exp_extent, extent); \
    } \
    if (extent - hsize < exp_h_free) { \
        nerrs++; \
        printf("Error at line %d in %s: header free expecting "OFFFMT" but got "OFFFMT"\n", \
               __LINE__,__FILE__, exp_h_free, extent - hsize); \
    } \
    if (has_fix_vars && v_free != exp_v_free) { \
        nerrs++; \
        printf("Error at line %d in %s: v_free expecting "OFFFMT" but got "OFFFMT"\n", \
               __LINE__,__FILE__, exp_v_free, v_free); \
    } \
    if (r_begin != exp_r_begin) { \
        nerrs++; \
        printf("Error at line %d in %s: record begin expecting "OFFFMT" but got "OFFFMT"\n", \
               __LINE__,__FILE__, exp_r_begin, r_begin); \
    } \
    /* read variables back and check contents */ \
    nerrs += check_vars(comm, ncid, varid, coll_io); \
    if (nerrs > 0) { \
        printf("Error at line %d in %s: check_vars failed\n", \
               __LINE__,__FILE__); \
        goto err_out; \
    } \
}

#define GROW_METADATA(growth) { \
    MPI_Offset len; \
    err = ncmpi_inq_attlen(ncid, NC_GLOBAL, "attr", &len); \
    len += growth; \
    char *attr = (char*) calloc(len + growth, 1); \
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "attr", len, attr);\
    CHECK_ERR \
    free(attr); \
    if (verbose && rank == 0) \
        printf("Line %d: grow header size from "OFFFMT" to "OFFFMT"\n", \
               __LINE__,hsize, hsize+growth); \
}

#define CHECK_ALIGNMENTS { \
    /* hints set in MPI info precede ncmpi__enddef */ \
    v_align = (env_v_align) ? env_v_align : \
              (info_v_align) ? info_v_align : \
              (v_align > 0) ? v_align : 512; \
    r_align = (env_r_align) ? env_r_align : \
              (info_r_align) ? info_r_align : \
              (r_align > 0) ? r_align : \
              (has_fix_vars) ? 4 : 512;\
    if (h_minfree == -1) h_minfree = 0; \
    if (v_minfree == -1) v_minfree = 0; \
    exp_hsize  = old_hsize + increment; \
    exp_extent = exp_hsize + h_minfree; \
    if (has_fix_vars) { \
        exp_extent = MAX(exp_extent, old_extent); \
        exp_extent = RNDUP(exp_extent, v_align); \
        exp_r_begin = exp_extent + fix_v_size + v_minfree; \
        exp_r_begin = MAX(exp_r_begin, old_r_begin); \
        exp_r_begin = RNDUP(exp_r_begin, r_align); \
    } else { \
        exp_r_begin = MAX(exp_extent, old_r_begin); \
        exp_r_begin = RNDUP(exp_r_begin, r_align); \
        exp_extent = exp_r_begin; \
    } \
    exp_h_free = exp_extent - exp_hsize; \
    exp_v_free = exp_r_begin - (exp_extent + fix_v_size); \
    CHECK_HEADER_SIZE \
}

#define PRINT_HINTS \
    if (verbose && rank == 0) { \
        printf("\n========================================\n"); \
        printf("  Line %d hsize %lld extent %lld r_begin %lld has_fix_vars %d\n", \
               __LINE__,hsize,extent,r_begin, has_fix_vars); \
        printf("  Line %d ncmpi__enddef() increment %lld h_minfree %lld v_align %lld v_minfree %lld r_align %lld\n", \
               __LINE__, increment, h_minfree, v_align, v_minfree, r_align); \
    }

/* test alignments hints
 * 1. set in environment variable PNETCDF_HINTS,
 * 2. set in ncmpi__enddef()
 * 3. set in MPI File info, which is passed into ncmpi_create/ncmpi_open
 * Note precedence of hints: PNETCDF_HINTS > ncmpi__enddef() > MPI info.
 */
static int
tst_fmt(const char *out_path,
        int         coll_io,
        MPI_Info    global_info,
        int         has_fix_vars,
        MPI_Offset *env_align,  /* [3] 0 means unset in PNETCDF_HINTS */
        MPI_Offset *info_align) /* [3] 0 means unset in MPI info */
{
    int i, rank, nprocs, ncid, err, nerrs=0;
    int *buf[4], dimid[3], varid[4];
    MPI_Info info=MPI_INFO_NULL;
    MPI_Offset bufLen[4], start[2], count[2], increment, fix_v_size;

    MPI_Offset hsize=0, old_hsize=-1, exp_hsize=-1;
    MPI_Offset extent=0, old_extent=-1, exp_extent=-1;
    MPI_Offset h_free=0, old_h_free, exp_h_free;
    MPI_Offset v_free=0, old_v_free, exp_v_free;
    MPI_Offset r_begin=0, old_r_begin=-1, exp_r_begin=-1;
    MPI_Offset h_minfree, v_align, v_minfree, r_align;
    MPI_Offset env_h_align=0, env_v_align=0, env_r_align=0;
    MPI_Offset info_h_align=0, info_v_align=0, info_r_align=0;

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Info_dup(global_info, &info);

    if (verbose && rank == 0)
        printf("---- has_fix_vars=%d env_align=%s info_align=%s\n",
               has_fix_vars,(env_align==NULL)?"NULL":"SET",
               (info_align==NULL)?"NULL":"SET");

    if (env_align != NULL) {
        env_h_align = env_align[0];  /* 0 means unset in PNETCDF_HINTS */
        env_v_align = env_align[1];  /* 0 means unset in PNETCDF_HINTS */
        env_r_align = env_align[2];  /* 0 means unset in PNETCDF_HINTS */
        if (env_v_align == 0) env_v_align = env_h_align;
    }
    if (info_align != NULL) {
        char str[16];
        info_h_align = info_align[0]; /* 0 means unset in MPI info */
        info_v_align = info_align[1]; /* 0 means unset in MPI info */
        info_r_align = info_align[2]; /* 0 means unset in MPI info */
        if (info_h_align) {
            sprintf(str, OFFFMT, info_h_align);
            MPI_Info_set(info, "nc_header_align_size", str);
        }
        if (info_v_align) {
            sprintf(str, OFFFMT, info_v_align);
            MPI_Info_set(info, "nc_var_align_size", str);
        }
        if (info_r_align) {
            sprintf(str, OFFFMT, info_r_align);
            MPI_Info_set(info, "nc_record_align_size", str);
        }
        if (info_v_align == 0) info_v_align = info_h_align;
    }
    if (verbose && rank == 0)
        printf("---- has_fix_vars=%d env_align="OFFFMT" "OFFFMT" "OFFFMT" info_align="OFFFMT" "OFFFMT" "OFFFMT"\n",
               has_fix_vars,env_h_align,env_v_align,env_r_align,
               info_h_align,info_v_align,info_r_align);

    /* create a new file */
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim", LEN*nprocs, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "ta", NC_INT, 2, dimid,   &varid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "tb", NC_INT, 2, dimid,   &varid[1]); CHECK_ERR
    if (has_fix_vars) {
        err = ncmpi_def_var(ncid, "fa", NC_INT, 1, dimid+1, &varid[2]);
        CHECK_ERR
        err = ncmpi_def_var(ncid, "fb", NC_INT, 1, dimid+1, &varid[3]);
        CHECK_ERR
    }
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "attr", 0, NULL); CHECK_ERR

    increment = 0; h_minfree = v_minfree = v_align = r_align = -1;
    PRINT_HINTS
    err = ncmpi_enddef(ncid); CHECK_ERR

    buf[0] = (int*) malloc(sizeof(int) * 2 * LEN * 4);

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* write to all variables, 2 records */
    start[0] = 0; start[1] = rank * LEN;
    count[0] = 2; count[1] = LEN;

    bufLen[0] = count[0]*count[1];
    for (i=0; i<bufLen[0]; i++) buf[0][i] = rank + i + 1;

#ifdef USE_BLOCKING_APIS
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf[0]);
    else
        err = ncmpi_put_vara_int(ncid, varid[0], start, count, buf[0]);
#else
    err = ncmpi_iput_vara_int(ncid, varid[0], start, count, buf[0], NULL);
#endif
    CHECK_ERR

    buf[1] = buf[0] + bufLen[0];
    bufLen[1] = count[0]*count[1];
    for (i=0; i<bufLen[1]; i++) buf[1][i] = rank + i + 10;

#ifdef USE_BLOCKING_APIS
    if (coll_io)
        err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf[1]);
    else
        err = ncmpi_put_vara_int(ncid, varid[1], start, count, buf[1]);
#else
    err = ncmpi_iput_vara_int(ncid, varid[1], start, count, buf[1], NULL);
#endif
    CHECK_ERR

    if (has_fix_vars) {
        buf[2] = buf[1] + bufLen[1];
        bufLen[2] = count[1];
        for (i=0; i<bufLen[2]; i++) buf[2][i] = rank + i + 100;

#ifdef USE_BLOCKING_APIS
        if (coll_io)
            err = ncmpi_put_vara_int_all(ncid, varid[2], start+1, count+1, buf[2]);
        else
            err = ncmpi_put_vara_int(ncid, varid[2], start+1, count+1, buf[2]);
#else
        err = ncmpi_iput_vara_int(ncid, varid[2], start+1, count+1, buf[2], NULL);
#endif
        CHECK_ERR

        buf[3] = buf[2] + bufLen[2];
        bufLen[3] = count[1];
        for (i=0; i<bufLen[3]; i++) buf[3][i] = rank + i + 1000;

#ifdef USE_BLOCKING_APIS
        if (coll_io)
            err = ncmpi_put_vara_int_all(ncid, varid[3], start+1, count+1, buf[3]);
        else
            err = ncmpi_put_vara_int(ncid, varid[3], start+1, count+1, buf[3]);
#else
        err = ncmpi_iput_vara_int(ncid, varid[3], start+1, count+1, buf[3], NULL);
#endif
        CHECK_ERR
        fix_v_size = sizeof(int) * LEN * nprocs * 2;
    }
    else
        fix_v_size = 0;

#ifndef USE_BLOCKING_APIS
    if (coll_io)
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    else
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR
#endif

    GET_HEADER_SIZE

    old_hsize   = hsize;
    old_extent  = extent;
    old_r_begin = r_begin;

    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR
    increment = 0;

    /* explicitly disable header alignment */
    h_minfree = 0;  /* header free space */
    v_align   = 4;  /* alignment for variable section (also header extent) */
    v_minfree = 0;  /* free space between fixed and record variable sections */
    r_align   = 4;  /* alignment for record variable section */

    PRINT_HINTS
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE

    /* expect nothing changed */
    exp_hsize   = old_hsize;
    exp_extent  = old_extent;
    exp_h_free  = old_h_free;
    exp_v_free  = old_v_free;
    exp_r_begin = old_r_begin;
    CHECK_HEADER_SIZE

    /* file sync before reading */
    if (!coll_io) {
        err = ncmpi_sync(ncid);
        CHECK_ERR
    }
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* reopen the file and check file header size and extent */
    err = ncmpi_open(comm, out_path, NC_WRITE, info, &ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_varid(ncid, "ta", &varid[0]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "tb", &varid[1]); CHECK_ERR
    if (has_fix_vars) {
        err = ncmpi_inq_varid(ncid, "fa", &varid[2]); CHECK_ERR
        err = ncmpi_inq_varid(ncid, "fb", &varid[3]); CHECK_ERR
    }

    GET_HEADER_SIZE

    /* expect nothing changed */
    exp_hsize   = old_hsize;
    exp_extent  = old_extent;
    exp_h_free  = old_h_free;
    exp_v_free  = old_v_free;
    exp_r_begin = old_r_begin;
    CHECK_HEADER_SIZE

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    /* grow header size */
    increment = extent - hsize + 8;
    GROW_METADATA(increment)

    /* exit define mode */
    h_minfree = v_minfree = v_align = r_align = -1;
    PRINT_HINTS
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE
    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    /* add a header free space, disable header alignment */
    increment = 0;
    h_minfree = 76;
    v_minfree = 44;
    v_align = 4;
    r_align = 4;
    PRINT_HINTS
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE
    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    /* grow header size */
    increment = extent - hsize - 16;
    GROW_METADATA(increment);

    /* disable header alignment */
    h_minfree = 0;
    v_minfree = 31;
    v_align = 4;
    r_align = 4;
    PRINT_HINTS
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE
    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    /* grow header size */
    increment = extent - hsize - 4;
    GROW_METADATA(increment);

    /* exit define mode */
    h_minfree = v_minfree = v_align = r_align = -1;
    PRINT_HINTS
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE
    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    increment = 0;
    h_minfree = 0;
    v_minfree = 31;
    r_align = 4;

    /* increase v_align */
    v_align = extent + 4;
    PRINT_HINTS
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE

    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    increment = 0;
    h_minfree = 0;
    v_minfree = 31;
    v_align = 0;

    /* increase r_align */
    r_align = r_begin + 4;
    PRINT_HINTS
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE

    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    /* add nothing */
    increment = 0; h_minfree = v_minfree = v_align = r_align = -1;
    PRINT_HINTS
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE

    CHECK_ALIGNMENTS

err_out:
    err = ncmpi_close(ncid); CHECK_ERR
    free(buf[0]);
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char str[256], *saved_env;
    int rank, err, nerrs=0, has_fix_vars;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Offset env_align[3], info_align[3];

    MPI_Comm_rank(comm, &rank);

    verbose = 0;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* retrieve value of environment variable PNETCDF_HINTS */
    saved_env = getenv("PNETCDF_HINTS");
    if (verbose && rank == 0 && saved_env != NULL)
        printf("PNETCDF_HINTS=%s\n",saved_env);

    /* No alignment hints should be set in the environment variable
     * PNETCDF_HINTS (there can be other kinds) and the info object passed into
     * this subroutine (there can be other hints) before running this test
     * program. Check seq_runs.sh and parallel_run.sh first before running this
     * test program.
     */
    for (has_fix_vars=1; has_fix_vars>=0; has_fix_vars--) {

        /* Test when there is no alignment hints set at all */
        nerrs += tst_fmt(out_path, coll_io, info, has_fix_vars, NULL, NULL);
        if (nerrs > 0) goto err_out;

        /* Test when there is no alignment hints set in environment variable
         * PNETCDF_HINTS and set alignment hints in MPI Info object.
         */
        info_align[0] = 28;  /*  7 x 4 */
        info_align[1] = 44;  /* 11 x 4 */
        info_align[2] = 52;  /* 13 x 4 */

        nerrs += tst_fmt(out_path, coll_io, info, has_fix_vars, NULL, info_align);
        if (nerrs > 0) goto err_out;

        /* Set hints in environment variable PNETCDF_HINTS, but no hints set in
         * MPI Info object.
         */
        env_align[0] = 68;  /* 17 x 4 */
        env_align[1] = 76;  /* 19 x 4 */
        env_align[2] = 92;  /* 23 x 4 */
        sprintf(str, "nc_header_align_size="OFFFMT";nc_var_align_size="OFFFMT";nc_record_align_size="OFFFMT"\n",
                env_align[0], env_align[1], env_align[2]);
        setenv("PNETCDF_HINTS", str, 1);

        nerrs += tst_fmt(out_path, coll_io, info, has_fix_vars, env_align, NULL);
        if (nerrs > 0) goto err_out;

        /* Test if the alignment hints set in environment variable
         * PNETCDF_HINTS take precedence over hints set in MPI Info object, when
         * Hints are both set in environment variable PNETCDF_HINTS and in MPI
         * Info object.
         */
        nerrs += tst_fmt(out_path, coll_io, info, has_fix_vars, env_align, info_align);
        if (nerrs > 0) goto err_out;

        /* Test a different set of alignment hints set in environment variable
         * PNETCDF_HINTS.
         */
        env_align[0] = 68;  /* 17 x 4 */
        env_align[1] = 0;
        env_align[2] = 92;  /* 23 x 4 */
        sprintf(str, "nc_header_align_size="OFFFMT";nc_record_align_size="OFFFMT"\n",
                env_align[0], env_align[2]);
        setenv("PNETCDF_HINTS", str, 1);

        nerrs += tst_fmt(out_path, coll_io, info, has_fix_vars, env_align, info_align);
        if (nerrs > 0) goto err_out;

        /* Test a different set of alignment hints set in environment variable
         * PNETCDF_HINTS.
         */
        env_align[0] = 0;
        env_align[1] = 76;  /* 19 x 4 */
        env_align[2] = 92;  /* 23 x 4 */
        sprintf(str, "nc_var_align_size="OFFFMT";nc_record_align_size="OFFFMT"\n",
                env_align[1], env_align[2]);
        setenv("PNETCDF_HINTS", str, 1);

        nerrs += tst_fmt(out_path, coll_io, info, has_fix_vars, env_align, info_align);
        if (nerrs > 0) goto err_out;

        /* Test a different set of alignment hints set in environment variable
         * PNETCDF_HINTS.
         */
        env_align[0] = 0;
        env_align[1] = 76;  /* 19 x 4 */
        env_align[2] = 0;
        sprintf(str, "nc_var_align_size="OFFFMT"\n", env_align[1]);
        setenv("PNETCDF_HINTS", str, 1);

        nerrs += tst_fmt(out_path, coll_io, info, has_fix_vars, env_align, info_align);
        if (nerrs > 0) goto err_out;

        /* Test a different set of alignment hints set in environment variable
         * PNETCDF_HINTS.
         */
        env_align[0] = 0;  /* 17 x 4 */
        env_align[1] = 0;
        env_align[2] = 92;  /* 23 x 4 */
        sprintf(str, "nc_record_align_size="OFFFMT"\n", env_align[2]);
        setenv("PNETCDF_HINTS", str, 1);

        nerrs += tst_fmt(out_path, coll_io, info, has_fix_vars, env_align, info_align);
        if (nerrs > 0) goto err_out;

        /* restore the original value set in environment variable PNETCDF_HINTS */
        if (saved_env != NULL) setenv("PNETCDF_HINTS", saved_env, 1);
        else                   unsetenv("PNETCDF_HINTS");
    }

err_out:
    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1;   /* test intra-node aggregation */
    opt.drv      = 1;   /* test PNCIO driver */
    opt.ind      = 1;   /* test hint romio_no_indep_rw */
    opt.chk      = 300; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1;   /* test burst-buffering feature */
    opt.mod      = 1;   /* test independent data mode */
    opt.hdr_diff = 1;   /* run ncmpidiff for file header only */
    opt.var_diff = 1;   /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "header alignment", opt, test_io);

    MPI_Finalize();

    return err;
}
