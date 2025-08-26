/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests all alignment features available from PnetCDF:
 * 1. v_align: header extent alignment (starting file offset of data section)
 * 2. h_minfree: header free space
 * 3. r_align: record variable section alignment
 * 4. v_minfree: free space between record variable section and the end of last
 *    fix-sized variable.
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
#include <pnetcdf.h>

#include <testutils.h>

#define LEN 101

#define RNDUP(x, unit) ((((x) + (unit) - 1) / (unit)) * (unit))

static int verbose;

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

/*----< check_vars() >-------------------------------------------------------*/
/* read back variables from file and check their contents */
static int
check_vars(MPI_Comm comm, int ncid, int *varid)
{
    int i, nerrs=0, err, rank, *buf=NULL, nvars;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);

    err = ncmpi_inq_nvars(ncid, &nvars); CHECK_ERROUT

    /* check fix-sized variables */
    start[0] = 0;
    start[1] = rank * LEN;
    count[1] = LEN;

    buf = (int*) malloc(sizeof(int) * 2 * LEN);

    /* check record variables */
    count[0] = 2;

    for (i=0; i<count[0]*count[1]; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[0], start,   count,   buf);
    CHECK_ERROUT
    for (i=0; i<2*LEN; i++)
        CHECK_VAL(ncid, varid[0], i, buf[i], rank+i+1000)

    for (i=0; i<count[0]*count[1]; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[1], start,   count,   buf);
    CHECK_ERROUT
    for (i=0; i<2*LEN; i++)
        CHECK_VAL(ncid, varid[1], i, buf[i], rank+i+10000)

    if (nvars == 2) goto err_out;

    /* check fix-sized variables */
    count[0] = 1;

    for (i=0; i<count[0]*count[1]; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[2], start+1, count+1, buf);
    CHECK_ERROUT
    for (i=0; i<LEN; i++)
        CHECK_VAL(ncid, varid[2], i, buf[i], rank+i)

    for (i=0; i<count[0]*count[1]; i++) buf[i] = -1;
    err = ncmpi_get_vara_int_all(ncid, varid[3], start+1, count+1, buf);
    CHECK_ERROUT
    for (i=0; i<LEN; i++)
        CHECK_VAL(ncid, varid[3], i, buf[i], rank+i+100)

err_out:
    if (buf != NULL) free(buf);
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
    nerrs += check_vars(comm, ncid, varid); \
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
              (!has_fix_vars && env_r_align) ? env_r_align : \
              (info_v_align) ? info_v_align : v_align; \
    r_align = (env_r_align) ? env_r_align : \
              (info_r_align) ? info_r_align : r_align; \
    exp_hsize  = old_hsize + increment; \
    exp_extent = RNDUP(exp_hsize + h_minfree, v_align); \
    old_extent = RNDUP(old_extent, v_align); \
    exp_extent = (exp_extent < old_extent) ? old_extent : exp_extent; \
    if (has_fix_vars) { \
        exp_r_begin = RNDUP(exp_extent + fix_v_size + v_minfree, r_align); \
        old_r_begin = RNDUP(old_r_begin, r_align); \
    } \
    else /* v_minfree and r_align are ignored */ \
        exp_r_begin = exp_extent; \
    exp_r_begin = (exp_r_begin < old_r_begin) ? old_r_begin : exp_r_begin; \
    exp_h_free  = exp_extent - exp_hsize; \
    exp_v_free  = exp_r_begin - (exp_extent + fix_v_size); \
    CHECK_HEADER_SIZE \
}

/* test alignments hints
 * 1. set in environment variable PNETCDF_HINTS,
 * 2. set in ncmpi__enddef()
 * 3. set in MPI File info, which is passed into ncmpi_create/ncmpi_open
 * Note precedence of hints: PNETCDF_HINTS > ncmpi__enddef() > MPI info.
 */
static int
tst_fmt(char       *filename,
        int         cmode,
        int         has_fix_vars,
        MPI_Offset *env_align,  /* [3] 0 means unset in PNETCDF_HINTS */
        MPI_Offset *info_align) /* [3] 0 means unset in MPI info */
{
    int i, rank, nprocs, ncid, err, nerrs=0;
    int *buf, dimid[3], varid[4];
    MPI_Info info=MPI_INFO_NULL;
    MPI_Offset start[2], count[2], increment, fix_v_size;

    MPI_Offset hsize=0, old_hsize, exp_hsize;
    MPI_Offset extent=0, old_extent, exp_extent;
    MPI_Offset h_free=0, old_h_free, exp_h_free;
    MPI_Offset v_free=0, old_v_free, exp_v_free;
    MPI_Offset r_begin=0, old_r_begin, exp_r_begin;
    MPI_Offset h_minfree, v_align, v_minfree, r_align;
    MPI_Offset env_h_align=0, env_v_align=0, env_r_align=0;
    MPI_Offset info_h_align=0, info_v_align=0, info_r_align=0;

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (verbose && rank == 0)
        printf("---- cmode=%d has_fix_vars=%d env_align=%s info_align=%s\n",
               cmode,has_fix_vars,(env_align==NULL)?"NULL":"SET",
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
        MPI_Info_create(&info);
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
        printf("---- cmode=%d has_fix_vars=%d env_align="OFFFMT" "OFFFMT" "OFFFMT" info_align="OFFFMT" "OFFFMT" "OFFFMT"\n",
               cmode,has_fix_vars,env_h_align,env_v_align,env_r_align,
               info_h_align,info_v_align,info_r_align);


    /* create a new file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR

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

    err = ncmpi_enddef(ncid); CHECK_ERR

    /* write to all variables, 2 records */
    start[0] = 0; start[1] = rank * LEN;
    count[0] = 2; count[1] = LEN;

    buf = (int*) malloc(sizeof(int) * count[0] * count[1]);

    for (i=0; i<count[0] * count[1]; i++) buf[i] = rank + i + 1000;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, buf);
    CHECK_ERR
    for (i=0; i<count[0] * count[1]; i++) buf[i] = rank + i + 10000;
    err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, buf);
    CHECK_ERR

    if (has_fix_vars) {
        for (i=0; i<count[1]; i++) buf[i] = rank + i;
        err = ncmpi_put_vara_int_all(ncid, varid[2], start+1, count+1, buf);
        CHECK_ERR
        for (i=0; i<count[1]; i++) buf[i] = rank + i + 100;
        err = ncmpi_put_vara_int_all(ncid, varid[3], start+1, count+1, buf);
        CHECK_ERR
        fix_v_size = sizeof(int) * LEN * nprocs * 2;
    }
    else
        fix_v_size = 0;

    GET_HEADER_SIZE

    /* PnetCDF default v_align is 512 when no hints given */
    v_align = (env_v_align) ? env_v_align :
              (!has_fix_vars && env_r_align) ? env_r_align :
              (info_v_align) ? info_v_align : 512;
    r_align = (env_r_align) ? env_r_align : (info_r_align) ? info_r_align : 4;
    increment = 0;
    h_minfree = 0;
    v_minfree = 0;

    exp_hsize = hsize;
    exp_extent = RNDUP(hsize, v_align);
    if (has_fix_vars)
        exp_r_begin = RNDUP(exp_extent + fix_v_size, r_align);
    else /* r_align and v_minfree are ignored */
        exp_r_begin = exp_extent;
    exp_h_free = exp_extent - exp_hsize;
    exp_v_free = v_free;
    CHECK_HEADER_SIZE

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR
    increment = 0;

    /* explicitly disable header alignment */
    h_minfree = 0;  /* header free space */
    v_align   = 4;  /* alignment for variable section (also header extent) */
    v_minfree = 0;  /* free space between fixed and record variable sections */
    r_align   = 4;  /* alignment for record variable section */

    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    GET_HEADER_SIZE

    /* expect nothing changed */
    exp_hsize   = old_hsize;
    exp_extent  = old_extent;
    exp_h_free  = old_h_free;
    exp_v_free  = old_v_free;
    exp_r_begin = old_r_begin;
    CHECK_HEADER_SIZE

    err = ncmpi_close(ncid); CHECK_ERR

    /* reopen the file and check file header size and extent */
    err = ncmpi_open(comm, filename, NC_WRITE, info, &ncid); CHECK_ERR

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
    h_minfree = 0;
    v_minfree = 0;
    v_align = 4;
    r_align = 4;
    err = ncmpi_enddef(ncid); CHECK_ERR

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
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

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
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    GET_HEADER_SIZE
    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    /* grow header size */
    increment = extent - hsize - 4;
    GROW_METADATA(increment);

    /* exit define mode */
    h_minfree = 0;
    v_minfree = 0;
    v_align = 4;
    r_align = 4;
    err = ncmpi_enddef(ncid); CHECK_ERR

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
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    GET_HEADER_SIZE

    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    increment = 0;
    h_minfree = 0;
    v_minfree = 31;
    v_align = 4;

    /* increase r_align */
    r_align = r_begin + 4;
    err = ncmpi__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
    CHECK_ERR

    GET_HEADER_SIZE

    CHECK_ALIGNMENTS

    /* enter redefine mode -------------------------------------------*/
    err = ncmpi_redef(ncid); CHECK_ERR

    increment = 0;
    h_minfree = 0;
    v_minfree = 0;
    v_align = 4;
    r_align = 4;

    /* add nothing */
    err = ncmpi_enddef(ncid); CHECK_ERR

    GET_HEADER_SIZE

    CHECK_ALIGNMENTS

err_out:
    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    return nerrs;
}

int main(int argc, char** argv)
{
    char filename[256], str[256];
    int i, rank, err, nerrs=0, cmode[3], has_fix_vars;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Offset env_align[3], info_align[3];

    MPI_Init(&argc, &argv);
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
        if (verbose) printf("\n");
    }
    cmode[0] = 0;
    cmode[1] = NC_64BIT_OFFSET;
    cmode[2] = NC_64BIT_DATA;

    for (has_fix_vars=1; has_fix_vars>=0; has_fix_vars--) {
        /* No hints set in environment variable PNETCDF_HINTS.
         * No hints set in MPI Info object.
         */
        unsetenv("PNETCDF_HINTS");
        for (i=0; i<3; i++) {
            nerrs += tst_fmt(filename, cmode[i], has_fix_vars, NULL, NULL);
            if (nerrs > 0) goto main_exit;
        }

        info_align[0] = 28;  /*  7 x 4 */
        info_align[1] = 44;  /* 11 x 4 */
        info_align[2] = 52;  /* 13 x 4 */

        /* No hints set in environment variable PNETCDF_HINTS.
         * Hints set in MPI Info object.
         */
        for (i=0; i<3; i++) {
            nerrs += tst_fmt(filename, cmode[i], has_fix_vars, NULL,
                             info_align);
            if (nerrs > 0) goto main_exit;
        }

        env_align[0] = 68;  /* 17 x 4 */
        env_align[1] = 76;  /* 19 x 4 */
        env_align[2] = 92;  /* 23 x 4 */
        sprintf(str, "nc_header_align_size="OFFFMT";nc_var_align_size="OFFFMT";nc_record_align_size="OFFFMT"\n",
                env_align[0], env_align[1], env_align[2]);
        setenv("PNETCDF_HINTS", str, 1);

        /* Set hints in environment variable PNETCDF_HINTS.
         * No hints set in MPI Info object.
         */
        for (i=0; i<3; i++) {
            nerrs += tst_fmt(filename, cmode[i], has_fix_vars, env_align, NULL);
            if (nerrs > 0) goto main_exit;
        }

        /* Set hints in environment variable PNETCDF_HINTS.
         * Set hints in MPI Info object.
         */
        for (i=0; i<3; i++) {
            nerrs += tst_fmt(filename, cmode[i], has_fix_vars, env_align,
                             info_align);
            if (nerrs > 0) goto main_exit;
        }

        env_align[0] = 68;  /* 17 x 4 */
        env_align[1] = 0;
        env_align[2] = 92;  /* 23 x 4 */
        sprintf(str, "nc_header_align_size="OFFFMT";nc_record_align_size="OFFFMT"\n",
                env_align[0], env_align[2]);
        setenv("PNETCDF_HINTS", str, 1);

        /* Set hints in environment variable PNETCDF_HINTS.
         * No hints set in MPI Info object.
         */
        for (i=0; i<3; i++) {
            nerrs += tst_fmt(filename, cmode[i], has_fix_vars, env_align,
                             info_align);
            if (nerrs > 0) goto main_exit;
        }

        env_align[0] = 0;
        env_align[1] = 76;  /* 19 x 4 */
        env_align[2] = 92;  /* 23 x 4 */
        sprintf(str, "nc_var_align_size="OFFFMT";nc_record_align_size="OFFFMT"\n",
                env_align[1], env_align[2]);
        setenv("PNETCDF_HINTS", str, 1);

        /* Set hints in environment variable PNETCDF_HINTS.
         * No hints set in MPI Info object.
         */
        for (i=0; i<3; i++) {
            nerrs += tst_fmt(filename, cmode[i], has_fix_vars, env_align,
                             info_align);
            if (nerrs > 0) goto main_exit;
        }

        env_align[0] = 0;
        env_align[1] = 76;  /* 19 x 4 */
        env_align[2] = 0;
        sprintf(str, "nc_var_align_size="OFFFMT"\n", env_align[1]);
        setenv("PNETCDF_HINTS", str, 1);

        /* Set hints in environment variable PNETCDF_HINTS.
         * No hints set in MPI Info object.
         */
        for (i=0; i<3; i++) {
            nerrs += tst_fmt(filename, cmode[i], has_fix_vars, env_align,
                             info_align);
            if (nerrs > 0) goto main_exit;
        }

        env_align[0] = 0;  /* 17 x 4 */
        env_align[1] = 0;
        env_align[2] = 92;  /* 23 x 4 */
        sprintf(str, "nc_record_align_size="OFFFMT"\n", env_align[2]);
        setenv("PNETCDF_HINTS", str, 1);

        /* Set hints in environment variable PNETCDF_HINTS.
         * No hints set in MPI Info object.
         */
        for (i=0; i<3; i++) {
            nerrs += tst_fmt(filename, cmode[i], has_fix_vars, env_align,
                             info_align);
            if (nerrs > 0) goto main_exit;
        }
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

main_exit:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

