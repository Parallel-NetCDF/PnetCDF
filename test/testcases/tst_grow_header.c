/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests header extent growth by re-entering the define mode
 * multiple times and add more fix-sized and record variables.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_grow_header tst_grow_header.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./tst_grow_header testfile.nc
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

#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

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
check_vars(MPI_Comm    comm,
           int         ncid,
           int         coll_io,
           MPI_Offset *start,
           MPI_Offset *count)
{
    int id, j, r, nerrs=0, err, rank, *buf=NULL, rec_dim, dimids[3], nvars;
    MPI_Offset nrecs;

    MPI_Comm_rank(comm, &rank);

    err = ncmpi_inq_nvars(ncid, &nvars); CHECK_ERROUT
    err = ncmpi_inq_dimid(ncid, "time", &rec_dim); CHECK_ERROUT
    err = ncmpi_inq_dimlen(ncid, rec_dim, &nrecs); CHECK_ERROUT

    buf = (int*) malloc(sizeof(int) * count[1] * count[2]);

    for (id=0; id<nvars; id++) { /* check fix-sized variables first */
        err = ncmpi_inq_vardimid(ncid, id, dimids); CHECK_ERROUT
        if (dimids[0] == rec_dim) continue;

        if (coll_io)
            err = ncmpi_get_vara_int_all(ncid, id, start+1, count+1, buf);
        else
            err = ncmpi_get_vara_int(ncid, id, start+1, count+1, buf);
        CHECK_ERR

        for (j=0; j<count[1] * count[2]; j++) {
            int expect = rank + j + id * 10;
            if (buf[j] != expect) {
                char name[64];
                err = ncmpi_inq_varname(ncid, id, name);
                printf("Error: fix var %s [%d] expect %d but got %d\n",
                       name, j, expect, buf[j]);
                nerrs++;
                break;
            }
        }
        MPI_Allreduce(&nerrs, &err, 1, MPI_INT, MPI_MAX, comm);
        if (err > 0) goto err_out;
    }

    for (id=0; id<nvars; id++) { /* check record variables only */
        err = ncmpi_inq_vardimid(ncid, id, dimids); CHECK_ERROUT
        if (dimids[0] != rec_dim) continue;

        for (r=0; r<nrecs; r++) {
            start[0] = r;
            if (coll_io)
                err = ncmpi_get_vara_int_all(ncid, id, start, count, buf);
            else
                err = ncmpi_get_vara_int(ncid, id, start, count, buf);
            CHECK_ERR
            for (j=0; j<count[1] * count[2]; j++) {
                int expect = rank + j + id * 10 + r;
                if (buf[j] != expect) {
                    char name[64];
                    err = ncmpi_inq_varname(ncid, id, name);
                    printf("Error: rec var %s [%d][%d] expect %d but got %d\n",
                        name, r, j, expect, buf[j]);
                    nerrs++;
                    break;
                }
            }
            MPI_Allreduce(&nerrs, &err, 1, MPI_INT, MPI_MAX, comm);
            if (err > 0) goto err_out;
        }
    }

err_out:
    if (buf != NULL) free(buf);
    return err;
}

#define GET_HEADER_SIZE { \
    old_hsize   = hsize; \
    old_extent  = extent; \
    err = ncmpi_inq_header_size(ncid, &hsize); CHECK_ERR \
    err = ncmpi_inq_header_extent(ncid, &extent); CHECK_ERR \
    if (verbose && rank == 0) { \
        printf("Line %d: header size   old = "OFFFMT" new = "OFFFMT"\n", \
               __LINE__,old_hsize, hsize); \
        printf("Line %d: header extent old = "OFFFMT" new = "OFFFMT"\n", \
               __LINE__,old_extent, extent); \
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
    /* read variables back and check contents */ \
    nerrs += check_vars(comm, ncid, coll_io, start, count); \
    if (nerrs > 0) { \
        printf("Error at line %d in %s: variable contents unexpected\n", \
               __LINE__,__FILE__ ); \
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

#define ADD_FIX_VAR(id) { \
    char name[64]; \
    sprintf(name, "v%d", id); \
    err = ncmpi_def_var(ncid, name, NC_INT, 2, dimid+1, varid+id); \
    CHECK_ERR \
    if (verbose && rank == 0) \
        printf("Line %d: add fix var %s of size %6zd\n", \
               __LINE__,name, sizeof(int)*10*nprocs*LEN); \
}
#define ADD_REC_VAR(id) { \
    char name[64]; \
    sprintf(name, "v%d", id); \
    err = ncmpi_def_var(ncid, name, NC_INT, 3, dimid, varid+id); \
    CHECK_ERR \
    if (verbose && rank == 0) \
        printf("Line %d: add rec var %s of size %6zd\n", \
               __LINE__,name, sizeof(int)*10*nprocs*LEN); \
}
#define WRITE_FIX_VAR(id) { \
    for (i=0; i<count[1] * count[2]; i++) buf[i] = rank + i + id * 10; \
    err = ncmpi_iput_vara_int(ncid, varid[id], start+1, count+1, buf, NULL); \
    CHECK_ERR \
}
#define WRITE_REC_VAR(id) { \
    for (i=0; i<count[1] * count[2]; i++) \
        buf[i] = rank + i + id * 10 + (int)start[0]; \
    err = ncmpi_iput_vara_int(ncid, varid[id], start, count, buf, NULL); \
    CHECK_ERR \
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    global_info)
{
    int i, rank, nprocs, ncid, err, nerrs=0;
    int *buf=NULL, dimid[3], varid[16];
    MPI_Info info=MPI_INFO_NULL;
    MPI_Offset start[3], count[3], increment;

    MPI_Offset hsize=0, old_hsize, exp_hsize;
    MPI_Offset extent=0, old_extent, exp_extent;

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Info_dup(global_info, &info);
    MPI_Info_set(info, "nc_var_align_size", "1024");

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file */
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", 10*nprocs, &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", LEN,       &dimid[2]); CHECK_ERR

    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "attr", 0, NULL); CHECK_ERR

    ADD_REC_VAR(0)
    ADD_REC_VAR(1)
    ADD_FIX_VAR(2)
    ADD_FIX_VAR(3)

    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* write to all variables, 1 record */
    start[0] = 0; start[1] = rank * 10; start[2] = 0;
    count[0] = 1; count[1] = 10;        count[2] = LEN;

    buf = (int*) malloc(sizeof(int) * count[1] * count[2]);

    /* write record variables */
    start[0] = 0; WRITE_REC_VAR(0)
    start[0] = 0; WRITE_REC_VAR(1)

    /* write fix-sized variables */
    WRITE_FIX_VAR(2)
    WRITE_FIX_VAR(3)

    if (coll_io)
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    else
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* reopen the file */
    err = ncmpi_open(comm, out_path, NC_WRITE, info, &ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    err = ncmpi_inq_varid(ncid, "v0", &varid[0]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "v1", &varid[1]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "v2", &varid[2]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "v3", &varid[3]); CHECK_ERR

    GET_HEADER_SIZE

    /* enter redefine mode to grow header size without growing header extent,
     * adding no new variable
     */
    err = ncmpi_redef(ncid); CHECK_ERR
    increment = MIN(8, extent - hsize);
    GROW_METADATA(increment)
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE

    exp_hsize  = old_hsize + increment;
    exp_extent = old_extent;
    CHECK_HEADER_SIZE

    /* enter redefine mode to grow header size without growing header extent,
     * adding a new fix-sized variable
     */
    err = ncmpi_redef(ncid); CHECK_ERR
    ADD_FIX_VAR(4)
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    WRITE_FIX_VAR(4)

    if (coll_io)
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    else
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    GET_HEADER_SIZE

    exp_hsize  = hsize;
    exp_extent = old_extent;
    CHECK_HEADER_SIZE

    /* enter redefine mode to grow header size without growing header extent,
     * adding a new record variable
     */
    err = ncmpi_redef(ncid); CHECK_ERR
    ADD_REC_VAR(5)
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    start[0] = 0; WRITE_REC_VAR(5) /* write 1st record */
    start[0] = 1; WRITE_REC_VAR(5) /* write 2nd record */
    start[0] = 1; WRITE_REC_VAR(0) /* write 2nd record */
    start[0] = 1; WRITE_REC_VAR(1) /* write 2nd record */

    if (coll_io)
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    else
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    GET_HEADER_SIZE

    exp_hsize  = hsize;
    exp_extent = old_extent;
    CHECK_HEADER_SIZE

    /* enter redefine mode to grow header extent
     * and adding a new fix-sized variable
     */
    err = ncmpi_redef(ncid); CHECK_ERR
    increment = extent - hsize;
    GROW_METADATA(increment)
    ADD_FIX_VAR(6)
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    WRITE_FIX_VAR(6)

    if (coll_io)
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    else
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    GET_HEADER_SIZE

    exp_hsize  = hsize;
    exp_extent = 1024 * 2;
    CHECK_HEADER_SIZE

    /* enter redefine mode to grow header extent
     * and adding a new record variable
     */
    err = ncmpi_redef(ncid); CHECK_ERR
    increment = extent - hsize;
    GROW_METADATA(increment)
    ADD_REC_VAR(7)
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    start[0] = 0; WRITE_REC_VAR(7) /* write 1st record */
    start[0] = 1; WRITE_REC_VAR(7) /* write 2nd record */
    start[0] = 2; WRITE_REC_VAR(7) /* write 3rd record */
    start[0] = 2; WRITE_REC_VAR(0) /* write 3rd record */
    start[0] = 2; WRITE_REC_VAR(1) /* write 3rd record */
    start[0] = 2; WRITE_REC_VAR(5) /* write 3rd record */

    if (coll_io)
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    else
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    GET_HEADER_SIZE

    exp_hsize  = hsize;
    exp_extent = 1024 * 3;
    CHECK_HEADER_SIZE

    /* enter redefine mode to grow header extent
     * and adding a new fix-sized variable and a record variable
     */
    err = ncmpi_redef(ncid); CHECK_ERR
    increment = extent - hsize;
    GROW_METADATA(increment)
    ADD_REC_VAR(8)
    ADD_FIX_VAR(9)
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    WRITE_FIX_VAR(9)

    start[0] = 0; WRITE_REC_VAR(8) /* write 1st record */
    start[0] = 1; WRITE_REC_VAR(8) /* write 2nd record */
    start[0] = 2; WRITE_REC_VAR(8) /* write 3rd record */

    if (coll_io)
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    else
        err = ncmpi_wait(ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_ERR

    GET_HEADER_SIZE

    exp_hsize  = hsize;
    exp_extent = 1024 * 4;
    CHECK_HEADER_SIZE

    /* enter redefine mode to grow header extent
     * but add nothing else
     */
    err = ncmpi_redef(ncid); CHECK_ERR
    increment = extent - hsize + 8;
    GROW_METADATA(increment)
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    GET_HEADER_SIZE

    exp_hsize  = old_hsize + increment;
    exp_extent = 1024 * 5;
    CHECK_HEADER_SIZE

    err = ncmpi_close(ncid); CHECK_ERR

err_out:
    if (buf != NULL) free(buf);
    MPI_Info_free(&info);

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
    opt.chk      = 4096; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1;    /* test burst-buffering feature */
    opt.mod      = 1;    /* test independent data mode */
    opt.hdr_diff = 1;    /* run ncmpidiff for file header only */
    opt.var_diff = 1;    /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "header grow", opt, test_io);

    MPI_Finalize();

    return err;
}
