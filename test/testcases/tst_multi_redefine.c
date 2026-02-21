/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests entering define modes multiple times, growing header
 * extension causing moving data section to plances of higher offsets, and
 * checking the contents of all variables defined.
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

static int verbose;

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef long long longlong;
typedef unsigned long long ulonglong;

#define END_DEF { \
    err = ncmpi_enddef(ncid); CHECK_ERROUT \
    if (!coll_io) { \
        err = ncmpi_begin_indep_data(ncid); \
        CHECK_ERR \
    } \
    err = ncmpi_inq_header_size(ncid, &new_hdr_size); CHECK_ERROUT \
    err = ncmpi_inq_header_extent(ncid, &new_hdr_ext); CHECK_ERROUT \
    if (verbose && rank == 0) { \
        printf("Add var %2d: header      size grows from %4lld to %4lld\n", \
               k, old_hdr_size, new_hdr_size); \
        if (new_hdr_ext > old_hdr_ext) \
            printf("Add var %2d: header extension grows from %4lld to %4lld\n", \
                k, old_hdr_ext, new_hdr_ext); \
    } \
    old_hdr_size = new_hdr_size; \
    old_hdr_ext  = new_hdr_ext; \
}

#define EXP_VAL(nn) ((j + nn + rank) % NC_MAX_BYTE)

#define DEF_VAR(xtype) { \
    err = ncmpi_redef(ncid); CHECK_ERROUT \
    sprintf(str, "fix_var_%d", k); \
    err = ncmpi_def_var(ncid, str, xtype, 2, dimids+1, &fix_varids[k]); \
    CHECK_ERROUT \
    sprintf(str, "attribute of fix-sized variable %d", k); \
    err = ncmpi_put_att_text(ncid, fix_varids[k], "attr", strlen(str), str); \
    CHECK_ERROUT \
    sprintf(str, "rec_var_%d", k); \
    err = ncmpi_def_var(ncid, str, xtype, 3, dimids, &rec_varids[k]); \
    CHECK_ERROUT \
    sprintf(str, "attribute of record variable %d", k); \
    err = ncmpi_put_att_text(ncid, rec_varids[k], "attr", strlen(str), str); \
    CHECK_ERROUT \
    END_DEF \
}

#define PUT_BUF_CHAR { \
    char *buf = (char*) malloc(sizeof(char) * nelems); \
    for (j=0; j<nelems; j++) buf[j] = EXP_VAL(k); \
    if (coll_io) \
        err = ncmpi_put_vara_text_all(ncid, fix_varids[k], start+1, count+1, buf); \
    else \
        err = ncmpi_put_vara_text(ncid, fix_varids[k], start+1, count+1, buf); \
    CHECK_ERROUT \
    if (coll_io) \
        err = ncmpi_put_vara_text_all(ncid, rec_varids[k], start, count, buf); \
    else \
        err = ncmpi_put_vara_text(ncid, rec_varids[k], start, count, buf); \
    CHECK_ERROUT \
    free(buf); \
}

#define PUT_BUF(itype) { \
    itype *buf = (itype*) malloc(sizeof(itype) * nelems); \
    for (j=0; j<nelems; j++) buf[j] = EXP_VAL(k); \
    if (coll_io) \
        err = ncmpi_put_vara_##itype##_all(ncid, fix_varids[k], start+1, count+1, buf); \
    else \
        err = ncmpi_put_vara_##itype(ncid, fix_varids[k], start+1, count+1, buf); \
    CHECK_ERROUT \
    if (coll_io) \
        err = ncmpi_put_vara_##itype##_all(ncid, rec_varids[k], start, count, buf); \
    else \
        err = ncmpi_put_vara_##itype(ncid, rec_varids[k], start, count, buf); \
    CHECK_ERROUT \
    free(buf); \
}

#define CHECK_BUF_CHAR { \
    /* file sync before reading */ \
    err = ncmpi_sync(ncid); \
    CHECK_ERR \
    MPI_Barrier(MPI_COMM_WORLD); \
    char *buf = (char*) malloc(sizeof(char) * nelems); \
    if (coll_io) \
        err = ncmpi_get_vara_text_all(ncid, fix_varids[j], start+1, count+1, buf); \
    else \
        err = ncmpi_get_vara_text(ncid, fix_varids[j], start+1, count+1, buf); \
    CHECK_ERROUT \
    for (m=0; m<nelems; m++) { \
        char exp = EXP_VAL(m); \
        if (buf[m] != exp) { \
            printf("Error: var %d expect buf[%d] = %d but got %d\n", \
                   j, m, (int)exp, (int)buf[m]); \
            nerrs++; \
            goto err_out; \
        } \
    } \
    memset(buf, 0, sizeof(char) * nelems); \
    if (coll_io) \
        err = ncmpi_get_vara_text_all(ncid, rec_varids[j], start, count, buf); \
    else \
        err = ncmpi_get_vara_text(ncid, rec_varids[j], start, count, buf); \
    CHECK_ERROUT \
    for (m=0; m<nelems; m++) { \
        char exp = EXP_VAL(m); \
        if (buf[m] != exp) { \
            printf("Error: var %d expect buf[%d] = %d but got %d\n", \
                   j, m, (int)exp, (int)buf[m]); \
            nerrs++; \
            goto err_out; \
        } \
    } \
    free(buf); \
}

#define CHECK_BUF \
    /* file sync before reading */ \
    err = ncmpi_sync(ncid); \
    CHECK_ERR \
    MPI_Barrier(MPI_COMM_WORLD); \
    for (j=0; j<k; j++) { \
        if (j % ntypes == 0) { \
            CHECK_BUF_CHAR \
            continue; \
        } \
        double *buf = (double*) malloc(sizeof(double) * nelems); \
        if (coll_io) \
            err = ncmpi_get_vara_double_all(ncid, fix_varids[j], start+1, count+1, buf); \
        else \
            err = ncmpi_get_vara_double(ncid, fix_varids[j], start+1, count+1, buf); \
        CHECK_ERROUT \
        for (m=0; m<nelems; m++) { \
            double exp = EXP_VAL(m); \
            if (buf[m] != exp) { \
                printf("Error: var %d expect buf[%d] = %.f but got %.f\n", \
                       j, m, exp, buf[m]); \
                nerrs++; \
                goto err_out; \
            } \
        } \
        memset(buf, 0, sizeof(double) * nelems); \
        if (coll_io) \
            err = ncmpi_get_vara_double_all(ncid, rec_varids[j], start, count, buf); \
        else \
            err = ncmpi_get_vara_double(ncid, rec_varids[j], start, count, buf); \
        CHECK_ERROUT \
        for (m=0; m<nelems; m++) { \
            double exp = EXP_VAL(m); \
            if (buf[m] != exp) { \
                printf("Error: var %d expect buf[%d] = %.f but got %.f\n", \
                       j, m, exp, buf[m]); \
                nerrs++; \
                goto err_out; \
            } \
        } \
        free(buf); \
    }

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char str[64];
    int i, j, k, m, rank, nprocs, ncid, err, nerrs=0, dimids[3];
    int *fix_varids=NULL, *rec_varids=NULL, ntypes, nRounds;
    MPI_Offset old_hdr_size, new_hdr_size;
    MPI_Offset old_hdr_ext, new_hdr_ext;
    MPI_Offset len_y, len_x, nelems, start[3], count[3];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    nRounds = NROUNDS;
    len_y = NY;
    len_x = NX;

    /* CDF-1 and 2 have 6 external data types, CDF-5 has 11 */
    ntypes = (format == NC_FORMAT_64BIT_DATA) ? 11 : 6;

    fix_varids = (int*) malloc(sizeof(int) * nRounds * ntypes);
    rec_varids = (int*) malloc(sizeof(int) * nRounds * ntypes);

    start[0] = 0;
    start[1] = 0;
    start[2] = rank * len_x;
    count[0] = 1;
    count[1] = len_y;
    count[2] = len_x;
    nelems = count[1] * count[2];

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERROUT

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
        /* add 2 new variables of type NC_CHAR */
        DEF_VAR(NC_CHAR)
        PUT_BUF_CHAR
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_BYTE */
        DEF_VAR(NC_BYTE)
        PUT_BUF(schar)
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_SHORT */
        DEF_VAR(NC_SHORT)
        PUT_BUF(short)
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_INT */
        DEF_VAR(NC_INT)
        PUT_BUF(int)
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_FLOAT */
        DEF_VAR(NC_FLOAT)
        PUT_BUF(float)
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_DOUBLE */
        DEF_VAR(NC_DOUBLE)
        PUT_BUF(double)
        k++;
        CHECK_BUF

        if (format == NC_FORMAT_64BIT_DATA) {
            /* add 2 new variables of type NC_UBYTE */
            DEF_VAR(NC_UBYTE)
            PUT_BUF(uchar)
            k++;
            CHECK_BUF

            /* add 2 new variables of type NC_USHORT */
            DEF_VAR(NC_USHORT)
            PUT_BUF(ushort)
            k++;
            CHECK_BUF

            /* add 2 new variables of type NC_UINT */
            DEF_VAR(NC_UINT)
            PUT_BUF(uint)
            k++;
            CHECK_BUF

            /* add 2 new variables of type NC_INT64 */
            DEF_VAR(NC_INT64)
            PUT_BUF(longlong)
            k++;
            CHECK_BUF

            /* add 2 new variables of type NC_UINT64 */
            DEF_VAR(NC_UINT64)
            PUT_BUF(ulonglong)
            k++;
            CHECK_BUF
        }
    }

    err = ncmpi_close(ncid); CHECK_ERROUT

err_out:
    if (fix_varids != NULL) free(fix_varids);
    if (rec_varids != NULL) free(rec_varids);

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
    opt.chk      = 1024; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1;    /* test burst-buffering feature */
    opt.mod      = 1;    /* test independent data mode */
    opt.hdr_diff = 1;    /* run ncmpidiff for file header only */
    opt.var_diff = 1;    /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "re-enter define mode", opt, test_io);

    MPI_Finalize();

    return err;
}
