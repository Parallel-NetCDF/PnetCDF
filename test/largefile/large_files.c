/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  This program tests writing 2 record variables into a file of size larger
 *  than 4 GiB. First variable is a 4D array of total size > 4 GiB and the
 *  second is a small 2D array. The contents of both variables are read back
 *  to check whether they are written correctly.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NUMRECS 1
#define I_LEN 4104
#define J_LEN 1023
#define K_LEN 1023
#define N_LEN 2

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int err, nerrs=0, rank, ncid, varid, x_id;
    int n, rec, i, j, k, dims[4];
    MPI_Offset start[4] = {0, 0, 0, 0};
    MPI_Offset count[4] = {1, 1, J_LEN, K_LEN};

    /* I/O buffers */
    signed char *buf;

    // printf("\n*** Testing large files, slowly.\n");
    // printf("*** Creating large file %s...", out_path);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank > 0) goto err_out;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_SELF, out_path, NC_CLOBBER, info, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line %d in %s: (%s)\n", __LINE__,__FILE__,ncmpi_strerrno(err));
        return 1;
    }

    buf = (signed char*) malloc(J_LEN * K_LEN);

    /* define a 4D record variable */
    err = ncmpi_def_dim(ncid, "rec", NC_UNLIMITED, &dims[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "i",   I_LEN,        &dims[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "j",   J_LEN,        &dims[2]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "k",   K_LEN,        &dims[3]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_BYTE, 4, dims, &varid);
    CHECK_ERR

    /* define a 2D record variable */
    err = ncmpi_def_dim(ncid, "n", N_LEN, &dims[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "x", NC_BYTE, 2, dims, &x_id); CHECK_ERR

    /* don't initialize variables with fill values */
    err = ncmpi_set_fill(ncid, NC_NOFILL, 0); CHECK_ERR

    /* leave define mode */
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* write var */
    n = 0;
    for (rec=0; rec<NUMRECS; rec++) {
        start[0] = rec;
        for (i=0; i<I_LEN; i++) { /* initialize write buf */
            for (j=0; j<J_LEN; j++) {
                for (k=0; k<K_LEN; k++) {
                    buf[j*K_LEN+k] = (signed char)n;
                    n = (n == 127) ? 0 : (n+1);
                }
            }
            start[1] = i;
            count[1] = 1;
            err = ncmpi_put_vara_schar_all(ncid, varid, start, count, buf);
            CHECK_ERR
        }
        buf[0] = 42; /* initialize write buf */
        buf[1] = 21;
        start[1] = 0;
        count[1] = N_LEN;
        err = ncmpi_put_vara_schar_all(ncid, x_id, start, count, buf);
        CHECK_ERR
    }
    err = ncmpi_close(ncid); CHECK_ERR

    // printf("ok\n");
    // printf("*** Reading large file %s...", out_path);

    err = ncmpi_open(MPI_COMM_SELF, out_path, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* read variables and check their contents */
    n = 0;
    for (rec=0; rec<NUMRECS; rec++) {
        start[0] = rec;
        for (i=0; i<I_LEN; i++) {
            for (j=0; j<J_LEN; j++) /* set read buf to all -1 */
                for (k=0; k<K_LEN; k++)
                    buf[j*K_LEN+k] = (signed char) -1;

            start[1] = i;
            count[1] = 1;
            err = ncmpi_get_vara_schar_all(ncid, varid, start, count, buf);
            CHECK_ERR

            for (j=0; j<J_LEN; j++) {
                for (k=0; k<K_LEN; k++) {
                    if (buf[j*K_LEN+k] != (signed char) n % 128) {
                        printf("Error on read, var[%d, %d, %d, %d] = %d wrong, should be %d !\n",
                               rec, i, j, k, buf[j*K_LEN+k], (signed char)n);
                        nerrs++;
                        goto err_out;
                    }
                    n = (n == 127) ? 0 : (n+1);
                }
            }
        }

        /* read variable x and check its contents */
        buf[0] = buf[1] = -1;
        start[1] = 0;
        count[1] = N_LEN;
        ncmpi_get_vara_schar_all(ncid, x_id, start, count, buf); CHECK_ERR
        if (buf[0] != 42 || buf[1] != 21) {
            printf("Error on read, x[] = %d, %d\n", buf[0], buf[1]);
            nerrs++;
        }
    }
    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);

err_out:
    MPI_Barrier(MPI_COMM_WORLD);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 0; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "> 4 GiB file", opt, test_io);

    MPI_Finalize();

    return err;
}
