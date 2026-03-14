/*********************************************************************
 *
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example tests a single call of ncmpi_put_varn_int_all() to write a
 * sequence of requests with arbitrary array indices, all with length == 1.
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicc -O2 -o tst_varn_var1 tst_varn_var1.c -lpnetcdf
 *    % mpiexec -n 4 ./tst_varn_var1 /pvfs2/wkliao/testfile.nc
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *             Y = 4 ;
 *             X = 10 ;
 *             time = UNLIMITED ; // (4 currently)
 *    variables:
 *             int fix_var(Y, X) ;
 *             int rec_var(time, X) ;
 *    data:
 *
 *     fix_var =
 *       0, _, -1, _, -2, _, -3, _, -4, _, 0, _, -1, _, -2, _, -3, _, -4, _,
 *       0, _, -1, _, -2, _, -3, _, -4, _, 0, _, -1, _, -2, _, -3, _, -4, _,
 *       0, _, -1, _, -2, _, -3, _, -4, _, 0, _, -1, _, -2, _, -3, _, -4, _,
 *       0, _, -1, _, -2, _, -3, _, -4, _, 0, _, -1, _, -2, _, -3, _, -4, _ ;
 *
 *     rec_var =
 *       0, _, -1, _, -2, _, -3, _, -4, _, 0, _, -1, _, -2, _, -3, _, -4, _,
 *       0, _, -1, _, -2, _, -3, _, -4, _, 0, _, -1, _, -2, _, -3, _, -4, _,
 *       0, _, -1, _, -2, _, -3, _, -4, _, 0, _, -1, _, -2, _, -3, _, -4, _,
 *       0, _, -1, _, -2, _, -3, _, -4, _, 0, _, -1, _, -2, _, -3, _, -4, _ ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), memset() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NY 40
#define NX 40
#define NDIMS 2

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int i, j, k, rank, nprocs, err, nerrs=0;
    int ncid, varid[2], dimid[2], nreqs, req, *buf;
    MPI_Offset **starts=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    buf = (int*) malloc(sizeof(int) * NY * NX);

    nreqs = NY * NX * nprocs;
    starts    = (MPI_Offset**) malloc(sizeof(MPI_Offset*) * nreqs);
    starts[0] = (MPI_Offset*)  calloc(nreqs * NDIMS, sizeof(MPI_Offset));
    for (i=1; i<nreqs; i++)
        starts[i] = starts[i-1] + NDIMS;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* create a global array of size NY * NX */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX * nprocs, &dimid[1]);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "fix_var", NC_INT, NDIMS, dimid, &varid[0]);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "rec_var", NC_INT, NDIMS, dimid, &varid[1]);
    CHECK_ERR

    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    for (i=0; i<NY; i++) {
        err = ncmpi_fill_var_rec(ncid, varid[1], i);
        CHECK_ERR
    }

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* test using blocking varn API with counts == NULL */

    for (k=0; k<2; k++) {
        /* write using put_varn API equivalent to multiple put_var1 */
        nreqs = 0;
        for (i=0; i<NY; i++) {
            for (j=0; j<2; j++) {
                starts[nreqs][0] = i;
                starts[nreqs][1] = j * 2 * nprocs + 2 * rank;
                buf[nreqs] = -rank;
                nreqs++;
            }
        }
        if (coll_io)
            err = ncmpi_put_varn_int_all(ncid, varid[k], nreqs, starts, NULL, buf);
        else
            err = ncmpi_put_varn_int(ncid, varid[k], nreqs, starts, NULL, buf);
        CHECK_ERR

        /* file sync before reading */
        err = ncmpi_sync(ncid);
        CHECK_ERR
        MPI_Barrier(MPI_COMM_WORLD);


        /* read back */
        for (i=0; i<nreqs; i++) buf[i] = -99;

        /* read using varn API equivalent to multiple get_var1 */
        if (coll_io)
            err = ncmpi_get_varn_int_all(ncid, varid[k], nreqs, starts, NULL, buf);
        else
            err = ncmpi_get_varn_int(ncid, varid[k], nreqs, starts, NULL, buf);
        CHECK_ERR

        /* check read buf contents */
        for (i=0; i<nreqs; i++)
            if (buf[i] != -rank) {
                printf("Error at line %d in %s: expecting buf[%d]=%d but got %d\n",
                        __LINE__,__FILE__,i,-rank,buf[i]);
                nerrs++;
                goto err_out;
            }
    }

    /* test using nonblocking varn API with counts == NULL */

    for (k=0; k<2; k++) {
        /* write using iput_varn API equivalent to multiple iput_var1 */
        nreqs = 0;
        for (i=0; i<NY; i++) {
            for (j=0; j<2; j++) {
                starts[nreqs][0] = i;
                starts[nreqs][1] = j * 2 * nprocs + 2 * rank;
                buf[nreqs] = -rank;
                nreqs++;
            }
        }
        err = ncmpi_iput_varn_int(ncid, varid[k], nreqs, starts, NULL, buf, &req);
        CHECK_ERR

        if (coll_io)
            err = ncmpi_wait_all(ncid, 1, &req, NULL);
        else
            err = ncmpi_wait(ncid, 1, &req, NULL);
        CHECK_ERR

        /* file sync before reading */
        err = ncmpi_sync(ncid);
        CHECK_ERR
        MPI_Barrier(MPI_COMM_WORLD);

        /* read back */
        for (i=0; i<nreqs; i++) buf[i] = -99;

        /* read using varn API equivalent to multiple get_var1 */
        err = ncmpi_iget_varn_int(ncid, varid[k], nreqs, starts, NULL, buf, &req);
        CHECK_ERR

        if (coll_io)
            err = ncmpi_wait_all(ncid, 1, &req, NULL);
        else
            err = ncmpi_wait(ncid, 1, &req, NULL);
        CHECK_ERR

        /* check read buf contents */
        for (i=0; i<nreqs; i++)
            if (buf[i] != -rank) {
                printf("Error at line %d in %s: expecting buf[%d]=%d but got %d\n",
                        __LINE__,__FILE__,i,-rank,buf[i]);
                nerrs++;
                goto err_out;
            }
    }

err_out:
    err = ncmpi_close(ncid);
    CHECK_ERR

    free(buf);
    free(starts[0]);
    free(starts);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "ncmpi_put_varn_int_all()", opt, test_io);

    MPI_Finalize();

    return err;
}
