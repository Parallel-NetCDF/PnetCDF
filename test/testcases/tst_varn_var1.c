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

#define NY 4
#define NX 4
#define NDIMS 2

int main(int argc, char** argv)
{
    char filename[256];
    int i, j, k, rank, nprocs, err, nerrs=0;
    int ncid, cmode, varid[2], dimid[2], nreqs, req, *buf;
    MPI_Offset **starts=NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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
        sprintf(cmd_str, "*** TESTING C   %s for ncmpi_put_varn_int_all() ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    buf = (int*) malloc(sizeof(int) * NY * NX);

    nreqs = NY * NX * nprocs;
    starts    = (MPI_Offset**) malloc(sizeof(MPI_Offset*) * nreqs);
    starts[0] = (MPI_Offset*)  calloc(nreqs * NDIMS, sizeof(MPI_Offset));
    for (i=1; i<nreqs; i++)
        starts[i] = starts[i-1] + NDIMS;

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
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
        err = ncmpi_put_varn_int_all(ncid, varid[k], nreqs, starts, NULL, buf);
        CHECK_ERR

        /* read back */
        for (i=0; i<nreqs; i++) buf[i] = -99;

        /* read using varn API equivalent to multiple get_var1 */
        err = ncmpi_get_varn_int_all(ncid, varid[k], nreqs, starts, NULL, buf);
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

        err = ncmpi_wait_all(ncid, 1, &req, NULL);
        CHECK_ERR

        /* read back */
        for (i=0; i<nreqs; i++) buf[i] = -99;

        /* read using varn API equivalent to multiple get_var1 */
        err = ncmpi_iget_varn_int(ncid, varid[k], nreqs, starts, NULL, buf, &req);
        CHECK_ERR

        err = ncmpi_wait_all(ncid, 1, &req, NULL);
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

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

