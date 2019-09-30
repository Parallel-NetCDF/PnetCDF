/*
   This program test reading a compressed netcdf 4 file
   The data is written using NetCDF and read by PnetCDF
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>

#include <netcdf.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NZ 4
#define NY 4
#define NX 4

static int create_nc4(char *filename)
{
    int i, err, nerrs=0;
    int ncid, varid, dimids[3], *buf;
    int shuffle, deflate, deflate_level;

    /* create a new file, if already exists, clobber it */
    err = nc_create(filename, NC_NETCDF4|NC_CLOBBER, &ncid); CHECK_ERR

    /* define dimensions */
    err = nc_def_dim(ncid, "z", NZ, &dimids[0]); CHECK_ERR
    err = nc_def_dim(ncid, "y", NY, &dimids[1]); CHECK_ERR
    err = nc_def_dim(ncid, "x", NX, &dimids[2]); CHECK_ERR

    /* define variable */
    err = nc_def_var(ncid, "var", NC_INT, 3, dimids, &varid); CHECK_ERR

    /* Set the compression settings for the variable */
    shuffle = NC_SHUFFLE;
    deflate = 1;
    deflate_level = 5;
    err = nc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level); CHECK_ERR

    /* exit define mode */
    err = nc_enddef(ncid); CHECK_ERR

    /* initialize buffer contents */
    buf = (int*) malloc(NZ*NY*NZ * sizeof(int));
    for (i=0; i<NZ*NY*NZ; i++) buf[i] = i;

    /* write the entire variable */
    err = nc_put_var_int(ncid, varid, buf); CHECK_ERR

    free(buf);

    /* flush write data to file system */
    err = nc_sync(ncid); CHECK_ERR

    err = nc_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char **argv) {
    char filename[512];
    int i, err, nerrs=0, rank, np;
    int ncid, ndims, nvars, varid, dimids[3], *buf;
    MPI_Offset bufLen;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for reading compressed NetCDF4 file", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    /* rank 0 creates a NETCDF4 file */
    if (rank == 0) {
        /* remove the file system type prefix name if there is any.
         * For example, when filename = "lustre:/home/foo/testfile.nc", remove
         * "lustre:" to make path = "/home/foo/testfile.nc" in open() below
         */
        char *path = strchr(filename, ':');
        if (path == NULL) path = filename; /* no prefix */
        else              path++;

        err = create_nc4(path);
        if (err) {
            nerrs++;
            goto fn_exit;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* Read the NetCDF4 file with PnetCDF APIs */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    err = ncmpi_inq(ncid, &ndims, &nvars, NULL, NULL); CHECK_ERR

    /* check if number of dimensions defined in the file is expected */
    if (ndims != 3) {
        printf("Error: expect ndims=3 but got %d\n", ndims);
        nerrs++;
        goto fn_exit;
    }
    /* check if number of variables defined in the file is expected */
    if (nvars != 1) {
        printf("Error: expect nvars=1 but got %d\n", nvars);
        nerrs++;
        goto fn_exit;
    }

    /* obtain variable ID */
    err = ncmpi_inq_varid(ncid, "var", &varid); CHECK_ERR

    /* check if number of dimensions of the variable is expected */
    err = ncmpi_inq_varndims(ncid, varid, &ndims); CHECK_ERR
    if (ndims != 3) {
        printf("Error: expect variable's ndims=3 but got %d\n", ndims);
        nerrs++;
        goto fn_exit;
    }

    /* obtain dimension IDs of the variable and their lengths */
    err = ncmpi_inq_vardimid(ncid, varid, dimids); CHECK_ERR
    bufLen = 1;
    for (i=0; i<3; i++) {
        MPI_Offset dimLen;
        err = ncmpi_inq_dimlen(ncid, dimids[i], &dimLen); CHECK_ERR
        bufLen *= dimLen;
    }

    buf = (int*) malloc(bufLen * sizeof(int));

    /* read the entire variable */
    err = ncmpi_get_var_int_all(ncid, varid, buf); CHECK_ERR

    /* check the contents of variable */
    for (i=0; i<bufLen; i++)
        if (buf[i] != i) {
            printf("Error: expect buf[%d]=%d, but got %d\n",i,i,buf[i]);
            nerrs++;
            break;
        }

    free(buf);

    /* Close file */
    ncmpi_close(ncid);

    /* check if there is any malloc residue */
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

fn_exit:
   MPI_Finalize();
   return (nerrs > 0);
}
